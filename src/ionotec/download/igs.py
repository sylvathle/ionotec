from pathlib import Path
from datetime import timedelta

import re
import json
import time
import tempfile
import subprocess
import requests
from concurrent.futures import ThreadPoolExecutor

from ionotec.download.extract import extract_and_remove_gz

# NASA CDDIS archive of daily IGS RINEX files.
#   .../daily/YYYY/DDD/YYd/   RINEX 3/4 long names (.crx.gz / .rnx.gz) and
#                             RINEX 2 Hatanaka (ssssDDD0.YYd.gz or .Z)
#   .../daily/YYYY/DDD/YYo/   RINEX 2 plain observation (ssssDDD0.YYo.gz or .Z)
# Authentication: Earthdata credentials in ~/.netrc
#   machine urs.earthdata.nasa.gov login <user> password <pass>
BASE_URL = "https://cddis.nasa.gov/archive/gnss/data/daily"

# True  -> use curl (needs curl and a working ~/.netrc); requests is only
#          used if curl is not installed. Default: plain requests gets
#          "Connection reset by peer" on the CDDIS daily listings.
# False -> try requests first and fall back to curl if it fails.
USE_CURL = True

# Successful daily listings are saved here so each day is fetched only once
# (None -> no disk cache). The file name contains the directory kind (d/o).
CACHE_DIR = Path.home() / ".cache" / "ionotec" / "igs"

DEFAULT_WORKERS = 4   # parallel connections to CDDIS

# RINEX 3/4 long name, daily observation, Hatanaka or plain:
#   ABMF00GLP_R_20240960000_01D_30S_MO.crx.gz
_RE_R3 = re.compile(
    r"^([A-Z0-9]{9})_[A-Z]_\d{11}_\d{2}[A-Z]_\d{2}[A-Z]_MO\.(?:crx|rnx)\.gz$",
    re.IGNORECASE,
)
# RINEX 2 short name: Hatanaka (d) or plain (o), gzip or compress:
#   abmf0960.24d.gz   abmf0960.24d.Z   abmf0960.24o.gz   abmf0960.24o.Z
_RE_R2 = re.compile(
    r"^([A-Z0-9]{4})\d{3}[0-9A-X]\.\d{2}([do])\.(?:gz|Z)$",
    re.IGNORECASE,
)

# Same patterns, to find file names anywhere in a page
_RE_FIND_R3 = re.compile(
    r"(?<![A-Za-z0-9])[A-Za-z0-9]{9}_[A-Za-z]_\d{11}_\d{2}[A-Za-z]_\d{2}[A-Za-z]"
    r"_MO\.(?:crx|rnx)\.gz(?![A-Za-z0-9_.])",
    re.IGNORECASE,
)
_RE_FIND_R2 = re.compile(
    r"(?<![A-Za-z0-9])[A-Za-z0-9]{4}\d{3}[0-9A-Xa-x]\.\d{2}[do]\.(?:gz|Z)"
    r"(?![A-Za-z0-9_.])",
    re.IGNORECASE,
)

# One listing per day and directory kind is enough for the whole process
_LISTING_CACHE = {}


# ----------------------------------------------------------------------
# Transport
# ----------------------------------------------------------------------

def _day_url(d, kind="d"):
    """Directory holding the daily files of date d: 'd' = YYd (RINEX 3/4
    long names + RINEX 2 Hatanaka), 'o' = YYo (RINEX 2 plain)."""
    yy = d.year % 100
    return f"{BASE_URL}/{d.year}/{d.timetuple().tm_yday:03d}/{yy:02d}{kind}/"


def _new_session():
    # Plain session with default headers (this is what works with CDDIS)
    return requests.Session()


def _fetch_requests(session, url, dest=None, max_retry_time=60):
    """
    GET url like get_dcb_from_cddis does: default headers, stream=True,
    timeout=30, retry every second for up to max_retry_time s, on any
    exception.

    dest=None -> returns the body (bytes); otherwise writes to dest and
    returns True. Returns None on 404.
    """
    start_time = time.time()

    while True:
        try:
            with session.get(url, stream=True, timeout=30) as r:

                if r.status_code == 404:
                    return None

                r.raise_for_status()

                if dest is None:
                    return r.content

                with open(dest, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                return True

        except Exception as e:
            if dest is not None and Path(dest).exists():
                Path(dest).unlink()

            elapsed = time.time() - start_time

            if elapsed >= max_retry_time:
                print(f"Failed {url} after {elapsed:.1f} s: {e}")
                raise

            print(f"Request failed ({e}). Retrying in 1 second...")
            time.sleep(1)


def _fetch_curl(url, dest=None, max_attempts=6):
    """
    Same as _fetch_requests but through curl (-n: ~/.netrc, -L: follow the
    Earthdata redirects). The cookie engine is on (-b "") but no cookie
    file is read or written.
    Network failures (e.g. 'Connection reset by peer') are retried here,
    with a growing pause, because curl's own --retry ignores them.
    Files are resumed after a reset (-C -).
    Returns bytes (dest=None) / True, or None on 404.
    """
    if dest is None:
        tmp = tempfile.NamedTemporaryFile(delete=False)
        tmp.close()
        out = Path(tmp.name)
    else:
        out = Path(dest)

    cmd = [
        "curl", "-sS", "-n", "-L", "-b", "",
        "-o", str(out), "-w", "%{http_code}", url,
    ]

    if dest is not None:
        cmd[1:1] = ["-C", "-"]

    try:
        for attempt in range(max_attempts):

            r = subprocess.run(cmd, capture_output=True, text=True,
                               timeout=900)

            code = int(r.stdout.strip()[-3:] or 0) if r.returncode == 0 else 0

            if r.returncode == 0 and code == 404:
                return None

            # 416: the file was already complete
            if dest is not None and r.returncode == 0 and code == 416:
                return True

            # (33) server cannot resume -> restart from scratch
            if r.returncode == 33 and dest is not None:
                out.unlink(missing_ok=True)
                continue

            if r.returncode == 0 and 0 < code < 400:
                if dest is None:
                    return out.read_bytes()
                return True

            msg = r.stderr.strip() or f"HTTP {code}"

            if attempt == max_attempts - 1:
                raise ConnectionError(f"curl: {msg}")

            wait = 10 * (attempt + 1)
            print(f"WARNING {url}: {msg} -> retry in {wait} s")
            time.sleep(wait)

    finally:
        if dest is None:
            out.unlink(missing_ok=True)
        elif out.exists() and out.stat().st_size == 0:
            out.unlink(missing_ok=True)


def _fetch(session, url, dest=None):
    """curl first (USE_CURL), requests as fallback, or the reverse."""
    if USE_CURL:
        try:
            return _fetch_curl(url, dest)
        except FileNotFoundError:
            print("WARNING curl not found, using requests")
            return _fetch_requests(session, url, dest)

    try:
        return _fetch_requests(session, url, dest)
    except Exception as e:
        print(f"WARNING requests failed for {url}: {e}")
        print("        falling back to curl")

    return _fetch_curl(url, dest)


# ----------------------------------------------------------------------
# Listings
# ----------------------------------------------------------------------

def _list_day(session, d, kind="d"):
    """
    Observation file names available on CDDIS for date d in the YYd
    (kind="d": long names + RINEX 2 Hatanaka) or YYo (kind="o": RINEX 2
    plain) directory. Cached in memory and on disk.
    Returns None if the directory does not exist.
    """
    day = d.date() if hasattr(d, "date") else d
    key = (day, kind)

    if key in _LISTING_CACHE:
        return _LISTING_CACHE[key]

    cache_file = None
    if CACHE_DIR is not None:
        cache_file = Path(CACHE_DIR) / f"{day:%Y%j}_{kind}.json"
        if cache_file.exists():
            try:
                _LISTING_CACHE[key] = json.loads(cache_file.read_text())
                return _LISTING_CACHE[key]
            except Exception:
                cache_file.unlink(missing_ok=True)

    body = _fetch(session, _day_url(d, kind))

    if body is None:
        return None

    text = body.decode("utf-8", "replace")

    # Scan the whole page for file names, whatever the HTML looks like
    found = set(_RE_FIND_R2.findall(text))
    if kind == "d":
        found |= set(_RE_FIND_R3.findall(text))

    # keep only what belongs in this directory kind
    names = {
        n for n in found
        if _RE_R3.match(n) or (_RE_R2.match(n) and _RE_R2.match(n).group(2).lower() == kind)
    }

    if not names:
        debug = Path(tempfile.gettempdir()) / f"igs_listing_{day:%Y%j}_{kind}.html"
        debug.write_text(text)
        print(
            f"WARNING no observation file found in {_day_url(d, kind)} "
            f"({len(body)} bytes). Page saved to {debug}"
        )
        print("        first 300 characters:", repr(text[:300]))

    _LISTING_CACHE[key] = sorted(names)

    if cache_file is not None and names:
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        cache_file.write_text(json.dumps(_LISTING_CACHE[key]))

    return _LISTING_CACHE[key]


def _days(start, end):
    days = []
    d = start
    while d <= end:
        days.append(d)
        d += timedelta(days=1)
    return days


def _list_days(days, workers=DEFAULT_WORKERS, kind="d"):
    """{day: list of file names, or None if missing/failed}, in parallel."""

    def one(d):
        try:
            with _new_session() as session:
                return _list_day(session, d, kind)
        except Exception as e:
            print(f"ERROR {d}: {e}")
            return None

    with ThreadPoolExecutor(max_workers=max(1, workers)) as ex:
        return dict(zip(days, ex.map(one, days)))


def _station_ids(fname):
    """(full_id, site4) in lower case; full_id is None for RINEX 2 names."""
    m = _RE_R3.match(fname)
    if m:
        full = m.group(1).lower()
        return full, full[:4]

    m = _RE_R2.match(fname)
    if m:
        return None, m.group(1).lower()

    return None, None


def _is_short(fname, kind):
    m = _RE_R2.match(fname)
    return bool(m) and m.group(2).lower() == kind


# ----------------------------------------------------------------------
# Which file to take: RINEX 3/4 first, RINEX 2 only if there is none
# ----------------------------------------------------------------------

def _pick(station, names_d, names_o=None):
    """
    Priority list [(kind, file name)] for a 4-character station, from the
    listings of the YYd (names_d) and YYo (names_o) directories.
    ONE file is chosen:
        1. RINEX 3/4 long name (YYd)
        2. RINEX 2 Hatanaka    (YYd)
        3. RINEX 2 plain       (YYo)
    Empty list if the station has none of them.
    """
    longs = sorted(
        n for n in names_d or []
        if _station_ids(n)[0] and _station_ids(n)[1] == station
    )
    if longs:
        return [("d", longs[0])]

    shorts_d = sorted(
        n for n in names_d or []
        if _is_short(n, "d") and _station_ids(n)[1] == station
    )
    if shorts_d:
        return [("d", shorts_d[0])]

    shorts_o = sorted(
        n for n in names_o or []
        if _is_short(n, "o") and _station_ids(n)[1] == station
    )
    if shorts_o:
        return [("o", shorts_o[0])]

    return []


def _blind_candidates(station, d, long_name=None):
    """Names built without any listing, in priority order. The RINEX 3/4
    long name is only known for a 9-character id (the sampling is assumed
    to be 30 s)."""
    site = station[:4]
    doy = d.timetuple().tm_yday
    yy = d.year % 100

    out = []

    if long_name:
        out.append(("d", long_name))

    for kind in ("d", "o"):
        for ext in ("gz", "Z"):
            out.append((kind, f"{site}{doy:03d}0.{yy:02d}{kind}.{ext}"))

    return out


def stations_available_igs(start, end, full_id=False,
                           workers=DEFAULT_WORKERS, include_plain=False):
    """
    Return the union of stations having daily RINEX observation files
    on the NASA CDDIS archive during the requested period: RINEX 3/4
    (long names) and, for stations that only have RINEX 2, RINEX 2
    (ssssDDD0.YYd.gz / .Z).

    Archive:
        https://cddis.nasa.gov/archive/gnss/data/daily/YYYY/DDD/YYd/

    Parameters
    ----------
    full_id : bool
        False -> 4-character lower-case names (like stations_available_cors).
        True  -> 9-character upper-case ids (e.g. 'ABMF00GLP') for stations
                 with a long name; RINEX 2-only stations have only 4
                 characters (e.g. 'PGC5').
    workers : int
        Daily listings fetched in parallel.
    include_plain : bool
        Also read the YYo directory (plain RINEX 2, ssssDDD0.YYo.gz/.Z) to
        find stations that have neither RINEX 3/4 nor Hatanaka RINEX 2.
        One more (slow) listing per day, so off by default.
    """

    days = _days(start, end)

    list_stations = set()
    short_only = set()

    def add(names):
        for fname in names:
            full, site = _station_ids(fname)

            if full_id and full:
                list_stations.add(full.upper())
            elif full_id:
                short_only.add(site.upper())
            else:
                list_stations.add(site)

    for d, names in _list_days(days, workers, "d").items():
        if names is None:
            print(f"NONE  no listing for {d}")
            continue
        add(names)

    if include_plain:
        for d, names in _list_days(days, workers, "o").items():
            if names is None:
                print(f"NONE  no YYo listing for {d}")
                continue
            add(names)

    if full_id:
        # a station seen with a 9-character id does not need its 4-character entry
        with_long = {s[:4] for s in list_stations}
        list_stations |= {s for s in short_only if s not in with_long}

    return sorted(list_stations)


# ----------------------------------------------------------------------
# Download
# ----------------------------------------------------------------------

def _extract(path):
    """Extract a .gz (extract_and_remove_gz) or a Unix-compressed .Z
    (gzip -dc reads it); the archive is removed. Returns the extracted Path."""
    path = Path(path)

    if path.suffix != ".Z":
        return extract_and_remove_gz(path)

    out = path.with_suffix("")

    try:
        with open(out, "wb") as f:
            subprocess.run(["gzip", "-dc", str(path)], stdout=f, check=True,
                           stderr=subprocess.PIPE)

    except FileNotFoundError:
        out.unlink(missing_ok=True)
        raise RuntimeError("gzip was not found in the PATH")

    except subprocess.CalledProcessError as e:
        out.unlink(missing_ok=True)
        raise RuntimeError(f"gzip failed on {path.name}: {e.stderr!r}")

    path.unlink()

    return out


def _rnx2crx(rnx_path):
    """
    Hatanaka-compress a RINEX file: `rnx2crx -f file.rnx` -> file.crx.
    The .rnx file is removed afterwards. Returns the .crx Path.
    """
    rnx_path = Path(rnx_path)
    crx_path = rnx_path.with_suffix(".crx")

    try:
        subprocess.run(
            ["rnx2crx", "-f", str(rnx_path)],
            check=True,
            capture_output=True,
            text=True,
        )

    except FileNotFoundError:
        raise RuntimeError("rnx2crx was not found in the PATH")

    except subprocess.CalledProcessError as e:
        crx_path.unlink(missing_ok=True)
        raise RuntimeError(f"rnx2crx failed on {rnx_path.name}: {e.stderr}")

    if not crx_path.exists():
        raise RuntimeError(f"rnx2crx did not create {crx_path.name}")

    rnx_path.unlink()

    return crx_path


def _crx2rnx(crx_path):
    """
    Un-compress a RINEX 2 Hatanaka file: `crx2rnx -f file.YYd` ->
    file.YYo. The .YYd file is removed afterwards. Returns the .YYo Path.
    """
    crx_path = Path(crx_path)
    rnx_path = crx_path.with_suffix(crx_path.suffix[:-1] + "o")   # .24d -> .24o

    try:
        subprocess.run(
            ["crx2rnx", "-f", str(crx_path)],
            check=True,
            capture_output=True,
            text=True,
        )

    except FileNotFoundError:
        raise RuntimeError("crx2rnx was not found in the PATH")

    except subprocess.CalledProcessError as e:
        rnx_path.unlink(missing_ok=True)
        raise RuntimeError(f"crx2rnx failed on {crx_path.name}: {e.stderr}")

    if not rnx_path.exists():
        raise RuntimeError(f"crx2rnx did not create {rnx_path.name}")

    crx_path.unlink()

    return rnx_path


def _is_r2_hatanaka(path):
    """True for a RINEX 2 Hatanaka file name (ssssDDD0.YYd)."""
    return re.fullmatch(r"\.\d{2}d", Path(path).suffix, re.IGNORECASE) is not None


def _final_path(gz_path):
    """Where the downloaded archive ends up:
        X.crx.gz -> X.crx           (RINEX 3/4 Hatanaka, kept)
        X.rnx.gz -> X.crx           (after rnx2crx)
        X.24d.gz / .Z -> X.24o      (RINEX 2 Hatanaka, after crx2rnx)
        X.24o.gz / .Z -> X.24o      (RINEX 2 plain, kept)"""
    p = Path(gz_path).with_suffix("")

    if p.suffix.lower() == ".rnx":
        return p.with_suffix(".crx")

    if _is_r2_hatanaka(p):
        return p.with_suffix(p.suffix[:-1] + "o")

    return p


def _download_one(d, candidates, output_dir):
    """
    Download + extract ONE file: the first candidate that exists
    (404 -> next one). A plain RINEX 3/4 .rnx is converted to Hatanaka
    (.crx) with rnx2crx; a RINEX 2 Hatanaka file (.YYd) is converted to
    plain RINEX 2 (.YYo) with crx2rnx. A network/extraction error stops
    the search (no silent fallback to an older RINEX version).
    Returns the final Path, or None.
    """
    for kind, fname in candidates:

        gz_path = output_dir / fname
        extracted_path = gz_path.with_suffix("")
        final_path = _final_path(gz_path)

        if final_path.exists():
            print(f"SKIP  {final_path.name}")
            return final_path

        try:
            print(f"GET   {fname}")

            with _new_session() as session:
                ok = _fetch(session, _day_url(d, kind) + fname, dest=gz_path)

            if ok is None:                         # 404: try the next one
                print(f"NONE  {fname}")
                continue

            extracted = _extract(gz_path)

            if extracted.suffix.lower() == ".rnx":
                extracted = _rnx2crx(extracted)          # RINEX 3/4 -> .crx

            #elif _is_r2_hatanaka(extracted):
            #    extracted = _crx2rnx(extracted)          # RINEX 2 .YYd -> .YYo

            return extracted

        except Exception as e:
            # leave nothing half-done (.gz/.Z, extracted file, converted file)
            for p in (gz_path, extracted_path, final_path):
                p.unlink(missing_ok=True)

            print(f"ERROR {fname}: {e}")
            return None

    return None


def download_igs_stations(stations, start, end, output_dir=".",
                          workers=DEFAULT_WORKERS):
    """
    Download several stations at once, all files in parallel.

    stations: iterable of 4- or 9-character names.
    One file per station and day: RINEX 3/4 if the station has it,
    RINEX 2 (.YYd, then .YYo) only if it does not.
    Returns {station: [Path, ...]}.
    """

    stations = [s.lower() for s in stations]

    for s in stations:
        if len(s) not in (4, 9):
            raise ValueError(f"bad station name: {s!r}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    days = _days(start, end)

    # --- listings: YYd for every day, YYo only where some station needs it
    listings_d, listings_o = {}, {}

    if any(len(s) == 4 for s in stations):
        listings_d = _list_days(days, workers, "d")

        need_o = {
            d for d in days
            for s in stations
            if len(s) == 4
            and listings_d.get(d) is not None
            and not _pick(s, listings_d[d])
        }

        if need_o:
            listings_o = _list_days(sorted(need_o), workers, "o")

    # --- plan: (station, day, candidates)
    tasks = []

    for st in stations:
        for d in days:

            doy = d.timetuple().tm_yday

            if len(st) == 9:
                # No listing: the long name, then the RINEX 2 names
                long_name = f"{st.upper()}_R_{d.year}{doy:03d}0000_01D_30S_MO.crx.gz"
                tasks.append((st, d, _blind_candidates(st, d, long_name)))
                continue

            names_d = listings_d.get(d)

            if names_d is None:
                print(f"NONE  {st} on {d} (no listing, DOY {doy:03d})")
                continue

            cands = _pick(st, names_d, listings_o.get(d))

            if cands:
                tasks.append((st, d, cands))
            else:
                print(f"NONE  {st} on {d} (DOY {doy:03d})")

    with ThreadPoolExecutor(max_workers=max(1, workers)) as ex:
        results = list(ex.map(
            lambda t: _download_one(t[1], t[2], output_dir), tasks
        ))

    out = {st: [] for st in stations}
    for (st, _, _), path in zip(tasks, results):
        if path is not None:
            out[st].append(path)

    print()
    print(f"Files available locally: {sum(len(v) for v in out.values())}")

    return out


def download_igs(station, start, end, output_dir=".",
                 workers=DEFAULT_WORKERS):
    """
    Download IGS daily RINEX observation files from NASA CDDIS.

    station: 4-character site name or 9-character id (case-insensitive).
    One file per station and day, in this priority:
        1. RINEX 3/4 long name (SSSS00CCC_R_..._MO.crx.gz / .rnx.gz)
        2. RINEX 2 Hatanaka (ssssDDD0.YYd.gz / .Z)
        3. RINEX 2 plain    (ssssDDD0.YYo.gz / .Z)
    A plain RINEX 3/4 .rnx is converted to Hatanaka (.crx) with `rnx2crx`;
    a RINEX 2 Hatanaka file (.YYd) is stored uncompressed (.YYo), converted
    with `crx2rnx`. Both tools must be in the PATH.
    Days are downloaded in parallel (`workers` connections).

    Returns
    -------
    list of Path
        Locally available RINEX files.
    """

    station = station.lower()

    return download_igs_stations(
        [station], start, end, output_dir, workers
    )[station]


def download_file(station, date, output_dir="."):
    """
    Download ONE daily RINEX observation file of a station from CDDIS.

    station : 4-character site name (case-insensitive)
    date    : datetime / date of the day wanted

    Priority (the first one that exists is downloaded, and only that one):
        1. RINEX 3/4 long name  SSSS00CCC_R_YYYYDDD0000_01D_30S_MO.crx.gz
           (or .rnx.gz), YYd directory
        2. RINEX 2 Hatanaka     ssssDDD0.YYd.gz / .Z, YYd directory
        3. RINEX 2 plain        ssssDDD0.YYo.gz / .Z, YYo directory

    The next version is only tried when the previous one does not exist
    (404 / not in the listing). A network error stops the search instead
    of silently falling back to an older RINEX version.

    The archive is extracted (.gz with extract_and_remove_gz, .Z with
    gzip -dc) and removed. A plain RINEX 3/4 .rnx is converted to Hatanaka
    (.crx) with `rnx2crx`; a RINEX 2 Hatanaka file (.YYd) is converted to
    plain RINEX 2 (.YYo) with `crx2rnx`. Both tools must be in the PATH and
    the intermediate file is removed.

    Returns
    -------
    Path of the final file: .../PGC500CAN_R_20240010000_01D_30S_MO.crx for
    RINEX 3/4, .../pgc50010.24o for RINEX 2, or None if no file was found /
    on error.
    """

    station = station.lower()

    if len(station) != 4:
        raise ValueError("station must be a 4-character station name")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    doy = date.timetuple().tm_yday

    with _new_session() as session:

        try:
            names_d = _list_day(session, date, "d")
            names_o = None

            if names_d is not None and not _pick(station, names_d):
                names_o = _list_day(session, date, "o")

            candidates = _pick(station, names_d, names_o)

        except Exception as e:
            print(f"WARNING listing {date}: {e} -> trying the file names directly")
            candidates = _blind_candidates(station, date)

    if not candidates:
        print(f"NONE  {station} on {date} (DOY {doy:03d})")
        return None

    path = _download_one(date, candidates, output_dir)

    if path is None:
        print(f"NONE  {station} on {date} (DOY {doy:03d})")

    return path
