from pathlib import Path
from datetime import timedelta

import re
import json
import time
import tempfile
import threading
import subprocess
import requests
from concurrent.futures import ThreadPoolExecutor

from ionotec.download.extract import extract_and_remove_gz

# NASA CDDIS archive of daily IGS RINEX files.
# Authentication: Earthdata credentials in ~/.netrc
#   machine urs.earthdata.nasa.gov login <user> password <pass>
# (requests picks them up automatically, also after the redirect to urs.)
BASE_URL = "https://cddis.nasa.gov/archive/gnss/data/daily"

# True  -> use curl (needs curl and a working ~/.netrc); requests is only
#          used if curl is not installed. Default: plain requests gets
#          "Connection reset by peer" on the CDDIS daily listings.
# False -> try requests first and fall back to curl if it fails.
USE_CURL = True
_COOKIES = Path(tempfile.gettempdir()) / "igs_cddis_cookies.txt"

# Successful daily listings are saved here so each day is fetched only once
# (None -> no disk cache).
CACHE_DIR = Path.home() / ".cache" / "ionotec" / "igs"

# RINEX 3 long name, daily observation, Hatanaka or plain:
#   ABMF00GLP_R_20240960000_01D_30S_MO.crx.gz
_RE_R3 = re.compile(
    r"^([A-Z0-9]{9})_[A-Z]_\d{11}_\d{2}[A-Z]_\d{2}[A-Z]_MO\.(?:crx|rnx)\.gz$",
    re.IGNORECASE,
)
# RINEX 2 Hatanaka:  abmf0960.24d.gz
_RE_R2 = re.compile(
    r"^([A-Z0-9]{4})\d{3}[0-9A-X]\.\d{2}d\.gz$",
    re.IGNORECASE,
)

# Same patterns, to find file names anywhere in a page
_RE_FIND_R3 = re.compile(
    r"(?<![A-Za-z0-9])[A-Za-z0-9]{9}_[A-Za-z]_\d{11}_\d{2}[A-Za-z]_\d{2}[A-Za-z]"
    r"_MO\.(?:crx|rnx)\.gz",
    re.IGNORECASE,
)
_RE_FIND_R2 = re.compile(
    r"(?<![A-Za-z0-9])[A-Za-z0-9]{4}\d{3}[0-9A-Xa-x]\.\d{2}d\.gz",
    re.IGNORECASE,
)

# One listing per day is enough for the whole process
_LISTING_CACHE = {}


def _day_url(d):
    """Directory holding the daily Hatanaka files of date d."""
    yy = d.year % 100
    return f"{BASE_URL}/{d.year}/{d.timetuple().tm_yday:03d}/{yy:02d}d/"


def _cookie_file():
    """Cookie jar of the current thread (parallel curls can't share one)."""
    return _COOKIES.with_name(
        f"{_COOKIES.stem}_{threading.get_ident()}{_COOKIES.suffix}"
    )


def _new_session():
    # Plain session with default headers (this is what works with CDDIS)
    return requests.Session()


def _fetch_requests(session, url, dest=None, max_retry_time=60):
    """
    GET url exactly like get_dcb_from_cddis does: default headers,
    stream=True, timeout=30, retry every second for up to max_retry_time s,
    on any exception.

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
    Earthdata redirects, cookie jar kept between calls).
    Network failures (e.g. 'Connection reset by peer') are retried here,
    with a growing pause, because curl's own --retry ignores them.
    Returns bytes (dest=None) / True, or None on 404.
    """
    if dest is None:
        tmp = tempfile.NamedTemporaryFile(delete=False)
        tmp.close()
        out = Path(tmp.name)
    else:
        out = Path(dest)

    cmd = [
        "curl", "-sS", "-n", "-L",
        "-b", "",
        "-o", str(out), "-w", "%{http_code}", url,
    ]

    # Files: after a connection reset, continue from the bytes already
    # received instead of restarting from 0 (-C -). Not for listings.
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


def _list_day(session, d):
    """
    Observation file names (.gz) available on CDDIS for date d (cached).
    Returns None if the directory does not exist.
    """
    key = d.date() if hasattr(d, "date") else d

    if key in _LISTING_CACHE:
        return _LISTING_CACHE[key]

    cache_file = None
    if CACHE_DIR is not None:
        cache_file = Path(CACHE_DIR) / f"{key:%Y%j}.json"
        if cache_file.exists():
            try:
                _LISTING_CACHE[key] = json.loads(cache_file.read_text())
                return _LISTING_CACHE[key]
            except Exception:
                cache_file.unlink(missing_ok=True)

    body = _fetch(session, _day_url(d))

    if body is None:
        return None

    text = body.decode("utf-8", "replace")

    # Scan the whole page for file names, whatever the HTML looks like
    names = set(_RE_FIND_R3.findall(text)) | set(_RE_FIND_R2.findall(text))

    if not names:
        debug = Path(tempfile.gettempdir()) / f"igs_listing_{key:%Y%j}.html"
        debug.write_text(text)
        print(
            f"WARNING no observation file found in {_day_url(d)} "
            f"({len(body)} bytes). Page saved to {debug}"
        )
        print("        first 300 characters:", repr(text[:300]))

    _LISTING_CACHE[key] = sorted(names)

    if cache_file is not None and names:
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        cache_file.write_text(json.dumps(_LISTING_CACHE[key]))

    return _LISTING_CACHE[key]


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


DEFAULT_WORKERS = 4   # parallel connections to CDDIS


def _days(start, end):
    days = []
    d = start
    while d <= end:
        days.append(d)
        d += timedelta(days=1)
    return days


def _list_days(days, workers=DEFAULT_WORKERS):
    """{day: list of file names, or None if missing/failed}, in parallel."""

    def one(d):
        try:
            with _new_session() as session:
                return _list_day(session, d)
        except Exception as e:
            print(f"ERROR {d}: {e}")
            return None

    with ThreadPoolExecutor(max_workers=max(1, workers)) as ex:
        return dict(zip(days, ex.map(one, days)))


def stations_available_igs(start, end, full_id=False,
                           workers=DEFAULT_WORKERS):
    """
    Return the union of stations having daily RINEX observation files
    on the NASA CDDIS archive during the requested period.

    Archive:
        https://cddis.nasa.gov/archive/gnss/data/daily/YYYY/DDD/YYd/

    File:
        SSSS00CCC_R_YYYYDDD0000_01D_30S_MO.crx.gz

    Parameters
    ----------
    full_id : bool
        False -> 4-character lower-case names (like stations_available_cors).
        True  -> 9-character upper-case ids (e.g. 'ABMF00GLP').
    workers : int
        Daily listings fetched in parallel.
    """

    list_stations = set()

    listings = _list_days(_days(start, end), workers)

    for d, names in listings.items():

        if names is None:
            print(f"NONE  no listing for {d}")
            continue

        for fname in names:

            full, site = _station_ids(fname)

            if full_id:
                list_stations.add((full or site).upper())
            else:
                list_stations.add(site)

    return sorted(list_stations)


def _plan(station, days, listings):
    """
    [(day, file name, fallback file name or None)] to download for one
    station: ONE file per day.

    Long names (SSSS00CCC_R_..._MO.crx.gz, RINEX 3/4 naming) are preferred;
    the RINEX 2 short name (ssssDDD0.YYd.gz) is only used when the day has
    no long-name file.
    """
    tasks = []

    for d in days:

        doy = d.timetuple().tm_yday

        if len(station) == 9:
            # No listing needed: build the long name; if it does not
            # exist, the RINEX 2 file is the fallback.
            long_name = (
                f"{station.upper()}_R_{d.year}{doy:03d}0000_01D_30S_MO.crx.gz"
            )
            short_name = f"{station[:4]}{doy:03d}0.{d.year % 100:02d}d.gz"
            tasks.append((d, long_name, short_name))
            continue

        long_names, short_names = [], []

        for fname in listings.get(d) or []:
            full, site = _station_ids(fname)

            if site != station:
                continue

            (long_names if full else short_names).append(fname)

        if long_names:
            tasks.append((d, sorted(long_names)[0], None))
        elif short_names:
            tasks.append((d, sorted(short_names)[0], None))
        else:
            print(f"NONE  {station} on {d} (DOY {doy:03d})")

    return tasks


def _download_one(d, fname, output_dir, fallback=None):
    """Download + extract one file. Returns the Path, or None.
    If the file does not exist on the server (404) and `fallback` is given,
    that file is downloaded instead."""
    gz_path = output_dir / fname
    rinex_path = gz_path.with_suffix("")

    # Already extracted
    if rinex_path.exists():
        print(f"SKIP  {rinex_path.name}")
        return rinex_path

    try:
        print(f"GET   {fname}")

        with _new_session() as session:
            ok = _fetch(session, _day_url(d) + fname, dest=gz_path)

        if ok is None:
            print(f"NONE  {fname}")

            if fallback is not None:
                return _download_one(d, fallback, output_dir)

            return None

        return extract_and_remove_gz(gz_path)

    except Exception as e:

        # keep a partial file after a network failure: the next run
        # resumes it; anything else is removed
        if gz_path.exists() and not isinstance(e, ConnectionError):
            gz_path.unlink()

        print(f"ERROR {fname}: {e}")
        return None


def download_igs_stations(stations, start, end, output_dir=".",
                          workers=DEFAULT_WORKERS):
    """
    Download several stations at once, all files in parallel.

    stations: iterable of 4- or 9-character names.
    Returns {station: [Path, ...]}.
    """

    stations = [s.lower() for s in stations]

    for s in stations:
        if len(s) not in (4, 9):
            raise ValueError(f"bad station name: {s!r}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    days = _days(start, end)

    listings = {}
    if any(len(s) == 4 for s in stations):
        listings = _list_days(days, workers)

    tasks = []                      # (station, day, fname, fallback)
    for st in stations:
        tasks += [(st, d, f, fb) for d, f, fb in _plan(st, days, listings)]

    with ThreadPoolExecutor(max_workers=max(1, workers)) as ex:
        results = list(ex.map(
            lambda t: _download_one(t[1], t[2], output_dir, t[3]), tasks
        ))

    out = {st: [] for st in stations}
    for (st, _, _, _), path in zip(tasks, results):
        if path is not None:
            out[st].append(path)

    print()
    print(f"Files available locally: {sum(len(v) for v in out.values())}")

    return out


def download_igs(station, start, end, output_dir=".",
                 workers=DEFAULT_WORKERS):
    """
    Download IGS daily RINEX observation files from NASA CDDIS.

    Archive structure:
        https://cddis.nasa.gov/archive/gnss/data/daily/YYYY/DDD/YYd/

    station: 4-character site name or 9-character id (case-insensitive).
    With a 9-character id the file URL is built directly (no directory
    listing); with a 4-character name the daily listing is read (cached).
    Days are downloaded in parallel (`workers` connections).

    One file per station and day: the long-name file (RINEX 3/4 naming,
    SSSS00CCC_R_..._MO.crx.gz) when it exists, the RINEX 2 file
    (ssssDDD0.YYd.gz) only if it does not.

    The downloaded .gz files are extracted using extract_and_remove_gz().

    Returns
    -------
    list of Path
        Locally available RINEX files.
    """

    station = station.lower()

    return download_igs_stations(
        [station], start, end, output_dir, workers
    )[station]
