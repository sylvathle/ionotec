from pathlib import Path
from datetime import timedelta

import re
import time
import subprocess
import requests

from ionotec.download.extract import extract_and_remove_gz

BASE_URL = "https://data.earthscope.org/archive/gnss/data/daily/rinex4/obs"

# RINEX 4 daily observation, any sampling interval (15S, 30S, ...):
#   AB0900USA_R_20240600000_01D_15S_MO.rnx.gz
_RE_OBS = re.compile(
    r"^([A-Z0-9]{9})_[A-Z]_\d{11}_\d{2}[A-Z]_\d{2}[A-Z]_MO\.(?:rnx|crx)\.gz$",
    re.IGNORECASE,
)
# Same pattern, to find file names anywhere in the ?list response
_RE_FIND = re.compile(
    r"(?<![A-Za-z0-9])[A-Za-z0-9]{9}_[A-Za-z]_\d{11}_\d{2}[A-Za-z]_\d{2}[A-Za-z]"
    r"_MO\.(?:rnx|crx)\.gz(?![A-Za-z0-9_.])",
    re.IGNORECASE,
)

# One listing per day is enough for the whole process
_LISTING_CACHE = {}


def _get_token():
    """
    Access token of the EarthScope CLI (`es login` must have been run once).
    """
    try:
        return subprocess.check_output(
            ["es", "user", "get-access-token"],
            text=True,
        ).strip()

    except FileNotFoundError:
        raise RuntimeError(
            "The EarthScope CLI 'es' was not found in the PATH"
        )

    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            "Could not get an EarthScope access token "
            f"(is 'es login' done?): {e}"
        )


def _new_session():
    s = requests.Session()
    s.headers["Authorization"] = f"Bearer {_get_token()}"
    return s


def _get(session, url, max_retries=4, **kwargs):
    """
    GET with retries and token refresh (the token expires after a while).
    Returns the Response, or None on 404. The caller must close it
    (use 'with').
    """
    refreshed = False

    for attempt in range(max_retries):
        try:
            r = session.get(url, timeout=(10, 120), **kwargs)

            if r.status_code == 404:
                r.close()
                return None

            if r.status_code == 401 and not refreshed:
                r.close()
                session.headers["Authorization"] = f"Bearer {_get_token()}"
                refreshed = True
                continue

            r.raise_for_status()
            return r

        except (
            requests.exceptions.ConnectionError,
            requests.exceptions.ChunkedEncodingError,
            requests.exceptions.Timeout,
        ) as e:
            if attempt == max_retries - 1:
                raise

            wait = 2 * (attempt + 1)
            print(f"WARNING {url}: {e} -> retry in {wait} s")
            time.sleep(wait)


def _rnx2crx(rnx_path):
    """
    Hatanaka-compress a RINEX file: `rnx2crx file.rnx` -> file.crx.
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


def _day_url(d):
    """Directory holding the daily files of date d."""
    return f"{BASE_URL}/{d.year}/{d.timetuple().tm_yday:03d}/"


def _list_day(session, d):
    """
    File names available on EarthScope for date d (cached).
    Returns None if the directory does not exist.
    """
    key = d.date() if hasattr(d, "date") else d

    if key in _LISTING_CACHE:
        return _LISTING_CACHE[key]

    r = _get(session, _day_url(d) + "?list")

    if r is None:
        return None

    with r:
        names = sorted(set(_RE_FIND.findall(r.text)))

    _LISTING_CACHE[key] = names

    return names


def _station_ids(fname):
    """(full_id, site4) in lower case, or (None, None)."""
    m = _RE_OBS.match(fname)
    if m:
        full = m.group(1).lower()
        return full, full[:4]

    return None, None


def stations_available_earthscope(start, end, full_id=False):
    """
    Return the union of stations having RINEX 4 daily observation files
    on EarthScope during the requested period.

    Archive:
        https://data.earthscope.org/archive/gnss/data/daily/rinex4/obs/YYYY/DDD/

    File:
        SSSS00CCC_R_YYYYDDD0000_01D_15S_MO.rnx.gz

    Example:
        AB0900USA_R_20240600000_01D_15S_MO.rnx.gz

    Parameters
    ----------
    full_id : bool
        False -> 4-character lower-case names (like the other sources).
        True  -> 9-character upper-case ids (e.g. 'AB0900USA').

    Needs the EarthScope CLI ('es') logged in.
    """

    list_stations = set()

    d = start

    with _new_session() as session:

        while d < end:

            try:
                names = _list_day(session, d)

            except Exception as e:
                print(f"ERROR {d}: {e}")
                d += timedelta(days=1)
                continue

            if names is None:
                print(f"NONE  no directory for {d}")
                d += timedelta(days=1)
                continue

            for fname in names:

                full, site = _station_ids(fname)

                if full_id:
                    list_stations.add(full.upper())
                else:
                    list_stations.add(site)

            d += timedelta(days=1)

    return sorted(list_stations)


def download_earthscope(station, start, end, output_dir="."):
    """
    Download EarthScope daily RINEX 4 observation files.

    Archive structure:
        https://data.earthscope.org/archive/gnss/data/daily/rinex4/obs/YYYY/DDD/

    Filename:
        SSSS00CCC_R_YYYYDDD0000_01D_15S_MO.rnx.gz

    Example:
        https://data.earthscope.org/archive/gnss/data/daily/rinex4/obs/2024/060/AB0900USA_R_20240600000_01D_15S_MO.rnx.gz

    station: 4-character site name or 9-character id (case-insensitive).
    The sampling interval is whatever the archive has (15S, 30S, ...), so
    the daily listing is read (once per day, cached). If several files
    match the same station and day, the first one in the listing is used.

    The downloaded .gz files are extracted using extract_and_remove_gz(),
    and the resulting plain .rnx files are converted to Hatanaka (.crx)
    with `rnx2crx` (which must be in the PATH); the .rnx is then removed.
    Needs the EarthScope CLI ('es') logged in.

    Returns
    -------
    list of Path
        Locally available RINEX files.
    """

    station = station.lower()

    if len(station) not in (4, 9):
        raise ValueError("station must be a 4 or 9 character name")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    list_files = []

    d = start

    with _new_session() as session:

        while d <= end:

            doy = d.timetuple().tm_yday

            try:
                names = _list_day(session, d)

            except Exception as e:
                print(f"ERROR listing {d}: {e}")
                d += timedelta(days=1)
                continue

            fname = None

            for name in names or []:

                full, site = _station_ids(name)

                if (len(station) == 9 and full == station) or \
                   (len(station) == 4 and site == station):
                    fname = name
                    break

            if fname is None:
                print(f"NONE  {station} on {d} (DOY {doy:03d})")
                d += timedelta(days=1)
                continue

            gz_path = output_dir / fname
            # X.rnx.gz -> X.rnx -> (rnx2crx) -> X.crx
            rinex_path = gz_path.with_suffix("").with_suffix(".crx")

            # Already extracted
            if rinex_path.exists():
                print(f"SKIP  {rinex_path.name}")
                list_files.append(rinex_path)
                d += timedelta(days=1)
                continue

            try:
                print(f"GET   {fname}")

                r = _get(session, _day_url(d) + fname, stream=True)

                if r is None:
                    print(f"NONE  {fname}")
                    d += timedelta(days=1)
                    continue

                with r:
                    with open(gz_path, "wb") as f:
                        for chunk in r.iter_content(
                            chunk_size=1024 * 1024
                        ):
                            if chunk:
                                f.write(chunk)

                extracted = extract_and_remove_gz(gz_path)

                # Plain RINEX -> Hatanaka (.crx); .crx.gz files are
                # already Hatanaka and need no conversion.
                if extracted.suffix.lower() == ".rnx":
                    extracted = _rnx2crx(extracted)

                list_files.append(extracted)

            except Exception as e:

                # leave nothing half-done (.gz, extracted .rnx, .crx)
                for p in (gz_path, gz_path.with_suffix(""), rinex_path):
                    p.unlink(missing_ok=True)

                print(f"ERROR {fname}: {e}")

            d += timedelta(days=1)

    print()
    print(f"Files available locally: {len(list_files)}")

    return list_files


# RINEX 2 daily observation (plain), used when RINEX 4 is not available:
#   https://data.earthscope.org/archive/gnss/data/daily/rinex2/obs/YYYY/DDD/CHNR0970.24o.gz
BASE_URL_RINEX2 = "https://data.earthscope.org/archive/gnss/data/daily/rinex2/obs"


def download_file(station, date, output_dir=".", rinex2_fallback=True):
    """
    Download ONE daily RINEX observation file of a station from EarthScope.

    station : 4-character site name (case-insensitive)
    date    : datetime / date of the day wanted
    rinex2_fallback : if False, only RINEX 4 is looked for

    Priority (the first one that exists is downloaded, and only that one):
        1. RINEX 4    SSSS00CCC_R_YYYYDDD0000_01D_15S_MO.rnx.gz
                      (rinex4/obs/YYYY/DDD/, any sampling interval)
        2. RINEX 2    SSSSDDD0.YYo.gz  (rinex2/obs/YYYY/DDD/)

    RINEX 2 is only tried when RINEX 4 does not exist (not in the daily
    listing / 404). A network error stops the search instead of silently
    falling back to an older RINEX version.

    The .gz is extracted with extract_and_remove_gz(). A plain RINEX 4
    .rnx is converted to Hatanaka (.crx) with `rnx2crx` (must be in the
    PATH) and the .rnx removed; a RINEX 2 file stays plain (.YYo).
    Needs the EarthScope CLI ('es') logged in.

    Returns
    -------
    Path of the final file (e.g. .../AB0900USA_R_20240600000_01D_15S_MO.crx
    for RINEX 4, .../CHNR0970.24o for RINEX 2), or None if no file was
    found / on error.
    """

    station = station.lower()

    if len(station) != 4:
        raise ValueError("station must be a 4-character station name")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    doy = date.timetuple().tm_yday
    yy = date.year % 100

    with _new_session() as session:

        # --- candidates in priority order: (url, file name)
        candidates = []

        try:
            names = _list_day(session, date) or []

            for name in names:
                if _station_ids(name)[1] == station:
                    candidates.append((_day_url(date) + name, name))
                    break                        # first one in the listing

        except Exception as e:
            print(f"WARNING listing {date}: {e} -> RINEX 4 name unknown")

        if rinex2_fallback:
            dir2 = f"{BASE_URL_RINEX2}/{date.year}/{doy:03d}/"

            # the archive names are upper case (CHNR0970.24o.gz); lower
            # case is tried too in case a file is stored that way
            for st in (station.upper(), station):
                name = f"{st}{doy:03d}0.{yy:02d}o.gz"
                candidates.append((dir2 + name, name))

        # --- take the first one that exists
        for url, fname in candidates:

            gz_path = output_dir / fname
            extracted_path = gz_path.with_suffix("")

            # X.rnx.gz -> X.rnx -> (rnx2crx) -> X.crx ; X.YYo.gz -> X.YYo
            final_path = (
                extracted_path.with_suffix(".crx")
                if extracted_path.suffix.lower() == ".rnx"
                else extracted_path
            )

            if final_path.exists():
                print(f"SKIP  {final_path.name}")
                return final_path

            try:
                print(f"GET   {fname}")

                r = _get(session, url, stream=True)

                if r is None:                     # 404: try the next one
                    print(f"NONE  {fname}")
                    continue

                with r:
                    with open(gz_path, "wb") as f:
                        for chunk in r.iter_content(chunk_size=1024 * 1024):
                            if chunk:
                                f.write(chunk)

                extracted = extract_and_remove_gz(gz_path)

                # Plain RINEX 4 -> Hatanaka; .crx.gz and RINEX 2 are kept
                if extracted.suffix.lower() == ".rnx":
                    extracted = _rnx2crx(extracted)

                return extracted

            except Exception as e:
                # leave nothing half-done (.gz, extracted .rnx, .crx)
                for p in (gz_path, extracted_path, final_path):
                    p.unlink(missing_ok=True)

                print(f"ERROR {fname}: {e}")
                return None

    print(f"NONE  {station} on {date} (DOY {doy:03d})")
    return None
