from pathlib import Path
from datetime import timedelta

import re
import time
import requests

from ionotec.download.extract import extract_and_remove_gz

BASE_URL = "http://garner.ucsd.edu/pub/rinex"

# RINEX 3 daily Hatanaka, any sampling interval (15S, 30S, 01H, ...):
#   ABMF00GLP_R_20240010000_01D_15S_MO.crx.gz
_RE_R3 = re.compile(
    r"^([A-Z0-9]{9})_[A-Z]_\d{11}_\d{2}[A-Z]_\d{2}[A-Z]_MO\.crx\.gz$",
    re.IGNORECASE,
)
# Same pattern, to find file names anywhere in the directory page
_RE_FIND = re.compile(
    r"(?<![A-Za-z0-9])[A-Za-z0-9]{9}_[A-Za-z]_\d{11}_\d{2}[A-Za-z]_\d{2}[A-Za-z]"
    r"_MO\.crx\.gz",
    re.IGNORECASE,
)

# One directory listing per day is enough for the whole process
_LISTING_CACHE = {}


def _day_url(d):
    """Directory holding the daily files of date d."""
    return f"{BASE_URL}/{d.year}/{d.timetuple().tm_yday:03d}/"


def _get(session, url, max_retries=4, **kwargs):
    """
    GET with retries. Returns the Response, or None on 404.
    The caller must close it (use 'with').
    """
    for attempt in range(max_retries):
        try:
            r = session.get(url, timeout=(10, 60), **kwargs)

            if r.status_code == 404:
                r.close()
                return None

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


def _list_day(session, d):
    """
    File names available on Garner for date d (cached).
    Returns None if the directory does not exist.
    """
    key = d.date() if hasattr(d, "date") else d

    if key in _LISTING_CACHE:
        return _LISTING_CACHE[key]

    r = _get(session, _day_url(d))

    if r is None:
        return None

    with r:
        names = sorted(set(_RE_FIND.findall(r.text)))

    _LISTING_CACHE[key] = names

    return names


def _station_ids(fname):
    """(full_id, site4) in lower case, or (None, None)."""
    m = _RE_R3.match(fname)
    if m:
        full = m.group(1).lower()
        return full, full[:4]

    return None, None


def stations_available_garner(start, end, full_id=False):
    """
    Return the union of stations having daily RINEX observation files
    on the Garner (SOPAC/UCSD) archive during the requested period.

    Archive:
        http://garner.ucsd.edu/pub/rinex/YYYY/DDD/

    File:
        SSSS00CCC_R_YYYYDDD0000_01D_15S_MO.crx.gz

    Parameters
    ----------
    full_id : bool
        False -> 4-character lower-case names (like stations_available_cors).
        True  -> 9-character upper-case ids (e.g. 'ABMF00GLP').
    """

    list_stations = set()

    d = start

    with requests.Session() as session:

        while d <= end:

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


def download_garner(station, start, end, output_dir="."):
    """
    Download Garner daily RINEX 3 Hatanaka files.

    Archive structure:
        http://garner.ucsd.edu/pub/rinex/YYYY/DDD/

    station: 4-character site name or 9-character id (case-insensitive).
    The sampling interval is whatever the archive has (15S, 30S, ...);
    if several files match the same station and day, the first one in the
    listing is used.

    The downloaded .gz files are extracted using extract_and_remove_gz().

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

    with requests.Session() as session:

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
            rinex_path = gz_path.with_suffix("")

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
                        for chunk in r.iter_content(chunk_size=1024 * 1024):
                            if chunk:
                                f.write(chunk)

                rinex_path = extract_and_remove_gz(gz_path)
                list_files.append(rinex_path)

            except Exception as e:

                if gz_path.exists():
                    gz_path.unlink()

                print(f"ERROR {fname}: {e}")

            d += timedelta(days=1)

    print()
    print(f"Files available locally: {len(list_files)}")

    return list_files



def download_file(station, date, output_dir="."):
    """
    Download ONE daily RINEX 3 Hatanaka observation file of a station
    from the Garner (SOPAC/UCSD) archive.

    Parameters
    ----------
    station : str
        4-character site name or 9-character id (case-insensitive).
    date : datetime or date
        The target date for the observation file.
    output_dir : str or Path, optional
        Directory where the file will be saved. Default is current directory.

    Returns
    -------
    Path or None
        Path of the final extracted file, or None if no file was found or on error.
    """
    station = station.lower()

    if len(station) not in (4, 9):
        raise ValueError("station must be a 4 or 9 character name")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    doy = date.timetuple().tm_yday

    with requests.Session() as session:
        try:
            names = _list_day(session, date)
        except Exception as e:
            print(f"WARNING listing {date}: {e}")
            names = None

        if names is None:
            print(f"NONE  no directory for {date}")
            return None

        fname = None
        for name in names:
            full, site = _station_ids(name)
            if (len(station) == 9 and full == station) or \
               (len(station) == 4 and site == station):
                fname = name
                break

    if fname is None:
        print(f"NONE  {station} on {date} (DOY {doy:03d})")
        return None

    gz_path = output_dir / fname
    rinex_path = gz_path.with_suffix("")

    # Already extracted locally
    if rinex_path.exists():
        print(f"SKIP  {rinex_path.name}")
        return rinex_path

    try:
        print(f"GET   {fname}")
        
        with requests.Session() as session:
            r = _get(session, _day_url(date) + fname, stream=True)

            if r is None:
                print(f"NONE  {fname}")
                return None

            with r:
                with open(gz_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            f.write(chunk)

        # Extract .gz -> .crx
        rinex_path = extract_and_remove_gz(gz_path)
        return rinex_path

    except Exception as e:
        if gz_path.exists():
            gz_path.unlink()
        rinex_path.unlink(missing_ok=True)

        print(f"ERROR {fname}: {e}")
        return None
