from pathlib import Path
from datetime import timedelta
from urllib.parse import urljoin

import re
import time
import subprocess
import requests

from ionotec.download.extract import extract_and_remove_gz

BASE_URL = "https://data.geonet.org.nz/v1/data/gnss/rinex"

# RINEX 3 daily observation, any sampling interval (15S, 30S, ...):
#   ARTA00NZL_R_20241300000_01D_30S_MO.rnx.gz
_RE_OBS = re.compile(
    r"^([A-Z0-9]{9})_[A-Z]_\d{11}_\d{2}[A-Z]_\d{2}[A-Z]_MO\.(?:rnx|crx)\.gz$",
    re.IGNORECASE,
)
# Same pattern, to find file names anywhere in the directory page
_RE_FIND = re.compile(
    r"(?<![A-Za-z0-9])[A-Za-z0-9]{9}_[A-Za-z]_\d{11}_\d{2}[A-Za-z]_\d{2}[A-Za-z]"
    r"_MO\.(?:rnx|crx)\.gz(?![A-Za-z0-9_.])",
    re.IGNORECASE,
)

# One listing per day is enough for the whole process: {day: {name: url}}
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
            r = session.get(url, timeout=(10, 120), **kwargs)

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
    {file name: url} of the RINEX 3 observation files GeoNet has for
    date d (cached). Returns None if the directory does not exist.
    """
    key = d.date() if hasattr(d, "date") else d

    if key in _LISTING_CACHE:
        return _LISTING_CACHE[key]

    day_url = _day_url(d)

    r = _get(session, day_url)

    if r is None:
        return None

    with r:
        text = r.text

    files = {}

    # Links of the HTML listing (absolute or relative)
    for href in re.findall(r'href="([^"]+)"', text, re.IGNORECASE):
        fname = href.split("?")[0].rstrip("/").split("/")[-1]

        if _RE_OBS.match(fname):
            files[fname] = urljoin(day_url, href)

    # Names that only appear as plain text
    for fname in _RE_FIND.findall(text):
        files.setdefault(fname, day_url + fname)

    _LISTING_CACHE[key] = files

    return files


def _station_ids(fname):
    """(full_id, site4) in lower case, or (None, None)."""
    m = _RE_OBS.match(fname)
    if m:
        full = m.group(1).lower()
        return full, full[:4]

    return None, None


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


def stations_available_geonet(start, end, full_id=False):
    """
    Return the union of stations having RINEX 3 daily observation files
    on GeoNet (New Zealand) during the requested period.

    Archive:
        https://data.geonet.org.nz/v1/data/gnss/rinex/YYYY/DDD/

    File:
        SSSS00CCC_R_YYYYDDD0000_01D_30S_MO.rnx.gz

    Example:
        ARTA00NZL_R_20241300000_01D_30S_MO.rnx.gz

    Parameters
    ----------
    full_id : bool
        False -> 4-character lower-case names (like the other sources).
        True  -> 9-character upper-case ids (e.g. 'ARTA00NZL').
    """

    list_stations = set()

    d = start

    with requests.Session() as session:

        while d < end:

            try:
                files = _list_day(session, d)

            except Exception as e:
                print(f"ERROR {d}: {e}")
                d += timedelta(days=1)
                continue

            if files is None:
                print(f"NONE  no directory for {d}")
                d += timedelta(days=1)
                continue

            for fname in files:

                full, site = _station_ids(fname)

                if full_id:
                    list_stations.add(full.upper())
                else:
                    list_stations.add(site)

            d += timedelta(days=1)

    return sorted(list_stations)


def download_geonet(station, start, end, output_dir="."):
    """
    Download GeoNet RINEX 3 daily observation files for a station.

    Archive structure:
        https://data.geonet.org.nz/v1/data/gnss/rinex/YYYY/DDD/

    Filename:
        SSSS00CCC_R_YYYYDDD0000_01D_30S_MO.rnx.gz

    station: 4-character site name or 9-character id (case-insensitive).
    If several files match the same station and day (e.g. different
    sampling intervals), the first one in the listing is used.

    The downloaded .gz files are extracted using extract_and_remove_gz(),
    and the resulting plain .rnx files are converted to Hatanaka (.crx)
    with `rnx2crx` (which must be in the PATH); the .rnx is then removed.

    Returns
    -------
    list of Path
        Locally available Hatanaka RINEX (.crx) files.
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
                files = _list_day(session, d)

            except Exception as e:
                print(f"ERROR listing {d}: {e}")
                d += timedelta(days=1)
                continue

            fname = None

            for name in files or {}:

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

                r = _get(session, files[fname], stream=True)

                if r is None:
                    print(f"NONE  {fname}")
                    d += timedelta(days=1)
                    continue

                with r:
                    with open(gz_path, "wb") as f:
                        for chunk in r.iter_content(chunk_size=1024 * 1024):
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


def download_file(station, date, output_dir="."):
    """
    Download ONE daily GeoNet RINEX 3 observation file for a station.

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
        Path of the final extracted Hatanaka (.crx) file, or None if no file
        was found or on error.
    """
    station = station.lower()

    if len(station) not in (4, 9):
        raise ValueError("station must be a 4 or 9 character name")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    doy = date.timetuple().tm_yday

    with requests.Session() as session:
        try:
            files = _list_day(session, date)
        except Exception as e:
            print(f"WARNING listing {date}: {e}")
            files = None

        if files is None:
            print(f"NONE  no directory for {date}")
            return None

        fname = None
        for name in files:
            full, site = _station_ids(name)
            if (len(station) == 9 and full == station) or \
               (len(station) == 4 and site == station):
                fname = name
                break

    if fname is None:
        print(f"NONE  {station} on {date} (DOY {doy:03d})")
        return None

    gz_path = output_dir / fname
    # X.rnx.gz -> X.rnx -> (rnx2crx) -> X.crx
    rinex_path = gz_path.with_suffix("").with_suffix(".crx")

    # Already extracted locally
    if rinex_path.exists():
        print(f"SKIP  {rinex_path.name}")
        return rinex_path

    try:
        print(f"GET   {fname}")

        with requests.Session() as session:
            r = _get(session, files[fname], stream=True)

            if r is None:
                print(f"NONE  {fname}")
                return None

            with r:
                with open(gz_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            f.write(chunk)

        extracted = extract_and_remove_gz(gz_path)

        # Plain RINEX -> Hatanaka (.crx); .crx.gz files are
        # already Hatanaka and need no conversion.
        if extracted.suffix.lower() == ".rnx":
            extracted = _rnx2crx(extracted)

        return extracted

    except Exception as e:
        # Leave nothing half-done (.gz, extracted .rnx, .crx)
        for p in (gz_path, gz_path.with_suffix(""), rinex_path):
            p.unlink(missing_ok=True)

        print(f"ERROR {fname}: {e}")
        return None
