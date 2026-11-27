"""
Download daily RINEX observation files from PANGA (Pacific Northwest
Geodetic Array, Central Washington University).

Archive (anonymous HTTPS, no login):
    https://www.geodesy.org/pub/data/{yyyy}/{doy}/rinex3/   RINEX 3/4 long names
        ADEL00USA_R_20241250000_01D_30S_MO.crx.bz2     (01S, 05S or 30S sampling)
    https://www.geodesy.org/pub/data/{yyyy}/{doy}/30sec/    RINEX 2 Hatanaka, 30 s
        adel1250.24d.bz2

Priority for each station and day (highest RINEX version first):
    1. RINEX 3/4, 30 s
    2. RINEX 3/4 at a higher rate (5 s, then 1 s) - only if highrate=True;
       these daily files are much larger (1 s: ~25 MB compressed vs ~2 MB at 30 s)
    3. RINEX 2, 30 s
Files are bz2-compressed; .gz and .Z are also handled.
Same public API as sonel.py / nrcan.py.
"""

import bz2
import re
import shutil
import subprocess
import time
from datetime import timedelta
from pathlib import Path

import requests

from ionotec.download.extract import extract_and_remove_gz

BASE_URL = "https://www.geodesy.org/pub/data"
R3_DIR = "rinex3"
R2_DIR = "30sec"

# RINEX 3/4 long name:  ADEL00USA_R_20241250000_01D_30S_MO.crx.bz2
_RE_R3 = re.compile(r"^([A-Z0-9]{9})_[A-Z]_\d{11}_\d{2}[DHM]_(\d{2})([SMHD])_MO\.", re.IGNORECASE)
# RINEX 2 Hatanaka:     adel1250.24d.bz2
_RE_R2 = re.compile(r"^([A-Z0-9]{4})\d{3}[0-9A-X]\.\d{2}d\.(bz2|gz|Z)$", re.IGNORECASE)
_COMPRESSED = (".bz2", ".gz", ".z")
# File links in the Apache directory listing
_RE_HREF = re.compile(r'href="([^"?/][^"/]*)"', re.IGNORECASE)
_UNIT_S = {"S": 1, "M": 60, "H": 3600, "D": 86400}


# ----------------------------------------------------------------------
# HTTPS helpers
# ----------------------------------------------------------------------

class _PangaHTTP:
    """HTTPS session to the PANGA archive, with retries."""

    def __init__(self, timeout=120):
        self.timeout = timeout
        self.session = requests.Session()

    def close(self):
        self.session.close()

    @staticmethod
    def day_url(year, doy, sub):
        return f"{BASE_URL}/{year}/{doy:03d}/{sub}"

    def list_dir(self, url, retries=3):
        """Names of the files in one directory, None if it does not exist."""
        url = url + "/"
        for attempt in range(retries):
            try:
                r = self.session.get(url, timeout=self.timeout)
                if r.status_code == 404:
                    return None
                r.raise_for_status()
                return [n for n in _RE_HREF.findall(r.text) if not n.startswith("http")]
            except requests.RequestException as e:
                print(f"WARNING {url} (attempt {attempt + 1}): {e}")
                time.sleep(2 ** attempt)
        raise ConnectionError(f"Could not list {url}")

    def list_day(self, year, doy):
        """
        Observation files of one day as {fname: remote_dir}, from rinex3/ and 30sec/.
        Returns None if neither directory exists.
        """
        out, found = {}, False
        for sub in (R3_DIR, R2_DIR):
            remote_dir = self.day_url(year, doy, sub)
            names = self.list_dir(remote_dir)
            if names is None:
                continue
            found = True
            for n in names:
                if _is_obs_file(n):
                    out[n] = remote_dir
        return out if found else None

    def retrieve(self, remote_dir, fname, local_path, retries=3):
        url = f"{remote_dir}/{fname}"
        tmp_path = Path(f"{local_path}.part")
        for attempt in range(retries):
            try:
                with self.session.get(url, stream=True, timeout=self.timeout) as r:
                    r.raise_for_status()
                    with open(tmp_path, "wb") as f:
                        for chunk in r.iter_content(1 << 20):
                            f.write(chunk)
                tmp_path.replace(local_path)
                return
            except Exception as e:
                print(f"WARNING {fname} (attempt {attempt + 1}): {e}")
                tmp_path.unlink(missing_ok=True)
                time.sleep(2 ** attempt)
        raise ConnectionError(f"Could not download {fname}")


def _is_obs_file(fname):
    """True for compressed daily observation files (R3/R4 MO or R2 Hatanaka)."""
    if not fname.lower().endswith(_COMPRESSED):
        return False
    return bool(_RE_R3.match(fname) or _RE_R2.match(fname))


def _station_ids(fname):
    """
    Return (full_id, site4) in lower case from a file name.
    full_id is the 9-character id (RINEX 3/4) or None (RINEX 2).
    """
    m = _RE_R3.match(fname)
    if m:
        full = m.group(1).lower()
        return full, full[:4]
    m = _RE_R2.match(fname)
    if m:
        return None, m.group(1).lower()
    return None, None


def _sampling_s(fname):
    """Sampling interval in seconds of a RINEX 3/4 long name (30 for RINEX 2)."""
    m = _RE_R3.match(fname)
    if not m:
        return 30
    return int(m.group(2)) * _UNIT_S[m.group(3).upper()]


def _local_name(fname):
    """Name of the extracted file (compression suffix removed)."""
    return re.sub(r"\.(bz2|gz|Z)$", "", fname)


def _extract(path):
    """
    Decompress a downloaded file and remove the archive.
    .bz2 -> bz2 module ; .gz -> extract_and_remove_gz ; .Z -> gzip -d.
    """
    path = Path(path)
    out = path.with_suffix("")

    if path.suffix == ".bz2":
        with bz2.open(path, "rb") as fi, open(out, "wb") as fo:
            shutil.copyfileobj(fi, fo, 1 << 20)
        path.unlink()
        return out

    if path.suffix.lower() == ".gz":
        return Path(extract_and_remove_gz(path))

    if path.suffix == ".Z":
        gzip_bin = shutil.which("gzip")
        if gzip_bin is None:
            raise RuntimeError("gzip is needed to decompress .Z files")
        with open(out, "wb") as f:
            subprocess.run([gzip_bin, "-dc", str(path)], stdout=f, check=True)
        path.unlink()
        return out

    raise ValueError(f"Unknown compression: {path.name}")


def _rinex_version(path):
    """RINEX version (float) read from the header of an extracted file, or None."""
    try:
        with open(path, "r", errors="replace") as f:
            for _ in range(3):
                line = f.readline()
                if "RINEX VERSION / TYPE" in line:
                    return float(line[:9])
    except Exception:
        pass
    return None


def _download(conn, remote_dir, fname, output_dir):
    """Download + extract one file. Returns the extracted Path, or None."""
    archive = output_dir / fname
    try:
        print(f"GET   {fname}")
        conn.retrieve(remote_dir, fname, archive)
        return _extract(archive)
    except Exception as e:
        archive.unlink(missing_ok=True)
        (output_dir / _local_name(fname)).unlink(missing_ok=True)
        print(f"ERROR {fname}: {e}")
        return None


def _ordered_candidates(r3, r2, highrate):
    """
    Candidates in priority order: RINEX 3/4 30 s (or coarser), then RINEX 3/4
    high rate (5 s before 1 s) if highrate, then RINEX 2.
    """
    r3_std = sorted(f for f in r3 if _sampling_s(f) >= 30)
    r3_hr = sorted((f for f in r3 if _sampling_s(f) < 30), key=_sampling_s, reverse=True)
    return r3_std + (r3_hr if highrate else []) + sorted(r2)


def _best_file(conn, dirs, r3, r2, output_dir, highrate):
    """First candidate (in priority order) that is cached or downloads correctly."""
    for fname in _ordered_candidates(r3, r2, highrate):
        local = output_dir / _local_name(fname)
        if local.exists():
            print(f"SKIP  {local.name}")
            return local
        path = _download(conn, dirs[fname], fname, output_dir)
        if path is not None:
            return path
    return None


def _iter_days(start, end, inclusive=True):
    d = start
    while d <= end if inclusive else d < end:
        yield d, d.year, d.timetuple().tm_yday
        d += timedelta(days=1)


def _split_candidates(names, station):
    """(r3, r2) observation files of `station` (4 or 9 characters, lower case)."""
    r3, r2 = [], []
    for fname in names or []:
        full, site = _station_ids(fname)
        if len(station) == 9:
            if full == station:
                r3.append(fname)
        elif site == station:
            (r3 if full else r2).append(fname)
    return r3, r2


def _check_station(station):
    station = station.lower()
    if len(station) not in (4, 9):
        raise ValueError("station must be a 4 or 9 character name")
    return station


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------

def stations_available_panga(start, end, full_id=False):
    """
    Return the union of stations having daily RINEX observation files
    (rinex3/ or 30sec/) on the PANGA archive during [start, end).

    Parameters
    ----------
    full_id : bool
        False -> 4-character names (lower case), like stations_available_cors.
        True  -> 9-character ids (e.g. 'ADEL00USA'). Stations that only
                 have RINEX 2 files are returned with their 4-character name.

    Returns
    -------
    sorted list of str
    """
    list_stations = set()
    conn = _PangaHTTP()

    try:
        for d, year, doy in _iter_days(start, end, inclusive=False):
            try:
                names = conn.list_day(year, doy)
            except Exception as e:
                print(f"ERROR {d}: {e}")
                continue

            if names is None:
                print(f"NONE  no directory for {d} (DOY {doy:03d})")
                continue

            r3_sites, r2_sites = set(), set()
            for fname in names:
                full, site = _station_ids(fname)
                if full:
                    r3_sites.add(site)
                    list_stations.add(full.upper() if full_id else site)
                else:
                    r2_sites.add(site)

            for site in r2_sites:
                if not full_id or site not in r3_sites:
                    list_stations.add(site.upper() if full_id else site)

    finally:
        conn.close()

    return sorted(list_stations)


def list_files_panga(station, start, end, highrate=True):
    """
    List the observation files of one station available on the PANGA
    archive between start and end (inclusive), in download priority order.

    station: 4-character site name or 9-character id (case-insensitive).

    Returns
    -------
    list of (date, filename)
    """
    station = _check_station(station)
    found = []
    conn = _PangaHTTP()

    try:
        for d, year, doy in _iter_days(start, end):
            try:
                names = conn.list_day(year, doy)
            except Exception as e:
                print(f"ERROR {d}: {e}")
                continue

            r3, r2 = _split_candidates(names, station)
            found += [(d, f) for f in _ordered_candidates(r3, r2, highrate)]

    finally:
        conn.close()

    return found


def download_panga(station, start, end, output_dir=".", highrate=True):
    """
    Download PANGA daily RINEX observation files, one per day, with the
    highest RINEX version available (see module docstring for the priority).

    station: 4-character site name or 9-character id (case-insensitive).
    highrate: allow 1 s / 5 s RINEX 3 daily files before falling back to
              RINEX 2 30 s (large downloads).

    Returns
    -------
    list of Path
        Locally available RINEX files.
    """
    station = _check_station(station)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    list_files = []
    conn = _PangaHTTP()

    try:
        for d, year, doy in _iter_days(start, end):
            try:
                names = conn.list_day(year, doy)
            except Exception as e:
                print(f"ERROR accessing {BASE_URL}/{year}/{doy:03d}: {e}")
                continue

            r3, r2 = _split_candidates(names, station)
            if not _ordered_candidates(r3, r2, highrate):
                print(f"NONE  {station} on {d} (DOY {doy:03d})")
                continue

            path = _best_file(conn, names, r3, r2, output_dir, highrate)
            if path is not None:
                list_files.append(path)

    finally:
        conn.close()

    print()
    print(f"Files available locally: {len(list_files)}")

    return list_files


def download_file(station, date, output_dir=".", highrate=True):
    """
    Download ONE daily RINEX observation file for a station from PANGA,
    matching by the first 4 characters and prioritizing RINEX 3/4
    (9-character station ID) over RINEX 2 if available.

    Parameters
    ----------
    station : str
        4-character site name (case-insensitive).
    date : datetime or date
        The target date for the observation file.
    output_dir : str or Path, optional
        Directory where the file will be saved. Default is current directory.
    highrate : bool, optional
        Allow 1 s / 5 s RINEX 3 daily files before RINEX 2 30 s. Default True.

    Returns
    -------
    Path or None
        Path of the final extracted file, or None if no file was found or on error.
    """
    station = station.lower()
    if len(station) != 4:
        raise ValueError("station must be a 4-character station name")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    year = date.year
    doy = date.timetuple().tm_yday

    conn = _PangaHTTP()

    try:
        try:
            names = conn.list_day(year, doy)
        except Exception as e:
            print(f"ERROR accessing {BASE_URL}/{year}/{doy:03d}: {e}")
            return None

        if names is None:
            print(f"NONE  no directory for {date} (DOY {doy:03d})")
            return None

        r3, r2 = _split_candidates(names, station)
        if not _ordered_candidates(r3, r2, highrate):
            print(f"NONE  {station} on {date} (DOY {doy:03d})")
            return None

        return _best_file(conn, names, r3, r2, output_dir, highrate)

    finally:
        conn.close()
