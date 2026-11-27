"""
Download daily RINEX observation files from the Canadian Geodetic Survey
(CGS, Natural Resources Canada) public archive.

Archive (anonymous HTTPS):
    https://cacs.nrcan.gc.ca/gps/data/gpsdata/{yy}{doy}/{yy}d/

The daily directory holds two kinds of observation files (30 s):
    RINEX 3/4 long name:  ALGO00CAN_R_20241250000_01D_30S_MO.crx.gz
    RINEX 2 Hatanaka:     algo1250.24d.Z          (Unix compress, .Z)

From 1 March 2027 CGS only publishes RINEX >= 3.
Same public API as sonel.py; RINEX 3/4 is always preferred over RINEX 2.
"""

import re
import shutil
import subprocess
import time
from datetime import timedelta
from pathlib import Path

import requests

from ionotec.download.extract import extract_and_remove_gz

BASE_URL = "https://cacs.nrcan.gc.ca/gps/data/gpsdata"

# RINEX 3/4 long name:  ALGO00CAN_R_20241250000_01D_30S_MO.crx.gz
_RE_R3 = re.compile(r"^([A-Z0-9]{9})_[A-Z]_\d{11}_", re.IGNORECASE)
# RINEX 2 Hatanaka:     algo1250.24d.Z  (also accept .gz)
_RE_R2 = re.compile(r"^([A-Z0-9]{4})\d{3}[0-9A-X]\.\d{2}d\.(Z|gz)$", re.IGNORECASE)
# Navigation files (R3: _MN. / _GN. ..., R2: .24n, .24g ...)
_RE_NAV = re.compile(r"_[A-Z]N\.|\.\d{2}[nglhpq]\.", re.IGNORECASE)
# File links in the Apache directory listing
_RE_HREF = re.compile(r'href="([^"?/][^"/]*)"', re.IGNORECASE)


# ----------------------------------------------------------------------
# HTTPS helpers
# ----------------------------------------------------------------------

class _NrcanHTTP:
    """HTTPS session to the CGS archive, with retries."""

    def __init__(self, timeout=60):
        self.timeout = timeout
        self.session = requests.Session()

    def close(self):
        self.session.close()

    @staticmethod
    def day_url(year, doy):
        yy = year % 100
        return f"{BASE_URL}/{yy:02d}{doy:03d}/{yy:02d}d"

    def list_day(self, year, doy, retries=3):
        """
        Names of the files in gpsdata/YYDDD/YYd/.
        Returns None if the directory does not exist, [] if it is empty.
        """
        url = self.day_url(year, doy) + "/"

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
    """True for compressed daily observation files (R3/R4 .gz or R2 Hatanaka .Z/.gz)."""
    if not fname.lower().endswith((".gz", ".z")):
        return False
    if _RE_NAV.search(fname):
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


def _extract(path):
    """
    Decompress a downloaded file and remove the archive.
    .gz -> extract_and_remove_gz ; .Z (Unix compress, RINEX 2) -> gzip -d.
    Returns the Path of the extracted file.
    """
    path = Path(path)

    if path.suffix.lower() == ".gz":
        return Path(extract_and_remove_gz(path))

    if path.suffix == ".Z":
        out = path.with_suffix("")
        gzip_bin = shutil.which("gzip")
        if gzip_bin is None:
            raise RuntimeError("gzip is needed to decompress .Z files")
        with open(out, "wb") as f:
            subprocess.run([gzip_bin, "-dc", str(path)], stdout=f, check=True)
        path.unlink()
        return out

    raise ValueError(f"Unknown compression: {path.name}")


def _local_name(fname):
    """Name of the extracted file (archive suffix removed)."""
    return re.sub(r"\.(gz|Z)$", "", fname)


def _rinex_version(path):
    """RINEX version (float) read from the header of an extracted file
    (plain or Hatanaka), or None."""
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


def _best_file(conn, remote_dir, r3_candidates, r2_candidates, output_dir):
    """
    Get the highest RINEX version available for one station and day.
    RINEX 3/4 long names come first; when several long-name files exist
    their headers decide (highest version kept, the others deleted).
    Falls back to RINEX 2 only if no long-name file could be obtained.
    """
    for group in (sorted(r3_candidates), sorted(r2_candidates)):
        got = []
        for fname in group:
            local = output_dir / _local_name(fname)
            if local.exists():
                print(f"SKIP  {local.name}")
                path = local
            else:
                path = _download(conn, remote_dir, fname, output_dir)
            if path is not None:
                got.append(path)

        if got:
            if len(got) == 1:
                return got[0]
            best = max(got, key=lambda p: _rinex_version(p) or 0.0)
            for p in got:
                if p != best:
                    p.unlink(missing_ok=True)
            return best

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
        if not _is_obs_file(fname):
            continue
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

def stations_available_nrcan(start, end, full_id=False):
    """
    Return the union of stations having daily RINEX observation files
    on the NRCan archive during the requested period [start, end).

    Archive:
        https://cacs.nrcan.gc.ca/gps/data/gpsdata/YYDDD/YYd/

    Parameters
    ----------
    full_id : bool
        False -> 4-character names (lower case), like stations_available_cors.
        True  -> 9-character ids (e.g. 'ALGO00CAN'). Stations that only
                 have RINEX 2 files are returned with their 4-character name.

    Returns
    -------
    sorted list of str
    """
    list_stations = set()
    conn = _NrcanHTTP()

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

            r3_sites = set()
            r2_sites = set()
            for fname in names:
                if not _is_obs_file(fname):
                    continue
                full, site = _station_ids(fname)
                if full:
                    r3_sites.add(site)
                    list_stations.add(full.upper() if full_id else site)
                else:
                    r2_sites.add(site)

            for site in r2_sites:
                # with full_id, only add the 4-char name when no long name exists
                if not full_id or site not in r3_sites:
                    list_stations.add(site.upper() if full_id else site)

    finally:
        conn.close()

    return sorted(list_stations)


def list_files_nrcan(station, start, end):
    """
    List the observation files of one station available on the NRCan
    archive between start and end (inclusive).

    station: 4-character site name or 9-character id (case-insensitive).

    Returns
    -------
    list of (date, filename)   (RINEX 3/4 first for each day)
    """
    station = _check_station(station)
    found = []
    conn = _NrcanHTTP()

    try:
        for d, year, doy in _iter_days(start, end):
            try:
                names = conn.list_day(year, doy)
            except Exception as e:
                print(f"ERROR {d}: {e}")
                continue

            r3, r2 = _split_candidates(names, station)
            found += [(d, f) for f in sorted(r3) + sorted(r2)]

    finally:
        conn.close()

    return found


def download_nrcan(station, start, end, output_dir="."):
    """
    Download NRCan daily RINEX observation files, one per day,
    with the highest RINEX version available (RINEX 3/4 over RINEX 2).

    Archive structure:
        https://cacs.nrcan.gc.ca/gps/data/gpsdata/YYDDD/YYd/

    station: 4-character site name or 9-character id (case-insensitive).
    .gz files are extracted with extract_and_remove_gz(), .Z files with gzip.

    Returns
    -------
    list of Path
        Locally available RINEX files.
    """
    station = _check_station(station)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    list_files = []
    conn = _NrcanHTTP()

    try:
        for d, year, doy in _iter_days(start, end):
            remote_dir = conn.day_url(year, doy)

            try:
                names = conn.list_day(year, doy)
            except Exception as e:
                print(f"ERROR accessing {remote_dir}: {e}")
                continue

            r3, r2 = _split_candidates(names, station)
            if not (r3 or r2):
                print(f"NONE  {station} on {d} (DOY {doy:03d})")
                continue

            path = _best_file(conn, remote_dir, r3, r2, output_dir)
            if path is not None:
                list_files.append(path)

    finally:
        conn.close()

    print()
    print(f"Files available locally: {len(list_files)}")

    return list_files


def download_file(station, date, output_dir="."):
    """
    Download ONE daily RINEX observation file for a station from NRCan,
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

    conn = _NrcanHTTP()
    remote_dir = conn.day_url(year, doy)

    try:
        try:
            names = conn.list_day(year, doy)
        except Exception as e:
            print(f"ERROR accessing {remote_dir}: {e}")
            return None

        if names is None:
            print(f"NONE  no directory for {date} (DOY {doy:03d})")
            return None

        r3, r2 = _split_candidates(names, station)
        if not (r3 or r2):
            print(f"NONE  {station} on {date} (DOY {doy:03d})")
            return None

        return _best_file(conn, remote_dir, r3, r2, output_dir)

    finally:
        conn.close()
