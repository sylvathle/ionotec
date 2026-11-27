"""
Download daily RINEX observation files from the IGS Global Data Center
at IGN France (IGS network plus EUREF, RGP and other networks: ~1,900
stations per day).

Archive (anonymous FTP; HTTPS not reachable):
    ftp://igs.ign.fr/pub/igs/data/{yyyy}/{doy}/       (flat: obs + nav together)

    LPGS00ARG_R_20241250000_01D_30S_MO.crx.gz      RINEX 3/4, 30 s
    AB0700USA_R_20241250000_01D_15S_MO.crx.gz      RINEX 3/4, 15 s
    lpgs1250.24d.Z / .gz                           RINEX 2 Hatanaka

Priority for each station and day:
    1. RINEX 3/4 at 30 s
    2. RINEX 3/4 at 15 s (only when no 30 s file exists)
    3. RINEX 3/4 at any other rate (coarser first)
    4. RINEX 2 (Hatanaka; only one of .gz / .Z is downloaded)
Same public API as sonel.py.
"""

import re
import shutil
import subprocess
import time
from datetime import timedelta
from ftplib import FTP, error_perm, error_temp
from pathlib import Path

from ionotec.download.extract import extract_and_remove_gz

FTP_HOST = "igs.ign.fr"
FTP_ROOT = "/pub/igs/data"

# RINEX 3/4 daily observation file:  LPGS00ARG_R_20241250000_01D_30S_MO.crx.gz
_RE_R3 = re.compile(r"^([A-Z0-9]{9})_[A-Z]_\d{11}_01D_(\d{2})([SMHD])_MO\.", re.IGNORECASE)
_UNIT_S = {"S": 1, "M": 60, "H": 3600, "D": 86400}
# RINEX 2 Hatanaka:     lpgs1250.24d.gz, lpgs1250.24d.Z
_RE_R2 = re.compile(r"^([A-Z0-9]{4})\d{3}[0-9A-X]\.\d{2}d\.(Z|gz)$", re.IGNORECASE)
# Navigation files (R3: _MN. / _GN. ..., R2: .24n, .24g ...)
_RE_NAV = re.compile(r"_[A-Z]N\.|\.\d{2}[nglhpq]\.", re.IGNORECASE)


# ----------------------------------------------------------------------
# HTTPS helpers
# ----------------------------------------------------------------------

class _IgnFTP:
    """FTP connection to IGN that reconnects when the server drops it."""

    def __init__(self, timeout=60):
        self.timeout = timeout
        self.ftp = None
        self.connect()

    def connect(self, retries=6):
        """Connect + anonymous login, retried with backoff (421 / network errors)."""
        self.close()
        for attempt in range(retries):
            try:
                self.ftp = FTP(FTP_HOST, timeout=self.timeout)
                self.ftp.login()
                return
            except (error_temp, OSError) as e:
                self.close()
                wait = min(5 * 2 ** attempt, 120)
                print(f"WARNING connect (attempt {attempt + 1}): {e} - retry in {wait}s")
                time.sleep(wait)
        raise ConnectionError(f"Could not connect to {FTP_HOST}")

    def close(self):
        if self.ftp is not None:
            try:
                self.ftp.quit()
            except Exception:
                try:
                    self.ftp.close()
                except Exception:
                    pass
            self.ftp = None

    @staticmethod
    def day_dir(year, doy):
        return f"{FTP_ROOT}/{year}/{doy:03d}"

    def list_day(self, year, doy, retries=3):
        """
        Names of the files in YYYY/DDD.
        Returns None if the directory does not exist, [] if it is empty.
        """
        remote_dir = self.day_dir(year, doy)

        for attempt in range(retries):
            try:
                self.ftp.cwd(remote_dir)
                names = self.ftp.nlst()
                return [n.split("/")[-1] for n in names]

            except error_perm as e:
                # 550: no such directory / no files
                if str(e).startswith("550"):
                    return None
                raise

            except Exception as e:
                print(f"WARNING {remote_dir} (attempt {attempt + 1}): {e}")
                self.connect()

        raise ConnectionError(f"Could not list {remote_dir}")

    def retrieve(self, remote_dir, fname, local_path, retries=3):
        tmp_path = Path(f"{local_path}.part")
        for attempt in range(retries):
            try:
                self.ftp.cwd(remote_dir)
                with open(tmp_path, "wb") as f:
                    self.ftp.retrbinary(f"RETR {fname}", f.write)
                tmp_path.replace(local_path)
                return

            except Exception as e:
                print(f"WARNING {fname} (attempt {attempt + 1}): {e}")
                tmp_path.unlink(missing_ok=True)
                self.connect()

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


def _sampling_s(fname):
    """Sampling interval (s) of a RINEX 3/4 name; None for RINEX 2."""
    m = _RE_R3.match(fname)
    if not m:
        return None
    return int(m.group(2)) * _UNIT_S[m.group(3).upper()]


def _r3_rank(fname):
    """Sort key: 30 s first, then 15 s, then coarser rates, then finer ones."""
    s = _sampling_s(fname)
    if s == 30:
        return (0, 0, fname)
    if s == 15:
        return (1, 0, fname)
    if s > 30:
        return (2, s, fname)
    return (3, -s, fname)


def _ordered_candidates(r3, r2):
    """Candidates in priority order (see module docstring); duplicate
    RINEX 2 copies (.gz and .Z of the same file) are kept only once."""
    out, seen = [], set()
    for fname in sorted(r3, key=_r3_rank) + sorted(r2, key=lambda f: not f.endswith(".gz")):
        if _local_name(fname) in seen:
            continue
        seen.add(_local_name(fname))
        out.append(fname)
    return out


def _best_file(conn, remote_dir, r3_candidates, r2_candidates, output_dir):
    """First candidate, in priority order, that is cached or downloads correctly."""
    for fname in _ordered_candidates(r3_candidates, r2_candidates):
        local = output_dir / _local_name(fname)
        if local.exists():
            print(f"SKIP  {local.name}")
            return local
        path = _download(conn, remote_dir, fname, output_dir)
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

def stations_available_ign(start, end, full_id=False):
    """
    Return the union of stations having daily RINEX observation files
    on the IGN archive during the requested period [start, end).

    Archive:
        ftp://igs.ign.fr/pub/igs/data/YYYY/DDD/

    Parameters
    ----------
    full_id : bool
        False -> 4-character names (lower case), like stations_available_cors.
        True  -> 9-character ids (e.g. 'LPGS00ARG'). Stations that only
                 have RINEX 2 files are returned with their 4-character name.

    Returns
    -------
    sorted list of str
    """
    list_stations = set()
    conn = _IgnFTP()

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


def list_files_ign(station, start, end):
    """
    List the observation files of one station available on the IGN
    archive between start and end (inclusive).

    station: 4-character site name or 9-character id (case-insensitive).

    Returns
    -------
    list of (date, filename)   (RINEX 3/4 first for each day)
    """
    station = _check_station(station)
    found = []
    conn = _IgnFTP()

    try:
        for d, year, doy in _iter_days(start, end):
            try:
                names = conn.list_day(year, doy)
            except Exception as e:
                print(f"ERROR {d}: {e}")
                continue

            r3, r2 = _split_candidates(names, station)
            found += [(d, f) for f in _ordered_candidates(r3, r2)]

    finally:
        conn.close()

    return found


def download_ign(station, start, end, output_dir="."):
    """
    Download IGN daily RINEX observation files, one per day,
    with the highest RINEX version available (RINEX 3/4 over RINEX 2).

    Archive structure:
        ftp://igs.ign.fr/pub/igs/data/YYYY/DDD/

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
    conn = _IgnFTP()

    try:
        for d, year, doy in _iter_days(start, end):
            remote_dir = conn.day_dir(year, doy)

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
    Download ONE daily RINEX observation file for a station from IGN,
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

    conn = _IgnFTP()
    remote_dir = conn.day_dir(year, doy)

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
