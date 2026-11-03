from ftplib import FTP, error_perm
from pathlib import Path
from datetime import timedelta

import re

from ionotec.download.extract import extract_and_remove_gz

FTP_HOST = "ftp.sonel.org"
FTP_ROOT = "/gps/data"

# RINEX 3 long name:  ALBO00ESP_R_20240100000_01D_30S_MO.crx.gz
_RE_R3 = re.compile(r"^([A-Z0-9]{9})_[A-Z]_\d{11}_", re.IGNORECASE)
# RINEX 2 Hatanaka:   albo0100.24d.gz
_RE_R2 = re.compile(r"^([A-Z0-9]{4})\d{3}[0-9A-X]\.\d{2}d\.gz$", re.IGNORECASE)
# Navigation files (R3: _MN. / _GN. ..., R2: .24n, .24g ...)
_RE_NAV = re.compile(r"_[A-Z]N\.|\.\d{2}[nglhpq]\.", re.IGNORECASE)


# ----------------------------------------------------------------------
# FTP helpers
# ----------------------------------------------------------------------

class _SonelFTP:
    """FTP connection to SONEL that reconnects when the server drops it."""

    def __init__(self, timeout=60):
        self.timeout = timeout
        self.ftp = None
        self.connect()

    def connect(self):
        self.close()
        self.ftp = FTP(FTP_HOST, timeout=self.timeout)
        self.ftp.login()

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

    def list_day(self, year, doy, retries=3):
        """
        Names of the files in /gps/data/YYYY/DDD.
        Returns None if the directory does not exist, [] if it is empty.
        """
        remote_dir = f"{FTP_ROOT}/{year}/{doy:03d}"

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
        for attempt in range(retries):
            try:
                self.ftp.cwd(remote_dir)
                with open(local_path, "wb") as f:
                    self.ftp.retrbinary(f"RETR {fname}", f.write)
                return

            except Exception as e:
                print(f"WARNING {fname} (attempt {attempt + 1}): {e}")
                Path(local_path).unlink(missing_ok=True)
                self.connect()

        raise ConnectionError(f"Could not download {fname}")


def _is_obs_file(fname):
    """True for compressed daily observation files (R3 or R2 Hatanaka)."""
    if not fname.lower().endswith(".gz"):
        return False
    if _RE_NAV.search(fname):
        return False
    return bool(_RE_R3.match(fname) or _RE_R2.match(fname))


def _station_ids(fname):
    """
    Return (full_id, site4) in lower case from a file name.
    full_id is the 9-character id (RINEX 3) or None (RINEX 2).
    """
    m = _RE_R3.match(fname)
    if m:
        full = m.group(1).lower()
        return full, full[:4]

    m = _RE_R2.match(fname)
    if m:
        return None, m.group(1).lower()

    return None, None


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------

def stations_available_sonel(start, end, full_id=False):
    """
    Return the union of stations having daily RINEX observation files
    on the SONEL FTP server during the requested period.

    Archive:
        ftp://ftp.sonel.org/gps/data/YYYY/DDD/

    Parameters
    ----------
    full_id : bool
        False -> 4-character names (lower case), like stations_available_cors.
        True  -> 9-character ids (e.g. 'ALBO00ESP') as used in
                 process_sonel.list_stations. RINEX 2 files only have
                 the 4-character name, so they are returned as is.

    Returns
    -------
    sorted list of str
    """

    list_stations = set()

    conn = _SonelFTP()

    try:
        d = start

        while d <= end:

            year = d.year
            doy = d.timetuple().tm_yday

            try:
                names = conn.list_day(year, doy)
            except Exception as e:
                print(f"ERROR {d}: {e}")
                d += timedelta(days=1)
                continue

            if names is None:
                print(f"NONE  no directory for {d} (DOY {doy:03d})")
                d += timedelta(days=1)
                continue

            for fname in names:

                if not _is_obs_file(fname):
                    continue

                full, site = _station_ids(fname)

                if full_id:
                    list_stations.add((full or site).upper())
                else:
                    list_stations.add(site)

            d += timedelta(days=1)

    finally:
        conn.close()

    return sorted(list_stations)


def list_files_sonel(station, start, end):
    """
    List the observation files of one station available on the SONEL FTP
    server between start and end (inclusive).

    station: 4-character site name or 9-character id (case-insensitive).

    Returns
    -------
    list of (date, filename)
    """

    station = station.lower()

    if len(station) not in (4, 9):
        raise ValueError("station must be a 4 or 9 character name")

    found = []

    conn = _SonelFTP()

    try:
        d = start

        while d <= end:

            year = d.year
            doy = d.timetuple().tm_yday

            try:
                names = conn.list_day(year, doy)
            except Exception as e:
                print(f"ERROR {d}: {e}")
                d += timedelta(days=1)
                continue

            for fname in names or []:

                if not _is_obs_file(fname):
                    continue

                full, site = _station_ids(fname)

                if (len(station) == 9 and full == station) or \
                   (len(station) == 4 and site == station):
                    found.append((d, fname))

            d += timedelta(days=1)

    finally:
        conn.close()

    return found


def download_sonel(station, start, end, output_dir="."):
    """
    Download SONEL daily RINEX observation files.

    Archive structure:
        ftp://ftp.sonel.org/gps/data/YYYY/DDD/

    station: 4-character site name or 9-character id (case-insensitive).
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

    conn = _SonelFTP()

    try:
        d = start

        while d <= end:

            year = d.year
            doy = d.timetuple().tm_yday

            remote_dir = f"{FTP_ROOT}/{year}/{doy:03d}"

            try:
                names = conn.list_day(year, doy)
            except Exception as e:
                print(f"ERROR accessing {remote_dir}: {e}")
                d += timedelta(days=1)
                continue

            found = False

            for fname in names or []:

                if not _is_obs_file(fname):
                    continue

                full, site = _station_ids(fname)

                if len(station) == 9:
                    if full != station:
                        continue
                elif site != station:
                    continue

                found = True

                gz_path = output_dir / fname
                rinex_path = gz_path.with_suffix("")

                # Already extracted
                if rinex_path.exists():
                    print(f"SKIP  {rinex_path.name}")
                    list_files.append(rinex_path)
                    continue

                try:
                    print(f"GET   {fname}")

                    conn.retrieve(remote_dir, fname, gz_path)

                    rinex_path = extract_and_remove_gz(gz_path)
                    list_files.append(rinex_path)

                except Exception as e:

                    if gz_path.exists():
                        gz_path.unlink()

                    print(f"ERROR {fname}: {e}")

            if not found:
                print(f"NONE  {station} on {d} (DOY {doy:03d})")

            d += timedelta(days=1)

    finally:
        conn.close()

    print()
    print(f"Files available locally: {len(list_files)}")

    return list_files

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
    gz_path = output_dir / fname
    try:
        print(f"GET   {fname}")
        conn.retrieve(remote_dir, fname, gz_path)
        return extract_and_remove_gz(gz_path)
    except Exception as e:
        gz_path.unlink(missing_ok=True)
        gz_path.with_suffix("").unlink(missing_ok=True)
        print(f"ERROR {fname}: {e}")
        return None


def download_file(station, date, output_dir="."):
    """
    Download ONE daily RINEX observation file for a station from SONEL,
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
    remote_dir = f"{FTP_ROOT}/{year}/{doy:03d}"

    conn = _SonelFTP()

    try:
        try:
            names = conn.list_day(year, doy)
        except Exception as e:
            print(f"ERROR accessing {remote_dir}: {e}")
            return None

        if names is None:
            print(f"NONE  no directory for {date} (DOY {doy:03d})")
            return None

        r3_candidates = []
        r2_candidates = []

        for fname in names:
            if not _is_obs_file(fname):
                continue

            full, site = _station_ids(fname)
            if site != station:
                continue

            if full:  # RINEX 3/4 long name (9-character ID)
                r3_candidates.append(fname)
            else:     # RINEX 2 short name
                r2_candidates.append(fname)

        # Prioritize RINEX 3/4 first, fall back to RINEX 2 if none available
        candidates = r3_candidates + r2_candidates

        if not candidates:
            print(f"NONE  {station} on {date} (DOY {doy:03d})")
            return None

        # Take the top priority matching file
        fname = candidates[0]
        gz_path = output_dir / fname
        rinex_path = gz_path.with_suffix("")

        # Already extracted locally
        if rinex_path.exists():
            print(f"SKIP  {rinex_path.name}")
            return rinex_path

        try:
            print(f"GET   {fname}")
            conn.retrieve(remote_dir, fname, gz_path)

            rinex_path = extract_and_remove_gz(gz_path)
            return rinex_path

        except Exception as e:
            if gz_path.exists():
                gz_path.unlink()
            rinex_path.unlink(missing_ok=True)

            print(f"ERROR {fname}: {e}")
            return None

    finally:
        conn.close()
