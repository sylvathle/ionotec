from ftplib import FTP
from pathlib import Path

import subprocess

from datetime import timedelta
from ionotec.download.extract import extract_and_remove_gz

FTP_HOST = "ftp.earthobservatory.sg"

# Compression extensions found on the server, in order of preference
COMPRESSIONS = (".gz", ".Z")


def stations_available_sugar(start, end):

    ftp = FTP(FTP_HOST, timeout=60)
    ftp.login()

    list_stations = set()

    for year in range(start.year, end.year + 1):

        remote_dir = f"/SugarData/{year}"

        try:
            ftp.cwd(remote_dir)
            stations = ftp.nlst()

        except Exception as e:
            print(f"ERROR {year}: {e}")
            continue

        for station in stations:

            station = station.split("/")[-1]

            # SuGAr stations are 4-character RINEX v2 codes
            if len(station) != 4:
                continue

            list_stations.add(station)

    ftp.quit()

    return sorted(list_stations)


def extract_and_remove_z(z_path):
    """Decompress a Unix-compress (.Z, LZW) file and remove the archive.

    Uses `gzip -dc` (gzip reads .Z natively); falls back to the `unlzw3`
    package if gzip is not available.

    Returns
    -------
    pathlib.Path
        Path of the decompressed file (archive name without the .Z suffix).
    """
    z_path = Path(z_path)
    out_path = z_path.with_suffix("")

    try:
        try:
            with out_path.open("wb") as f:
                subprocess.run(
                    ["gzip", "-dc", str(z_path)],
                    stdout=f, stderr=subprocess.PIPE, check=True,
                )
        except (FileNotFoundError, subprocess.CalledProcessError):
            # gzip missing or refused the file -> pure-Python fallback
            import unlzw3  # pip install unlzw3
            out_path.write_bytes(unlzw3.unlzw(z_path.read_bytes()))
    except Exception:
        out_path.unlink(missing_ok=True)
        raise

    z_path.unlink()
    return out_path


def _extract(path):
    """Extract a downloaded archive according to its extension."""
    if path.suffix == ".Z":
        return extract_and_remove_z(path)
    return extract_and_remove_gz(path)


def _fetch_day(ftp, station, date, output_dir):
    """Download + extract one daily file on an already-logged-in FTP session.

    Tries each extension in COMPRESSIONS (.gz, then .Z).
    Returns the extracted path, or None if unavailable.
    """
    year = date.year
    doy = date.timetuple().tm_yday
    base = f"{station}{doy:03d}0.{year % 100:02d}d"
    output_path = output_dir / base

    if output_path.exists():
        print(f"SKIP  {output_path.name}")
        return output_path

    remote_dir = f"/SugarData/{year}/{station}/RINEX2-HATANAKA"
    try:
        ftp.cwd(remote_dir)
    except Exception as e:
        print(f"ERROR accessing {remote_dir}: {e}")
        return None

    last_err = None
    for ext in COMPRESSIONS:
        fname = base + ext
        local = output_dir / fname

        try:
            print(f"GET   {fname}")
            with local.open("wb") as f:
                ftp.retrbinary(f"RETR {fname}", f.write)
        except Exception as e:
            # not on server with this extension (or transfer failed): try next
            local.unlink(missing_ok=True)
            last_err = e
            continue

        try:
            out = _extract(local)
            print(f"OK    {out.name}")
            return out
        except Exception as e:
            local.unlink(missing_ok=True)
            print(f"ERROR extracting {fname}: {e}")
            return None

    print(f"NONE  {station} on {date} (DOY {doy:03d}): {last_err}")
    return None


def download_file(station, date, output_dir="."):
    """Download one daily SuGAr RINEX 2 Hatanaka observation file.

    The remote file may be compressed as .gz or .Z; both are handled.

    Parameters
    ----------
    station : str
        Four-character SuGAr station code (case-insensitive).
    date : datetime.date or datetime.datetime
        Date of the daily file to download.
    output_dir : str or pathlib.Path
        Directory where the extracted file is saved.

    Returns
    -------
    pathlib.Path or None
        Path to the extracted .YYd file, or None if unavailable or if the
        download/extraction fails.
    """
    station = station.lower()

    if len(station) != 4:
        raise ValueError(
            f"SuGAr station must be a 4-character RINEX v2 code: {station}"
        )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ftp = None
    try:
        ftp = FTP(FTP_HOST, timeout=60)
        ftp.login()
        return _fetch_day(ftp, station, date, output_dir)

    except Exception as e:
        print(f"NONE  {station} on {date}: {e}")
        return None

    finally:
        if ftp is not None:
            try:
                ftp.quit()
            except Exception:
                ftp.close()


def download_sugar(station, start, end, output_dir="."):
    """Download daily files for [start, end) (end exclusive) on one session.

    Returns the list of extracted file paths.
    """
    station = station.lower()

    if len(station) != 4:
        raise ValueError(
            f"SuGAr station must be a 4-character RINEX v2 code: {station}"
        )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    list_files = []

    ftp = FTP(FTP_HOST, timeout=60)
    ftp.login()

    try:
        d = start
        while d < end:
            path = _fetch_day(ftp, station, d, output_dir)
            if path is not None:
                list_files.append(path)
            d += timedelta(days=1)
    finally:
        try:
            ftp.quit()
        except Exception:
            ftp.close()

    return list_files
