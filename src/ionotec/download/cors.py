from ftplib import FTP
from pathlib import Path

import subprocess
import requests
import re
import time

from datetime import timedelta
from ionotec.download.extract import extract_and_remove_gz

def stations_available_cors(start, end):
    """
    Return the union of 4-character NOAA CORS stations
    having daily Hatanaka RINEX files during the requested period.

    Archive:
        https://noaa-cors-pds.s3.amazonaws.com/rinex/YYYY/DDD/

    File:
        ssssDDD0.YYd.gz
    """

    list_stations = set()

    d = start

    with requests.Session() as session:

        while d <= end:

            year = d.year
            doy = d.timetuple().tm_yday

            prefix = f"rinex/{year}/{doy:03d}/"

            continuation_token = None

            while True:

                params = {
                    "list-type": "2",
                    "prefix": prefix,
                }

                if continuation_token is not None:
                    params["continuation-token"] = continuation_token

                try:
                    r = session.get(
                        "https://noaa-cors-pds.s3.amazonaws.com/",
                        params=params,
                        timeout=60
                    )
                    r.raise_for_status()

                except Exception as e:
                    print(f"ERROR {d}: {e}")
                    break

                # Extract object keys
                keys = re.findall(
                    r"<Key>([^<]+)</Key>",
                    r.text
                )

                # Expected:
                # rinex/YYYY/DDD/ssss/ssssDDD0.YYd.gz
                pattern = re.compile(
                    rf"^rinex/{year}/{doy:03d}/"
                    rf"([A-Za-z0-9]{{4}})/"
                    rf"[A-Za-z0-9]{{4}}{doy:03d}0\.\d{{2}}d\.gz$",
                    re.IGNORECASE
                )

                for key in keys:

                    match = pattern.match(key)

                    if match:
                        list_stations.add(
                            match.group(1).lower()
                        )

                # Is there another page?
                if "<IsTruncated>true</IsTruncated>" not in r.text:
                    break

                match = re.search(
                    r"<NextContinuationToken>([^<]+)</NextContinuationToken>",
                    r.text
                )

                if not match:
                    print(
                        f"WARNING: no continuation token for {d}"
                    )
                    break

                continuation_token = match.group(1)

            d += timedelta(days=1)

    return sorted(list_stations)


def download_cors(station, start, end, output_dir="."):
    """
    Download NOAA CORS RINEX 2 daily Hatanaka files.

    Archive structure:
        https://noaa-cors-pds.s3.amazonaws.com/rinex/YYYY/DDD/ssss/

    Filename:
        ssssDDD0.YYd.gz

    Example:
        https://noaa-cors-pds.s3.amazonaws.com/rinex/2024/096/akrd/akrd0960.24d.gz

    The downloaded .gz files are extracted using
    extract_and_remove_gz().

    Returns
    -------
    list of Path
        Locally available RINEX files.
    """

    station = station.lower()

    if len(station) != 4:
        raise ValueError("station must be a 4-character station name")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    list_files = []

    d = start

    with requests.Session() as session:

        while d <= end:

            year = d.year
            doy = d.timetuple().tm_yday
            yy = year % 100

            fname = f"{station}{doy:03d}0.{yy:02d}d.gz"

            url = (
                f"https://noaa-cors-pds.s3.amazonaws.com/"
                f"rinex/{year}/{doy:03d}/{station}/{fname}"
            )

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

                with session.get(
                    url,
                    stream=True,
                    timeout=120
                ) as response:

                    if response.status_code == 404:
                        print(
                            f"NONE  {station} on {d} "
                            f"(DOY {doy:03d})"
                        )
                        d += timedelta(days=1)
                        continue

                    response.raise_for_status()

                    with open(gz_path, "wb") as f:
                        for chunk in response.iter_content(
                            chunk_size=1024 * 1024
                        ):
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


# True  -> a downloaded RINEX 2 Hatanaka file (.YYd) is stored uncompressed
#          (.YYo) with `crx2rnx`, as in igs.py.
# False -> keep the Hatanaka file (.YYd).
CONVERT_RINEX2_TO_PLAIN = True


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


def download_file(station, date, output_dir="."):
    """
    Download ONE NOAA CORS daily RINEX 2 observation file of a station.

    station : 4-character site name (case-insensitive)
    date    : datetime / date of the day wanted

    Archive:
        https://noaa-cors-pds.s3.amazonaws.com/rinex/YYYY/DDD/ssss/ssssDDD0.YYd.gz

    This archive only provides RINEX 2 Hatanaka files, so there is no
    version to choose from. The .gz is extracted with extract_and_remove_gz();
    if CONVERT_RINEX2_TO_PLAIN is True (default) the Hatanaka file (.YYd) is
    converted to plain RINEX 2 (.YYo) with `crx2rnx` (must be in the PATH)
    and the .YYd removed.

    Returns
    -------
    Path of the final file (e.g. .../akrd0960.24o, or .24d if
    CONVERT_RINEX2_TO_PLAIN is False), or None if the file does not exist
    or on error.
    """

    station = station.lower()

    if len(station) != 4:
        raise ValueError("station must be a 4-character station name")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    year = date.year
    doy = date.timetuple().tm_yday
    yy = year % 100

    fname = f"{station}{doy:03d}0.{yy:02d}d.gz"

    url = (
        f"https://noaa-cors-pds.s3.amazonaws.com/"
        f"rinex/{year}/{doy:03d}/{station}/{fname}"
    )

    gz_path = output_dir / fname
    extracted_path = gz_path.with_suffix("")                  # ssssDDD0.YYd

    final_path = (
        extracted_path.with_suffix(f".{yy:02d}o")             # ssssDDD0.YYo
        if CONVERT_RINEX2_TO_PLAIN
        else extracted_path
    )

    # Already downloaded
    if final_path.exists():
        print(f"SKIP  {final_path.name}")
        return final_path

    try:
        print(f"GET   {fname}")

        with requests.Session() as session:

            for attempt in range(3):
                try:
                    response = session.get(url, stream=True, timeout=120)
                    break

                except (
                    requests.exceptions.ConnectionError,
                    requests.exceptions.Timeout,
                ) as e:
                    if attempt == 2:
                        raise

                    print(f"WARNING {fname}: {e} -> retry in {2 * (attempt + 1)} s")
                    time.sleep(2 * (attempt + 1))

            with response:

                if response.status_code == 404:
                    print(f"NONE  {station} on {date} (DOY {doy:03d})")
                    return None

                response.raise_for_status()

                with open(gz_path, "wb") as f:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            f.write(chunk)

        extracted = extract_and_remove_gz(gz_path)

        if CONVERT_RINEX2_TO_PLAIN:
            extracted = _crx2rnx(extracted)

        return extracted

    except Exception as e:

        # leave nothing half-done (.gz, extracted .YYd, .YYo)
        for p in (gz_path, extracted_path, final_path):
            p.unlink(missing_ok=True)

        print(f"ERROR {fname}: {e}")
        return None
