from pathlib import Path
from datetime import timedelta
from urllib.parse import urljoin

import re
import time
import requests

from ionotec.download.extract import extract_and_remove_gz

BASE_URL = "https://rinex.geodetic.gov.hk"

# RINEX 3 daily Hatanaka, 30 s:
#   HKQT00HKG_R_20241160000_01D_30S_MO.crx.gz
_RE_OBS = re.compile(
    r"^([A-Z0-9]{9})_[A-Z]_\d{11}_\d{2}[A-Z]_\d{2}[A-Z]_MO\.crx\.gz$",
    re.IGNORECASE,
)
# Same pattern, to find file names anywhere in the directory page
_RE_FIND = re.compile(
    r"(?<![A-Za-z0-9])[A-Za-z0-9]{9}_[A-Za-z]_\d{11}_\d{2}[A-Za-z]_\d{2}[A-Za-z]"
    r"_MO\.crx\.gz(?![A-Za-z0-9_.])",
    re.IGNORECASE,
)


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


def stations_available_hongkong(start, end):
    """
    Return the union of 4-character Hong Kong GNSS stations
    available during the requested period.

    Archive structure:
        https://rinex.geodetic.gov.hk/rinex3/YYYY/DOY/ssss/30s/

    The station names are the sub-directories of
        https://rinex.geodetic.gov.hk/rinex3/YYYY/DOY/

    (a station directory may exist without a file for that day;
    download_hongkong prints NONE in that case).
    """

    list_stations = set()

    d = start

    with requests.Session() as session:

        while d < end:

            year = d.year
            doy = d.timetuple().tm_yday

            remote_dir = f"{BASE_URL}/rinex3/{year}/{doy:03d}/"

            try:
                r = _get(session, remote_dir)

            except Exception as e:
                print(f"ERROR {d}: {e}")
                d += timedelta(days=1)
                continue

            if r is None:
                print(f"NONE  no directory for {d}")
                d += timedelta(days=1)
                continue

            with r:
                text = r.text

            # Station sub-directories, e.g.
            # HREF="/rinex3/2024/100/hkpc/"
            stations = re.findall(
                rf'HREF="[^"]*/rinex3/{year}/{doy:03d}/([A-Za-z0-9]{{4}})/?"',
                text,
                re.IGNORECASE,
            )

            list_stations.update(s.lower() for s in stations)

            d += timedelta(days=1)

    return sorted(list_stations)


def download_hongkong(station, start, end, output_dir="."):
    """
    Download Hong Kong daily RINEX 3 Hatanaka files (30 s).

    Archive structure:
        https://rinex.geodetic.gov.hk/rinex3/YYYY/DOY/ssss/30s/

    Filename:
        SSSS00CCC_R_YYYYDDD0000_01D_30S_MO.crx.gz

    Example:
        HKQT00HKG_R_20241160000_01D_30S_MO.crx.gz

    station: 4-character site name or 9-character id (case-insensitive).

    The downloaded .gz files are extracted using extract_and_remove_gz();
    the files are kept as Hatanaka (.crx).

    Returns
    -------
    list of Path
        Locally available RINEX files.
    """

    station = station.lower()

    if len(station) not in (4, 9):
        raise ValueError("station must be a 4 or 9 character name")

    site = station[:4]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    list_files = []

    d = start

    with requests.Session() as session:

        while d <= end:

            year = d.year
            doy = d.timetuple().tm_yday

            remote_dir = f"{BASE_URL}/rinex3/{year}/{doy:03d}/{site}/30s/"

            try:
                r = _get(session, remote_dir)

            except Exception as e:
                print(f"ERROR listing {remote_dir}: {e}")
                d += timedelta(days=1)
                continue

            if r is None:
                print(f"NONE  {station} on {d} (DOY {doy:03d})")
                d += timedelta(days=1)
                continue

            with r:
                text = r.text

            # {file name: url}, from the exact hrefs of the server
            files = {}

            for href in re.findall(r'HREF="([^"]+)"', text, re.IGNORECASE):
                fname = href.split("?")[0].rstrip("/").split("/")[-1]

                if _RE_OBS.match(fname):
                    files[fname] = urljoin(remote_dir, href)

            for fname in _RE_FIND.findall(text):
                files.setdefault(fname, remote_dir + fname)

            # Keep the requested station (and day)
            names = [
                f for f in files
                if f.lower().startswith(station)
                and f"_{year}{doy:03d}" in f
            ]

            if not names:
                print(f"NONE  {station} on {d} (DOY {doy:03d})")
                d += timedelta(days=1)
                continue

            for fname in names:

                gz_path = output_dir / fname
                rinex_path = gz_path.with_suffix("")

                # Already extracted
                if rinex_path.exists():
                    print(f"SKIP  {rinex_path.name}")
                    list_files.append(rinex_path)
                    continue

                try:
                    print(f"GET   {fname}")

                    r = _get(session, files[fname], stream=True)

                    if r is None:
                        print(f"NONE  {fname}")
                        continue

                    with r:
                        with open(gz_path, "wb") as f:
                            for chunk in r.iter_content(
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



def download_file(station, date, output_dir="."):
    """
    Télécharge UN SEUL fichier d'observation RINEX 3 Hatanaka journalier (30 s) de Hong Kong.

    Parameters
    ----------
    station : str
        Nom de station à 4 caractères ou identifiant à 9 caractères (insensible à la casse).
    date : datetime ou date
        La date cible pour le fichier d'observation.
    output_dir : str ou Path, optional
        Répertoire où le fichier sera enregistré. Par défaut, le répertoire courant.

    Returns
    -------
    Path ou None
        Chemin du fichier final extrait, ou None si aucun fichier n'a été trouvé ou en cas d'erreur.
    """
    station = station.lower()

    if len(station) not in (4, 9):
        raise ValueError("station must be a 4 or 9 character name")

    site = station[:4]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    year = date.year
    doy = date.timetuple().tm_yday

    remote_dir = f"{BASE_URL}/rinex3/{year}/{doy:03d}/{site}/30s/"

    with requests.Session() as session:
        try:
            r = _get(session, remote_dir)
        except Exception as e:
            print(f"ERROR listing {remote_dir}: {e}")
            return None

        if r is None:
            print(f"NONE  {station} on {date} (DOY {doy:03d})")
            return None

        with r:
            text = r.text

        # {nom_de_fichier : url}, récupéré depuis les liens de la page HTML
        files = {}

        for href in re.findall(r'HREF="([^"]+)"', text, re.IGNORECASE):
            fname = href.split("?")[0].rstrip("/").split("/")[-1]

            if _RE_OBS.match(fname):
                files[fname] = urljoin(remote_dir, href)

        for fname in _RE_FIND.findall(text):
            files.setdefault(fname, remote_dir + fname)

        # Filtrer pour ne garder que la station et la date requises
        names = [
            f for f in files
            if f.lower().startswith(station)
            and f"_{year}{doy:03d}" in f
        ]

        if not names:
            print(f"NONE  {station} on {date} (DOY {doy:03d})")
            return None

        # On prend le premier fichier correspondant
        fname = names[0]
        gz_path = output_dir / fname
        rinex_path = gz_path.with_suffix("")

        # Déjà extrait localement
        if rinex_path.exists():
            print(f"SKIP  {rinex_path.name}")
            return rinex_path

        try:
            print(f"GET   {fname}")

            r = _get(session, files[fname], stream=True)

            if r is None:
                print(f"NONE  {fname}")
                return None

            with r:
                with open(gz_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            f.write(chunk)

            # Extraction du fichier .gz vers .crx
            rinex_path = extract_and_remove_gz(gz_path)
            return rinex_path

        except Exception as e:
            if gz_path.exists():
                gz_path.unlink()
            rinex_path.unlink(missing_ok=True)

            print(f"ERROR {fname}: {e}")
            return None
