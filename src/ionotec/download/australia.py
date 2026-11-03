from pathlib import Path
from datetime import timedelta

import re
import time
import requests

from ionotec.download.extract import extract_and_remove_gz

# Geoscience Australia GNSS archive (public S3 bucket)
BUCKET_URL = "https://ga-gnss-data-rinex-v1.s3.amazonaws.com"
DAILY_PREFIX = "public/daily"

# RINEX 3/4 long name, daily observation, Hatanaka or plain:
#   ABMF00GLP_R_20160160000_01D_30S_MO.crx.gz
_RE_LONG = re.compile(
    r"^([A-Z0-9]{9})_[A-Z]_\d{11}_\d{2}[A-Z]_\d{2}[A-Z]_MO\.(?:crx|rnx)\.gz$",
    re.IGNORECASE,
)
# RINEX 2 short name, Hatanaka:
#   abmf0160.16d.gz
_RE_SHORT = re.compile(
    r"^([A-Z0-9]{4})\d{3}0\.\d{2}d\.gz$",
    re.IGNORECASE,
)


def _day_prefix(d):
    return f"{DAILY_PREFIX}/{d.year}/{d.timetuple().tm_yday:03d}/"


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


def _list_keys(session, prefix):
    """All object keys of the bucket starting with `prefix` (S3 listing)."""
    keys = []
    token = None

    while True:
        params = {"list-type": "2", "prefix": prefix}

        if token is not None:
            params["continuation-token"] = token

        r = _get(session, BUCKET_URL + "/", params=params)

        if r is None:
            return keys

        with r:
            text = r.text

        keys += re.findall(r"<Key>([^<]+)</Key>", text)

        if "<IsTruncated>true</IsTruncated>" not in text:
            return keys

        m = re.search(
            r"<NextContinuationToken>([^<]+)</NextContinuationToken>", text
        )

        if not m:
            print(f"WARNING no continuation token for {prefix}")
            return keys

        token = m.group(1)


def _station_ids(fname):
    """(full_id, site4) in lower case, or (None, None). full_id is None
    for RINEX 2 short names."""
    m = _RE_LONG.match(fname)
    if m:
        full = m.group(1).lower()
        return full, full[:4]

    m = _RE_SHORT.match(fname)
    if m:
        return None, m.group(1).lower()

    return None, None


def _is_obs(fname):
    return bool(_RE_LONG.match(fname) or _RE_SHORT.match(fname))


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


def stations_available_australia(start, end, full_id=False):
    """
    Return the union of stations having daily RINEX observation files
    in the Geoscience Australia archive during the requested period
    (RINEX 4/3 long names and RINEX 2 short names).

    Archive:
        https://ga-gnss-data-rinex-v1.s3.amazonaws.com/public/daily/YYYY/DDD/

    File:
        SSSS00CCC_R_YYYYDDD0000_01D_30S_MO.crx.gz   (RINEX 3/4)
        ssssDDD0.YYd.gz                             (RINEX 2)

    Parameters
    ----------
    full_id : bool
        False -> 4-character lower-case names (like the other sources).
        True  -> 9-character upper-case ids when the station has a long
                 name; RINEX 2-only stations are returned with 4 characters.
    """

    list_stations = set()

    d = start

    with requests.Session() as session:

        while d < end:

            try:
                keys = _list_keys(session, _day_prefix(d))

            except Exception as e:
                print(f"ERROR {d}: {e}")
                d += timedelta(days=1)
                continue

            if not keys:
                print(f"NONE  no files for {d}")

            for key in keys:

                fname = key.split("/")[-1]

                if not _is_obs(fname):
                    continue

                full, site = _station_ids(fname)

                if full_id:
                    list_stations.add((full or site).upper())
                else:
                    list_stations.add(site)

            d += timedelta(days=1)

    if full_id:
        # a station seen with a 9-character id does not also need its
        # 4-character RINEX 2 entry
        longs = {s[:4] for s in list_stations if len(s) == 9}
        list_stations = {
            s for s in list_stations if len(s) == 9 or s[:4].lower() not in
            {l.lower() for l in longs}
        }

    return sorted(list_stations)


def _download(session, d, fname, output_dir):
    """Download + extract one file. Returns the extracted Path, or None."""
    gz_path = output_dir / fname
    url = f"{BUCKET_URL}/{_day_prefix(d)}{fname}"

    try:
        print(f"GET   {fname}")

        r = _get(session, url, stream=True)

        if r is None:
            print(f"NONE  {fname}")
            return None

        with r:
            with open(gz_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        f.write(chunk)

        return extract_and_remove_gz(gz_path)

    except Exception as e:

        gz_path.unlink(missing_ok=True)
        gz_path.with_suffix("").unlink(missing_ok=True)

        print(f"ERROR {fname}: {e}")
        return None


def download_australia(station, start, end, output_dir="."):
    """
    Download Geoscience Australia daily RINEX observation files.

    Archive structure:
        https://ga-gnss-data-rinex-v1.s3.amazonaws.com/public/daily/YYYY/DDD/

    station: 4-character site name or 9-character id (case-insensitive).
    (RINEX 2 short names only exist for 4-character names.)

    Version preference, per station and day:
        1. RINEX 4, then RINEX 3  (long names SSSS00CCC_R_..._MO.crx.gz).
           The file name does not tell v4 from v3: if several long-name
           files exist for the day they are all downloaded, the version
           is read from the header, the highest is kept and the others
           are deleted. A single long-name file is taken as is.
        2. RINEX 2 (ssssDDD0.YYd.gz) only if no long-name file could be
           obtained.

    The downloaded .gz files are extracted using extract_and_remove_gz();
    the files are kept as Hatanaka (.crx / .YYd).

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

            doy = d.timetuple().tm_yday

            # List only this station's files (names may be upper or lower
            # case, and S3 prefixes are case-sensitive)
            names = set()

            try:
                for case in {site.upper(), site.lower()}:
                    for key in _list_keys(session, _day_prefix(d) + case):
                        names.add(key.split("/")[-1])

            except Exception as e:
                print(f"ERROR listing {d}: {e}")
                d += timedelta(days=1)
                continue

            long_names, short_names = [], []

            for name in sorted(names):

                full, s = _station_ids(name)

                if s != site:
                    continue

                if _RE_LONG.match(name):
                    if len(station) == 4 or full == station:
                        long_names.append(name)
                else:
                    if len(station) == 4:
                        short_names.append(name)

            if not long_names and not short_names:
                print(f"NONE  {station} on {d} (DOY {doy:03d})")
                d += timedelta(days=1)
                continue

            # Already extracted (an earlier run keeps only the best one)
            local = None
            for name in long_names + short_names:
                p = (output_dir / name).with_suffix("")
                if p.exists():
                    local = p
                    break

            if local is not None:
                print(f"SKIP  {local.name}")
                list_files.append(local)
                d += timedelta(days=1)
                continue

            chosen = None

            # Tier 1: RINEX 4 / 3 (long names)
            if long_names:

                got = []

                for name in long_names:
                    p = _download(session, d, name, output_dir)

                    if p is not None:
                        got.append(p)

                    # a single candidate needs no version check
                    if len(long_names) == 1:
                        break

                if len(got) == 1:
                    chosen = got[0]

                elif len(got) > 1:
                    # highest version wins (ties: first in name order)
                    versions = [_rinex_version(p) or 0.0 for p in got]
                    best = max(range(len(got)), key=lambda i: versions[i])
                    chosen = got[best]

                    print(
                        f"      versions {dict(zip([p.name for p in got], versions))}"
                        f" -> keeping {chosen.name}"
                    )

                    for p in got:
                        if p != chosen:
                            p.unlink(missing_ok=True)

            # Tier 2: RINEX 2
            if chosen is None and short_names:
                chosen = _download(session, d, short_names[0], output_dir)

            if chosen is not None:
                list_files.append(chosen)

            d += timedelta(days=1)

    print()
    print(f"Files available locally: {len(list_files)}")

    return list_files


def download_file(station, date, output_dir="."):
    """
    Download ONE daily Geoscience Australia RINEX observation file for a station.

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

    site = station[:4]
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    doy = date.timetuple().tm_yday

    with requests.Session() as session:
        # List only this station's files (S3 prefixes are case-sensitive)
        names = set()

        try:
            for case in {site.upper(), site.lower()}:
                for key in _list_keys(session, _day_prefix(date) + case):
                    names.add(key.split("/")[-1])

        except Exception as e:
            print(f"ERROR listing {date}: {e}")
            return None

        long_names, short_names = [], []

        for name in sorted(names):
            full, s = _station_ids(name)

            if s != site:
                continue

            if _RE_LONG.match(name):
                if len(station) == 4 or full == station:
                    long_names.append(name)
            else:
                if len(station) == 4:
                    short_names.append(name)

        if not long_names and not short_names:
            print(f"NONE  {station} on {date} (DOY {doy:03d})")
            return None

        # Check if already locally extracted (best candidate first)
        for name in long_names + short_names:
            p = (output_dir / name).with_suffix("")
            if p.exists():
                print(f"SKIP  {p.name}")
                return p

        chosen = None

        # Tier 1: RINEX 4 / 3 (long names)
        if long_names:
            got = []

            for name in long_names:
                p = _download(session, date, name, output_dir)
                if p is not None:
                    got.append(p)

                # A single candidate needs no version check
                if len(long_names) == 1:
                    break

            if len(got) == 1:
                chosen = got[0]

            elif len(got) > 1:
                # Highest version wins (ties: first in name order)
                versions = [_rinex_version(p) or 0.0 for p in got]
                best = max(range(len(got)), key=lambda i: versions[i])
                chosen = got[best]

                print(
                    f"      versions {dict(zip([p.name for p in got], versions))}"
                    f" -> keeping {chosen.name}"
                )

                for p in got:
                    if p != chosen:
                        p.unlink(missing_ok=True)

        # Tier 2: RINEX 2
        if chosen is None and short_names:
            chosen = _download(session, date, short_names[0], output_dir)

        if chosen is None:
            print(f"NONE  {station} on {date} (DOY {doy:03d})")

        return chosen
