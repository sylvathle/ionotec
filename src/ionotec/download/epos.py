"""
Download daily RINEX observation files through the EPOS-GNSS data gateway.

EPOS-GNSS (European Plate Observing System) indexes, in one place, the daily
RINEX files of ~1,500 European stations hosted by many national data centres,
for example (2024-01-02):

    RENAG / RGP (France)          renag.resif.fr
    RING (INGV, Italy)            gnssdata-epos-ring.ingv.it
    FReDNet (OGS / CRS, Italy)    frednet.crs.ogs.it
    EUREF EPN                     epncb.oma.be
    SWEPOS, ESTPOS, OS Net, ...   datacenter.gnss-epos.eu
    Norway / Sweden / Spain tide-gauge stations   data.sonel.org
    + Spanish regional networks, Belgium, Portugal, Romania, Greece,
      French overseas (IPGP), Iceland, ...

Index (GLASS API, anonymous):
    https://gnssdata-epos.oca.eu/GlassFramework/webresources/files/combination/
        dateRange=YYYY-MM-DD,YYYY-MM-DD/json
Each entry gives the station marker (9 characters), the RINEX format
(RINEX2/3/4), sampling, size, md5 and the direct URL at the data centre.
The daily index is cached locally (~1 MB per day) so it is fetched only once.

Only 24 h / 30 s observation files are used. For each station and day:
    RINEX 4 > RINEX 3 > RINEX 2; when the same file is offered by several
    data centres, the next copy is tried if a download fails.

Same public API as sonel.py / nrcan.py.
"""

import json
import re
import shutil
import subprocess
import time
from datetime import timedelta
from pathlib import Path

import requests

from ionotec.download.extract import extract_and_remove_gz

GLASS_URL = "https://gnssdata-epos.oca.eu/GlassFramework/webresources/files/combination"
INDEX_CACHE = Path.home() / ".cache" / "ionotec" / "epos"

_FORMAT_RANK = {"RINEX4": 0, "RINEX3": 1, "RINEX2": 2}


# ----------------------------------------------------------------------
# GLASS index
# ----------------------------------------------------------------------

def _get(session, url, retries=3, timeout=(15, 300), **kw):
    for attempt in range(retries):
        try:
            r = session.get(url, timeout=timeout, **kw)
            r.raise_for_status()
            return r
        except requests.RequestException as e:
            if attempt == retries - 1:
                raise
            print(f"WARNING {url} (attempt {attempt + 1}): {e}")
            time.sleep(2 ** attempt)


def day_index(session, day, cache_dir=INDEX_CACHE):
    """
    All EPOS daily (24 h, 30 s) observation files of `day`, as a list of dicts
    {marker, format, url, name, size, md5}. Cached in cache_dir (one JSON per day).
    Days older than 3 days are cached permanently; recent days are re-queried.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{day:%Y%m%d}.json"

    recent = (time.time() - time.mktime(day.timetuple())) < 3 * 86400
    if path.exists() and not recent:
        return json.loads(path.read_text())

    url = f"{GLASS_URL}/dateRange={day:%Y-%m-%d},{day:%Y-%m-%d}/json"
    r = _get(session, url)
    entries = r.json() if r.content.strip() else []

    files = []
    for e in entries:
        ft = e.get("file_type") or {}
        if ft.get("sampling_window") != "24h" or ft.get("sampling_frequency") != "30s":
            continue
        if ft.get("format") not in _FORMAT_RANK:
            continue
        if (e.get("reference_date") or "")[:10] != f"{day:%Y-%m-%d}":
            continue
        files.append({
            "marker": e["station_marker"].upper(),
            "format": ft["format"],
            "url": e["url"],
            "name": e["name"],
            "size": e.get("file_size"),
            "md5": e.get("md5_checksum"),
        })

    tmp = path.with_suffix(".part")
    tmp.write_text(json.dumps(files))
    tmp.replace(path)
    return files


# ----------------------------------------------------------------------
# File helpers
# ----------------------------------------------------------------------

def _local_name(fname):
    """Name of the extracted file (compression suffix removed)."""
    return re.sub(r"\.(gz|Z)$", "", fname)


def _extract(path):
    """
    Decompress a downloaded file and remove the archive.
    .gz -> extract_and_remove_gz ; .Z -> gzip -d ; anything else is kept as is
    (some data centres serve uncompressed Hatanaka .crx files).
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

    return path


def _download(session, entry, output_dir):
    """Download + extract one index entry. Returns the extracted Path, or None."""
    archive = output_dir / entry["name"]
    tmp = Path(f"{archive}.part")
    try:
        print(f"GET   {entry['name']}  ({entry['url'].split('/')[2]})")
        with session.get(entry["url"], stream=True, timeout=(15, 300)) as r:
            r.raise_for_status()
            with open(tmp, "wb") as f:
                for chunk in r.iter_content(1 << 20):
                    f.write(chunk)
        if entry.get("size") and tmp.stat().st_size != int(entry["size"]):
            raise IOError(f"size {tmp.stat().st_size} != {entry['size']}")
        tmp.replace(archive)
        return _extract(archive)
    except Exception as e:
        tmp.unlink(missing_ok=True)
        archive.unlink(missing_ok=True)
        (output_dir / _local_name(entry["name"])).unlink(missing_ok=True)
        print(f"ERROR {entry['name']}: {e}")
        return None


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


def _candidates(files, station):
    """Index entries of `station` (4 or 9 characters) in priority order."""
    station = station.upper()
    sel = [e for e in files
           if (e["marker"] == station if len(station) == 9 else e["marker"][:4] == station)]
    return sorted(sel, key=lambda e: _FORMAT_RANK[e["format"]])


def _best_file(session, entries, output_dir):
    """First entry (highest RINEX version, then next data centre) cached or downloaded."""
    for e in entries:
        local = output_dir / _local_name(e["name"])
        if local.exists():
            print(f"SKIP  {local.name}")
            return local
    for e in entries:
        path = _download(session, e, output_dir)
        if path is not None:
            return path
    return None


def _iter_days(start, end, inclusive=True):
    d = start
    while d <= end if inclusive else d < end:
        yield d
        d += timedelta(days=1)


def _check_station(station):
    station = station.upper()
    if len(station) not in (4, 9):
        raise ValueError("station must be a 4 or 9 character name")
    return station


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------

def stations_available_epos(start, end, full_id=False):
    """
    Union of the stations having daily 30 s RINEX files in EPOS during [start, end).

    full_id : False -> 4-character names (lower case), like stations_available_cors.
              True  -> 9-character markers (e.g. 'AGDE00FRA').
    """
    stations = set()
    with requests.Session() as session:
        for d in _iter_days(start, end, inclusive=False):
            try:
                files = day_index(session, d)
            except Exception as e:
                print(f"ERROR index {d}: {e}")
                continue
            for e in files:
                stations.add(e["marker"] if full_id else e["marker"][:4].lower())
    return sorted(stations)


def list_files_epos(station, start, end):
    """
    Files of one station between start and end (inclusive), in priority order.

    Returns list of (date, url).
    """
    station = _check_station(station)
    found = []
    with requests.Session() as session:
        for d in _iter_days(start, end):
            try:
                files = day_index(session, d)
            except Exception as e:
                print(f"ERROR index {d}: {e}")
                continue
            found += [(d, e["url"]) for e in _candidates(files, station)]
    return found


def download_epos(station, start, end, output_dir="."):
    """
    Download EPOS daily RINEX files of one station, one per day (start and end
    included), highest RINEX version first.

    Returns list of Path (locally available RINEX files).
    """
    station = _check_station(station)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    list_files = []
    with requests.Session() as session:
        for d in _iter_days(start, end):
            try:
                files = day_index(session, d)
            except Exception as e:
                print(f"ERROR index {d}: {e}")
                continue
            entries = _candidates(files, station)
            if not entries:
                print(f"NONE  {station} on {d}")
                continue
            path = _best_file(session, entries, output_dir)
            if path is not None:
                list_files.append(path)

    print()
    print(f"Files available locally: {len(list_files)}")
    return list_files


def download_file(station, date, output_dir="."):
    """
    Download ONE daily RINEX observation file of a station from EPOS
    (4-character name, case-insensitive), RINEX 4/3 preferred over RINEX 2.

    Returns the Path of the extracted file, or None.
    """
    station = station.upper()
    if len(station) != 4:
        raise ValueError("station must be a 4-character station name")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with requests.Session() as session:
        try:
            files = day_index(session, date)
        except Exception as e:
            print(f"ERROR index {date}: {e}")
            return None
        entries = _candidates(files, station)
        if not entries:
            print(f"NONE  {station} on {date}")
            return None
        return _best_file(session, entries, output_dir)
