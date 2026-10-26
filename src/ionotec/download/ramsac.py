from pathlib import Path
from datetime import timedelta

import ssl
import json
import time
import datetime as _dt
import subprocess
import requests

# RAMSAC (IGN Argentina) daily RINEX 2 Hatanaka files, Unix-compressed (.Z):
#   https://ramsac.ign.gob.ar/api/v1/rinex/download/5/fosa0010.24d.Z
# The server does not send its intermediate certificate (Let's Encrypt YR1),
# so requests needs a CA bundle that includes it. It is built automatically
# on the first TLS error (build_ca_bundle). By hand it would be:
#   curl -sS -o yr1.der http://yr1.i.lencr.org/
#   openssl x509 -inform DER -in yr1.der -out yr1.pem
#   cat /etc/ssl/certs/ca-certificates.crt yr1.pem > ~/.cache/ionotec/ramsac-bundle.pem
BASE_URL = "https://ramsac.ign.gob.ar/api/v1/rinex/download"
CA_BUNDLE = Path.home() / ".cache" / "ionotec" / "ramsac-bundle.pem"

# Selector that appears in the download URL (5 works for the daily files;
# 1 gives a 500 error and unknown values give 403). Tried in this order.
IDS = (5,)

# Status codes that mean "this file does not exist"
_NONE_CODES = (403, 404, 500)

# Station codes of the RAMSAC status page (lower case)
STATIONS = [s.lower() for s in """
25MA 3ARO ABRA AGGO ALTA ALUM AUTF AZUL BATE BCAR BCH1 BOLS BORC BSEN CAEP CAQI
CATA CFAG CGSJ CHAC CHIM CHLT CHMA CHOY CORD COYQ CSJ1 CSLO CSOM DINO DORE EBY1
ECAR EPGZ EPSF ESQU EZZA FEDE FMAT FOSA FUTF GGUA GRCA GRLS GUAY GVIL IGM1 ITAI
JBAL JCBC JCHL JVCH JVGO LARJ LHCL LMHS LNQM LPGS LROS MBIO MECO MGLO MGNO MGUE
MPL2 MZAC MZAE MZAL MZAR MZAU MZGA MZRF NEQN NESA NOYA NUGA OAFA OBRA OJOA ORAN
PACO PATA PBCA PDE3 PEBA PEJ1 PELU PIAG PIBV PORT PRCV PRNA PUMA PWIL PWRO PY01
RAF1 RCRO RDCM RECO RIO2 RIO4 RODE RSAL RSAU RSCL RUFI RWSN SAIS SAJA SANT SARM
SBAL SICO SJSF SL01 SMAN SMDM SPRZ SRLP SUAR SUR1 SVIC TATU TAVA TERO TGTA TILC
TOLH TOSF TRN1 TUC1 TUCU UCOR UNPA UNRO UNSA UNSJ UYBU UYCO UYPA UYSA UYSO UYTD
VBCA VCON VCOP VDMA VIMA VIME VIRA VLAN YCBA YEMA ZPLA
""".split()]

# Last file published by each station that is no longer current, from the
# RAMSAC status page (snapshot of 2026-10-22). A station cannot have a file
# after its last one, so later days are not even requested. Stations that
# are current when the snapshot was taken are not listed (never skipped).
# Set USE_LAST_FILE_SKIP = False to ignore this table.
USE_LAST_FILE_SKIP = True
LAST_FILE = {k: _dt.date.fromisoformat(v) for k, v in {
    "alta": "2025-10-09", "borc": "2025-03-18", "coyq": "2020-01-26",
    "csom": "2026-09-08", "eby1": "2026-06-28", "ecar": "2025-09-11",
    "grls": "2026-08-03", "meco": "2024-01-22", "prcv": "2024-09-09",
    "rdcm": "2026-07-18", "rsal": "2025-11-29", "sprz": "2026-02-10",
    "vcop": "2025-09-17", "vira": "2023-02-25", "vlan": "2025-06-09",
    "ycba": "2025-12-12",
}.items()}

# Results of the existence checks are kept on disk, so a station-day is
# never requested twice (None -> no disk cache). "Found" is kept forever;
# "not found" only for days older than NEGATIVE_CACHE_AFTER_DAYS, because
# recent files may still be published.
CACHE_FILE = Path.home() / ".cache" / "ionotec" / "ramsac_exists.json"
NEGATIVE_CACHE_AFTER_DAYS = 7
_CACHE = None


_BUNDLE_TRIED = False


def build_ca_bundle(path=None, host="ramsac.ign.gob.ar"):
    """
    Build a CA bundle that works for RAMSAC: the usual trusted roots
    (certifi) plus the intermediate certificate(s) the server forgets to
    send. They are found by following the "CA Issuers" links of the
    server's certificate (Let's Encrypt YR1 -> its issuer -> ...).
    Needs the `cryptography` package (installed with paramiko).
    """
    import certifi
    from cryptography import x509
    from cryptography.hazmat.primitives import serialization
    from cryptography.x509.oid import AuthorityInformationAccessOID, ExtensionOID

    path = Path(path or CA_BUNDLE)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Only reads the server's own certificate (no verification needed).
    leaf_pem = ssl.get_server_certificate((host, 443))
    cert = x509.load_pem_x509_certificate(leaf_pem.encode())

    extra = []

    for _ in range(5):

        if cert.subject == cert.issuer:          # self-signed root reached
            break

        try:
            aia = cert.extensions.get_extension_for_oid(
                ExtensionOID.AUTHORITY_INFORMATION_ACCESS
            ).value
        except x509.ExtensionNotFound:
            break

        urls = [
            d.access_location.value for d in aia
            if d.access_method == AuthorityInformationAccessOID.CA_ISSUERS
        ]

        if not urls:
            break

        data = requests.get(urls[0], timeout=30).content

        try:
            cert = x509.load_der_x509_certificate(data)
        except ValueError:
            cert = x509.load_pem_x509_certificate(data)

        extra.append(cert.public_bytes(serialization.Encoding.PEM))

    path.write_bytes(
        Path(certifi.where()).read_bytes() + b"\n" + b"\n".join(extra)
    )

    return path


def _ensure_bundle(session):
    """Build the CA bundle once per process; True if the caller should retry."""
    global _BUNDLE_TRIED

    if _BUNDLE_TRIED or CA_BUNDLE is None:
        return False

    _BUNDLE_TRIED = True

    try:
        print("INFO  building the CA bundle for ramsac.ign.gob.ar ...")
        build_ca_bundle()

    except Exception as e:
        print(f"WARNING could not build the CA bundle: {e}")
        return False

    session.verify = str(CA_BUNDLE)
    return True


def _new_session():
    s = requests.Session()
    s.verify = str(CA_BUNDLE) if CA_BUNDLE is not None and Path(CA_BUNDLE).exists() else True
    return s


def _filename(station, d):
    """ssssDDD0.YYd.Z  (daily file, day of year, 2-digit year)"""
    return f"{station}{d.timetuple().tm_yday:03d}0.{d.year % 100:02d}d.Z"


def _call(session, method, url, max_retries=4, **kwargs):
    """Request with retries on connection problems. Returns the Response."""
    for attempt in range(max_retries):
        try:
            return session.request(method, url, timeout=(10, 120), **kwargs)

        except requests.exceptions.SSLError as e:
            if _ensure_bundle(session):
                continue                      # retry with the new bundle

            raise RuntimeError(
                "TLS certificate error: ramsac.ign.gob.ar does not send its "
                "intermediate certificate. Build the CA bundle described at "
                f"the top of ramsac.py (expected at {CA_BUNDLE}). {e}"
            )

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


def _load_cache():
    global _CACHE

    if _CACHE is None:
        _CACHE = {}

        if CACHE_FILE is not None and Path(CACHE_FILE).exists():
            try:
                _CACHE = json.loads(Path(CACHE_FILE).read_text())
            except Exception:
                _CACHE = {}

    return _CACHE


def _save_cache():
    if CACHE_FILE is None or _CACHE is None:
        return

    path = Path(CACHE_FILE)
    path.parent.mkdir(parents=True, exist_ok=True)

    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(_CACHE))
    tmp.replace(path)


def _remember(station, d, found):
    """Store an existence result (negative ones only for old days)."""
    day = d.date() if hasattr(d, "date") else d

    if not found:
        limit = _dt.date.today() - _dt.timedelta(days=NEGATIVE_CACHE_AFTER_DAYS)

        if day > limit:
            return

    _load_cache()[f"{station}:{day:%Y%j}"] = bool(found)


def _exists_cached(session, station, d):
    """_exists() with the disk cache; returns True / False."""
    day = d.date() if hasattr(d, "date") else d
    key = f"{station}:{day:%Y%j}"

    cache = _load_cache()

    if key in cache:
        return cache[key]

    found = _exists(session, station, d) is not None
    _remember(station, d, found)

    return found


def _after_last_file(station, d):
    day = d.date() if hasattr(d, "date") else d
    last = LAST_FILE.get(station)

    return USE_LAST_FILE_SKIP and last is not None and day > last


def _exists(session, station, d):
    """The selector ID of the first URL that has the day's file, or None.
    Uses HEAD (no download); falls back to a streamed GET if HEAD is refused."""
    fname = _filename(station, d)

    for sid in IDS:
        url = f"{BASE_URL}/{sid}/{fname}"

        r = _call(session, "HEAD", url, allow_redirects=True)

        if r.status_code in (405, 501):
            r = _call(session, "GET", url, stream=True)
            r.close()

        if r.status_code == 200:
            return sid

        if r.status_code not in _NONE_CODES:
            r.raise_for_status()

    return None


def _uncompress(z_path):
    """Extract a Unix-compressed .Z file (gzip -dc also reads .Z); remove
    the .Z afterwards. Returns the extracted Path."""
    z_path = Path(z_path)
    out_path = z_path.with_suffix("")

    for cmd in (["gzip", "-dc", str(z_path)], ["uncompress", "-c", str(z_path)]):
        try:
            with open(out_path, "wb") as f:
                subprocess.run(cmd, stdout=f, check=True, stderr=subprocess.PIPE)

            z_path.unlink()
            return out_path

        except FileNotFoundError:
            continue

        except subprocess.CalledProcessError as e:
            out_path.unlink(missing_ok=True)
            raise RuntimeError(f"{cmd[0]} failed on {z_path.name}: {e.stderr!r}")

    out_path.unlink(missing_ok=True)
    raise RuntimeError("neither gzip nor uncompress was found in the PATH")


def stations_available_ramsac(start, end, stations=None):
    """
    Return the union of RAMSAC stations having a daily RINEX file during
    the requested period (4-character lower-case names, like the other
    sources).

    There is no listing endpoint that I know of, so every candidate station
    is checked for every day with a HEAD request (no download), except:
      - stations already found are not checked again on later days;
      - days after a station's last published file (LAST_FILE) are skipped;
      - results are cached on disk (CACHE_FILE), so a station-day is never
        requested twice, even across runs.
    Pass `stations=[...]` to restrict the candidates.

    Archive:
        https://ramsac.ign.gob.ar/api/v1/rinex/download/5/ssssDDD0.YYd.Z
    """

    candidates = [s.lower() for s in (stations or STATIONS)]

    list_stations = set()

    n_new = 0

    try:
        with _new_session() as session:

            d = start

            while d <= end:

                for st in candidates:

                    if st in list_stations or _after_last_file(st, d):
                        continue

                    day = d.date() if hasattr(d, "date") else d
                    was_cached = f"{st}:{day:%Y%j}" in _load_cache()

                    try:
                        if _exists_cached(session, st, d):
                            list_stations.add(st)

                    except RuntimeError:
                        raise

                    except Exception as e:
                        print(f"ERROR {st} {d}: {e}")
                        continue

                    if not was_cached:
                        n_new += 1

                        if n_new % 100 == 0:
                            _save_cache()           # survive interruptions

                d += timedelta(days=1)

    finally:
        _save_cache()

    return sorted(list_stations)


def download_ramsac(station, start, end, output_dir="."):
    """
    Download RAMSAC (IGN Argentina) daily RINEX 2 Hatanaka files.

    Archive structure:
        https://ramsac.ign.gob.ar/api/v1/rinex/download/5/ssssDDD0.YYd.Z

    station: 4-character name (case-insensitive; 9-character ids are cut
    to their first 4 characters).

    The downloaded .Z files are extracted with gzip (or uncompress) and
    removed; the files are kept as Hatanaka RINEX 2 (.YYd).

    Returns
    -------
    list of Path
        Locally available RINEX files.
    """

    station = station.lower()

    if len(station) == 9:
        station = station[:4]

    if len(station) != 4:
        raise ValueError("station must be a 4-character station name")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    list_files = []

    with _new_session() as session:

        d = start

        while d <= end:

            doy = d.timetuple().tm_yday
            fname = _filename(station, d)

            z_path = output_dir / fname
            rinex_path = z_path.with_suffix("")

            # Already extracted
            if rinex_path.exists():
                print(f"SKIP  {rinex_path.name}")
                list_files.append(rinex_path)
                d += timedelta(days=1)
                continue

            try:
                done = False

                for sid in IDS:

                    print(f"GET   {fname}")

                    r = _call(
                        session, "GET", f"{BASE_URL}/{sid}/{fname}", stream=True
                    )

                    with r:
                        if r.status_code in _NONE_CODES:
                            continue

                        r.raise_for_status()

                        with open(z_path, "wb") as f:
                            for chunk in r.iter_content(chunk_size=1024 * 1024):
                                if chunk:
                                    f.write(chunk)

                    list_files.append(_uncompress(z_path))
                    _remember(station, d, True)
                    done = True
                    break

                if not done:
                    _remember(station, d, False)
                    print(f"NONE  {station} on {d} (DOY {doy:03d})")

            except RuntimeError:
                raise

            except Exception as e:

                z_path.unlink(missing_ok=True)
                rinex_path.unlink(missing_ok=True)

                print(f"ERROR {fname}: {e}")

            d += timedelta(days=1)

    _save_cache()

    print()
    print(f"Files available locally: {len(list_files)}")

    return list_files
