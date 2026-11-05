import os
import io
import gzip
import shutil
import requests
from bs4 import BeautifulSoup
import os, sys, time
import numpy as np
import pandas as pd
from datetime import datetime, timedelta

#from requests.auth import HTTPDigestAuth
 
IGS_SNX_URL = "https://files.igs.org/pub/station/general/igs.snx"

 
def extract_gz_files(source_file,dest_file):
    # Create destination folder if it doesn't exist
    if source_file.endswith(".gz"):
        #gz_path = os.path.join(source_folder, filename)

        # Output file name (remove .gz)
        #output_filename = os.path.splitext(filename)[0]
        #output_path = os.path.join(destination_folder, output_filename)

        #print(f"Extracting {source_file} -> {dest_file}")

        # Decompress
        with gzip.open(source_file, 'rb') as f_in:
            with open(dest_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)


def get_page(session, url, max_retries=5):
    for attempt in range(max_retries):
        try:
            print(f"GET {url} (attempt {attempt + 1}/{max_retries})")

            with session.get(
                url,
                timeout=(10, 60),
                headers={"Connection": "close"},
            ) as r:
                r.raise_for_status()
                content = r.content

            return BeautifulSoup(content, "html.parser")

        except (
            requests.exceptions.ConnectionError,
            requests.exceptions.ChunkedEncodingError,
            requests.exceptions.Timeout,
        ) as exc:
            if attempt == max_retries - 1:
                raise

            wait = min(5 * (2 ** attempt), 60)
            print(f"Request failed: {exc}")
            print(f"Retrying in {wait} seconds...")
            time.sleep(wait)
                


def get_rinex_from_cddis(year,doy,suff,DEST_DIR_BASE,list_stations=None,nfirst=-1):
    
    print ('getRinex',suff,list_stations)
    str_doy = str(doy)
    if doy<10: str_doy = '00'+str_doy
    elif doy<100: str_doy = '0'+str_doy
    BASE_URL = "https://cddis.nasa.gov/archive/gnss/data/daily/"+str(year)+"/"+str_doy+"/"+suff+"/"
    #DEST_DIR = "/home/sylvain/Documents/jupyter_project/TEC/receiver_dcb/rinex/IGS/"+str(year)+"/"+str(doy)+"/"+suff+"/"
    DEST_DIR = DEST_DIR_BASE+str(year)+"/"+str_doy+"/"+suff+"/"
    #DEST_DIR = "./downloads/TEC/RINEX/IGS/"+str(year)+"/"+str(doy)+"/"+suff+"/"

    os.makedirs(DEST_DIR, exist_ok=True)

    session = requests.Session()

    max_retry_time = 60  # seconds
    start_time = time.time()
    
    print (BASE_URL)
    while True:
        #try:
        if True:
            with session.get(BASE_URL, stream=True, timeout=30) as r:
                print ('get BASE URL')
                print("trust_env:", session.trust_env)
                print("netrc authentication:", session.auth)

                settings = session.merge_environment_settings(
                    BASE_URL, {}, None, None, None
                )
                print("trust_env:", session.trust_env)
                print("Session auth:", session.auth)
                print("Environment settings:", settings)
                r.raise_for_status()
                #r.raise_for_status()
                #soup = BeautifulSoup(r.text, "html.parser")
                soup = get_page(session, BASE_URL)
                print (soup)
            break
    
        #except Exception as e:
        #    elapsed = time.time() - start_time
    #
    #        if elapsed >= max_retry_time:
    #            print(f"Failed to contact {BASE_URL} after {elapsed:.1f} s: {e}")
        #        return
    
            #print(f"Concat failed ({e}). Retrying in 1 second...")
    #        time.sleep(1)
        

    # Step 2: extract file links
    links = []
    for a in soup.find_all("a"):
        href = a.get("href")

        if not href:
            continue
    
        if (".rnx" not in href and ".crx" not in href):
            continue
    
        #print (station,href)
        #print (href.strip())
        if list_stations!=None:
            for station in list_stations:
                if station==href.lower()[:4]: links.append(href.strip())
        else: links.append(href.strip())

   

        splthref = href.split('.')
        suff = splthref[0]
        ext = splthref[1]

        file_in_dest = suff+'.'+ext
        if os.path.exists(DEST_DIR+file_in_dest): 
            #print (file_in_dest,"exists")
            continue
            
    #print(f"Found {len(links)} files")
    if len(links)==0: 
        print ("No file for ",station, doy, year)
        return []
    #headers = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.10; rv:39.0)'}
    if nfirst==-1: nfirst =len(links)

    list_rinex_downloaded = []
    
    ilink = 0
    # Step 3: download each file
    for file in links:
        print (file)

        url = BASE_URL + file

        local_path = os.path.join(DEST_DIR, file)
        
        # Skip if the compressed file already exists
        if os.path.isfile(local_path.replace('.gz','')):
            #print(f"Skipping {file} (already downloaded)")
            list_rinex_downloaded.append(local_path.replace('.gz',''))
            continue

        #print(f"Downloading {file}...")
        max_retry_time = 60  # seconds
        start_time = time.time()
        
        while True:
            try:
                with session.get(url, stream=True, timeout=30) as r:
                    r.raise_for_status()
        
                    with open(local_path, "wb") as f:
                        for chunk in r.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
        
                print(f"Successfully downloaded {file}")
                break
        
            except Exception as e:
                elapsed = time.time() - start_time
        
                if elapsed >= max_retry_time:
                    print(f"Failed to download {file} after {elapsed:.1f} s: {e}")
                    raise
        
                print(f"Download failed ({e}). Retrying in 1 second...")
                time.sleep(1)


        extract_gz_files(local_path,local_path.replace('.gz',''))
        os.remove(local_path)
        list_rinex_downloaded.append(local_path.replace('.gz',''))
        ilink+=1
        if ilink>=nfirst: 
            return list_rinex_downloaded
            break
    return list_rinex_downloaded
        

def down_obs_files(datemin,datemax,station,dest):
    
    interval = timedelta(days=1)
    d = datemin
    downloaded_files = []
    while d<datemax:
        year = d.year
        doy = d.timetuple().tm_yday
        yy = d.strftime("%y")
        list_download = get_rinex_from_cddis(year,doy,str(yy)+'d',dest,list_stations=[station],nfirst=-1)
        downloaded_files += list_download
        d += interval

    return downloaded_files


def get_dcb_from_cddis(year, doy, DEST_DIR_BASE):
    BASE_URL = "https://cddis.nasa.gov/archive/gnss/products/bias/" + str(year) + "/"
    DEST_DIR = DEST_DIR_BASE + str(year) + "/"
    os.makedirs(DEST_DIR, exist_ok=True)
    session = requests.Session()

    str_doy = str(doy)
    if doy < 10:
        str_doy = "00" + str_doy
    elif doy < 100:
        str_doy = "0" + str_doy

    # Define the target patterns in strict priority order
    target_patterns = [
        f"CAS0MGXRAP_{year}{str_doy}0000_01D_01D_DCB.BSX",
        f"CAS0OPSRAP_{year}{str_doy}0000_01D_01D_DCB.BIA",
        f"CAS1OPSRAP_{year}{str_doy}0000_01D_01D_DCB.BIA"
    ]

    # Step 1: Check if any uncompressed file already exists locally in priority order
    for pattern in target_patterns:
        local_path = os.path.join(DEST_DIR, pattern)
        if os.path.exists(local_path):
            return local_path

    # Step 2: Fetch the directory page from CDDIS if not found locally
    max_retry_time = 600  # seconds
    start_time = time.time()
    soup = None

    while True:
        try:
            with session.get(BASE_URL, stream=True, timeout=30) as r:
                r.raise_for_status()
                soup = BeautifulSoup(r.text, "html.parser")
            break

        except Exception as e:
            elapsed = time.time() - start_time
            if elapsed >= max_retry_time:
                print(f"Failed to contact {BASE_URL} after {elapsed:.1f} s: {e}")
                return
            print(f"Connection failed ({e}). Retrying in 1 second...")
            time.sleep(1)

    # Collect all available links and strip whitespace to prevent trailing spaces/%20 issues
    available_links = {a.get("href").strip() for a in soup.find_all("a") if a.get("href")}

    # Step 3: Check patterns on the server in priority order and download only the first match found
    for pattern in target_patterns:
        matching_link = None
        for href in available_links:
            if pattern in href:
                matching_link = href
                break

        if not matching_link:
            continue  # Try the next lower-priority pattern

        local_file_path = os.path.join(DEST_DIR, matching_link)
        uncompressed_path = local_file_path.replace(".gz", "")

        # If the uncompressed file already exists from a previous run
        if os.path.exists(uncompressed_path):
            return uncompressed_path

        # Download the highest priority file found
        print(f"Downloading {matching_link}...")
        download_retry_time = 60
        download_start = time.time()

        while True:
            try:
                with session.get(BASE_URL + matching_link, stream=True, timeout=30) as r:
                    r.raise_for_status()
                    with open(local_file_path, "wb") as f:
                        for chunk in r.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)

                print(f"Successfully downloaded {matching_link}")
                break

            except Exception as e:
                elapsed = time.time() - download_start
                if elapsed >= download_retry_time:
                    print(f"Failed to download {matching_link} after {elapsed:.1f} s: {e}")
                    raise
                print(f"Download failed ({e}). Retrying in 1 second...")
                time.sleep(1)

        # Extract if compressed, clean up, and return immediately
        if matching_link.endswith(".gz"):
            extract_gz_files(local_file_path, uncompressed_path)
            os.remove(local_file_path)

        return uncompressed_path

    # If none of the patterns exist locally or on the server
    print(f"No valid DCB file found for year {year}, DOY {str_doy}")
    return None



def read_igs_snx(filename=None, url=IGS_SNX_URL):
    """
    Read IGS SINEX station site information from the +SITE/ID block.

    Parameters
    ----------
    filename : str or None
        Local SINEX file. If None, download the current IGS igs.snx.
    url : str
        URL used when filename is None.

    Returns
    -------
    pandas.DataFrame
        Columns:
            station     : 4-character station code
            point       : point code
            domes       : DOMES number
            description : Station description
            longitude   : Geocentric longitude [decimal degrees]
            latitude    : Geocentric latitude [decimal degrees]
            height      : Height [m]
    """

    # Helper function to convert Degrees, Minutes, Seconds to Decimal Degrees
    def dms_to_deg(d_str, m_str, s_str):
        try:
            deg = float(d_str)
            min = float(m_str)
            sec = float(s_str)
            sign = -1 if deg < 0 else 1
            return sign * (abs(deg) + min/60 + sec/3600)
        except ValueError:
            return np.nan # Use NaN for invalid conversions

    # ------------------------------------------------------------------
    # Open the file either locally or remotely
    # ------------------------------------------------------------------
    if filename is None:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        file = io.StringIO(response.text)
    else:
        file = open(filename, "r")

    records = []

    try:
        in_site_id = False

        for line_num, line in enumerate(file):

            # Start of SITE/ID block
            if line.startswith("+SITE/ID"):
                in_site_id = True
                continue

            # End of SITE/ID block
            if line.startswith("-SITE/ID"):
                in_site_id = False
                continue

            if not in_site_id:
                continue

            # Skip comments/header within the block (lines starting with '*')
            if line.startswith("*") or not line.strip():
                continue

            # Example line from igs.snx:
            # abmf  A 97103M001 P Les Abymes,Guadeloupe  298 28 20.9  16 15 44.3   -25.7
            # Based on IGS SINEX (1-based) column specs and example:
            # CODE PT __DOMES__ T _STATION DESCRIPTION__ _LONGITUDE_ _LATITUDE__ HEIGHT_
            # 1-4  6  8-16      18 20-40                  42-45 (D)   55-58 (D)   68-73
            #                                            47-48 (M)   60-61 (M)
            #                                            50-53 (S)   63-66 (S)
            try:
                station_code = line[1:5].strip()
                point_code = line[8:9].strip()
                domes = line[9:18].strip()
                # line[17:18] is 'T' (Type), which we skip for now
                description = line[21:43].strip()

                # Parse longitude DMS (1-based: 42-45, 47-48, 50-53)
                lon_d_str = line[44:47].strip()
                lon_m_str = line[48:50].strip()
                lon_s_str = line[51:55].strip()

                # Parse latitude DMS (1-based: 55-58, 60-61, 63-66)
                lat_d_str = line[56:59].strip()
                lat_m_str = line[60:62].strip()
                lat_s_str = line[63:67].strip()

                # Parse height (1-based: 68-73)
                height_str = line[68:75].strip()

                longitude = dms_to_deg(lon_d_str, lon_m_str, lon_s_str)
                latitude = dms_to_deg(lat_d_str, lat_m_str, lat_s_str)
                height = float(height_str) if height_str else np.nan

                # Normalize longitude to -180 to 180
                if not np.isnan(longitude):
                    longitude = (longitude + 180) % 360 - 180

            except (ValueError, IndexError):
                # If parsing fails for a line, skip it
                # For debugging: print(f"Skipping line {line_num+1} due to parsing error: {line.strip()}")
                continue

            records.append({
                "station": station_code,
                "point": point_code,
                "domes": domes,
                "description": description,
                "longitude": longitude,
                "latitude": latitude,
                "height": height,
            })

    finally:
        file.close()

    # Return DataFrame directly, no pivoting needed for this block
    df = pd.DataFrame(records)

    return df
