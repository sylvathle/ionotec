import os
import gzip
import shutil
import requests
from bs4 import BeautifulSoup
import os, sys, time

#from requests.auth import HTTPDigestAuth
 

 
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
                


def get_rinex_from_cddis(year,doy,suff,DEST_DIR_BASE,list_stations=None,nfirst=-1):
    
    BASE_URL = "https://cddis.nasa.gov/archive/gnss/data/daily/"+str(year)+"/"+str(doy)+"/"+suff+"/"
    #DEST_DIR = "/home/sylvain/Documents/jupyter_project/TEC/receiver_dcb/rinex/IGS/"+str(year)+"/"+str(doy)+"/"+suff+"/"
    DEST_DIR = DEST_DIR_BASE+str(year)+"/"+str(doy)+"/"+suff+"/"
    #DEST_DIR = "./downloads/TEC/RINEX/IGS/"+str(year)+"/"+str(doy)+"/"+suff+"/"

    os.makedirs(DEST_DIR, exist_ok=True)

    session = requests.Session()

    max_retry_time = 60  # seconds
    start_time = time.time()
    
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
    
            #print(f"Concat failed ({e}). Retrying in 1 second...")
            time.sleep(1)
        

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
                if station in href.lower(): links.append(href.strip())
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
        print ("No file for ",station,suff)
        return []
    #headers = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.10; rv:39.0)'}
    if nfirst==-1: nfirst =len(links)

    list_rinex_downloaded = []
    
    ilink = 0
    # Step 3: download each file
    for file in links:

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
        




def get_dcb_from_cddis(year,doy,DEST_DIR_BASE):

    BASE_URL = "https://cddis.nasa.gov/archive/gnss/products/bias/"+str(year)+"/"
    DEST_DIR = DEST_DIR_BASE + str(year) + "/"
    os.makedirs(DEST_DIR, exist_ok=True)
    session = requests.Session()

    if os.path.exists(DEST_DIR+"CAS0MGXRAP_"+str(year)+str(doy)+"0000_01D_01D_DCB.BSX"): 
        return DEST_DIR+"CAS0MGXRAP_"+str(year)+str(doy)+"0000_01D_01D_DCB.BSX"
    if os.path.exists(DEST_DIR+"CAS0OPSRAP_"+str(year)+str(doy)+"0000_01D_01D_DCB.BIA"): 
        return DEST_DIR+"CAS0OPSRAP_"+str(year)+str(doy)+"0000_01D_01D_DCB.BIA"
    
    #print (DEST_DIR)

    max_retry_time = 60  # seconds
    start_time = time.time()

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
    
            print(f"Concat failed ({e}). Retrying in 1 second...")
            time.sleep(1)

    # Step 2: extract file links
    links = []
    for a in soup.find_all("a"):
        href = a.get("href")

        if not href:
            continue
    
        if ("CAS0MGXRAP_"+str(year)+str(doy) not in href) and ("CAS0OPSRAP_"+str(year)+str(doy) not in href):
            continue
        print (href)
    
        splthref = href.split('.')
        suff = splthref[0]
        ext = splthref[1]

        file_in_dest = suff+'.'+ext
        if os.path.exists(DEST_DIR+file_in_dest): 
            #print (file_in_dest,"exists")
            continue
        links.append(href.strip())
    
    #print(f"Found {len(links)} files")
    #print (links)
    if len(links)==0: 
        return
    #headers = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.10; rv:39.0)'}

    list_dcb_downloaded = []

    # Step 3: download each file
    for file in links:

        url = BASE_URL + file

        local_path = os.path.join(DEST_DIR, file)
        
        # Skip if the compressed file already exists
        if os.path.isfile(local_path.replace('.gz','')):
            #print(f"Skipping {file} (already downloaded)")
            return local_path.replace('.gz','')
            #continue

        print(f"Downloading {file}...")
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


        
        extract_gz_files(local_path,DEST_DIR+file.replace('.gz',''))
        os.remove(local_path)
        return local_path.replace('.gz','')

