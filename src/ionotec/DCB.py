"""
Download P1P2YYMM.DCB files from the CODE FTP server at:
    ftp://ftp.aiub.unibe.ch/CODE/YYYY/P1P2YYMM.DCB
 
Uses Python's built-in ftplib so no third-party packages are needed.
 
Usage
-----
    python download_p1p2_dcb.py                        # current month
    python download_p1p2_dcb.py 2024-03                # specific month
    python download_p1p2_dcb.py 2023-01 2023-06        # date range (inclusive)
 
The file for a given date is determined by its year+month (YYMM).
"""
 
import sys
import os
import argparse
from datetime import datetime, date, timedelta
import time
import re
from calendar import monthrange
from ftplib import FTP, error_perm, all_errors
from pathlib import Path
#import gzip
#import shutil

import scipy.constants as csts
import math

#import requests
#from bs4 import BeautifulSoup

import pandas as pd
from . import decompress
from . import stations as st
from . import freq

from . import igs
 
# ── Constants ────────────────────────────────────────────────────────────────
 
FTP_HOST        = "ftp.aiub.unibe.ch"
FTP_TIMEOUT     = 60   # seconds per connection
RETRY_ATTEMPTS  = 3    # how many times to try each file
RETRY_DELAY     = 5    # seconds to wait between retries
INTER_FILE_DELAY = 2   # seconds to wait between successive files
#DCB_dir = "C:\\Users\\User\\Documents\\GitHub\\pytec\\tests\\DCB\\"



 
# ── Helpers ───────────────────────────────────────────────────────────────────
 
def _dcb_path(
    dt: datetime,
) -> tuple[str, str]:
    """
    Return (remote_path, filename) for the P1P2 DCB file covering *dt*.
 
    Remote path is relative to the FTP root, e.g. 'CODE/2024/P1P22403.DCB'.
    """
    yy       = dt.strftime("%y")    # two-digit year
    mm       = dt.strftime("%m")    # zero-padded month
    yyyy     = dt.strftime("%Y")    # four-digit year (directory)
    filename = f"P1P2{yy}{mm}.DCB.Z"
    remote   = f"CODE/{yyyy}/{filename}"
    return remote, filename



def _dcb_local_path(
    dt: datetime,
) -> tuple[str, str]:
    """
    Return (remote_path, filename) for the P1P2 DCB file covering *dt*.
 
    Remote path is relative to the FTP root, e.g. 'CODE/2024/P1P22403.DCB'.
    """
    yy       = dt.strftime("%y")    # two-digit year
    mm       = dt.strftime("%m")    # zero-padded month
    yyyy     = dt.strftime("%Y")    # four-digit year (directory)
    filename = f"P1P2{yy}{mm}.DCB"
    return os.path.join(DCB_dir,filename)

    
# ── Core download ─────────────────────────────────────────────────────────────
 
def _cleanup(path: str) -> None:
    """Remove a file if it exists and is empty (leftover from a failed transfer)."""
    if os.path.exists(path) and os.path.getsize(path) == 0:
        os.remove(path)
 
 
def _ftp_download_file(
    remote: str,
    filename: str,
    dest_dir: str,
) -> str:
    """
    Download each (remote_path, filename) pair using a **fresh FTP connection
    per file** with retries and a small pause between files.
 
    Returns a list of successfully saved local paths.
    """
    os.makedirs(dest_dir, exist_ok=True)
    saved: str = ""
 
    #for idx, (remote, filename) in enumerate(remote_paths):
    if True:
        local_path = os.path.join(dest_dir, filename)
 
        # Pause between files to avoid hammering the server
        #if idx > 0:
        #    time.sleep(INTER_FILE_DELAY)
 
        print(f" {filename} … ", flush=True)
 
        success = False
        for attempt in range(1, RETRY_ATTEMPTS + 1):
            try: 
                #Open a *fresh* FTP connection, download one file, then close.
                #Raises on any error so the caller can retry.
                with FTP(FTP_HOST, timeout=FTP_TIMEOUT) as ftp:
                    ftp.login()           # anonymous
                    ftp.set_pasv(True)    # passive mode — avoids NAT/firewall issues on Windows
                    with open(local_path, "wb") as fh:
                        ftp.retrbinary(f"RETR /{remote}", fh.write)
                        
                size = os.path.getsize(local_path)
                print(f"  OK  ({size:,} bytes)  →  {local_path}")
                saved = local_path
                success = True
                break
 
            except error_perm as exc:
                # 550 = file not found / permission denied — no point retrying
                print(f"  SKIPPED  (server: {exc})")
                _cleanup(local_path)
                break
 
            except Exception as exc:
                _cleanup(local_path)
                if attempt < RETRY_ATTEMPTS:
                    print(f"  attempt {attempt} failed ({exc}); retrying in {RETRY_DELAY}s …", flush=True)
                    time.sleep(RETRY_DELAY)
                else:
                    print(f"  FAILED after {RETRY_ATTEMPTS} attempts: {exc}")
 
    return saved
    
 
# ── Public API ────────────────────────────────────────────────────────────────
 
def download_dcb(
    dt: datetime,
    dest_dir: str = "."
) -> str:
    """
    Download the P1P2 DCB file for the month that contains *dt*.
 
    Parameters
    ----------
    dt       : any ``datetime.date`` (or ``datetime.datetime``) object
    dest_dir : directory where the file will be saved (created if needed)
 
    Returns
    -------
    Local path of the downloaded file.
    """
    #if isinstance(dt, datetime):
    #    dt = dt.date()
    url, filename = _dcb_path(dt)
    dcb_file = dest_dir+filename.replace(".Z","")
    if os.path.exists(dcb_file): return dcb_file 
    #dcb_Z = download_file(url, filename, dest_dir)
    dcb_Z = _ftp_download_file(url, filename, dest_dir)

    # Decompress the file
    with open(dcb_Z, "rb") as f:
        decompressed = decompress.decompress_z(f.read())   # returns bytes
    output_path = dcb_Z.replace('.Z','')
    with open(output_path, "wb") as f:
            f.write(decompressed)
    os.remove(dcb_Z)
    return output_path

    

    
    
def getSatDCB(date,target_sv):
    dcb_path = _dcb_local_path(date)
    if not os.path.exists(dcb_path):
        download_dcb(date, dest_dir=DCB_dir)
        print ("Downloaded path", dcb_path)
        
    """
    Parses a standard GNSS P1-P2 DCB file and returns the float bias value 
    for a specific satellite vehicle (e.g., 'G17').
    
    Returns None if the satellite is not found.
    """
    if not os.path.exists(dcb_path):
        raise FileNotFoundError(f"The file {dcb_path} does not exist.")
        
    target_sv = target_sv.strip().upper()
    
    with open(dcb_path, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip empty lines and header comments
            if not line or line.startswith('%') or line.startswith('#'):
                continue
                
            parts = line.split()
            
            # Check if the row starts with a valid satellite PRN identifier (e.g., G01, R12)
            if re.match(r'^[GREC]\d{2}$', parts[0]):
                current_sv = parts[0].upper()
                
                # If this row matches the satellite of interest
                if current_sv == target_sv:
                    try:
                        # Extract the numeric bias value (usually the last element)
                        return float(parts[-2])
                    except ValueError:
                        continue
                        
    # Return None if the satellite isn't listed in the file
    return None

# Regex that matches a data line, e.g.:
#   "G01                          -6.791       0.010"
#   "R26                          -0.819       0.048"
_DATA_LINE = re.compile(
    r"^([A-Z]\d{2})\s+"          # PRN  e.g. G01, R07, E11
    r"(-?\d+\.\d+)\s+"           # VALUE (ns), possibly negative
    r"(-?\d+\.\d+)"              # RMS   (ns)
)
 
# Filename pattern  P1P2YYMM.DCB  (case-insensitive)
_FNAME_PAT = re.compile(r"P1P2(\d{2})(\d{2})\.DCB", re.IGNORECASE)
 
 
 
def read_p1p2_dcb(date) -> pd.DataFrame:
    """
    Parse a CODE P1P2YYMM.DCB file into a tidy DataFrame.
 
    Parameters
    ----------
    path : path-like
        Path to the .DCB file.
 
    Returns
    -------
    pd.DataFrame with columns:
        sv    (str)              – satellite PRN, e.g. 'G01'
        time  (datetime, UTC)   – first day of the solution month
        DCB   (float)           – P1-P2 differential code bias  [ns]
        STD   (float)           – RMS uncertainty of the DCB    [ns]
    """

    dcb_path = _dcb_local_path(date)
    if not os.path.exists(dcb_path):
        download_dcb(date, dest_dir=DCB_dir)
        print ("Downloaded path", dcb_path)
        
    path = Path(dcb_path)
    epoch = date
 
    rows: list[dict] = []
    with path.open("r", encoding="ascii", errors="replace") as fh:
        for line in fh:
            m = _DATA_LINE.match(line)
            if m:
                rows.append(
                    {
                        "sv":   m.group(1),
                        "time": epoch,
                        "DCB":  float(m.group(2)),
                        "STD":  float(m.group(3)),
                    }
                )
 
    if not rows:
        raise ValueError(f"No data rows found in '{path}'.")
 
    df = pd.DataFrame(rows, columns=["sv", "time", "DCB", "STD"])
    return df

def read_dcb(dcb_path) -> pd.DataFrame:
    """
    Parse a CODE P1P2YYMM.DCB file into a tidy DataFrame.
 
    Parameters
    ----------
    path : path-like
        Path to the .DCB file.
 
    Returns
    -------
    pd.DataFrame with columns:
        sv    (str)              – satellite PRN, e.g. 'G01'
        time  (datetime, UTC)   – first day of the solution month
        DCB   (float)           – P1-P2 differential code bias  [ns]
        STD   (float)           – RMS uncertainty of the DCB    [ns]
    """
        
    path = Path(dcb_path)

    #print (path.name)
    epoch = datetime.strptime(path.name, "P1P2%y%m.DCB").date()
    #print (epoch)
    #sys.exit()
    #epoch = date
 
    rows: list[dict] = []
    with path.open("r", encoding="ascii", errors="replace") as fh:
        for line in fh:
            m = _DATA_LINE.match(line)
            if m:
                rows.append(
                    {
                        "sv":   m.group(1),
                        "time": epoch,
                        "C1": "P1",
                        "C2": "P2",
                        "dcb":  float(m.group(2)),
                        "std":  float(m.group(3)),
                    }
                )
 
    if not rows:
        raise ValueError(f"No data rows found in '{path}'.")
 
    df = pd.DataFrame(rows, columns=["sv", "time", "C1","C2", "dcb", "std"])
    return df


def getBias_fromfile(sat,file,chanel1, chanel2):
    fbias = open(file,'r')
    #if sat[0]=="G":
    for line in fbias.readlines():
        splt_line = line.split()
        if splt_line[0] == "DSB":
            if splt_line[2]==sat:
                if splt_line[3]==chanel1:
                    if splt_line[4]==chanel2:
                        #print (sat,chanel1,chanel2,splt_line[	8])
                        return float(splt_line[8])
    return float('NaN')

def load_dcb(datemin=None,datemax=None):

    DCB_dir = st.root_dir + 'DCB/'
    os.makedirs(DCB_dir,exist_ok=True)

    list_used_dcb = []
    
    if (datemin is None) or (datemax is None): return pd.DataFrame()
        
    directory_DCB_path = Path(DCB_dir)
    list_sorted_DCB = []
    for file in directory_DCB_path.rglob("*"): 
        if file.is_file():
            list_sorted_DCB.append(str(file.resolve()))
            
    for i in range((datemax.date() - datemin.date()).days + 1):
        d = datemin.date() + timedelta(days=i)
        year = d.year
        doy = (d - date(year,1,1)).days + 1
        DCB_file = DCB_dir +str(year)+ "/CAS0MGXRAP_"+str(year)+str(doy)+"0000_01D_01D_DCB.BSX"
        DCB_file = igs.get_dcb_from_cddis(year,doy,DCB_dir)
        list_used_dcb.append(DCB_file)
    
    df_dcb = pd.DataFrame()
    
    for f in list_used_dcb:
        dict_dcb = {"time":[],"sv":[],"C1":[],"C2":[],"dcb":[],"std":[]}
        if f[-3:]=="BIA" or f[-3:]=="BSX":
            #print (f)
            fsplit = f.split("_")
            strt = fsplit[-4]
            #year = int(strt[:4])
            #doy = int(strt[4:7])
            date_obj = datetime.strptime(strt[:7], "%Y%j").date()
            fbias = open(f,'r')
            #if sat[0]=="G":
            for line in fbias.readlines():
                splt_line = line.split()
                if splt_line[0] == "DSB":
                    sv = splt_line[2]
                    if len(sv)!=3: continue    
                    ns_to_tecu = freq.getAlpha(sv,splt_line[3],splt_line[4]) * csts.c * 1e-9  / 1e16
                    
                    if math.isnan(ns_to_tecu): continue
                    #ns_to_tecu = 1
                    dict_dcb["time"].append(date_obj)
                    dict_dcb["sv"].append(splt_line[2])
                    dict_dcb["C1"].append(splt_line[3])
                    dict_dcb["C2"].append(splt_line[4])

                    
                    dict_dcb["dcb"].append(float(splt_line[-2])*ns_to_tecu)
                    dict_dcb["std"].append(float(splt_line[-1])*ns_to_tecu)

            df_dcb = pd.concat([df_dcb,pd.DataFrame(dict_dcb)])
        else:
            df = read_dcb(f)
            df_dcb = pd.concat([df_dcb,df])
            #print (df)
            #print ("Need to implement load of ",f)
    
    df_dcb['time'] = pd.to_datetime(df_dcb['time'])
    #df_dcb.set_index('time',inplace=True)
    #print (df_dcb)
    #sys.exit()
    return df_dcb
