#############################################################################3
#
# MIT License
#
# Copyright (c) 2021 Sylvain Blunier
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  SOFTWARE.
#
#############################################################################3

#import georinex as gr
import pandas as pd
import numpy as np

import pymap3d as pm

from os import listdir, getenv, path, makedirs
from os.path import isfile, join

import sys

from . import rinex

#import os
#root_dir = os.environ.get('PYTEC_PATH')

from pathlib import Path
# create output dir
base_dir = Path(sys.argv[0]).resolve().parent
target_dir = base_dir / "output/"
#target_dir.mkdir(parents=True, exist_ok=True)
#print ("Output directory is :",target_dir)
root_dir = str(target_dir)+"/"

#print ("root_dir", root_dir)

#def create_output(folder_name: str) -> Path:
#    """
#    Create a folder in the same directory as the script that started execution
###    (e.g. main.py).
#    """
#    main_file = Path(sys.argv[0]).resolve()
#    base_dir = main_file.parent

#    target_dir = base_dir / folder_name
#    target_dir.mkdir(parents=True, exist_ok=True)

#   return target_dir

#def create_output(base_dir: Path, folder_name: str) -> Path:
#    base_dir = Path(base_dir).resolve()
#    target = base_dir / folder_name
#    target.mkdir(parents=True, exist_ok=True)
#    return target

csv_stations = root_dir+"stations.csv"

def resume_station(folder):
    '''  Function not coordinated by now with only code
    Creates stations.csv that contains information about the navigation
    stations
    '''

    search_dir = Path(folder)
    files_o = [f for f in search_dir.rglob("*o") if f.is_file()]
    files_d = [f for f in search_dir.rglob("*d") if f.is_file()]
    files_crx = [f for f in search_dir.rglob("*crx") if f.is_file()]
    files_rnx = [f for f in search_dir.rglob("*rnx") if f.is_file()]

    files = files_o+files_d+files_crx+files_rnx
    
    d = {"station":[],"X":[],"Y":[],"Z":[],"resolution(s)":[]}
    for f in files:
        fname = f.name
        f_path = f.resolve()
        #print (f.replace(folder,""))
        #continue
        #if fname[-1]!="o": continue


        #try: header = gr.rinexheader(f_path)
        try: 
            rx = rinex.rinex(f_path) 
            header = rx.read_header()
        except ValueError:
            print ("Error in file",f_path)
        if "INTERVAL" in header.keys(): interval = float(header["INTERVAL"].replace(" ",""))
        else: interval = 1.0
        
        if header['type']!='O': continue

        station = header['name_station']
        #print (fname,station)
        if station in d["station"]: continue
        
        pos_antena = header['position']
        d["station"].append(station)
        d["X"].append(pos_antena[0])
        d["Y"].append(pos_antena[1])
        d["Z"].append(pos_antena[2])
        d["resolution(s)"].append(interval)
        #d["br"].append(float("nan"))
    

    df = pd.DataFrame(d)
    df["lat"],df["lon"],df["alt"] = pm.ecef2geodetic(df["X"],df["Y"],df["Z"])
    df.sort_values(by="station",inplace=True)
    try:
        df.to_csv(csv_stations,index=False)
        return df.set_index("station")
    except:
        return False

#def resume_station(year_folder,force=False):
#    '''  Function not coordinated by now with only code
#    Creates stations.csv that contains information about the navigation
#    stations
#    '''

def get_closest_stations(p):
    ''' Input: list of three float corresponding to the (X,Y,Z) coordinates
        Output: a list of strings corresponding to the station ordered by the closest to farthest '''
    df = pd.read_csv(csv_stations)
    df.set_index("station",inplace=True)

    dfdist = df.assign(distance = lambda x: np.sqrt((x["X"]-p[0])**2+(x["Y"]-p[1])**2+(x["Z"]-p[2])**2))
    dfdist.sort_values(by="distance",inplace=True)
    #print (dfdist.head())
    return dfdist.index.values.tolist()

def get_station_pos(station):
    '''	Input: string, name of a station
    Output: [X,Y,Z] np.ndarray corresponding to its position in the ECEF reference system '''

    df = pd.read_csv(csv_stations).set_index("station")
    pos = df[["X","Y","Z"]].loc[station].values
    return pos

def get_station_interval(station):

    df = pd.read_csv(csv_stations).set_index("station")
    pos = float(df[["interval"]].loc[station].values)
    return pos
