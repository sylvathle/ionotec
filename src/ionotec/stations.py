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

import georinex as gr
import pandas as pd
import numpy as np

from os import listdir, getenv, path, makedirs
from os.path import isfile, join

import sys

#import os
#root_dir = os.environ.get('PYTEC_PATH')

from pathlib import Path
# create output dir
base_dir = Path(sys.argv[0]).resolve().parent
target_dir = base_dir / "output/"
target_dir.mkdir(parents=True, exist_ok=True)
root_dir = str(target_dir)+"/"
print ("root_dir", root_dir)

def create_output(folder_name: str) -> Path:
    """
    Create a folder in the same directory as the script that started execution
    (e.g. main.py).
    """
    main_file = Path(sys.argv[0]).resolve()
    base_dir = main_file.parent

    target_dir = base_dir / folder_name
    target_dir.mkdir(parents=True, exist_ok=True)

    return target_dir

#def create_output(base_dir: Path, folder_name: str) -> Path:
#    base_dir = Path(base_dir).resolve()
#    target = base_dir / folder_name
#    target.mkdir(parents=True, exist_ok=True)
#    return target

csv_stations = root_dir+"stations.csv"

def resume_station(year,force=False):
    '''  Function not coordinated by now with only code
    Creates stations.csv that contains information about the navigation
    stations
    '''
    
    print ("LIST stations")
    
    if path.exists(csv_stations) and not force: return
    year_folder = root_dir + str(year)+ "/"
    for f in listdir(year_folder):
        doy = f.replace(".zip","")
        day_folder = year_folder+ doy + "/"
        print (day_folder)

        files = [f for f in listdir(day_folder) if isfile(join(day_folder, f))]
        d = {"station":[],"X":[],"Y":[],"Z":[],"interval":[]}
        for f in files:
            print (f)
            if f[-1]!="o": continue
            station = f[:4]
            if station in d["station"]: continue
            if ((station+str(doy)+"0."+str(year)[2:]+"o" not in files) and
                (station+str(doy)+"1."+str(year)[2:]+"o" not in files)): continue

            rinex = station + str(doy) + "0." + str(year)[2:]

            try: header = gr.rinexheader(day_folder+rinex+"o")
            except ValueError:
                print ("Error in file",f)
            if "INTERVAL" in header.keys(): interval = float(header["INTERVAL"].replace(" ",""))
            else: interval = 1.0
            pos_antena = header['position']
            d["station"].append(station)
            d["X"].append(pos_antena[0])
            d["Y"].append(pos_antena[1])
            d["Z"].append(pos_antena[2])
            d["interval"].append(interval)
            #d["br"].append(float("nan"))
    df = pd.DataFrame(d)
    df.sort_values(by="station",inplace=True)
    df.to_csv(csv_stations,index=False)

def get_closest_stations(p):
    ''' Input: list of three float corresponding to the (X,Y,Z) coordinates
        Output: a list of strings corresponding to the station ordered by the closest to farthest '''
    df = pd.read_csv(csv_stations)
    df.set_index("station",inplace=True)
    #p = df[["X","Y","Z"]].iloc[1]
    #print (p)
    dfdist = df.assign(distance = lambda x: ((x["X"]-p[0])**2+(x["Y"]-p[1])**2+(x["Z"]-p[2])**2))
    dfdist.sort_values(by="distance",inplace=True)
    return dfdist.index.values

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

#class BaseClass:
#    def base_method(self) -> str:
#        """
#        Base method.
#        """
#        return "hello from BaseClass"

#    def __call__(self) -> str:
#        return self.base_method()


#def base_function() -> str:
#    """
#    Base function.
#    """
#    return "hello from base function"
