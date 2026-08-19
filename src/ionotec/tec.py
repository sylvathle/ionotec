#############################################################################
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
#############################################################################
##import os,sys
#import psutil

#import gnss
import georinex as gr

#import stations as st
import pandas as pd
import scipy.constants as csts
import datetime
import numpy as np
import sys,os,subprocess
import matplotlib.pyplot as plt
#import matplotlib.dates as md
import pymap3d as pm
#from os import listdir,path,mkdir
#from os.path import isfile, join
#import julian
import math	
import time
import re


import psutil
#import gc
#from pympler import asizeof

from pathlib import Path

#from .stations import root_dir
from . import stations as st
from . import gnss
from . import rinex as rx
from . import freq
from . import DCB


import warnings

from pathlib import Path



pd.options.mode.chained_assignment = None

# Earth Radius
R_E = 6371000





def tleft(border):
    return border["t_left"]


'''
def filter_slope_leap(list_borders,series,diffs):
    list_borders.sort(key=tleft)
    #return list_borders
    for i in range(len(list_borders)-1):
        if list_borders[i+1]["left"]<=list_borders[i]["right"]:
            list_borders[i+1]["left"]=list_borders[i]["right"]+1
            list_borders[i+1]["t_left"]=series.index[list_borders[i]["right"]+1]
    i=0
    while i<len(list_borders):
        #short segment
        t_short_segment = abs((list_borders[i]["t_right"]-list_borders[i]["t_left"]).seconds)<3*60
        short_segment = list_borders[i]["right"]-list_borders[i]["left"]<6
        negative_segment = list_borders[i]["t_right"]<list_borders[i]["t_left"]
        if t_short_segment or negative_segment or short_segment:
            list_borders.pop(i)
        else: i+=1
    i=1
    while i<len(list_borders)-1:
        if (list_borders[i]["t_right"]-list_borders[i]["t_left"]).seconds>30*60:
            i+=1
            continue
        A,B,M,m=fit_lin(diffs[list_borders[i]["left"]:list_borders[i]["right"]],series[list_borders[i]["left"]:list_borders[i]["right"]].values)
        if not (list_borders[i-1]["right_A"] is None or list_borders[i+1]["left_A"] is None):
            if abs(A-list_borders[i-1]["right_A"])>0.02 and abs(list_borders[i+1]["left_A"]-A)>0.02:
                list_borders.pop(i)
            else: i+=1
        else:
            i+=1
    return list_borders
'''

def filter_slope_leap(list_borders, series, diffs):
    # --- Step 1: order segments chronologically ---
    list_borders.sort(key=tleft)

    # --- Step 2: resolve overlaps between consecutive segments ---
    # If segment i+1 starts before (or exactly where) segment i ends,
    # push its left edge to just after segment i's right edge.
    for i in range(len(list_borders) - 1):
        if list_borders[i + 1]["left"] <= list_borders[i]["right"]:
            list_borders[i + 1]["left"] = list_borders[i]["right"] + 1
            list_borders[i + 1]["t_left"] = series.index[list_borders[i]["right"] + 1]


    # --- Step 3: drop degenerate / too-short segments ---
    i = 0
    while i < len(list_borders):
        # duration under 3 minutes
        t_short_segment = abs((list_borders[i]["t_right"] - list_borders[i]["t_left"]).seconds) < 3 * 60
        # fewer than 6 samples
        short_segment = list_borders[i]["right"] - list_borders[i]["left"] < 6
        # t_right earlier than t_left (invalid/negative duration, can happen after Step 2's shift)
        negative_segment = list_borders[i]["t_right"] < list_borders[i]["t_left"]

        if t_short_segment or negative_segment or short_segment:
            list_borders.pop(i)          # drop it, don't advance i (next element shifts into this slot)
        else:
            i += 1

    # --- Step 4: drop segments whose slope doesn't stand out from both neighbors ---
    # Skips the first and last segment (i starts at 1, stops before len-1),
    # since those have no left/right neighbor respectively to compare against.
    i = 1
    while i < len(list_borders) - 1:
        # only re-examine segments short enough to plausibly be a spurious leap (<30 min)
        if (list_borders[i]["t_right"] - list_borders[i]["t_left"]).seconds > 30 * 60:
            i += 1
            continue

        # fit a line to diffs vs. series values over this segment
        # A = slope, B = intercept, M/m = max/min (unused here, but returned by fit_lin)
        A, B, M, m = fit_lin(
            diffs[list_borders[i]["left"]:list_borders[i]["right"]],
            series[list_borders[i]["left"]:list_borders[i]["right"]].values,
        )

        # only compare if both neighbors already have a computed slope
        #if not (list_borders[i - 1]["right_A"] is None or list_borders[i + 1]["left_A"] is None):
        #    # if this segment's slope differs meaningfully from BOTH neighbors,
        #    # it's likely a spurious leap -> drop it
        #    if abs(A - list_borders[i - 1]["right_A"]) > 0.02 and abs(list_borders[i + 1]["left_A"] - A) > 0.02:
        #        print ('Second pop', list_borders[i - 1]["right_A"] , list_borders[i + 1]["left_A"]  , A - list_borders[i - 1]["right_A"],list_borders[i + 1]["left_A"] - A)
        #        list_borders.pop(i)
        #    else:
        #        i += 1
        #else:
        #    i += 1
        i += 1
    return list_borders


'''def plot_leap(diffs,series,s,A,B,N,borders,title=""):
    fig, ax = plt.subplots(1,figsize=(10,7))
    xx = np.array(diffs[s:s+N])
    yy = A*xx+B
    if s<len(diffs):
        xright=diffs[s+N]
        tx_right = series.index[s+N]
        yright=A*xright+B
        ax.plot([tx_right],yright,'rx')
        arrow_yrightmargin = (series[s+N]-yright)*0.012
        x0 = md.date2num(series.index[s+N])
        ax.arrow(x0,yright+arrow_yrightmargin,x0-x0,series[s+N]-yright-2*arrow_yrightmargin,width=0.00003,head_length=0.07*(series[s+N]-yright-2*arrow_yrightmargin),length_includes_head=True,color="black")
    if s>0:
        xleft=diffs[s-1]
        tx_left = series.index[s-1]
        yleft=A*xleft+B
        ax.plot([tx_left],yleft,'rx')
        arrow_yleftmargin = (series[s-1]-yleft)*0.012
        x1 = md.date2num(series.index[s-1])
        ax.arrow(x1,yleft+arrow_yleftmargin,x1-x1,series[s-1]-yleft-2*arrow_yleftmargin,width=0.00003,head_length=0.07*(series[s-1]-yleft-2*arrow_yleftmargin),length_includes_head=True,color="black")

    txx = [series.index[s+i] for i in range(N)]

    ax.plot(txx,yy,'k')
    plt.plot(series.iloc[max(0,s-4):min(s+N+6,len(series))],'x')
    ax.plot(series,'bx')
    for b in borders:
        if b["left"]!=None:
            ax.axvline(x=series.index[b["left"]],color='r')
    for b in borders:
        if b["right"]!=None:
            ax.axvline(x=series.index[b["right"]],color='b')

    plt.title(title)
    mng = plt.get_current_fig_manager()
    mng.full_screen_toggle()

    plt.show()
    plt.close()
'''

def process_br(df,h):
    # Coefficients of the quadratic error function that will be computed
    a,b=0,0

    # Compute cos\chi for full serie
    df["ci"] = np.cos(np.arcsin(R_E*np.cos(df["elevation"])/(R_E+h)))
    df.dropna(subset=["STEC_l","elevation"],inplace=True)

    if len(df)<2: return float('NaN')

    # Create list of time of the station removing duplicates
    time_series = df.groupby(level="time").size()
    time_series = df.index


    # Compute sums
    for t in time_series:
    #for t in df.index:
        sum_ci = 0 # sum of cos\chi_i for a coefficient (to be squared)
        sum_ci2 = 0 # sum of squared of cos\chi_i for a and b coefficient
        sum_sici = 0 # sum of STEC_i*cos\chi_i for b coefficient
        sum_sici2 = 0 # sum of STEC_i*cos\chi_i^2 for b coefficient
        N=0 # Number of satellites with data at time t
        # Get subserie containing data at time t
        df_t = df.loc[t]
        # If only one satellite has data, df_t is a Series, not a Dataframe, nothing to compute
        if isinstance(df_t,pd.Series): continue
        # Compute sums over each satellite
        for index,row in df_t.iterrows():
            si = row["STEC_l"]
            ci = row["ci"]
            sum_ci += ci
            sum_ci2 += ci**2
            sum_sici += si*ci
            sum_sici2 += si*ci*ci
            N+=1
        if N==0: continue
        # Update a and b over time serie, we ignore the main 1/N factor since it cancels for bias
        # We also ignore "2" factor for a since is cancels in -b/(2a) root calculation
        a=a+(sum_ci2-(1/N)*sum_ci**2)/N
        b=b+(sum_sici2-(1/N)*sum_sici*sum_ci)/N
        if math.isnan(a) or math.isnan(b): return float('NaN')

    if a==0: return float("nan")
    # root of error function = receiver bias.
    br = b/a
    return br




class tec:

    ## Definition of the constants for the equations of STEC
    gps_f1, gps_f2, gps_f5 = 1575.42 * 1e6, 1227.60 * 1e6, 1176.45e6
    gps_lambda1, gps_lambda2, gps_lambda5 = csts.c/gps_f1, csts.c/gps_f2, csts.c/gps_f5
    gps_alpha = gps_f1**2*gps_f2**2/(gps_f1**2-gps_f2**2)/40.318

    channels = {}

    is_gnss_processed = False

    dcb_extensions = ['BSX','BIA','DCB']
    nav_extension = ['n', 'g', 'h', 'rnx']
    obs_extension = ['crx', 'd', 'o']

    folder_debug = 'debug/'

    dict_f_obs = {}
    list_f_nav = []
    list_f_dcb = []


    source_data_folder = ""

    def __init__(self):
        self.station = ""
        self.coord = ""


        #Resolution that will be used for the process (should by 60 seconds)
        self.resolution = 60

        # Resolution of the rinex station, in seconds.
        #self.obs_resolution = ""

        # DataFrame containing observation data
        self.df_obs = pd.DataFrame()

        # List of satellites that have been observed
        self.sv = []

        # List of dictionnaries containing borders of signal for each satellite
        self.borders={}

        ## MAIN OUTPUT
        self.list_df = {}

        self.list_obs_stations = []
        
        self.df_sat_DCB = pd.DataFrame()
        
        # Time and year of file
        self.year = None
        self.doy = None
        self.t_min = {}
        self.t_max = {}

        # Bias of the antenna, calculated by method "compute_reveiver_bias"
        # Value stored in csv "stations.csv"
        self.br = {'station':[],'constellation':[],'time_i':[],'time_f':[],'br':[]}
        
        if not os.path.exists(st.root_dir + "TEC/"):
            os.mkdir(st.root_dir + "TEC/")
        
    def load_rinex_folder(self,folder,date_min=None, date_max=None,h=350000):
        self.h = h

        self.datemin = date_min
        self.datemax = date_max
        self.source_data_folder = rinex_folder
        self.prepare_files()

        self.gnss = gnss.gnss(self.list_f_nav)

        self.sat_dcb = DCB.load_dcb(self.list_f_dcb,datemin=self.datemin,datemax=self.datemax)
        #

    def load_rinex_lists(self,list_f_obs,list_f_nav=None,list_f_dcb=None,datemin=None, datemax=None, h=350000):
        self.h = h

        self.list_f_obs = list_f_obs 
        
        # List of navigation files
        self.list_f_nav = list_f_nav if not list_f_nav is None else []
        
        # File containing satellite bias
        self.list_f_dcb = list_f_dcb if not list_f_dcb is None else []
        
        self.list_df = {}
                
        
        for f_obs in list_f_obs:
            rfile = rx.rinex(f_obs)
            header = rfile.read_header()
            ## Avoid running files that are not in the requested date range
            datelim_cond = (
		(datemin is None or header["t_first_obs"].date() >= datemin.date()) and
		(datemax is None or header["t_first_obs"].date() <= datemax.date())
            )
            if not datelim_cond: continue
            
            self.station = header["name_station"].replace(" ","").lower()
            self.coord = header['position']
            list_df_f_obs = rfile.read_data()

            for const, df in list_df_f_obs.items():
                df.set_index('time',inplace=True)
                if const in self.list_df.keys():
                    self.list_df[const] = pd.concat([self.list_df[const],df])
                else:
                    self.list_df[const] = df
                #print (self.list_df[const])
            
        self.datemin = None
         
        for constellation in self.list_df.keys():
            if self.datemin is None: self.datemin = min(self.list_df[constellation].index)
            else: self.datemin = min(self.datemin, min(self.list_df[constellation].index))
        if not datemin is None: self.datemin = max(self.datemin,datemin)
        #else: self.datemin = datemin

        self.datemax= None
        #if datemax is None: 
        for constellation in self.list_df.keys():
            if self.datemax is None: self.datemax = max(self.list_df[constellation].index)
            else: self.datemax = max(self.datemax, max(self.list_df[constellation].index))
        if not datemax is None: self.datemax = min(self.datemax,datemax)
        #else: self.datemax = datemax      
        

        const_to_del = []
        for k in self.list_df.keys():
            if k!="S": continue
            const_to_del.append(k)
        for const in const_to_del:
            del self.list_df[const]
            
        list_all_sv = []
        for const in self.list_df.keys():
            list_all_sv += self.list_df[const]['sv'].unique().tolist()            

            
            
        self.gnss = gnss.gnss(datemin=self.datemin,datemax=self.datemax,list_satellites=list_all_sv)
        self.gnss.load_all_sats()        

        self.sat_dcb = DCB.load_dcb(self.list_f_dcb,datemin=self.datemin,datemax=self.datemax)
        #print (self.sat_dcb)



    def prepare_files(self):
        directory_path = Path(self.source_data_folder)
        
        # List all files recursively
        files = sorted([file for file in directory_path.rglob("*") if file.is_file()])
        list_day_f_obs = {}
        for file in files:
            ext = file.name.split(".")[-1]
            end_ext = ''.join([char for char in ext if not char.isdigit()])

            if ext=="crx":
                year = int(file.name[-26:-22])
                doy = int(file.name[-22:-19])
            elif ext=="rnx":
                year = int(file.name[-22:-18])
                doy = int(file.name[-18:-15])
            elif end_ext in ["d","o","n","g","h"]:
                year = int(file.name[-3:-1])+2000
                doy = int(file.name[-8:-5])
        
            if end_ext in self.obs_extension:
                date = datetime.datetime(year, 1, 1, 0, 0, 0, 0) + datetime.timedelta(days=doy - 1)

                if self.datemin is not None:
                    if (date < self.datemin): continue
                if self.datemax is not None:
                    if (date > self.datemax): continue              
                station = file.name[:4].lower()
                if station not in self.list_obs_stations:
                    self.list_obs_stations.append(station)
                

                if not station in list_day_f_obs.keys(): 
                    list_day_f_obs[station] = {date:[file]}
                    #self.dict_f_obs[station] = [file]
                else: 
                    if not date in list_day_f_obs[station].keys():
                        list_day_f_obs[station][date] = [file]
                    else:
                        list_day_f_obs[station][date].append(file)# = {"day":date,"files":[]}
                        
                    
            if end_ext in self.nav_extension:
                date = datetime.datetime(year, 1, 1, 0, 0, 0, 0) + datetime.timedelta(days=doy - 1)
                if self.datemin is not None:
                    if (date < self.datemin): continue
                if self.datemax is not None:
                    if (date > self.datemax): continue                
                self.list_f_nav.append(file)

            if end_ext in self.dcb_extensions:
                if (end_ext=='BSX') or (end_ext=='BIA'):
                    year = int(file.name[-27:-23])
                    doy = int(file.name[-23:-20])
                    date = datetime.datetime(year, 1, 1, 0, 0, 0, 0) + datetime.timedelta(days=doy - 1)
                    if self.datemin is not None:
                        if (date < self.datemin): continue
                    if self.datemax is not None:
                        if (date > self.datemax): continue             
                elif (end_ext=='DCB'):
                    year = int(file.name[-8:-6])+2000
                    month = int(file.name[-6:-4])
                    if self.datemin is not None:
                        if year<self.datemin.year: continue
                        if year==self.datemin.year and month<self.datemin.month: continue
                    if self.datemax is not None:
                        if year>self.datemax.year: continue
                        if year==self.datemax.year and month>self.datemax.month: continue                
                self.list_f_dcb.append(str(file))


        for station in list_day_f_obs.keys():
            self.dict_f_obs[station] = []
            for date in list_day_f_obs[station].keys():              
                for file in list_day_f_obs[station][date]:
                    if '.crx' in file.name: 
                        self.dict_f_obs[station].append(str(file))
                        break
                    if file.name[-1]=='d': 
                        self.dict_f_obs[station].append(str(file))
                        break
                    if file.name[-1]=='o': 
                        self.dict_f_obs[station].append(str(file))
                        break

    def get_observation_station_list(self):
        return self.list_obs_stations

    def rinex_to_stec(self,station):
        ''' Extract the relevant data for tec calculation from observation and navigation
            Compute STEC of pseudo range and code phase
        '''                   
        
        ### GPS
        if 'G' in self.list_df.keys():
            const = 'G'
            list_cols = self.list_df[const].columns

            chan={}

            if ("C1C" in list_cols) and ("C2W" in list_cols) and ("L1C" in list_cols) and ("L2W" in list_cols) and ("S1C" in list_cols) and ("S2W" in list_cols):
                C1,C2,L1,L2,S1,S2 = "C1C","C2W","L1C","L2W","S1C","S2W"
                chan = {"C1":C1,"C2":C2,"L1":L1,"L2":L2,"S1":S1,"S2":S2}
            elif ("C1" in list_cols) and ("P2" in list_cols) and ("L1" in list_cols) and ("L2" in list_cols):
                C1,C2,L1,L2,S1,S2 = "C1","P2","L1","L2","S1","S2"
                chan = {"C1":"P1","C2":"P2","L1":"L1","L2":"L2","S1":"S1","S2":"S2"}
            elif ("P1" in list_cols) and ("P2" in list_cols) and ("L1" in list_cols) and ("L2" in list_cols):
                C1,C2,L1,L2,S1,S2 = "P1","P2","L1","L2","S1","S2"
                chan = {"C1":"P1","C2":"P2","L1":"L1","L2":"L2","S1":"S1","S2":"S2"}

            if chan:

                self.channels[const] = []
                self.channels[const].append(chan)
                #self.list_df[const].set_index("time",inplace=True)
                self.t_min[const] = min(self.list_df[const].index)
                self.t_max[const] = max(self.list_df[const].index)  
                
                self.list_df[const].rename(columns={S1:'S1', S2:'S2'},inplace=True)      
                
                self.list_df[const]['STEC_l'] = (self.list_df[const][L1]*self.gps_lambda1-self.list_df[const][L2]*self.gps_lambda2)*self.gps_alpha/1e16
                self.list_df[const]['STEC_p'] = (self.list_df[const][C2]-self.list_df[const][C1])*self.gps_alpha/1e16

                #self.list_df[const] = self.list_df[const][['sv',"STEC_l","STEC_p"]]
                self.list_df[const] = self.list_df[const].dropna(subset=["STEC_l","STEC_p"])
                
                self.list_df[const]["C1"] = chan["C1"]
                self.list_df[const]["C2"] = chan["C2"]
                

            else: del self.list[const]

            
        ### GLONASS
        if 'R' in self.list_df.keys():

            const = 'R'
            list_sv = self.list_df['R']['sv'].unique().tolist()
            if len(list_sv)==0: return   

            list_cols = self.list_df[const].columns
            #print (const,list_cols)
            
            chan = {}
            if ("C1P" in list_cols) and ("C2P" in list_cols) and ("L1P" in list_cols) and ("L2P" in list_cols) and ("S1P" in list_cols) and ("S2P" in list_cols):
                C1,C2,L1,L2,S1,S2 = "C1P","C2P","L1P","L2P","S1P","S2P"
                for c in ['P','C']:
                    varc1 = 'C1'+c
                    if varc1 in list_cols:
                        C1 = varc1
                        S1 = 'S1'+c
                        break
                for c in ['P','C']:
                    varc2 = 'C2'+c
                    if varc2 in list_cols:
                        C2 = varc2
                        S2 = 'S2'+c
                        break
                for c in ['P','C']:
                    varc1 = 'L1'+c
                    if varc1 in list_cols:
                        L1 = varc1
                        break
                for c in ['P','C']:
                    varc2 = 'L2'+c
                    if varc2 in list_cols:
                        L2 = varc2
                        break  
                chan = {"C1":C1,"C2":C2,"L1":L1,"L2":L2,"S1":S1,"S2":S2,"S1":"S1","S2":"S2"}
            elif ("L1" in list_cols) and ("L2" in list_cols):
                L1,L2,S1,S2 = "L1","L2","S1","S2"
                if ("P2" in list_cols):
                    C2="P2"
                    if ("P1" in list_cols): C1="P1"
                    elif ("C1" in list_cols): C1="C1"
                elif ("C2" in list_cols):
                    C2="C2"
                    if ("P1" in list_cols): C1="P1"
                    elif ("C1" in list_cols): C1="C1"
                chan = {"C1":"P1","C2":"P2","L1":L1,"L2":L2,"S1":"S1","S2":"S2"}

            if chan:
                self.channels[const] = []
                self.channels[const].append(chan)
                
                #self.list_df[const].set_index("time",inplace=True)
                self.t_min[const] = min(self.list_df[const].index)
                self.t_max[const] = max(self.list_df[const].index)      
    
                df_glonass = pd.DataFrame()
                for sv in list_sv:
                    df_sat_glo  = self.list_df[const][self.list_df[const]["sv"]==sv]
                    df_sat_glo["f1"] = freq.getFrequency(sv,chan['C1'])
                    df_sat_glo["f2"] = freq.getFrequency(sv,chan['C2'])
                    df_sat_glo["alpha"] = freq.getAlpha(sv,chan['C1'],chan['C2'])
                    df_glonass = pd.concat([df_glonass,df_sat_glo])

                self.list_df[const] = df_glonass.dropna(subset=[C1,C2,L1,L2])
                self.list_df[const]["lambda1"] = csts.c/self.list_df[const]["f1"]
                self.list_df[const]["lambda2"] = csts.c/self.list_df[const]["f2"]
                
                self.list_df[const].rename(columns={S1:'S1', S2:'S2'},inplace=True)
                
                #print (self.list_df[const].columns)
                
                self.list_df[const]["STEC_p"] = (self.list_df[const][C2] - self.list_df[const][C1])*self.list_df[const]["alpha"]/1e16
                self.list_df[const]["STEC_l"] = (self.list_df[const]["lambda1"]*self.list_df[const][L1] - \
                                               self.list_df[const]["lambda2"]*self.list_df[const][L2])*self.list_df[const]["alpha"]/1e16
                #print (self.list_df[const].columns)
                self.list_df[const].dropna(subset=["STEC_p","STEC_l"],inplace=True)
                
                if len(self.list_df[const])!=0:
                    self.list_df[const]["C1"] = chan["C1"]
                    self.list_df[const]["C2"] = chan["C2"]
                    self.list_df[const] = self.list_df[const][['sv',"C1","C2",'S1','S2',"STEC_l","STEC_p"]]
                    self.t_min[const] = min(self.list_df[const].index)
                    self.t_max[const] = max(self.list_df[const].index)
            else: del self.list_df[const]

 
        ## GALILEO
        if "E" in self.list_df.keys():

            const = 'E'
            self.list_df[const] = self.list_df[const].dropna(axis=1, how='all')
            list_cols = self.list_df[const].columns
            
            #print (const,self.list_df[const].columns)

            
            C1, C2 = '', ''
            L1, L2 = '', ''
            S1, S2 = '', ''
            chan = {}
            for c in ['C','X','A','B','Z']:
                varc1 = 'C1'+c
                if varc1 in list_cols:
                    C1 = varc1
                    S1 = 'S1'+c
                    break
            for c in ['Q','X','I']:
                varc2 = 'C5'+c
                if varc2 in list_cols:
                    C2 = varc2
                    S2 = 'S5'+c
                    break
            for c in ['C','X','A','B','Z']:
                varc1 = 'L1'+c
                if varc1 in list_cols:
                    L1 = varc1
                    break
            for c in ['Q','X','I']:
                varc2 = 'L5'+c
                if varc2 in list_cols:
                    L2 = varc2
                    break
            if (C1!='') and (C2!='') and (L1!='') and (L2!='') and (S1!='') and (S2!=''): 
                chan = {"C1":C1,"C2":C2,"L1":L1,"L2":L2,"S1":S1,"S2":S2}

            if chan:
                self.channels[const] = []
                self.channels[const].append(chan)
                
                #self.list_df[const].set_index("time",inplace=True)
                self.list_df[const]["STEC_p"] = (self.list_df[const][C2] - self.list_df[const][C1])*self.gps_alpha/1e16
                self.list_df[const]["STEC_l"] = (self.gps_lambda1*self.list_df[const][L1] - self.gps_lambda5*self.list_df[const][L2])*self.gps_alpha/1e16
                self.list_df[const]["C1"] = chan["C1"]
                self.list_df[const]["C2"] = chan["C2"]
                self.list_df[const].rename(columns={S1:'S1', S2:'S2'},inplace=True)
                self.list_df[const] = self.list_df[const][['sv',"C1","C2",'S1','S2',"STEC_l","STEC_p"]]
                #self.list_df['E'].dropna(inplace=True)
                self.list_df[const] = self.list_df[const].dropna(subset=["STEC_l","STEC_p"])
                
                self.t_min[const] = min(self.list_df[const].index)
                self.t_max[const] = max(self.list_df[const].index)
            else:
                del self.list_df[const]




        #### Beidu
    
        if "C" in self.list_df.keys():
            const = 'C'
            bds_f1 = 1561.098 * 1e6
            bds_f2 = 1207.14 * 1e6
            bds_f3 = 1268.52 * 1e6
            
            bds_lambda1 = csts.c/bds_f1
            bds_lambda2 = csts.c/bds_f2
            bds_lambda3 = csts.c/bds_f3
           
            beidu_columns = self.list_df[const].columns
            beidu_alpha_1 = (bds_f1**2*bds_f2**2/(bds_f1**2-bds_f2**2)/40.318)/1e16
            beidu_alpha_2 = (bds_f1**2*bds_f3**2/(bds_f1**2-bds_f3**2)/40.318)/1e16
            beidu_alpha_3 = (bds_f2**2*bds_f3**2/(bds_f2**2-bds_f3**2)/40.318)/1e16
            list_cols = ['time','sv']
            
            
            self.channels[const] = []
            
            df_beidu_2_6 = pd.DataFrame()
            df_beidu_2_7 = pd.DataFrame()
            if ("L2I" in beidu_columns) and ("L7I" in beidu_columns) and ("C7I" in beidu_columns) and ("C2I" in beidu_columns):
                C1,C2 = 'C2I','C7I'
                L1,L2 = 'L2I','L7I'
                df_beidu_2_7 = pd.DataFrame()
                chan = {"C1":C1,"C2":C2,"L1":L1,"L2":L2}
                
                df_beidu_2_7 = self.list_df[const].copy()
                df_beidu_2_7["C1"] = chan["C1"]
                df_beidu_2_7["C2"] = chan["C2"]
                df_beidu_2_7.rename(columns={'S2I':'S1', 'S7I':'S2'},inplace=True)
                df_beidu_2_7["STEC_l"] = (df_beidu_2_7[L1]*bds_lambda1-df_beidu_2_7[L2]*bds_lambda2)*beidu_alpha_1
                df_beidu_2_7["STEC_p"] = (df_beidu_2_7[C2]-df_beidu_2_7[C1])*beidu_alpha_1
                df_beidu_2_7.dropna(subset="STEC_l",inplace=True)
                if len(df_beidu_2_7)!=0: 
                    self.channels[const].append(chan)
            if ("L2I" in beidu_columns) and ("L6I" in beidu_columns) and ("C6I" in beidu_columns) and ("C2I" in beidu_columns):
                C1,C2 = 'C2I','C6I'
                L1,L2 = 'L2I','L6I'
                df_beidu_2_6 = pd.DataFrame()
                chan = {"C1":C1,"C2":C2,"L1":L1,"L2":L2}
                df_beidu_2_6 = self.list_df[const].copy()
                df_beidu_2_6["C1"] = chan["C1"]
                df_beidu_2_6["C2"] = chan["C2"]
                df_beidu_2_6.rename(columns={'S2I':'S1', 'S6I':'S2'},inplace=True)
                df_beidu_2_6["STEC_l"] = (df_beidu_2_6[L1]*bds_lambda1-df_beidu_2_6[L2]*bds_lambda3)*beidu_alpha_2
                df_beidu_2_6["STEC_p"] = (df_beidu_2_6[C2]-df_beidu_2_6[C1])*beidu_alpha_2
                df_beidu_2_6.dropna(subset="STEC_l",inplace=True)
                if len(df_beidu_2_6)!=0: 
                    self.channels[const].append(chan)

            self.list_df[const] = pd.DataFrame()
            if len(df_beidu_2_6)!=0: self.list_df[const] = pd.concat([self.list_df[const],df_beidu_2_6])
            if len(df_beidu_2_7)!=0: self.list_df[const] = pd.concat([self.list_df[const],df_beidu_2_7])


            if len(self.list_df[const])!=0:
                #self.list_df[const].set_index("time",inplace=True)
                #print (self.list_df[const])
                self.list_df[const] = self.list_df[const][['sv',"C1","C2",'S1','S2',"STEC_l","STEC_p"]]
                self.list_df[const].dropna(subset="STEC_l",inplace=True)         
                self.t_min[const] = min(self.list_df[const].index)
                self.t_max[const] = max(self.list_df[const].index)
            else:
                del self.list_df[const]
        #### QZSS
        if "J" in self.list_df.keys():
            const = 'J'
            
            self.list_df[const] = self.list_df[const].dropna(axis=1, how='all')

            C1, C2, L1, L2, S1, S2 = '', '', '', '', '', ''
            for c in ['C','X','S','L','Z']:
                varc1 = 'C1'+c
                if varc1 in self.list_df[const].columns:
                    C1 = varc1
                    S1 = 'S1'+c
                    break
            for c in ['Q','X','I']:
                varc2 = 'C5'+c
                if varc2 in self.list_df[const].columns:
                    C2 = varc2
                    S2 = 'S2'+c
                    break
            for c in ['C','X','S','L','Z']:
                varc1 = 'L1'+c
                if varc1 in self.list_df[const].columns:
                    L1 = varc1
                    break
            for c in ['Q','X','I']:
                varc2 = 'L5'+c
                if varc2 in self.list_df[const].columns:
                    L2 = varc2
                    break
            if (C1!='') and (C2!='') and (L1!='') and (L2!='') and (S1!='') and (S2!=''): 
                chan = {"C1":C1,"C2":C2,"L1":L1,"L2":L2,"S1":S1,"S2":S2}
                    
            #chan = {"C1":C1,"C2":C2,"L1":L1,"L2":L2}

            if chan:
                self.channels[const] = []
                self.channels[const].append(chan)
                
                qzss_alpha = self.gps_f1**2*self.gps_f5**2/(self.gps_f1**2-self.gps_f5**2)/40.318
                #self.list_df[const].set_index("time",inplace=True)
                self.list_df[const]["STEC_l"] = (self.list_df[const][L1]*self.gps_lambda1-self.list_df[const][L2]*self.gps_lambda5)*qzss_alpha/1e16
                self.list_df[const]["STEC_p"] = (self.list_df[const][C2]-self.list_df[const][C1])*qzss_alpha/1e16
                self.list_df[const]["C1"] = chan["C1"]
                self.list_df[const]["C2"] = chan["C2"]
                self.list_df[const].rename(columns={S1:'S1', S2:'S2'},inplace=True)
                self.list_df[const].dropna(subset=["STEC_l","STEC_p"],inplace=True)
                self.list_df[const] = self.list_df[const][['sv',"C1","C2",'S1','S2',"STEC_l","STEC_p"]]

                self.t_min[const] = min(self.list_df[const].index)
                self.t_max[const] = max(self.list_df[const].index)

            else:
                del self.list_df[const]



        if "S" in self.list_df.keys():

            const = "S"
            list_columns=self.list_df[const]
            C1, C2, L1, L2 = 'C1C', 'C5I', 'L1C', 'L5I'
            chan = {}
            if (C1 in list_columns) \
                    and (C2 in list_columns)\
                    and (L1 in list_columns)\
                    and (L2 in list_columns):

                chan = {"C1":C1,"C2":C2,"L1":L1,"L2":L2}

            if chan:

                self.channels[const] = []
                self.channels[const].append(chan)
                #self.list_df[const].set_index("time",inplace=True)
                self.t_min[const] = min(self.list_df[const].index)
                self.t_max[const] = max(self.list_df[const].index)        

                sbas_alpha = self.gps_f1**2*self.gps_f5**2/(self.gps_f1**2-self.gps_f5**2)/40.318
                self.list_df[const]['STEC_l'] = (self.list_df[const][L1]*self.gps_lambda1-self.list_df[const][L2]*self.gps_lambda5)*sbas_alpha/1e16
                self.list_df[const]['STEC_p'] = (self.list_df[const][C2]-self.list_df[const][C1])*sbas_alpha/1e16
    
                self.list_df[const]["C1"] = chan["C1"]
                self.list_df[const]["C2"] = chan["C2"]
                
                self.list_df[const] = self.list_df[const][['sv',"C1","C2","STEC_l","STEC_p"]]

            else: del self.list_df[const]

        #print ("Just loaded all obs rinex")
        #process = psutil.Process(os.getpid())
        #mem = process.memory_info()
        #print(f"RSS: {mem.rss / 1024**2:.2f} MB")
        #print(f"VMS: {mem.vms / 1024**2:.2f} MB")

        #for const in self.list_df.keys(): print(f"list_df[{const}]: {sys.getsizeof(self.list_df[const])/ 1024**2:.2f} MB")

        # List satellites seen by the station and prepare dict of list_borders
        const_to_del = []
        
        for const in self.list_df.keys():
            if len(self.list_df[const])==0:
                const_to_del.append(const)
                continue

            for s in self.list_df[const]["sv"].values:
                if s not in self.sv:
                    self.sv.append(s)
                    self.borders[s]=[]
                            

        # Remove constellation that have no data
        for const in const_to_del:
            del self.list_df[const]
            del self.channels[const]

        dict_oper = {"S1": "mean", "S2": "mean", "STEC_l": "mean","STEC_p": "mean"}
        for const in self.list_df.keys():
            #self.list_df[const]=self.list_df[const].groupby(["sv","C1","C2"]).resample("1min").agg({"STEC_l": "mean","STEC_p": "mean"}).reset_index().set_index('time').dropna(subset=['STEC_l','STEC_p'])
            #print (const)
            self.list_df[const]=self.list_df[const].groupby(["sv","C1","C2"]).resample("1min").agg(dict_oper).reset_index().set_index('time').dropna(subset=['STEC_l','STEC_p'])
            
            #self.list_df[const] = (
            #        self.list_df[const]
            #            .groupby(["sv","C1","C2"])
            #            .resample("1min")
            #            .mean()              # aggregation for numeric columns
            #            .reset_index()
            #        ).set_index('time')
            #print (self.list_df[const])

        #const_to_del = []
        #for const in self.list_df.keys():
        #    if const=='C': continue
        #    const_to_del.append(const)
            
        #for const in const_to_del:
        #    del self.list_df[const]
        #    del self.channels[const]
 
        return True

    def add_satellite_pos(self):
        const_without_pos = []
        
        for const in self.list_df.keys():

            self.list_df[const] = self.gnss.getElevation(self.list_df[const],self.coord)
            
            #self.list_df[const].dropna(subset="elevation",inplace=True)
            self.list_df[const] = self.list_df[const][self.list_df[const]["elevation"]>0]
            
            if len(self.list_df[const])==0: 
                const_without_pos.append(const)
                continue

            self.list_df[const] = self.gnss.getPiercingPoint(self.list_df[const],self.coord,self.h)            

           # print ("dropna elevation")
           # self.list_df[const].dropna(subset="elevation",inplace=True)
           # process = psutil.Process(os.getpid())
           # mem = process.memory_info()
           # print(f"RSS: {mem.rss / 1024**2:.2f} MB")
           # print(f"VMS: {mem.vms / 1024**2:.2f} MB")
           # print(f"list_df[{const}]: {sys.getsizeof(self.list_df[const])/ 1024**2:.2f} MB")
            #print (self.list_df[const])

            self.list_df[const] = self.list_df[const][["sv","C1","C2","elevation","lat","lon","alt","S1","S2","STEC_l","STEC_p"]]
                        
        for const in const_without_pos:
            del self.list_df[const]
            if const in self.channels.keys(): del self.channels[const]
        return True
        
    def list_leaps_series(self,series,tol_dev=0.2,tol_sig=10,N=5):
    
        series.dropna(inplace=True)
	
        indices = series.index
        # Discard too short series
        if len(series)<=N: return []        

        # List time in seconds from first index of series
        diffs = [0]
        for i in range(len(series)-1):
            diffs.append(diffs[-1]+(series.index[i+1]-series.index[i]).seconds)

        list_series = series.tolist()
        
        tol_deviation=tol_dev*self.resolution

        # List time in seconds from first index of series
        diffs = [0]
        for i in range(len(series)-1):
            diffs.append(diffs[-1]+(series.index[i+1]-series.index[i]).seconds)
        
        s=0
        
        list_left_borders = []
        list_right_borders = []

        border = {'left':None, 'right':None}
        
        while s<len(series):
            A,B,max_dev,mean_dev = fit_lin(diffs[s:s+N],series[s:s+N].values)
            # Compute distance of point s and s+N+1 with fit
            left_dev=float(abs(list_series[s-1]-A*diffs[s-1]-B)) if s>0 else 0
            right_dev=float(abs(list_series[s+N]-A*diffs[s+N]-B)) if s+N<len(series) else 0
            #list_fit_params[s]={"A":float(A),"B":float(B),"max_dev":float(max_dev),"left_dev":float(left_dev),"right_dev":float(right_dev)}
 
            

            
            #Looking for left borders
            # beginning of serie
            if s==0: 
                list_left_borders.append(0)
                border['left'] = indices[0]
            else:
                left_dev=float(abs(list_series[s-1]-A*diffs[s-1]-B))
                if left_dev>tol_deviation: 
                    if len(list_left_borders)>0:
                        t_segment_size = (indices[s]-indices[list_left_borders[-1]]).seconds
                        n_segment_size = s - list_left_borders[-1]
                        if (t_segment_size<=N*self.resolution) or (n_segment_size<=N):
                            list_left_borders[-1]=s
                            border['left']=indices[s]
                            #else: 
                        else: # Big enough to be considered as a segment
                            list_left_borders.append(s)
                            if border['left']==None or border['left']>indices[s-1]: 
                                border['left']=indices[s]
                                
                            else:
                                border['right'] = indices[s-1]
                                print ('completed 1 ---------> ',border)
                                border = {'left':indices[s], 'right':None}
                    else:
                        list_left_borders.append(s)
                        border['left']=indices[s]
                else:
                    if border['left']==None: print ('PROBLEM in leap series, no left border found while running towards right border 11')
                    else: border['right'] = indices[s-1]
            
            # Looking for right borders
            if s==len(series)-1: 
                list_right_borders.append(len(series)-1)
                if border['left']!=None: 
                    border['right']=indices[len(series)-1]
                    print ('completed 2 ---------> ',border)
                # else: last segment is useless because no left border, nothing to do here
                print ('OUT')
                return
            else:
                if s+N<len(series):
                    right_dev=float(abs(list_series[s+N]-A*diffs[s+N]-B))
                    if right_dev>tol_deviation: 
                        if len(list_right_borders)>0:
                            t_segment_size = (indices[s+N-1]-indices[list_right_borders[-1]]).seconds
                            n_segment_size = s+N-1 - list_right_borders[-1]
                            if (t_segment_size<=N*self.resolution) or (n_segment_size<=N):
                                list_right_borders[-1]=list_right_borders[-1]#s+N-1
                                #border['left']=indices[s+N-1]
                                border = {'left':indices[s+N], 'right':None}
                            else: 
                                list_right_borders.append(s+N-1)
                                if border['left']==None:
                                    border = {'left':indices[s+N], 'right':None}
                                else:
                                    border['right'] = indices[s+N-1]
                                    print ('completed 3---------> ',border)
                                    border = {'left':indices[s+N], 'right':None}
                            
                        else: 
                            list_right_borders.append(s+N-1)
                            if border['left']==None:
                                border = {'left':indices[s+N], 'right':None}
                            else:
                                border['right'] = indices[s+N-1]
                                print ('completed 4---------> ',border)
                                border = {'left':indices[s+N], 'right':None}
                    else:
                        if border['left']==None: print ('PROBLEM in leap series, no left border found while running towards right border 22')
                        else: border['right'] = indices[s+N-1]                        
                    
            s+=1
        
        print ('left')
        for l in list_left_borders: print (l,indices[l])
        print ('right')
        for r in list_right_borders: print (r,indices[r])
        
        sys.exit()
        
        border = {'left':None,'right':None}
        s_left = 0
        s_right = 0
        while s_left<len(list_left_borders):
            
            if s_left+1<len(list_left_borders):
                if list_left_borders[s_left+1]<list_right_borders[s_right]:
                    right = list_left_borders[s_left+1]
                else:
                    right = list_right_borders[s_right]
                    s_right = s_right+2                               
            else:
                right = list_right_borders[s_right]
            border = {'left':list_left_borders[s_left] ,'right':right }
            s_left = s_left+1
            
            #while list_left_borders[s_left]<list_right_borders[s_right]:
            #    s_left+=1
                
            
            #s_right+=1    
            print (s_left,s_right,list_left_borders[s_left],list_left_borders[s_right],border,len(list_left_borders))
            #s_left+=1
        
        
        sys.exit()
            
                    


            


    '''
    def list_leaps_series(self,series,tol_dev=0.2,tol_sig=10,N_=5,debug=False):
       

        #Remove nan data in series
        series.dropna(inplace=True)
        indices = series.index
        # Discard too short series
        if len(series)<=N_: return []

        # Final list that will contain all the borders
        list_borders = []
        N=N_

        tol_deviation=tol_dev*self.resolution
        

        # List time in seconds from first index of series
        diffs = [0]
        for i in range(len(series)-1):
            diffs.append(diffs[-1]+(series.index[i+1]-series.index[i]).seconds)

        

        list_fit_params = {}

        unstable_left=True

        list_series = series.tolist()
                
        s=0
        while s<len(list_series)-N:
            A,B,max_dev,mean_dev = fit_lin(diffs[s:s+N],series[s:s+N].values)
            # Compute distance of point s and s+N+1 with fit
            left_dev=float(abs(list_series[s-1]-A*diffs[s-1]-B)) if s>0 else 0
            right_dev=float(abs(list_series[s+N]-A*diffs[s+N]-B)) #if s+N<len(series) else 0
            #list_fit_params[s]={"A":float(A),"B":float(B),"max_dev":float(max_dev),"left_dev":float(left_dev),"right_dev":float(right_dev)}
            list_fit_params[s]={"A":round(float(A),2),"B":round(float(B),2),"max_dev":round(float(max_dev),2),"left_dev":round(float(left_dev),2),"right_dev":round(float(right_dev),2)}
            s+=1
        

        border={"left":None,"t_left":None,"right":None,"t_right":None,"left_A":None,"left_B":None,"right_A":None,"right_B":None}


        retroactive_borders = []

        s1=0
        s2=0
        oldest_left_border = 0
        while s1<len(series):
            # Looking for an right border
            #print (s1,s2,indices[s1].strftime('%y-%m-%d %H:%M'),indices[s2].strftime('%y-%m-%d %H:%M'),oldest_left_border)
            #print ('\t',border)
            if border["right"] is None:
                if s1>=len(series)-N+1:
                    print (indices[s1].strftime('%H:%M'),indices[s2].strftime('%H:%M'),'RETURN list_borders',border)
                    #list_borders = filter_slope_leap(list_borders,series,diffs)
                    return list_borders

                if s1==len(series)-N:
                    
                    border["right"]=s1+N-1
                    border["t_right"]=series.index[s1+N-1]
                    border["right_A"]=list_fit_params[s1-1]["A"]
                    border["right_B"]=list_fit_params[s1-1]["B"]
                    if border['left']==None:
                        border['left']=oldest_left_border
                        border['t_left']=series.index[oldest_left_border]
                        border['left_A']=list_fit_params[oldest_left_border]["A"]
                        border['left_B']=list_fit_params[oldest_left_border]["B"]
                        list_borders.append(border)
                    print (indices[s1].strftime('%H:%M'),indices[s2].strftime('%H:%M'),'s1==len(series)-N',border)
                    s2=s1-1
                elif list_fit_params[s1]["right_dev"]>tol_deviation or (list_fit_params[s1]["right_dev"]>tol_sig*list_fit_params[s1]["max_dev"] and list_fit_params[s1]["right_dev"]>0.04*self.resolution):
                    
                    
                    border["right"]=s1+N-1
                    border["t_right"]=series.index[s1+N-1]
                    if s1!=0:
                        border["right_A"]=list_fit_params[s1]["A"]
                        border["right_B"]=list_fit_params[s1]["B"]
                    print (indices[s1].strftime('%H:%M'),indices[s2].strftime('%H:%M'), 'found gap, while border right is None',border)
                    s2=s1
                else:
                    print (indices[s1].strftime('%H:%M'),indices[s2].strftime('%H:%M'),'s1+1',border)
                    s1+=1
            else:


                print (s1)
                ## Conditions for left border
                # s2 cannot be negative, if at 0 this is a begining of segment
                s2_null = s2==0
                # If s2 reached right border of segment found at its left, must stop
                s2_reached_last_segment=False
                if len(list_borders)>0: s2_reached_last_segment = s2==list_borders[-1]["right"]
                # Conditions based on tolerance threshold
                s2_above_tolerance = False
                s2_std = False
                if list_fit_params[s2]["left_dev"]!=None:
                    # Next point is clearly above tolerance
                    s2_above_tolerance = list_fit_params[s2]["left_dev"]>tol_deviation
                    # Filter by standard deviation of linear left segment
                    s2_std = list_fit_params[s2]["left_dev"]>tol_sig*list_fit_params[s2]["max_dev"] and list_fit_params[s2]["left_dev"]>0.04*self.resolution

                if s2_null or s2_reached_last_segment or s2_above_tolerance or s2_std:

                    if border["right"]-s2<=N:
                        if len(retroactive_borders)!=0:
                            for b in retroactive_borders:
                                #print ('   ',b['t_left'],b['t_right'])
                                print ('Increment list border A',b)
                                list_borders.append(b)
                            retroactive_borders=[]
                            #if debug: plot_leap(diffs,series,s2,list_fit_params[s2]["A"],list_fit_params[s2]["B"],N,list_borders,"forward")

                        border={"left":None,"t_left":None,"right":None,"t_right":None,"left_A":None,"left_B":None,"right_A":None,"right_B":None}
                        print (indices[s1].strftime('%H:%M'),indices[s2].strftime('%H:%M'),'borderright-s2<=N',border)
                        s1=s1+N
                        continue
                    
                    ## Patch, for the case a left border was detected and leave a big segment behind that was healthy
                    
                    
                    
                    border["left"]=s2
                    border["t_left"]=series.index[s2]
                    border["left_A"]=list_fit_params[s2]["A"]
                    border["left_B"]=list_fit_params[s2]["B"]

                    print (indices[s1].strftime('%H:%M'),indices[s2].strftime('%H:%M'),'assing s2 to left',border)
                    
                    retroactive_borders.append(border)
                    if oldest_left_border<s2-N:
                         
                        border_left = {"left":oldest_left_border,"t_left":indices[oldest_left_border],"right":s2-2,
                        		"t_right":indices[s2-2],"left_A":list_fit_params[oldest_left_border]["A"],"left_B":list_fit_params[oldest_left_border]["B"],
                        		"right_A":list_fit_params[s2-2]["A"],"right_B":list_fit_params[s2-2]["B"]}
                        retroactive_borders.append(border_left)
                        print (indices[s1].strftime('%H:%M'),indices[s2].strftime('%H:%M'), 'Assign olders_left to left:', oldest_left_border,border)
                    oldest_left_border = border["right"]+1
                    
                    # If right border of last segment is not reached, must keep
                    # searching backward for a new segment
                    if len(list_borders)>0 and border["left"]-list_borders[-1]["right"]>N:
                        if s2-1>=N: 
                            border={"left":None,"t_left":None,"right":s2-1,"t_right":series.index[s2-1],"left_A":None,"left_B":None,"right_A":list_fit_params[s2-N]["A"],"right_B":list_fit_params[s2-N]["B"]}
                            print (indices[s1].strftime('%H:%M'),indices[s2].strftime('%H:%M'),'Assign s2-N to right',border)
                        else:
                            border={"left":None,"t_left":None,"right":s2-1,"t_right":series.index[s2-1],"left_A":None,"left_B":None,"right_A":None,"right_B":None}
                            print (indices[s1].strftime('%H:%M'),indices[s2].strftime('%H:%M'),'Assign s2-1 to right',border)
                        s2-=N
                        continue
                    retroactive_borders.reverse()
                    for b in retroactive_borders:
                        #print ('   ',b['t_left'],b['t_right'])
                        print ('Increment list border B',b)
                        list_borders.append(b)
                    retroactive_borders=[]

                    if s1>=len(series)-N-1:
                        #print ('HERE1')
                        
                        #list_borders = filter_slope_leap(list_borders,series,diffs)
                        #print (list_borders)
                        print ('return list_borders because right border has been found',indices[s1],indices[s2],border)
                        for b in list_borders: print (b['t_left'],b['t_right'])
                        return list_borders
                    border={"left":None,"t_left":None,"right":None,"t_right":None,"left_A":None,"left_B":None,"right_A":None,"right_B":None}
                    s1=s1+N
                else: s2-=1
        
        #print ('HERE2')
        #list_borders = filter_slope_leap(list_borders,series,diffs)
        print ('return, End of function')
        for b in list_borders: print (b['t_left'],b['t_right'])
        return list_borders
        '''



    def list_leaps(self,df):
        #borders_l,borders_p=[],[]
        #stec_l = ''
        #stec_p = ''
        #for c in df.columns:
        #    if "STEC_l" in c: 
        #        stec_l=c
        #    if "STEC_p" in c: 
        #        stec_p=c
                #break
                
        #print ('sl')
        borders_sl = self.list_leaps_series(df['STEC_l'],1/60,3,3)
        #for b in borders_sl: print (b['t_left'],b['t_right'])

        #borders_slp = self.list_leaps_series(df['STEC_p'],tol_dev=2.5/6.,tol_sig=6,N=8)
        
        #print (df)
        #print ('*******************************************')
        #print ('sl')
        #for bord in borders_sl: print (bord)
        #print ('*******************************************') 	
        #print ('slp')
        #for bord in borders_slp: print (bord)

        sys.exit()

        # We look at the borders found for slp and sll leaps, and make the
        # match
        # For each border found for sll
        #if len(borders_sl)!=0:
        #    for bsll in borders_sl:
        #        # For each border found for slp
        #        for bslp in borders_slp:
        #            # Test if left border of sll in inside slp borders item
        #            if bsll["t_left"]>=bslp["t_left"] and bsll["t_left"]<bslp["t_right"]:
        #                # Left of new border must be the one of bsll
        #                newb = bsll.copy()
        #                newb["t_left"]=bsll["t_left"]
        #                # We asign to the new border the earliest right border
        #                if bsll["t_right"]<bslp["t_right"]: newb["t_right"]=bsll["t_right"]
        #                else: newb["t_right"]=bslp["t_right"]
        #                borders_l.append(newb)
        #                break
        #            # If last condition not fulfilled, test in the right sll border
        #            # is inside the slp border item
        #            elif bsll["t_right"]>bslp["t_left"] and bsll["t_right"]<=bslp["t_left"]:
        #                newb = bsll.copy()
        #                newb["t_right"]=bsll["t_right"]
        #                #Left border must be outside blsp, otherwise it would
        #                #have been seen by the previous condition
        #                newb["t_left"]=bslp["t_left"]
        #                borders_l.append(newb)
        #                break
        #            # Case the bsll fully contains bslp, take bslp borders
        #            elif bsll["t_left"]<bslp["t_left"] and bsll["t_right"]>bslp["t_right"]:
        #                borders_l.append(bslp)
        #                break
                    


        return borders_sl,borders_slp

    def add_baseline(self):
    
        N_min_stec_p = 10
        signal_strength_threshold = 38
    
        for const, df_data in self.list_df.items():

            if len(df_data)==0: continue

            df_data = df_data.dropna(subset=["elevation"])

            t_begin = df_data.index[0]
            t_end = df_data.index[-1]

            df_filtered = pd.DataFrame()

            list_sv = df_data["sv"].unique().tolist()

            for sat in list_sv:
            
                if sat!='C14': continue
                print (sat)


                df_sat = df_data[df_data["sv"]==sat]
                df_sat_filter = pd.DataFrame()

                for channel in self.channels[const]:

                    if channel['C2']!='C6I': continue
                    chan_filter = (df_sat["C1"]==channel["C1"]) & (df_sat["C2"]==channel["C2"])
                    df_sat_chanel = df_sat[chan_filter]
                        
                    if len(df_sat_chanel["STEC_l"].dropna())==0: continue

                    # Make list of arcs of satellite using elevation information
                    elevations = df_sat_chanel["elevation"]
                    list_arcs = gnss.get_arcs(elevations,t_begin,t_end)

                    # Squared of sin of elevation for baseline (Brs) pondering
                    df_sat_chanel["sin2_ele"] = np.sin(df_sat_chanel['elevation'])**2
                    

                    # Individual values of the sum for future baseline calculation
                    #df_sat_chanel['BRs'] = (df_sat_chanel["STEC_p"]-df_sat_chanel["STEC_l"])*df_sat_chanel["sin2_ele"]
                    
                    # cos(chi) to calculate VTEC from STEC
                    df_sat_chanel['cos_chi'] = np.cos(np.arcsin(R_E*np.cos(df_sat_chanel["elevation"])/(R_E+self.h)))
                    n_arc=0

                    df_arcs = pd.DataFrame()
                    iarc = 0
                    for arc in list_arcs:
                        

                    
                        #if iarc==0: 
                        #    iarc+=1
                        #    continue
                        arc_as_baseline = False
                        # Must gets borders that are inside the arc segment
                        ##arc_borders = None
                        #arc_borders = []

                        # Extract data from satellite arc
                        date_filter = (df_sat_chanel.index>=arc["start"]) & (df_sat_chanel.index<=arc["end"])
                        df_arc = df_sat_chanel[date_filter]

                        if len(df_arc)<=8:continue      
                                     
                        arc_borders_l, arc_borders_p = self.list_leaps(df_arc)

                        
                        if (len(arc_borders_l)==0) or (len(arc_borders_p)==0): continue
                        
                        list_diff_S1, list_diff_S2  = [], []
                        
                        date_border_filter  = (df_arc.index>=arc_borders_p[0]["t_left"]) & (df_arc.index<=arc_borders_p[0]["t_right"])
                        df_border = df_arc[date_border_filter]
                    
                        signal_strength1 = (df_border['S1'] * df_border['sin2_ele']).sum() / df_border['sin2_ele'].sum()
                        signal_strength2 = (df_border['S2'] * df_border['sin2_ele']).sum() / df_border['sin2_ele'].sum()
                        
                        was_signal_cancelled = False if (signal_strength1>signal_strength_threshold) and (signal_strength2>signal_strength_threshold) else True
                        if was_signal_cancelled:
                            df_arc.loc[arc_borders_p[0]["t_left"]:arc_borders_p[0]["t_right"], 'STEC_p'] = np.nan
                        #print (df_arc.loc[arc_borders_p[0]["t_left"]:arc_borders_p[0]["t_right"], 'STEC_p'].tolist())
                        
                        has_baseline = True if was_signal_cancelled else False
                        
                        for i in range(1,len(arc_borders_p)):
                            
                            date_border_filter  = (df_arc.index>=arc_borders_p[i]["t_left"]) & (df_arc.index<=arc_borders_p[i]["t_right"])
                            df_border = df_arc[date_border_filter]
                            
                            signal_strength1 = (df_border['S1'] * df_border['sin2_ele']).sum() / df_border['sin2_ele'].sum()
                            signal_strength2 = (df_border['S2'] * df_border['sin2_ele']).sum() / df_border['sin2_ele'].sum()
                            
                            is_strong_signal_1 = True if (signal_strength1>signal_strength_threshold) and (signal_strength2>signal_strength_threshold) else False

                            #dt = (df_arc.index[arc_borders_p[i+1]['left']] - df_arc.index[arc_borders_p[i]['right']]).seconds
                            diff_S1 = df_arc['S1'].iloc[arc_borders_p[i]['left']] - df_arc['S1'].iloc[arc_borders_p[i-1]['right']]
                            diff_S2 = df_arc['S2'].iloc[arc_borders_p[i]['left']] - df_arc['S2'].iloc[arc_borders_p[i-1]['right']]
                            
                            # If there is a drop in signal strength, regarless the signal strength, the segment that dropped is unreliable 
                            if diff_S1<-5 or diff_S2<-5:
                                iright = i+1 if i<len(arc_borders_p)-1 else i
                                ileft = i-1 if i>0 else 0
                                df_arc.loc[arc_borders_p[ileft]["t_right"]:arc_borders_p[iright]["t_left"], 'STEC_p'] = np.nan
                                was_signal_cancelled=True

                                
                            # If there was no notorious drop, we look if the previous signal was cancelled
                            #   if in addition there was no recovery, segment i is cancelled
                            elif was_signal_cancelled:
                                if is_strong_signal_1: 
                                    was_signal_cancelled=False
                                    has_baseline=True
                                else:
                                    iright = i+1 if i<len(arc_borders_p)-1 else i
                                    ileft = i-1 if i>0 else 0
                                    df_arc.loc[arc_borders_p[ileft]["t_right"]:arc_borders_p[iright]["t_left"], 'STEC_p'] = np.nan
                                    was_signal_cancelled=True
                            else:
                                has_baseline=True



                        print (has_baseline)
                        if not has_baseline: continue
                            
                            #print (diff_S1,diff_S2,was_signal_cancelled,is_strong_signal_1,signal_strength1,signal_strength2)
                            #print (arc_borders_p[i]["t_left"],arc_borders_p[i]["t_right"])
                            #print (df_arc.loc[arc_borders_p[i]["t_left"]:arc_borders_p[i]["t_right"], 'STEC_p'].tolist()[:10])
                            
                            #is_next_strong_signal = is_strong_signal

                       #

                            
                            #print (signal_strength1,signal_strength2,border["t_left"],border["t_right"])
                            
                            #if (signal_strength1<40) or (signal_strength2<40):
                            #    df_arc.loc[border["t_left"]:border["t_right"], 'STEC_p'] = np.nan
                                
                                
                        #### Make a fit of STEC_p to detect deviations of STEC_l
                        x = (df_arc.index - df_arc.index[0]).total_seconds()
                        y = df_arc["STEC_p"].to_numpy()
                        w = df_arc["sin2_ele"].to_numpy()

                        # 1. Masque booléen pour exclure les NaN de y (et s'assurer que x et w sont valides)
                        valid_mask = ~np.isnan(y) & ~np.isnan(x) & ~np.isnan(w)

                        # 2. Filtrage des données pour le fit
                        x_fit = x[valid_mask]
                        y_fit = y[valid_mask]
                        w_fit = w[valid_mask]

                        # 3. Ajustement dynamique du degré basé sur le nombre de points VALIDES
                        degree_fit = min(len(x_fit) - 1, 10)

                        # 4. Calcul du fit si on a assez de points valides
                        if degree_fit >= 0:
                            try: 
                                coeffs = np.polyfit(x_fit, y_fit, deg=degree_fit, w=w_fit)
                            except: 
                                print ('STEC_p fit did not work\n')
                                continue
                        else:
                            print ('Not enough valid data points for fit for degree '+str(degree_fit)+' '+str(len(x))+' '+str(len(x_fit))+'\n')
                            continue

                        # 5. Évaluation sur l'axe 'x' complet (conserve les dimensions d'origine du DataFrame)
                        poly = np.poly1d(coeffs)
                        #df_arc["STEC_p_fit"] = poly(x)      
                        df_arc.loc[valid_mask, "STEC_p_fit"] = poly(x_fit)                          

                        sigma = (df_arc["STEC_p"] - df_arc["STEC_p_fit"]).std()

                        if len(arc_borders_l)==0: continue

        
                        #### Clasify segments #####
                        has_sane_parts = False
                        sane_indices = []
                        #list_size_segs = []
                        d_df_seg = []
                        segment_anchored = []

                        list_pop_border = []


                        i = 0
                        while i<len(arc_borders_l):
                            print (arc_borders_l[i])
  
                            date_border_filter  = (df_arc.index>=arc_borders_l[i]["t_left"]) & (df_arc.index<=arc_borders_l[i]["t_right"])
                            df_seg = df_arc[date_border_filter]

                            brs = None
                            N_valid_STEC_p = len(df_seg['STEC_p'].dropna())
                            #print (df_seg.index[0])
                            if N_valid_STEC_p>N_min_stec_p:
                                df_seg_br = df_seg.dropna(subset=['STEC_p'])
                                df_seg_br['BRs'] = (df_seg_br["STEC_p"]-df_seg_br["STEC_l"])*df_seg_br["sin2_ele"]
                                brs=df_seg_br['BRs'].sum(skipna=True)/df_seg_br["sin2_ele"].sum(skipna=True)
                                #print (df_seg_br)
                                df_seg['STEC_l'] = df_seg['STEC_l']+brs
                                segment_anchored.append(True)
                                #if (N_valid_STEC_p == len(df_seg['STEC_p'])) and (max(df_seg['STEC_l']-df_seg['STEC_p_fit'])>5*sigma): 
                                #    df_seg["STEC_l"] = df_seg['STEC_p_fit']
                            else:
                                segment_anchored.append(False)

                            if i>0 and segment_anchored[i-1] and segment_anchored[i]:
                                d_df_seg[i-1] = pd.concat([d_df_seg[i-1],df_seg])
                                arc_borders_l[i-1]['t_right'] = arc_borders_l[i]['t_right'] 
                                arc_borders_l[i-1]['right_A'] = arc_borders_l[i]['right_A'] 
                                arc_borders_l[i-1]['right_B'] = arc_borders_l[i]['right_B'] 
                                arc_borders_l.pop(i)
                                segment_anchored.pop(i)
                                i = i-1
                            else:
                                d_df_seg.append(df_seg)
                            i = i+1

                        #i=1
                        #while i<len(d_df_seg):
                        #    d_df_seg[0] = pd.concat([d_df_seg[0],d_df_seg[i]])
                        #    print (d_df_seg[0])
                        #    d_df_seg.pop(i)
                        #continue
                        
                        print (segment_anchored)
                        for i in range(len(segment_anchored)):
                            print (segment_anchored[i],arc_borders_l[i])
                        # If no STEC_l segment could be anchored in STEC_p then there is nothing useful, go to next arc.
                        if True not in segment_anchored: continue


                        #if len(segment_anchored)<=1: continue 

                        while len(d_df_seg)!=1:

                            
                        #for i in range(1,len(d_df_seg)):
                            #print ('------------------------------------------')
                            #print (i,segment_anchored[0],segment_anchored[1])
                            #print (1,len(d_df_seg),len(segment_anchored))
                            #print (d_df_seg[0][['STEC_l','STEC_p']])
                            #print (d_df_seg[1][['STEC_l','STEC_p']])

                            # Case segment i is anchored: we see if left neightbour need rescue
                            if segment_anchored[1]:
                                # Case left segment is NOT anchored (otherwise nothing to do)
                                if not segment_anchored[0]:
                                    #t_dist=(arc_borders_l[1]["t_right"]-arc_borders_l[0]["t_left"]).seconds
                                    #slope = (arc_borders_l[1]["right_A"]+arc_borders_l[0]["left_A"])/2        
                                    t_dist=(arc_borders_l[1]["t_left"]-arc_borders_l[0]["t_right"]).seconds
                                    slope = (arc_borders_l[1]["left_A"]+arc_borders_l[0]["right_A"])/2        
                                    glue=d_df_seg[1]["STEC_l"].iloc[0]-d_df_seg[0]["STEC_l"].iloc[-1]-slope*t_dist
                                    d_df_seg[0]["STEC_l"] = d_df_seg[0]["STEC_l"] + glue
                                    
                                    #print ('0 1',t_dist,slope,glue)
                                    #print (d_df_seg[0])
                                    #print (d_df_seg[1])
                                    #d_df_seg[0] = pd.concat([d_df_seg[0],d_df_seg[1]])
                                    #d_df_seg.pop(i)
                                    #arc_borders_l[i-1]['t_right'] = arc_borders_l[i]['t_right'] 
                                    #arc_borders_l[i-1]['right_A'] = arc_borders_l[i]['right_A'] 
                                    #arc_borders_l[i-1]['right_B'] = arc_borders_l[i]['right_B'] 
                                    #arc_borders_l.pop(i)
                                # if both segment are anchored they will be merged in one, then we remove 
                                #   the one that says false which by definition is the first one
                                segment_anchored.pop(0)


                            # Case segment i is not anchored, we seek for help at the left, we paste anyway
                            else:
                                t_dist=(arc_borders_l[1]["t_left"]-arc_borders_l[0]["t_right"]).seconds
                                slope = (arc_borders_l[1]["left_A"]+arc_borders_l[0]["right_A"])/2        
                                glue=float(d_df_seg[1]["STEC_l"].iloc[0]-d_df_seg[0]["STEC_l"].iloc[-1]-slope*t_dist)
                                #print ('d_df_seg[1]["STEC_l"].iloc[0]=  ',d_df_seg[1]["STEC_l"].iloc[0])
                                #print ('d_df_seg[0]["STEC_l"].iloc[-1]=  ',d_df_seg[0]["STEC_l"].iloc[-1])
                                #print ('slope = ',slope,'t_dist = ',t_dist,'glue = ',glue)
                                
                                # Case the left segment is anchored then the i segment must be shifted
                                if segment_anchored[0]:
                                    #print ('1 0',t_dist,slope,glue,type(glue))
                                    #print (d_df_seg[0])
                                    d_df_seg[1]["STEC_l"] = d_df_seg[1]["STEC_l"] - glue
                                    #print (d_df_seg[1])
                                else:
                                    #print ('0 0',t_dist,slope,glue)
                                    d_df_seg[0]["STEC_l"] = d_df_seg[0]["STEC_l"] + glue
                                ## Here we remove the second bool, if both are not anchored we want false, 
                                #   else we want to keep the true value from the first one
                                segment_anchored.pop(1)

                            d_df_seg[0] = pd.concat([d_df_seg[0],d_df_seg[1]])
                            d_df_seg.pop(1)
                            arc_borders_l[0]['t_right'] = arc_borders_l[1]['t_right'] 
                            arc_borders_l[0]['right_A'] = arc_borders_l[1]['right_A'] 
                            arc_borders_l[0]['right_B'] = arc_borders_l[1]['right_B'] 
                            arc_borders_l.pop(1)



                            #i = i-1
                                
                            #print (d_df_seg[0])

                        
                                




                            



                        '''
                        for i,b in enumerate(arc_borders_l):
                        
                            size_seg = (arc_borders_l[i]["t_right"]-arc_borders_l[i]["t_left"]).seconds
                            #list_size_segs.append(size_seg)

        
                            if size_seg>50*60:
                                df_seg = None
                                date_border_filter  = (df_arc.index>=arc_borders_l[i]["t_left"]) & (df_arc.index<=arc_borders_l[i]["t_right"])
                                df_seg = df_arc[date_border_filter]

                                #df_br = df_seg.dropna(subset=['BRs',"sin2_ele"])

                                #We take only long enough segments of STEC_p to trust the signal 
                                # if not then we take STEC_l as it is 
                                if len(df_seg['STEC_p'].dropna())>N_min_stec_p:
                                    print ('Using baseline here 11',b['t_left'],b['t_right'])
                                    arc_has_baseline = True
                                    df_seg_br = df_seg.dropna(subset=['STEC_p'])
                                    df_seg_br['BRs'] = (df_seg_br["STEC_p"]-df_seg_br["STEC_l"])*df_seg_br["sin2_ele"]
                                    brs=df_seg_br['BRs'].sum(skipna=True)/df_seg_br["sin2_ele"].sum(skipna=True)
                                    df_seg["STEC_l"] = df_seg["STEC_l"] +  brs
                                    df_seg['diffSTEC'] = df_seg['STEC_l']-df_seg['STEC_p_fit']
                                    ## In case STEC_l corrected with is going too far away from the true value we the part that was adjusted with STEC_p.
                                    if max(df_seg['diffSTEC'])>5*sigma: df_seg["STEC_l"] = df_seg['STEC_p_fit']
                                    segment_anchored.append(True)
                                else: segment_anchored.append(False)
                                     
                                d_df_seg.append(df_seg)
                                has_sane_parts = True
                                sane_indices.append(i)
                            else:
                                segment_anchored.append(False)
                                d_df_seg.append(None)
                        if not has_sane_parts: 
                            continue

                        print (sane_indices,len(arc_borders_l))
                        
                        # paste left part of arc
                        for i in range(sane_indices[0]-1,-1,-1):
                            df_seg = None
                            date_border_filter  = (df_arc.index>=arc_borders_l[i]["t_left"]) & (df_arc.index<=arc_borders_l[i]["t_right"])
                            df_seg = df_arc[date_border_filter]
                            print (i,arc_borders_l[i]["t_left"],arc_borders_l[i]["t_right"])
                            
                            brs = None
                            N_valid_STEC_p = len(df_seg['STEC_p'].dropna())
                            if N_valid_STEC_p>N_min_stec_p:
                                df_seg_br = df_seg.dropna(subset=['STEC_p'])
                                df_seg_br['BRs'] = (df_seg_br["STEC_p"]-df_seg_br["STEC_l"])*df_seg_br["sin2_ele"]
                                brs=df_seg_br['BRs'].sum(skipna=True)/df_seg_br["sin2_ele"].sum(skipna=True)
                                
                            

                            t_dist=(arc_borders_l[i+1]["t_left"]-arc_borders_l[i]["t_right"]).seconds
                            slope = (arc_borders_l[i+1]["left_A"]+arc_borders_l[i]["right_A"])/2        
                            glue=d_df_seg[i+1]["STEC_l"].iloc[0]-df_seg["STEC_l"].iloc[-1]-slope*t_dist
                            #df_br = df_seg.dropna(subset=['BRs',"sin2_ele"])
                            
                            # Case no baseline could be established
                            if brs is None:
                                df_seg["STEC_l"] = df_seg["STEC_l"] + glue
                                print ('Using glue 11',arc_borders_l[i]['t_left'],arc_borders_l[i]['t_right'],glue)
                            else:
                                arc_has_baseline = True
                                df_seg['diffSTEC'] = df_seg['STEC_l']-df_seg['STEC_p_fit']
                                if (N_valid_STEC_p == len(df_seg['STEC_p'])) and (max(df_seg['diffSTEC'])>5*sigma): 
                                    df_seg["STEC_l"] = df_seg['STEC_p_fit']
                                    segment_anchored[i]=True
                                elif abs(glue-brs)>20:
                                    print ('Using baseline here 22',arc_borders_l[i]['t_left'],arc_borders_l[i]['t_right'])
                                    df_seg["STEC_l"] = df_seg["STEC_l"] + brs
                                    segment_anchored[i]=True
                                else: 
                                    print ('Using glue 11.2',arc_borders_l[i]['t_left'],arc_borders_l[i]['t_right'],'glue=',glue,'brs=',brs)
                                    df_seg["STEC_l"] = df_seg["STEC_l"] + glue
                            
                            #df_seg["STEC_l"] = df_seg["STEC_l"] + brs
                            #df_seg['diffSTEC'] = df_seg['STEC_l']-df_seg['STEC_p_fit']
                            #if max(df_seg['diffSTEC'])>5*sigma: df_seg["STEC_l"] = df_seg['STEC_p_fit']
                            d_df_seg[i] = df_seg

                        # paste right part of arc
                        for i in range(sane_indices[-1]+1,len(arc_borders_l)):
                            df_seg = None
                            date_border_filter  = (df_arc.index>=arc_borders_l[i]["t_left"]) & (df_arc.index<=arc_borders_l[i]["t_right"])
                            df_seg = df_arc[date_border_filter]
                            
                            brs = None
                            N_valid_STEC_p = len(df_seg['STEC_p'].dropna())
                            if N_valid_STEC_p>N_min_stec_p:
                                df_seg_br = df_seg.dropna(subset=['STEC_p'])
                                df_seg_br['BRs'] = (df_seg_br["STEC_p"]-df_seg_br["STEC_l"])*df_seg_br["sin2_ele"]
                                brs=df_seg_br['BRs'].sum(skipna=True)/df_seg_br["sin2_ele"].sum(skipna=True)
                                
                            #df_seg['BRs'] = (df_seg["STEC_p"]-df_seg["STEC_l"])*df_seg["sin2_ele"]
                            
                            t_dist=(arc_borders_l[i]["t_left"]-arc_borders_l[i-1]["t_right"]).seconds
                            slope = arc_borders_l[i-1]["right_A"]       
                            glue=d_df_seg[i-1]["STEC_l"].iloc[-1]-df_seg["STEC_l"].iloc[0]+slope*t_dist
                            
                            #df_br = df_seg.dropna(subset=['BRs',"sin2_ele"])
                            #brs=df_seg['BRs'].sum(skipna=True)/df_seg["sin2_ele"].sum(skipna=True)
                            # Case no baseline could be established
                            if brs is None:
                                df_seg["STEC_l"] = df_seg["STEC_l"] + glue
                                print ('Using glue 22',arc_borders_l[i]['t_left'],arc_borders_l[i]['t_right'])
                            else:
                                arc_has_baseline = True
                                df_seg['diffSTEC'] = df_seg['STEC_l']-df_seg['STEC_p_fit']
                                if (N_valid_STEC_p == len(df_seg['STEC_p'])) and (max(df_seg['diffSTEC'])>5*sigma): 
                                    df_seg["STEC_l"] = df_seg['STEC_p_fit']
                                elif abs(glue-brs)>20:
                                    print ('Using baseline here 33')
                                    df_seg["STEC_l"] = df_seg["STEC_l"] + brs
                                else: 
                                    df_seg["STEC_l"] = df_seg["STEC_l"] + glue
                            
                            #if not (brs is None) and (abs(glue-brs)>20): correction = brs
                            #else: correction = glue
                            #df_seg["STEC_l"] = df_seg["STEC_l"] + correction
                            ##df_seg["STEC_l"] = df_seg["STEC_l"] + brs
                            #df_seg['diffSTEC'] = df_seg['STEC_l']-df_seg['STEC_p_fit']
                            #if max(df_seg['diffSTEC'])>5*sigma: df_seg["STEC_l"] = df_seg['STEC_p_fit']
                            d_df_seg[i] = df_seg

                        # Look at intermediate segments
                        i = 0
                        if len(sane_indices)>=2:
                            for i,s in enumerate(sane_indices[0:len(sane_indices)-1]):
                                ileft = s+1
                                iright = sane_indices[i+1]-1
                                only_right = False
                                only_left = False
                                while ileft<=iright:
                                    if not only_left:
                                        t_left = (arc_borders_l[ileft]["t_left"]-arc_borders_l[ileft-1]["t_right"]).seconds
                                        t_right = (arc_borders_l[iright]["t_left"]-arc_borders_l[ileft]["t_right"]).seconds
                                        if t_left<=t_right or only_right:
                                            df_seg = None
                                            date_border_filter  = (df_arc.index>=arc_borders_l[ileft]["t_left"]) & (df_arc.index<=arc_borders_l[ileft]["t_right"])
                                            df_seg = df_arc[date_border_filter]

                                            brs = None
                                            N_valid_STEC_p = len(df_seg['STEC_p'].dropna())
                                            if N_valid_STEC_p>N_min_stec_p:
                                                df_seg_br = df_seg.dropna(subset=['STEC_p'])
                                                df_seg_br['BRs'] = (df_seg_br["STEC_p"]-df_seg_br["STEC_l"])*df_seg_br["sin2_ele"]
                                                brs=df_seg_br['BRs'].sum(skipna=True)/df_seg_br["sin2_ele"].sum(skipna=True)
                                            
                                            #df_seg['BRs'] = (df_seg["STEC_p"]-df_seg["STEC_l"])*df_seg["sin2_ele"]
                                            
                                            #df_seg = df_arc.loc[arc_borders[ileft]["t_left"]:arc_borders[ileft]["t_right"]]
                                            t_dist=(arc_borders_l[ileft]["t_left"]-arc_borders_l[ileft-1]["t_right"]).seconds        
                                            slope = arc_borders_l[ileft-1]["right_A"]        
                                            glue=d_df_seg[ileft-1]["STEC_l"].iloc[-1]-df_seg["STEC_l"].iloc[0]+slope*t_dist

                                            if brs is None:
                                                print ('Using glue 33',arc_borders_l[i]['t_left'],arc_borders_l[i]['t_right'])
                                                df_seg["STEC_l"] = df_seg["STEC_l"] + glue
                                            else:
                                                arc_has_baseline = True
                                                df_seg['diffSTEC'] = df_seg['STEC_l']-df_seg['STEC_p_fit']
                                                if (N_valid_STEC_p == len(df_seg['STEC_p'])) and (max(df_seg['diffSTEC'])>5*sigma):  
                                                    df_seg["STEC_l"] = df_seg['STEC_p_fit']
                                                elif abs(glue-brs)>20:
                                                    print ('Using baseline here 44')
                                                    df_seg["STEC_l"] = df_seg["STEC_l"] + brs
                                                else: 
                                                    df_seg["STEC_l"] = df_seg["STEC_l"] + glue                                            
                                            #df_br = df_seg.dropna(subset=['BRs',"sin2_ele"])
                                            #brs=df_seg['BRs'].sum(skipna=True)/df_seg["sin2_ele"].sum(skipna=True)
                                            #if abs(glue-brs)>20: correction = brs
                                            #else: correction = glue
                                            #df_seg["STEC_l"] = df_seg["STEC_l"] + correction
                                            #df_seg["STEC_l"] = df_seg["STEC_l"] + brs
                                            #df_seg['diffSTEC'] = df_seg['STEC_l']-df_seg['STEC_p_fit']
                                            #if max(df_seg['diffSTEC'])>5*sigma: df_seg["STEC_l"] = df_seg['STEC_p_fit']
        
                                            d_df_seg[ileft] = df_seg
                                            ileft+=1
                                        elif not only_right:
                                            only_left=True
                                        if ileft>iright: break
                                    if not only_right:
                                        t_left = (arc_borders_l[iright]["t_left"]-arc_borders_l[ileft]["t_right"]).seconds
                                        t_right = (arc_borders_l[iright+1]["t_left"]-arc_borders_l[iright]["t_right"]).seconds
                                        if t_left>=t_right or only_left:
                                            df_seg = None
                                            date_border_filter  = (df_arc.index>=arc_borders_l[iright]["t_left"]) & (df_arc.index<=arc_borders_l[iright]["t_right"])
                                            df_seg = df_arc[date_border_filter]

                                            brs = None
                                            N_valid_STEC_p = len(df_seg['STEC_p'].dropna())
                                            if N_valid_STEC_p>N_min_stec_p:
                                                df_seg_br = df_seg.dropna(subset=['STEC_p'])
                                                df_seg_br['BRs'] = (df_seg_br["STEC_p"]-df_seg_br["STEC_l"])*df_seg_br["sin2_ele"]
                                                brs=df_seg_br['BRs'].sum(skipna=True)/df_seg_br["sin2_ele"].sum(skipna=True)
                                                                                            
                                            #df_seg['BRs'] = (df_seg["STEC_p"]-df_seg["STEC_l"])*df_seg["sin2_ele"]
                                            #df_seg = df_arc.loc[arc_borders[iright]["t_left"]:arc_borders[iright]["t_right"]]
                                            t_dist=(arc_borders_l[iright+1]["t_left"]-arc_borders_l[iright]["t_right"]).seconds        
                                            slope = arc_borders_l[iright+1]["right_A"]        
                                            glue=d_df_seg[iright+1]["STEC_l"].iloc[0]-df_seg["STEC_l"].iloc[-1]+slope*t_dist
                                            #df_br = df_seg.dropna(subset=['BRs',"sin2_ele"])
                                            #brs=df_br['BRs'].sum(skipna=True)/df_br["sin2_ele"].sum(skipna=True)
                                            #if abs(glue-brs)>20: correction = brs
                                            #else: correction = glue
                                            ##df_seg["STEC_l"] = df_seg["STEC_l"] + correction
                                            #df_seg["STEC_l"] = df_seg["STEC_l"] + brs
                                            #df_seg['diffSTEC'] = df_seg['STEC_l']-df_seg['STEC_p_fit']
                                            #if max(df_seg['diffSTEC'])>5*sigma: df_seg["STEC_l"] = df_seg['STEC_p_fit'] 
                                            if brs is None:
                                                print ('Using glue 44',arc_borders_l[i]['t_left'],arc_borders_l[i]['t_right'])
                                                df_seg["STEC_l"] = df_seg["STEC_l"] + glue
                                            else:
                                                arc_has_baseline = True
                                                df_seg['diffSTEC'] = df_seg['STEC_l']-df_seg['STEC_p_fit']
                                                if (N_valid_STEC_p == len(df_seg['STEC_p'])) and (max(df_seg['diffSTEC'])>5*sigma):  
                                                    df_seg["STEC_l"] = df_seg['STEC_p_fit']
                                                elif abs(glue-brs)>20:
                                                    print ('Using baseline here 55')
                                                    df_seg["STEC_l"] = df_seg["STEC_l"] + brs
                                                else: 
                                                    df_seg["STEC_l"] = df_seg["STEC_l"] + glue
                                                           
                                            d_df_seg[iright] = df_seg
                                            iright-=1
                                        elif not only_left:
                                            only_right=True
                        '''

                        #df_filter_arc = pd.DataFrame()
                        #for df in d_df_seg:
                        #    if df is None:continue
                        #    df_filter_arc = pd.concat([df_filter_arc,df])                        
                        
                        if len(d_df_seg)!=1: continue
                        df_arcs = pd.concat([df_arcs,d_df_seg[0]])    
        
                    if len(df_arcs)==0: continue
                    if len(df_sat_filter)==0: df_sat_filter = df_arcs
                    else:                      
                        df_sat_filter["STEC_l"] = df_arcs["STEC_l"]
                        df_sat_filter["STEC_p"] = df_arcs["STEC_p"]

                    df_sat_filter['C1'] = channel['C1']
                    df_sat_filter['C2'] = channel['C2']
                 
 
                    df_filtered = pd.concat([df_filtered,df_sat_filter])


            if len(df_filtered)==0: continue
            self.list_df[const] = df_filtered

            self.list_df[const]['time_day'] = self.list_df[const].index.normalize()
            self.list_df[const].reset_index(inplace=True)   # 'time' becomes a column

            p1p2_mask = (self.list_df[const]['C1'] == 'P1') & (self.list_df[const]['C2'] == 'P2')
            self.list_df[const].loc[p1p2_mask, 'time_day'] = self.list_df[const].loc[p1p2_mask, 'time_day'].dt.to_period('M').dt.to_timestamp()
    
            # Merge on daily time + sv + C1 + C2
            self.list_df[const] = self.list_df[const].merge(
                self.sat_dcb[['time', 'sv', 'C1', 'C2', 'dcb']],
                left_on=['time_day', 'sv', 'C1', 'C2'],
                right_on=['time','sv', 'C1', 'C2'],
                how='inner',
                suffixes=('', '_df2')
            )

            if len(self.list_df[const])==0: continue 
                
            # Restore original index and drop helper columns
            
            self.list_df[const] = self.list_df[const].set_index('time').drop(columns=['time_day', 'time_df2'])

            
            if const!='S':
                #self.list_df[const]["STEC_l"] += self.list_df[const]["dcb"]
                self.list_df[const]["VTEC"]=self.list_df[const]["STEC_l"]*self.list_df[const]['cos_chi']

    
    def compute_receiver_bias(self,resolution=datetime.timedelta(days=1)):
        ''' Sylvain Blunier 06/2021 v1.0.0
                Function that computes and returns the biais of some station.
                * Input: string of the name of the feather files that contains data for compute
                This file must be in directory "./feather/"
                'station_feather' feather file must at least possess columns:
                    time,elevation,STEC_sl
                * Output: float representing receiver bias
                Algorithm:
                This algorithm minimize the sum of variances between satellites at each time
                *section 4.1 of https://hal.archives-ouvertes.fr/hal-00317176/file/angeo-21-2083-2003.pdf
                *algebra: https://colab.research.google.com/drive/1UCZHR0t-9jyyjAnLuMN3N0Z2NB6tgI_l?usp=sharing
        '''
 
        for const, df in self.list_df.items():

            #list_stec = [s for s in df.columns if 'STEC_l' in s]
            if len(df)==0: continue

            br = process_br(df[["elevation",'STEC_l']],self.h)
        
            self.br['station'].append(self.station)
            self.br['time_i'].append(min(df.index))
            self.br['time_f'].append(max(df.index))
            self.br['constellation'].append(const)
            #self.br['C1'].append(channel["C1"])
            #self.br['C2'].append(channel["C2"])
            self.br['br'].append(br)

            df_br_station = pd.DataFrame(self.br)
            
            
            self.list_df[const]['time_day'] = self.list_df[const].index.normalize()
            self.list_df[const].reset_index(inplace=True)   # 'time' becomes a column

            self.list_df[const]['br'] = br
            
            df_br_station['time_day'] = df_br_station['time_i'].dt.normalize()
            df_br_station_unique = df_br_station[df_br_station['station']==self.station]
            df_br_station_unique = df_br_station_unique[df_br_station_unique['constellation']==const]

            
        self.df_br = pd.DataFrame(self.br)#.set_index("time_i")
        f_br = st.root_dir + "TEC/receiver_bias.csv"

        self.df_br.set_index("station",inplace=True)
        if os.path.exists(f_br):
            df_br_stored = pd.read_csv(f_br).set_index("station")
            df_br_stored = df_br_stored[df_br_stored.index!=self.station]
            self.df_br = pd.concat([df_br_stored,self.df_br])
        self.df_br.to_csv(f_br)


    def add_receiver_bias(self):
        self.compute_receiver_bias()
        for const in self.list_df.keys():

            if len(self.list_df[const])==0: continue 


            if const!="S":
                self.list_df[const]["VTEC"] = self.list_df[const]["STEC_l"]-self.list_df[const]["br"]
                self.list_df[const]["VTEC"] *= np.cos(np.arcsin(R_E*np.cos(self.list_df[const]["elevation"])/(R_E+self.h)))

                self.list_df[const].reset_index(inplace=True)
                self.list_df[const] = self.list_df[const][["time","sv","lat","lon","elevation","STEC_l","VTEC"]]
            else:
                self.list_df[const].reset_index(inplace=True)
                self.list_df[const] = self.list_df[const][["time","sv","lat","lon","elevation","STEC_l"]]

            self.list_df[const] = self.list_df[const].groupby(by=["time","sv"],as_index=False).mean()

    #def load_sat_DCB(self):

            

    #def load_sat_DCB(self):
    #
    #    # Case no DCB file is provided, we try downloading them
    #    if len(self.list_f_dcb) == 0:
    #    
    #        directory_DCB_path = Path(st.root_dir+'DCB/')
    #        list_sorted_DCB = []
    #        for file in directory_DCB_path.rglob("*"): 
    #            if file.is_file():
    #                list_sorted_DCB.append(str(file.resolve()))
    #        #list_sorted_DCB = sorted([file.resolve() for file in directory_DCB_path.rglob("*") if file.is_file()])
    #        
    #        for i in range((self.datemax.date() - self.datemin.date()).days + 1):
    #            d = self.datemin.date() + datetime.timedelta(days=i)
    #            year = d.year
    #            doy = (d - datetime.date(year,1,1)).days + 1
    #            DCB_file = st.root_dir+'DCB/' +str(year)+ "/CAS0MGXRAP_"+str(year)+str(doy)+"0000_01D_01D_DCB.BSX"
    #            print (DCB_file)
    #            print (list_sorted_DCB)
    #            if (DCB_file not in list_sorted_DCB):
    #                #print ("Should be here")
    #                DCB.get_dcb_from_cddis(year,doy)
    #            self.list_f_dcb.append(DCB_file)
    #    #else: 
    #        
    #
    #
    #    if len(self.list_f_dcb)==0: return pd.DataFrame()
    #    else:
    #        self.sat_dcb = DCB.load_dcb(self.list_f_dcb)
    #        return
        
    def to_feather(self):
        df_obs = pd.DataFrame()

        for const in self.list_df.keys():
            df_obs = pd.concat([
                df_obs,
                self.list_df[const]#[["time","sv","lat","lon","elevation","STEC_l","VTEC"]]
            ])
            

        feather_path = st.root_dir + "TEC/" + self.station
        df_obs.to_feather(feather_path+".feather")
    
    def compute_vtec(self,station):	
        self.list_df = {}
        self.channels = {}
     
        print ("RINEX to STEC")
        self.rinex_to_stec(station)
        
        print ("Add satellite position")
        self.add_satellite_pos()
        
        print ("Calculating baseline to correct Slant TEC")
        self.add_baseline()

        print ("Calculating receiver bias, correct Slant TEC, compute VTEC")
        self.add_receiver_bias()

        df_obs = pd.DataFrame()

        for const in self.list_df.keys():
            df_obs = pd.concat([
                df_obs,
                self.list_df[const]
            ])

        return df_obs
        
        #self.to_feather(st.root_dir + "TEC/" + str(self.year) + "/" + self.station + ".feather" )
        #self.to_feather()
    
