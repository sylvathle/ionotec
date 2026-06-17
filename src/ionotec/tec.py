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

from pathlib import Path

#from .stations import root_dir
from . import stations as st
from . import gnss
from . import rinex as rx
from . import freq
from . import DCB


import warnings

pd.options.mode.chained_assignment = None

# Earth Radius
R_E = 6371000



def fit_lin(t,sig):
    N=len(t)
    # Coefficients of paraboloid of error function (N*mse)
    a,b,c,d,e=0,N,0,0,0 # b=N for B^2 coef
    # Iterate over subseries to calculate coefficients
    for i in range(N):
        a+=t[i]**2 # A^2 coef
        c+=2*t[i] # A*B coef
        d-=2*t[i]*sig[i] # A coef
        e-=2*sig[i] # B coef

    # Forward A and B parameters of linear fit (solve the minimum of mse)
    A=-(2*b*d-c*e)/(4*a*b-c**2)
    B=-(c*d-2*a*e)/(c**2-4*a*b)

    
    sigma=0
    sigmas = []
    # Calculate value of err function
    for i in range(N):
        s = abs(sig[i]-A*t[i]-B)
        sigma+=s
        sigmas.append(s)
    max_deviation = max(sigmas)
    mean_deviation = np.mean(sigmas)

    return A,B,max_deviation,mean_deviation

def tleft(border):
    return border["t_left"]


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


def plot_leap(diffs,series,s,A,B,N,borders,title=""):
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
        if math.isnan(a) or math.isnan(b): sys.exit()

    if a==0: return float("nan")
    # root of error function = receiver bias.
    br = b/a
    return br

from pathlib import Path



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
    

    dict_f_obs = {}
    list_f_nav = []
    list_f_dcb = []


    source_data_folder = ""

    def __init__(self,rinex_folder,
                 date_min=None,date_max=None,
                 min_lat=None,min_lon=None,
                 h=350000):

        self.station = ""
        self.coord = ""

        # Estimated high of ionosphere
        self.h = h

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
            os.mkdir(st.root_dir + "TEC/", exist_ok=True)
        
        self.datemin = date_min
        self.datemax = date_max
        self.source_data_folder = rinex_folder
        self.prepare_files()

        self.gnss = gnss.gnss(self.list_f_nav)

        self.load_sat_DCB()
        self.sat_dcb.reset_index(inplace=True)

    ''' Function that prepare the obs, nav and bias list '''
    def prepare_files(self):

        # Define the target directory
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


    '''
    def __init__(self,list_f_obs,f_nav,list_f_sat_DCB=[],outfolder="output",resolution=60,h=350000):
    
        # List of observation files
        self.list_f_obs = list_f_obs 
        # List of navigation files
        self.f_nav = f_nav
        # File containing satellite bias
        self.list_f_sat_DCB = list_f_sat_DCB
        #Output feather file
        #if f_out == "": 
        #    self.f_feather = self.list_f_obs[min(len(self.list_f_obs)-1,1)][:-4]+"tec.feather"
        #else:
        #    if f_out[-8:] ==".feather": self.f_feather = f_out
        #    else: self.f_feather = f_out + ".feather"
    
        self.sv_gps = [str(i+1) for i in range(32)]
        for i,sv in enumerate(self.sv_gps):
            if i<9: self.sv_gps[i] = "G0"+self.sv_gps[i]
            else: self.sv_gps[i] = "G"+self.sv_gps[i]

        self.sv_glonass = [str(i+1) for i in range(24)]
        for i,sv in enumerate(self.sv_glonass):
            if i<9: self.sv_glonass[i] = "R0"+self.sv_glonass[i]
            else: self.sv_glonass[i] = "R"+self.sv_glonass[i]

        self.sv_galileo = [str(i+1) for i in range(36)]
        for i,sv in enumerate(self.sv_galileo):
            if i<9: self.sv_galileo[i] = "E0"+self.sv_galileo[i]
            else: self.sv_galileo[i] = "E"+self.sv_galileo[i]
        
        ##if not os.path.isfile(self.f_sat_bias): 
    
        # Test if each file in list of stations exists    
        for f_obs in self.list_f_obs:
            if not os.path.exists(f_obs):
                return       
        
        self.station = ""
        self.coord = ""

        # Estimated high of ionosphere
        self.h = h

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

        
        self.df_sat_DCB = pd.DataFrame()
        
        # Time and year of file
        self.year = None
        self.doy = None
        self.t_min = {}
        self.t_max = {}

        # Bias of the antenna, calculated by method "compute_reveiver_bias"
        # Value stored in csv "stations.csv"
        self.br = {'station':[],'constellation':[],'time_i':[],'time_f':[],'chanel1':[],'chanel2':[],'br':[]}
        
        if not os.path.exists(st.root_dir + "TEC/"):
            os.mkdir(st.root_dir + "TEC/")

    '''

   
    def rinex_to_stec(self,station):
        ''' Extract the relevant data for tec calculation from observation and navigation
            Compute STEC of pseudo range and code phase
        '''       

        
        list_f_obs = self.dict_f_obs[station]


        for f_obs in list_f_obs:
            print (f_obs)
            try:
                rfile = rx.rinex(f_obs)
                header = rfile.read_header()
            except ValueError: continue

            self.station =header["name_station"].replace(" ","").lower()
            self.coord = header['position']
            list_df_f_obs = rfile.read_data()
            #try: list_df_f_obs = rfile.read_data()
            #except:
            #    print ("PROBLEM READING ",f_obs)
            #    continue
        
            for const, df in list_df_f_obs.items():
                #if const not in ['G']: continue
                if const in self.list_df.keys():
                    self.list_df[const] = pd.concat([self.list_df[const],df])
                else:
                    self.list_df[const] = df.copy()

        #const_to_del = []
        #for k in self.list_df.keys():
        #    if k=="S": continue
        #    const_to_del.append(k)
        #for const in const_to_del:
        #    del self.list_df[const]
        
        ### GPS
        if 'G' in self.list_df.keys():
            const = 'G'
            list_cols = self.list_df[const].columns

            chan={}

            if ("C1C" in list_cols) and ("C2W" in list_cols) and ("L1C" in list_cols) and ("L2W" in list_cols):
                C1,C2,L1,L2 = "C1C","C2W","L1C","L2W"
                chan = {"C1":C1,"C2":C2,"L1":L1,"L2":L2}
            elif ("C1" in list_cols) and ("P2" in list_cols) and ("L1" in list_cols) and ("L2" in list_cols):
                C1,C2,L1,L2 = "C1","P2","L1","L2"
                chan = {"C1":"P1","C2":"P2","L1":"L1","L2":"L2"}
            elif ("P1" in list_cols) and ("P2" in list_cols) and ("L1" in list_cols) and ("L2" in list_cols):
                C1,C2,L1,L2 = "P1","P2","L1","L2"
                chan = {"C1":"P1","C2":"P2","L1":"L1","L2":"L2"}

            if chan:

                self.channels[const] = []
                self.channels[const].append(chan)
                self.list_df[const].set_index("time",inplace=True)
                self.t_min[const] = min(self.list_df[const].index)
                self.t_max[const] = max(self.list_df[const].index)        
                
                self.list_df[const]['STEC_l'] = (self.list_df[const][L1]*self.gps_lambda1-self.list_df[const][L2]*self.gps_lambda2)*self.gps_alpha/1e16
                self.list_df[const]['STEC_p'] = (self.list_df[const][C2]-self.list_df[const][C1])*self.gps_alpha/1e16
    
                self.list_df[const]["C1"] = chan["C1"]
                self.list_df[const]["C2"] = chan["C2"]
                
                self.list_df[const] = self.list_df[const][['sv',"C1","C2","STEC_l","STEC_p"]]

            else: del self.list[const]

            
        ### GLONASS
        if 'R' in self.list_df.keys():
            const = 'R'
            list_sv = self.list_df['R']['sv'].unique().tolist()
            if len(list_sv)==0: 
                del self.list_df[const]
                return   

            list_cols = self.list_df[const].columns

            chan = {}
            if ("C1P" in list_cols) and ("C2P" in list_cols) and ("L1P" in list_cols) and ("L2P" in list_cols):
                C1,C2,L1,L2 = "C1P","C2P","L1P","L2P"
                for c in ['P','C']:
                    varc1 = 'C1'+c
                    if varc1 in list_cols:
                        C1 = varc1
                        break
                for c in ['P','C']:
                    varc2 = 'C2'+c
                    if varc2 in list_cols:
                        C2 = varc2
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
                chan = {"C1":C1,"C2":C2,"L1":L1,"L2":L2}
            elif ("L1" in list_cols) and ("L2" in list_cols):
                L1,L2 = "L1","L2"
                if ("P2" in list_cols):
                    C2="P2"
                    if ("P1" in list_cols): C1="P1"
                    elif ("C1" in list_cols): C1="C1"
                elif ("C2" in list_cols):
                    C2="C2"
                    if ("P1" in list_cols): C1="P1"
                    elif ("C1" in list_cols): C1="C1"
                chan = {"C1":"P1","C2":"P2","L1":L1,"L2":L2}

            if chan:
                self.channels[const] = []
                self.channels[const].append(chan)
                
                self.list_df[const].set_index("time",inplace=True)
                #self.t_min[const] = min(self.list_df[const].index)
                #self.t_max[const] = max(self.list_df[const].index)      
    
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

       
                self.list_df[const]["STEC_p"] = (self.list_df[const][C2] - self.list_df[const][C1])*self.list_df[const]["alpha"]/1e16
                self.list_df[const]["STEC_l"] = (self.list_df[const]["lambda1"]*self.list_df[const][L1] - \
                                               self.list_df[const]["lambda2"]*self.list_df[const][L2])*self.list_df[const]["alpha"]/1e16

                
                self.list_df[const].dropna(subset=["STEC_p","STEC_l"],inplace=True)

                
                if len(self.list_df[const])!=0:
                    self.list_df[const]["C1"] = chan["C1"]
                    self.list_df[const]["C2"] = chan["C2"]
                    self.list_df[const] = self.list_df[const][['sv',"C1","C2","STEC_l","STEC_p"]]
                    self.t_min[const] = min(self.list_df[const].index)
                    self.t_max[const] = max(self.list_df[const].index)
                else: 
                    del self.list_df[const]
                    del self.channels[const]
            else: del self.list_df[const]
 
        ## GALILEO
        if "E" in self.list_df.keys():

            const = 'E'
            self.list_df[const] = self.list_df[const].dropna(axis=1, how='all')
            list_cols = self.list_df[const].columns
            
            C1, C2 = '', ''
            L1, L2 = '', ''
            chan = {}
            for c in ['C','X','A','B','Z']:
                varc1 = 'C1'+c
                if varc1 in list_cols:
                    C1 = varc1
                    break
            for c in ['Q','X','I']:
                varc2 = 'C5'+c
                if varc2 in list_cols:
                    C2 = varc2
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
            if (C1!='') and (C2!='') and (L1!='') and (L2!=''): 
                chan = {"C1":C1,"C2":C2,"L1":L1,"L2":L2}

            if chan:
                self.channels[const] = []
                self.channels[const].append(chan)
                
                self.list_df[const].set_index("time",inplace=True)
                self.list_df[const]["STEC_p"] = (self.list_df[const][C2] - self.list_df[const][C1])*self.gps_alpha/1e16
                self.list_df[const]["STEC_l"] = (self.gps_lambda1*self.list_df[const][L1] - self.gps_lambda5*self.list_df[const][L2])*self.gps_alpha/1e16
                self.list_df[const]["C1"] = chan["C1"]
                self.list_df[const]["C2"] = chan["C2"]
                self.list_df[const] = self.list_df[const][['sv',"C1","C2","STEC_l","STEC_p"]]
                #self.list_df['E'].dropna(inplace=True)
                self.list_df[const].dropna(subset=["STEC_l"],inplace=True)
                
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
           
            beidu_columns = self.list_df["C"].columns
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
                df_beidu_2_6["STEC_l"] = (df_beidu_2_6[L1]*bds_lambda1-df_beidu_2_6[L2]*bds_lambda3)*beidu_alpha_2
                df_beidu_2_6["STEC_p"] = (df_beidu_2_6[C2]-df_beidu_2_6[C1])*beidu_alpha_2
                df_beidu_2_6.dropna(subset="STEC_l",inplace=True)
                if len(df_beidu_2_6)!=0: 
                    self.channels[const].append(chan)

            self.list_df[const] = pd.DataFrame()
            if len(df_beidu_2_6)!=0: self.list_df[const] = pd.concat([self.list_df[const],df_beidu_2_6])
            if len(df_beidu_2_7)!=0: self.list_df[const] = pd.concat([self.list_df[const],df_beidu_2_7])


            if len(self.list_df[const])!=0:
                self.list_df[const].set_index("time",inplace=True)
                self.list_df[const] = self.list_df[const][['sv',"C1","C2","STEC_l","STEC_p"]]
                self.list_df[const].dropna(subset="STEC_l",inplace=True)
            
                self.t_min[const] = min(self.list_df[const].index)
                self.t_max[const] = max(self.list_df[const].index)
            else:
                del self.list_df[const]
        #### QZSS
        if "J" in self.list_df.keys():
            const = 'J'
            
            self.list_df[const] = self.list_df[const].dropna(axis=1, how='all')

            C1, C2, L1, L2 = '', '', '', ''
            for c in ['C','X','S','L','Z']:
                varc1 = 'C1'+c
                if varc1 in self.list_df[const].columns:
                    C1 = varc1
                    break
            for c in ['Q','X','I']:
                varc2 = 'C5'+c
                if varc2 in self.list_df[const].columns:
                    C2 = varc2
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
            if (C1!='') and (C2!='') and (L1!='') and (L2!=''): 
                chan = {"C1":C1,"C2":C2,"L1":L1,"L2":L2}
                    
            chan = {"C1":C1,"C2":C2,"L1":L1,"L2":L2}

            if chan:
                self.channels[const] = []
                self.channels[const].append(chan)
                
                qzss_alpha = self.gps_f1**2*self.gps_f5**2/(self.gps_f1**2-self.gps_f5**2)/40.318
                self.list_df[const].set_index("time",inplace=True)
                self.list_df[const]["STEC_l"] = (self.list_df[const][L1]*self.gps_lambda1-self.list_df[const][L2]*self.gps_lambda5)*qzss_alpha/1e16
                self.list_df[const]["STEC_p"] = (self.list_df[const][C2]-self.list_df[const][C1])*qzss_alpha/1e16
                self.list_df[const]["C1"] = chan["C1"]
                self.list_df[const]["C2"] = chan["C2"]
                self.list_df[const] = self.list_df[const][['sv',"C1","C2","STEC_l","STEC_p"]]

                self.t_min[const] = min(self.list_df[const].index)
                self.t_max[const] = max(self.list_df[const].index)

            else:
                del self.list_df[const]



        if "S" in self.list_df.keys():

            const = "S"
            C1, C2, L1, L2 = 'C1C', 'C5I', 'L1C', 'L5I'

            chan = {"C1":C1,"C2":C2,"L1":L1,"L2":L2}

            if chan:

                self.channels[const] = []
                self.channels[const].append(chan)
                self.list_df[const].set_index("time",inplace=True)
                self.t_min[const] = min(self.list_df[const].index)
                self.t_max[const] = max(self.list_df[const].index)        

                sbas_alpha = self.gps_f1**2*self.gps_f5**2/(self.gps_f1**2-self.gps_f5**2)/40.318
                self.list_df[const]['STEC_l'] = (self.list_df[const][L1]*self.gps_lambda1-self.list_df[const][L2]*self.gps_lambda5)*sbas_alpha/1e16
                self.list_df[const]['STEC_p'] = (self.list_df[const][C2]-self.list_df[const][C1])*sbas_alpha/1e16
    
                self.list_df[const]["C1"] = chan["C1"]
                self.list_df[const]["C2"] = chan["C2"]
                
                self.list_df[const] = self.list_df[const][['sv',"C1","C2","STEC_l","STEC_p"]]

            else: del self.list_df[const]

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
            self.gnss.load_sats(const,self.t_min[const],self.t_max[const])            
            self.list_df[const] = self.gnss.getElevation(self.list_df[const],self.coord)
            self.list_df[const].dropna(subset="elevation",inplace=True)
            self.list_df[const] = self.gnss.getPiercingPoint(self.list_df[const],self.coord,self.h)            
            self.list_df[const].dropna(subset="elevation",inplace=True)

            self.list_df[const] = self.list_df[const][["sv","C1","C2","elevation","lat","lon","alt","STEC_l","STEC_p"]]
                        
        for const in const_without_pos:
            del self.list_df[const]
            del self.channels[const]
        return True


    def list_leaps_series(self,series,tol_dev=0.2,tol_sig=10,N_=5,debug=False):
        ''' detects leaps in the STEC_sl time series
            input: Data pandas Series indexes with time
            output: list of dictionnaries
                Each dictionary contains first and last date of time series without leap
        '''

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

        s=0

        list_fit_params = {}

        unstable_left=True

        list_series = series.tolist()
        
        while s<len(list_series)-N:

            #Going for N next point without big jump
            if False:
                has_big_time_jump=True
                while has_big_time_jump and s<len(series)-N:
                    for i in range(N-1):
                        if diffs[s+1+i]-diffs[s+i]>2*60:
                            print ("TIME JUMP")
                            if not border["left"] is None:
                                border["right"]=s+i
                                list_borders.append(border)
                                border={"left":None,"right":None,"left_A":None,"left_B":None,"right_A":None,"right_B":None}
                                unstable_left=True
                            s=s+1
                            break
                        has_big_time_jump=False

            #Get linear fit parameter of N right points
            A,B,max_dev,mean_dev = fit_lin(diffs[s:s+N],series[s:s+N].values)
            # Compute distance of point s and s+N+1 with fit
            left_dev=abs(list_series[s-1]-A*diffs[s-1]-B) if s>0 else None

            right_dev=abs(list_series[s+N]-A*diffs[s+N]-B) if s+N<len(series) else None
            list_fit_params[s]={"A":A,"B":B,"max_dev":max_dev,"left_dev":left_dev,"right_dev":right_dev}
            s+=1


        border={"left":None,"t_left":None,"right":None,"t_right":None,"left_A":None,"left_B":None,"right_A":None,"right_B":None}

        retroactive_borders = []

        s1=0
        s2=0
        while s1<len(series):
            # Looking for an right border
            if border["right"] is None:
                if s1>=len(series)-N+1:
                    if debug:
                        print ("RIGHT",s1,series.index[s1],len(series),list_fit_params[s1]["left_dev"],list_fit_params[s1]["right_dev"],tol_dev)
                        plot_leap(diffs,series,s1-N,list_fit_params[s1-N]["A"],list_fit_params[s1-N]["B"],N,list_borders,"forward")
                    list_borders = filter_slope_leap(list_borders,series,diffs)
                    return list_borders
                    #return list_borders
                if s1==len(series)-N:
                    border["right"]=s1+N-1
                    border["t_right"]=series.index[s1+N-1]
                    border["right_A"]=list_fit_params[s1-1]["A"]
                    border["right_B"]=list_fit_params[s1-1]["B"]

                    s2=s1-1
                elif list_fit_params[s1]["right_dev"]>tol_deviation or (list_fit_params[s1]["right_dev"]>tol_sig*list_fit_params[s1]["max_dev"] and list_fit_params[s1]["right_dev"]>0.04*self.resolution):
                    if debug:
                        if s1>=0:
                            print ("RIGHT",s1,series.index[s1],len(series),list_fit_params[s1]["left_dev"],list_fit_params[s1]["right_dev"],tol_dev)
                            plot_leap(diffs,series,s1,list_fit_params[s1]["A"],list_fit_params[s1]["B"],N,list_borders)
                    border["right"]=s1+N-1
                    border["t_right"]=series.index[s1+N-1]
                    if s1!=0:
                        border["right_A"]=list_fit_params[s1]["A"]
                        border["right_B"]=list_fit_params[s1]["B"]
                    s2=s1
                else:
                    if debug:
                        if s1>=0: print ("RIGHT",s1-N,series.index[s1],len(series),list_fit_params[s1]["left_dev"],list_fit_params[s1]["right_dev"],tol_dev)
                    s1+=1
            else:

                if debug:
                    print ("LEFT",s2,series.index[s2],len(series),list_fit_params[s2]["left_dev"],list_fit_params[s2]["right_dev"],tol_dev)
                    plot_leap(diffs,series,s2,list_fit_params[s2]["A"],list_fit_params[s2]["B"],N,list_borders)
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
                                list_borders.append(b)
                            retroactive_borders=[]
                            #if debug: plot_leap(diffs,series,s2,list_fit_params[s2]["A"],list_fit_params[s2]["B"],N,list_borders,"forward")

                        border={"left":None,"t_left":None,"right":None,"t_right":None,"left_A":None,"left_B":None,"right_A":None,"right_B":None}
                        s1=s1+N
                        continue
                    border["left"]=s2
                    border["t_left"]=series.index[s2]
                    border["left_A"]=list_fit_params[s2]["A"]
                    border["left_B"]=list_fit_params[s2]["B"]
                    retroactive_borders.append(border)
                    # If right border of last segment is not reached, must keep
                    # searching backward for a new segment
                    if len(list_borders)>0 and border["left"]-list_borders[-1]["right"]>N:
                        if s2-1>=N:
                            border={"left":None,"t_left":None,"right":s2-1,"t_right":series.index[s2-1],"left_A":None,"left_B":None,"right_A":list_fit_params[s2-N]["A"],"right_B":list_fit_params[s2-N]["B"]}
                        else:
                            border={"left":None,"t_left":None,"right":s2-1,"t_right":series.index[s2-1],"left_A":None,"left_B":None,"right_A":None,"right_B":None}
                        s2-=N
                        continue
                    retroactive_borders.reverse()
                    for b in retroactive_borders:
                        list_borders.append(b)
                    retroactive_borders=[]

                    if s1>=len(series)-N-1:

                        if debug:
                            plot_leap(diffs,series,s2,list_fit_params[s2]["A"],list_fit_params[s2]["B"],N,list_borders)
                        list_borders = filter_slope_leap(list_borders,series,diffs)
                        return list_borders
                    border={"left":None,"t_left":None,"right":None,"t_right":None,"left_A":None,"left_B":None,"right_A":None,"right_B":None}
                    s1=s1+N
                else: s2-=1
        list_borders = filter_slope_leap(list_borders,series,diffs)
        return list_borders

    def list_leaps(self,df):
        borders=[]
        stec_l = ''
        stec_p = ''
        for c in df.columns:
            if "STEC_l" in c: 
                stec_l=c
            if "STEC_p" in c: 
                stec_p=c
                #break
        borders_sl = self.list_leaps_series(df[stec_l],0.4,3,3,debug=False)
        borders_slp = self.list_leaps_series(df[stec_p],100,150,5,debug=False)

        # We look at the borders found for slp and sll leaps, and make the
        # match
        # For each border found for sll
        if len(borders_sl)==0: return borders
        for bsll in borders_sl:
            # For each border found for slp
            for bslp in borders_slp:
                # Test if left border of sll in inside slp borders item
                if bsll["t_left"]>=bslp["t_left"] and bsll["t_left"]<bslp["t_right"]:
                    # Left of new border must be the one of bsll
                    newb = bsll.copy()
                    newb["t_left"]=bsll["t_left"]
                    # We asign to the new border the earliest right border
                    if bsll["t_right"]<bslp["t_right"]: newb["t_right"]=bsll["t_right"]
                    else: newb["t_right"]=bslp["t_right"]
                    borders.append(newb)
                    break
                # If last condition not fulfilled, test in the right sll border
                # is inside the slp border item
                elif bsll["t_right"]>bslp["t_left"] and bsll["t_right"]<=bslp["t_left"]:
                    newb = bsll.copy()
                    newb["t_right"]=bsll["t_right"]
                    #Left border must be outside blsp, otherwise it would
                    #have been seen by the previous condition
                    newb["t_left"]=bslp["t_left"]
                    borders.append(newb)
                    break
                # Case the bsll fully contains bslp, take bslp borders
                elif bsll["t_left"]<bslp["t_left"] and bsll["t_right"]>bslp["t_right"]:
                    borders.append(bslp)
                    break

        return borders

    def add_baseline(self):

        for const, df_data in self.list_df.items():

            df_data = df_data.dropna(subset=["elevation"])
            
            if len(df_data)==0: continue
            

            t_begin = df_data.index[0]
            t_end = df_data.index[-1]

            df_filtered = pd.DataFrame()

            list_sv = df_data["sv"].unique().tolist()

            for sat in list_sv:

                df_sat = df_data[df_data["sv"]==sat]
                df_sat_filter = pd.DataFrame()

                for channel in self.channels[const]:

                    
                    chan_filter = (df_sat["C1"]==channel["C1"]) & (df_sat["C2"]==channel["C2"])
                    df_sat_chanel = df_sat[chan_filter]

                        
                    if len(df_sat_chanel["STEC_l"].dropna())==0: continue
                    

                    # Make list of arcs of satellite using elevation information
                    elevations = df_sat_chanel["elevation"]
                    list_arcs = gnss.get_arcs(elevations,t_begin,t_end)

                    # Squared of sin of elevation for baseline (Brs) pondering
                    df_sat_chanel["sin2_ele"] = np.sin(df_sat_chanel['elevation'])**2

                    # Individual values of the sum for future baseline calculation
                    df_sat_chanel['BRs'] = (df_sat_chanel["STEC_p"]-df_sat_chanel["STEC_l"])*df_sat_chanel["sin2_ele"]
                    
                    # cos(chi) to calculate VTEC from STEC
                    df_sat_chanel['cos_chi'] = np.cos(np.arcsin(R_E*np.cos(df_sat_chanel["elevation"])/(R_E+self.h)))
                    n_arc=0

                    df_arcs = pd.DataFrame()
                    for arc in list_arcs:
                      
                        # Must gets borders that are inside the arc segment
                        arc_borders = None
                        arc_borders = []

                        # Extract data from satellite arc
                        date_filter = (df_sat_chanel.index>=arc["start"]) & (df_sat_chanel.index<=arc["end"])
                        df_arc = df_sat_chanel[date_filter]

                        if len(df_arc)<=8:continue                   
                        arc_borders = self.list_leaps(df_arc)

                        
                        
                        if len(arc_borders)==0: continue
        
                        # Eliminate segments with too low elevation satellites
                        list_max_ele = []
                        i=0
                        while i<len(arc_borders):
                            t_dist=(arc_borders[i]["t_right"]-arc_borders[i]["t_left"]).seconds
                            date_border_filter  = (df_arc.index>=arc_borders[i]["t_left"]) & (df_arc.index<=arc_borders[i]["t_right"])
                            max_ele = max(df_arc[date_border_filter]["elevation"])
                            #max_ele = max(df_arc.loc[arc_borders[i]["t_left"]:arc_borders[i]["t_right"]]["elevation"])
                            if max_ele<30*np.pi/180 and t_dist<20*60:
                                arc_borders.pop(i)
                            else:
                                i+=1
                                list_max_ele.append(max_ele)

                        #### Clasify segments #####
                        d_df_seg = []
                        has_sane_parts = False
                        sane_indices = []
                        list_size_segs = []

                        for i,b in enumerate(arc_borders):
                            size_seg = (arc_borders[i]["t_right"]-arc_borders[i]["t_left"]).seconds
                            list_size_segs.append(size_seg)

        
                            if size_seg>50*60:
                                df_seg = None
                                date_border_filter  = (df_arc.index>=arc_borders[i]["t_left"]) & (df_arc.index<=arc_borders[i]["t_right"])
                                df_seg = df_arc[date_border_filter]

                                df_br = df_seg.dropna(subset=['BRs',"sin2_ele"])


                                brs=df_br['BRs'].sum(skipna=True)/df_br["sin2_ele"].sum(skipna=True)
                                df_seg["STEC_l"] = df_seg["STEC_l"] + brs
                                d_df_seg.append(df_seg)
                                has_sane_parts = True
                                sane_indices.append(i)
                            else:
                                d_df_seg.append(None)
                        if not has_sane_parts: 
                            continue

                        # paste left part of arc
                        for i in range(sane_indices[0]-1,-1,-1):
                            df_seg = None
                            date_border_filter  = (df_arc.index>=arc_borders[i]["t_left"]) & (df_arc.index<=arc_borders[i]["t_right"])
                            df_seg = df_arc[date_border_filter]

                            t_dist=(arc_borders[i+1]["t_left"]-arc_borders[i]["t_right"]).seconds
                            slope = (arc_borders[i+1]["left_A"]+arc_borders[i]["right_A"])/2
        
                            glue=d_df_seg[i+1]["STEC_l"].iloc[0]-df_seg["STEC_l"].iloc[-1]-slope*t_dist
                            df_br = df_seg.dropna(subset=['BRs',"sin2_ele"])
                            brs=df_br['BRs'].sum(skipna=True)/df_br["sin2_ele"].sum(skipna=True)
        
                            if abs(glue-brs)>20: correction = brs
                            else: correction = glue
                            df_seg["STEC_l"] = df_seg["STEC_l"] + correction
                            d_df_seg[i] = df_seg

                        # paste right part of arc
                        for i in range(sane_indices[-1]+1,len(arc_borders)):
                            df_seg = None
                            date_border_filter  = (df_arc.index>=arc_borders[i]["t_left"]) & (df_arc.index<=arc_borders[i]["t_right"])
                            df_seg = df_arc[date_border_filter]
                            t_dist=(arc_borders[i]["t_left"]-arc_borders[i-1]["t_right"]).seconds
        
                            slope = arc_borders[i-1]["right_A"]
        
                            glue=d_df_seg[i-1]["STEC_l"].iloc[-1]-df_seg["STEC_l"].iloc[0]+slope*t_dist
                            df_br = df_seg.dropna(subset=['BRs',"sin2_ele"])
                            brs=df_br['BRs'].sum(skipna=True)/df_br["sin2_ele"].sum(skipna=True)
                            if abs(glue-brs)>20: correction = brs
                            else: correction = glue
                            df_seg["STEC_l"] = df_seg["STEC_l"] + correction
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
                                        t_left = (arc_borders[ileft]["t_left"]-arc_borders[ileft-1]["t_right"]).seconds
                                        t_right = (arc_borders[iright]["t_left"]-arc_borders[ileft]["t_right"]).seconds
                                        if t_left<=t_right or only_right:
                                            df_seg = None
                                            date_border_filter  = (df_arc.index>=arc_borders[ileft]["t_left"]) & (df_arc.index<=arc_borders[ileft]["t_right"])
                                            df_seg = df_arc[date_border_filter]
                                            #df_seg = df_arc.loc[arc_borders[ileft]["t_left"]:arc_borders[ileft]["t_right"]]
                                            t_dist=(arc_borders[ileft]["t_left"]-arc_borders[ileft-1]["t_right"]).seconds
        
                                            slope = arc_borders[ileft-1]["right_A"]
        
                                            glue=d_df_seg[ileft-1]["STEC_l"].iloc[-1]-df_seg["STEC_l"].iloc[0]+slope*t_dist
                                            df_br = df_seg.dropna(subset=['BRs',"sin2_ele"])
                                            brs=df_br['BRs'].sum(skipna=True)/df_br["sin2_ele"].sum(skipna=True)
                                            if abs(glue-brs)>20: correction = brs
                                            else: correction = glue
                                            df_seg["STEC_l"] = df_seg["STEC_l"] + correction
        
                                            d_df_seg[ileft] = df_seg
                                            ileft+=1
                                        elif not only_right:
                                            only_left=True
                                        if ileft>iright: break
                                    if not only_right:
                                        t_left = (arc_borders[iright]["t_left"]-arc_borders[ileft]["t_right"]).seconds
                                        t_right = (arc_borders[iright+1]["t_left"]-arc_borders[iright]["t_right"]).seconds
                                        if t_left>=t_right or only_left:
                                            df_seg = None
                                            date_border_filter  = (df_arc.index>=arc_borders[iright]["t_left"]) & (df_arc.index<=arc_borders[iright]["t_right"])
                                            df_seg = df_arc[date_border_filter]
                                            #df_seg = df_arc.loc[arc_borders[iright]["t_left"]:arc_borders[iright]["t_right"]]
                                            t_dist=(arc_borders[iright+1]["t_left"]-arc_borders[iright]["t_right"]).seconds
        
                                            slope = arc_borders[iright+1]["right_A"]
        
                                            glue=d_df_seg[iright+1]["STEC_l"].iloc[0]-df_seg["STEC_l"].iloc[-1]+slope*t_dist
                                            df_br = df_seg.dropna(subset=['BRs',"sin2_ele"])
                                            brs=df_br['BRs'].sum(skipna=True)/df_br["sin2_ele"].sum(skipna=True)
                                            if abs(glue-brs)>20: correction = brs
                                            else: correction = glue
                                            df_seg["STEC_l"] = df_seg["STEC_l"] + correction
        
                                            d_df_seg[iright] = df_seg
                                            iright-=1
                                        elif not only_left:
                                            only_right=True

                        df_filter_arc = pd.DataFrame()
                        for df in d_df_seg:
                            if df is None:continue
                            df_filter_arc = pd.concat([df_filter_arc,df])                        
                        
                        
                        df_arcs = pd.concat([df_arcs,df_filter_arc])    
        
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
                self.list_df[const]["STEC_l"] += self.list_df[const]["dcb"]
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

        

    def load_sat_DCB(self):

        if len(self.list_f_dcb)==0: return pd.DataFrame()
        else:
            self.sat_dcb = DCB.load_dcb(self.list_f_dcb)
            return
        

    def to_feather(self):

        df_obs = pd.DataFrame()

        for const in self.list_df.keys():
            df_obs = pd.concat([
                df_obs,
                self.list_df[const]#[["time","sv","lat","lon","elevation","STEC_l","VTEC"]]
            ])
            

        feather_path = st.root_dir + "TEC/" + self.station
        df_obs.to_csv(feather_path+".csv")
        
    
    def compute_vtec(self,station):	

        self.list_df = {}
        self.channels = {}
     
        self.rinex_to_stec(station)
        
        self.add_satellite_pos()
        
        print ("Calculating baseline to correct Slant TEC")
        self.add_baseline()

        print ("Calculating receiver bias, correct Slant TEC, compute VTEC")
        self.add_receiver_bias()
        
        #self.to_feather(st.root_dir + "TEC/" + str(self.year) + "/" + self.station + ".feather" )
        self.to_feather()
    
