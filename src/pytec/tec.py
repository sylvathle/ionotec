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
#import matplotlib.pyplot as plt
#import matplotlib.dates as md
import pymap3d as pm
#from os import listdir,path,mkdir
#from os.path import isfile, join
import julian
import math	
import time

#from .stations import root_dir
from . import stations as st
from . import gnss



pd.options.mode.chained_assignment = None

# Earth Radius
R_E = 6371000

## Definition of the constants for the equations of STEC
f1,f2 = 1575.42 * 1e6, 1227.60 * 1e6
lambda1,lambda2=csts.c/f1,csts.c/f2
alpha = f1**2*f2**2/(f1**2-f2**2)/40.318

# List of possible indices in observation RINEX files, P1 not available
#full_vars_list = ['P1','P2','C1','C2','L1','L2','S1','S2']
full_vars_list = ['P2','C1','C2','L1','L2','S1','S2']

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

#def filter_overlap_leap(list_borders,series):
#    return list_borders

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
    print ("plot_leap")
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
    #if deviation>tol_deviation: print ("deviation<tol_deviation",deviation,tol_deviation)
    #if deviation>10*max_deviation: print ("deviation>10*max_deviation")
    plt.title(title)
    mng = plt.get_current_fig_manager()
    mng.full_screen_toggle()

    #plt.savefig("cycle_slips_png/almaG01_"+str(n_frame)+".png")
    plt.show()
    plt.close()

from pathlib import Path

class tec:

    def __init__(self,list_f_obs,f_nav,f_sat_bias,outfolder="output",resolution=60,h=400000):
    
        # List of observation files
        self.list_f_obs = list_f_obs        
        # List of navigation files
        self.f_nav = f_nav
        # File containing satellite bias
        self.f_sat_bias = f_sat_bias
        #Output feather file
        #if f_out == "": 
        #    self.f_feather = self.list_f_obs[min(len(self.list_f_obs)-1,1)][:-4]+"tec.feather"
        #else:
        #    if f_out[-8:] ==".feather": self.f_feather = f_out
        #    else: self.f_feather = f_out + ".feather"
        
        if not os.path.isfile(self.f_sat_bias): 
            print ("bias file location not found, will take 0, might affect strongly the results")
    
        # Test if each file in list of stations exists    
        for f_obs in self.list_f_obs:
            if not os.path.exists(f_obs):
                print ("File",f_obs,"does not exist, nothing to process")
                return       
        
        self.station = ""
        self.coord = ""

        # Estimated high of ionosphere
        self.h = h

        #Resolution that will be used for the process (should by 60 seconds)
        self.resolution = 60

        # Resolution of the rinex station, in seconds.
        self.obs_resolution = ""

        # DataFrame containing observation data
        self.df_obs = pd.DataFrame()

        # List of satellites that have been observed
        self.sv = []

        # List of dictionnaries containing borders of signal for each satellite
        self.borders={}

        # Time and year of file
        self.year = None
        self.doy = None

        # Bias of the antenna, calculated by method "compute_reveiver_bias"
        # Value stored in csv "stations.csv"
        self.br = 0
        
        if not os.path.exists(st.root_dir + "TEC/"):
            try: os.mkdir(st.root_dir + "TEC/")
            except OSError as e:
                    if e.errno!=17: print ("FAIL creation of directory "+st.root_dir + "TEC/", e )
            else: print ("Successfully created the directory "+st.root_dir + "TEC/")

    def rinex_to_stec(self):
        ''' Extract the relevant data for tec calculation from observation and navigation
            Compute STEC of pseudo range and code phase
        '''

        # Header for observation
        try: 
            # Header of the central file, the one we want to save the information
            header = gr.rinexheader(self.list_f_obs[min(len(self.list_f_obs)-1,1)])
            # Header of the first file, to extract 
        except ValueError: return False
        self.station =header["MARKER NAME"].replace(" ","").lower()
        self.coord = header['position']
        if "INTERVAL" in header.keys():
            self.obs_resolution = float(header["INTERVAL"].replace(" ",""))
        else:
            self.obs_resolution = 1.0

        t_obs_splt = header["TIME OF FIRST OBS"].split(" ")
        t_obs_list = []
        for t in t_obs_splt:
            if t not in ["","GPS"]:
                t_obs_list.append(float(t))

        date_observation = datetime.datetime(int(t_obs_list[0]),int(t_obs_list[1]),int(t_obs_list[2]),int(t_obs_list[3]),int(t_obs_list[4]),int(t_obs_list[5]))
        self.year = date_observation.year
        self.doy = str(date_observation.timetuple().tm_yday)
        #sys.exit()

        list_vars = []

        # lista de campos accesibles en el header
        fields = header['fields']

        # Revisamos si las variables están disponibles
        for l in full_vars_list:
            if l in fields: list_vars.append(l)

        # Concat all data in one file
        first=True
        self.df_obs = pd.DataFrame()
        for f_obs in self.list_f_obs:
            # Read data of observation file
            obs=gr.load(f_obs,meas=list_vars,use="G")
 
            if first: 
                # Rinex a DataFrame
                self.df_obs = obs.to_dataframe()
                self.df_obs.index.set_names(["time","sv"],inplace=True)
                self.df_obs.reset_index(level=["sv"],inplace=True)
                self.df_obs.index = pd.to_datetime(self.df_obs.index)
                first=False
            else:
                # Rinex a DataFrame
                df_day = obs.to_dataframe()
                df_day.index.set_names(["time","sv"],inplace=True)
                df_day.reset_index(level=["sv"],inplace=True)
                df_day.index = pd.to_datetime(df_day.index)
                self.df_obs = pd.concat([df_day,self.df_obs])
            
            

        # List satellites seen by the station and prepare dict of list_borders
        for s in self.df_obs["sv"].values:
            if s not in self.sv:
                self.sv.append(s)
                self.borders[s]=[]

        # Process Slant TEC for pseudo range (slp)
        self.df_obs['STEC_slp'] = (self.df_obs['P2'] - self.df_obs['C1'])*alpha/1e16
        # Process Slant TEC for code phase (sll)
        self.df_obs['STEC_sll'] = (lambda1*self.df_obs['L1'] - lambda2*self.df_obs['L2'])*alpha/1e16

        # Changing resolution
        self.df_obs = self.df_obs.groupby("sv").resample(str(self.resolution)+"S").mean()
        self.df_obs.reset_index(level=["sv"],inplace=True)

        # Remove nan values (may be optimized)
        self.df_obs.dropna(subset=["STEC_slp","STEC_sll"],inplace=True)

        ## When intermediate feather not used anymore .reset_index() should be removed
        self.df_obs = self.df_obs[["sv","STEC_slp","STEC_sll"]]#.reset_index()
        return True


    def add_satellite_pos(self):
        #time_list =  self.df_obs.index
        gps = gnss.gnss(self.f_nav)
        
        d_in, d_out = min(self.df_obs.index), max(self.df_obs.index)
        gps.load_sats(d_in,d_out)
        df_data = pd.DataFrame()

        t_tot = 0
        t_full = 0
        tf = time.time()
        
        
        self.df_obs = gps.getElevation(self.df_obs,self.coord)
        self.df_obs = gps.getPiercingPoint(self.df_obs,self.coord,self.h)
        
        return True


    def list_leaps_series(self,series,tol_dev=0.2,tol_sig=10,N_=5,debug=False):
        ''' detects leaps in the STEC_sl time series
            input: Data pandas Series indexes with time
            output: list of dictionnaries
                Each dictionary contains first and last date of time series without leap
        '''
        #series = self.df_obs[self.df_obs["sv"]==sat][signal].loc[start:end]

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

        while s<len(series)-N:

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
            left_dev=abs(series[s-1]-A*diffs[s-1]-B) if s>0 else None
            right_dev=abs(series[s+N]-A*diffs[s+N]-B) if s+N<len(series) else None
            list_fit_params[s]={"A":A,"B":B,"max_dev":max_dev,"left_dev":left_dev,"right_dev":right_dev}
            s+=1


        border={"left":None,"t_left":None,"right":None,"t_right":None,"left_A":None,"left_B":None,"right_A":None,"right_B":None}
        #print (list_fit_params)

        retroactive_borders = []

        s1=0
        s2=0
        while s1<len(series):
            # Looking for an right border
            if border["right"] is None:
                if s1>=len(series)-N+1:
                    #print ("plot 0")
                    #list_borders = remove_overlaps(list_borders)
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
                    #print ("plot 2")
                    s2=s1
                else:
                    if debug:
                        if s1>=0: print ("RIGHT",s1-N,series.index[s1],len(series),list_fit_params[s1]["left_dev"],list_fit_params[s1]["right_dev"],tol_dev)
                    s1+=1
            else:
                #print (s1,border)
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
                    #print (border["right"]-s2)
                    if border["right"]-s2<=N:
                        if len(retroactive_borders)!=0:
                            for b in retroactive_borders:
                                #print ("Adding border to list_borders 1",b)
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
                    #print ("new", border)
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
                    #print ("plot 3")
                    #if debug: plot_leap(diffs,series,s2,list_fit_params[s2]["A"],list_fit_params[s2]["B"],N,list_borders,"forward")

                    if s1>=len(series)-N-1:
                        #print (list_borders)
                        #list_borders = remove_overlaps(list_borders)
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
        borders_sl = self.list_leaps_series(df["STEC_sll"],0.4,3,3,debug=False)
        borders_slp = self.list_leaps_series(df["STEC_slp"],100,150,5,debug=False)

        # We look at the borders found for slp and sll leaps, and make the
        # match
        # For each border found for sll
        if len(borders_sl)==0 or len(borders_slp)==0: return borders
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

        df_data = self.df_obs.copy()
        if len(df_data)==0: return

        t_begin = df_data.index[0]
        t_end = df_data.index[-1]

        df_filtered = pd.DataFrame()

        for sat in self.sv:

            df_sat = df_data[df_data["sv"]==sat]
            if len(df_sat)==0: continue

            # Make list of arcs of satellite using elevation information
            elevations = df_sat["elevation"]
            list_arcs = gnss.get_arcs(elevations,t_begin,t_end)

            # Get satellite bias and correct STEC_sl with it
            sat_bias = gnss.getBias_fromfile(sat,self.f_sat_bias) * alpha * csts.c * 1e-9  / 1e16
            df_sat['STEC_sl']=df_sat['STEC_sll']

            # Squared of sin of elevation for baseline (Brs) pondering
            df_sat["sin2_ele"] = np.sin(df_sat['elevation'])**2

            # Individual values of the sum for future baseline calculation
            df_sat['BRs'] = (df_sat['STEC_slp']-df_sat['STEC_sll'])*df_sat["sin2_ele"]

            # cos(chi) to calculate VTEC from STEC
            df_sat['cos_chi'] = np.cos(np.arcsin(R_E*np.cos(df_sat["elevation"])/(R_E+self.h)))

            n_arc=0
            # Run over each arc
            for arc in list_arcs:
                # Must gets borders that are inside the arc segment
                arc_borders = None
                arc_borders = []

                # Extract data from satellite arc
                df_arc = df_sat.loc[arc["start"]:arc["end"]]
                if len(df_arc)<=8:continue

                arc_borders = self.list_leaps(df_arc)

                if len(arc_borders)==0: continue

                # Eliminate segments with too low elevation satellites
                list_max_ele = []
                #for i,b in enumerate(arc_borders):
                i=0
                while i<len(arc_borders):
                    t_dist=(arc_borders[i]["t_right"]-arc_borders[i]["t_left"]).seconds
                    max_ele = max(df_arc.loc[arc_borders[i]["t_left"]:arc_borders[i]["t_right"]]["elevation"])
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
                        df_seg = df_arc.loc[b["t_left"]:b["t_right"]]
                        df_br = df_seg.dropna(subset=["BRs","sin2_ele"])
                        brs=df_br['BRs'].sum(skipna=True)/df_br["sin2_ele"].sum(skipna=True)
                        df_seg["STEC_sl"] = df_seg["STEC_sl"] + brs
                        d_df_seg.append(df_seg)
                        has_sane_parts = True
                        sane_indices.append(i)
                    else:
                        d_df_seg.append(None)
                if not has_sane_parts: continue

                # paste left part of arc
                for i in range(sane_indices[0]-1,-1,-1):
                    df_seg = None
                    df_seg = df_arc.loc[arc_borders[i]["t_left"]:arc_borders[i]["t_right"]]
                    t_dist=(arc_borders[i+1]["t_left"]-arc_borders[i]["t_right"]).seconds
                    slope = (arc_borders[i+1]["left_A"]+arc_borders[i]["right_A"])/2

                    #df_filter_arc = df_arc.loc[arc_borders[i]["t_left"]:arc_borders[i]["t_right"]]
                    glue=d_df_seg[i+1]["STEC_sl"].iloc[0]-df_seg["STEC_sl"].iloc[-1]-slope*t_dist
                    df_br = df_seg.dropna(subset=["BRs","sin2_ele"])
                    brs=df_br['BRs'].sum(skipna=True)/df_br["sin2_ele"].sum(skipna=True)

                    if abs(glue-brs)>20: correction = brs
                    else: correction = glue
                    df_seg["STEC_sl"] = df_seg["STEC_sl"] + correction
                    d_df_seg[i] = df_seg

                # paste right part of arc
                for i in range(sane_indices[-1]+1,len(arc_borders)):
                    df_seg = None
                    df_seg = df_arc.loc[arc_borders[i]["t_left"]:arc_borders[i]["t_right"]]
                    t_dist=(arc_borders[i]["t_left"]-arc_borders[i-1]["t_right"]).seconds

                    slope = arc_borders[i-1]["right_A"]

                    glue=d_df_seg[i-1]["STEC_sl"].iloc[-1]-df_seg["STEC_sl"].iloc[0]+slope*t_dist
                    df_br = df_seg.dropna(subset=["BRs","sin2_ele"])
                    brs=df_br['BRs'].sum(skipna=True)/df_br["sin2_ele"].sum(skipna=True)
                    if abs(glue-brs)>20: correction = brs
                    else: correction = glue
                    df_seg["STEC_sl"] = df_seg["STEC_sl"] + correction
                    d_df_seg[i] = df_seg


                # Look at intermediate segments
                i = 0
                if len(sane_indices)>=2:
                #if False:
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
                                    df_seg = df_arc.loc[arc_borders[ileft]["t_left"]:arc_borders[ileft]["t_right"]]
                                    t_dist=(arc_borders[ileft]["t_left"]-arc_borders[ileft-1]["t_right"]).seconds

                                    slope = arc_borders[ileft-1]["right_A"]

                                    glue=d_df_seg[ileft-1]["STEC_sl"].iloc[-1]-df_seg["STEC_sl"].iloc[0]+slope*t_dist
                                    df_br = df_seg.dropna(subset=["BRs","sin2_ele"])
                                    brs=df_br['BRs'].sum(skipna=True)/df_br["sin2_ele"].sum(skipna=True)
                                    if abs(glue-brs)>20: correction = brs
                                    else: correction = glue
                                    df_seg["STEC_sl"] = df_seg["STEC_sl"] + correction

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
                                    df_seg = df_arc.loc[arc_borders[iright]["t_left"]:arc_borders[iright]["t_right"]]
                                    t_dist=(arc_borders[iright+1]["t_left"]-arc_borders[iright]["t_right"]).seconds

                                    slope = arc_borders[iright+1]["right_A"]

                                    glue=d_df_seg[iright+1]["STEC_sl"].iloc[0]-df_seg["STEC_sl"].iloc[-1]+slope*t_dist
                                    df_br = df_seg.dropna(subset=["BRs","sin2_ele"])
                                    brs=df_br['BRs'].sum(skipna=True)/df_br["sin2_ele"].sum(skipna=True)
                                    if abs(glue-brs)>20: correction = brs
                                    else: correction = glue
                                    df_seg["STEC_sl"] = df_seg["STEC_sl"] + correction

                                    d_df_seg[iright] = df_seg
                                    iright-=1
                                elif not only_left:
                                    only_right=True


                df_filter_arc = pd.DataFrame()
                for df in d_df_seg:
                    if df is None:continue
                    df_filter_arc = pd.concat([df_filter_arc,df])

                # Calculate VTEC (still not take into account receiver bias)
                df_filter_arc["STEC_sl"] = df_filter_arc["STEC_sl"]+sat_bias
                df_filter_arc['VTEC']=df_filter_arc['STEC_sl']*df_filter_arc['cos_chi']
                df_filter_arc=df_filter_arc[["sv","lat","lon","elevation","cos_chi","STEC_sll","STEC_slp","STEC_sl","VTEC"]]
                df_filtered = pd.concat([df_filtered,df_filter_arc])

        self.df_obs = df_filtered.copy()



    def compute_receiver_bias(self):
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

        df = self.df_obs.copy()

        # Remove rows with nan STEC_sl
        df.dropna(subset=["STEC_sl","elevation"],inplace=True)
        df = df[df["elevation"]>30*np.pi/180]

        # Compute cos\chi for full serie
        df["ci"] = np.cos(np.arcsin(R_E*np.cos(df["elevation"])/(R_E+self.h)))

        # Coefficients of the cuadratic error function that will be computed
        a,b=0,0

        # Create list of time of the station removing duplicates
        time_series = df.groupby(level="time").size()
        time_series = time_series.index

        #print (len(time_series))
        # Compute sums
        for t in time_series:
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
                si = row["STEC_sl"]
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
            #print (a,b,si,ci,sum_ci,sum_ci2,sum_sici,sum_sici2,N)
            if math.isnan(a) or math.isnan(b): sys.exit()

        if a==0: return float("nan")
        # root of error function = receiver bias.
        br = b/a

        self.br = br

        d = {"station":[self.station],"br":[self.br]}

        f_br = st.root_dir + "TEC/" + str(self.year) + "/receiver_bias.csv"
        if not os.path.exists(st.root_dir + "TEC/"+ str(self.year)):
            try: os.mkdir(st.root_dir + "TEC/"+ str(self.year))
            except OSError as e:
                    if e.errno!=17: print ("FAIL creation of directory "+st.root_dir + "TEC/"+ str(self.year), e )
            else: print ("Successfully created the directory "+st.root_dir + "TEC/"+ str(self.year))
                    
        if os.path.exists(f_br):
            df_br = pd.read_csv(f_br).set_index("station")
            df_br = df_br[df_br.index!=self.station]
            df_br = pd.concat([df_br,pd.DataFrame(d).set_index("station")])
            df_br.to_csv(f_br)
        else:
            df_br = pd.DataFrame(d).set_index("station")
            df_br.to_csv(f_br)


    def get_receiver_bias(self,force_compute=False):
        ''' Function that returns the bias of the station given in argument
                If the information is not in stations.csv just return 0
                    (with warning, receiver bias should not be ignored)
            * Input: - station: string. station name
                    - compute_if_nan: bool, False by default, if True, calls compute_receiver_bias
            * Output: float. bias value in TECU
        '''
        receiver_bias = pd.read_csv(st.root_dir+"stations.csv").set_index("station")

        if not self.station in receiver_bias.index:
            print (f"WARNING must update stations.csv, with {station} not in csv, return br=0")
            return 0

        br = receiver_bias.loc[self.station]["br"]
        if force_compute:
            self.br = self.compute_receiver_bias()
            return self.br
        self.br=br
        return br

    def add_receiver_bias(self):
        if len(self.df_obs)==0: return
        self.compute_receiver_bias()

        if math.isnan(self.br):
            print (f"Warning: will null receiver bias")
            self.br=0

        # Remove rows with nan STEC_sl
        self.df_obs.dropna(subset=["STEC_sl"],inplace=True)

        # Compute VTEC with receiver bias
        self.df_obs["VTEC"]=(self.df_obs["STEC_sl"]-self.br)*np.cos(np.arcsin(R_E*np.cos(self.df_obs["elevation"])/(R_E+self.h)))

    def to_feather(self, f_feather):
        print (f_feather)
        self.df_obs[["sv","lat","lon","elevation","STEC_slp","STEC_sl","VTEC"]].reset_index().to_feather(f_feather)


    def compute_vtec(self):	
 
        print ("Converting code and pseudorange from observation rinex to slant tec")
        if not self.rinex_to_stec(): return
    	
        print ("Calculating satellite position from navigation rinex")
        self.add_satellite_pos()
    
        print ("Calculating baseline to correct Slant TEC")
        self.add_baseline()
    
        print ("Calculating receiver bias, correct Slant TEC, compute VTEC")
        self.add_receiver_bias()
        
        print ("Save to feather:",st.root_dir + "TEC/" + str(self.year) + "/" + self.station + str(self.doy) + "tec.feather")
        self.to_feather(st.root_dir + "TEC/" + str(self.year) + "/" + self.station + str(self.doy) + "tec.feather" )
        #self.to_feather(self.f_feather)
    
