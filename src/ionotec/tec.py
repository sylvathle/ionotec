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
#import georinex as gr

#import stations as st
import pandas as pd
import scipy.constants as csts
from datetime import datetime, timedelta
import numpy as np
import sys,os,subprocess
#import matplotlib.pyplot as plt
#import matplotlib.dates as md
import pymap3d as pm
#from os import listdir,path,mkdir
#from os.path import isfile, join
#import julian
#import math	
#import time
#import re


#import psutil
#import gc
#from pympler import asizeof

from pathlib import Path

#from .stations import root_dir
from . import stations as st
from . import rinex as rx
from . import freq
from . import DCB
from . import gnss
from . import reconstruction as reco


#import warnings




#pd.options.mode.chained_assignment = None

# Earth Radius
R_E = 6371000


class tec:

    def __init__(self):

        #Resolution that will be used for the process (should by 60 seconds)
        #self.resolution = 60
        self.h = 350000
        self.rDCB_interval = timedelta(days=3)

        self.list_obs_stations = []
        self.list_tec_stations = {}
        self.list_station_df_obs = {}

    def set_root_dir(self,root_dir):
        st.root_dir = root_dir
        Path(st.root_dir).mkdir(parents=True, exist_ok=True)

    def set_rDCB_interval(self,rdcb_interval): self.rDCB_interval = rdcb_interval

    def set_h(self,h): self.h = h

    def get_observation_station_list(self):
        return self.list_obs_stations

    '''
    def load_rinex_folder(self,folder,date_min=None, date_max=None,h=350000):
        self.h = h

        self.datemin = date_min
        self.datemax = date_max
        self.source_data_folder = rinex_folder
        self.prepare_files()

        self.gnss = gnss.gnss(self.list_f_nav)

        self.sat_dcb = DCB.load_dcb(self.list_f_dcb,datemin=self.datemin,datemax=self.datemax)
        #
    '''



    def run(self,list_f_obs):
        #self.h = h

        #self.list_f_obs = list_f_obs 
        
        # List of navigation files
        #self.list_f_nav = list_f_nav if not list_f_nav is None else []
        
        # File containing satellite bias
        #self.list_f_dcb = list_f_dcb if not list_f_dcb is None else []



        if not os.path.exists(st.root_dir + "TEC/"):
            os.mkdir(st.root_dir + "TEC/")
        
        if not isinstance(list_f_obs,list): list_f_obs = [list_f_obs]
        
        for f_obs in list_f_obs:

            if isinstance(f_obs,str): f_obs = Path(f_obs)
            if not f_obs.is_file(): 
                print (f_obs.resolve(),"doesn't exist, Skipping it")
                continue
            station = f_obs.name[:4].lower()

            if station not in self.list_tec_stations.keys():
                self.list_tec_stations[station] = tec_station([f_obs])
                self.list_obs_stations.append(station)
            else:
                self.list_tec_stations[station].add_f_obs(f_obs)

        Path(st.root_dir).mkdir(parents=True, exist_ok=True)
        print ("Output directory is :",st.root_dir)


        for station in self.list_tec_stations.keys():
            print ("Running tec for station: "+station)
            self.list_station_df_obs[station] = self.list_tec_stations[station].run(self.h,self.rDCB_interval)
            #print (self.list_station_df_obs[station])


        return self.list_station_df_obs




class tec_station:

    def __init__(self,list_f_obs):
        
        self.list_f_obs = list_f_obs

        self.channels = {}
        self.list_df = {}

        self.station = ""
        self.datemin = None
        self.datemax = None
        

    def add_f_obs(self,f_obs):
        self.list_f_obs.append(f_obs)


    def load_files(self):

        for f_obs in sorted(self.list_f_obs):

            rfile = rx.rinex(f_obs)
            header = rfile.read_header()

            ## Avoid running files that are not in the requested date range
            #datelim_cond = (
		    #    (self.datemin is None or header["t_first_obs"].date() >= datemin.date()) and
		    #    (self.datemax is None or header["t_first_obs"].date() <= datemax.date())
            #)
            #if not datelim_cond: continue
            
            self.station = header["name_station"].replace(" ","").lower()
            self.coord = header['position']
            list_df_f_obs = rfile.read_data()

            for const, df in list_df_f_obs.items():
                df.set_index('time',inplace=True)
                if const in self.list_df.keys():
                    self.list_df[const] = pd.concat([self.list_df[const],df])
                else:
                    self.list_df[const] = df

            

        self.datemin = None
        for constellation in self.list_df.keys():
            if self.datemin is None: self.datemin = min(self.list_df[constellation].index)
            else: self.datemin = min(self.datemin, min(self.list_df[constellation].index))


        self.datemax= None
        for constellation in self.list_df.keys():
            if self.datemax is None: self.datemax = max(self.list_df[constellation].index)
            else: self.datemax = max(self.datemax, max(self.list_df[constellation].index))


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

        self.sat_dcb = DCB.load_dcb(datemin=self.datemin,datemax=self.datemax)


    '''
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

    '''


    def rinex_to_stec(self):
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
                chan = {"C1":"C1W","C2":"C2W","L1":"L1","L2":"L2","S1":"S1","S2":"S2"}
            elif ("P1" in list_cols) and ("P2" in list_cols) and ("L1" in list_cols) and ("L2" in list_cols):
                C1,C2,L1,L2,S1,S2 = "P1","P2","L1","L2","S1","S2"
                chan = {"C1":"C1W","C2":"C2W","L1":"L1","L2":"L2","S1":"S1","S2":"S2"}

            if chan:

                self.channels[const] = []
                self.channels[const].append(chan)
                #self.list_df[const].set_index("time",inplace=True)
                #self.t_min[const] = min(self.list_df[const].index)
                #self.t_max[const] = max(self.list_df[const].index)  
                
                self.list_df[const].rename(columns={S1:'S1', S2:'S2'},inplace=True)      
                
                self.list_df[const]['STEC_l'] = (self.list_df[const][L1]*freq.gps_lambda1-self.list_df[const][L2]*freq.gps_lambda2)*freq.gps_alpha/1e16
                self.list_df[const]['STEC_p'] = (self.list_df[const][C2]-self.list_df[const][C1])*freq.gps_alpha/1e16

                #self.list_df[const] = self.list_df[const][['sv',"STEC_l","STEC_p"]]
                self.list_df[const] = self.list_df[const].dropna(subset=["STEC_l","STEC_p"])
                
                self.list_df[const]["C1"] = chan["C1"]
                self.list_df[const]["C2"] = chan["C2"]
                

            else: del self.list_df[const]

            
        ### GLONASS
        if 'R' in self.list_df.keys():

            const = 'R'
            list_sv = self.list_df['R']['sv'].unique().tolist()
            if len(list_sv)==0: return   

            list_cols = self.list_df[const].columns
            
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
                chan = {"C1":"","C2":"","L1":L1,"L2":L2,"S1":"S1","S2":"S2"}
                if ("P2" in list_cols):
                    C2="P2"
                    chan["C2"] = "C2P"
                    if ("P1" in list_cols): 
                        C1="P1"
                        chan["C1"]="C1P"
                    elif ("C1" in list_cols): 
                        C1="C1"
                        chan["C1"]="C1C"
                elif ("C2" in list_cols):
                    C2="C2"
                    chan["C2"] = "C2C"
                    if ("P1" in list_cols): 
                        C1="P1"
                        chan["C1"]="C1P"
                    elif ("C1" in list_cols): 
                        C1="C1"
                        chan["C1"]="C1C"
                    

            if chan:
                self.channels[const] = []
                self.channels[const].append(chan)
                
                #self.list_df[const].set_index("time",inplace=True)
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
                
                self.list_df[const].rename(columns={S1:'S1', S2:'S2'},inplace=True)
                
                
                self.list_df[const]["STEC_p"] = (self.list_df[const][C2] - self.list_df[const][C1])*self.list_df[const]["alpha"]/1e16
                self.list_df[const]["STEC_l"] = (self.list_df[const]["lambda1"]*self.list_df[const][L1] - \
                                               self.list_df[const]["lambda2"]*self.list_df[const][L2])*self.list_df[const]["alpha"]/1e16
                self.list_df[const].dropna(subset=["STEC_p","STEC_l"],inplace=True)
                
                if len(self.list_df[const])!=0:
                    self.list_df[const]["C1"] = chan["C1"]
                    self.list_df[const]["C2"] = chan["C2"]
                    self.list_df[const] = self.list_df[const][['sv',"C1","C2",'S1','S2',"STEC_l","STEC_p"]]
                    #self.t_min[const] = min(self.list_df[const].index)
                    #self.t_max[const] = max(self.list_df[const].index)
            else: del self.list_df[const]

 
        ## GALILEO
        if "E" in self.list_df.keys():

            const = 'E'
            self.list_df[const] = self.list_df[const].dropna(axis=1, how='all')
            list_cols = self.list_df[const].columns
            

            
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
            if not ((S1 in list_cols) or (S1 in list_cols)): chan = {}

            if chan:
                self.channels[const] = []
                self.channels[const].append(chan)
                
                #self.list_df[const].set_index("time",inplace=True)
                self.list_df[const]["STEC_p"] = (self.list_df[const][C2] - self.list_df[const][C1])*freq.gps_alpha/1e16
                self.list_df[const]["STEC_l"] = (freq.gps_lambda1*self.list_df[const][L1] - freq.gps_lambda5*self.list_df[const][L2])*freq.gps_alpha/1e16
                self.list_df[const]["C1"] = chan["C1"]
                self.list_df[const]["C2"] = chan["C2"]
                self.list_df[const].rename(columns={S1:'S1', S2:'S2'},inplace=True)
                self.list_df[const] = self.list_df[const][['sv',"C1","C2",'S1','S2',"STEC_l","STEC_p"]]
                #self.list_df['E'].dropna(inplace=True)
                self.list_df[const] = self.list_df[const].dropna(subset=["STEC_l","STEC_p"],how='any')

                
                #self.t_min[const] = min(self.list_df[const].index)
                #self.t_max[const] = max(self.list_df[const].index)
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
            if ("L2I" in beidu_columns) and ("L7I" in beidu_columns) and ("C7I" in beidu_columns) and ("C2I" in beidu_columns) and ('S2I' in beidu_columns) and ('S7I' in beidu_columns):
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

            if ("L2I" in beidu_columns) and ("L6I" in beidu_columns) and ("C6I" in beidu_columns) and ("C2I" in beidu_columns) and ('S2I' in beidu_columns) and ('S6I' in beidu_columns):
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
                self.list_df[const] = self.list_df[const][['sv',"C1","C2",'S1','S2',"STEC_l","STEC_p"]]
                self.list_df[const].dropna(subset="STEC_l",inplace=True)         
                #self.t_min[const] = min(self.list_df[const].index)
                #self.t_max[const] = max(self.list_df[const].index)
            else:
                del self.list_df[const]
                
                
        #### QZSS
        if "J" in self.list_df.keys():
            const = 'J'
            
            self.list_df[const] = self.list_df[const].dropna(axis=1, how='all')
            chan = {}

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
                    S2 = 'S5'+c
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
                

                
                qzss_alpha = freq.gps_f1**2*freq.gps_f5**2/(freq.gps_f1**2-freq.gps_f5**2)/40.318
                #self.list_df[const].set_index("time",inplace=True)
                self.list_df[const]["STEC_l"] = (self.list_df[const][L1]*freq.gps_lambda1-self.list_df[const][L2]*freq.gps_lambda5)*qzss_alpha/1e16
                self.list_df[const]["STEC_p"] = (self.list_df[const][C2]-self.list_df[const][C1])*qzss_alpha/1e16
                self.list_df[const]["C1"] = chan["C1"]
                self.list_df[const]["C2"] = chan["C2"]
                self.list_df[const].rename(columns={S1:'S1', S2:'S2'},inplace=True)

                self.list_df[const].dropna(subset=["STEC_l","STEC_p"],inplace=True)
                self.list_df[const] = self.list_df[const][['sv',"C1","C2",'S1','S2',"STEC_l","STEC_p"]]

                #self.t_min[const] = min(self.list_df[const].index)
                #self.t_max[const] = max(self.list_df[const].index)

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
                #self.t_min[const] = min(self.list_df[const].index)
                #self.t_max[const] = max(self.list_df[const].index)        

                sbas_alpha = freq.gps_f1**2*freq.gps_f5**2/(freq.gps_f1**2-freq.gps_f5**2)/40.318
                self.list_df[const]['STEC_l'] = (self.list_df[const][L1]*freq.gps_lambda1-self.list_df[const][L2]*freq.gps_lambda5)*sbas_alpha/1e16
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

            #for s in self.list_df[const]["sv"].values:
            #    if s not in self.sv:
            #        self.sv.append(s)
            #        self.borders[s]=[]
                            

        # Remove constellation that have no data
        for const in const_to_del:
            del self.list_df[const]
            del self.channels[const]

        dict_oper = {"S1": "mean", "S2": "mean", "STEC_l": "mean","STEC_p": "mean"}
        for const in self.list_df.keys():
            self.list_df[const]=self.list_df[const].groupby(["sv","C1","C2"]).resample("1min").agg(dict_oper).reset_index().set_index('time').dropna(subset=['STEC_l','STEC_p'])
            
 
        return True
        
        
        

    def add_satellite_pos(self):
        const_without_pos = []
        
        for const in self.list_df.keys():

            self.list_df[const] = self.gnss.getElevation(self.list_df[const],self.coord)
            self.list_df[const] = self.list_df[const][self.list_df[const]["elevation"]>0]
            
            if len(self.list_df[const])==0: 
                const_without_pos.append(const)
                continue

            self.list_df[const] = self.gnss.getPiercingPoint(self.list_df[const],self.coord,self.h)
            # cos(chi) to calculate VTEC from STEC
            self.list_df[const]['cos_chi'] = np.cos(np.arcsin(R_E*np.cos(self.list_df[const]["elevation"])/(R_E+self.h)))   
            self.list_df[const] = self.list_df[const][["sv","C1","C2","elevation","cos_chi","lat","lon","alt","S1","S2","STEC_l","STEC_p"]]
                        
        for const in const_without_pos:
            del self.list_df[const]
            if const in self.channels.keys(): del self.channels[const]
        return True
        
        
        
    def add_baseline(self):
    
        N_min_stec_p = 10
        signal_strength_threshold = 38

        const_to_del = []
    
        for const, df_data in self.list_df.items():

            if len(df_data)==0: continue

            df_data = df_data.dropna(subset=["elevation"])

            t_begin = df_data.index[0]
            t_end = df_data.index[-1]

            df_filtered = pd.DataFrame()

            list_sv = df_data["sv"].unique().tolist()
            
            df_sats_corrected = pd.DataFrame()

            for sat in list_sv:
            
                #if sat!='G04': continue
                #print (sat)


                df_sat = df_data[df_data["sv"]==sat]
                #df_sat_filter = pd.DataFrame()
                
                df_sat_corrected = pd.DataFrame()

                for channel in self.channels[const]:

                    #if channel['C2']!='C6I': continue
                    chan_filter = (df_sat["C1"]==channel["C1"]) & (df_sat["C2"]==channel["C2"])
                    df_sat_chanel = df_sat[chan_filter]
                        
                    if len(df_sat_chanel["STEC_l"].dropna())==0: continue

                    # Make list of arcs of satellite using elevation information
                    elevations = df_sat_chanel["elevation"]
                    list_arcs = gnss.get_arcs(elevations,t_begin,t_end)

                    # Squared of sin of elevation for baseline (Brs) pondering
                    df_sat_chanel["sin2_ele"] = np.sin(df_sat_chanel['elevation'])**2                    

                    n_arc=0
                    
                    for iarc, arc in enumerate(list_arcs):

                        # Extract data from satellite arc
                        date_filter = (df_sat_chanel.index>=arc["start"]) & (df_sat_chanel.index<=arc["end"])
                        df_arc = df_sat_chanel[date_filter]

                        if len(df_arc)<=8:continue      
                                     
                        df_arc.dropna(subset=['STEC_l'],inplace=True)

                        df_arc = reco.correct_signal(df_arc)

                        df_sats_corrected = pd.concat([df_sats_corrected,df_arc])    
        
                    if len(df_sats_corrected)==0: continue
            
            self.list_df[const] = df_sats_corrected

            if len(self.list_df[const])==0: 
                const_to_del.append(const)
                continue

            self.list_df[const]['time_day'] = self.list_df[const].index.normalize()
            self.list_df[const].reset_index(inplace=True)   # 'time' becomes a column

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

            self.list_df[const]["STEC_l"] += self.list_df[const]["dcb"]
            self.list_df[const]["VTEC"]=self.list_df[const]["STEC_l"]*self.list_df[const]['cos_chi']

        for const in const_to_del:
            del self.list_df[const]



    def correct_receiver_DCB(self):
    
        ''' Sylvain Blunier 08/2026 v2.0.0
                Function that computes and returns the biais of some station and correct STEC and VTEC
                The algorithm minimize the sum of variances between satellites at each time
                *section 4.1 of https://hal.archives-ouvertes.fr/hal-00317176/file/angeo-21-2083-2003.pdf
                *algebra: https://colab.research.google.com/drive/1UCZHR0t-9jyyjAnLuMN3N0Z2NB6tgI_l?usp=sharing
        '''
        
        # Bias of the antenna, calculated by method "compute_reveiver_bias"
        # Value stored in csv "stations.csv"
        dict_rDCB = {'station':[],'constellation':[],'C1':[],'C2':[],'time_i':[],'time_f':[],'DCB':[]}
           
        for const in self.list_df.keys():

            if len(self.list_df[const])==0: continue
            
            d = self.datemin
            
            self.list_df[const]["rDCB"] = float("NaN")
            
            while d<self.datemax:

                interval_process_mask = (self.list_df[const].index>=d) & (self.list_df[const].index<d+self.rDCB_interval)
                
                for chanel in self.channels[const]:
                    chanel_process_mask = interval_process_mask & (self.list_df[const]["C1"]==chanel["C1"]) & (self.list_df[const]["C2"]==chanel["C2"])
                    df_interval_chanel = self.list_df[const][chanel_process_mask]      
                    
                    if len(df_interval_chanel)==0: continue              

                    rDCB = reco.process_receiver_dcb(df_interval_chanel,elevation_filter = 10)

                    self.list_df[const].loc[chanel_process_mask, "rDCB"] = rDCB

                    dict_rDCB['station'].append(self.station)
                    dict_rDCB['time_i'].append(min(df_interval_chanel.index))
                    dict_rDCB['time_f'].append(max(df_interval_chanel.index))
                    dict_rDCB['constellation'].append(const)
                    dict_rDCB['C1'].append(chanel["C1"])
                    dict_rDCB['C2'].append(chanel["C2"])
                    dict_rDCB['DCB'].append(rDCB)
                    
                d += self.rDCB_interval

        
        ## Stored DCB information    
        self.df_br = pd.DataFrame(dict_rDCB)#.set_index("time_i")
        self.df_br.set_index("station",inplace=True)
        
        f_br = st.root_dir + "TEC/DCB_receiver.csv"
        if os.path.exists(f_br):
            df_br_stored = pd.read_csv(f_br).set_index("station")
            df_br_stored = df_br_stored[df_br_stored.index!=self.station]
            self.df_br = pd.concat([df_br_stored,self.df_br])
        self.df_br.to_csv(f_br)

        const_to_del = []
        
        ## Correct STEC and VTEC with rDCB
        for const in self.list_df.keys():

            if len(self.list_df[const])==0: continue 

            self.list_df[const]["STEC"] = self.list_df[const]["STEC_l"]-self.list_df[const]["rDCB"]
            self.list_df[const]["VTEC"] = self.list_df[const]["STEC"] * self.list_df[const]["cos_chi"] 
                #np.cos(np.arcsin(R_E*np.cos(self.list_df[const]["elevation"])/(R_E+self.h)))

            self.list_df[const].reset_index(inplace=True)
            self.list_df[const].dropna('STEC',inplace=True)
            if len(self_df[const])==0:
                const_to_del.append(const)
                continue
            self.list_df[const] = self.list_df[const][["time","sv","lat","lon","elevation","cos_chi","STEC_l","STEC","VTEC","rDCB","dcb"]]
            #else:
            #    self.list_df[const].reset_index(inplace=True)
            #    self.list_df[const] = self.list_df[const][["time","sv","lat","lon","elevation","STEC"]]

            self.list_df[const] = self.list_df[const].groupby(by=["time","sv"],as_index=False).mean()
            self.list_df[const].set_index('time',inplace=True)

       for const in const_to_del:
           del self.list_df[const]

        
    def to_feather(self):
        df_obs = pd.DataFrame()

        for const in self.list_df.keys():
            df_obs = pd.concat([
                df_obs,
                self.list_df[const]#[["time","sv","lat","lon","elevation","STEC_l","VTEC"]]
            ])
            

        feather_path = st.root_dir + "TEC/" + self.station
        df_obs.to_feather(feather_path+".feather")
    
    #def compute_vtec(self):	
    def run(self, h, rdcb_interval ,store = True):
        #self.list_df = {}
        #self.channels = {}
        self.h = h
        self.rDCB_interval = rdcb_interval

        print ("Load files")
        self.load_files()

        print ("RINEX to STEC")
        self.rinex_to_stec()
        
        print ("Add GNSS Navigation Information (IPP & Elevation)")
        self.add_satellite_pos()
        
        print ("Correct Slant TEC discontinuities and baseline")
        self.add_baseline()

        print ("Calculating receiver DCB, correct Slant TEC, compute VTEC")
        self.correct_receiver_DCB()

        df_obs = pd.DataFrame()

        for const in self.list_df.keys():
            df_obs = pd.concat([
               df_obs,
               self.list_df[const][["sv","lat","lon","elevation","cos_chi","STEC","VTEC"]]
            ])

        
            

        if len(df_obs)>0:
            if store:
                for year in df_obs.index.year.unique():
                    # Filter the dataframe for the current year
                    df_year = df_obs[df_obs.index.year == year]
                    folder = st.root_dir + "TEC/" + str(year) + "/"
                    Path(folder).mkdir(parents=True, exist_ok=True)
                    feather_path = folder + self.station
                    df_year.to_feather(feather_path+".feather")
            return df_obs.dropna(subset=["STEC","VTEC"])
        else:
            return None
        

