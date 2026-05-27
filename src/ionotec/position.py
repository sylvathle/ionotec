import sys
import datetime
import time

import shutil
#shutil.rmtree('output')

import pymap3d as pm

from . import stations
from . import gnss
#from . import tec
from . import rinex
from . import freq


from os import listdir, makedirs
from os.path import isfile, join, exists
from pathlib import Path

#import georinex as gr
import pandas as pd

import pymap3d as pm
import numpy as np

import scipy.constants as csts


def in_large_gap(regular_idx, raw_idx, threshold):
    """Returns a boolean mask: True where a grid point is inside a gap >= threshold."""
    raw_sorted = raw_idx.sort_values()
    mask = np.zeros(len(regular_idx), dtype=bool)

    for i, t in enumerate(regular_idx):
        # last raw point at or before t
        before = raw_sorted[raw_sorted <= t]
        # first raw point at or after t
        after = raw_sorted[raw_sorted >= t]

        if before.empty or after.empty:
            mask[i] = True   # outside data range → NaN
            continue

        gap = after[0] - before[-1]
        if gap >= threshold:
            mask[i] = True

    return mask




class position:
    
    def __init__(self,f_obs,f_nav,f_bias):

        self.f_obs = f_obs
        self.f_nav = f_nav
        self.f_bias = f_bias

        self.df_pos = pd.DataFrame()

        self.C1 = 'C1C'
        self.C2 = 'C2W'
        self.C1 = 'P1'
        self.C2 = 'P2'


    def computePosition(self):

        self.computeObs()
        self.computeNav()
        print (self.f_nav)
        
    def computeNav(self):

        # First, load the needed data from reported navigation rinex files
        list_nav = {"time":[],"sv":[],"SVclockBias":[],"X":[],"Y":[],"Z":[]}
        
        for rinex_nav in self.f_nav:
            nav = rinex.rinex(rinex_nav)
            head_nav = nav.read_header()
            df_nav = nav.read_data()
            #print (df_nav)
            #list_t = df_nav['time'].unique()

            list_sv = self.df_pos['sv'].unique()
            #print ('Here')
            for sv in list_sv:
                #print (sv)
                #if sv!="G14": continue
                df_pos_sv = self.df_pos[self.df_pos['sv']==sv]
                list_t = df_pos_sv['time'].unique()
                
                df_nav_sv = df_nav[df_nav['sv']==sv]
                #list_t = df_nav_sv['time'].unique()
                #print (df_nav_sv)
                irow=0
                row = df_nav_sv.iloc[0]
                #print (row)
                #print (list_t)
                for t in list_t:
                    #if t > datetime.datetime(2024,5,7,20,1,0,0): continue
                    #if t < datetime.datetime(2024,5,7,19,59,0,0): continue
                    #df_t = df_nav_sv[df_nav_sv['time']==t]
                    #print ('\t',t)
                    while (row['time']<t) or (irow>=len(df_nav_sv)):
                        irow+=1
                        row = df_nav_sv.iloc[irow]
                    #print ('\t',row)
                    X,Y,Z = gnss.gps_nav_to_XYZ(row,t)
                    #print (sv,t,row['time'],X,Y,Z)
                    list_nav['time'].append(t)
                    list_nav['sv'].append(sv)
                    list_nav['SVclockBias'].append(row['SVclockBias'])
                    list_nav['X'].append(X)
                    list_nav['Y'].append(Y)
                    list_nav['Z'].append(Z)


                    #print (t,row['time'],X,Y,Z)
                #sys.exit()
                    #print (t,len(df_t))
                #for i,row in df_t.iterrows(): 
                    #print (row)
                 #   X,Y,Z = gnss.gps_nav_to_XYZ(row,t)
                 #   list_nav['time'].append(row['time'])
                 #   list_nav['sv'].append(row['sv'])
                 #   list_nav['SVclockBias'].append(row['SVclockBias'])
                 #   list_nav['X'].append(X)
                 #   list_nav['Y'].append(Y)
                 #   list_nav['Z'].append(Z)
                    
        df_nav = pd.DataFrame(list_nav)
        #df_nav = df_nav[~df_nav.index.duplicated(keep="last")]

        #print (df_nav[df_nav['sv']=='G03'])
        
        self.df_pos = self.df_pos.merge(df_nav, on=['time', 'sv'], how='left')
        self.df_pos.dropna()
        #print (self.df_pos)
        #merged = merged.set_index('time')

        #cols_to_interpolate = ['SVclockBias', 'X', 'Y', 'Z']

        #self.df_pos = (
        #    merged
        #        .sort_index()  # important after outer join — rows may be out of order
        #        .groupby('sv', group_keys=True)
        #        .apply(lambda g: g[cols_to_interpolate]
        #        .interpolate(method='polynomial',order=3)
        #        .join(g.drop(columns=cols_to_interpolate)))
        #).reset_index().dropna()

        df_br = pd.read_csv("output/TEC/None/receiver_bias.csv")
        br = float(df_br[df_br['constellation']=='G']['br'][0])*1e16

        list_sv = self.df_pos['sv'].unique()
        self.df_pos['alpha'] = float('NaN')
        self.df_pos['f1'] = float('NaN')
        self.df_pos['f2'] = float('NaN')
        self.df_pos['dcb(ns)'] = float('NaN')
        for sv in list_sv :
            #print (sv)
            row_sv = (self.df_pos['sv']==sv)
            alpha = freq.getAlpha(sv,self.C1, self.C2)
            #print (alpha)
            self.df_pos.loc[row_sv,'alpha'] = alpha
            self.df_pos.loc[row_sv,'f1'] = freq.getFrequency(sv,self.C1)
            self.df_pos.loc[row_sv,'f2'] = freq.getFrequency(sv,self.C2)
            #print (self.df_pos)
            dcb = gnss.getBias_fromfile(sv,self.f_bias[0],self.C1, self.C2)
            self.df_pos.loc[row_sv,'dcb(ns)'] = dcb
            
        self.df_pos['br_station(ns)'] = br/self.df_pos['alpha']

        #Correct with satellite bias
        #self.df_pos[self.C1] = self.df_pos[self.C1] + csts.c*self.df_pos["SVclockBias"]
        #self.df_pos[self.C2] = self.df_pos[self.C2] + csts.c*self.df_pos["SVclockBias"]

        # eliminate ionospheric effects
        self.df_pos['a'] = self.df_pos['f1']**2 / (self.df_pos['f1']**2 - self.df_pos['f2']**2)
        self.df_pos['b'] = self.df_pos['f2']**2 / (self.df_pos['f1']**2 - self.df_pos['f2']**2)
        self.df_pos['d'] = self.df_pos[self.C1]
        #self.df_pos['d'] = self.df_pos[self.C1]*self.df_pos['a'] - self.df_pos[self.C2]*self.df_pos['b']
        #self.df_pos['d'] += csts.c*self.df_pos['dcb(ns)']*1e-9
        #self.df_pos['d'] -= csts.c*self.df_pos['br_station(ns)']*1e-9

        self.df_pos['d'] += csts.c*self.df_pos["SVclockBias"]

        #print ("clock Bias (meters):",csts.c*self.df_pos["SVclockBias"])


        
        self.df_pos.to_csv("df_pos.csv",index=False)

            
    
    def computeObs(self):
        
        for rinex_obs in self.f_obs:
            obs = rinex.rinex(rinex_obs)
            head_obs = obs.read_header()
            list_df_obs =  obs.read_data()
            df_obs_gps = list_df_obs['G']
            #df_obs_gps = df_obs_gps[['time', 'sv', 'C1C', 'C2W', 'L1C', 'L2W' ]]
            #df_obs_gps.dropna(subset=['C1C', 'C2W', 'L1C', 'L2W'],inplace=True)
            df_obs_gps = df_obs_gps[['time', 'sv', self.C1, self.C2 ]]
            df_obs_gps.dropna(subset=[self.C1, self.C2],inplace=True)
            self.df_pos = pd.concat([self.df_pos,df_obs_gps])
        self.df_pos= self.df_pos.sort_values('time')
        #print (self.df_pos)
    