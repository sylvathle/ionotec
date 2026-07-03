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
#############################################################################

import os,sys,shutil

#import sys
import numpy as np
import georinex as gr
import datetime
import math
import numpy as np
#import urllib.request
import pymap3d as pm
#from unlzw import unlzw

from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
from pathlib import Path

import psutil

pd.options.mode.chained_assignment = None

#import random
#import time

#import matplotlib.pyplot as plt
from . import stations as st
from . import rinex as rx

from os import listdir
from os.path import isfile, join

#Ionosphere altitude
Alt_i = 300000
# Earth Radius
R_E = 6370000
omega_e = 7.2921151467e-5
GM = 3.986005e14
sqrt_GM = math.sqrt(GM)
C20 = -1082.63e-6
a_major_ellipsoid_axis = 6378136

#WG84_to_PZ90=np.matrix([[1,-0.33/3600.0,0],[0.33/3600.0,1,0],[0,0,1]])


GPS_EPOCH  = datetime.datetime(1980, 1, 6)
BDS_EPOCH  = datetime.datetime(2006, 1, 1)

def getDateFromSat(Week, Toe, system='G'):
    epoch = BDS_EPOCH if system == 'C' else GPS_EPOCH
    return epoch + datetime.timedelta(weeks=Week, seconds=Toe)

def removeOutsiders(df,thres=datetime.timedelta(days=1)):
    #### Remove outsiders
    T = df.index
    diff_t = []
    for it in range(len(T)-1):
        diff_t.append(T[it+1]-T[it])
    
    outsider=[]
    istail = False
    thres_t = datetime.timedelta(days=1)
    #hasValidData = False
    for d in range(len(diff_t)):
        if diff_t[d]>thres_t: 
            outsider.append(True)
        else:
            if not istail: 
                istail=True
                outsider.append(False)
            outsider.append(False)
            #hasValidData = True


    #if not hasValidData:
    #   return pd.DataFrame()

            
    df["out"] = outsider
    df = df[df["out"]==False]
    df.drop(columns=["out"],inplace=True)    
    return df

def ephemeris_to_XYZ(data,date):

    Toe = data['Toe']
    #TGD = data['TGD']
    IDOT = data['IDOT']
    #IODC = data['IODC']
    Week = data['Week']
	#Square root of the semi-major axis m^{1/2}
    sqrtA =  data['sqrtA']
	#Eccentricity
    e =  data['Eccentricity']
	#Inclination angle at reference time
    Io =  data['Io']
	#Longitude of ascending node at reference time
    Omega0 =  data['Omega0']
	#Argument of perigee (semicircles)
    omega = data['omega']
	#Mean anomaly at reference time (semicircles)
    M0 = data['M0']
	#Mean motion difference from computed value
    DeltaN = data['DeltaN']
	#Rate of change of right ascension
    OmegaDot = data['OmegaDot']
	#Amplitude of the sine harmonic correction term to the argument of latitude (rad)
    Cus = data['Cus']
	#Amplitude of the cosine harmonic correction term to the argument of latitude (rad)
    Cuc = data['Cuc']
	#Amplitude of the sine harmonic correction term to the angle of inclination (rad)
    Cis = data['Cis']
	#Amplitude of the cosine harmonic correction term to the angle of inclination (rad)
    Cic = data['Cic']
	#Amplitude of the sine harmonic correction term to the orbit radius (m)
    Crs = data['Crs']
	#Amplitude of the cosine harmonic correction term to the orbit radius (m)
    Crc = data['Crc']
    d0 = getDateFromSat(Week,Toe,data['sv'][0])

    
    tk = (date-d0).total_seconds()
    #if tk > 302400: 
    #    print ("big tk",tk)
    #    return [float('NaN'),float('NaN'),float('NaN')]
    #    #tk -= 604800
    #elif tk < -302400: 
    #    print ("big tk",tk)
    #    return [float('NaN'),float('NaN'),float('NaN')]
    #    #tk += 604800


    n0 = sqrt_GM/sqrtA**3
    n = n0 + DeltaN
    
    Mk = M0 + n*tk
    Ek = Mk
    for i in range(3):
        Ek = Ek - (Ek-e*math.sin(Ek)-Mk)/(1-e*math.cos(Ek))

    vk = 2*math.atan(math.sqrt((1+e)/(1-e))*math.tan(Ek/2))
    
    Phik = vk + omega
    delta_uk = Cuc * math.cos(2*Phik) + Cus * math.sin(2*Phik)
    delta_rk = Crc * math.cos(2*Phik) + Crs * math.sin(2*Phik)
    delta_ik = Cic * math.cos(2*Phik) + Cis * math.sin(2*Phik)

    uk = Phik + delta_uk
    rk = sqrtA**2 * (1 - e * math.cos(Ek)) + delta_rk


    Xk_p = rk * math.cos(uk)
    Yk_p = rk * math.sin(uk)
    
    ik = Io + IDOT*tk + delta_ik
    Omegak = Omega0 + (OmegaDot - omega_e)*tk - omega_e * Toe
    Xk_orb = Xk_p * math.cos(Omegak) - Yk_p * math.cos(ik) * math.sin(Omegak)
    Yk_orb = Xk_p * math.sin(Omegak) + Yk_p * math.cos(ik) * math.cos(Omegak)
    Zk_orb = Yk_p * math.sin(ik)   



    sv = data['sv']  # e.g. 'C01', 'C03', ...
    prn = int(sv[1:])
    is_geo = sv[0] == 'C' and (1 <= prn <= 5 or 59 <= prn <= 63)
    
    if is_geo:
        # BeiDou GEO: extra rotation by 5 degrees around X, then by GAST
        phi = -5.0 * math.pi / 180.0          # 5 deg tilt correction
        GAST     = omega_e * tk   
        Omegak = Omega0 + OmegaDot * tk - omega_e * Toe
        
        # Recompute all three orbital coordinates with corrected ik
        Xk_orb = Xk_p * math.cos(Omegak) - Yk_p * math.cos(ik) * math.sin(Omegak)
        Yk_orb = Xk_p * math.sin(Omegak) + Yk_p * math.cos(ik) * math.cos(Omegak)
        Zk_orb = Yk_p * math.sin(ik)
        
        
        # Rotation matrix Rz(-GAST) * Rx(-5deg)
        Rx = np.array([
            [1,           0,            0],
            [0,  math.cos(phi),  math.sin(phi)],
            [0, -math.sin(phi),  math.cos(phi)]
        ])
        Rz = np.array([
            [ math.cos(GAST), math.sin(GAST), 0],
            [-math.sin(GAST), math.cos(GAST), 0],
            [0,               0,              1]
        ])
        xyz = Rz @ Rx @ np.array([Xk_orb, Yk_orb, Zk_orb])
        Xk, Yk, Zk = xyz

    else:
        # MEO / IGSO: same as GPS
        Xk, Yk, Zk = Xk_orb, Yk_orb, Zk_orb
    
    return [Xk,Yk,Zk]

def split_interpolate_concat(df, gap_threshold='3h'):
    """
    Split a dataframe on gaps larger than gap_threshold,
    interpolate each chunk, then concatenate.
    
    Parameters
    ----------
    df         : DataFrame with a DatetimeIndex
    columns    : list of columns to interpolate
    gap_threshold : str or pd.Timedelta, max allowed gap
    """
    threshold = pd.Timedelta(gap_threshold)

    # 1. Find where real data exists (at least one non-NaN value)
    has_data = df.notna().any(axis=1)
    real_times = df.index[has_data]

    # 2. Compute gaps between consecutive real data points
    gaps = real_times[1:] - real_times[:-1]

    # 3. Find split points: indices (in real_times) where gap exceeds threshold
    split_positions = np.where(gaps > threshold)[0] + 1  # +1 → start of new chunk

    # 4. Build index boundaries for slicing the original df
    #    Each chunk spans from one real timestamp to the next split point
    boundaries = [real_times[0]]
    for pos in split_positions:
        boundaries.append(real_times[pos - 1])  # end of previous chunk
        boundaries.append(real_times[pos])       # start of next chunk
    boundaries.append(real_times[-1])

    # Pair up into (start, end) slices
    slices = [(boundaries[i], boundaries[i + 1]) for i in range(0, len(boundaries), 2)]

    # 5. Slice, interpolate, collect
    chunks = []
    for start, end in slices:
        chunk = df.loc[start:end].copy()
        if len(chunk.dropna())<4: continue
        chunk = chunk.interpolate(method='polynomial', order=3)
        chunks.append(chunk)

    # 6. Concatenate all chunks
    return pd.concat(chunks)


def remove_duplicated_dates(file_path):
    f = open(file_path,'r')
  
    is_header = True
    i = 0
    data_line_batch = ""
    dict_prn_dates = {}
  
    new_file_lines = ""
  
    for line in f:
  
      # Check where is the header, just to copy in the new file
      if "END OF HEADER" in line:
        is_header = False
        new_file_lines += line
        continue
      # If we are in the header copy and next line
      if is_header:
        new_file_lines += line
        continue
  
      # if we are in the first line of the data for a PRN and date, read date
      if i%8==0:
        # First time there is nothing to write yet.
        if i!=0:
          # If the satellite was already found before, 
          # check that the date was not already seen
          if PRN in dict_prn_dates.keys():
            # If date not seen, add it to seen dates, and add data_batch to new file 
            #    else do nothing
            if date not in dict_prn_dates[PRN]:
              dict_prn_dates[PRN].append(date)
              new_file_lines += data_line_batch
          # If satellite not found, then create key, put date in dict
          #    and copy data_batch to new file
          else:
            dict_prn_dates[PRN] = [date]
            new_file_lines += data_line_batch
  
        PRN = line[:3]
        # Get date
        date = datetime.datetime.strptime(line[4:23],"%Y %m %d %H %M %S")
  
        # Empty data_line_batch
        data_line_batch = ""
      
      data_line_batch += line
      i += 1
  
    ## Do the last data since it is not processed at the end of the loop.
    if PRN in dict_prn_dates.keys():
      if date not in dict_prn_dates[PRN]:
        new_file_lines += data_line_batch
    else: 
      new_file_lines += data_line_batch
    f.close()
  
    f = open(file_path,'w')
    f.write(new_file_lines)
    f.close()


class gnss:
    '''
        gnss Class to perform eficiently the position and elevation of the
        satellites
        for now it only works on GPS constellation, future version should deal
        with GLONASS
    '''

    file_format = 'feather' ## csv format not implemented

    #n_sat_gps = 32
    #n_sat_glonass = 24
    #n_sat_galileo = 25
    resolution = 60

    list_f_rinex_nav = []

    rinex_doy = 1
    rinex_year = 2020
    f_doy_reported = ""
    
    
    def __init__(self,f_nav=[],datemin=None, datemax=None, form='feather'):
        self.gnss_dir = st.root_dir + "GNSS/"
        
        if not os.path.exists(self.gnss_dir):
            try: os.mkdir(self.gnss_dir)
            except OSError as e:
                if e.errno!=17: print ("FAIL creation of directory "+self.gnss_dir, e )
            else: print ("Successfully created the directory "+self.gnss_dir)

        self.datemin = datemin
        self.datemax = datemax

        self.file_format = form
        self.list_f_rinex_nav = f_nav
        self.dict_df_pos = {}
        self.df_pos = pd.DataFrame()

        self.list_constellation = ['G','E','C','J','R','S']
        self.df_nav_gnss = {}
        
        #self.df_doy_processed = pd.DataFrame()
        self.dict_doy_processed = {}# {"date":[],"G_res":[],"R_res":[]}
        self.csv_record_processing = self.gnss_dir+"/record_processing.csv"
        
        self.resolution = 1
      
        self.compute_position()


    def inform_date_processed(self,list_date,constellation):
        # Make sure there is data before 10:00
        cutoff_I = datetime.time(4, 0)
        cutoff_F = datetime.time(20, 0)
        dates = {
            dt.date()
            for dt in list_date
                if (dt.time() < cutoff_F) and (cutoff_I > dt.time())
        }


        for date in dates:  
            if date in self.dict_doy_processed.keys():
                if (self.resolution<self.dict_doy_processed[date][constellation+"_res"]) or\
                    (self.dict_doy_processed[date][constellation+"_res"]==0):
                    self.dict_doy_processed[date][constellation+"_res"] = self.resolution
            else:
                self.dict_doy_processed[date] = {"G_res":0,"R_res":0,"E_res":0,"C_res":0,"J_res":0,"S_res":0}
                self.dict_doy_processed[date][constellation+"_res"]=self.resolution
    
    # Function that compute the position of the satellites from the navigation file 
    #   with the resolution informed when instanciating the gnss object
    def compute_position(self):

        directory_path = Path(self.gnss_dir + "/")

        for const in self.list_constellation:
            self.df_nav_gnss[const] = pd.DataFrame()

        dict_file_processed = {'files':[]}
        
        if os.path.exists(self.csv_record_processing):
            df_doy_processed = pd.read_csv(self.csv_record_processing)
            dict_file_processed['files'] = df_doy_processed['files'].to_list()


        for f_rinex_nav in self.list_f_rinex_nav:
            
            if f_rinex_nav.name in dict_file_processed['files']: continue

            try: nav = rx.rinex(str(f_rinex_nav))
            except: continue
            head = nav.read_header()
            svtype = head['constellation']
           
            try: df_nav = nav.read_data()
            except:
                print ("WARNING: Could not read",str(f_rinex_nav))
                continue
            
            if len(df_nav)==0: continue


            df_nav.index.set_names(["time","sv"],inplace=True)
            df_nav.reset_index(level=["sv"],inplace=True)

            for c in df_nav.columns:
                if 'spare' in c:
                    df_nav.drop(columns=[c],inplace=True)
            df_nav.dropna(inplace=True)
            if (svtype=='G') or (svtype=='J'): df_nav.rename(columns={"GPSWeek":"Week"},inplace=True)
            if svtype=='C': df_nav.rename(columns={"BDTWeek":"Week"},inplace=True)
            if svtype=='E': df_nav.rename(columns={"GALWeek":"Week"},inplace=True)

            df_nav.index = pd.to_datetime(df_nav.index)

            dict_file_processed['files'].append(f_rinex_nav.name)

            if len(df_nav)==0: 
                continue

            date = df_nav.index[int(len(df_nav)/2)].date()
            
            df_nav.sort_index(ascending=True,inplace=True)
            self.df_nav_gnss[svtype] = pd.concat([self.df_nav_gnss[svtype],df_nav])

        min_date = self.datemin.replace(second=0, microsecond=0)
        max_date = self.datemax.replace(second=0, microsecond=0)

        time_list = []
        t = min_date

        self.resolution = 60
                            
        while t<max_date:
            time_list.append(t)
            t = t+datetime.timedelta(seconds=self.resolution)

        for const in self.df_nav_gnss.keys():
            if (len(self.df_nav_gnss[const])<=2): continue
            list_sv = self.df_nav_gnss[const]["sv"].unique().tolist()
            self.df_nav_gnss[const] = self.df_nav_gnss[const].groupby(["time","sv"]).mean()
            self.df_nav_gnss[const].reset_index(level=["sv"],inplace=True)


            for sv in list_sv:
                print (sv)

                if len(sv)>3: continue
                df_sat = self.df_nav_gnss[const][self.df_nav_gnss[const]["sv"]==sv]
                if len(df_sat)<3: continue

                #if sv=="G14": print (df_sat)

                #if os.path.exists("output/GNSS/"+sv+".feather"): continue


                #df_sat = removeOutsiders(df_sat)
                #df_sat = df_sat[ (df_sat.index>self.datemin) & (df_sat.index<self.datemax) ]
                
                #min_date = min(df_sat.index)
                #max_date = max(df_sat.index)



                if sv[0] in ['G','E','C','J']:
                    i_time_list = 0
                    date = time_list[i_time_list]
                    # Intermediate dictionary intended to contain data of the satellite under process 
                    dict_sat_pos = {"time":[],"X":[],"Y":[],"Z":[]}
                    df_sat = df_sat[df_sat['sqrtA']!=0]
                    for t_nav, row in df_sat.iterrows():
                        #if sv=="G14": print (t_nav)
                        # Calculate position for each time in time_list that are before the next available position information
                        while date < t_nav and i_time_list<len(time_list):
                            #if sv=="G14": print ("\t",date)
                            sat_pos = ephemeris_to_XYZ(row,date)
                            dict_sat_pos["time"].append(date)
                            dict_sat_pos["X"].append(sat_pos[0])
                            dict_sat_pos["Y"].append(sat_pos[1])
                            dict_sat_pos["Z"].append(sat_pos[2])
                            i_time_list += 1
                            if i_time_list==len(time_list): break
                            date = time_list[i_time_list]
    
                    # Case navigation data does not provide position until 00:00, use last available date 
                    while i_time_list<len(time_list):
                        sat_pos = ephemeris_to_XYZ(row,date)
                        dict_sat_pos["time"].append(date)
                        dict_sat_pos["X"].append(sat_pos[0])
                        dict_sat_pos["Y"].append(sat_pos[1])
                        dict_sat_pos["Z"].append(sat_pos[2])
                        i_time_list += 1
                        if i_time_list==len(time_list): break
                        date = time_list[i_time_list]
        
                    df_sat = pd.DataFrame(dict_sat_pos)
                    df_sat["time"] = pd.to_datetime(df_sat["time"])
                    df_sat.set_index("time",inplace=True)
                else:
                    new_index = df_sat.index.union(time_list)
                    df_sat = df_sat.reindex(new_index).sort_index()
                    df_sat = split_interpolate_concat(df_sat[['X','Y','Z']], gap_threshold='3h')
                    df_sat['sv'] = sv

                csv_nav_sat = self.gnss_dir + "/" + sv +".feather"
    
                if os.path.exists(csv_nav_sat):
                    df_former_sat = pd.read_feather(csv_nav_sat)
                    df_former_sat = df_former_sat[~df_former_sat.index.isin(df_sat.index)]
                    df_sat = pd.concat([df_former_sat,df_sat])
                    df_sat = df_sat.sort_index()


                df_sat.drop_duplicates().to_feather(csv_nav_sat)
        
        df_doy_processed = pd.DataFrame(dict_file_processed)
        df_doy_processed.to_csv(self.csv_record_processing,index=True)

    def load_all_sats(self):

        '''
            Constellation: G: GPS
                            R: GLONASS
                            E: GALILEO
        '''
    
        ##year = d_in.year
        first = True

        gnss_data_dir = Path(self.gnss_dir )

        folder = Path(gnss_data_dir)
        list_gnss = [f for f in folder.glob('*feather') if f.is_file()]

        tot_mem_gnss_loaded = 0
        tot_mem_gnss = 0

        for fgnss in list_gnss:

            sat_path = fgnss.resolve()
            sat_file = fgnss.name
            split_sat = sat_file.split(".")

            if split_sat[-1]!="feather": continue
            sat = split_sat[0]
            
            self.dict_df_pos[sat] = pd.read_feather(sat_path)
    
            self.dict_df_pos[sat].index = pd.to_datetime(self.dict_df_pos[sat].index)
            #print(f"dict_df_pos[{sat}]: {sys.getsizeof(self.dict_df_pos[sat])/ 1024**2:.2f} MB, len(dict_df_pos[{sat}])={len(self.dict_df_pos[sat])}")
            tot_mem_gnss_loaded += sys.getsizeof(self.dict_df_pos[sat])

            mask = (self.dict_df_pos[sat].index>=self.datemin) & (self.dict_df_pos[sat].index<=self.datemax)
            self.dict_df_pos[sat]=self.dict_df_pos[sat].loc[mask]

            tot_mem_gnss += sys.getsizeof(self.dict_df_pos[sat])


            if first:
                self.df_pos = self.dict_df_pos[sat]
                self.df_pos["sv"] = sat
                first = False
            else:
                df_inter = self.dict_df_pos[sat]
                df_inter["sv"] = sat
                self.df_pos = pd.concat([self.df_pos,df_inter])

        self.df_pos.drop_duplicates(inplace=True)


    
    # Function that loads the satellite of a specific year, files must exist, otherwise load empty
    def load_sats(self,const,d_in,d_out):

        '''
            Constellation: G: GPS
                            R: GLONASS
                            E: GALILEO
        '''
    
        ##year = d_in.year
        first = True

        gnss_data_dir = Path(self.gnss_dir )

        folder = Path(gnss_data_dir)
        list_gnss = [f for f in folder.glob(const+'*feather') if f.is_file()]
        

        for fgnss in list_gnss:

            sat_path = fgnss.resolve()
            sat_file = fgnss.name
            split_sat = sat_file.split(".")

            if split_sat[-1]!="feather": continue
            sat = split_sat[0]
            
            self.dict_df_pos[sat] = pd.read_feather(sat_path)


            self.dict_df_pos[sat].index = pd.to_datetime(self.dict_df_pos[sat].index)
            mask = (self.dict_df_pos[sat].index>=d_in) & (self.dict_df_pos[sat].index<=d_out)
            self.dict_df_pos[sat]=self.dict_df_pos[sat].loc[mask]

            if first:
                self.df_pos = self.dict_df_pos[sat]
                self.df_pos["sv"] = sat
                first = False
            else:
                df_inter = self.dict_df_pos[sat]
                df_inter["sv"] = sat
                self.df_pos = pd.concat([self.df_pos,df_inter])

        self.df_pos.drop_duplicates(inplace=True)


    
    def getElevation(self,df,pos_antena):

        norm_antena = np.linalg.norm(pos_antena)


        
        dfdif = pd.DataFrame() 
        if len(self.df_pos)==0: 
            df["X"] = float('NaN')
            df["Y"] = float('NaN')
            df["Z"] = float('NaN')
            df["elevation"] = float('NaN')
            return df
        dfdif['Xdif'] = self.df_pos["X"]-pos_antena[0]
        dfdif['Ydif'] = self.df_pos["Y"]-pos_antena[1]
        dfdif['Zdif'] = self.df_pos["Z"]-pos_antena[2]

        nv = pos_antena[0]*dfdif["Xdif"] + pos_antena[1]*dfdif["Ydif"] + pos_antena[2]*dfdif["Zdif"]
        norm = norm_antena * np.linalg.norm(dfdif[['Xdif','Ydif','Zdif']].values,axis=1) 
        
        sin_phi = nv/norm
        self.df_pos["elevation"] = np.arcsin(sin_phi)


        return df.merge(self.df_pos,how='left',on=['sv','time']) 


    def getPiercingPoint(self,df,pos_antena,h=350000):
        R_I = R_E + h
        df_inter = pd.DataFrame()

        #### Calculation of the position of the intersection of Satellite-Antena line with ionosphere which is a second degree equation
        df_inter["Xas"] = df["X"]-pos_antena[0]
        df_inter["Yas"] = df["Y"]-pos_antena[1]
        df_inter["Zas"] = df["Z"]-pos_antena[2]
        
        # return elevation of satellite given the ID os the satellite, the position of the antenna.
        df_inter["AS_squared"] = df_inter["Xas"]**2 + df_inter["Yas"]**2 + df_inter["Zas"]**2
        #Parameters of the polynomial to solve
        df_inter["b_poly"] = 2*df_inter["Xas"]*(pos_antena[0]*df_inter["Xas"] + pos_antena[1]*df_inter["Yas"] + pos_antena[2]*df_inter["Zas"])
        df_inter["c_poly"] = df_inter["Xas"]**2 * ( pos_antena[0]**2 + pos_antena[1]**2+ pos_antena[2]**2 - R_I**2 )
        df_inter["Delta"] = df_inter["b_poly"]**2 -4*df_inter["AS_squared"]*df_inter["c_poly"]
        
        df_inter["x_ion_1"] = (-df_inter["b_poly"]+np.sqrt(df_inter["Delta"]))/(2*df_inter["AS_squared"]) + pos_antena[0]
        df_inter["y_ion_1"] = df_inter["Yas"]/df_inter["Xas"]*(df_inter["x_ion_1"]-pos_antena[0])+pos_antena[1]
        df_inter["z_ion_1"] = df_inter["Zas"]/df_inter["Xas"]*(df_inter["x_ion_1"]-pos_antena[0])+pos_antena[2] 
        
        df_inter["x_ion_2"] = (-df_inter["b_poly"]-np.sqrt(df_inter["Delta"]))/(2*df_inter["AS_squared"]) + pos_antena[0]
        df_inter["y_ion_2"] = df_inter["Yas"]/df_inter["Xas"]*(df_inter["x_ion_2"]-pos_antena[0])+pos_antena[1]
        df_inter["z_ion_2"] = df_inter["Zas"]/df_inter["Xas"]*(df_inter["x_ion_2"]-pos_antena[0])+pos_antena[2]
        
        df_inter["X1xposAntena"] = (df_inter["x_ion_1"]-pos_antena[0])*df_inter["Xas"] + (df_inter["y_ion_1"]-pos_antena[1])*df_inter["Yas"] + (df_inter["z_ion_1"]-pos_antena[2])*df_inter["Zas"]
        df_inter["X2xposAntena"] = (df_inter["x_ion_2"]-pos_antena[0])*df_inter["Xas"] + (df_inter["y_ion_2"]-pos_antena[1])*df_inter["Yas"] + (df_inter["z_ion_2"]-pos_antena[2])*df_inter["Zas"]

        
        condition = (df_inter["X1xposAntena"]<0) | ( (df_inter["X1xposAntena"]>0) & (df_inter["X2xposAntena"]>0) &  (df_inter["X2xposAntena"]<df_inter["X1xposAntena"]))
        
        df_inter["pos_ion_X"] = np.where(condition,df_inter["x_ion_2"],df_inter["x_ion_1"]) 
        df_inter["pos_ion_Y"] = np.where(condition,df_inter["y_ion_2"],df_inter["y_ion_1"]) 
        df_inter["pos_ion_Z"] = np.where(condition,df_inter["z_ion_2"],df_inter["z_ion_1"]) 

        df["lat"],df["lon"],df["alt"] = pm.ecef2geodetic(df_inter["pos_ion_X"],df_inter["pos_ion_Y"],df_inter["pos_ion_Z"])
        df.drop(columns=['X','Y','Z'],inplace=True)
        return df


def get_arcs(elevations,t_begin=None,t_end=None):
    elevations = elevations.dropna()
    values = elevations.values
    index = elevations.index

    if t_begin==None: t_begin = index[0]
    if t_end==None: t_end = index[-1]

    # List of dictionaries containing first and last date of arc and if arc is complete in the time serie
    list_arcs=[]

    # Declare first dict, first arc.
    arc={}

    # Variable to track max elevation during the arc.
    max_elevation = 0
    iin=0

    # If first date has data, the first arc is already running
    # First date of arc is first non nan value, and arc not fully available in time serie
    if len(index)==0 or len(values)==0: return list_arcs
    ioldold,voldold = index[0],values[0]
    iold,vold = index[0],values[0]
    i,v = index[0],values[0]
    #isnan1 = math.isnan(v)
    is_in_arc = False
    if v>0:
        if len(index)>1: inew,vnew=index[1],values[1]
        else:  inew,vnew=index[0],values[0]
        arc["start"]=i
        iin=0
        max_elevation=v
        arc["end"]=None
        if len(values)>1:
            if values[1]<values[0]:
                arc["max_ele"]=values[0]
                arc["tmax"]=index[0]
                arc["imax"]=0
            else:
                arc["max_ele"]=values[1]
                arc["tmax"]=index[1]
                arc["imax"]=0

        if (i-t_begin).seconds<300 and v<15.0*math.pi/180: arc["full"]=False
        else: arc["full"]=True
        is_in_arc = True



	# Get elevation sense, rising (positive) or falling (negative)
    last_last_rising=False
    last_rising=False
    rising=False

    if len(elevations)==1:
        arc["tmax"]=index[0]
        arc["imax"]=0
        arc["max_ele"]=values[0]
    for k in range(len(elevations)-1):
		# Get new value of elevation
        inew,vnew = index[k+1],values[k+1]
        rising=vnew>v
        if not is_in_arc and vnew<0:
            ioldold,voldol=iold,vold
            iold,vold=i,v
            i,v=inew,vnew
            last_rising = rising
            continue

        if not is_in_arc and vnew>0:
            is_in_arc=True
            arc["start"]=inew
            iin=k
            arc["end"]=None
            arc["tmax"]=i
            arc["imax"]=0
            arc["full"]=True
            max_elevation = v
            arc["max_ele"]=max_elevation

        # if elevation is changing from up to down, max reached
        if k>0 and last_rising and not rising and max_elevation<vnew:
            max_elevation=vnew
            arc["tmax"]=i
            arc["imax"]=k-iin-1
            arc["max_ele"]=max_elevation

		# Case arc is finished when satellite goes down skyline
        if not rising and (vnew<0):
            arc["end"]=i
            list_arcs.append(arc)
            is_in_arc=False
            arc={} # empty dictionnary waiting for next arc

		# Case time serie between falling and next rising is truncated
        if k>1 and  not last_last_rising and rising and (i-iold>datetime.timedelta(hours=1)):
            arc["end"]=iold#inew
            list_arcs.append(arc)
            arc={}
            arc["start"]=i
            iin=k
            arc["end"]=None
            arc["tmax"]=inew
            arc["imax"]=0
            arc["full"]=True
            max_elevation = vnew
            arc["max_ele"]=max_elevation
            is_in_arc=True

        if k>0 and not last_rising and rising and (v<0):
            arc["end"]=i
            list_arcs.append(arc)
            arc={}
            arc["start"]=i
            iin=k
            arc["end"]=None
            arc["tmax"]=i
            arc["imax"]=0
            arc["full"]=True
            max_elevation = v
            arc["max_ele"]=max_elevation
            is_in_arc=True

        iold,vold=i,v
        i,v=inew,vnew

        last_last_rising = last_rising
        last_rising = rising


	# Special treatment for last data of time series, if an arc was running end it and set no "full" arc
    if is_in_arc:
        arc["end"]=inew
        if rising:
            arc["tmax"]=inew
            arc["imax"]=k-iin
        else: arc["full"]=True
        arc["max_ele"]=elevations[arc["tmax"]]
        list_arcs.append(arc)

    return list_arcs
