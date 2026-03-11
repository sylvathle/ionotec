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

pd.options.mode.chained_assignment = None

#import random
#import time

#import matplotlib.pyplot as plt
from . import stations as st

from os import listdir
from os.path import isfile, join

#Ionosphere altitude
Alt_i = 300000
# Earth Radius
R_E = 6370000
omega_e = 7.292115e-5
GM = 3.986005e14
C20 = -1082.63e-6
a_major_ellipsoid_axis = 6378136
WG84_to_PZ90=np.matrix([[1,-0.33/3600.0,0],[0.33/3600.0,1,0],[0,0,1]])


#root_dir = "/media/sblunier/CSN_TEC/"
#root_dir = os.getenv("ROOT_RINEX")
#root_dir = "/home/sblunier/Work/spaceweather/tec/"
#dir_rinex_nav = "/home/sblunier/Work/spaceweather/data/GNSS/"
#dir_rinex_nav = root_dir+"2020/349/"
#dir_sat_bias = "/home/sblunier/Work/spaceweather/data/GNSS/CODE/"
#dir_sat_bias = dir_rinex_nav

def F1(pos,v):
	return np.array(v)

def F2(pos,v,a):
	r=np.linalg.norm(pos)
	xpp = -GM*pos[0]/r**3 + 1.5*C20*GM*a_major_ellipsoid_axis**2 * pos[0] * abs(1-5*pos[2]**2/r**2) / r**5  + omega_e**2*pos[0] + 2*omega_e*v[1]+a[0]
	ypp = -GM*pos[1]/r**3 + 1.5*C20*GM*a_major_ellipsoid_axis**2 * pos[1] * abs(1-5*pos[2]**2/r**2) / r**5 + omega_e**2*pos[1] - 2*omega_e*v[0]+a[1]
	zpp = -GM*pos[2]/r**3 + 1.5*C20*GM*a_major_ellipsoid_axis**2 * pos[2] * abs(3-5*pos[2]**2/r**2) / r**5
	return np.array([xpp,ypp,zpp])

def getDateFromSat(weekNumber,seconds):
	date_ref = datetime.datetime(1980,1,6,0,0,0)
	days = weekNumber*7
	sec = int(seconds) #- 604800
	while sec > 86399:
		days = days + 1
		sec -= 86400
	date = date_ref + datetime.timedelta(days=days) + datetime.timedelta(seconds=sec)
	return date

def removeOutsiders(df,thres=datetime.timedelta(days=1)):
    #### Remove outsiders
    T = df.index
    diff_t = []
    for it in range(len(T)-1):
        diff_t.append(T[it+1]-T[it])
    
    outsider=[]
    istail = False
    thres_t = datetime.timedelta(days=1)
    for d in range(len(diff_t)):
        if diff_t[d]>thres_t: 
            outsider.append(True)
        else:
            if not istail: 
                istail=True
                outsider.append(False)
            outsider.append(False)
            
    #print (outsider)    
    #print (len(outsider))
    #print (len(diff_t))

    df["out"] = outsider
    df = df[df["out"]==False]
    df.drop(columns=["out"],inplace=True)    
    return df

def gps_nav_to_XYZ(data,date):
    Toe = data['Toe']
    TGD = data['TGD']
    IDOT = data['IDOT']
    IODC = data['IODC']
    GPSWeek = data['GPSWeek']
    TransTime = data['TransTime']
    #clock correction biais
    clockBiais = data['SVclockBias']
    #clock correction drift
    clockDrift = data['SVclockDrift']
	#clock correction drift rate
    clockDriftRate = data['SVclockDriftRate']
	#Square root of the semi-major axis m^{1/2}
    sqrtA =  data['sqrtA']
	#Eccentricity
    Eccentricity =  data['Eccentricity']
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
	#Rate of change of inclination
    IDOT = data['IDOT']
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

    d0 = getDateFromSat(GPSWeek,Toe)

    tk = (date-d0).seconds

    n0 = math.sqrt(GM)/sqrtA**3
    n = n0 + DeltaN

    Mk = M0 + n*tk
    Ek = Mk
    for i in range(3): Ek = Mk + Eccentricity*math.sin(Ek)

    sin_vk = math.sqrt(1-Eccentricity**2)*math.sin(Ek) / (1.0-Eccentricity*math.cos(Ek))
    cos_vk = (math.cos(Ek)-Eccentricity)/(1.0-Eccentricity*math.cos(Ek))
    vk=0
    if sin_vk>0: vk = math.acos(cos_vk)
    else: vk = -math.acos(cos_vk)

    Phik = vk + omega
    delta_uk = Cuc * math.cos(2*Phik) + Cus * math.sin(2*Phik)
    delta_rk = Crc * math.cos(2*Phik) + Crs * math.sin(2*Phik)
    delta_ik = Cic * math.cos(2*Phik) + Cis * math.sin(2*Phik)

    uk = Phik + delta_uk
    rk = sqrtA**2 * (1 - Eccentricity * math.cos(Ek)) + delta_rk
    ik = Io + IDOT*tk + delta_ik

    Xk_p = rk * math.cos(uk)
    Yk_p = rk * math.sin(uk)

    Omegak = Omega0 + (OmegaDot-omega_e)*tk - omega_e * Toe

    Xk = Xk_p * math.cos(Omegak) - Yk_p * math.sin(Omegak)* math.cos(ik)
    Yk = Xk_p * math.sin(Omegak) + Yk_p * math.cos(Omegak)* math.cos(ik)
    Zk = Yk_p * math.sin(ik)
    return [Xk,Yk,Zk]

#dir_sats = "sats/"

class gnss:
    '''
        gnss Class to perform eficiently the position and elevation of the
        satellites
        for now it only works on GPS constellation, future version should deal
        with GLONASS
    '''

    file_format = 'feather' ## csv format not implemented

    n_sat_gps = 32
    n_sat_glonass = 24
    n_sat_galileo = 25
    resolution = 60

    list_f_rinex_nav = []

    #list_cols_gps = ["Toe","TGD","IDOT","IODC","GPSWeek","TransTime","SVclockBias","SVclockDrift",\
    #    "SVclockDriftRate","sqrtA","Eccentricity","Io","Omega0","omega","M0","DeltaN","OmegaDot","Cus","Cuc",\
    #    "Cis","Crs","Crc","Cic"] 
    #list_cols_glonass = ["sv","X","Y","Z","dX","dY","dZ","dX2","dY2","dZ2"]


    rinex_doy = 1
    rinex_year = 2020
    f_doy_reported = ""
    
    def __init__(self,f_nav=[],resolution=60,form='feather'):
        
        self.gnss_dir = st.root_dir + "GNSS/"
        
        if not os.path.exists(self.gnss_dir):
            try: os.mkdir(self.gnss_dir)
            except OSError as e:
                if e.errno!=17: print ("FAIL creation of directory "+self.gnss_dir, e )
            else: print ("Successfully created the directory "+self.gnss_dir)

        self.file_format = form
        self.list_f_rinex_nav = f_nav
        self.dict_df_pos = {}
        self.df_pos = pd.DataFrame()
        
        #self.df_doy_processed = pd.DataFrame()
        self.dict_doy_processed = {}# {"date":[],"G_res":[],"R_res":[]}
        self.csv_record_processing = self.gnss_dir+"/record_processing.csv"
        
        self.resolution = resolution
      
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
        #print (list_date)
        #print (constellation)
        #print (dates)

        for date in dates:  
            if date in self.dict_doy_processed.keys():
                if (self.resolution<self.dict_doy_processed[date][constellation+"_res"]) or\
                    (self.dict_doy_processed[date][constellation+"_res"]==0):
                    self.dict_doy_processed[date][constellation+"_res"] = self.resolution
            else:
                self.dict_doy_processed[date] = {"G_res":0,"R_res":0}
                #self.dict_doy_processed["date"].append(date)
                #self.dict_doy_processed["G_res"].append(0)
                #self.dict_doy_processed["R_res"].append(0)
                self.dict_doy_processed[date][constellation+"_res"]=self.resolution
        #print (self.dict_doy_processed)
    
    # Function that compute the position of the satellites from the navigation file 
    #   with the resolution informed when instanciating the gnss object
    def compute_position(self):

        list_nav_gps = []
        self.df_nav_glonass = pd.DataFrame()
        self.df_nav_gps = pd.DataFrame()      
        
        if os.path.exists(self.csv_record_processing):
            df_doy_processed = pd.read_csv(self.csv_record_processing,index_col="time")
            #df_doy_processed = df_doy_processed.T
            #print (df_doy_processed)
            df_doy_processed.index = pd.to_datetime(df_doy_processed.index).date
            df_doy_processed = df_doy_processed.T
            self.dict_doy_processed = df_doy_processed.to_dict()
        #self.dict_dict_processed.index = pd.to_datetime(df_doy_processed["time"]).dt.date
        #print (self.dict_doy_processed)
        #sys.exit()
        
        for f_rinex_nav in self.list_f_rinex_nav:
            #print (f_rinex_nav)
            nav = gr.load(f_rinex_nav)
            head = nav.head()
            svtype = head.svtype[0]
            
            df_nav = nav.to_dataframe()
            #print (df_nav)
            df_nav.index.set_names(["time","sv"],inplace=True)
            df_nav.reset_index(level=["sv"],inplace=True)
            df_nav.dropna(inplace=True)

            df_nav.index = pd.to_datetime(df_nav.index)
            if len(df_nav)==0: 
                #print ("----------------------",f_rinex_nav)
                continue
            #print (f_rinex_nav)
            
            #print (df_nav)

            date = df_nav.index[int(len(df_nav)/2)].date()
            #print ("date",f_rinex_nav,date)

            
            df_nav.sort_index(ascending=True,inplace=True)
            #print (df_nav)
            #filename = head.filename
            #df_nav["sta"]  = filename[:4]
            if svtype=='R':
                if date in self.dict_doy_processed.keys():
                    #if self.dict_doy_processed[date]["R_res"]==0: continue
                    if self.dict_doy_processed[date]["R_res"]>=self.resolution: continue
                self.df_nav_glonass = pd.concat([self.df_nav_glonass,df_nav[["sv","X","Y","Z"]]])
            elif svtype=='G':
                if date in self.dict_doy_processed.keys():
                    #if self.dict_doy_processed[date]["G_res"]==0: continue
                    if self.dict_doy_processed[date]["G_res"]>=self.resolution: continue
                self.df_nav_gps = pd.concat([self.df_nav_gps,df_nav])

        self.compute_gps_pos()
        self.compute_glonass_pos()

        df_doy_processed = pd.DataFrame(self.dict_doy_processed)
        df_doy_processed = df_doy_processed.transpose()
        #print (df_doy_processed)
        df_doy_processed.index.set_names(["time"],inplace=True)
        df_doy_processed.sort_index(inplace=True)
        #print (df_doy_processed)
        df_doy_processed.to_csv(self.csv_record_processing,index=True)

    def compute_gps_pos(self):
        if (len(self.df_nav_gps)<=2): return

        print ("GPS POS")
        #self.df_nav_gps.index.set_names(["time","sv"],inplace=True)
        #self.df_nav_gps.reset_index(level=["sv"],inplace=True)
        #self.df_nav_gps.dropna(inplace=True)

        #self.df_nav_gps.index = pd.to_datetime(self.df_nav_gps.index)
        #self.df_nav_gps.sort_index(ascending=True,inplace=True)
        
        list_sv = self.df_nav_gps["sv"].unique().tolist()
        self.df_nav_gps = self.df_nav_gps.groupby(["time","sv"]).mean()
        self.df_nav_gps.reset_index(level=["sv"],inplace=True)
        list_times = []#self.df_nav_gps.index
        
        for sv in list_sv:
            df_sat = self.df_nav_gps[self.df_nav_gps["sv"]==sv]
            if len(df_sat)<3: continue
            #df_sat_new = df_sat_new.groupby(["time","sv"]).mean()
            #df_sat_new.reset_index(level=["sv"],inplace=True)
            df_sat = removeOutsiders(df_sat)
            list_times = list_times+df_sat.index.tolist()

            # If no data has been found for this satellite, go for next one without creating the feather
            #if len(df_sat)==0:
            #    print ("No position for satellite",sv)
            #    continue
                
            min_date = min(df_sat.index)
            max_date = max(df_sat.index)
            
            gps_time_list = []
            t = min_date
                        
            while t<=max_date:
                gps_time_list.append(t)
                t = t+datetime.timedelta(seconds=self.resolution)

            #print (self.gps_time_list)
            i_time_list = 0
            date = gps_time_list[i_time_list]
            # Intermediate dictionnary intended to contain data of the satellite under process 
            dict_sat_pos = {"time":[],"X":[],"Y":[],"Z":[]}
            for t_nav, row in df_sat.iterrows():
                # Calculate position for each time in time_list that are before the next available position information
                while date < t_nav and i_time_list<len(gps_time_list):
                    sat_pos = gps_nav_to_XYZ(row,date)
                    dict_sat_pos["time"].append(date)
                    dict_sat_pos["X"].append(sat_pos[0])
                    dict_sat_pos["Y"].append(sat_pos[1])
                    dict_sat_pos["Z"].append(sat_pos[2])
                    i_time_list += 1
                    if i_time_list==len(gps_time_list): break
                    date = gps_time_list[i_time_list]

            # Case navigation data does not provide position until 00:00, use last available date 
            while i_time_list<len(gps_time_list):
                sat_pos = gps_nav_to_XYZ(row,date)
                dict_sat_pos["time"].append(date)
                dict_sat_pos["X"].append(sat_pos[0])
                dict_sat_pos["Y"].append(sat_pos[1])
                dict_sat_pos["Z"].append(sat_pos[2])
                i_time_list += 1
                if i_time_list==len(gps_time_list): break
                date = gps_time_list[i_time_list]

            df_sat = pd.DataFrame(dict_sat_pos)
            df_sat["time"] = pd.to_datetime(df_sat["time"])
            df_sat.set_index("time",inplace=True)
            #df_sat.reset_index(inplace=True)
            #print (df_sat)
            list_years = list({d.year for d in df_sat.index.tolist()})
            for y in list_years:
                df_year = df_sat[df_sat.index.year == y].dropna()
                csv_nav_sat = self.gnss_dir + str(y) + "/" + sv + ".feather"
                if not os.path.exists(self.gnss_dir + str(y)): os.mkdir(self.gnss_dir + str(y))
                #df_sat.to_feather(csv_nav_sat)
                #csv_nav_sat = self.gnss_dir + str(y) + "/" + sv + ".csv"
                #print (csv_nav_sat)
                if os.path.exists(csv_nav_sat):
                    df_former_sat = pd.read_feather(csv_nav_sat)
                    df_sat = pd.concat([df_sat,df_former_sat])
                    #print (df_sat)
                df_sat.to_feather(csv_nav_sat)

        list_times = sorted(list(set(list_times)))
        #list_times = sorted(list_times)
        self.inform_date_processed(list_times,'G')

    
    def compute_glonass_pos(self):
        if (len(self.df_nav_glonass)<=2): return

        print ("GLONASS POS")
        #self.df_nav_glonass.index.set_names(["time","sv"],inplace=True)
        #self.df_nav_glonass.reset_index(level=["sv"],inplace=True)
        #self.df_nav_glonass.dropna(inplace=True)
        
        #self.df_nav_glonass.index = pd.to_datetime(self.df_nav_glonass.index)
        #self.df_nav_glonass.sort_index(ascending=True,inplace=True)

        #min_date = min(self.df_nav_glonass.index)
        #max_date = max(self.df_nav_glonass.index)

        #self.glonass_time_list = []
        #t = min_date
        #while t<=max_date:
        #    self.glonass_time_list.append(t)
        #    t = t+datetime.timedelta(seconds=self.resolution)
        #print (self.glonass_time_list)
        #self.df_nav_glonass = self.df_nav_glonass[["sv","X","Y","Z","dX","dY","dZ","dX2","dY2","dZ2"]]#.iloc[8:38]
        #self.df_nav_glonass = self.df_nav_glonass[["sv","X","Y","Z"]]#.iloc[8:38]
        #pd.set_option('display.precision', 15)
        self.df_nav_glonass = self.df_nav_glonass.groupby(["time","sv"]).mean()
        self.df_nav_glonass.reset_index(level=["sv"],inplace=True)
        list_times = self.df_nav_glonass.index
        
        list_sv = self.df_nav_glonass["sv"].unique().tolist()
        for sv in list_sv:
            df_sat = self.df_nav_glonass[self.df_nav_glonass["sv"]==sv]
            if len(df_sat)<3: continue

            #print (df_sat)
            df_sat = removeOutsiders(df_sat)
            if len(df_sat)<3: continue

            min_date = min(df_sat.index)
            max_date = max(df_sat.index)
            glonass_time_list = []
            t = min_date
            while t<=max_date:
                glonass_time_list.append(t)
                t = t+datetime.timedelta(seconds=self.resolution)            
            
            new_index = df_sat.index.union(glonass_time_list)
            df_interp = df_sat.reindex(new_index).sort_index()
            #print (self.time_list)
            df_interp = df_interp[['X','Y','Z']].interpolate(method='quadratic',axis=0)
            list_years = list({d.year for d in df_interp.index.tolist()})
            for y in list_years:
                df_year = df_interp[df_interp.index.year == y].dropna()
                csv_nav_sat = self.gnss_dir + str(y) + "/" + sv + ".feather"
                if not os.path.exists(self.gnss_dir + str(y)): os.mkdir(self.gnss_dir + str(y))
                if os.path.exists(csv_nav_sat):
                    df_former_sat = pd.read_feather(csv_nav_sat)
                    df_year = pd.concat([df_year,df_former_sat])
                df_year.to_feather(csv_nav_sat)  
                
        self.inform_date_processed(list_times,'R')

        
            
    # Function that loads the satellite of a specific year, files must exist, otherwise load empty
    def load_sats(self,d_in,d_out):

        '''
            Constellation: G: GPS
                            R: GLONASS
                            E: GALILEO
        '''
    
        year = d_in.year
        first = True

        nsat = 36
         
        for n_sat,g in [(self.n_sat_gps,"G"),(self.n_sat_glonass,"R"),(self.n_sat_galileo,"E")]:
            for i in range(1,n_sat+1):
                #Name of the satellite
                sat = ""
                if i<10: sat = g+"0"+str(i)
                else: sat = g+str(i)
    
                sat_file = self.gnss_dir + str(year) + "/" + sat + ".feather"
                #sat_file = self.gps_dir + str(year) + "/" + sat + ".csv"
    
                if os.path.exists(sat_file):
                    self.dict_df_pos[sat] = pd.read_feather(sat_file)
                    #self.dict_df_pos[sat] = pd.read_csv(sat_file)
                    #print (sat)
                    #print (self.dict_df_pos[sat])
                    #self.dict_df_pos[sat].set_index("time",inplace=True)
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
                else:
                    self.dict_df_pos[sat] = pd.DataFrame()

    
    def getElevation(self,df,pos_antena):

        norm_antena = np.linalg.norm(pos_antena)
        
        dfdif = pd.DataFrame()   
        dfdif['Xdif'] = self.df_pos["X"]-pos_antena[0]
        dfdif['Ydif'] = self.df_pos["Y"]-pos_antena[1]
        dfdif['Zdif'] = self.df_pos["Z"]-pos_antena[2]

        nv = pos_antena[0]*dfdif["Xdif"] + pos_antena[1]*dfdif["Ydif"] + pos_antena[2]*dfdif["Zdif"]
        norm = norm_antena * np.linalg.norm(dfdif[['Xdif','Ydif','Zdif']].values,axis=1) 
        
        sin_phi = nv/norm
        self.df_pos["elevation"] = np.arcsin(sin_phi)
        return pd.merge(df,self.df_pos,how='left',on=['sv','time']) 


    def getPiercingPoint(self,df,pos_antena,h=400000):
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
        return df

'''
def getBias_fromfile(sat,f_bias):
    path_local_file = f_bias
    f = open(path_local_file)
    for lin in f:
         splt_lin = lin.split()
         if len(splt_lin)==0: continue
         if splt_lin[0]!=sat: continue
         else: return float(splt_lin[1])
    print ("No bias found for "+sat)
    return 0.0
'''

def getBias_fromfile(sat,file):
    fbias = open(file,'r')
    if sat[0]=="G":
        for line in fbias.readlines():
            splt_line = line.split()
            if splt_line[0] == "DSB":
                if splt_line[2]==sat:
                    if splt_line[3]=="C1W":
                        if splt_line[4]=="C2W":
                            return float(splt_line[8])
    if sat[0]=="R":
        for line in fbias.readlines():
            splt_line = line.split()
            if splt_line[0] == "DSB":
                if splt_line[2]==sat:
                    if splt_line[3]=="C1P":
                        if splt_line[4]=="C2P":
                            return float(splt_line[8])
    return float('NaN')

'''
	return list of dictionnaries, each correspond to information of the arcs in he time series
	each arc is described with keys:
		- "start" : First point of satellite arc
		- "end" : Last point of satellite arc
		- "tmax": time of maximum elevation
		- "full": True if begining and end of arc is contained is the time serie, if truncated, False
	Input: pandas Serie of elevations with time index.
		t_begin and t_end are to indicate time borders of the time serie.

'''

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

        if i-t_begin<datetime.timedelta(minutes=5) and v<15.0*math.pi/180: arc["full"]=False
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
		#print (inew,vnew,last_rising,rising,is_in_arc)
		#isnan2 = math.isnan(vnew)
        if not is_in_arc and vnew<0:
            ioldold,voldol=iold,vold
            iold,vold=i,v
            i,v=inew,vnew
            last_rising = rising
            continue

        if not is_in_arc and vnew>0:
            #print ("not in arc and vnew>0")
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
            #print ("going down, skyline passed")
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
			#print ("GAP")
			#print (i,iold,i-iold)
			#print (list_arcs)

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
			#print ("CHANGING ARC v<0")
			#print (i,iold,i-iold)
			#print (list_arcs)

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

	#for arc in list_arcs: print (arc)
    return list_arcs
