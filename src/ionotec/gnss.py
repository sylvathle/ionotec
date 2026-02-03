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

dir_sats = "sats/"

class gnss:
    '''
        gnss Class to perform eficiently the position and elevation of the
        satellites
        for now it only works on GPS constellation, future version should deal
        with GLONASS
    '''
    def __init__(self,f_nav=[],resolution=60):
     
        self.gps_dir = st.root_dir + "GPS/"
        if not os.path.exists(self.gps_dir):
            try: os.mkdir(self.gps_dir)
            except OSError as e:
                if e.errno!=17: print ("FAIL creation of directory "+self.gps_dir, e )
            else: print ("Successfully created the directory "+self.gps_dir)

        self.list_f_rinex_nav = f_nav
        self.dict_df_pos = {}
        self.df_pos = pd.DataFrame()
        self.resolution = resolution
        self.compute_position(self.list_f_rinex_nav)


    # Function that loads the satellite of a specific year, files must exist, otherwise load empty
    def load_sats(self,d_in,d_out):
    
        year = d_in.year
        first = True
         
        for i in range(1,33):
            #Name of the satellite
            sat = ""
            if i<10: sat = "G0"+str(i)
            else: sat = "G"+str(i)

            sat_file = self.gps_dir + str(year) + "/" + sat + ".feather"

            if os.path.exists(sat_file):
                self.dict_df_pos[sat] = pd.read_feather(sat_file)
                self.dict_df_pos[sat].set_index("time",inplace=True)
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

       
    # Function that compute the position of the satellites from the navigation file 
    #   with the resolution informed when instanciating the gnss object
    def compute_position(self,form='feather'):
    
        # Columns that need to be used
        list_cols = ["Toe","TGD","IDOT","IODC","GPSWeek","TransTime","SVclockBias","SVclockDrift",\
        "SVclockDriftRate","sqrtA","Eccentricity","Io","Omega0","omega","M0","DeltaN","OmegaDot","Cus","Cuc",\
        "Cis","Crs","Crc","Cic"] 

        for f_rinex_nav in self.list_f_rinex_nav:

            rinex_doy = int(f_rinex_nav[-8:-5])
            rinex_year = int("20"+f_rinex_nav[-3:-1])

            # Construct path of satellite files and create directory if does not exists
            year_gpsdir = self.gps_dir + str(rinex_year)
            if not os.path.exists(year_gpsdir):
                try: os.mkdir(year_gpsdir)
                except OSError as e:
                    if e.errno!=17: print ("FAIL creation of directory "+year_gpsdir, e )
                    else: print ("Successfully created the directory "+year_gpsdir)    
            
            #Annex file reporting doy that have been considered for each satellite
            f_doy_reported = year_gpsdir + "/reported_doy.feather"
            
            
            if not os.path.exists(f_doy_reported):
                dict_doy_reported = {"sat":[],"res":[]}
                for doy in range(366):
                    dict_doy_reported[str(doy)] = []
                for i in range(1,33):
                    #Name of the satellite
                    sat = ""
                    if i<10: sat = "G0"+str(i)
                    else: sat = "G"+str(i)
                    dict_doy_reported["sat"].append(sat)
                    dict_doy_reported["res"].append(self.resolution )
                    for doy in range(366):
                        dict_doy_reported[str(doy)].append(0)
                pd.DataFrame(dict_doy_reported).to_feather(f_doy_reported)      

            # Load report of already processed files in order to know which doy have been processed 
            #   for each satellite and with which resolution
            df_doy_reported = pd.read_feather(f_doy_reported)
            df_doy_reported.set_index("sat",inplace=True)
           
            # Prepare time list that will be used
            t_base = datetime.datetime(rinex_year,1,1,0,0,0,0)+datetime.timedelta(days=rinex_doy-1)
            t=t_base
            time_list = []
            while t<t_base+datetime.timedelta(days=1):
                time_list.append(t)
                t += datetime.timedelta(seconds=self.resolution)

        
            # Load provided navigation file
            try:
                nav = gr.load(f_rinex_nav)
                df_nav = nav.to_dataframe()
            except ValueError:
                print ("COuld not load file", f_rinex_nav)
                continue
                

            # Make if friendly to use with dataframe
            df_nav.index.set_names(["time","sv"],inplace=True)
            df_nav.reset_index(level=["sv"],inplace=True)
            df_nav.dropna(inplace=True)
       
            # Iter over all the satellites
            for i in range(1,33):
                #Name of the satellite
                sat = ""
                if i<10: sat = "G0"+str(i)
                else: sat = "G"+str(i)

                # flag to know if position of satellite has already been computed for this doy
                is_doy_reported = df_doy_reported[str(rinex_doy)].loc[sat]==1
            
                # flag to know if the resolution is fine
                is_res_enough = df_doy_reported["res"].loc[sat]<=self.resolution

                # If this doy has been reported with enough time-resolution, got for next satellite
                if is_doy_reported and is_res_enough: continue
            
                # import data of this satellite
                df_sat_new = df_nav[df_nav["sv"]==sat][list_cols]          
            
                # If no data has been found for this satellite, go for next one without creating the feather
                if len(df_sat_new)==0:
                    print ("No position for satellite",sat)
                    continue
            
                # File correponding to satellite of interest
                csv_nav_sat = year_gpsdir + "/" + sat + ".feather"
            
                # Intermediate dictionnary intended to contain data of the satellite under process 
                dict_sat_pos = {"time":[],"X":[],"Y":[],"Z":[]}
            
                # Iterator of time list of positions to be computed
                i_time_list = 0
                # Get first time
                date = time_list[i_time_list]
            
                # Iter over the times available in the navigation file
                for t_nav, row in df_sat_new.iterrows():
                    # Calculate position for each time in time_list that are before the next available position information
                    while date < t_nav and i_time_list<len(time_list):
                        sat_pos = gps_nav_to_XYZ(row,date)
                        dict_sat_pos["time"].append(date)
                        dict_sat_pos["X"].append(sat_pos[0])
                        dict_sat_pos["Y"].append(sat_pos[1])
                        dict_sat_pos["Z"].append(sat_pos[2])
                        i_time_list += 1
                        if i_time_list==len(time_list): break
                        date = time_list[i_time_list]
                last_row = df_sat_new.iloc[-1]
            
                # Case novigation data does not provide position until 00:00, use last available date
                while i_time_list<len(time_list):
                    sat_pos = gps_nav_to_XYZ(last_row,date)
                    dict_sat_pos["time"].append(date)
                    dict_sat_pos["X"].append(sat_pos[0])
                    dict_sat_pos["Y"].append(sat_pos[1])
                    dict_sat_pos["Z"].append(sat_pos[2])
                    i_time_list += 1
                    if i_time_list==len(time_list): break
                    date = time_list[i_time_list]
            
                df_sat_nav = pd.DataFrame(dict_sat_pos)
                df_sat_nav.set_index("time",inplace=True)
            
                # Check if feather file is not already there in order not to use navigation file
                if os.path.exists(csv_nav_sat):
                    #print (csv_nav_sat, "file exist, concat with new results")
                    df_sat_old = pd.read_feather(csv_nav_sat) 
                    df_sat_old.set_index("time",inplace=True)
                    df_sat_old.index = pd.to_datetime(df_sat_old.index)
                    df_sat_nav = pd.concat([df_sat_nav,df_sat_old])
                          
                df_sat_nav.reset_index(inplace=True)
                df_sat_nav.to_feather(csv_nav_sat)
            
                # Report that this doy has been treated with the corresponding resolution to avoid reprocessing
                df_doy_reported[str(rinex_doy)].loc[sat]=1
                df_doy_reported["res"].loc[sat]=self.resolution
        
            # Update report file
            df_doy_reported.reset_index(inplace=True)
            df_doy_reported.to_feather(f_doy_reported)
            

        
 
    
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

        #return 1.0, [20000000,2000000,-1000000]
        #pos_sat=self.getPosdf(df)
    def getElevation_old(self,sat,pos_antena,date,mode="RAD"):   
        n = np.array([pos_antena[0],pos_antena[1],pos_antena[2]])
        v = np.array([pos_sat[0]-pos_antena[0],pos_sat[1]-pos_antena[1],pos_sat[2]-pos_antena[2]])
        nv = np.dot(n,v)
        sin_phi = nv/(np.linalg.norm(n)*np.linalg.norm(v))
        return math.asin(sin_phi),pos_sat


    # Read position from saved files
    # Not working but not used
    def getPos(self,df):
        
        print (self.df_pos)
        print (df)
        print (pd.merge(df,self.df_pos,how='left',on=['sv','time']))
        sys.exit()
        
        if len(self.dict_df_pos[sat])==0: 
            return [float("NaN"),float("NaN"),float("NaN")]

        if date not in self.dict_df_pos[sat].index: 
            return [float("NaN"),float("NaN"),float("NaN")]

        pos = self.dict_df_pos[sat][["X","Y","Z"]].loc[date].tolist()
       
        return pos  

          
    #### Old function, should not be used anymore
    def to_CSV(self,date,feather=True):
        #resetCSV(date,feather)
        #print ("In resetcsv",gps_dir)
                
        year_gpsdir = self.gps_dir + str(date.year)
        if not os.path.exists(year_gpsdir):
            try: os.mkdir(year_gpsdir)
            except OSError as e:
                if e.errno!=17: print ("FAIL creation of directory "+year_gpsdir, e )
                else: print ("Successfully created the directory "+year_gpsdir)
                
        #Annex file reporting doy that have been considered for each satellite
        f_doy_reported = self.gps_dir + str(date.year) + "/reported_doy.feather"
        if not os.path.exists(f_doy_reported):
            dict_doy_reported = {"sat":[]}
            for doy in range(366):
                dict_doy_reported[str(doy)] = []
            for i in range(1,33):
                #Name of the satellite
                sat = ""
                if i<10: sat = "G0"+str(i)
                else: sat = "G"+str(i)
                dict_doy_reported["sat"].append(sat)
                
                for doy in range(366):
                    dict_doy_reported[str(doy)].append(0)
            pd.DataFrame(dict_doy_reported).to_feather(f_doy_reported)
        
        
        df_doy_reported = pd.read_feather(f_doy_reported)
        df_doy_reported.set_index("sat",inplace=True)
	
        print (df_doy_reported)

        try:
            nav = gr.load(self.f_rinex_nav)
            df_nav = nav.to_dataframe()
        except ValueError:
            print ("COuld not load file", self.f_rinex_nav)
        
        df_nav.index.set_names(["time","sv"],inplace=True)
        df_nav.reset_index(level=["sv"],inplace=True)
        df_nav.dropna(inplace=True)
         
        list_cols = ["Toe","TGD","IDOT","IODC","GPSWeek","TransTime","SVclockBias","SVclockDrift",\
        "SVclockDriftRate","sqrtA","Eccentricity","Io","Omega0","omega","M0","DeltaN","OmegaDot","Cus","Cuc",\
        "Cis","Crs","Crc","Cic"]
         
        for i in range(1,33):
            #Name of the satellite
            sat = ""
            if i<10: sat = "G0"+str(i)
            else: sat = "G"+str(i)

            is_doy_reported = df_doy_reported[str(self.doy)].loc[sat]==1

            # If this doy for this sat has been reported skip it
            if is_doy_reported: continue
                 
            # import data of this satellite
            df_sat_new = df_nav[df_nav["sv"]==sat][list_cols]
            # If other doys have been reported, we want to add the new doy to them
            if os.path.exists(year_gpsdir+"/"+sat+".feather"):
                df_sat = pd.read_feather(year_gpsdir+"/"+sat+".feather")
                df_sat.index = df_sat["time"]  
                df_sat_new = pd.concat(df_sat_new,df_sat)
            df_sat_new.reset_index(inplace=True)
            df_sat_new.to_feather(year_gpsdir+"/"+sat+".feather")
            
            # Report that this doy has been treated to no repeat
            df_doy_reported[str(self.doy)].loc[sat]=1
        df_doy_reported.reset_index(inplace=True)
        df_doy_reported.to_feather(f_doy_reported)

    #### Old function, should not be used anymore
    def getPos2(self,sat,date):
        doy = str(date.timetuple().tm_yday)
        year = str(date.year)
        
        # File correponding to satellite of interest
        csv_nav_sat = st.root_dir + "GPS/" + str(date.year) + "/" + sat + ".feather"
        
        # Check if feather file is not already there in order not to use navigation file
        if not os.path.exists(csv_nav_sat):
            print (csv_nav_sat, "does not exist")
            self.to_CSV(date)
        df_sat_nav = pd.read_feather(csv_nav_sat) 
        df_sat_nav.set_index("time",inplace=True)
        df_sat_nav.index = pd.to_datetime(df_sat_nav.index)
         
        dt=datetime.timedelta(days=2)
        closest_t = date
        
        if len(df_sat_nav)==0:
            print ("Position of sat ", sat, "unknown at",date )
            return [float("NaN"),float("NaN"),float("NaN")]
        
        for t in df_sat_nav.index:
            # Test for dates before time seeked
            #if d<t:
            if (t-date)<dt:
                closest_t = t
                dt = t-date
                
        #print (sat,t)
        data = df_sat_nav.loc[t]    

        return gps_nav_to_XYZ(data,date)    
        



    # return elevation of satellite given the ID os the satellite, the position of the antenna.
    def getElevation_old(self,sat,pos_antena,date,mode="RAD"):
        #return 1.0, [20000000,2000000,-1000000]
        pos_sat=self.getPos(sat,date)
        n = np.array([pos_antena[0],pos_antena[1],pos_antena[2]])
        v = np.array([pos_sat[0]-pos_antena[0],pos_sat[1]-pos_antena[1],pos_sat[2]-pos_antena[2]])
        nv = np.dot(n,v)
        sin_phi = nv/(np.linalg.norm(n)*np.linalg.norm(v))
        return math.asin(sin_phi),pos_sat


def get_pos_antenna(station,year,doy):
    df_antennas = pd.read_csv(st.root_dir+str(year)+"/"+str(doy)+"/antennas.csv")
    df_antennas.set_index("station",inplace=True)
    pos_antena = []
    pos_antena.append(df_antennas.loc[station,"x"])
    pos_antena.append(df_antennas.loc[station,"y"])
    pos_antena.append(df_antennas.loc[station,"z"])
    return pos_antena

def getPiercingPoint_old(pos_antena,pos_sat,h=400000):
    R_I = R_E + h
    #### Calculation of the position of the intersection of Satellite-Antena line with ionosphere which is a second degree equation
    Xas = pos_sat[0]-pos_antena[0]
    # return elevation of satellite given the ID os the satellite, the position of the antenna.
    Yas = pos_sat[1]-pos_antena[1]
    Zas = pos_sat[2]-pos_antena[2]
    AS_squared = Xas**2 + Yas**2 + Zas**2
    #Parameters of the polynomial to solve
    a_poly = AS_squared
    b_poly = 2*Xas*(pos_antena[0]*Xas + pos_antena[1]*Yas + pos_antena[2]*Zas)
    c_poly = Xas**2 * ( pos_antena[0]**2 + pos_antena[1]**2+ pos_antena[2]**2 - R_I**2 )
    Delta = b_poly**2 - 4*a_poly * c_poly
    x_ion, y_ion, z_ion = float('nan'),float('nan'),float('nan')
    pos_ion = [float('nan'),float('nan'),float('nan')]
    if Delta > 0:
        x_ion_1 = (-b_poly+np.sqrt(Delta))/(2*a_poly) + pos_antena[0]
        y_ion_1 = Yas/Xas*(x_ion_1-pos_antena[0])+pos_antena[1]
        z_ion_1 = Zas/Xas*(x_ion_1-pos_antena[0])+pos_antena[2]

        x_ion_2 = (-b_poly-np.sqrt(Delta))/(2*a_poly) + pos_antena[0]
        y_ion_2 = Yas/Xas*(x_ion_2-pos_antena[0])+pos_antena[1]
        z_ion_2 = Zas/Xas*(x_ion_2-pos_antena[0])+pos_antena[2]

        # point must be between satellite and antena, este must be the other solution
        X1xposAntena=(x_ion_1-pos_antena[0])*Xas + (y_ion_1-pos_antena[1])*Yas + (z_ion_1-pos_antena[2])*Zas
        X2xposAntena=(x_ion_2-pos_antena[0])*Xas + (y_ion_2-pos_antena[1])*Yas + (z_ion_2-pos_antena[2])*Zas
        if X1xposAntena<0 or (X1xposAntena>0 and X2xposAntena>0 and X2xposAntena<X1xposAntena): pos_ion = [x_ion_2,y_ion_2,z_ion_2]
        else: pos_ion = [x_ion_1,y_ion_1,z_ion_1]
    elif Delta==0:
        x_ion = -b_poly/(2*a_poly)
        pos_ion = [x_ion,Yas/Xas*(x_ion-pos_antena[0])+pos_antena[1],Zas/Xas*(x_ion-pos_antena[0])+pos_antena[2]]
    return pos_ion

def resetCSV(date,feather=True):
    gps_dir = st.root_dir + "GPS/" #str(date.year)+"/"+ str(date.timetuple().tm_yday) +"/GPS/"
    print ("In resetcsv",gps_dir)
    try: os.mkdir(gps_dir)
    except OSError as e:
        if e.errno!=17: print ("FAIL creation of directory "+gps_dir, e )
    else: print ("Successfully created the directory "+gps_dir)
    year_gpsdir = gps_dir + str(date.year)
    try: os.mkdir(year_gpsdir)
    except OSError as e:
        if e.errno!=17: print ("FAIL creation of directory "+year_gpsdir, e )
    else: print ("Successfully created the directory "+year_gpsdir)
    #for n_sat in range(1,33):
    #    if n_sat<10: sat = "G0" + str(n_sat)
    #    else: sat = "G" + str(n_sat)
    #    sat_dir = gps_dir+sat
    #    try: os.mkdir(sat_dir)
    #    except OSError as e:
    #        if e.errno!=17: print ("FAIL creation of directory "+gps_dir, err)
    #    else: print ("Successfully created the directory "+sat_dir)
    #    for filename in os.listdir(sat_dir):
    #        file_path = os.path.join(sat_dir, filename)
    #        try:
    #            if os.path.isfile(file_path) or os.path.islink(file_path): os.unlink(file_path)
    #            elif os.path.isdir(file_path): shutil.rmtree(file_path)
    #        except Exception as e: print('Failed to delete %s. Reason: %s' % (file_path, e))




def getGPSSatPosition(sat,d,pos_antena=None):
	#csv_nav_dir = dir_rinex_nav + "GPS/" + sat + "/"
	csv_nav_dir = st.root_dir + str(d.year)+"/"+ str(d.timetuple().tm_yday) +"/" + "GPS/" + sat + "/"
	if not os.path.isdir(csv_nav_dir):
		print (f"Calculating RINEX navigation to csv for year:{d.year}, doy:{d.timetuple().tm_yday}")
		makeCSV(d)
		print ("...Done")
	nav_csv_files = [f for f in listdir(csv_nav_dir) if isfile(join(csv_nav_dir, f))]
	if len(nav_csv_files)==0:
		print (f"Calculating RINEX navigation to csv for year:{d.year}, doy:{d.timetuple().tm_yday}")
		makeCSV(d)
		print ("...Done")
		nav_csv_files = [f for f in listdir(csv_nav_dir) if isfile(join(csv_nav_dir, f))]
	csv_nav = csv_nav_dir+nav_csv_files[0]

	df_data = pd.read_csv(csv_nav)
	df_data["time"]=pd.to_datetime(df_data["time"])
	df_data.set_index("time",inplace=True)
	dt=datetime.timedelta(days=2)
	closest_t = d
	for t in df_data.index:
		# Test for dates before time seeked
		#if d<t:
		if (t-d)<dt:
			closest_t = t
			dt = t-d

	data = df_data.loc[t]
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

	tk = (d-d0).seconds

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



def getElevation(sat,pos_antena,date,mode="RAD"):
	#return 1.0, [20000000,2000000,-1000000]
	pos_sat=getGPSSatPosition(sat,date,pos_antena)
	n = np.array([pos_antena[0],pos_antena[1],pos_antena[2]])
	v = np.array([pos_sat[0]-pos_antena[0],pos_sat[1]-pos_antena[1],pos_sat[2]-pos_antena[2]])
	nv = np.dot(n,v)
	sin_phi = nv/(np.linalg.norm(n)*np.linalg.norm(v))
	if mode=="DEG": return math.asin(sin_phi)*180.0/math.pi, pos_sat
	return math.asin(sin_phi),pos_sat

def getBias(sat,date,mode="P1P2"):
	if mode=="P1P2": filename = mode+date.strftime("%y%m")+".DCB"
	elif mode=="P1C1": filename = mode+date.strftime("%y%m")+".DCB"
	path_local_file = st.root_dir + str(date.year)+"/"+ str(date.timetuple().tm_yday) +"/" + filename
	if not os.path.exists(path_local_file):
		compressed_file = filename + ".Z"
		url_compressed_file = "http://ftp.aiub.unibe.ch/CODE/" + str(date.year) + "/" + compressed_file
		path_compressed_file = dir_sat_bias + compressed_file
		#print (path_local_file)
		try:
			print ("Download " + compressed_file)
			urllib.request.urlretrieve(url_compressed_file,path_compressed_file)
		except urllib.error.URLError as e:
			print ("Problem in retrieving: " + url_compressed_file +"\n" + e.reason)
			return 0
		with open(path_compressed_file,'rb') as fh:
			f=open(path_local_file,"w")
			f.write(unlzw(fh.read()).decode("utf-8"))
			f.close()
		os.remove(path_compressed_file)
	f = open(path_local_file)
	for lin in f:
		splt_lin = lin.split()
		if len(splt_lin)==0: continue
		if splt_lin[0]!=sat: continue
		else: return float(splt_lin[1])
	print ("No biais found")
	return 0.0

def getBias_fromfile(sat,f_bias):
    path_local_file = f_bias
    f = open(path_local_file)
    for lin in f:
         splt_lin = lin.split()
         if len(splt_lin)==0: continue
         if splt_lin[0]!=sat: continue
         else: return float(splt_lin[1])
    print ("No biais found")
    return 0.0


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
