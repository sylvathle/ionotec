import time
import datetime
import sys
import pandas as pd
import georinex as gr
import pymap3d as pm

def merge_str(a,b):
    new_str = ""
    i = 0
    while i < len(a):
        if i < len(b) and b[i] not in [" ","\n"]:
            new_str += b[i]
        else:
            new_str += a[i]
        i += 1
    while i < len(b):
        new_str += b[i]
        i += 1
        
    return new_str


def process_head_data_line(line):
    
    spline = line.split()
    year = int(spline[1])
    month = int(spline[2])
    day = int(spline[3])
    hour = int(spline[4])
    minute = int(spline[5])
    second = int(spline[6].split(".")[0])
    millisecond = int(spline[6][3:6])
    
    d = datetime.datetime(year,month,day,hour,minute,second,millisecond)
    
    nsat = int(spline[8])

    ### Get Satellites list
    istr = 0
    str_sats=spline[9]
    str_sats = line[41:].replace(" ","0")
    list_sats = []
    while istr<len(str_sats):
        list_sats.append(str_sats[istr:istr+3])
        istr=istr+3

    return d,nsat,list_sats

def getVal(strVal):
    #float(self.line[il:il+19].replace("D","E"))
    str_  = strVal.replace(" ","").replace("D","E").replace("\n","")
    try: return float(str_)
    except: return float('NaN')

def process_head_data_line_rinexv2(line):

    spline = line.split()
    year = int("20"+line[1:3])
    month = int(line[4:6])

    day = int(line[7:9])
    hour = int(line[10:12])
    minute = int(line[13:15])
    second = int(line[16:18])
    millisecond = int(line[19:22])
    d = datetime.datetime(year,month,day,hour,minute,second,millisecond)

    nsat = int(line[30:32])

    ### Get Satellites list
    istr = 0
    str_sats = line[32:]
    list_sats = []
    while istr<len(str_sats):
        list_sats.append(str_sats[istr:istr+3])
        istr=istr+3

    return d,nsat,list_sats



def update_hatanaka(list_item,sv,iobs,value):

    i = list_item[sv][iobs][-1]
    if i < 3: list_item[sv][iobs][-1] += 1
    list_item[sv][iobs][i] = value
    
    while i>0:
        i -= 1
        list_item[sv][iobs][i] += list_item[sv][iobs][i+1]


def reset_hatanaka(list_item):

    for sv in list_item.keys():
        for iobs in range(len(list_item[sv])):
            list_item[sv][iobs] = [float('NaN'),float('NaN'),float('NaN'),float('NaN'),0]
    


class rinex:

    def __init__(self,rinex_path_):
        
        self.rinex_path = rinex_path_
        self.file = open(self.rinex_path, 'r')
        self.dict_obs = {}
        self.header = {}
        self.list_df = {}
        self.def_nav = pd.DataFrame()

        self.line = ""
        self.isHeader = True
        
        self.dict_obs_name = {}
        self.dict_obs_shift = {}
        
        self.vars_list = ['P2','P1','C2','C1','L1','L2']
        
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file.close()

    def read_header(self):

        constellation = ''
        
        n_obs=0
        self.header['type_rinex'] = "N/A" 
        nline = 0
        while self.isHeader:
            self.line = self.file.readline()
            nline+=1
            label = self.line[60:]
            data = self.line[:60]
            if 'SYS / # / OBS TYPES' in label:
                if data[0]!=' ':
                    constellation = data[0]
                    n_obs = int(data[4:6])
                    self.dict_obs[constellation] = []
                    #{"time":[],"sv":[]}
                    self.dict_obs_name[constellation] = []
                    split_data = data[6:].split()
                    for obs in split_data:
                        #dict_obs[constellation][obs] = []
                        self.dict_obs_name[constellation].append(obs)
                        self.dict_obs_shift[constellation] = [0 for n in range(n_obs)]
                    continue
                else:
                    split_data = data[6:].split()
                    for obs in split_data:
                        #dict_obs[constellation][obs] = []
                        self.dict_obs_name[constellation].append(obs)
                        self.dict_obs_shift[constellation] = [0 for n in range(n_obs)]     
                    continue
            ## Case for Rinex v2
            if '# / TYPES OF OBSERV' in label:
                constellation = 'all'
                if data[4:6]!='  ':
                    n_obs = int(data[4:6])
                    self.dict_obs[constellation] = []
                    split_data = data[6:].split()
                    self.dict_obs_name[constellation] = []
                    for obs in split_data:
                        #dict_obs[constellation][obs] = []
                        self.dict_obs_name[constellation].append(obs)
                        self.dict_obs_shift[constellation] = [0 for n in range(n_obs)]
                    continue
                else:
                    split_data = data[6:].split()
                    for obs in split_data:
                        #dict_obs[constellation][obs] = []
                        self.dict_obs_name[constellation].append(obs)
                        self.dict_obs_shift[constellation] = [0 for n in range(n_obs)]  
                    continue
            if self.line[0]==">": 
                self.isHeader = False
                continue
            if 'END OF HEADER' in label:
                self.isHeader = False
                self.line = self.file.readline()
                continue
            if 'APPROX POSITION XYZ' in label:
                str_pos_antena = data.split()
                X = float(str_pos_antena[0])
                Y = float(str_pos_antena[1])
                Z = float(str_pos_antena[2])
                self.header["position"] = [X,Y,Z]
                lat,lon,alt = pm.ecef2geodetic(X,Y,Z)
                self.header["coord"] = [float(lat),float(lon),float(alt)]
                continue
            if 'MARKER NAME' in label:
                self.header["name_station"] = data.replace(" ","")[:4].lower()
                #self.header["name_station"] = self.header["name_station"][:4]
                continue
            if 'RINEX VERSION / TYPE' in label:
                spline = data.split()
                self.header["version"] = float(spline[0])
                if data[20]=='H': 
                    self.header['type']='N'
                    self.header['constellation'] = 'S'
                elif (data[20]=='N') and ('GPS' in data): 
                    self.header['type']='N'
                    self.header['constellation'] = 'G'
                elif data[20]=='G': 
                    self.header['type']='N'
                    self.header['constellation'] = 'R'
                else:
                    self.header['type']=data[20] 
                    self.header['constellation'] = data[40]
                continue
            if 'CRINEX VERS   / TYPE' in label:
                crinex_version = data[0:20].replace(" ","")
                type_rinex = data[20:40]
                if "COMPACT" in type_rinex:
                    self.header['type_rinex'] = type_rinex
            if 'INTERVAL' in label:
                self.header['resolution'] = float(data.split()[0]) 
                continue
            if 'TIME OF FIRST OBS' in label:
                dlist = data.split()
                y = int(dlist[0])
                m = int(dlist[1])
                d = int(dlist[2])
                H = int(dlist[3])
                M = int(dlist[4])
                S = int(dlist[5].split('.')[0])
                d = datetime.datetime(y,m,d,H,M,S,0)
                self.header['t_first_obs'] = d
                continue
        return self.header

    def read_data(self):

        if int(self.header["version"])==2:
            if "COMPACT" in self.header['type_rinex']:
                #self.read_obs_rinexv2()
                for constellation in ['G','R']:
                    obs = gr.load(self.rinex_path, meas=self.vars_list, use=constellation)
                    df_temp = obs.to_dataframe()
                    df_temp.reset_index(inplace=True)
                    self.list_df[constellation] = df_temp                    
               # for const in self.list_df.keys():
                #    print (self.list_df[const])
                #sys.exit()
                return self.list_df
            if self.header['type']=='O':     
                for constellation in ['G','R']:
                #for constellation in ['G']:
                    obs = gr.load(self.rinex_path, meas=self.vars_list, use=constellation)
                    df_temp = obs.to_dataframe()
                    df_temp.reset_index(inplace=True)
                    self.list_df[constellation] = df_temp
                return self.list_df
            if self.header['type']=='N':
                if self.header['constellation']=='S':
                    self.read_nav()
                    return self.df_nav
                else:
                    nav = gr.load(self.rinex_path)
                    self.df_nav = nav.to_dataframe()
                    #self.df_nav.reset_index(inplace=True)
                    self.df_nav.dropna(inplace=True)
                    #self.df_nav.set_index('time',inplace=True)
                    #self.df_nav.index = pd.to_datetime(self.df_nav.index)
                    return self.df_nav
               

        elif int(self.header["version"])==3: 
            if self.header['type']=='O':
                self.read_obs()
                return self.list_df
            if self.header['type']=='N':
                #nav = gr.load(self.rinex_path)
                #self.df_nav = nav.to_dataframe()
                if self.header['constellation'] in ['G','J','C','E']:
                    nav = self.read_nav3_ephemerids()
                else:
                    nav = self.read_nav3_xyz()
                return self.df_nav
                #nav = self.read_nav4()
                #return self.df_nav

        elif (int(self.header["version"])==4):
            if (self.header['type']=='N'):
                if self.header['constellation'] in ['G','J','C','E']:
                    nav = self.read_nav4_ephemerids()
                else:
                    nav = self.read_nav4_xyz()
                return self.df_nav
            elif "COMPACT" in self.header['type_rinex']:
                #print (self.header)
                self.read_obs()
                return self.list_df

        
        else:
            print ("Warning, file",self.rinex_path)
            print ("Rinex version is: ", self.header["version"], " which is currently no supported")
            return pd.DataFrame()
            
    
    def read_obs(self):
        
        data_head_line = ""
        
        ### First datatime
        data_head_line = self.line    
        d,nsat,list_sats = process_head_data_line(data_head_line)

        
        last_item = {}
        
        self.file.readline()
        isat = 0
        while isat<nsat:
            self.line = self.file.readline()
            spline = self.line.split(" ")
            sv = list_sats[isat]
            constellation = sv[0]
            last_item[sv]=[]
            obs_item = {"time":d, "sv":sv}
            iobs = 0
            for s in spline[:-1]:
                obs_name = self.dict_obs_name[constellation][iobs]
                if s=='': 
                    v=float('NaN')
                    obs_item[obs_name] = v
                    last_item[sv].append([v,float('NaN'),float('NaN'),float('NaN'),0])
                else: 
                    ss = s.split('&')
                    shift = int(ss[0])
                    self.dict_obs_shift[constellation][iobs] = shift
                    v = float(ss[-1])/10**shift
                    obs_item[obs_name] = v
                    last_item[sv].append([v,float('NaN'),float('NaN'),float('NaN'),1])
                iobs = iobs+1
            isat=isat+1
            
            self.dict_obs[constellation].append(obs_item)

        
        while self.line:
            self.line = self.file.readline() 
            
            ## Do not support when hatanaka is reseted
            ##   See /bodega1/TEC/RINEX/IGS/2024/130/24d/MAO000USA_R_20241300000_01D_30S_MO.crx
            data_head_line = merge_str(data_head_line,self.line.replace("&","0"))
            d,nsat,list_sats = process_head_data_line(data_head_line)

            #print (d)
            #print (self.line.replace("\n",""))
            #if self.line[0]=='>': break
        
            
            #line = file.readline()
            if not self.file.readline(): 
                break
            if self.line[0]=='>': 
                print ("LEAVE FOR FILE for restarting hatanaka '>' at beginning of file")
                #self.file.close()
                break

            isat = 0
            while isat<nsat:
                self.line = self.file.readline()
                spline = self.line[:-1].split(" ")
                sv = list_sats[isat]
                constellation = sv[0]
                obs_item = {"time":d, "sv":sv}
                iobs = 0
    
                
                if sv not in last_item.keys(): 
                    last_item[sv]=[]
        
                    while iobs<len(self.dict_obs_name[constellation]):
                        obs_name = self.dict_obs_name[constellation][iobs]
                        if iobs<len(spline): s = spline[iobs]
                        ### Case we see fewer variables reported than those announced: consider the value is NaN
                        else:
                            obs_item[obs_name] = float('NaN')
                            iobs = iobs+1
                            last_item[sv].append([float('NaN'),float('NaN'),float('NaN'),float('NaN'),0])
                            continue
                        if s=='': 
                            v=float('NaN')
                            last_item[sv].append([v,float('NaN'),float('NaN'),float('NaN'),0])
                            obs_item[obs_name] = v
                        else: 
                            ss = s.split('&')
                            str_v = ss[-1]
                            # Data is new
                            if len(ss)==2:
                                shift = int(ss[0])
                                self.dict_obs_shift[constellation][iobs] = shift
                                v = float(str_v)/10**shift
                                obs_item[obs_name] = v
                                last_item[sv].append([v,float('NaN'),float('NaN'),float('NaN'),1])
                            # Data is known
                            elif len(ss)==1:
                                shift = self.dict_obs_shift[constellation][iobs]
                                v = float(str_v)/10**shift
                                update_hatanaka(last_item,sv,iobs,v)                
                                obs_item[obs_name] = last_item[sv][iobs][0]
                        iobs = iobs+1
                    isat=isat+1
                    self.dict_obs[constellation].append(obs_item)
    
                else:
                    
                    while iobs<len(self.dict_obs_name[constellation]):
                        obs_name = self.dict_obs_name[constellation][iobs]
                        if iobs<len(spline): s = spline[iobs]
                        ### Case we see fewer variables reported than those announced: consider the value is NaN
                        else:
                            obs_item[obs_name] = float('NaN')
                            last_item[sv][iobs]=[float('NaN'),float('NaN'),float('NaN'),float('NaN'),0]
                            iobs = iobs+1
                            continue
                        if s=='': 
                            v=float('NaN')
                            last_item[sv][iobs]=[v,float('NaN'),float('NaN'),float('NaN'),0]
                            obs_item[obs_name] = v
                        else: 
                            ss = s.split('&')
                            str_v = ss[-1]
                            # Data is new
                            if len(ss)==2:
                                shift = int(ss[0])
                                v = float(str_v)/10**shift
                                obs_item[obs_name] = v
                                last_item[sv][iobs]=[v,float('NaN'),float('NaN'),float('NaN'),1]
                            # Data is known
                            elif len(ss)==1:
                                shift = self.dict_obs_shift[constellation][iobs]
                                v = float(str_v)/10**shift
                                update_hatanaka(last_item,sv,iobs,v)                
                                obs_item[obs_name] = last_item[sv][iobs][0]
                        iobs = iobs+1
                    isat=isat+1
                    self.dict_obs[constellation].append(obs_item)
        
        #self.file.readline()

        
        for constellation, observations in self.dict_obs.items():
            
            if len(observations)==0: continue
            ## Discard Indian and SBAS constellations
                # Indian constellation because only single frequency
            if constellation in ['I']: continue
            #if not constellation in ['J']: continue
            self.list_df[constellation] = pd.DataFrame(observations)
        
    def read_obs_rinexv2(self):

        
        ### First datatime
        data_head_line = self.line   
        d,nsat,list_sats = process_head_data_line_rinexv2(data_head_line)

        
         # ------------------------------------------------------------------ header
        #obs_types: list[str] = []
        #total_obs: int = 0
        #i = 0       
        #sys.exit()
        
        
        last_item = {}
        constellation = 'all'
        nline = 0
        self.line = self.file.readline()
        nline+=1
        isat = 0
        while isat<nsat:
            self.line = self.file.readline()
            self.line = self.line.replace("\n","")

            nline+=1
            spline = self.line.split(" ")
            sv = list_sats[isat]
            
            last_item[sv]=[]
            obs_item = {"time":d, "sv":sv}
            iobs = 0
            n_empty = 0

            while iobs<len(self.dict_obs_name[constellation]):
                obs_name = self.dict_obs_name[constellation][iobs]
                if iobs<len(spline): s = spline[iobs]
                ### Case we see fewer variables reported than those announced: consider the value is NaN
                else:
                    obs_item[obs_name] = float('NaN')
                    iobs = iobs+1
                    last_item[sv].append([float('NaN'),float('NaN'),float('NaN'),float('NaN'),0])
                    continue
                if s=='': 
                    v=float('NaN')
                    last_item[sv].append([v,float('NaN'),float('NaN'),float('NaN'),0])
                    obs_item[obs_name] = v
                else: 
                    ss = s.split('&')
                    str_v = ss[-1]
                    # Data is new
                    if len(ss)==2:
                        shift = int(ss[0])
                        v = float(str_v)/10**shift
                        obs_item[obs_name] = v
                        last_item[sv].append([v,float('NaN'),float('NaN'),float('NaN'),1])
                    # Data is known
                    else:
                        shift = self.dict_obs_shift[constellation][iobs]
                        v = float(str_v)/10**shift
                        update_hatanaka(last_item,sv,iobs,v)                
                        obs_item[obs_name] = last_item[sv][iobs][0]
                iobs = iobs+1
            isat=isat+1
            self.dict_obs[constellation].append(obs_item)

            
        
        while self.line:
            
            self.line = self.file.readline()

            nline+=1
            data_head_line = merge_str(data_head_line,self.line.replace("&","0"))
            d,nsat,list_sats = process_head_data_line_rinexv2(data_head_line)

            
            self.line = self.file.readline()
            if "COMMENT" in self.line:
                while "COMMENT" in self.line:
                    self.line = self.file.readline()
                    continue
                #data_head_line = merge_str(data_head_line,self.line.replace("&","0"))
                data_head_line = self.line 
                d,nsat,list_sats = process_head_data_line_rinexv2(data_head_line)
                self.line = self.file.readline()

            if not self.line:
                nline+=1
                continue
            
            
            isat = 0
            while isat<nsat:
                self.line = self.file.readline()

                nline+=1
                spline = self.line[:-1].split(" ")
                sv = list_sats[isat]
                constellation = 'all'
                obs_item = {"time":d, "sv":sv}
                iobs = 0
    
                
                if sv not in last_item.keys(): 
                    last_item[sv]=[]
        
                    while iobs<len(self.dict_obs_name[constellation]):
                        obs_name = self.dict_obs_name[constellation][iobs]
                        if iobs<len(spline): s = spline[iobs]
                        ### Case we see fewer variables reported than those announced: consider the value is NaN
                        else:
                            obs_item[obs_name] = float('NaN')
                            iobs = iobs+1
                            last_item[sv].append([float('NaN'),float('NaN'),float('NaN'),float('NaN'),0])
                            continue
                        if s=='': 
                            v=float('NaN')
                            last_item[sv].append([v,float('NaN'),float('NaN'),float('NaN'),0])
                            obs_item[obs_name] = v
                        else: 
                            ss = s.split('&')
                            str_v = ss[-1]
                            # Data is new
                            if len(ss)==2:
                                shift = int(ss[0])
                                v = float(str_v)/10**shift
                                obs_item[obs_name] = v
                                last_item[sv].append([v,float('NaN'),float('NaN'),float('NaN'),1])
                            # Data is known
                            else:
                                shift = self.dict_obs_shift[constellation][iobs]
                                v = float(str_v)/10**shift
                                update_hatanaka(last_item,sv,iobs,v)                
                                obs_item[obs_name] = last_item[sv][iobs][0]
                        iobs = iobs+1
                    isat=isat+1
                    self.dict_obs[constellation].append(obs_item)
    
                else:

                    while iobs<len(self.dict_obs_name[constellation]):
                        obs_name = self.dict_obs_name[constellation][iobs]
                        if iobs<len(spline): s = spline[iobs]
                        ### Case we see fewer variables reported than those announced: consider the value is NaN
                        else:
                            obs_item[obs_name] = float('NaN')
                            last_item[sv][iobs]=[float('NaN'),float('NaN'),float('NaN'),float('NaN'),0]
                            iobs = iobs+1
                            continue
                        if s=='': 
                            v=float('NaN')
                            last_item[sv][iobs]=[v,float('NaN'),float('NaN'),float('NaN'),0]
                            obs_item[obs_name] = v
                        else: 
                            ss = s.split('&')
                            str_v = ss[-1]
                            # Data is new
                            if len(ss)==2:
                                shift = int(ss[0])
                                v = float(str_v)/10**shift
                                obs_item[obs_name] = v
                                last_item[sv][iobs]=[v,float('NaN'),float('NaN'),float('NaN'),1]
                            # Data is known
                            else:
                                shift = self.dict_obs_shift[constellation][iobs]
                                v = float(str_v)/10**shift
                                update_hatanaka(last_item,sv,iobs,v)                
                                obs_item[obs_name] = last_item[sv][iobs][0]
                        iobs = iobs+1
                    isat=isat+1
                    self.dict_obs[constellation].append(obs_item)
        

        nline+=1

        constellations = set()

        # 2. Iterate through the observations in 'all'
        for obs_item in self.dict_obs['all']:
            # Get the 'sv' value (e.g., 'G01', 'R12', 'E05')
            sv_str = obs_item.get('sv')
            constellation = sv_str[0]  # Extracts 'G', 'R', 'E', etc.
            
            # Ensure the key exists in our dictionary structure
            if constellation not in self.dict_obs:
                self.dict_obs[constellation] = []
                
            # Append the observation item to its respective constellation list
            self.dict_obs[constellation].append(obs_item)


        for constellation, observations in self.dict_obs.items():
            if constellation in ['I']: continue
            self.list_df[constellation] = pd.DataFrame(observations)

        
        if "I" in self.list_df.keys():
            del self.list_df["I"]  
        if "all" in self.list_df.keys():
            del self.list_df["all"]   

    def read_nav3_ephemerids(self):

        i = 0
        dict_nav_data = {
            "sv":[], "time":[], "SVclockBias":[], "SVclockDrift": [], "SVclockDrift2":[],
            "IODE":[], "Crs":[], "DeltaN":[], "M0":[],
            "Cuc":[], "Eccentricity":[], "Cus":[], "sqrtA":[],
            "Toe":[], "Cic":[], "Omega0":[], "Cis":[],
            "Io":[], "Crc":[], "omega":[], "OmegaDot":[],
            "IDOT":[], "L2":[], "Week":[], "L2Pflag":[],
            "svAcc":[], "health":[], "TGD":[], "IODC":[],
            "transTime":[], "BNK":[]
        }

        
        while self.line:
            
            i += 1    
            if self.line[0]!=" ":
                #self.line = self.file.readline()
                sv=self.line[:3].replace(" ","0")
                y = int(self.line[4:8])
                m = int(self.line[9:11])
                d = int(self.line[12:14])
                H = int(self.line[15:17])
                M = int(self.line[18:20])
                S = int(self.line[21:23])
                d = datetime.datetime(y,m,d,H,M,S,0)
                dict_nav_data["sv"].append(sv)
                dict_nav_data["time"].append(d)
                il=23
                if len(self.line)>=il+19: clock_bias = getVal(self.line[il:il+19])
                else: clock_bias = float('NaN')
                dict_nav_data["SVclockBias"].append(clock_bias)
                il=il+19
                if len(self.line)>=il+19: SVclockDrift = getVal(self.line[il:il+19])
                else: SVclockDrift = float('NaN')
                dict_nav_data["SVclockDrift"].append(SVclockDrift)
                il=il+19
                if len(self.line)>=il+19: SVclockDrift2 = getVal(self.line[il:il+19])
                else: SVclockDrift2 = float('NaN')
                dict_nav_data["SVclockDrift2"].append(SVclockDrift2)

                self.line = self.file.readline()
                il=4
                if len(self.line)>=il+19: IODE = getVal(self.line[il:il+19])
                else: IODE = float('NaN')
                dict_nav_data["IODE"].append(IODE)
                il=il+19
                if len(self.line)>=il+19: Crs = getVal(self.line[il:il+19])
                else: Crs = float('NaN')
                dict_nav_data["Crs"].append(Crs)
                il=il+19                
                if len(self.line)>=il+19: DeltaN = getVal(self.line[il:il+19])
                else: DeltaN = float('NaN')
                dict_nav_data["DeltaN"].append(DeltaN)                
                il=il+19
                if len(self.line)>=il+19: M0 = getVal(self.line[il:il+19])
                else: M0 = float('NaN')
                dict_nav_data["M0"].append(M0)    

                self.line = self.file.readline()
                il=4
                if len(self.line)>=il+19: Cuc = getVal(self.line[il:il+19])
                else: Cuc = float('NaN')
                dict_nav_data["Cuc"].append(Cuc)
                il=il+19
                if len(self.line)>=il+19: Eccentricity = getVal(self.line[il:il+19])
                else: Eccentricity = float('NaN')
                dict_nav_data["Eccentricity"].append(Eccentricity)
                il=il+19
                if len(self.line)>=il+19: Cus = getVal(self.line[il:il+19])
                else: Cus = float('NaN')
                dict_nav_data["Cus"].append(Cus)                
                il=il+19
                if len(self.line)>=il+19: sqrtA = getVal(self.line[il:il+19])
                else: sqrtA = float('NaN')
                dict_nav_data["sqrtA"].append(sqrtA) 

                self.line = self.file.readline()
                il=4
                if len(self.line)>=il+19: Toe = getVal(self.line[il:il+19])
                else: Toe = float('NaN')
                dict_nav_data["Toe"].append(Toe)
                il=il+19
                if len(self.line)>=il+19: Cic = getVal(self.line[il:il+19])
                else: Cic = float('NaN')
                dict_nav_data["Cic"].append(Cic)
                il=il+19
                if len(self.line)>=il+19: Omega0 = getVal(self.line[il:il+19])
                else: Omega0 = float('NaN')
                dict_nav_data["Omega0"].append(Omega0)                
                il=il+19
                if len(self.line)>=il+19: Cis = getVal(self.line[il:il+19])
                else: Cis = float('NaN')
                dict_nav_data["Cis"].append(Cis) 

                self.line = self.file.readline()
                il=4
                if len(self.line)>=il+19: Io = getVal(self.line[il:il+19])
                else: Io = float('NaN')
                dict_nav_data["Io"].append(Io)
                il=il+19
                if len(self.line)>=il+19: Crc = getVal(self.line[il:il+19])
                else: Crc = float('NaN')
                dict_nav_data["Crc"].append(Crc)
                il=il+19
                if len(self.line)>=il+19: omega = getVal(self.line[il:il+19])
                else: omega = float('NaN')
                dict_nav_data["omega"].append(omega)                
                il=il+19
                if len(self.line)>=il+19: OmegaDot = getVal(self.line[il:il+19])
                else: OmegaDot = float('NaN')
                dict_nav_data["OmegaDot"].append(OmegaDot) 


                self.line = self.file.readline()
                il=4
                if len(self.line)>=il+19: IDOT = getVal(self.line[il:il+19])
                else: IDOT = float('NaN')
                dict_nav_data["IDOT"].append(IDOT)
                il=il+19
                if len(self.line)>=il+19: L2 = getVal(self.line[il:il+19])
                else: L2 = float('NaN')
                dict_nav_data["L2"].append(L2)
                il=il+19
                if len(self.line)>=il+19: Week = getVal(self.line[il:il+19])
                else: Week = float('NaN')
                dict_nav_data["Week"].append(Week)                
                il=il+19
                if len(self.line)>=il+19: L2Pflag = getVal(self.line[il:il+19])
                else: L2Pflag = float('NaN')
                dict_nav_data["L2Pflag"].append(L2Pflag) 
                
                self.line = self.file.readline()
                il=4
                if len(self.line)>=il+19: svAcc = getVal(self.line[il:il+19])
                else: svAcc = float('NaN')
                dict_nav_data["svAcc"].append(svAcc)
                il=il+19
                if len(self.line)>=il+19: health = getVal(self.line[il:il+19])
                else: health = float('NaN')
                dict_nav_data["health"].append(health)
                il=il+19
                if len(self.line)>=il+19: TGD = getVal(self.line[il:il+19])
                else: TGD = float('NaN')
                dict_nav_data["TGD"].append(TGD)                
                il=il+19
                if len(self.line)>=il+19: IODC = getVal(self.line[il:il+19])
                else: IODC = float('NaN')
                dict_nav_data["IODC"].append(IODC) 

                self.line = self.file.readline().replace("\n","")
                il=4
                transTime = float(self.line[il:il+19].replace("D","E"))
                if len(self.line)>=il+19: transTime = getVal(self.line[il:il+19])
                else: transTime = float('NaN')
                dict_nav_data["transTime"].append(transTime)
                il=il+19
                if len(self.line)>=il+19: BNK = getVal(self.line[il:il+19])
                else: BNK = float('NaN')
                dict_nav_data["BNK"].append(BNK)
            self.line = self.file.readline()

        self.df_nav = pd.DataFrame(dict_nav_data)
        self.df_nav.set_index(["time","sv"],inplace=True)
 

    def read_nav3_xyz(self):
        
        i = 0
        dict_nav_data = {
            "sv":[], "time":[], "SVclockBias":[], "SVrelFreqBias": [], "MessageFrameTime":[],
            "X":[], "dX":[], "dX2":[], "health":[],
            "Y":[], "dY":[], "dY2":[], "URA":[],
            "Z":[], "dZ":[], "dZ2":[], "IODN":[]           
        }
        
        while self.line:
            
            i += 1    
            if self.line[0]!=" ":

                sv=self.line[:3]
                sv = sv.replace(" ","0")
                y = int(self.line[4:8])
                m = int(self.line[9:11])
                d = int(self.line[12:14])
                H = int(self.line[15:17])
                M = int(self.line[18:20])
                S = int(self.line[21:23])
                d = datetime.datetime(y,m,d,H,M,S,0)
                dict_nav_data["sv"].append(sv)
                dict_nav_data["time"].append(d)
                il=23
                clock_bias = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["SVclockBias"].append(clock_bias)
                il=il+19
                rel_freq_bias = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["SVrelFreqBias"].append(rel_freq_bias)
                il=il+19
                transmission_time = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["MessageFrameTime"].append(transmission_time)

                self.line = self.file.readline()
                il=4
                X = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["X"].append(X*1e3)
                il=il+19
                dX = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["dX"].append(dX*1e3)
                il=il+19
                dX2 = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["dX2"].append(dX2*1e3)                
                il=il+19
                health = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["health"].append(health)    

                self.line = self.file.readline()
                il=4
                Y = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["Y"].append(Y*1e3)
                il=il+19
                dY = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["dY"].append(dY*1e3)
                il=il+19
                dY2 = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["dY2"].append(dY2*1e3)                
                il=il+19
                acc_code = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["URA"].append(acc_code) 

                self.line = self.file.readline()
                il=4
                Z = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["Z"].append(Z*1e3)
                il=il+19
                dZ = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["dZ"].append(dZ*1e3)
                il=il+19
                dZ2 = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["dZ2"].append(dZ2*1e3)                
                il=il+19
                iodn = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["IODN"].append(iodn) 

            self.line = self.file.readline()

        self.df_nav = pd.DataFrame(dict_nav_data)
        self.df_nav.set_index(["time","sv"],inplace=True)
    
            #if i>20: break

    def read_nav4_ephemerids(self):


        i = 0
        dict_nav_data = {
            "sv":[], "time":[], "SVclockBias":[], "SVclockDrift": [], "SVclockDrift2":[],
            "IODE":[], "Crs":[], "DeltaN":[], "M0":[],
            "Cuc":[], "Eccentricity":[], "Cus":[], "sqrtA":[],
            "Toe":[], "Cic":[], "Omega0":[], "Cis":[],
            "Io":[], "Crc":[], "omega":[], "OmegaDot":[],
            "IDOT":[], "L2":[], "Week":[], "L2Pflag":[],
            "svAcc":[], "health":[], "TGD":[], "IODC":[],
            "transTime":[], "BNK":[]
        }
        
        while self.line:
           
            i += 1    
            if self.line[0:5]=="> EPH":
                self.line = self.file.readline()
                sv=self.line[:3]
                sv = sv.replace(" ","0")
                y = int(self.line[4:8])
                m = int(self.line[9:11])
                d = int(self.line[12:14])
                H = int(self.line[15:17])
                M = int(self.line[18:20])
                S = int(self.line[21:23])
                d = datetime.datetime(y,m,d,H,M,S,0)
                dict_nav_data["sv"].append(sv)
                dict_nav_data["time"].append(d)
                #il=23
                il=23
                if len(self.line)>=il+19: clock_bias = getVal(self.line[il:il+19])
                else: clock_bias = float('NaN')
                dict_nav_data["SVclockBias"].append(clock_bias)
                il=il+19
                if len(self.line)>=il+19: SVclockDrift = getVal(self.line[il:il+19])
                else: SVclockDrift = float('NaN')
                dict_nav_data["SVclockDrift"].append(SVclockDrift)
                il=il+19
                if len(self.line)>=il+19: SVclockDrift2 = getVal(self.line[il:il+19])
                else: SVclockDrift2 = float('NaN')
                dict_nav_data["SVclockDrift2"].append(SVclockDrift2)

                self.line = self.file.readline()
                il=4
                if len(self.line)>=il+19: IODE = getVal(self.line[il:il+19])
                else: IODE = float('NaN')
                dict_nav_data["IODE"].append(IODE)
                il=il+19
                if len(self.line)>=il+19: Crs = getVal(self.line[il:il+19])
                else: Crs = float('NaN')
                dict_nav_data["Crs"].append(Crs)
                il=il+19                
                if len(self.line)>=il+19: DeltaN = getVal(self.line[il:il+19])
                else: DeltaN = float('NaN')
                dict_nav_data["DeltaN"].append(DeltaN)                
                il=il+19
                if len(self.line)>=il+19: M0 = getVal(self.line[il:il+19])
                else: M0 = float('NaN')
                dict_nav_data["M0"].append(M0)    

                self.line = self.file.readline()
                il=4
                if len(self.line)>=il+19: Cuc = getVal(self.line[il:il+19])
                else: Cuc = float('NaN')
                dict_nav_data["Cuc"].append(Cuc)
                il=il+19
                if len(self.line)>=il+19: Eccentricity = getVal(self.line[il:il+19])
                else: Eccentricity = float('NaN')
                dict_nav_data["Eccentricity"].append(Eccentricity)
                il=il+19
                if len(self.line)>=il+19: Cus = getVal(self.line[il:il+19])
                else: Cus = float('NaN')
                dict_nav_data["Cus"].append(Cus)                
                il=il+19
                if len(self.line)>=il+19: sqrtA = getVal(self.line[il:il+19])
                else: sqrtA = float('NaN')
                dict_nav_data["sqrtA"].append(sqrtA) 

                self.line = self.file.readline()
                il=4
                if len(self.line)>=il+19: Toe = getVal(self.line[il:il+19])
                else: Toe = float('NaN')
                dict_nav_data["Toe"].append(Toe)
                il=il+19
                if len(self.line)>=il+19: Cic = getVal(self.line[il:il+19])
                else: Cic = float('NaN')
                dict_nav_data["Cic"].append(Cic)
                il=il+19
                if len(self.line)>=il+19: Omega0 = getVal(self.line[il:il+19])
                else: Omega0 = float('NaN')
                dict_nav_data["Omega0"].append(Omega0)                
                il=il+19
                if len(self.line)>=il+19: Cis = getVal(self.line[il:il+19])
                else: Cis = float('NaN')
                dict_nav_data["Cis"].append(Cis) 

                self.line = self.file.readline()
                il=4
                if len(self.line)>=il+19: Io = getVal(self.line[il:il+19])
                else: Io = float('NaN')
                dict_nav_data["Io"].append(Io)
                il=il+19
                if len(self.line)>=il+19: Crc = getVal(self.line[il:il+19])
                else: Crc = float('NaN')
                dict_nav_data["Crc"].append(Crc)
                il=il+19
                if len(self.line)>=il+19: omega = getVal(self.line[il:il+19])
                else: omega = float('NaN')
                dict_nav_data["omega"].append(omega)                
                il=il+19
                if len(self.line)>=il+19: OmegaDot = getVal(self.line[il:il+19])
                else: OmegaDot = float('NaN')
                dict_nav_data["OmegaDot"].append(OmegaDot) 


                self.line = self.file.readline()
                il=4
                if len(self.line)>=il+19: IDOT = getVal(self.line[il:il+19])
                else: IDOT = float('NaN')
                dict_nav_data["IDOT"].append(IDOT)
                il=il+19
                if len(self.line)>=il+19: L2 = getVal(self.line[il:il+19])
                else: L2 = float('NaN')
                dict_nav_data["L2"].append(L2)
                il=il+19
                if len(self.line)>=il+19: Week = getVal(self.line[il:il+19])
                else: Week = float('NaN')
                dict_nav_data["Week"].append(Week)                
                il=il+19
                if len(self.line)>=il+19: L2Pflag = getVal(self.line[il:il+19])
                else: L2Pflag = float('NaN')
                dict_nav_data["L2Pflag"].append(L2Pflag) 
                
                self.line = self.file.readline()
                il=4
                if len(self.line)>=il+19: svAcc = getVal(self.line[il:il+19])
                else: svAcc = float('NaN')
                dict_nav_data["svAcc"].append(svAcc)
                il=il+19
                if len(self.line)>=il+19: health = getVal(self.line[il:il+19])
                else: health = float('NaN')
                dict_nav_data["health"].append(health)
                il=il+19
                if len(self.line)>=il+19: TGD = getVal(self.line[il:il+19])
                else: TGD = float('NaN')
                dict_nav_data["TGD"].append(TGD)                
                il=il+19
                if len(self.line)>=il+19: IODC = getVal(self.line[il:il+19])
                else: IODC = float('NaN')
                dict_nav_data["IODC"].append(IODC) 

                self.line = self.file.readline().replace("\n","")
                il=4
                transTime = float(self.line[il:il+19].replace("D","E"))
                if len(self.line)>=il+19: transTime = getVal(self.line[il:il+19])
                else: transTime = float('NaN')
                dict_nav_data["transTime"].append(transTime)
                il=il+19
                if len(self.line)>=il+19: BNK = getVal(self.line[il:il+19])
                else: BNK = float('NaN')
                dict_nav_data["BNK"].append(BNK)
                

                '''
                clock_bias = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["SVclockBias"].append(clock_bias)
                il=il+19
                rel_freq_bias = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["SVclockDrift"].append(rel_freq_bias)
                il=il+19
                transmission_time = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["SVclockDrift2"].append(transmission_time)

                self.line = self.file.readline()
                il=4
                IODE = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["IODE"].append(IODE)
                il=il+19
                Crs = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["Crs"].append(Crs)
                il=il+19
                DeltaN = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["DeltaN"].append(DeltaN)                
                il=il+19
                M0 = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["M0"].append(M0)    

                self.line = self.file.readline()
                il=4
                Cuc = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["Cuc"].append(Cuc)
                il=il+19
                Eccentricity = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["Eccentricity"].append(Eccentricity)
                il=il+19
                Cus = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["Cus"].append(Cus)                
                il=il+19
                sqrtA = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["sqrtA"].append(sqrtA) 

                self.line = self.file.readline()
                il=4
                Toe = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["Toe"].append(Toe)
                il=il+19
                Cic = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["Cic"].append(Cic)
                il=il+19
                Omega0 = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["Omega0"].append(Omega0)                
                il=il+19
                Cis = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["Cis"].append(Cis) 

                self.line = self.file.readline()
                il=4
                Io = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["Io"].append(Io)
                il=il+19
                Crc = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["Crc"].append(Crc)
                il=il+19
                omega = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["omega"].append(omega)                
                il=il+19
                OmegaDot = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["OmegaDot"].append(OmegaDot) 


                self.line = self.file.readline()
                il=4
                IDOT = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["IDOT"].append(IDOT)
                il=il+19
                L2 = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["L2"].append(L2)
                il=il+19
                Week = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["Week"].append(Week)                
                il=il+19
                if len(self.line)>=il+19: L2Pflag = getVal(self.line[il:il+19])
                else: L2Pflag = float('NaN')
                dict_nav_data["L2Pflag"].append(L2Pflag) 
                
                self.line = self.file.readline()
                il=4
                if len(self.line)>=il+19: svAcc = getVal(self.line[il:il+19])
                else: svAcc = float('NaN')
                dict_nav_data["svAcc"].append(svAcc)
                il=il+19
                if len(self.line)>=il+19: health = getVal(self.line[il:il+19])
                else: health = float('NaN')
                dict_nav_data["health"].append(health)
                il=il+19
                if len(self.line)>=il+19: TGD = getVal(self.line[il:il+19])
                else: TGD = float('NaN')
                dict_nav_data["TGD"].append(TGD)                
                il=il+19
                if len(self.line)>=il+19: IODC = getVal(self.line[il:il+19])
                else: IODC = float('NaN')
                dict_nav_data["IODC"].append(IODC) 

                self.line = self.file.readline()
                il=4
                if len(self.line)>=il+19: transTime = getVal(self.line[il:il+19])
                else: transTime = float('NaN')
                dict_nav_data["transTime"].append(transTime)
                il=il+19
                if len(self.line)>=il+19: BNK = getVal(self.line[il:il+19])
                else: BNK = float('NaN')
                dict_nav_data["BNK"].append(BNK)

                dict_nav_data
                '''
            self.line = self.file.readline()

        self.df_nav = pd.DataFrame(dict_nav_data)
        self.df_nav.set_index(["time","sv"],inplace=True)


    def read_nav4_xyz(self):
        
        i = 0
        dict_nav_data = {
            "sv":[], "time":[], "SVclockBias":[], "SVrelFreqBias": [], "MessageFrameTime":[],
            "X":[], "dX":[], "dX2":[], "health":[],
            "Y":[], "dY":[], "dY2":[], "URA":[],
            "Z":[], "dZ":[], "dZ2":[], "IODN":[]           
        }
        
        while self.line:
            i += 1    
            if self.line[0:5]=="> EPH":
                self.line = self.file.readline()
                sv=self.line[:3]
                sv = sv.replace(" ","0")
                y = int(self.line[4:8])
                m = int(self.line[9:11])
                d = int(self.line[12:14])
                H = int(self.line[15:17])
                M = int(self.line[18:20])
                S = int(self.line[21:23])
                d = datetime.datetime(y,m,d,H,M,S,0)
                dict_nav_data["sv"].append(sv)
                dict_nav_data["time"].append(d)
                il=23
                clock_bias = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["SVclockBias"].append(clock_bias)
                il=il+19
                rel_freq_bias = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["SVrelFreqBias"].append(rel_freq_bias)
                il=il+19
                transmission_time = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["MessageFrameTime"].append(transmission_time)

                self.line = self.file.readline()
                il=4
                X = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["X"].append(X*1e3)
                il=il+19
                dX = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["dX"].append(dX*1e3)
                il=il+19
                dX2 = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["dX2"].append(dX2*1e3)                
                il=il+19
                health = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["health"].append(health)    

                self.line = self.file.readline()
                il=4
                Y = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["Y"].append(Y*1e3)
                il=il+19
                dY = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["dY"].append(dY*1e3)
                il=il+19
                dY2 = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["dY2"].append(dY2*1e3)                
                il=il+19
                acc_code = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["URA"].append(acc_code) 

                self.line = self.file.readline()
                il=4
                Z = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["Z"].append(Z*1e3)
                il=il+19
                dZ = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["dZ"].append(dZ*1e3)
                il=il+19
                dZ2 = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["dZ2"].append(dZ2*1e3)                
                il=il+19
                iodn = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["IODN"].append(iodn) 

            self.line = self.file.readline()

        self.df_nav = pd.DataFrame(dict_nav_data)
        self.df_nav.set_index(["time","sv"],inplace=True)


