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
    #print (d)
    
    nsat = int(spline[8])

    ### Get Satellites list
    istr = 0
    str_sats=spline[9]
    list_sats = []
    while istr<len(spline[9]):
        list_sats.append(str_sats[istr:istr+3])
        istr=istr+3

    #print ()
        
    return d,nsat,list_sats

def process_head_data_line_rinexv2(line):

    spline = line.split()
    year = int("20"+line[1:3])
    month = int(line[4:6])

    day = int(line[7:9])
    #print (spline)
    #print (day)
    #hour = int(spline[4])
    hour = int(line[10:12])
    #print (hour)
    #minute = int(spline[5])
    minute = int(line[13:15])
    #print (minute)
    
    #second = int(spline[6].split(".")[0])
    second = int(line[16:18])
    #print (second)
    
    #millisecond = int(spline[6][3:6])
    millisecond = int(line[19:22])
    #print (millisecond)
    d = datetime.datetime(year,month,day,hour,minute,second,millisecond)
    #print (d)

    nsat = int(line[30:32])
    #nsat = int(spline[8])

    ### Get Satellites list
    istr = 0
    #str_sats=spline[9]
    str_sats = line[32:]
    #print (str_sats)
    list_sats = []
    while istr<len(str_sats):
        list_sats.append(str_sats[istr:istr+3])
        istr=istr+3

    #print ()
        
    return d,nsat,list_sats

#last_item = {}
#constellations_df = {}


def update_hatanaka(list_item,sv,iobs,value):

    i = list_item[sv][iobs][-1]
    #if sv=="C01":
    #    if iobs==1:
    #        print (sv,iobs,i,list_item[sv][iobs])
    if i < 3: list_item[sv][iobs][-1] += 1
    list_item[sv][iobs][i] = value
    
    while i>0:
        i -= 1
        list_item[sv][iobs][i] += list_item[sv][iobs][i+1]


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

        #print (self.header)
        if int(self.header["version"])==2:
            if "COMPACT" in self.header['type_rinex']:
                #print ("not using georinex")
                #self.read_obs_rinexv2()
                for constellation in ['G','R']:
                    obs = gr.load(self.rinex_path, meas=self.vars_list, use=constellation)
                    df_temp = obs.to_dataframe()
                    df_temp.reset_index(inplace=True)
                    self.list_df[constellation] = df_temp                    
                return self.list_df
            if self.header['type']=='O':     
                #for constellation in ['G','R']:
                for constellation in ['G']:
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
                nav = gr.load(self.rinex_path)
                self.df_nav = nav.to_dataframe()
                return self.df_nav
                #nav = self.read_nav4()
                #return self.df_nav

        elif (int(self.header["version"])==4) and (self.header['type']=='N'):
            print ("Will proces version 4 navigation")
            nav = self.read_nav4()
            return self.df_nav
        
        else:
            print ("Warning, file",self.rinex_path)
            print ("Rinex version is: ", self.header["version"], " which is currently no supported")
            return pd.DataFrame()
            
    
    def read_obs(self):
        
        data_head_line = ""
        
        ### First datatime
        data_head_line = self.line    
        #print (data_head_line)
        d,nsat,list_sats = process_head_data_line(data_head_line)
        #print (d,nsat,list_sats)
        
        last_item = {}
        
        self.file.readline()
        isat = 0
        while isat<nsat:
            self.line = self.file.readline()
            spline = self.line.split(" ")
            #print (list_sats[isat],len(spline),spline)
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

        tt = 0
        
        while self.line:
            self.line = self.file.readline() 
            
            data_head_line = merge_str(data_head_line,self.line.replace("&","0"))
            d,nsat,list_sats = process_head_data_line(data_head_line)
        
            #print (d,data_head_line)
            
            #line = file.readline()
            if not self.file.readline(): 
                #print ("finished")
                continue
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
                                v = float(str_v)/10**shift
                                obs_item[obs_name] = v
                                last_item[sv].append([v,float('NaN'),float('NaN'),float('NaN'),1])
                            # Data is known
                            else:
                                shift = dict_obs_shift[constellation][iobs]
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
        
        self.file.readline()
    
        
        for constellation, observations in self.dict_obs.items():
            if len(observations)==0: continue
            if constellation == 'I': continue
            if not constellation in ['C']: continue
            self.list_df[constellation] = pd.DataFrame(observations)
        
    def read_obs_rinexv2(self):
            
        #data_head_line = ""
        
        ### First datatime
        data_head_line = self.line    
        #print (data_head_line)
        d,nsat,list_sats = process_head_data_line_rinexv2(data_head_line)
        #print (d,nsat,list_sats)
        #text = self.file.readlines()
        #lines = text.splitlines()
        #nlines = len(lines)

        #print (nlines)
        #print (lines)
        
         # ------------------------------------------------------------------ header
        #obs_types: list[str] = []
        #total_obs: int = 0
        #i = 0       
        #sys.exit()
        
        
        last_item = {}
        constellation = 'all'
        nline = 0
        self.line = self.file.readline()
        #print (self.line)
        nline+=1
        isat = 0
        while isat<nsat:
            self.line = self.file.readline()
            self.line = self.line.replace("\n","")
           # print (self.line)
            nline+=1
            spline = self.line.split(" ")
            #print (list_sats[isat],len(spline),spline)
            sv = list_sats[isat]
            
            last_item[sv]=[]
            obs_item = {"time":d, "sv":sv}
            iobs = 0
            #print (spline)
            n_empty = 0
            #print (d)
            #print (sv)
            #print (spline)

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
                        shift = dict_obs_shift[constellation][iobs]
                        v = float(str_v)/10**shift
                        update_hatanaka(last_item,sv,iobs,v)                
                        obs_item[obs_name] = last_item[sv][iobs][0]
                iobs = iobs+1
            isat=isat+1
            #print (obs_item)
            self.dict_obs[constellation].append(obs_item)
        #print (self.dict_obs)
        #for k,v in last_item.items():
        #    print ("----------")
        #    print (k)
        #    print (v)
            
            #for s in spline:
        '''
            for iobs, obs_name in enumerate(self.dict_obs_name[constellation]):
                print (iobs)
                print (s)
                
                if s==' ': 
                    n_empty+=1    
                    if n_empty>1:
                        continue
                    else:
                        iobs+=1
                if iobs>=len(self.dict_obs_name[constellation]): continue
              
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
            #print (obs_item)
            self.dict_obs[constellation].append(obs_item)
        '''
        #sys.exit()
        tt = 0
        
        while self.line:
            
            self.line = self.file.readline()
            #print (self.line)
            #if "COMMENT" in self.line: print ("Problem with comment")
            #print (d,nline, self.line)
            nline+=1
            data_head_line = merge_str(data_head_line,self.line.replace("&","0"))
            d,nsat,list_sats = process_head_data_line_rinexv2(data_head_line)

            #print (data_head_line)
            #print (self.line)
            
            #print (d)
            ##print (self.line[:-1])
            #print (data_head_line[:-1])
            
            self.line = self.file.readline()
            #print ("----"+self.line+"----")
            if "COMMENT" in self.line:
                while "COMMENT" in self.line:
                    self.line = self.file.readline()
                    continue
                #data_head_line = merge_str(data_head_line,self.line.replace("&","0"))
                data_head_line = self.line 
                d,nsat,list_sats = process_head_data_line_rinexv2(data_head_line)
                #print (list_sats)
                self.line = self.file.readline()
            #print (self.line)
            if not self.line:
                nline+=1
                #print ("finished")
                continue
            
            
            isat = 0
            while isat<nsat:
                self.line = self.file.readline()
                #if d>=datetime.datetime(2024,4,10,0,3,7,0) and nline<5650:
                #    print (nline,isat,self.line[:-1])
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
                                shift = dict_obs_shift[constellation][iobs]
                                v = float(str_v)/10**shift
                                update_hatanaka(last_item,sv,iobs,v)                
                                obs_item[obs_name] = last_item[sv][iobs][0]
                        iobs = iobs+1
                    isat=isat+1
                    #print (obs_item)
                    self.dict_obs[constellation].append(obs_item)
    
                else:
                    #print (d)
                    #print (isat,nsat)
                    #print (sv)
                    #print (spline)
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
                                #print ("\t",iobs)
                                #print ("\t",v)
                                #print ("\t",last_item[sv])
                                update_hatanaka(last_item,sv,iobs,v)                
                                obs_item[obs_name] = last_item[sv][iobs][0]
                        iobs = iobs+1
                    isat=isat+1
                    self.dict_obs[constellation].append(obs_item)
        
        #self.line = self.file.readline()
        #print (self.line)
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

        #del (self.dict_obs['all'])
        
        #print (pd.DataFrame(self.dict_obs['all']))
        for constellation, observations in self.dict_obs.items():
            if constellation in ['S','C','J','I']: continue
            self.list_df[constellation] = pd.DataFrame(observations)
        #for constellation, df in self.list_df.items():
            #self.list_df[constellation] = pd.DataFrame(observations)
        #    print (constellation)
        #    print (df)
        #sys.exit()
        
        #if "I" in self.list_df.keys():
        #    del self.list_df["I"]   


    '''
    def read_obs_rinexv2(self):
 
        data_head_line = self.line
        d, nsat, list_sats = process_head_data_line_rinexv2(data_head_line)
     
        last_item = {}
        constellation = 'all'
        nline = 0
        self.line = self.file.readline()
        nline += 1
        isat = 0
        while isat < nsat:
            self.line = self.file.readline()
            nline += 1
            spline = self.line.split(" ")
            sv = list_sats[isat]
     
            last_item[sv] = []
            obs_item = {"time": d, "sv": sv}
            iobs = 0
            n_empty = 0
            for s in spline[:-1]:
                if s == ' ':
                    n_empty += 1
                    if n_empty > 1:
                        continue
                    else:
                        iobs += 1
                if iobs >= len(self.dict_obs_name[constellation]): continue
     
                obs_name = self.dict_obs_name[constellation][iobs]
                if s == '':
                    v = float('NaN')
                    obs_item[obs_name] = v
                    last_item[sv].append([v, float('NaN'), float('NaN'), float('NaN'), 0])
                else:
                    ss = s.split('&')
                    shift = int(ss[0])
                    self.dict_obs_shift[constellation][iobs] = shift
                    v = float(ss[-1]) / 10 ** shift
                    obs_item[obs_name] = v
                    last_item[sv].append([v, float('NaN'), float('NaN'), float('NaN'), 1])
                iobs = iobs + 1
            isat = isat + 1
            self.dict_obs[constellation].append(obs_item)
     
        nline = 0
     
        while self.line:
     
            self.line = self.file.readline()
            if not self.line:
                break
            nline += 1
     
            # ----------------------------------------------------------------
            # Build the merged epoch header line and check for special events.
            # RINEX 2.11 / CRINEX uses epoch flag != 0 to embed header records
            # (COMMENTs, antenna moves, etc.) inline in the data body.
            # The epoch flag lives at character position 28; the record count
            # (n_special) spans positions 29-31 of the epoch line.
            # When flag != 0 we must skip n_special lines and then re-read the
            # real epoch line that follows.
            # ----------------------------------------------------------------
            data_head_line = merge_str(data_head_line, self.line.replace("&", "0"))
     
            # Extract epoch flag and n_special from the merged line
            line_padded = data_head_line.ljust(35)
            print (line_padded)
            try:
                epoch_flag = int(line_padded[28:30].strip() or 0)
                n_special  = int(line_padded[29:32].strip() or 0)
            except ValueError:
                epoch_flag = 0
                n_special  = 0
     
            # If this is a special-event epoch (flag 1-6), skip the embedded
            # header/comment records and then read the real data epoch header.
            if epoch_flag != 0:
                for _ in range(n_special):
                    self.file.readline()   # discard special record (COMMENT etc.)
                # The next line is the real epoch header with flag=0
                self.line = self.file.readline()
                if not self.line:
                    break
                print ("--------------",n_special)
                print (self.line)
                data_head_line = merge_str(data_head_line, self.line.replace("&", "0"))
                nline += 1

            print (self.line)
            print (data_head_line)
            d, nsat, list_sats = process_head_data_line_rinexv2(data_head_line)
     
            # Read the blank separator line that follows every epoch header
            self.line = self.file.readline()
            if not self.line:
                break
            nline += 1
     
            isat = 0
            while isat < nsat:
                self.line = self.file.readline()
                nline += 1
                spline = self.line[:-1].split(" ")
                sv = list_sats[isat]
                constellation = 'all'
                obs_item = {"time": d, "sv": sv}
                iobs = 0
     
                if sv not in last_item.keys():
                    last_item[sv] = []
     
                    while iobs < len(self.dict_obs_name[constellation]):
                        obs_name = self.dict_obs_name[constellation][iobs]
                        if iobs < len(spline):
                            s = spline[iobs]
                        else:
                            obs_item[obs_name] = float('NaN')
                            iobs += 1
                            last_item[sv].append([float('NaN'), float('NaN'), float('NaN'), float('NaN'), 0])
                            continue
                        if s == '':
                            v = float('NaN')
                            last_item[sv].append([v, float('NaN'), float('NaN'), float('NaN'), 0])
                            obs_item[obs_name] = v
                        else:
                            ss = s.split('&')
                            str_v = ss[-1]
                            if len(ss) == 2:
                                shift = int(ss[0])
                                v = float(str_v) / 10 ** shift
                                obs_item[obs_name] = v
                                last_item[sv].append([v, float('NaN'), float('NaN'), float('NaN'), 1])
                            else:
                                shift = self.dict_obs_shift[constellation][iobs]
                                v = float(str_v) / 10 ** shift
                                update_hatanaka(last_item, sv, iobs, v)
                                obs_item[obs_name] = last_item[sv][iobs][0]
                        iobs += 1
                    isat += 1
                    self.dict_obs[constellation].append(obs_item)
     
                else:
                    while iobs < len(self.dict_obs_name[constellation]):
                        obs_name = self.dict_obs_name[constellation][iobs]
                        if iobs < len(spline):
                            s = spline[iobs]
                        else:
                            obs_item[obs_name] = float('NaN')
                            last_item[sv][iobs] = [float('NaN'), float('NaN'), float('NaN'), float('NaN'), 0]
                            iobs += 1
                            continue
                        if s == '':
                            v = float('NaN')
                            last_item[sv][iobs] = [v, float('NaN'), float('NaN'), float('NaN'), 0]
                            obs_item[obs_name] = v
                        else:
                            ss = s.split('&')
                            str_v = ss[-1]
                            if len(ss) == 2:
                                shift = int(ss[0])
                                v = float(str_v) / 10 ** shift
                                obs_item[obs_name] = v
                                last_item[sv][iobs] = [v, float('NaN'), float('NaN'), float('NaN'), 1]
                            else:
                                shift = self.dict_obs_shift[constellation][iobs]
                                v = float(str_v) / 10 ** shift
                                update_hatanaka(last_item, sv, iobs, v)
                                obs_item[obs_name] = last_item[sv][iobs][0]
                        iobs += 1
                    isat += 1
                    self.dict_obs[constellation].append(obs_item)
     
        nline += 1
     
        constellations = set()
        for obs_item in self.dict_obs['all']:
            sv_str = obs_item.get('sv')
            constellation = sv_str[0]
            if constellation not in self.dict_obs:
                self.dict_obs[constellation] = []
            self.dict_obs[constellation].append(obs_item)
     
        del self.dict_obs['all']
     
        for constellation, observations in self.dict_obs.items():
            if constellation in ['E', 'S', 'C', 'J', 'I']: continue
            self.list_df[constellation] = pd.DataFrame(observations)

    ## Only reads Glonass or SBAS
    def read_nav(self):
        i = 0
        dict_nav_data = {
            "sv":[], "time":[], "SVclockBias":[], "SVrelFreqBias": [], "MessageFrameTime":[],
            "X":[], "dX":[], "dX2":[], "health":[],
            "Y":[], "dY":[], "dY2":[], "URA":[],
            "Z":[], "dZ":[], "dZ2":[], "IODN":[]           
        }
        
        while self.line:
            
            i += 1
            #print (self.line)
          
            if self.line[:3]!='   ':
                spline = self.line[:22].split()
                sv = spline[0]
                if len(sv)==2: sv = self.header['constellation']+sv
                if len(sv)==1: sv = self.header['constellation']+'0'+sv
                y = int(spline[1])
                if y<100: y = y+2000
                m = int(spline[2])
                d = int(spline[3])
                H = int(spline[4])
                M = int(spline[5])
                S = int(spline[6].split('.')[0])
                d = datetime.datetime(y,m,d,H,M,S,0)
                dict_nav_data["sv"].append(sv)
                dict_nav_data["time"].append(d)
                il=22
                clock_bias = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["SVclockBias"].append(clock_bias)
                il=il+19
                rel_freq_bias = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["SVrelFreqBias"].append(rel_freq_bias)
                il=il+19
                transmission_time = float(self.line[il:il+19].replace("D","E"))
                dict_nav_data["MessageFrameTime"].append(transmission_time)

                self.line = self.file.readline()
                il=3
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
                il=3
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
                il=3
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
        #sys.exit()
    '''

            
            #if i>20: break

    def read_nav4(self):
        
        i = 0
        dict_nav_data = {
            "sv":[], "time":[], "SVclockBias":[], "SVrelFreqBias": [], "MessageFrameTime":[],
            "X":[], "dX":[], "dX2":[], "health":[],
            "Y":[], "dY":[], "dY2":[], "URA":[],
            "Z":[], "dZ":[], "dZ2":[], "IODN":[]           
        }
        
        while self.line:
            
            i += 1    
            if self.line[0]==">":
                self.line = self.file.readline()
                sv=self.line[:3]
                y = int(self.line[4:8])
                m = int(self.line[9:11])
                d = int(self.line[12:14])
                H = int(self.line[15:17])
                M = int(self.line[18:20])
                S = int(self.line[21:23])
                d = datetime.datetime(y,m,d,H,M,S,0)
                dict_nav_data["sv"].append(sv)
                dict_nav_data["time"].append(d)
                il=24
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