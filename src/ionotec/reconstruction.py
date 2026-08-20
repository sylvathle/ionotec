import numpy as np
import pandas as pd

def fit_lin(t,sig):
    N=len(t)
    if N<2: return None, None
    # Coefficients of paraboloid of error function (N*mse)
    a,b,c,d,e=0,N,0,0,0 # b=N for B^2 coef
    # Iterate over subseries to calculate coefficients
    for i in range(N):
        a+=t[i]**2 # A^2 coef
        c+=2*t[i] # A*B coef
        d-=2*t[i]*sig[i] # A coef
        e-=2*sig[i] # B coef

    if (a==0) and (c==0):
        return float('NaN'),float('NaN'),float('NaN'),float('NaN')
        
    # Forward A and B parameters of linear fit (solve the minimum of mse)
    A=-(2*b*d-c*e)/(4*a*b-c**2)
    B=-(c*d-2*a*e)/(c**2-4*a*b)


    return A,B
    
    
    
    
    
def list_leaps_series(series,tol_dev=0.2,N=5,resolution = 60):
 
    # Discard too short series
    if len(series)<=N: return []   

    indices = series.index
    list_series = series.tolist()
    
    # List time in seconds from first index of series
    diffs = [0]
    delta_t = [0]
    for i in range(len(list_series)-1):
        diffs.append(diffs[-1]+(indices[i+1]-indices[i]).seconds)
        delta_t.append((indices[i+1]-indices[i]).seconds)
    delta_t.append(0) # add a fake delta_t to be able to finish the future iterations
    tol_deviation=tol_dev*resolution
    
    border_left, border_right = None, None

    list_borders = []

    A_mem_right = []
    dev_mem_right = []
    
    s=0
    while s<len(series):
        # Get fit of segment s:s+N, this allows to look if s-1 or s+N value are borders of segments
        ### For the cases where s+N>len(series) we are in principle outside the diffs list, 
        ### the fit_lin responds None in case the fit could not be performed, 
        ### this will affect the left border
        A,B = fit_lin(diffs[s:s+N],list_series[s:s+N])

        # Deviation of left point before the fit, 0 if it's the first (left border of series)
        left_dev=float(abs(list_series[s-1]-A*diffs[s-1]-B)) if (s>0) and (A!=None) and (B!=None) else 0
        
        # Deviation of right point after the fit, 0 if it's the first (left border of series)
        right_dev=float(abs(list_series[s+N]-A*diffs[s+N]-B)) if s+N<len(list_series) else 0
        A_mem_right.append(A)
        dev_mem_right.append(right_dev)

        condition_border = False
        
        if s<N: 
            if (left_dev>tol_deviation) or (delta_t[s]>resolution*N) or (s==0):
                border_left, border_right, A_left = s, None, A
                A_left = A
        elif s<len(series)-1:
            if ((left_dev>tol_deviation) and ((dev_mem_right[s-N]>tol_deviation) or (dev_mem_right[s-N-1]>tol_deviation) or (dev_mem_right[s-N+1]>tol_deviation))) or (delta_t[s]>resolution*N):
                t_left_right = (indices[s-1]-indices[border_left]).seconds
                n_left_right = s-1 - border_left
    
                if (t_left_right>N*resolution) and (n_left_right>N): # Big enough to be considered as a segment
                    border_right = s-1
                    list_borders.append([indices[border_left],indices[border_right],A_left,A_mem_right[s-N]])
                    border_left, border_right = s, None
                    A_left = A
                else: # if too short, we start a new segment
                    border_left, border_right = s, None
                    A_left = A
        else:

            if (delta_t[s]<resolution*N) and (dev_mem_right[s-N]<tol_deviation):
                t_left_right = (indices[s]-indices[border_left]).seconds
                n_left_right = s - border_left
    
                if (t_left_right>N*resolution) and (n_left_right>N): # Big enough to be considered as a segment
                    border_right = s
                    list_borders.append([indices[border_left],indices[border_right],A_left,A_mem_right[s-N+1]])
            else:
                t_left_right = (indices[s-1]-indices[border_left]).seconds
                n_left_right = s-1 - border_left
    
                if (t_left_right>N*resolution) and (n_left_right>N): # Big enough to be considered as a segment
                    border_right = s-1
                    list_borders.append([indices[border_left],indices[border_right],A_left,A_mem_right[s-N]])

        s+=1
    return list_borders
    
    
    
### Adds a column to df_arc containing a polynomial fit of STEC_p, for each segment of borders_stecp
#### Outside the segments or if the segment cannot satisfactorily adjusted, it is filled with NaN. 
def fit_STECP(df_arc,borders_stecp):



    df_arc["STEC_p_fit"] = np.nan

    for borders in borders_stecp:

        filter_border = (df_arc.index>=borders[0]) & (df_arc.index<=borders[1]) 
        #### Make a fit of STEC_p to detect deviations of STEC_l
        x = (df_arc.index - df_arc.index[0]).total_seconds()
        y = df_arc["STEC_p"].to_numpy()
        w = df_arc["sin2_ele"].to_numpy()
    
        # 1. Masque booléen pour exclure les NaN de y (et s'assurer que x et w sont valides)
        valid_mask = ~np.isnan(y) & ~np.isnan(x) & ~np.isnan(w)
        valid_mask = valid_mask & filter_border
    
        # 2. Filtrage des données pour le fit
        x_fit = x[valid_mask]
        y_fit = y[valid_mask]
        w_fit = w[valid_mask]
    
        # 3. Ajustement dynamique du degré basé sur le nombre de points VALIDES
        degree_fit = min(len(x_fit) - 3, 10)
        #print (degree_fit,len(x_fit))
    
        # 4. Calcul du fit si on a assez de points valides
        if degree_fit >= 0:
            try: 

                coeffs, residuals, rank, singular_values, rcond = np.polyfit(
                    x_fit, y_fit, deg=degree_fit, full=True
                )
                
                # Calculate the condition number
                condition_number = singular_values[0] / singular_values[-1]
                
                # A very high condition number (e.g., > 1e12) means a poor fit
                if condition_number > 1e12 or rank < (degree_fit + 1):
                    #print("Fit is poorly conditioned based on singular values.")
                    continue
            except: 
                print ('STEC_p fit did not work\n')
                continue
        else:
            #print ('Not enough valid data points for fit for degree '+str(degree_fit)+' '+str(len(x))+' '+str(len(x_fit))+'\n')
            continue
    
        # 5. Évaluation sur l'axe 'x' complet (conserve les dimensions d'origine du DataFrame)
        poly = np.poly1d(coeffs) 
        df_arc.loc[valid_mask, "STEC_p_fit"] = poly(x_fit)                          

    return df_arc
    


## Correct STEC_l, based on the values of STEC_p, STEC_p_fit and the strengh of the signal.
#def correct_signal(df_arc,borders_stecp,borders_stecl):
def correct_signal(df_arc):

    df_arc.dropna(subset=['STEC_l'],inplace=True)

    borders_stecl = list_leaps_series(df_arc['STEC_l'],tol_dev=1./60.,N=3) # tol_dev=1TECu/minute
    borders_stecp = list_leaps_series(df_arc['STEC_p'],tol_dev=25./60.,N=8) # tol_dev=25TECu/minute

    N_min_stec_p = 10
    
    healthy_segment = [True for border in borders_stecp]
    last_s1_healthy = 0
    last_s2_healthy = 0
    for iborder in range(1,len(borders_stecp)):
        
        diff_s1 = df_arc.loc[borders_stecp[last_s1_healthy][1]]['S1'] - df_arc.loc[borders_stecp[iborder][0]]['S1']
        diff_s2 = df_arc.loc[borders_stecp[last_s2_healthy][1]]['S2'] - df_arc.loc[borders_stecp[iborder][0]]['S2']
        if diff_s1>10: healthy_segment[iborder] = False
        else: last_s1_healthy=iborder
        if diff_s2>10: healthy_segment[iborder] = False
        else: last_s2_healthy=iborder

        if not healthy_segment[iborder]:
            ileft = iborder-1
            iright = iborder+1 if iborder<len(borders_stecp)-1 else iborder
            #print (borders_stecp[ileft][1],borders_stecp[iright][0])
            df_arc.loc[borders_stecp[ileft][1]:borders_stecp[iright][0], 'STEC_p'] = np.nan
        else:
            df_arc.loc[borders_stecp[iborder-1][1]:borders_stecp[iborder][0], 'STEC_p'] = np.nan

    df_arc = fit_STECP(df_arc,borders_stecp)

    arc_anchored = False


    ### Remove data between segments
    # if no border is found, indicate nothing healthy here, remove and continue
    if len(borders_stecl)==0:
        df_arc["STEC_l"] = np.nan
        return df_arc

    filter_erractic_segments = (df_arc.index<borders_stecl[0][0])
    df_arc[filter_erractic_segments] = np.nan    
    for iborder in range(1,len(borders_stecl)):
    #    #print (iborder,borders_stecl[iborder-1][1],borders_stecl[iborder][0])
        filter_erractic_segments = (df_arc.index>borders_stecl[iborder-1][1]) & (df_arc.index<borders_stecl[iborder][0])
        df_arc[filter_erractic_segments] = np.nan
    filter_erractic_segments = (df_arc.index>borders_stecl[-1][1])
    df_arc[filter_erractic_segments] = np.nan    

    iborder = 0
    while iborder<len(borders_stecl):
        
        filter_segment = (df_arc.index>=borders_stecl[iborder][0]) & (df_arc.index<=borders_stecl[iborder][1])
        df_seg = df_arc[filter_segment]#/df_arc.loc[borders_stecl[iborder][0]:borders_stecl[iborder][1]]
        brs = None
        N_valid_STEC_p = len(df_seg['STEC_p'].dropna())

        if N_valid_STEC_p>N_min_stec_p:
            df_seg_br = df_seg.dropna(subset=['STEC_p'])
            df_seg_br['BRs'] = (df_seg_br["STEC_p"]-df_seg_br["STEC_l"])*df_seg_br["sin2_ele"]
            brs=df_seg_br['BRs'].sum(skipna=True)/df_seg_br["sin2_ele"].sum(skipna=True)
            df_arc.loc[filter_segment,'STEC_l']+=brs

            ## Check case where the STEC_l is not completely out of the trend of STEC_l
            #  We use a comparison with a polynomial fit of STEC_p
            #  If the difference between the fit and the STEC_p is way below the difference between STEC_l and the fit,
            #  we consider STEC_l is corrupted and we replace it by the fit of STEC_p
            df_seg = df_arc[filter_segment]
            sigma = (df_seg["STEC_p"] - df_seg["STEC_p_fit"]).std()
            df_seg['diffSTEC'] = df_seg['STEC_l']-df_seg['STEC_p_fit']
            
            ## In case STEC_l corrected with is going too far away from the true value we the part that was adjusted with STEC_p.
            if max(df_seg['diffSTEC'])>15*sigma: 
                print ('correction with fit',df_seg['sv'].tolist()[0],df_seg.index[0],df_seg.index[-1],max(df_seg['diffSTEC']),15*sigma)
                df_arc.loc[filter_segment,'STEC_l'] = df_arc.loc[filter_segment,'STEC_p_fit']

            arc_anchored = True
            iborder+=1
            continue
            
        # If no baseline could be established and the segment under analysis is not the first one, 
        #        we glue it to the previous one
        elif len(borders_stecl)>1:    
            if iborder==0:
                t_dist = (borders_stecl[1][0]-borders_stecl[0][1]).seconds
                slope = (borders_stecl[0][3]+borders_stecl[1][2])/2.0
                glue = df_arc["STEC_l"].loc[borders_stecl[1][0]] - df_arc["STEC_l"].loc[borders_stecl[0][1]] - slope*t_dist
                df_arc.loc[filter_segment,'STEC_l']+=glue
                borders_stecl[0][1]=borders_stecl[1][1]
                borders_stecl[0][3]=borders_stecl[1][3]
                del (borders_stecl[1])
                iborder-=1
            elif (iborder>0) and (iborder<len(borders_stecl)-1): 
                #We first look for the closest segment 
                t_dist_left=(borders_stecl[iborder][0]-borders_stecl[iborder-1][1]).seconds
                t_dist_right=(borders_stecl[iborder+1][0]-borders_stecl[iborder][1]).seconds
    
                if t_dist_left<=t_dist_right:
                    slope = (borders_stecl[iborder][2]+borders_stecl[iborder-1][3])/2
                    glue = df_arc["STEC_l"].loc[borders_stecl[iborder][0]] - df_arc["STEC_l"].loc[borders_stecl[iborder-1][1]] - slope*t_dist_left
                    df_arc.loc[filter_segment,'STEC_l']-=glue  
                else:
                    slope = (borders_stecl[iborder+1][2]+borders_stecl[iborder][3])/2
                    glue = df_arc["STEC_l"].loc[borders_stecl[iborder+1][0]] - df_arc["STEC_l"].loc[borders_stecl[iborder][1]] - slope*t_dist_right
                    df_arc.loc[filter_segment,'STEC_l']+=glue
                    borders_stecl[iborder][1]=borders_stecl[iborder+1][1]
                    borders_stecl[iborder][3]=borders_stecl[iborder+1][3]
                    del (borders_stecl[iborder+1])
                    iborder-=1
            else:
                t_dist = (borders_stecl[-1][0]-borders_stecl[-2][1]).seconds
                slope = (borders_stecl[-1][2]+borders_stecl[-2][3])/2.0
                glue = df_arc["STEC_l"].loc[borders_stecl[-1][0]] - df_arc["STEC_l"].loc[borders_stecl[-2][1]] - slope*t_dist
                df_arc.loc[filter_segment,'STEC_l']-=glue
                borders_stecl[0][1]=borders_stecl[1][1]
                borders_stecl[0][3]=borders_stecl[1][3]
                del (borders_stecl[1])
                iborder-=1
        iborder+=1

    if not arc_anchored:
        df_arc['STEC_l']=np.nan
                     
    return df_arc
    
    
    
    
def process_receiver_dcb(df,elevation_filter = 30):
    # Coefficients of the quadratic error function that will be computed
    a,b=0,0

    # Compute cos\chi for full serie
    #df["ci"] = np.cos(np.arcsin(R_E*np.cos(df["elevation"])/(R_E+h)))
    df.dropna(subset=["STEC_l","elevation"],inplace=True)
    df = df[df['elevation']>elevation_filter/180*np.pi]

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
            ci = row["cos_chi"]
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
        if np.isnan(a) or np.isnan(b): return float('NaN')

    if a==0: return float("nan")
    # root of error function = receiver bias.
    br = b/a
    return br

