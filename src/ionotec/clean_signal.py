import pandas as pd
import numpy as np
import sys

from . import reconstruction as reco


def is_in_interval(t,interval):
    if (t<interval[1]) and (t>=interval[0]): return True
    else: return False
        
def filter_outsider_tracks(df,intervals_not_filtered):

    threshold = 15
    
    df_out = df.copy()
    mask_intervals = pd.Series(False, index=df_out.index)
    for interval in intervals_not_filtered:
        mask_intervals = mask_intervals | ((df_out.index>=interval[0]) & (df_out.index<interval[1]))
    df_not_filtered = df_out[mask_intervals]
    df_out = df_out[~mask_intervals]
    while True:

        df_stats = df_out.groupby(df_out.index)["VTEC"].agg(
            VTEC_mean="mean",
            #VTEC_std_lower=lambda x: x[x < x.mean()].std(),
            #VTEC_std_upper=lambda x: x[x > x.mean()].std(),
        )
        
        #dict_metrics = {'sv':[],'C1':[],'C2':[],'ti':[],'tf':[],'av_diff':[],'av_std':[]}
        dict_metrics = {'sv':[],'C1':[],'C2':[],'ti':[],'tf':[],'av_diff':[]}
        sv_channel = list(df_out[['sv', 'C1', 'C2']].drop_duplicates().itertuples(index=False, name=None))
        for (sv,C1,C2) in sv_channel:
            
            df_sv = df_out[ (df_out['sv']==sv) & (df_out['C1']==C1) & (df_out['C2']==C2) ]

            border_vtec = reco.list_leaps_series(df_sv['VTEC'],tol_dev=10.,N=60,resolution=60)
            mask_border = pd.Series(False, index=df_sv.index)
            for border in border_vtec:
                mask_border = mask_border | ((df_sv.index>=border[0]) & (df_sv.index<border[1]))
            df_sv = df_sv[mask_border]
            
            list_t,list_sigma_low,list_sigma_up,list_sigma_diff = [],[],[],[]
        
            av_diff_VTEC = []
            std_low_diff_VTEC = []
            std_up_diff_VTEC = []
            
            for iborder,border in enumerate(border_vtec):

                df_segment = df_sv.loc[border[0]:border[1]]
                df_stat_segment = df_stats.loc[border[0]:border[1]]
        
                list_t.append(df_stat_segment.index[0]+(df_stat_segment.index[-1]-df_stat_segment.index[0])/2)

                df_segment['VTEC_mean'] = df_stat_segment['VTEC_mean']
                #df_segment['VTEC_std_upper'] = df_stat_segment['VTEC_std_upper']
                #df_segment['VTEC_std_lower'] = df_stat_segment['VTEC_std_lower']
                df_segment['diff_VTEC'] = df_segment['VTEC']-df_stat_segment['VTEC_mean']
                #df_segment['distance_sigma'] = np.nan
                mask_over_mean = df_segment['diff_VTEC']>0
                #df_segment.loc[mask_over_mean,'distance_sigma'] = df_segment.loc[mask_over_mean,'diff_VTEC'] - df_segment.loc[mask_over_mean,'VTEC_std_upper']
                #df_segment.loc[~mask_over_mean,'distance_sigma'] = df_segment.loc[~mask_over_mean,'VTEC_std_lower'] - df_segment.loc[~mask_over_mean,'diff_VTEC']
                av_diff = abs(df_segment['diff_VTEC'].mean())
                #av_std = abs(df_segment['distance_sigma'].mean())
                if (av_diff<threshold): continue 
                dict_metrics['sv'].append(sv)
                dict_metrics['C1'].append(C1)
                dict_metrics['C2'].append(C2)
                dict_metrics['ti'].append(border[0])
                dict_metrics['tf'].append(border[1])
                dict_metrics['av_diff'].append(av_diff)
                #dict_metrics['av_std'].append(df_segment['distance_sigma'].mean())

        
        
        df_metrics = pd.DataFrame(dict_metrics)
        if len(df_metrics)==0: break
        df_metrics.sort_values('av_diff',inplace=True,ascending=False)
        
        list_t_intervals_affected = []
        for irow,row in df_metrics.iterrows():
            if row['av_diff']<threshold: break
            interval_already_affected = False
            for interval in list_t_intervals_affected:
                if (is_in_interval(row['ti'],interval) or is_in_interval(row['tf'],interval)):
                    interval_already_affected = True
                    break
            if interval_already_affected: continue
            mask = (df_out['sv']==row['sv']) & (df_out['C1']==row['C1']) & (df_out['C2']==row['C2'])
            df_out = df_out[~mask]
            list_t_intervals_affected.append([row['ti'],row['tf']])
        break
            
    return pd.concat([df_not_filtered,df_out])
