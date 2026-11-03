import pandas as pd
import numpy as np
import sys
import warnings
from scipy.spatial import cKDTree

from . import reconstruction as reco


def is_in_interval(t,interval):
    if (t<interval[1]) and (t>=interval[0]): return True
    else: return False
        
def filter_outsider_tracks(df,intervals_not_filtered=[]):

    threshold = 15
    
    df_out = df
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

def get_dcb_outsider_periods(df, threshold_std=2):

    print (df)
    # Convert time columns to datetime objects if not already
    df['time_i'] = pd.to_datetime(df['time_i'])
    df['time_f'] = pd.to_datetime(df['time_f'])

    print ('clean dcb 1')

    d_outsiders = {'station':[],'const':[],'time_i':[],'time_f':[]}

    for station in df['station'].unique():


      print (station)

      # Filter for the given station
      station_df = df[(df['station'] == station)].copy()

      if station_df.empty:
          print(f"No data found for station={station}")
          continue

      # Create a unique identifier for each combination (constellation, C1, C2)
      station_df['combo_id'] = station_df['constellation'] + '_' + station_df['C1'] + '_' + station_df['C2']

      

      #print(f"\n--- Leave-One-Out DCB Statistics for Station: {station} ---")
      for combo_id_val in station_df['combo_id'].unique():
          df_combo = station_df[station_df['combo_id'] == combo_id_val].copy()
          if len(df_combo) <= 2:
              #print(f"  {combo_id_val}: Not enough data points to compute leave-one-out statistics.")
              continue

          #print(f"\n  Combo ID: {combo_id_val}")
          df_combo['DCB'] = df_combo['DCB'].astype(float)
          std_tot = df_combo['DCB'].std()
          # Not enough data to filter or std_tot is low enough to trust
          if math.isnan(std_tot) or std_tot<threshold_std: continue 

          list_idx = []
          list_std = []
          last_std = 0
          # Looking for outsiders
          for idx, row in df_combo.iterrows():
              removed_dcb_value = row['DCB']
              # Create a temporary DataFrame excluding the current row
              temp_df_combo = df_combo.drop(index=idx)

              temp_mean = temp_df_combo['DCB'].mean()
              temp_std = temp_df_combo['DCB'].std()
              list_std.append(temp_std)
              list_idx.append(idx)
              #if math.isnan(temp_std) or temp_std<threshold_std:
              #  print(f"  {station},{combo_id_val}:  Removed DCB: {removed_dcb_value:.4f} -> Mean (w/o): {temp_mean:.4f}, Std (w/o): {temp_std:.4f}, Last Std (w/o): {last_std:.4f}")
              last_std = temp_std
          
          for i in range(len(list_std)):
            av_std = 0
            for j in range(len(list_std)):
              if i!=j: av_std+=list_std[j]
            av_std/=len(list_std)-1
            if 1.5*list_std[i]<av_std:

              time_i = df_combo.loc[list_idx[i]]['time_i']
              time_f = df_combo.loc[list_idx[i]]['time_f']
              const = df_combo.loc[list_idx[i]]['constellation']
              print (f"  {station},{combo_id_val}: {list_idx[i]} {list_std[i]} {av_std} {const} {time_i} {time_f}")
              d_outsiders['station'].append(station)
              d_outsiders['const'].append(const)
              d_outsiders['time_i'].append(time_i)
              d_outsiders['time_f'].append(time_f)

    df = pd.DataFrame(d_outsiders)
    return df




def _flag_epoch(lat, lon, v, k, radius_deg, k_min, nsigma, min_resid, ridge, n_iter):
    """Leave-one-out robust local-plane test for all points of one epoch.
    Returns a boolean array: True = outlier. Points with too few neighbors
    are never flagged (undecided)."""
    n = len(v)
    flags = np.zeros(n, bool)
    if n <= k_min:
        return flags

    # --- neighbors: kNN on the unit sphere, limited to radius_deg (great-circle)
    la, lo = np.deg2rad(lat), np.deg2rad(lon)
    xyz = np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    chord = 2 * np.sin(np.deg2rad(radius_deg) / 2)
    d, idx = cKDTree(xyz).query(xyz, k=min(k + 1, n), distance_upper_bound=chord)
    valid = (idx < n) & (idx != np.arange(n)[:, None])      # drops missing + the point itself
    idx = np.where(valid, idx, 0)
    nv = valid.sum(1)

    # --- local tangent coordinates (scaled to [-1, 1]) and distance weights
    dy = (lat[idx] - lat[:, None]) / radius_deg
    dx = ((lon[idx] - lon[:, None] + 180) % 360 - 180) * np.cos(la)[:, None] / radius_deg
    dv = v[idx]
    r = np.where(valid, d, chord) / chord
    w_dist = valid * np.clip(1 - r**3, 0, 1) ** 3            # tricube
    A = np.stack([np.ones_like(dx), dx, dy], axis=-1)        # [N, k, 3]
    eye = 1e-9 * np.eye(3)

    w = w_dist
    for it in range(n_iter + 1):
        # weighted ridge plane fit WITHOUT the point; intercept = prediction at the point
        Aw = A * w[..., None]
        M = np.einsum("nki,nkj->nij", Aw, A) + eye
        sw = w.sum(1)
        M[:, 1, 1] += ridge * sw                              # shrink poorly constrained slopes
        M[:, 2, 2] += ridge * sw
        b = np.einsum("nki,nk->ni", Aw, dv)
        coef = np.linalg.solve(M, b[..., None])[..., 0]
        res = np.where(valid, dv - np.einsum("nkj,nj->nk", A, coef), np.nan)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)   # all-NaN rows (no neighbors)
            s = 1.4826 * np.nanmedian(np.abs(res), axis=1)    # robust noise of the neighborhood
        s = s * np.sqrt(nv / np.maximum(nv - 3, 1))           # dof correction for 3 fitted params
        s_eff = np.maximum(s, min_resid / nsigma)             # absolute floor in TECU

        if it < n_iter:                                       # Tukey biweight reweighting
            u = res / (4.685 * s_eff[:, None])
            w = w_dist * np.where(np.abs(u) < 1, (1 - u**2) ** 2, 0.0)

    decided = (nv >= k_min) & ((w > 0).sum(1) >= k_min)
    return decided & (np.abs(v - coef[:, 0]) > nsigma * s_eff)


def filter_vtec_outliers(df, station_col="station", sat_col="sv",
                         k=20, radius_deg=5.0, k_min=6,
                         nsigma=4.0, min_resid=3.0, ridge=0.02, n_iter=2,
                         return_info=False):
    """
    Spatial outlier filter for VTEC (columns lat, lon, vtec, station, satellite;
    DatetimeIndex with many points per epoch). For each epoch and each point:
      1. take up to k nearest neighbors within radius_deg (great-circle), excluding the point;
      2. fit a robust (distance- and biweight-weighted) local plane to the neighbors;
      3. flag the point if |vtec - plane(point)| > max(nsigma * robust_sigma, min_resid).
    If a point is flagged, ALL rows of its (station, satellite) pair are removed from the
    whole DataFrame, whatever the epoch or pass.

    radius_deg : neighborhood radius [deg]       k_min : min neighbors to decide
    nsigma     : threshold in robust sigmas      min_resid : absolute floor [TECU]
    ridge      : slope shrinkage (collinear neighbors -> falls back to local mean)
    n_iter     : biweight reweighting iterations
    Returns the cleaned DataFrame (and, with return_info=True, a dict {'flagged' (points), 'removed' (rows)} of boolean arrays).
    """
    lat = df["lat"].to_numpy(float)
    lon = df["lon"].to_numpy(float)
    v = df["vtec"].to_numpy(float)
    flags = np.zeros(len(df), bool)

    ok = np.flatnonzero(np.isfinite(lat) & np.isfinite(lon) & np.isfinite(v))
    codes, _ = pd.factorize(df.index[ok])
    order = np.argsort(codes, kind="stable")
    for g in np.split(order, np.flatnonzero(np.diff(codes[order])) + 1):
        pos = ok[g]
        flags[pos] = _flag_epoch(lat[pos], lon[pos], v[pos], k, radius_deg, k_min,
                                 nsigma, min_resid, ridge, n_iter)

    # every row of a (station, satellite) pair that has at least one flagged point
    removed = (pd.Series(flags).groupby([df[station_col].to_numpy(), df[sat_col].to_numpy()])
               .transform("any").to_numpy())

    clean = df[~removed]
    return (clean, {"flagged": flags, "removed": removed}) if return_info else clean
