import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import RBFInterpolator
import sys


def plot_station(df_station,png_file_name,mozaic=False):

    list_sats = df_station["sv"].unique()

    if mozaic:
        fig, axs = plt.subplots(4,2,figsize=(14,10),sharex=True,sharey=True)

        i,j=0,0
        n_sat = 0
        n=0
        N=0
        for sat in list_sats:
            i = int(n/2)
            j = n%2
            df_sat = df_station[df_station["sv"]==sat]
            axs[i,j].plot(df_sat["VTEC"],'b.',markersize=1)
            axsx = axs[i,j].twinx()
            axsx.plot(df_sat["elevation"]*180/3.1415926535,'r.',markersize=1,alpha=0.5)
            axsx.set_ylim([0,90])
            axs[i,j].set_title(sat)
            axs[i,j].grid(True)
            axs[i,j].set_ylabel("VTEC")
            axsx.set_ylabel("elevation")
            axsx.spines['right'].set_color('red')
            axsx.tick_params(axis='y', colors='red')
            axsx.yaxis.label.set_color('red')

            n = n+1
            n_sat = n_sat+1

            if n==8:
                n=0
                plt.savefig(png_file_name+"-"+str(N)+".png",bbox_inches='tight')
                plt.close()
                fig, axs = plt.subplots(4,2,figsize=(14,10),sharex=True,sharey=True)
                N = N+1
        if n!=0:
            plt.savefig(png_file_name+"-"+str(N)+".png",bbox_inches='tight')
            plt.close()
            

    else:
        fig, axs = plt.subplots(1,figsize=(10,7),sharex=True,sharey=True)
        
        i,j=0,0
        n_sat = 0
        n=0
        N=0
        min_t, max_t = min(df_station.index), max(df_station.index)
        for sat in list_sats:
            df_sat = df_station[df_station["sv"]==sat]
            axs.scatter(df_sat.index,df_sat["VTEC"].values,s=1,c=df_sat["elevation"]*2/3.1415926535,cmap="YlGnBu")

            axs.set_title(png_file_name.split("/")[-1]+" "+min_t.strftime("%d/%m/%Y %H:%M:%S") + " to " + max_t.strftime("%d/%m/%Y	 %H:%M:%S") )
            axs.grid(True)
            axs.set_ylabel("VTEC")

        print (png_file_name)
        plt.savefig(png_file_name+".png",bbox_inches='tight')

#Function to make interpolations internally 
def bi_int(df_station, data_col='VTEC',lon_margin=1,lat_margin=1,lon_density=100,lat_density=100):
    """Interpola datos (VTEC por defecto) en una grilla regular definida"""
    LAT = np.linspace(df_station['lat'].min()-lat_margin,df_station['lat'].max()+lat_margin,num=lat_density)
    LON = np.linspace(df_station['lon'].min()-lon_margin,df_station['lon'].max()+lon_margin,num=lon_density)
    grid_points = np.array ([[lat,lon] for lat in LAT for lon in LON])
    
    # Verificar suficiencia de datos
    if df_station.shape[0] < 4:
        print('Time Stamp dosent have enough datapoints.')
        return None
    interp = RBFInterpolator(df_station[['lat','lon']].values, df_station[data_col].values,
                             smoothing=0.01, kernel='thin_plate_spline')
    interpolated_values = interp(grid_points)
    return pd.DataFrame({'lat':grid_points[:,0],'lon':grid_points[:,1],data_col: interpolated_values})

#Function to plot maps directly from ionotec, uses RBFInterpolator from scipy by default
def plot_map(df_station,png_file_name,timestamp,data_col='VTEC',lon_density=100,lat_density=100,show_html=False):
    timestamp = pd.to_datetime(timestamp)
    map_data = bi_int(df_station.loc[timestamp],lon_density=lon_density,lat_density=lat_density)
    fig = px.scatter_geo(map_data,
                         lat='lat',lon='lon',
                         color=data_col,
                         color_continuous_scale='jet',
                         title = f'{data_col} at {timestamp}.',
                         projection='natural earth',
                         opacity=0.09,
                         labels={data_col:f'{data_col} [TECu]'},
                         )
    fig.update_layout(geo=dict(landcolor='rgb(212,212,212)',
                               countrycolor='rgb(255,255,255)',
                               showcountries=True,
                               showland=True),
                      title_x = 0.02,
                      title_y = 0.96,
                      margin=dict(r=10,l=10,t=10,b=10)
                      )
    fig.update_geos(projection_scale=1.75, fitbounds='locations',
                    lataxis_showgrid=True, lonaxis_showgrid=True)

    print(f"Saved image {png_file_name}")
    if show_html: fig.show()
    fig.write_image(png_file_name,width=720,height=405,scale=2)
