import matplotlib.pyplot as plt
import pandas as pd
import sys
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import plotly.express as px
from ionotec import stations

from pathlib import Path
import cartopy.crs as ccrs
import cartopy. feature as cfeature
from cartopy.mpl.geoaxes import GeoAxes
#import cartopy.io.img_tiles as cimgt

from scipy.interpolate import RBFInterpolator

# create output dir
base_dir = Path(sys.argv[0]).resolve().parent
target_dir = base_dir / "output/TEC/Figs"
target_dir.mkdir(parents=True, exist_ok=True)
target_dir = str(target_dir)

def plot_station(df_station,station_name,mozaic=False):

    list_sats = df_station["sv"].unique()

    #print (station_name)
    #print (list_sats)

    if mozaic:
        fig, axs = plt.subplots(4,2,figsize=(28,20),sharex=True,sharey=True)

        i,j=0,0
        n_sat = 0
        n=0
        N=0
        for sat in list_sats:
            i = int(n/2)
            j = n%2
            df_sat = df_station[df_station["sv"]==sat]
            axs[i,j].plot(df_sat["VTEC"],'b.',markersize=2)
            axsx = axs[i,j].twinx()
            axsx.plot(df_sat["elevation"]*180/3.1415926535,'r.',markersize=2,alpha=0.5)
            axsx.set_ylim([0,90])
            axs[i,j].set_title(sat,fontsize=15)
            axs[i,j].grid(True)
            axs[i,j].set_ylabel("VTEC",fontsize=15)
            axsx.set_ylabel("elevation")
            axsx.spines['right'].set_color('red')
            axsx.tick_params(axis='y', colors='red')
            axsx.yaxis.label.set_color('red')

            n = n+1
            n_sat = n_sat+1

            if n==8:
                n=0
                plt.savefig(station_name+"-"+str(N)+".png",bbox_inches='tight')
                plt.close()
                fig, axs = plt.subplots(4,2,figsize=(28,20),sharex=True,sharey=True)
                N = N+1
        if n!=0:
            plt.savefig(station_name+"-"+str(N)+".png",bbox_inches='tight')
            plt.close()
            

    else:
        fig, axs = plt.subplots(1,figsize=(40,21),sharex=True,sharey=True)

        #station_name = station_name.split("/")[-1]
        print (station_name)
        i,j=0,0
        n_sat = 0
        n=0
        N=0
        min_t, max_t = min(df_station.index), max(df_station.index)
        first_GPS = True
        first_GLONASS = True
        for sat in list_sats:
            df_sat = df_station[df_station["sv"]==sat]
            if sat[0]=="G": 
                if first_GPS: 
                    sc1 = axs.scatter(df_sat.index,df_sat["VTEC"].values,s=25,c=df_sat["elevation"]*180/np.pi,cmap="YlGnBu",vmin=0,vmax=90)
                    first_GPS=False
                else:
                    axs.scatter(df_sat.index,df_sat["VTEC"].values,s=25,c=df_sat["elevation"]*180/np.pi,cmap="YlGnBu",vmin=0,vmax=90)
            if sat[0]=="R": 
                if first_GLONASS: 
                    sc2 = axs.scatter(df_sat.index,df_sat["VTEC"].values,s=25,c=df_sat["elevation"]*180/np.pi,cmap="PuRd",vmin=0,vmax=90)
                    first_GLONASS=False
                else:
                    axs.scatter(df_sat.index,df_sat["VTEC"].values,s=25,c=df_sat["elevation"]*180/np.pi,cmap="PuRd",vmin=0,vmax=90)
            axs.tick_params(axis='both', which='major', labelsize=25)
            axs.tick_params(axis='both', which='minor', labelsize=25)
            axs.set_title(station_name+" "+min_t.strftime("%d/%m/%Y %H:%M:%S") + " - " + max_t.strftime("%d/%m/%Y %H:%M:%S") ,fontsize=30)
            axs.grid(True)
            axs.set_ylabel("VTEC(TECu)",fontsize=30)

        ## Add station location on the plot
        #ax_map = inset_axes(axs, width="20%", height="25%", loc='upper right', borderpad=2)
        #ax_inset.set_extent([-90, -30, -60, 15]) 
        #ax_inset.set_xticks([])
        #ax_inset.set_yticks([])
 
        df_stations = pd.read_csv(stations.csv_stations).set_index("station")
        pos = df_stations.loc[station_name]

        #request = cimgt.QuadtreeTiles()
        #sys.exit()
        ax_map = inset_axes(
            axs,
            width="20%",      # very narrow
            height="35%",    # short bar
            loc="upper left",
            bbox_to_anchor=(0.02,-0.05,1.0,1.0),
            borderpad=0,
            bbox_transform=axs.transAxes,
            axes_class=GeoAxes,
            axes_kwargs=dict(projection=ccrs.PlateCarree())
            #axes_kwargs=dict(projection=request.crs)
        
        )

        ax_map.patch.set_alpha(0.0)
        ax_map.set_facecolor('none')
        
        #ax_map.add_image(request, 8)
        #ax_map.set_title("Station Location")
        #ax_map = plt.axes(projection=ccrs.PlateCarree())
    
        # 3. Add geographic features
        #ax_map.add_feature(cfeature.COASTLINE)
        #ax_map.add_feature(cfeature.BORDERS, linestyle=':')
        #ax_map.add_feature(cfeature.LAND, facecolor='lightgray', alpha=0.3)
        ax_map.add_feature(cfeature.OCEAN, facecolor='none',alpha=0.0) # Light blue
        ax_map.add_feature(cfeature.LAND, facecolor='#F5F5F5', edgecolor='none') # Soft gray-white
        ax_map.add_feature(cfeature.COASTLINE, linewidth=0.6, edgecolor='#555555')
        ax_map.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5, edgecolor='#999999')
        ax_map.set_extent([-90, -30, -60, 15]) 

        ax_map.spines['geo'].set_visible(False)
        #ax_map.stock_img()

        #ax_map.add_feature(cfeature.OCEAN, facecolor='none',alpha=0.0) # Light blue

        

        ax_map.scatter(pos['lon'], pos['lat'], color='darkred', 
               s=80, marker='X', transform=ccrs.PlateCarree(), 
               label='Observation points', zorder=5)

        ax_map.text(pos['lon'], pos['lat'] + 2, station_name,
                    transform=ccrs.PlateCarree(),
                    fontsize=24,
                    fontweight='bold',
                    ha='center',        # Horizontal alignment: center
                    va='bottom',        # Vertical alignment: bottom
                    color='black')
        
        if not first_GPS:
            cax1 = inset_axes(
                axs,
                width="15%",      # very narrow
                height="5%",    # short bar
                loc="upper right",
                bbox_to_anchor=(0,0,0.99,0.99),
                borderpad=0,
                bbox_transform=axs.transAxes
            )
            cbar1 = fig.colorbar(sc1, cax=cax1, orientation='horizontal',location='top')
            #cbar = fig.colorbar(sc, cax=cax)
            cbar1.ax.set_title("GPS elevation",fontsize=25,y=0.17)

            if not first_GLONASS:
                cbar1.ax.tick_params(length=0)
                cbar1.set_ticks([])
            #cbar1.ax.xaxis.set_label_position('left')

        if not first_GLONASS:
            cax2 = inset_axes(
                axs,
                width="15%",      # very narrow
                height="5%",    # short bar
                loc="upper right",
                bbox_to_anchor=(0,-0.055,0.99,0.99),
                borderpad=0,
                bbox_transform=axs.transAxes
            )
            
        


            cbar2 = fig.colorbar(sc2, cax=cax2, orientation='horizontal',location='bottom')
            #cbar = fig.colorbar(sc, cax=cax)
            cbar2.ax.set_title("GLONASS elevation",fontsize=25,y=0.17)
        

        #cbar2 = fig.colorbar(sc2, ax=cax,pad=0.04)
        #cbar2.set_label("GLONASS elevation",fontsize=30)

        #cbar1.ax.tick_params(labelsize=10, length=3)
        #cbar2.ax.tick_params(labelsize=10, length=3)
        
        # Increase tick size
        #cbar1.ax.tick_params(labelsize=25)
        #cbar2.ax.tick_params(labelsize=25)

        # Optional: increase tick width & length
        #cbar1.ax.tick_params(width=2, length=6)
        #cbar2.ax.tick_params(width=2, length=6)

        #print (station_name)
        plt.savefig(target_dir+"/"+station_name+".png",bbox_inches='tight')

