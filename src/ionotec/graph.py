import matplotlib.pyplot as plt
import pandas as pd
import sys
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def plot_station(df_station,png_file_name,mozaic=False):

    list_sats = df_station["sv"].unique()

    #print (png_file_name)
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
                plt.savefig(png_file_name+"-"+str(N)+".png",bbox_inches='tight')
                plt.close()
                fig, axs = plt.subplots(4,2,figsize=(28,20),sharex=True,sharey=True)
                N = N+1
        if n!=0:
            plt.savefig(png_file_name+"-"+str(N)+".png",bbox_inches='tight')
            plt.close()
            

    else:
        fig, axs = plt.subplots(1,figsize=(40,21),sharex=True,sharey=True)
        
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
            axs.set_title(png_file_name.split("/")[-1]+" "+min_t.strftime("%d/%m/%Y %H:%M:%S") + " - " + max_t.strftime("%d/%m/%Y %H:%M:%S") ,fontsize=30)
            axs.grid(True)
            axs.set_ylabel("VTEC(TECu)",fontsize=30)

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

        #print (png_file_name)
        plt.savefig(png_file_name+".png",bbox_inches='tight')
