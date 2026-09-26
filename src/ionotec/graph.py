
import matplotlib.pyplot as plt
import numpy as np
import math
import seaborn as sns
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
from matplotlib.colors import to_rgba

import colorsys

def contrast_palette(n, s=0.9, l_levels=(0.35, 0.55, 0.70)):
    phi = 0.618033988749895                    # golden ratio conjugate
    return [colorsys.hls_to_rgb((i * phi) % 1.0,
                                l_levels[i % len(l_levels)],
                                s)
            for i in range(n)]


def plot_tracks(df,dest_folder='',lat_station=None,lon_station=None,station=''):
    # Get unique 'sv' values to assign distinct colors
    unique_svs = df['sv'].unique()
    num_svs = len(unique_svs)

    datemin = min(df.index)
    datemax = max(df.index)

    str_d1 = datemin.strftime("%Y%m%d")
    str_d2 = datemax.strftime("%Y%m%d")


    # Generate a color palette with sufficient contrast for each 'sv'
    # Using sns.hls_palette directly to control saturation (s) and lightness (l)
    #colors = sns.hls_palette(num_svs, l=.5, s=.9)
    #sv_to_color = {sv: colors[i] for i, sv in enumerate(unique_svs)}

    colors = contrast_palette(num_svs)
    sv_to_color = {sv: colors[i] for i, sv in enumerate(unique_svs)}


    fig = plt.figure(figsize=(18, 9)) # Adjust figure size for better readability
    ax = fig.add_subplot(111) # Get the main axis


    # Remove plt.tight_layout to allow the inset map to overlap the main scatter plot

    if isinstance(lat_station, float) and isinstance(lon_station, float):
        # Add the world map inset (bottom-left, Mollweide projection, no ocean, no label, no borders)
        ax_map = fig.add_axes([0.14, 0.68, 0.18, 0.18], projection=ccrs.Mollweide(), zorder=1) # Adjusted position, map should be behind scatter plot

        ax.set_zorder(2)          # main axes drawn after (above) the inset
        ax.patch.set_visible(False)  # let the map show through

        ax_map.set_global() # Set global extent for Mollweide projection
        ax_map.add_feature(cfeature.LAND)
        ax_map.add_feature(cfeature.COASTLINE)
        ax_map.plot(lon_station, lat_station, 'r*', markersize=10, transform=ccrs.PlateCarree()) # Removed label



    for (sv, C1, C2), group in df.groupby(['sv', 'C1', 'C2']):
        alpha = np.clip(group['elevation'].to_numpy() / (180), 0, 1)

        r, g, b, _ = to_rgba(sv_to_color[sv])
        rgba = np.column_stack([
            np.full(len(alpha), r),
            np.full(len(alpha), g),
            np.full(len(alpha), b),
            alpha,
        ])

        ax.scatter(group.index, group['VTEC'],
               c=rgba,          # per-point RGBA, alpha included
               s=18,
               linewidths=0,    # avoid edge colors ignoring alpha
               zorder=1)

    # Iterate through each unique combination of (sv, C1, C2)
    #for (sv, C1, C2), group in df.groupby(['sv', 'C1', 'C2']):
    #    # Calculate alpha based on elevation (normalized from 0 to pi/2)
    #    group_alpha = group['elevation'] / (math.pi / 2)
    #    # Clip alpha values to ensure they are within the valid range [0, 1]
    #    group_alpha = np.clip(group_alpha, 0, 1)
#
#        # Get the color assigned to the current satellite (sv)
#        color = sv_to_color[sv]
#
#        # Plot VTEC vs. time for the current group with calculated alpha and color
#        ax.scatter(group.index, group['VTEC'],
#                c=[color], # 'c' expects a sequence of colors or a single color
#                alpha=group_alpha,
#                s=18,
#                zorder=1) # Scatter plot should be on top

    # Plot modifications
    ax.set_xlabel('') # No xlabel as requested
    ax.set_ylabel('VTEC (TECu)', fontsize=14) # Changed ylabel and increased font size
    ax.set_title(station, fontsize=16) # Changed title and increased font size
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.tick_params(axis='x', rotation=0, labelsize=12) # No rotation, increased font size
    ax.tick_params(axis='y', labelsize=12) # Increased font size

    # Use DateFormatter for x-axis ticks (only month and day)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    # Use AutoDateLocator to ensure appropriate tick spacing
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())



    plt.savefig(dest_folder+'/simple_'+station+'_'+str_d1+'_'+str_d2+'.png',bbox_inches='tight')
    plt.close()









def plot_tracks_individuals(df,dest_folder='',station=''):

    # Get unique 'sv' values to assign distinct colors
    unique_svs = df['sv'].unique()
    num_svs = len(unique_svs)

    datemin = min(df.index)
    datemax = max(df.index)

    str_d1 = datemin.strftime("%Y%m%d")
    str_d2 = datemax.strftime("%Y%m%d")



    colors = contrast_palette(num_svs)
    sv_to_color = {sv: colors[i] for i, sv in enumerate(unique_svs)}

    # --- Determine global x and y limits for shared axes ---
    xmin = df.index.min()
    xmax = df.index.max()
    x_margin = 0.02*(xmax-xmin)
    xmin = xmin - x_margin
    xmax = xmax + x_margin

    ymin = df['VTEC'].min()
    ymax = df['VTEC'].max()
    y_margin = 0.1*(ymax-ymin)
    ymin = ymin - y_margin
    ymax = ymax + y_margin

    # --- Determine layout for subplots ---
    svs_per_small_plot = 4
    num_small_plot_rows = math.ceil(num_svs / svs_per_small_plot)
    num_cols_for_small_plots = 2 # Let's use 2 columns for the smaller plots

    # Total rows for gridspec: 1 for main plot + num_small_plot_rows
    num_rows_total = 1 + num_small_plot_rows

    fig = plt.figure(figsize=(25, 5 * num_rows_total)) # Adjust figure size dynamically

    gs = gridspec.GridSpec(num_rows_total, num_cols_for_small_plots, figure=fig)

    # --- Main axis for all tracks (without map) ---
    ax_main = fig.add_subplot(gs[0, :]) # Spans all columns in the first row

    for (sv, C1, C2), group in df.groupby(['sv', 'C1', 'C2']):
        group_alpha = np.clip(group['elevation'] / (math.pi / 2), 0, 1)
        color = sv_to_color[sv]
        ax_main.scatter(group.index, group['VTEC'],
                    c=[color],
                    alpha=group_alpha,
                    s=10)

    ax_main.set_xlabel('')
    ax_main.set_ylabel('VTEC (TECu)', fontsize=20) # Increased font size
    ax_main.set_title('mdo1', fontsize=25) # Increased font size
    ax_main.grid(True, linestyle='--', alpha=0.7)
    ax_main.tick_params(axis='x', rotation=0, labelsize=20)
    ax_main.tick_params(axis='y', labelsize=20) # Increased font size
    ax_main.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    ax_main.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax_main.set_xlim(xmin, xmax) # Apply global x-limits
    ax_main.set_ylim(ymin, ymax) # Apply global y-limits


    # --- Create and populate smaller axes (4 svs each) ---
    row_idx = 1 # Start from the second row for small plots
    col_idx = 0

    for i in range(0, num_svs, svs_per_small_plot):
        current_svs = unique_svs[i : i + svs_per_small_plot]

        # Create subplot for this group of SVs
        ax_small = fig.add_subplot(gs[row_idx, col_idx])

        legend_handles = []

        for sv_to_plot in current_svs:
            # Filter for the specific sv_to_plot and iterate through its C1, C2 combinations
            sv_groups = df[df['sv'] == sv_to_plot].groupby(['sv', 'C1', 'C2'])
    
            for (sv, C1, C2), group in sv_groups:
                group_alpha = np.clip(group['elevation'] / (math.pi / 2), 0, 1)
                color = sv_to_color[sv]
    
                ax_small.scatter(group.index, group['VTEC'],
                                 c=[color],
                                 alpha=group_alpha,
                                 s=10)
    
            # Create a proxy artist for the legend entry for each SV
            color = sv_to_color[sv_to_plot]
            legend_handles.append(plt.Line2D([0], [0], marker='o', color='w', label=sv_to_plot,
                                             markerfacecolor=color, markersize=8))
    
        # Removed title='SV' from legend and added fontsize
        ax_small.legend(handles=legend_handles, loc='upper right', bbox_to_anchor=(1.0, 1.0), ncol=len(current_svs), fontsize=15)
    
        # No title for subaxes as requested
        # ax_small.set_title(f'SVs: {', '.join(current_svs)}', fontsize=12)
    
        # Set y-label conditionally
        if col_idx == 0: # Only for the left column
            ax_small.set_ylabel('VTEC (TECu)', fontsize=18) # Increased font size
            ax_small.tick_params(axis='y', labelsize=18) # Increased font size
        else:
            ax_small.set_ylabel('') # No ylabel for right column
            ax_small.tick_params(axis='y', labelleft=False) # No y-ticks for right column
    
        ax_small.grid(True, linestyle='--', alpha=0.7)
    
        # Set x-axis labels/ticks on all subplots as requested
        ax_small.tick_params(axis='x', rotation=0, labelsize=18) # Increased font size, always show
        ax_small.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
        ax_small.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax_small.set_xlabel('') # Ensure no xlabel even for bottom plots
    
        ax_small.set_xlim(xmin, xmax) # Apply global x-limits
        ax_small.set_ylim(ymin, ymax) # Apply global y-limits
    
        # Move to the next column/row for subplot placement
        col_idx += 1
        if col_idx >= num_cols_for_small_plots:
            col_idx = 0
            row_idx += 1
    
    plt.tight_layout() # Adjust layout to prevent overlaps
    plt.savefig(dest_folder+'/full_'+station+'_'+str_d1+'_'+str_d2+'.png',bbox_inches='tight')
    plt.close()
