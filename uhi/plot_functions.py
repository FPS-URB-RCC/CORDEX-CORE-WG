import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import pandas as pd
import xarray as xr
from itertools import product
from shapely.geometry import Point, Polygon

var_map = {
    'tasmin': 'TMIN',
    'tasmax': 'TMAX'
}

def plot_climatology(ds, ucdb_city, urban_vicinity, variable, URBAN, 
                     valid_stations = None, time_series = None, city = None):
    """
    Plot the climatological data.

    Parameters:
        ds (xr.Dataset): Dataset containing the climatological data.
        ucdb_city (gpd.GeoDataFrame): GeoDataFrame of the city boundaries.
        urban_vicinity (object): Object representing urban vicinity.
        obs (pd.DataFrame, optional): DataFrame containing observational data (default is None).

    Returns:
        matplotlib.figure.Figure: The generated figure.
    """
    ds_var_period_mean = ds.mean('time').compute()
    rural_mean = ds_var_period_mean[variable].where(urban_vicinity['urmask'] == 0).mean().compute()
    data = ds_var_period_mean[variable] - rural_mean

    proj = ccrs.PlateCarree()
    fig, ax = plt.subplots(subplot_kw={'projection': proj}, figsize=(12, 6))
    
    # Compute the maximum absolute value
    max_abs_value = abs(data).max().item()
    
    if valid_stations is not None:
        obs_urban, obs_vicinity = [], []
        for index, obs in valid_stations.iterrows():
            obs_lon = obs['lon']
            obs_lat = obs['lat']
            
            # Create a point with latitude and longuitude
            point = Point(obs_lon, obs_lat)
            is_inside = ucdb_city.contains(point)

            # Clasify in urban and vicinity
            if is_inside.any():
                obs_urban.append(obs['code'])
            else:
                obs_vicinity.append(obs['code'])
        if len(obs_vicinity)>0:
            #calculate the mean for the outside observations
            obs_vicinity_mean = 0 
            for item in time_series:
                if item['code'] in obs_vicinity:
                    obs_vicinity_mean+=item['data'].mean()[0]
            obs_vicinity_mean = obs_vicinity_mean / len(obs_vicinity)
            #plot each observation
            for item in time_series:
                temp_obs=item['data'].mean()[0]-obs_vicinity_mean
                if abs(temp_obs)>max_abs_value:
                    max_abs_value=abs(temp_obs)
            for index, valid_stations in valid_stations.iterrows():
                obs_lon = valid_stations['lon']
                obs_lat = valid_stations['lat']
                for item in time_series:
                    if item['code'] == valid_stations['code']:
                        temp_obs=item['data'].mean()[0]-obs_vicinity_mean
                        ax.scatter(obs_lon, obs_lat, c = temp_obs, marker='o', cmap='bwr', 
                                   s = 40, edgecolors = 'gray', vmin = -max_abs_value, vmax = max_abs_value,
                                   zorder = 10000) 
    
    im1 = ax.pcolormesh(ds.lon, ds.lat, data.values,
                    cmap='bwr', alpha = 0.7,
                    vmin = - max_abs_value, 
                    vmax = max_abs_value)
    
    cbar = fig.colorbar(im1, ax = ax)
    cbar.set_label('°C', rotation = 90, fontsize = 14)
    
    ucdb_city.plot(ax=ax, facecolor="none", transform=proj, edgecolor="Green", linewidth=2, zorder = 1000)
    
    ax.coastlines()
    if variable == 'tasmin':
        ax.set_title(f"Minimum temperature anomaly for {city}", fontsize = 14)

    # Overlay the cell borders and handle NaNs
    URBAN.plot_urban_borders(urban_vicinity, ax)
    
    plt.subplots_adjust(wspace=0.1, hspace=0.1)
    
    return fig

def plot_time_series(ds_var, variable, urban_vicinity, 
                     time_series = None, valid_stations = None, data_squares = False,
                     percentile = 100, var_map = var_map, ucdb_city = None, city = None, cache = ''):
    '''
    Plot time series data with optional urban area overlay and additional time series overlay.

    Parameters:
    ds_var (xarray.Dataset): Dataset containing the variable data.
    variable (str): Name of the variable of interest.
    urban_vicinity (xarray.Dataset): Dataset containing information about urban areas.
    time_series (list of pandas.DataFrame, optional): List of time series dataframes to overlay on the plot.
    data_squares (bool, optional): Flag indicating whether to plot individual data squares for urban and rural areas.

    Returns:
    matplotlib.figure.Figure: The plotted figure.
    '''
    urban_area_legend = False
    not_urban_area_legend = False
    is_rural = urban_vicinity['urmask'] == 0
    is_urban = urban_vicinity['urmask'] == 1
    rural_mean = (ds_var[variable]
        .where(is_rural)
        .groupby('time.month')
        .mean(dim = [ds_var.cf['Y'].name, ds_var.cf['X'].name, 'time'])
        .compute()
    )
    urban = ds_var[variable].where(is_urban).groupby('time.month')
    urban_mean = urban.mean(dim=[ds_var.cf['Y'].name, 
                                 ds_var.cf['X'].name,'time']).compute()
                         
    ds_var_period_mean = ds_var.groupby('time.month').mean('time')                  
    ds_annomaly = ds_var_period_mean[variable] - rural_mean
    rural_anomaly = ds_annomaly.where(is_rural)
    urban_anomaly = ds_annomaly.where(is_urban)
    xr.Dataset(dict(
        rural_anomaly = rural_anomaly,
        urban_anomaly = urban_anomaly,
        rural_mean = rural_mean,
        urban_mean = urban_mean
    )).to_netcdf(cache)
                         
    # Plot mean annual cycle (urban and rural)
    fig, ax = plt.subplots(figsize=(15, 7)) 
    (urban_mean-rural_mean).plot(ax=ax,  color = 'r', linestyle='-', 
                                     linewidth = 4, label='Urban mean')
                         
    (rural_mean-rural_mean).plot(ax=ax,  color = 'b', linestyle='-', 
                                 linewidth = 4, label='Vicinity mean')
    #
    # Plot individual data squares for urban and rural areas if requested
    #
    if data_squares:
        #
        # Fill within percentiles
        #
        axis = [ds_annomaly.get_axis_num(ds_annomaly.cf['X'].name),
                ds_annomaly.get_axis_num(ds_annomaly.cf['Y'].name)]
        colors = ['blue', 'red']
        for index, anom in enumerate([rural_anomaly, urban_anomaly]): 
            lower_percentile = np.nanpercentile(anom, percentile, axis=axis)
            upper_percentile = np.nanpercentile(anom, 100-percentile, axis=axis)
            ax.fill_between(
                ds_var_period_mean['month'],
                lower_percentile, upper_percentile,
                color=colors[index], alpha=0.1
            )
            for i, j in product(anom.cf['X'].values, anom.cf['Y'].values):
                anom_val = anom.sel({ds_var.cf['X'].name:i, ds_var.cf['Y'].name:j})
                if not np.isnan(anom_val[0]):
                    anom_val.plot(ax=ax, color=colors[index], linewidth=0.5)
                         
        #Add manually the legend
        #if urban_area_legend==True:
        #    ax.plot([], [], color='r', linewidth=0.5, label = 'Urban Area')
        #if not_urban_area_legend==True:
        #    ax.plot([], [], color='b', linewidth=0.5, label = 'Vicinity Area')
    
    #Plot the observation if requested
    obs_urban, obs_vicinity = [], []
    if time_series is not None:
        urban_obs_legend = False

        var = var_map.get(variable, None)
        #obs_monthly_change_mean = [0] * 12
        obs_monthly_change_mean_urban = [0] * 12
        obs_monthly_change_mean_vicinity = [0] * 12

        for index, obs in valid_stations.iterrows():
            obs_lon = obs['lon']
            obs_lat = obs['lat']
            
            # Create a point with latitude and longuitude
            point = Point(obs_lon, obs_lat)
            is_inside = ucdb_city.contains(point)
            
            # Clasify in urban and vicinity
            if is_inside.any():
                obs_urban.append(obs['code'])
            else:
                obs_vicinity.append(obs['code'])
        
        #get the mean and plot not urban observatios
        if len(obs_vicinity) > 0:
            for item in time_series:
                if item['code'] in obs_vicinity:
                    time_series_df = pd.DataFrame(item['data'])
                    time_series_df.index = pd.to_datetime(time_series_df.index)
                    time_series_df['month'] = time_series_df.index.month
                    for i in range(1, 13):
                        monthly_data = time_series_df[var].loc[time_series_df.index.month == i].mean()
                        obs_monthly_change_mean_vicinity[i-1] += monthly_data
            
            obs_monthly_change_mean_vicinity = [x / len(obs_vicinity) for x in obs_monthly_change_mean_vicinity]
            plt.plot(range(1, 13), [0] * 12 , 
                     color='g', linestyle='-', linewidth = 4, label='Vicinity obs. mean', zorder = 2000)

            #Plot each stations individually
            for item in time_series:
                time_series_df = pd.DataFrame(item['data'])
                time_series_df.index = pd.to_datetime(time_series_df.index)
                time_series_df['month'] = time_series_df.index.month
                obs_monthly_change = []
                for i in range(1, 13):
                    monthly_data = time_series_df.loc[time_series_df.index.month == i].mean()
                    monthly_change = monthly_data[var] - obs_monthly_change_mean_vicinity[i-1]
                    obs_monthly_change.append(monthly_change)
                    if item['code'] in obs_urban:
                        obs_monthly_change_mean_urban[i-1] += monthly_change
                        color_obs='g'
                        urban_obs_legend=True
                    else:
                        color_obs='k'                       
                        
                plt.plot(range(1, 13), obs_monthly_change, marker='o', color=color_obs, 
                             linestyle='--', linewidth = 2)
    
            
            if len(obs_urban) > 0:
                obs_monthly_change_mean = [x / len(obs_urban) for x in obs_monthly_change_mean_urban]
                plt.plot(range(1, 13), obs_monthly_change_mean, 
                         color='k', linestyle='-', linewidth = 4, label='Urban obs. mean', zorder = 2000)
                        
            #Add manually the legend
            ax.plot([], [], color='g', marker='o',  linestyle='--', linewidth=2, label='Vicinity obs.')
            if urban_obs_legend==True:
                ax.plot([], [], color='k', marker='o', linestyle='--',  linewidth=2, label='Urban obs.')
            
        else:
            print("Due to limited data availability outside urban areas, we are currently unable to present observations for any region")
 
    # Add legend to the plot
    ax.legend(fontsize = 14)
    
    # Customize the plot
    #ax.set_xlabel('Month', fontsize = 18)
    if variable == 'tasmin':
        ax.set_title(f"Minimum temperature anomaly for {city}", fontsize = 18)
        ax.set_ylabel(f"Minimum temperature anomaly (°C)", fontsize = 18)
    ax.set_xticks(np.arange(1, 13))
    ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 
                        'Aug', 'Sep', 'Oct', 'Nov', 'Dec'], fontsize = 18)
    ax.tick_params(axis='y', labelsize=18)
    #plt.grid(axis='y', linestyle='--', color='gray', alpha=0.7)
    
    return fig