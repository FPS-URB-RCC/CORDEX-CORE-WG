import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import pandas as pd
import xarray as xr
from itertools import product

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
    fig, ax = plt.subplots(subplot_kw={'projection': proj}, figsize=(20, 10))
    
    # Compute the maximum absolute value
    max_abs_value = abs(data).max().item()
    
    if valid_stations is not None:
        for item in time_series:
            temp_obs=item['data'].mean()[0]-rural_mean
            if abs(temp_obs)>max_abs_value:
                max_abs_value=abs(temp_obs)
        for index, valid_stations in valid_stations.iterrows():
            obs_lon = valid_stations['lon']
            obs_lat = valid_stations['lat']
            for item in time_series:
                if item['code'] == valid_stations['code']:
                    temp_obs=item['data'].mean()[0]-rural_mean
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
                     percentile = 100, var_map = var_map, city = None, cache = ''):
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
    fig, ax = plt.subplots(figsize=(20, 10)) 
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
    obs_urban, obs_vicinity, obs_not_sected = [], [], []
    if time_series is not None:
        urban_obs_legend = False
        not_urban_obs_legend = False
        not_obs_legend = False
        var = var_map.get(variable, None)
        #obs_monthly_change_mean = [0] * 12
        obs_monthly_change_mean_urban = [0] * 12
        obs_monthly_change_mean_vicinity = [0] * 12
        obs_monthly_change_mean_not_sected = [0] * 12

        for index, obs in valid_stations.iterrows():
            obs_lon = obs['lon']
            obs_lat = obs['lat']
            
            # Calculate the differences
            dist = np.sqrt((urban_vicinity['lon'] - obs_lon)**2 + 
                           (urban_vicinity['lat'] - obs_lat)**2)
            min_dist_idx = np.unravel_index(np.argmin(dist.values, axis=None), dist.shape)
            
            # Select the urban mask value at the nearest grid point
            selected_value = urban_vicinity['urmask'].isel({
                urban_vicinity.cf['Y'].name : min_dist_idx[0],
                urban_vicinity.cf['X'].name : min_dist_idx[1]}).values
            urban_mask_value = selected_value.item()
            if selected_value.item()==1:
                color_obs='k'
                urban_obs_legend=True
                obs_urban.append(obs)
            elif selected_value.item()==0:
                color_obs='g'
                not_urban_obs_legend=True
                obs_vicinity.append(obs)
            else:
                color_obs='grey'
                not_obs_legend=True
                obs_not_sected.append(obs)
            for item in time_series:
                if item['code'] == obs['code']:                   
                    time_series_df = pd.DataFrame(item['data'])
                    time_series_df.index = pd.to_datetime(time_series_df.index)
                    time_series_df['month'] = time_series_df.index.month
                    obs_monthly_change = []
                    for i in range(1, 13):
                        monthly_data = time_series_df.loc[time_series_df.index.month == i].mean()
                        rural_data = rural_mean[i-1].values
                        monthly_change = monthly_data[var] - rural_data
                        obs_monthly_change.append(monthly_change)
                        #obs_monthly_change_mean[i-1] += monthly_change
                        
                        if color_obs == 'k':
                            obs_monthly_change_mean_urban[i-1] += monthly_change
                        if color_obs == 'g':
                            obs_monthly_change_mean_vicinity[i-1] += monthly_change
                        if color_obs == 'grey':
                            obs_monthly_change_mean_not_sected[i-1] += monthly_change
                        
                    plt.plot(range(1, 13), obs_monthly_change, marker='o', color=color_obs, 
                             linestyle='--', linewidth = 2)
                    
                    
        #Add manually the legend
        if urban_obs_legend==True:
            ax.plot([], [], color='k', marker='o',  linewidth=0.5, label='Urban obs.')
        if not_urban_obs_legend==True:
            ax.plot([], [], color='g', marker='o',  linewidth=0.5, label='Vicinity obs.')
        if not_obs_legend==True:
            ax.plot([], [], color='grey', marker='o',  linewidth=0.5, 
                    label='Outside study area obs.')

        
        if len(obs_urban) > 0:
            obs_monthly_change_mean = [x / len(obs_urban) for x in obs_monthly_change_mean_urban]
            plt.plot(range(1, 13), obs_monthly_change_mean, 
                     color='k', linestyle='-', linewidth = 4, label='Urban obs. mean', zorder = 2000) 
        if len(obs_vicinity) > 0:
            obs_monthly_change_mean = [x / len(obs_vicinity) for x in obs_monthly_change_mean_vicinity]
            plt.plot(range(1, 13), obs_monthly_change_mean, 
                     color='g', linestyle='-', linewidth = 4, label='Vicinity obs. mean', zorder = 2000) 
    
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