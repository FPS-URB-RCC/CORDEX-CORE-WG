### This script aims to create a mask for urban and its surroundings of a city
### UF >0.1 to define cities
### CORDEX-CORE analysis // Gaby Langendijk 

import netCDF4 as nc
import matplotlib.pyplot as plt
import xarray as xr

def create_masks(city, domain, center_lat, center_lon):
    file_path = f"D:\\cordexcore22\\CORDEX-CORE-WG\\interp\\sfturf_{domain}-22n_ECMWF-ERAINT_evaluation_r1i1p1_GERICS-REMO2015_v1_fx.nc"
    file_path2 = f"D:\\cordexcore22\\CORDEX-CORE-WG\\interp\\sftlf_{domain}-22n_ECMWF-ERAINT_evaluation_r1i1p1_GERICS-REMO2015_v1_fx.nc"
    

    file_path3 = f"D:\\cordexcore22\\CORDEX-CORE-WG\\interp\\sfturf_{domain}-22n_ECMWF-ERAINT_evaluation_r1i1p1_ICTP-RegCM4_v0_fx.nc"
    file_path4 = f"D:\\cordexcore22\\CORDEX-CORE-WG\\interp\\sftlf_{domain}-22n_ECMWF-ERAINT_evaluation_r1i1p1_ICTP-RegCM4_v0_fx.nc"
    
    ## Orography
    file_path5 = f"D:\\cordexcore22\\CORDEX-CORE-WG\\interp\\orog_{domain}-22n_ECMWF-ERAINT_evaluation_r1i1p1_GERICS-REMO2015_v1_fx.nc"

    # Open the netCDF files
    ds = xr.open_dataset(file_path)
    ds2 = xr.open_dataset(file_path2)
    ds3 = xr.open_dataset(file_path3)
    ds4 = xr.open_dataset(file_path4)
    ds5 = xr.open_dataset(file_path5)

    # Select the orog variable (oroography)
    orog = ds5['orog']

    # Select the "sftlf" variable (land fraction)
    sftlf_REMO = ds2['sftlf']
    #sftlf_RegCM = ds4['sftlf']

    # Define the size of the area you want to select (in degrees)
    size = 1.8   # Convert from km to degrees given the resolution

    # Select the area around the center
    sftlf_area_REMO = sftlf_REMO.sel(lat=slice(center_lat - size/2, center_lat + size/2),
                           lon=slice(center_lon - size/2, center_lon + size/2))
    #sftlf_area_RegCM = sftlf_RegCM.sel(lat=slice(center_lat - size/2, center_lat + size/2),
                           #lon=slice(center_lon - size/2, center_lon + size/2))

    # Filter out the grid boxes where "sftlf" is less than 20% for large water bodies
    sftlf_area_filtered_REMO = sftlf_area_REMO.where(sftlf_area_REMO >= 20)
    #sftlf_area_filtered_RegCM = sftlf_area_RegCM.where(sftlf_area_RegCM >= 20)

    # Select the "sftuf" variable
    sftuf_REMO = ds['urban']
    #sftuf_RegCM = ds3['sftuf']

    # Since it's a 3D variable with a time dimension of 1, select the first time step
    sftuf_REMO = sftuf_REMO.isel(time=0)
    #sftuf_RegCM = sftuf_RegCM.isel(time=0)

    # Select the area around the center
    sftuf_area_REMO = sftuf_REMO.sel(lat=slice(center_lat - size/2, center_lat + size/2),
                           lon=slice(center_lon - size/2, center_lon + size/2))
    #sftuf_area_RegCM = sftuf_RegCM.sel(lat=slice(center_lat - size/2, center_lat + size/2),
    #                        lon=slice(center_lon - size/2, center_lon + size/2))

    # Filter out the grid boxes where "sftuf" is less than 0.1
    sftuf_area_filtered_REMO = sftuf_area_REMO.where(sftuf_area_REMO > 0.1)
    #sftuf_area_filtered_RegCM = sftuf_area_RegCM.where(sftuf_area_RegCM > 0.1)

    # Filter out the grid boxes where "sftuf_area_filtered" and "sftlf_area_filtered" are larger than 0
    sftuf_area_filtered_positive_REMO = sftuf_area_filtered_REMO.where(sftuf_area_filtered_REMO > 0)
    #sftuf_area_filtered_positive_RegCM = sftuf_area_filtered_RegCM.where(sftuf_area_filtered_RegCM > 0)
    
    sftlf_area_filtered_positive_REMO = sftlf_area_filtered_REMO.where(sftlf_area_filtered_REMO > 0)
    #sftlf_area_filtered_positive_RegCM = sftlf_area_filtered_RegCM.where(sftlf_area_filtered_RegCM > 0)


    # Create a new dataset with the remaining grid boxes
    positive_gridboxes_REMO = xr.Dataset({
        'sftuf': sftuf_area_filtered_positive_REMO,
        'sftlf': sftlf_area_filtered_positive_REMO
    })

    mask_urban_REMO = positive_gridboxes_REMO['sftuf']

    # Select the area around the center for 'orog'
    orog_area_REMO = orog.sel(lat=slice(center_lat - size/2, center_lat + size/2),
                              lon=slice(center_lon - size/2, center_lon + size/2))


    # Get the orography at the city location
    orog_city = orog_area_REMO.sel(lat=center_lat, lon=center_lon, method='nearest')
    # Calculate the difference in orography between the city location and the surrounding grid boxes
    orog_diff = abs(orog_area_REMO - orog_city)
    # Create a mask where the difference in orography is not larger than 600 meters between city location and suurounding grid boxes
    orog_diff_mask = xr.where(orog_diff <= 600, 1, 0)


   # positive_gridboxes_RegCM = xr.Dataset({
   #     'sftuf': sftuf_area_filtered_positive_RegCM,
   #     'sftlf': sftlf_area_filtered_positive_RegCM
   # })

    # mask_urban_RegCM = positive_gridboxes_RegCM['sftuf']

    # Create a mask where "sftuf_area" is less than 0.1
    mask_REMO = sftuf_area_REMO < 0.1
    #mask_RegCM = sftuf_area_RegCM < 0.1

    # Apply the mask to "sftlf_area_filtered" grid boxes
    #mask_surroundings_REMO = sftlf_area_filtered_REMO.where(mask_REMO)
    # Apply the mask to "sftlf_area_filtered" grid boxes + orography 
    mask_surroundings_REMO = sftlf_area_filtered_REMO.where(mask_REMO & (orog_diff_mask == 1))
    #mask_surroundings_RegCM = sftlf_area_filtered_RegCM.where(mask_RegCM)

    # Plot the mask_urban
    mask_urban_REMO.plot(cmap='viridis', vmin=0, vmax=1)
    plt.title(f'Mask urban REMO - {city} with UF > 0.1 and no sea, excl high mountains')
    plt.show()

    # Plot the mask_surroundings
    mask_surroundings_REMO.plot()
    plt.title(f'Mask surroundings REMO - {city} with UF <= 0.1, no sea, excl high mountains')
    plt.show()


    # Plot the mask_urban
    #mask_urban_RegCM.plot(cmap='viridis', vmin=0, vmax=1)
    #plt.title(f'Mask urban RegCM - {city} with UF > 0.1 and no sea')
    #plt.show()

    # Plot the mask_surroundings
    #mask_surroundings_RegCM.plot()
    #plt.title(f'Mask surroundings RegCM - {city} with UF <= 0.1 and no sea')
    #plt.show()


    return mask_urban_REMO, mask_surroundings_REMO #, mask_urban_RegCM, mask_surroundings_RegCM

#create_masks("Johannesburg", "AFR", -26.20, 28.04) # Johannesburg
#create_masks("Tokyo", "EAS", 35.65, 139.83) # Tokyo


## list of dictionaries to store the city information
cities = [
    {"name": "Tokyo", "domain": "EAS", "lat": 35.65, "lon": 139.83},
    {"name": "Beijing", "domain": "EAS", "lat": 39.90, "lon": 116.40},
    {"name": "Johannesburg", "domain": "AFR", "lat": -26.20, "lon": 28.04},
    {"name": "Mexico City", "domain": "CAM", "lat": 19.43, "lon": -99.13},
    {"name": "New York", "domain": "NAM", "lat": 40.71, "lon": -74.01},
    {"name": "Sydney", "domain": "AUS", "lat": -33.87, "lon": 151.21},
    {"name": "Riyadh", "domain": "WAS", "lat": 24.71, "lon": 46.67},
    {"name": "Montreal", "domain": "NAM", "lat": 45.50, "lon": -73.57},
    {"name": "Buenos aires", "domain": "SAM", "lat": -34.60, "lon": -58.38},
    {"name": "Los Angeles", "domain": "NAM", "lat": 34.05, "lon": -118.24},
    {"name": "Jakarta", "domain": "SEA", "lat": -6.21, "lon": 106.85}
]

for city in cities:
    create_masks(city["name"], city["domain"], city["lat"], city["lon"])

exit()


create_masks("Tokyo", "EAS", 35.65, 139.83) # Tokyo
create_masks("Beijing", "EAS", 39.90, 116.40) # Beijing
create_masks("Johannesburg", "AFR", -26.20, 28.04) # Johannesburg
create_masks("Mexico City", "CAM", 19.43, -99.13) # Mexico City
create_masks("New York", "NAM", 40.71, -74.01) # New York
create_masks("Sydney", "AUS", -33.87, 151.21) # Sydney
create_masks("Riyadh", "WAS", 24.71, 46.67) # Riyadh
create_masks("Montreal", "NAM", 45.50, -73.57) # Montreal
create_masks("Buenos aires", "SAM", -34.60, -58.38) # Buenos aires
create_masks("Los Angeles", "NAM", 34.05, -118.24) # Los Angeles
create_masks("Jakarta", "SEA", -6.21, 106.85) # Jakarta







