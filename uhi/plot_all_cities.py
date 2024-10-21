#!/usr/bin/env python
# coding: utf-8

import glob
import os
import papermill as pm
import sys

from icecream import ic
from utils import RCM_DICT
from utils import YAMLconfig

# Notebook paths
input_notebook = 'urban_area_selection.ipynb'
output_notebook = 'urban_area_selection__papermill.ipynb'

# Climate variable and expected output
variable = 'tasmax'
expected_figure_number = 7

# Load cities configuration from YAML
cities = YAMLconfig('selected_cities.yaml')

# Default thresholds and limits
default_urban_th = 10
default_urban_sur_th = 10
default_lon_lim = 1
default_lat_lim = 1
default_min_city_size = 0

# Iterate over cities and process data
for city in cities:
    abbr_city = city
    long_city = cities[city]['name']
    domain = cities[city]['domain']

    # Generate parameters for the current city
    parameters = {
        'abbr_city': city.split('_')[0],
        'lon_city': cities[city]['lon'],
        'lat_city': cities[city]['lat'],
        'domain': domain,     
        'variable': variable,
        'urban_var': city.split('_')[2],
        'model': city.split('_')[1],
        'urban_th': cities[city].get('urban_th', default_urban_th),
        'urban_sur_th': cities[city].get('urban_sur_th', default_urban_sur_th),
        'lon_lim': cities[city].get('lon_lim', default_lon_lim),
        'lat_lim': cities[city].get('lat_lim', default_lat_lim),
        'min_city_size': cities[city].get('min_city_size', default_min_city_size),
    }

    model = parameters['model']
    model_str = RCM_DICT[domain][model]
    ic(abbr_city, model, domain)

    # Update directory to include urban variable
    directory = f"results/{parameters['urban_var']}_{abbr_city}-{domain}_{model_str}"
    if len(glob.glob(f"{directory}/*.pdf")) == expected_figure_number:
        continue

    # Skip specific domain/model combination if conditions are met
    if domain == 'EUR-22' and model == 'RegCM':
        ic()
        continue
    else:
        try:
            # Execute notebook using Papermill
            pm.execute_notebook(
                input_path=input_notebook,
                output_path=output_notebook,
                parameters=parameters,
                kernel_name='python3'
            )
        except Exception as e:
            # Handle errors by saving a failed version of the output notebook
            output_notebook_failed = output_notebook.replace('.ipynb', f'_ERROR_{abbr_city}-{domain}_{model}.ipynb')
            os.system(f'cp {output_notebook} {output_notebook_failed}')
            print(f'Error executing notebook: {e}')

# Generate EUR-22 counterpart for REMO EUR-11 simulations (RegCM4 not available)
for urban_var in ["sfturf", "sftimf"]:
    for city in cities:
        if cities[city]['domain'] != 'EUR-11':
            continue
        
        abbr_city = city
        long_city = cities[city]['name']
        model = 'REMO'
        domain = 'EUR-22'

        # Generate parameters for the current city
        parameters = {
            'abbr_city': city.split('_')[0],
            'lon_city': cities[city]['lon'],
            'lat_city': cities[city]['lat'],
            'domain': domain,     
            'variable': variable,
            'urban_th': cities[city].get('urban_th', default_urban_th),
            'urban_sur_th': cities[city].get('urban_sur_th', default_urban_sur_th),
            'lon_lim': cities[city].get('lon_lim', default_lon_lim),
            'lat_lim': cities[city].get('lat_lim', default_lat_lim),
            'min_city_size': cities[city].get('min_city_size', default_min_city_size),
        }

        model_str = 'GERICS_REMO2015'
        ic(city, model, cities[city]['domain'])

        # Update directory for the results
        directory = f"results/{urban_var}_{abbr_city}-{domain}_{model_str}"
        if len(glob.glob(f"{directory}/*.pdf")) == expected_figure_number:
            continue
        else:
            try:
                # Execute notebook using Papermill
                pm.execute_notebook(
                    input_path=input_notebook,
                    output_path=output_notebook,
                    parameters=parameters,
                    kernel_name='python3'
                )
            except:
                # Handle errors by saving a failed version of the output notebook
                output_notebook_failed = output_notebook.replace('.ipynb', f'_ERROR_{abbr_city}-{domain}_{model}.ipynb')
                os.system(f'cp {output_notebook} {output_notebook_failed}')
                print('error')

