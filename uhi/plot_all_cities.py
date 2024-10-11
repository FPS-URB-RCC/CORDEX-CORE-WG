#!/usr/bin/env python
# coding: utf-8

import glob
import os
import papermill as pm
import sys

from icecream import ic
from utils import RCM_DICT
from utils import YAMLconfig

input_notebook = 'urban_area_selection.ipynb'
output_notebook = 'urban_area_selection__papermill.ipynb'
variable = 'tasmax'
expected_figure_number = 7
cities = YAMLconfig('selected_cities.yaml')

default_urban_th = 0.1
default_urban_sur_th = 0.1
default_lon_lim = 1
default_lat_lim = 1
for urban_var in ["sfturf", "sftimf"]:
    for city in cities:
        abbr_city = city
        long_city = cities[city]['name']
        domain = cities[city]['domain']
        parameters = {
            'abbr_city': city,
            'lon_city': cities[city]['lon'],
            'lat_city': cities[city]['lat'],
            'domain': domain,     
            'variable': variable,
            'urban_th': cities[city].get('urban_th', default_urban_th),
            'urban_sur_th': cities[city].get('urban_sur_th', default_urban_sur_th),
            'lon_lim': cities[city].get('lon_lim', default_lon_lim),
            'lat_lim': cities[city].get('lat_lim', default_lat_lim),
        }
        for model in ['REMO', 'RegCM']:
            model_str = RCM_DICT[domain][model]
            ic(abbr_city, model, domain);
            parameters['model'] = model
            directory = f"results/{urban_var}_{abbr_city}-{domain}_{model_str}"
            if len(glob.glob(f"{directory}/*.pdf")) == expected_figure_number:
                continue
            if domain == 'EUR-22' and model == 'RegCM':
                ic()
                continue
            else:
                try:
                    pm.execute_notebook(
                        input_path=input_notebook,
                        output_path=output_notebook,
                        parameters=parameters,
                        kernel_name='python3'
                    )
                except:
                    output_notebook_failed = output_notebook.replace('.ipynb', f'_ERROR_{abbr_city}-{domain}_{model}.ipynb')
                    os.system(f'cp {output_notebook} {output_notebook_failed}')
                    print('error')
#
# Generate EUR-22 counterpart for REMO EUR-11 simulations.
# EUR-22 runs are not available for RegCM4
#
for urban_var in ["sfturf", "sftimf"]:
    for city in cities:
        if cities[city]['domain'] != 'EUR-11':
            continue
        abbr_city = city
        long_city = cities[city]['name']
        model = 'REMO'
        domain = 'EUR-22'  
        parameters = {
            'abbr_city': abbr_city,
            'lon_city': cities[city]['lon'],
            'lat_city': cities[city]['lat'],
            'domain': domain,
            'variable': variable,     
            'model': model,    
            'urban_th': default_urban_th,
            'urban_sur_th': default_urban_sur_th,
            'lon_lim': cities[city].get('lon_lim', default_lon_lim),
            'lat_lim': cities[city].get('lat_lim', default_lat_lim),
        }
    
        model_str = 'GERICS_REMO2015'
        ic(city, model, cities[city]['domain'])
        directory = f"results/{urban_var}_{abbr_city}-{domain}_{model_str}"
        if len(glob.glob(f"{directory}/*.pdf")) == expected_figure_number:
            continue
        else:
            try:
                pm.execute_notebook(
                    input_path=input_notebook,
                    output_path=output_notebook,
                    parameters=parameters,
                    kernel_name='python3'
                )
            except:
                output_notebook_failed = output_notebook.replace('.ipynb', f'_ERROR_{abbr_city}-{domain}_{model}.ipynb')
                os.system(f'cp {output_notebook} {output_notebook_failed}')
                print('error')
