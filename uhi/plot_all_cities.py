#!/usr/bin/env python
# coding: utf-8

import glob
import os
import papermill as pm
import sys

from icecream import ic
from utils import cities, RCM_DICT

input_notebook = 'urban_area_selection.ipynb'
output_notebook = 'urban_area_selection__papermill.ipynb'
expected_figure_number = 4

exceptions = [
    "Manila_REMO_SEA-22", # Urban areas have more than 30 of sea
    "Dhaka_REMO_WAS-22", # Urban fractions > 0.1 is not available
    "Tokyo_RegCM_EAS-22", # tasmin not available 
    "Beijing_RegCM_EAS-22", # tasmin not available
    "Chengdu_RegCM_EAS-22", # tasmin not available
    "Seoul_RegCM_EAS-22", # tasmin not available
    "Shanghai_RegCM_EAS-22", # tasmin not available
    "Chengdu_REMO_EAS-22", # Urban fractions > 0.1 is not available
    "Lagos_REMO_AFR-22", # Urban fractions > 0.1 is not available
    "Luanda_REMO_AFR-22", # Urban fractions > 0.1 is not available
    "Luanda_RegCM_AFR-22", # Urban fractions > 0.1 is not available
    "Khartoum_REMO_AFR-22", # Urban fractions > 0.1 is not available
    "Naples_REMO_EUR-22",  # Urban fractions > 0.1 is not available
    "Barcelona_REMO_EUR-11", # Urban fractions > 0.1 is not available
]

for city in cities:
    abbr_city = cities[city]['city']
    domain = cities[city]['domain']
    parameters = {
        'city': city,
        'lon_city': cities[city]['lon'],
        'lat_city': cities[city]['lat'],
        'domain': domain,     
    }
    for model in ['REMO', 'RegCM']:
        model_str = RCM_DICT[domain][model]
        ic(abbr_city, model, domain);
        parameters['model'] = model
        directory = f"results/{abbr_city}-{domain}_{model_str}"
        if len(glob.glob(f"{directory}/*.pdf")) == expected_figure_number:
            continue
        if domain == 'EUR-22' and model == 'RegCM':
            continue
        if f"{abbr_city}_{model}_{domain}" in exceptions:
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
                print('error')
#
# Generate EUR-22 counterpart for REMO EUR-11 simulations.
# EUR-22 runs are not available for RegCM4
#
for city in cities:
    if cities[city]['domain'] != 'EUR-11':
        continue
    abbr_city = cities[city]['city']
    model = 'REMO'
    domain = 'EUR-22'  
    parameters = {
        'city': city,
        'lon_city': cities[city]['lon'],
        'lat_city': cities[city]['lat'],
        'domain': domain,
        'model': model
    }

    model_str = 'GERICS_REMO2015'
    ic(city, model, cities[city]['domain'])
    directory = f"results/{abbr_city}-{domain}_{model_str}"
    if len(glob.glob(f"{directory}/*.pdf")) == expected_figure_number:
        continue
    if f"{abbr_city}_{model}_{domain}" in exceptions:
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
            print('error')
