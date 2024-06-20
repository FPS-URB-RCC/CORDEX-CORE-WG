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
