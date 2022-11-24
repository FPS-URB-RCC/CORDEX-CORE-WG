[![MyBinder badge](https://img.shields.io/badge/MyBinder-Jupyter%20Lab-E66581.svg)](https://mybinder.org/v2/gh/jesusff/pyclimenv/main?urlpath=git-pull%3Frepo%3Dhttps%253A%252F%252Fgithub.com%252FFPS-URB-RCC%252FCORDEX-CORE-WG%26urlpath%3Dlab%252Ftree%252FCORDEX-CORE-WG%252F%26branch%3Dmain)

# CORDEX-CORE-WG

This repository is dedicated to the CORDEX-CORE data analysis WG

The WG started under the frame of the CORDEX FPS URB-RCC in 2022. 
The WG aims to investigate available CORDEX-CORE (0.22°) data for large urban areas across the globe. 

## Contents

Files:

| File | Description |
|----------|-------------|
| [selected_cities.yaml](./selected_cities.yaml) | List of selected cities for the analysis |
| [city_info.csv](./city_info.csv) | Information (area, population, density, gp022, coordinates, Köppen climate type, ...) on cities reasonably represented in 0.22 deg. resolution CORDEX-CORE simulations. This is created from GHSL data by [`city_population_area.ipynb`](./city_population_area.ipynb) |

Notebooks:

| Notebook | Description |
|----------|-------------|
| [city_population_area.ipynb](./city_population_area.ipynb) | Displays a map of selected cities by population and area |
| [city_scatter.ipynb](./city_scatter.ipynb) | Displays area vs population scatterplots by Köppen climate |
