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
