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


## Including the data

Since the datasets are too big to be stored on GitHub, they are currently stored on RUB's geocloud, see [here](https://geo-cloud.geographie.ruhr-uni-bochum.de/index.php/s/jfNy2xkPktKqzR6). This folder is password-protected. In case you are interested to have access, please reach out to [Matthias Demuzere via email](mailto:matthias.demuzere@rub.de?subject=[GitHub]%CORDEX%CORE%DataRequest).


1. Either download the data or setup a Nextcloud-sync
1. While in the repo run the following command to create a symbolic link to have the data folder
   available as `./data` in the repositories root. It works best if you use absolute paths.
   ```bash
   ln --symbolic /home/<user>/<some-other-folders>/CORDEX-CORE-WG /home/<user>/<some-other-folders>/CORDEX-CORE-WG/data
   ```
1. You now should have a symbolic link for `./data`, e.g.:
   ```console
   20:58:21-demuzmp4@rub:CORDEX-CORE-WG$ ls -l
   total 616
   -rw-r--r-- 1 demuzmp4 Domain Users  23688 Jan 17 20:38 city_info.csv
   -rw-r--r-- 1 demuzmp4 Domain Users 473558 Jan 17 20:38 city_population_area.ipynb
   -rw-r--r-- 1 demuzmp4 Domain Users 116289 Jan 17 20:38 city_scatter.ipynb
   lrwxrwxrwx 1 demuzmp4 Domain Users     44 Jan 17 20:57 data -> /home/demuzmp4/Nextcloud/data/CORDEX-CORE-WG
   drwxr-xr-x 2 demuzmp4 Domain Users   4096 Jan 17 20:38 indices
   -rw-r--r-- 1 demuzmp4 Domain Users   2666 Jan 17 20:59 README.md
   -rw-r--r-- 1 demuzmp4 Domain Users   1710 Jan 17 20:38 selected_cities.yaml
   ```