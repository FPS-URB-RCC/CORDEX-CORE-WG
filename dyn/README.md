# FPS URB RCC: Analysis of impact on atmospheric dynamics from representation of urban cities in CORDEX-CORE simulations

L. Fita, UBA/CIMA/IFAECI, C. A. Buenos Aires, Argentina

Here is analysed whether the presence of cities affect atmospheric dynamics of CORDEX-CORE simulations.

It is assume that a city will dynamically affect the atmosphere if in the vertical around its location, potential temperature, wind, humidity present differences respect their surroundings. The absence of such affect, implies that CORDEX-CORE simulations will not be suitable to analyse impacts on for example: landsea breeze, dynamis of convection, air quality and others.

These scripts only work when used with the PyNCplot tools (https://git.cima.fcen.uba.ar/lluis.fita/pyncplot/-/wikis/home). ESGF CORDEX 3-dimensional monthly mean data (ta, hurs, huss, ua, va and orog) has been downloaded into the CIMA's 'papa-deimos' server. These scripts are prepared to use ESGF-CORDEX downloaded data

# urban_3D.bash
This is the main script which controls all the workflow:
 - Retrieval of all the cities to analyze from 'CORDEX-CORE-WG/city_info.csv' repository file
 - For each city:
   - Create a unique file with all the data for each variable 6-grid point centered around the city center
   - Plot figures: 925 hPa horizontal map, SN/WE vertical sections of the direct values and its anomallies

# tamon_pressure_ver.py
Python script to perform the individual city plots for potential temperature
