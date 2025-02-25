import xarray as xr
import pandas as pd
import glob
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from icecream import ic
from utils import YAMLconfig

# Base directory and climate variable
directory = 'results_190225'
variable = 'tasmax'
urban_var = 'sftimf'  # sfturf/sftimf

cities = YAMLconfig('selected_cities.yaml')


# Function to extract city, domain, and model from the filename
def parse_filename(filename):
    parts = filename.split('/')[-1].split('_')
    city = parts[2].split('-')[0]
    domain = '-'.join(parts[2].split('-')[1:3])  # The domain can have two parts, so we take 1:3
    model = '_'.join(parts[6:7])  # Capture both parts to include 'ECMWF-ERAINT_evaluation' and 'GERICS_REMO2015'
    
    return city, domain, model
# Function to read acycle data from NetCDF files
def read_acycle_data(filelist):
    data = []
    for file in filelist:
        city, domain, model = parse_filename(file)
        ds = xr.open_dataset(file)
        diff = ds['urban_mean'] - ds['rural_mean']
        djf = diff.isel(month=[11, 0, 1]).mean()  # Winter average (DJF)
        jja = diff.isel(month=[6, 7, 8]).mean()  # Summer average (JJA)
        ann = diff.mean()  # Annual average
        data.append([city, domain, model, djf.values.item(), jja.values.item(), ann.values.item()])
    return(pd.DataFrame(data, columns=['City', 'Domain', 'Model', 'DJF', 'JJA', 'Ann']))

# Function to read urmask data (urban and rural mask)
def read_urmask_data(filelist):
    data = []
    for file in filelist:
        city, domain, model = parse_filename(file)
        ds = xr.open_dataset(file)
        nurban = np.sum(ds['urmask'] == 1)  # Count of urban cells
        nrural = np.sum(ds['urmask'] == 0)  # Count of rural cells
        data.append([city, domain, model, nurban.item(), nrural.item()])
    return(pd.DataFrame(data, columns=['City', 'Domain', 'Model', 'n_urban', 'n_rural']))

# Function to swap the levels of columns in a DataFrame
def swap_column_levels(df):
    df.columns = df.columns.swaplevel(0, 1)
    return(df.sort_index(axis=1))

# Function to split the index of the DataFrame into multi-index
def split_index(dframe):
    index_split = dframe.index.str.split('-', expand=True)
    multi_index = pd.MultiIndex.from_tuples(index_split, names=['City', 'Domain', 'Resolution'])
    dframe.index = multi_index
    return dframe

# Function to plot a heatmap from DataFrame
def plot_heatmap(heatmap_data, out='heatmap.pdf', title='', height=16, **kwargs):
    plt.figure(figsize=(10, height))
    ax = sns.heatmap(heatmap_data, annot=True, **kwargs)
    # Cross out missing values
    heatmap_data_mask = heatmap_data.isna()
    for i in range(heatmap_data.shape[0]):
        for j in range(heatmap_data.shape[1]):
            if heatmap_data_mask.iloc[i, j]:
                ax.add_line(plt.Line2D([j, j+1], [i, i+1], color='black', linewidth=2))
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out)

# Load cities information from YAML file
cities = YAMLconfig('selected_cities.yaml')

filelist_acycle = glob.glob(f"{directory}/*/{variable}_*_acycle-ur.nc")
filelist_urmask = glob.glob(f"{directory}/*/urmask_*_fx.nc")

cachefile = f'{directory}/{variable}_uhi_heatmap.csv'

# Create a cache to avoid reading data multiple times
cachefile = f'{directory}/{variable}_{urban_var}_uhi_heatmap.csv'
if os.path.exists(cachefile):
    df = pd.read_csv(cachefile)
else:
    df1 = read_acycle_data(filelist_acycle)
    df2 = read_urmask_data(filelist_urmask)
    df1.set_index(['City', 'Domain', 'Model'], inplace=True)
    df2.set_index(['City', 'Domain', 'Model'], inplace=True)
    df = df1.join(df2)
    df.reset_index(inplace=True)
    df.to_csv(cachefile)

# Calculate the rural-to-urban ratio
df['rur_to_urb_ratio'] = df['n_rural'] / df['n_urban']

cities = YAMLconfig('selected_cities_description.yaml')

is_coastal = {k: 'Coastal' if v['coastal'] else 'InLand' for k, v in cities.items()}
is_mountain = {k: 'Mountain' if v['mountain'] else 'NotMountain' for k, v in cities.items()}

df['is_coastal'] = df['City'].map(is_coastal)
df['is_mountain'] = df['City'].map(is_mountain)

# Create pivot tables for cell counts and UHI (Urban Heat Island) values
num_cells = df.pivot_table(
    index=['City', 'Domain', 'is_coastal', 'is_mountain'], columns='Model',
    values=['n_urban', 'n_rural', 'rur_to_urb_ratio']
)
num_cells = swap_column_levels(num_cells)

heatmap_data = df.pivot_table(index=['City','Domain' ,'is_coastal','is_mountain'], columns='Model', values=['DJF', 'JJA', 'Ann'])
heatmap_data = swap_column_levels(heatmap_data)
#sorted_cities = [city for city in num_cells.index if city in heatmap_data.index]
#heatmap_data = heatmap_data.loc[sorted_cities]
#num_cells = split_index(num_cells)

plot_heatmap(heatmap_data.sort_index(level=['is_coastal','is_mountain','City', 'Domain']),
    out=f'{directory}/{variable}_{urban_var}_uhi_heatmap.pdf', title=f"Seasonal UHI ({variable})", height = 16,
    center=0, cmap='RdBu_r', vmin=-2, vmax=2, fmt=".2f")

plot_heatmap(num_cells,
    out=f'{directory}/{urban_var}_ncells_heatmap.pdf', title="Number of urban / rural cells", height=16,
    cmap='BuPu', fmt=".0f", vmax=25)
