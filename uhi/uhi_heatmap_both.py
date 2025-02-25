import xarray as xr
import pandas as pd
import glob
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from icecream import ic
from utils import YAMLconfig

directory = 'results_190225'
variable = 'tasmax'
urban_var = 'sftimf'  # sfturf/sftimf


def swap_column_levels(df):
    df.columns = df.columns.swaplevel(0, 1)
    return(df.sort_index(axis=1))

def plot_heatmap(heatmap_data, out='heatmap.pdf', title='', height=16, **kwargs):
    plt.figure(figsize=(10, height))
    ax = sns.heatmap(heatmap_data, annot=True, **kwargs)
    # Cross out missing
    heatmap_data_mask = heatmap_data.isna()
    for i in range(heatmap_data.shape[0]):
        for j in range(heatmap_data.shape[1]):
            if heatmap_data_mask.iloc[i, j]:
                ax.add_line(plt.Line2D([j, j+1], [i, i+1], color='black', linewidth=2))
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out)

cities = YAMLconfig('selected_cities.yaml')

filelist_acycle = glob.glob(f"{directory}/*/{variable}_*_acycle-ur.nc")
filelist_obs = glob.glob(f"{directory}/*/{variable}_*acycle-ur-obs.nc")

cachefile = f'{directory}/{variable}_{urban_var}_uhi_heatmap.csv'
cachefile_obs = f'{directory}/{variable}_uhi_heatmap_obs.csv'

df1 = pd.read_csv(cachefile)
df2 = pd.read_csv(cachefile_obs)
df2['Model'] = 'observations'
df = pd.concat([df1, df2], axis=0)
df = df.replace(999999, np.nan)

cities = YAMLconfig('selected_cities_description.yaml')

is_coastal = {k: 'Coastal' if v['coastal'] else 'InLand' for k, v in cities.items()}
is_mountain = {k: 'Mountain' if v['mountain'] else 'NotMountain' for k, v in cities.items()}

df['is_coastal'] = df['City'].map(is_coastal)
df['is_mountain'] = df['City'].map(is_mountain)

heatmap_data = df.pivot_table(index=['City','Domain','is_coastal','is_mountain'], columns='Model', values=['DJF', 'JJA', 'Ann'])
heatmap_data = swap_column_levels(heatmap_data)

plot_heatmap(heatmap_data.sort_index(level=['is_coastal', 'is_mountain', 'City', 'Domain']),
    out=f'{directory}/{variable}_uhi_heatmap_with_obs.pdf', title=f"Seasonal UHI ({variable})", height = 16,
    center=0, cmap='RdBu_r', vmin=-4, vmax=4, fmt=".2f")