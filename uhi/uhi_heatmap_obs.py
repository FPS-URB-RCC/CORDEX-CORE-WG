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

def parse_filename(filename):
    parts = filename.split('/')[-1].split('_')
    city = parts[2].split('-')[0]
    domain = '-'.join(parts[2].split('-')[1:])
    return city, domain

def read_acycle_data(filelist):
    data = []
    cities=[]
    for file in filelist:
        city, domain = parse_filename(file)
        if city not in cities:
            cities.append(city)
            ds = xr.open_dataset(file)
            if 'urban_mean_anom' not in ds.variables:
                print(f"Warning: no observations on avaliable on {city}")
                djf = 999999
                jja = 999999
                ann = 999999
                data.append([city, domain,djf, jja, ann])
            else:
                djf = ds['urban_mean_anom'].sel(dim_0=[12, 1, 2]).mean(dim='dim_0')
                jja = ds['urban_mean_anom'].sel(dim_0=[6, 7, 8]).mean(dim='dim_0')
                ann = ds['urban_mean_anom'].mean(dim='dim_0')
                data.append([city, domain, djf.values.item(), jja.values.item(), ann.values.item()])
    return(pd.DataFrame(data, columns=['City','Domain', 'DJF', 'JJA', 'Ann']))

def split_index(dframe):
    index_split = dframe.index.str.split('-', expand=True)
    multi_index = pd.MultiIndex.from_tuples(index_split, names=['City', 'Domain', 'Resolution'])
    dframe.index = multi_index
    return(dframe)

def plot_heatmap(heatmap_data, out='heatmap.pdf', title='', height=16, **kwargs):
    plt.figure(figsize=(6, height))
    ax = sns.heatmap(heatmap_data, annot=True, **kwargs)
    ax.set(xlabel="", ylabel="")
    ax.xaxis.tick_top()
    # Cross out missing
    heatmap_data_mask = heatmap_data.isna()
    for i in range(heatmap_data.shape[0]):
        for j in range(heatmap_data.shape[1]):
            if heatmap_data_mask.iloc[i, j]:
                ax.add_line(plt.Line2D([j, j+1], [i, i+1], color='black', linewidth=2))
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out)

cities = YAMLconfig('selected_cities.yaml')
filelist_acycle = glob.glob(f"{directory}/*/{variable}_*acycle-ur-obs.nc")

cachefile = f'{directory}/{variable}_uhi_heatmap_obs.csv'
if os.path.exists(cachefile):
    df = pd.read_csv(cachefile)
else:
    df = read_acycle_data(filelist_acycle)
    df.set_index(['City', 'Domain'], inplace=True)
    df.reset_index(inplace=True)
    df.to_csv(cachefile, index=False)

cities = YAMLconfig('selected_cities_description.yaml')

is_coastal = {k: 'Coastal' if v['coastal'] else 'InLand' for k, v in cities.items()}
is_mountain = {k: 'Mountain' if v['mountain'] else 'NotMountain' for k, v in cities.items()}

df['is_coastal'] = df['City'].map(is_coastal)
df['is_mountain'] = df['City'].map(is_mountain)

heatmap_data = df.pivot_table(index=['City', 'is_coastal', 'is_mountain','Domain'], values=['DJF', 'JJA', 'Ann'])
heatmap_data = heatmap_data.replace(999999, np.nan)

plot_heatmap(heatmap_data.sort_index(level=['is_coastal', 'is_mountain', 'City', 'Domain']),
    out=f'{directory}/{variable}_uhi_heatmap_obs.pdf', title=f"Seasonal UHI observations", height = 16,
    center=0, cmap='RdBu_r', vmin=-4, vmax=4, fmt=".2f")