import xarray as xr
import pandas as pd
import glob
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from icecream import ic

directory = "results"

def parse_filename(filename):
    parts = filename.split('/')[-1].split('_')
    city = parts[1]#.split('-')[0]
    model = parts[6][:6]
    return city, model

def read_acycle_data(filelist):
    data = []
    for file in filelist:
        city, model = parse_filename(file)
        ds = xr.open_dataset(file)
        diff = ds['urban_mean'] - ds['rural_mean']
        djf = diff.sel(month=[12, 1, 2]).mean(dim='month')
        jja = diff.sel(month=[6, 7, 8]).mean(dim='month')
        ann = diff.mean(dim='month')
        data.append([city, model, djf.values.item(), jja.values.item(), ann.values.item()])
    return(pd.DataFrame(data, columns=['City', 'Model', 'DJF', 'JJA', 'Ann']))

def read_urmask_data(filelist):
    data = []
    for file in filelist:
        city, model = parse_filename(file)
        ds = xr.open_dataset(file)
        nurban = np.sum(ds['urmask'] == 1)
        nrural = np.sum(ds['urmask'] == 0)
        data.append([city, model, nurban.item(), nrural.item()])
    return(pd.DataFrame(data, columns=['City', 'Model', 'n_urban', 'n_rural']))

def swap_column_levels(df):
    df.columns = df.columns.swaplevel(0, 1)
    return(df.sort_index(axis=1))

def split_index(dframe):
    index_split = dframe.index.str.split('-', expand=True)
    multi_index = pd.MultiIndex.from_tuples(index_split, names=['City', 'Domain', 'Resolution'])
    dframe.index = multi_index
    return(dframe)

def plot_heatmap(heatmap_data, out='heatmap.pdf', title='', height=16, **kwargs):
    plt.figure(figsize=(6, height))
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

filelist_acycle = glob.glob(f"{directory}/*/tasmin_*_acycle-ur.nc")
filelist_urmask = glob.glob(f"{directory}/*/urmask_*_fx.nc")

cachefile = f'{directory}/uhi_heatmap.csv'
if os.path.exists(cachefile):
    df = pd.read_csv(cachefile)
else:
    df1 = read_acycle_data(filelist_acycle)
    df2 = read_urmask_data(filelist_urmask)
    df1.set_index(['City', 'Model'], inplace=True)
    df2.set_index(['City', 'Model'], inplace=True)
    df = df1.join(df2)
    df.reset_index(inplace=True)
    df.to_csv(cachefile)

df['rur_to_urb_ratio'] = df['n_rural'] / df['n_urban']
num_cells = df.pivot_table(
    index='City', columns='Model',
    values=['n_urban', 'n_rural', 'rur_to_urb_ratio']
)
num_cells = swap_column_levels(num_cells)
num_cells = num_cells.sort_values(by=('REMO20', 'n_urban'), ascending=False)
heatmap_data = df.pivot_table(index='City', columns='Model', values=['DJF', 'JJA', 'Ann'])
sorted_cities = [city for city in num_cells.index if city in heatmap_data.index]
heatmap_data = heatmap_data.loc[sorted_cities]
heatmap_data = swap_column_levels(heatmap_data)
heatmap_data = split_index(heatmap_data)
num_cells = split_index(num_cells)

plot_heatmap(heatmap_data.xs('11', level='Resolution'),
    out=f'{directory}/uhi_heatmap_EUR-11.pdf', title="Seasonal UHI (tasmin)", height = 6,
    center=0, cmap='RdBu_r', vmin=-2, vmax=2, fmt=".2f")
plot_heatmap(heatmap_data.xs('22', level='Resolution'),
    out=f'{directory}/uhi_heatmap_CDX-22.pdf', title="Seasonal UHI (tasmin)", height = 16,
    center=0, cmap='RdBu_r', vmin=-2, vmax=2, fmt=".2f")
plot_heatmap(num_cells.xs('11', level='Resolution'),
    out=f'{directory}/ncells_heatmap_EUR-11.pdf', title="Number of urban / rural cells", height = 6,
    cmap='BuPu', fmt=".0f", vmax=25)
plot_heatmap(num_cells.xs('22', level='Resolution'),
    out=f'{directory}/ncells_heatmap_CDX-22.pdf', title="Number of urban / rural cells", height = 16,
    cmap='BuPu', fmt=".0f", vmax=25)
