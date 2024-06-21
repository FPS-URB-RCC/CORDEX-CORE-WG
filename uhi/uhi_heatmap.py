import xarray as xr
import pandas as pd
import glob
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from icecream import ic

directory = "results"

filelist_acycle = glob.glob(f"{directory}/*/tasmin_*_acycle-ur.nc")
filelist_urmask = glob.glob(f"{directory}/*/urmask_*_fx.nc")

def parse_filename(filename):
    parts = filename.split('/')[-1].split('_')
    city = parts[1].split('-')[0]
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
        ic(nurban)
        nrural = np.sum(ds['urmask'] == 0)
        data.append([city, model, nurban.item(), nrural.item()])
    return(pd.DataFrame(data, columns=['City', 'Model', 'n_urban', 'n_rural']))

cachefile = 'uhi_heatmap.csv'
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

ic(df)
continentality = [
    "coastal",  # Dhaka
    "coastal",  # Jakarta
    "coastal",  # MumbaiBombay
    "inland",   # NewDelhi
    "inland",   # London
    "inland",   # Cairo
    "inland",   # NewDelhi
    "coastal",  # MumbaiBombay
    "coastal",  # NewYork
    "inland",   # Paris
    "coastal",  # Jakarta
    "mountain", # Tehran
    "coastal",  # Manila
    "coastal",  # NewYork
    "mountain"  # Tehran
]
#df['Continentality'] = continentality
#heatmap_data = df.pivot_table(index=['Continentality', 'City'], columns='Model', values=['DJF', 'JJA'])

heatmap_data = df.pivot_table(index='City', columns='Model', values=['DJF', 'JJA'])

# Create a heatmap
plt.figure(figsize=(4, 14))
ax = sns.heatmap(heatmap_data, annot=True, center=0, cmap='RdBu_r', fmt=".2f")

# Overlay crosses on NaN cells
heatmap_data_mask = heatmap_data.isna()
for i in range(heatmap_data.shape[0]):
    for j in range(heatmap_data.shape[1]):
        if heatmap_data_mask.iloc[i, j]:
            ax.add_line(plt.Line2D([j, j+1], [i, i+1], color='black', linewidth=2))

plt.title("Seasonal UHI (tasmin)")
plt.tight_layout()
plt.show()

