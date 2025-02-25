import xarray as xr
import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import os

# Directories containing the results
directories = {
    "Threshold 10": "results_10",
    "Threshold 20": "results_20",
    "Threshold 30": "results_30",
    "Threshold 40": "results_40",
    "Threshold 50": "results_50",
    "ucdb": "results_ucdb"
}

def count_cells(directory):
    """Counts urban cells from NetCDF files in the specified directory."""
    filelist = glob.glob(f"{directory}/*sftimf*/urmask*.nc", recursive=True)
    data = []

    for file in filelist:
        ds = xr.open_dataset(file, engine='h5netcdf')
        city, model = parse_filename(file)
        nurban = np.sum(ds['urmask'] == 1).item()
        data.append([city, model, nurban])

    return pd.DataFrame(data, columns=["City", "Model", "Urban"])

def parse_filename(filepath):
    """Extracts city and model information from the folder name."""
    folder_name = os.path.basename(os.path.dirname(filepath))
    parts = folder_name.split('_')
    city = parts[0]
    model = parts[1] if len(parts) > 1 else "Unknown"
    return city, model

# Combine data from all directories
result_dfs = []
for label, directory in directories.items():
    df = count_cells(directory)
    df["Threshold"] = label
    result_dfs.append(df)

combined_df = pd.concat(result_dfs)

# Pivot the DataFrame for easier plotting
pivot_df = combined_df.pivot_table(
    index='City', columns=['Model', 'Threshold'], values='Urban', fill_value=0
).reset_index()

# Reorganize the DataFrame for plotting
pivot_df = pivot_df[~pivot_df["City"].str.endswith("11")]

# Prepare a new DataFrame with thresholds and models
new_df = pd.DataFrame()

for city in pivot_df['City'].unique():
    city_data = {'City': city}
    for threshold in directories.keys():
        if threshold == 'ucdb':
            try:
                ucdb_value = pivot_df.loc[pivot_df['City'] == city, ('GERICS', threshold)].values[0]
                city_data['ucdb'] = ucdb_value
            except KeyError:
                city_data['ucdb'] = 0
        else:
            for model in ['ICTP', 'GERICS']:
                try:
                    value = pivot_df.loc[pivot_df['City'] == city, (model, threshold)].values[0]
                    city_data[f"{model}_{threshold}"] = value
                except KeyError:
                    city_data[f"{model}_{threshold}"] = 0
    new_df = pd.concat([new_df, pd.DataFrame([city_data])], ignore_index=True)
    new_df = new_df.drop_duplicates(ignore_index=True)

# Reorganize columns and drop duplicates
column_order = ['City']
for city in new_df['City'].unique():
    for model in ['ICTP', 'GERICS']:
        for threshold in directories.keys():
            if threshold != 'ucdb':
                column_name = f"{model}_{threshold}"
                if column_name not in column_order:
                    column_order.append(column_name)
    if 'ucdb' not in column_order:
        column_order.append('ucdb')

new_df = new_df[column_order]
new_df = new_df.loc[:, ~new_df.columns.duplicated()]

# Save the DataFrame to a CSV file
new_df.to_csv('urban_cells_reorganized.csv', index=False)

# Prepare the data for plotting
plot_data = pd.melt(
    new_df,
    id_vars='City',
    var_name='Category',
    value_name='Urban Cells'
)

# Extract Model and Threshold from the 'Category' column
plot_data[['Model', 'Threshold']] = plot_data['Category'].str.extract(r'(\w+)_Threshold (\d+)')

# Assign colors for ICTP and GERICS
threshold_colors = {
    "10": "#ff9999",
    "20": "#ff4d4d",
    "30": "#ff0000",
    "40": "#b30000",
    "50": "#660000",
}

ictp_colors = {
    "10": "#9999ff",
    "20": "#4d4dff",
    "30": "#0000ff",
    "40": "#0000b3",
    "50": "#000066",
}

ucdb_color = "#000000"  # Black for ucdb

# Plot data
plt.figure(figsize=(14, 6))

# Bar width for grouped bars
bar_width = 0.4
city_count = len(plot_data['City'].unique())

# Plot ICTP data
for threshold, color in ictp_colors.items():
    ictp_data = plot_data[(plot_data['Model'] == 'ICTP') & (plot_data['Threshold'] == threshold)]
    x_positions = np.arange(len(ictp_data)) - 0.2
    plt.bar(
        x_positions,
        ictp_data['Urban Cells'],
        color=color,
        label=f"ICTP {threshold}",
        width=bar_width,
    )

# Plot GERICS data
for threshold, color in threshold_colors.items():
    gerics_data = plot_data[(plot_data['Model'] == 'GERICS') & (plot_data['Threshold'] == threshold)]
    x_positions = np.arange(len(gerics_data)) + 0.2
    plt.bar(
        x_positions,
        gerics_data['Urban Cells'],
        color=color,
        label=f"GERICS {threshold}",
        width=bar_width,
    )

# Plot ucdb data
ucdb_data = plot_data[(plot_data['Category'] == 'ucdb')]
x_positions = np.arange(len(ucdb_data))
plt.bar(
    x_positions,
    ucdb_data['Urban Cells'],
    color=ucdb_color,
    label="Urban cells in polygon",
    width=0.2,
)

# Add labels, title, and legend
plt.xlabel('City')
plt.ylabel('Urban Cells')
plt.title('Urban Cells per Model and Threshold')
plt.xticks(np.arange(city_count), plot_data['City'].unique(), rotation=45, ha='right')
plt.legend(title='Threshold', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save the plot as a PNG file
plt.savefig('plot_ncells_urban_EUR-11.png', bbox_inches='tight')
plt.show()
