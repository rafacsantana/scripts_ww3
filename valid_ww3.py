# This script is a Python adaptation of the MATLAB code provided.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import netCDF4 as nc
from datetime import datetime, timedelta
import xarray as xr
from sys import exit

# Define a function to generate date range (equivalent to MATLAB datenum)
def generate_date_range(start_date, end_date, freq='D'):
    return pd.date_range(start=start_date, end=end_date, freq=freq)

# Define experiment and time range
experiment_name = 'NZWAVE-ERA5'
time_range = generate_date_range('2021-01-01', '2021-01-03')

lat_obs=-43.7567 # from ECAN website -43+(45/60);
lon_obs=173.3358 # 173+(20/60);

# Example function to load NetCDF data
def load_netcdf_data(file_path, variables):
    data = {}
    with nc.Dataset(file_path, mode='r') as ds:
        for var in variables:
            data[var] = ds.variables[var][:].data
    return data

# Example: Define paths (modify as needed)
base_path = '/esi/project/niwa03150/santanarc/hindcast/'

# Example variables to load
variables_to_load = ['lon', 'lat', 'hsig']

# Load sample data
file_example = f"{base_path}{experiment_name}/2021/01/01/00/ww3g_2021010100-utc_nzwave_era5.nc"
data_example = load_netcdf_data(file_example, variables_to_load)

# Extract relevant data
lon_mod = data_example['lon']
lat_mod = data_example['lat']



hs_mod = data_example['hsig'][0, :, :]
#tp = data_example['tpeak'][:, :, 0]
#depth = data_example['depth'][:, :, 0]

# Example: Plot significant wave height
plt.figure(figsize=(10, 8))
plt.pcolormesh(lon_mod, lat_mod, hs_mod, shading='auto', cmap='viridis')
plt.colorbar(label='Significant Wave Height (m)')
plt.title(f'Significant Wave Height on {time_range[0].strftime("%Y-%m-%d")}')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()

exit()

# Example: Time series plotting function
def plot_time_series(dates, values, variable_name, location):
    plt.figure(figsize=(12, 6))
    plt.plot(dates, values, label=f'{variable_name}')
    plt.title(f'{variable_name} Time Series at {location}')
    plt.xlabel('Date')
    plt.ylabel(variable_name)
    plt.grid(True)
    plt.legend()
    plt.show()

# Example use of plot_time_series
# (assuming hs values are extracted for a particular location)
sample_location_hs = hs[100, 100]  # Just an example location
plot_time_series(time_range, np.full(len(time_range), sample_location_hs), 'Hs (m)', 'Sample Location')

