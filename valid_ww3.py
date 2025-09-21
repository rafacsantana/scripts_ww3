# This script is a Python adaptation of the MATLAB code provided.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import netCDF4 as nc
from datetime import datetime, timedelta
import xarray as xr
import importlib
from sys import exit
import os
import wutils # Importing custom utility functions from santanarc_wave_utils
# Reload the santanarc_wave_utils module to ensure any changes are reflected
importlib.reload(wutils)

plt.close('all')  # Close all existing plots
plt.ion() # Enable interactive mode for plotting
# Set the default font size for plots
#plt.rcParams.update({'font.size': 14})

# computing elapsed time
start_time = datetime.now()
print(f"Script started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
# Define a function to compute elapsed time
def elapsed_time(start_time):
    elapsed = datetime.now() - start_time
    return f"Elapsed time: {elapsed.total_seconds()} seconds"
print("Starting script...")


# Define a function to generate date range (equivalent to MATLAB datenum)
def generate_date_range(start_date, end_date, freq='D'):
    return pd.date_range(start=start_date, end=end_date, freq=freq)

time_range = generate_date_range('2018-01-01', '2025-12-31')
#time_range = generate_date_range('2018-01-01', '2019-1-31')
time_range = generate_date_range('2021-10-20', '2021-10-22')
#time_range = generate_date_range('1983-01-01', '2023-12-31')
#time_range = generate_date_range('1983-01-01', '1983-1-3')
#time_range = generate_date_range('2014-12-30', '2030-1-2')
#time_range = generate_date_range('2005-1-1', '2100-1-2')

# Example: Define paths (modify as needed)
base_path = '/esi/project/niwa03150/santanarc/hindcast/'
fig_path = '/esi/project/niwa03150/santanarc/hindcast/figures/'

# listing the experiments available in base_path
def list_experiments(base_path):
    # List all directories in the base_path
    experiments = [d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d))]
    return experiments
# List available experiments
experiments = list_experiments(base_path)

# keep experiments that have 'WAVE' in the name
experiments = [exp for exp in experiments if 'WAVE-HR' in exp]

print("Available experiments:")
for exp in experiments:
    print(f"- {exp}")
# If you want to use a specific experiment, uncomment the line below and set the experiment name
# experiment_name = 'NZWAVE-ERA5'  # Example experiment name

# Define experiment and time range
#experiment_name = 'NZWAVE-ERA5'
experiments = ['NZWAVE-ERA5']
#experiments = ['NZWAVE-HR']

# loop through experiments and run the analysis on them
for experiment_name in experiments:
    print(f"Processing experiment: {experiment_name}")
    
    
    station_name = 'Banks_Peninsula'  # Example station name
    lat_obs=-43.7567 # from ECAN website -43+(45/60);
    lon_obs=173.3358 # 173+(20/60);
    
    
    # Example function to load NetCDF longitude and latitude data
    #def load_netcdf_data(file_path, variables):
    #    data = {}
    #    with nc.Dataset(file_path, mode='r') as ds:
    #        for var in variables:
    #            data[var] = ds.variables[var][:].data
    #    return data
    #
    #variables_to_load = ['lon', 'lat'] # , 'hsig']
    
    # Load sample data
    #file_example = f"{base_path}{experiment_name}/2021/01/01/00/ww3g_2021010100-utc_nzwave_nzlam.nc" # era5.nc"
    #data_example = load_netcdf_data(file_example, variables_to_load)
    
    # Extract relevant data
    #lon_mod = data_example['lon']
    #lat_mod = data_example['lat']
    
    
    
    # Load and concatenate data for significant wave height (hsig) over the specified time range
    hsig_variable = 'hsig'  # Example variable to load
    data_hsig = wutils.load_and_concatenate_data(base_path + experiment_name, time_range, hsig_variable, lon_obs, lat_obs)
    # Check if data was loaded successfully
    if data_hsig is None:
        print("No data loaded for the specified variable and time range.")
        exit()
    
    # Plot a figure of significant wave height data with time as the x-axis
    hs = data_hsig.values  # Convert to numpy array
    # Ensure the data is in the correct shape (time, lat, lon)
    if len(hs.shape) == 3:
        hs = hs.reshape(hs.shape[0], -1)  # Flatten spatial dimensions if needed
    # If data is 2D, reshape it to (time, lat, lon)
    if len(hs.shape) == 2:
        hs = hs.reshape(hs.shape[0], 1, -1)  # Add a singleton dimension for latitude
    # Ensure the data is in the correct shape (time, lat, lon)
    if len(hs.shape) == 1:
        hs = hs.reshape(-1, 1, 1)  # Add singleton dimensions for latitude and longitude
    # Check if the data has the expected dimensions
    if len(hs.shape) != 3:
        print(f"Unexpected data shape: {hs.shape}. Expected shape is (time, lat, lon).")
        exit()
    # Ensure the data is in the correct shape (time, lat, lon)
    if len(hs.shape) == 2:
        hs = hs.reshape(hs.shape[0], 1, -1)  # Add a singleton dimension for latitude
    # If data is 1D, reshape it to (time, lat, lon)
    if len(hs.shape) == 1:
        hs = hs.reshape(-1, 1, 1)  # Add singleton dimensions for latitude and longitude
    # Ensure the data is in the correct shape (time, lat, lon)
    if len(hs.shape) != 3:
        print(f"Unexpected data shape: {hs.shape}. Expected shape is (time, lat, lon).")
        exit()
    
    
    fig = plt.figure(figsize=(16, 8))    
    # Plot significant wave height data using just dots without a line
    plt.plot(data_hsig.time.values, hs[:, 0, 0], label='Significant Wave Height (m)', color='blue', marker='.', linestyle='None', markersize=5)
    
    sta_name= station_name.replace('_', ' ').replace('-', ' ')
    plt.title(f'Significant Wave Height at {sta_name}')
    plt.xlabel('Time')
    plt.ylabel('Significant Wave Height (m)')
    plt.legend()
    #plt.grid()
    
    
    plt.show()
    
    # Saving figure. Create a repo in fig_path if it doesn't exist with the experiment name
    fig_dir = f"{fig_path}{experiment_name}/"
    import os
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    # Save the figure with a descriptive name and include time range
    plt.savefig(f"{fig_dir}{station_name}_significant_wave_height_{time_range[0].strftime('%Y%m%d')}_{time_range[-1].strftime('%Y%m%d')}.png", bbox_inches='tight')



# Print elapsed time
print(elapsed_time(start_time))
# Exit the script
print("Script completed successfully.")


exit()


# Extract longitude and latitude from the dataset
lon_mod = data_hsig.lon.values
lat_mod = data_hsig.lat.values
# Convert time to pandas datetime for plotting
time_mod = pd.to_datetime(data_hsig.time.values)
# Print loaded data for debugging
print(f"Loaded data shape: {hs.shape}")
# Print the first few values for verification
print(f"First few significant wave height values: {hs.flatten()[:5]}")
# Print time range for debugging
print(f"Time range: {time_mod[0]} to {time_mod[-1]}")
# Print longitude and latitude for debugging
print(f"Longitude range: {lon_mod.min()} to {lon_mod.max()}")
# Print latitude range for debugging
print(f"Latitude range: {lat_mod.min()} to {lat_mod.max()}")
# Exit early if no data is loaded
if data_hsig is None:
    print("No significant wave height data loaded.")
    exit()
# Print the loaded data for debugging
print("Data loaded successfully.")
# Print the shape of the loaded data
print(f"Shape of loaded data: {data_hsig.shape}")
# Print the first few values of the loaded data
print(f"First few values of significant wave height: {data_hsig.values.flatten()[:5]}")
# Print the time range of the loaded data
print(f"Time range of loaded data: {time_mod[0]} to {time_mod[-1]}")
# If you want to load additional variables, you can modify the variables_to_load list
# and call the load_and_concatenate_data function again.



# Example variables to load

exit()

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

