import os
import re
import xarray as xr
from datetime import datetime, timedelta
# Function to list files in a directory that match a given pattern
# and provide files with full path.




def list_selected_files(directory_path, pattern):
    """
    Lists files in a directory that match a given pattern and provide files with full path.

    Args:
        directory_path (str): The path to the directory to search.
        pattern (str): The pattern to match against filenames (e.g., "*.txt", "image_*.png").

    Returns:
        list: A list of filenames that match the pattern.

    # Example usage:
    directory = "/path/to/your/directory"  # Replace with your actual directory path
    txt_files = list_selected_files(directory, "*.txt")
    print(f"Text files: {txt_files}")
    
    image_files = list_selected_files(directory, "image_*.png")
    print(f"Image files: {image_files}")

    """
    selected_files = []
    try:
        all_entries = os.listdir(directory_path)
        
        # Convert the wildcard pattern to a regular expression pattern
        # Replace '*' with '.*' to match any sequence of characters
        # Escape other special regex characters if present in the pattern
        regex_pattern = pattern.replace('.', r'\.').replace('*', '.*')
        
        for entry in all_entries:
            full_path = os.path.join(directory_path, entry)
            if os.path.isfile(full_path) and re.match(regex_pattern, entry):
                selected_files.append(entry)

    except FileNotFoundError:
        print(f"Directory not found: {directory_path}")
    except OSError as e:
        print(f"Error accessing directory: {e}")
    #make sure that there is only on file
    if len(selected_files) > 1:
        print(f"Warning: More than one file matches the pattern '{pattern}'. Only the first match will be returned.")
    elif len(selected_files) == 0:
        print(f"No files match the pattern '{pattern}' in the directory '{directory_path}'.")
        return directory_path+'ww3g.nc'

    elif len(selected_files) == 1:
        #print(f"File found: {selected_files[0]}")
        selected_files=directory_path + '/' + selected_files[0]
        return selected_files


# Function to load and concatenate data from multiple files and a given variable based on time range and lon_obs and lat_obs using xarray
def load_and_concatenate_data(experiment, time_range, variables, lon_obs, lat_obs):
    data_list = []
    for single_date in time_range:
        # Construct the file path based on the date and any available ww3g file
        file_path = list_selected_files(f"{experiment}/{single_date.strftime('%Y/%m/%d/%H')}/", f"ww3g_{single_date.strftime('%Y%m%d%H')}-utc_*.nc")
        #file_path = os.listdir(f"{experiment}/{single_date.strftime('%Y/%m/%d/%H')}/ww3g_{single_date.strftime('%Y%m%d%H')}-utc_*.nc") # nzwave+nzlam.nc" #_era5.nc"
        #file_path = !ls f"{experiment}/{single_date.strftime('%Y/%m/%d/%H')}/ww3g_{single_date.strftime('%Y%m%d%H')}-utc_*.nc" # nzwave+nzlam.nc" #_era5.nc"
        # display the file path being processed
        print(f"Processing file: {file_path}")
 
        ## checking that the file exists but and it is larger than 0 bytes
        if not os.path.exists(file_path) or not os.path.isfile(file_path) or not os.path.getsize(file_path) > 0:
        #if os.path.getsize(file_path) == 0:
           print(f"File is empty: {file_path}")
           continue

        else:

           try:
               #ds = xr.open_dataset(file_path)
               ds = xr.open_dataset(file_path, decode_times=True)  # decode_times=True is default
               # selecting  only 24h worth of data
               ds=ds.sel(time = slice(single_date, single_date + timedelta(hours=23.9)))
               ds=ds.sel(lon=lon_obs, lat=lat_obs, method='nearest')
               data_list.append(ds[variables])

           except FileNotFoundError:
               print(f"File not found: {file_path}")
               continue
           ## checking that the file exists but and it is larger than 0 bytes
           except OSError as e:
               print(f"Error reading file {file_path}: {e}")
               continue
        

     # select data based on lon_obs and lat_obs
    #if lon_obs is not None and lat_obs is not None:
    #    data_list = [ds.sel(lon=lon_obs, lat=lat_obs, method='nearest') for ds in data_list]

        # Concatenate the data along the time dimension
    if not data_list:
        print("No data found for the specified time range and location.")
        return None
    
    return xr.concat(data_list, dim='time')
    




