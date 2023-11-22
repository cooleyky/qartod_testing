"""
    @author Kylene Cooley
    @brief Functions used to evaluate statistics on OOI QC flags.
"""
# Import libraries used in this module
import os
import re
import warnings
from tqdm import tqdm
import io
import time
import sys
import ast
import stat

import netrc
import requests
from requests.adapters import HTTPAdapter
import numpy as np
import pandas as pd
import xarray as xr
import dask
from dask.diagnostics import ProgressBar         
from urllib3.util import Retry

from ooinet.M2M import get_thredds_url, get_thredds_catalog, get_deployments, \
    clean_catalog, URLS 
from ooinet.Instrument.common import process_file
from ooi_data_explorations.common import list_files, merge_frames, ENCODINGS

# Initialize Session object for M2M data requests
SESSION = requests.Session()
retry = Retry(connect=5, backoff_factor=0.5)
adapter = HTTPAdapter(max_retries=retry)
SESSION.mount('https://', adapter)


def build_data_path(refdes, method, stream, prefix='', folder='interim',
                    suffix='.nc'):
    """ Build a file name and path to folder from the data stream
    information to return a data path for opening or saving datasets.
    
    Parameters
    ----------
        refdes: string
            built from OOI site, node, and sensor for chosen dataset
        method: string
            'recovered_inst', 'recovered_host', or 'telemetered'(?)
        stream: string
            name of data stream
        prefix: string
            usually 'prod' or 'dev', but empty by default
        folder: string
            'interim' (default), 'processed', 'raw', or 'external'
    
    Returns
    -------
        ds_path: Path object
             absolute path to dataset
    """
     # build filename from dataset type and source
    filename = '-'.join((prefix,refdes,method,stream)) + suffix
    # path to data folder from notebook folder
    data_folder = os.path.relpath('../data')
    # Join folder path and file name to build full relative path
    ds_path=os.path.join(data_folder,folder,filename)
    return ds_path


def ooinet_gold_copy_request(refdes, method, stream, use_dask=False):
    """ Requests gold copy data via M2M, downloads requested datasets,
    and saves files containing a single deployment in the external data
    folder.
    
    Parameters
    ----------
      refdes: string
          Contains the site, node, and sensor ID of interest separated
          by '-'
      method: string
          Represents method of data retrieval, either 'recovered_inst',
          'recovered_host', or 'telemetered'
      stream: string
          Describes data product from given method of data retrieval,
          set by OOI
    
    Returns
    -------
      sensor_files: list
          Names of files that were requested and downloaded
    """
    # Use the gold copy THREDDs datasets
    thredds_url = get_thredds_url(refdes, method, stream, goldCopy=True)

    # Get the THREDDs catalog
    thredds_catalog = get_thredds_catalog(thredds_url)
    deployments = get_deployments(refdes)

    # Remove ancillary files from list of files from THREDDs catalog
    sensor_files = clean_catalog(thredds_catalog, stream, deployments)

    # Now build the url to access the data
    sensor_files = [re.sub("catalog.html\?dataset=",
                           URLS["goldCopy_fileServer"], file) for file in
                    sensor_files if "blank" not in file]

    # build path to folder where data will be saved
    folder_path = os.path.join(os.path.abspath('../data/external'), method,
                               stream, refdes)
    
    # make folder if it does not already exist
    os.makedirs(folder_path, exist_ok=True)

    # preprocess the data and save to disk
    for file in tqdm(sensor_files, 
                     desc='Downloading and Processing Data Files'):
        file_name = re.findall("deployment.*\.nc$", file)[0]
        response = SESSION.get(file, timeout=(3.05, 120))
        if response.ok:
            # load the data file
            if use_dask:
                ds = xr.open_dataset(io.BytesIO(response.content),
                                     decode_cf=False, chunks=10000)
            else:
                ds = xr.load_dataset(io.BytesIO(response.content),
                                     decode_cf=False)
            # Preprocess downloaded data
            ds = process_file(ds)
            ds = ds.drop_vars(['serial_number', 'dcl_controller_timestamp',
                               'date_time_string'], errors='ignore')
            file_path = os.path.join(folder_path, file_name)
            ds.to_netcdf(file_path)
            # ds.to_netcdf(file_path, mode='w', format='NETCDF4',
                         # engine='h5netcdf', encoding=ENCODINGS)
        else:
            print("Bad request: unable to download file %s" % file_name)
    return sensor_files


def dev1_data_request(site, node, sensor, method, stream, params):
    """ Make a request for a netcdf dataset from dev1 with data stream
    identifiers.
    
    Input:
    ------
        refdes: string containing the site, node, and sensor ID of
                interest separated by '-'
        method: string representing method of data retrieval, either
                'recovered_inst', 'recovered_host', or 'telemetered'
        stream: string for resulting data product from given method of
                data retrieval
    
    Returns:
    --------
        data: xarray Dataset containing the concatenated and preprocessed
              data files

    To-do: This function will also save the individual data files in
    the external data folder, organized by site, node, sensor from
    refdes.
    """
    # Initialize credentials for dev1 server
    # This process is borrowed from ooinet.M2M module
    try:
        nrc = netrc.netrc()
        AUTH = nrc.authenticators(
            'ooinet-dev1-west.intra.oceanobservatories.org')
        devlogin, devpassword = AUTH[0], AUTH[2]
        if AUTH is None:
            raise RuntimeError(
                'No entry found for machine \
                ``ooinet-dev1-west.oceanobservatories.org`` in the .netrc \
                file')
    except FileNotFoundError as e:
        raise OSError(e, os.strerror(e.errno), os.path.expanduser('~'))

    # Sub in ooinet-dev1-west.intra.oceanobservatories.org into the avaiable
    # API urls
    Dev01_urls = {}
    for key in URLS:
        url = URLS.get(key)
        if "opendap" in url:
            dev1_url = re.sub("opendap", "opendap-dev1-west.intra", url)
        else:
            dev1_url = re.sub("ooinet", "ooinet-dev1-west.intra", url)
        Dev01_urls[key] = dev1_url
    
    # Use the Dev1 data catalog URL for the request
    api_base_url = Dev01_urls['data']

    # Use the fileServer URL for downloading data files from the
    # THREDDS server
    tds_url = Dev01_urls['fileServer']

    # Create the request URL
    data_request_url ='/'.join((api_base_url, site, node, sensor, method,
                                stream))

    # Set parameters (optional, but not really)
    # We specify a date range to control the size of the dataset requested 
    params = params

    # Build and send the data request
    r = requests.get(data_request_url, params=params, auth=(devlogin, devpassword))
    dev1_request = r.json()

    # Wait until asynchronous request is completed before attempting to
    # download data
    print('Waiting for Dev 01 to process and prepare data request, this may \
            take up to 20 minutes.')
    url = [url for url in dev1_request['allURLs'] if
           re.match(r'.*async_results.*', url)][0]
    check_complete = url + '/status.txt'
    with tqdm(total=400, desc='Waiting', file=sys.stdout) as bar:
        for i in range(400):
            r = SESSION.get(check_complete)
            bar.update()
            bar.refresh()
            if r.status_code == requests.codes.ok:
                bar.n = 400
                bar.last_print_n = 400
                break
            else:
                time.sleep(3)


    # Download the NetCDF data files from the Dev01 thredds server and
    # load into xarray dataset
    url = dev1_request['outputURL']
    files = list_files(url)
    use_dask = False
    frames =[]

    for file in tqdm(files, desc='Downloading and Processing Data Files'):
        file_url = re.sub('catalog.html\?dataset=', tds_url, file)
        r = SESSION.get(file_url, timeout=(3.05, 120))
        if r.ok:
            # load the data file
            if use_dask:
                ds = xr.open_dataset(io.BytesIO(r.content), decode_cf=False,
                                     chunks=10000)
            else:
                ds = xr.load_dataset(io.BytesIO(r.content), decode_cf=False)

            # Perform preprocessing on opened dataset
            # In this workflow we are omitting preprocessing to the
            # *_qartod_executed variables as in
            # ooi_data_explorations.common.process_file()
            # Swap in 'time' for the primary dimension and chunk to 
            # optimize
            ds = ds.swap_dims({'obs': 'time'})
            ds = ds.chunk({'time': 100})
            
            ds = ds.reset_coords()
            keys = ['obs', 'id', 'provenance', 'driver_timestamp',
                    'ingestion_timestamp', 'port_timestamp',
                    'preferred_timestamp']
            for key in keys:
                if key in ds.variables:
                    ds = ds.drop_vars(key)

            # Since the CF decoding of the time is failing, explicitly
            # reset all instances where the units are seconds since
            # 1900-01-01 to the correct CF units and convert the values
            # to datetime64[ns] types
            time_pattern = re.compile(r'^seconds since 1900-01-01.*$')
            ntp_date = np.datetime64('1900-01-01')
            for v in ds.variables:
                if 'units' in ds[v].attrs.keys():
                    if isinstance(ds[v].attrs['units'], str):  
                        # because some units use non-standard characters...
                        if time_pattern.match(ds[v].attrs['units']):
                            del(ds[v].attrs['_FillValue'])  # no fill values for time!
                            ds[v].attrs['units'] = 'seconds since 1900-01-01T00:00:00.000Z'
                            np_time = ntp_date + (ds[v] * 1e9).astype('timedelta64[ns]')
                            ds[v] = np_time

            # Sort by time
            ds = ds.sortby('time') 
            # Clear-up some global attributes we will no longer be using
            keys = ['DODS.strlen', 'DODS.dimName',
                    'DODS_EXTRA.Unlimited_Dimension', '_NCProperties',
                    'feature_Type']
            for key in keys:
                if key in ds.attrs:
                    del(ds.attrs[key])

            try: 
                ds.encoding['unlimited_dims']
                del ds.encoding['unlimited_dims']
            except KeyError:
                pass

            # Resetting cdm_data_type from Point to Station and the
            # featureType from point to timeSeries

            ds.attrs['cdm_data_type'] = 'Station'
            ds.attrs['featureType'] = 'timeSeries'

            # Update some global attributes

            ds.attrs['acknowledgement'] = 'National Science Foundation'
            ds.attrs['comment'] = 'Data collected from the OOI Dev01 M2M API and reworked for use in locally stored NetCDF files.'
            frames.append(ds)
        else:
            failed_file = file.rpartition('/')
            warnings.warn('Failed to download %s' % failed_file[-1]) 

    if frames != []:
        # merge the data frames into a single data set
        data = merge_frames(frames)
    else:
        print('no data frames merged')
        data = None    
    return data


def timeseries_dict_to_xarray(dictionary, ds):
    """ Converts time series data in a dictionary to an Xarray Dataset.
    
    Input:
    ------
      dictionary
      ds
    
    Returns:
    --------
      dict_as_ds
    """

    # convert dict to data frame
    dict_as_df = pd.DataFrame.from_dict(dictionary)

    # Add time vector to df and set time as index
    dict_as_df = dict_as_df.assign(time=ds.time.values)
    dict_as_df = dict_as_df.set_index("time")
    
    # convert df to xarray
    dict_as_ds = dict_as_df.to_xarray() 

    return dict_as_ds


def parse_qartod_executed(ds, parameters):
    """ Parses the qartod tests for the given parameter into separate
    variables.
    
    Parameters:
    -----------
    ds: xarray.DataSet
        The dataset downloaded from OOI with the QARTOD flags applied.
    parameters: list[str]
        The name of the parameters in the dataset to parse the QARTOD flags
        
    Returns:
    --------
    ds: xarray.DataSet
        The dataset with the QARTOD test for the given parameters split
        out into new seperate data variables using the naming
        convention: {parameter}_qartod_{test_name}
    """
    # For the params into a list if only a string
    if type(parameters) is not list:
        parameters = list(parameters)
    
    # Iterate through each parameter
    for param in parameters:
        # Generate the qartod executed name
        qartod_name = f"{param}_qartod_executed"
        
        if qartod_name not in ds.variables:
            continue
    
        # Fix the test types
        ds[qartod_name] = ds[qartod_name].astype(str)
    
        # Get the test order
        test_order = ds[qartod_name].attrs["tests_executed"].split(",")
    
        # Iterate through the available tests and create separate
        # variables with the results
        for test in test_order:
            test_index = test_order.index(test)
            test_name = f"{param}_qartod_{test.strip()}"
            ds[test_name] = ds[qartod_name].str.get(test_index)
    return ds


def get_test_parameters(ds):
    """ Create a dictionary of key-value pairs of dataset 
    variable_name:alternate_parameter_name for parameters that have
    undergone QARTOD testing.
    
    Input:
    ------
    
    Returns:
    --------
    """
    
    test_parameters={}
    for var in ds.variables:
        if "qartod_results" in var:
            # Get the parameter name
            param = var.split("_qartod")[0]

            # Check if the parameter has an alternative ooinet_name
            if "alternate_parameter_name" in ds[param].attrs:
                ooinet_name = ds[param].attrs["alternate_parameter_name"]
            else:
                ooinet_name = param

            # Save the results in a dictionary
            test_parameters.update({
                param: ooinet_name
            })
    # Print out the results
    return test_parameters
    
    
def qartod_results_summary(ds, params, test):
    """ Calculate the statistics for parameter qartod flags. This
    function takes in a list of the parameters and the associated
    QARTOD tests to calculate the number of each flag and the percent
    of total flags that are in each category.
    
    Parameters
    ----------
    ds: xarray.DataSet
        An xarray dataset which contains the data
    params: list[strings]
        A list of the variables/parameters in the given dataset that
        have been tested with QARTOD
    tests: list[strings]
        A list of the QARTOD test names which to parse for the given
        parameters.
        
    Returns
    -------
    results: dict
        A dictionary which contains the number of each QARTOD flag and
        the percent of the total flags for each test applied to each
        parameter in the given dataset.
        
        results = {'parameter':
                        {'test_name':
                            {'total data points': int,
                            'good data points': (int, %),
                            'suspect data points': (int, %),
                            'bad data points': (int, %)}
                            },
                        }
    """
    # Check that the inputs are a list
    if type(params) is not list:
        params = [params]
            
    # Initialize the result dictionary and iterate 
    # through the parameters for each test
    results = {}
    for param in params:
        
        # Now iterate through each test
        test_results = {}
        
            
        # First, check that the test was applied
        test_name = f"{param}_qartod_{test}_test"
        if test_name not in ds.variables:
            continue
            
        # Count the total number of values
        n = ds[test_name].count().compute().values
        
        # First calculate the gross range results
        good = np.where(ds[test_name] == "1")[0]

        # Count the number of suspect/interesting
        suspect = np.where(ds[test_name] == "3")[0]

        # Count the number of fails
        bad = np.where(ds[test_name] == "4'")[0]

        test_results.update({"total": int(n),
                "good": (len(good), np.round(len(good)/n*100, 2)),
                "suspect": (len(suspect), np.round(len(suspect)/n*100, 2)),
                "fail": (len(bad), np.round(len(bad)/n*100, 2))
        })
        
        # Save the test results for each parameter
        results.update({
            param: test_results
        })    
    return results


def qartod_summary_expanded(ds, params, deployment, test):
    """ Calculate the statistics for parameter qartod flags. This
    function takes in a list of the parameters and the associated
    QARTOD tests to calculate the number of each flag and the percent
    of total flags that are in each category.
    
    Parameters
    ----------
    ds: xarray.DataSet
        An xarray dataset which contains the data
    params: list[strings]
        A list of the variables/parameters in the given dataset that
        have been tested with QARTOD
    tests: list[strings]
        A list of the QARTOD test names which to parse for the given
        parameters.
        
    Returns
    -------
    results: dict
        A dictionary which contains the number of each QARTOD flag and
        the percent of the total flags of each parameter for which the
        specified test was applied to the given dataset.
        
        results = {'deployment': int,
                    'test_name total flags': int,
                    'good flags': int,
                    'good %': int,
                    'suspect flags': int,
                    'suspect %': int,
                    'bad flags': int,
                    'bad %': int,...}
    
    Version 26 July 2023, Kylene M Cooley
    """
    # Check that the inputs are a list
    if type(params) is not list:
        params = [params]
            
    # Initialize the result dictionary 
    results = {}
    
    # add key for deployment to results dictionary
    results.update({"deployment" : f"{deployment}" })   
    
    # iterate through the parameters for each test
    for param in params:
            
        # First, check that the test was applied
        test_name = f"{param}_qartod_{test}_test"
        if test_name not in ds.variables:
            results.update({f"{param} total": "NaN",
                    })

            
        else:    
            # Count the total number of values
            n = ds[test_name].count().compute().values

            # First calculate the gross range results
            good = np.where(ds[test_name] == "1")[0]

            # Count the number of suspect/interesting
            suspect = np.where(ds[test_name] == "3")[0]

            # Count the number of fails
            bad = np.where(ds[test_name] == "4'")[0]

            results.update({f"{param} total": int(n),
                    f"{param.split('_')[-1]} good": len(good),
                    f"{param.split('_')[-1]} good %": np.round(len(good)/n*100, 2),
                    f"{param.split('_')[-1]} suspect": len(suspect),
                    f"{param.split('_')[-1]} suspect %": np.round(len(suspect)/n*100, 2),
                    f"{param.split('_')[-1]} fail": len(bad),
                    f"{param.split('_')[-1]} fail %": np.round(len(bad)/n*100, 2)
                    })
            
    return results

def get_deployment_ds(paths_copy):
    """ Loads multi-file dataset for files whose names contain the same
    deployment number as the first file in the list.
    """    
    # get deployment number from first file name
    file = paths_copy[0]    
    deployment = re.findall('deployment00[0-2][0-9]', file)[0][-2:]
    
    # Open each dataset and append to a set of data frames
    deployment_files = [x for x in paths_copy if
                        f'deployment00{deployment}' in x]
    if len(deployment_files)>1:
        deployment_ds = [xr.open_dataset(single_file) for single_file
                         in deployment_files]
        deployment_ds = merge_frames(deployment_ds)
    else:  
        deployment_ds = xr.open_dataset(file)

    [paths_copy.remove(x) for x in deployment_files]
    return deployment_ds, deployment, paths_copy

def collect_statistics(file_paths, test_name):
    """ Calls other functions to calculate statistics from a set of
    files and a name of a QARTOD test. The statistics are organized in
    a DataFrame.
    
    Parameters:
    -----------
        file_paths: list of paths to each file that will have
            statistics calculated. File names must include "deployment00##".
        test_name: string of QARTOD test name, i.e. "gross_range",
            "climatology".
        
    Returns:
    --------
        statistics: Pandas DataFrame containing statistics on each
            parameter with a QARTOD test in order of deployment number,
            then statistics of the full record.
        
    Version 23 Aug 2023, Kylene M Cooley    
    """
    
    # Initialize empty dictionary for statistics
    statistics = {}
    
    # Create a copy of list of file paths for individual deployment
    # statistics
    paths_copy = file_paths.copy()

    while len(paths_copy)>0:
        # Open a dataset with a single deployment
        file_ds, deployment, paths_copy = get_deployment_ds(paths_copy)

        # Get parameters that have QARTOD executed from expected test
        # dataset
        test_parameters = get_test_parameters(file_ds)
        parameters = list(test_parameters.keys())

        # Separate QARTOD test flags in expected test dataset by QARTOD
        # test name
        file_ds = parse_qartod_executed(file_ds, parameters)

        # Update summary statistics dictionary for each deployment,
        # then for all deployments
        print("Evaluating statistics on QARTOD flags for deployment \
            "f"{deployment}")
        summary_results = qartod_summary_expanded(file_ds, parameters,
                                                  deployment, test_name)
        statistics.update({f"{deployment}" : summary_results })

    # Open all data files and create merged full dataset
    merged_ds = [xr.open_dataset(single_file) for single_file in file_paths]
    merged_ds = merge_frames(merged_ds)
    deployment = "all"

    # Summary of flags from merged dataset for full data record
    print("Evaluating statistics on QARTOD flags for all deployments")
    merged_ds = parse_qartod_executed(merged_ds, parameters)
    summary_results = qartod_summary_expanded(merged_ds, parameters,
                                              deployment, test_name)
    statistics.update({ "all" : summary_results })

    # Create data frame from dictionary and check contents
    statistics = pd.DataFrame.from_dict(statistics, orient='index')
    statistics = statistics.set_index('deployment')
    return statistics