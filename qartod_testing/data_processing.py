import os
import re
import xarray as xr
import dask
from dask.diagnostics import ProgressBar
from ooinet import M2M
from ooinet.Instrument.common import process_file
import ooi_data_explorations.common as common
from ooi_data_explorations.common import ENCODINGS
import warnings
from tqdm import tqdm
import numpy as np
import netrc
import requests
from requests.adapters import HTTPAdapter
import io
import time
import sys
import pandas as pd
import ast
from ioos_qc.qartod import gross_range_test, climatology_test, ClimatologyConfig
from urllib3.util import Retry
import stat

# Initialize Session object for M2M data requests
SESSION = requests.Session()
retry = Retry(connect=5, backoff_factor=0.5)
adapter = HTTPAdapter(max_retries=retry)
SESSION.mount('https://', adapter)

def build_data_path(refdes,method,stream,prefix='',folder='interim',suffix='.nc'):
    """ Returns data path for opening or saving datasets for a particular sensor.
    
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
    
    filename = '-'.join((prefix,refdes,method,stream))+suffix           # build filename from dataset type and source

    data_folder = os.path.relpath('../data')                            # path to data folder from notebook folder

    ds_path=os.path.join(data_folder,folder,filename)                   # build full relative path 
    
    return ds_path

def ooinet_gold_copy_request(refdes, method, stream, use_dask=False):
    """ Requests gold copy data via M2M, downloads requested datasets, and saves files containing 
    a single deployment in the external data folder.
    
    Parameters
    ----------
      refdes: string 
          Contains the site, node, and sensor ID of interest separated by '-'
      method: string 
          Represents method of data retrieval, either 'recovered_inst', 'recovered_host', or 'telemetered'
      stream: string
          Describes data product from given method of data retrieval, set by OOI
    
    Returns
    -------
      sensor_files: list
          Names of files that were requested and downloaded
    """

    # Use the gold copy THREDDs datasets
    thredds_url = M2M.get_thredds_url(refdes, method, stream, goldCopy=True)

    # Get the THREDDs catalog
    thredds_catalog = M2M.get_thredds_catalog(thredds_url)
    deployments = M2M.get_deployments(refdes)

    # Remove ancillary files from list of files from THREDDs catalog
    sensor_files = M2M.clean_catalog(thredds_catalog, stream, deployments) 

    # Now build the url to access the data
    sensor_files = [re.sub("catalog.html\?dataset=", M2M.URLS["goldCopy_fileServer"], file) for file in sensor_files]

    # build path to folder where data will be saved
    folder_path = os.path.join(os.path.abspath('../data/external'), method, stream, refdes)
    
    # make folder if it does not already exist
    os.makedirs(folder_path, exist_ok=True)

    # preprocess the data and save to disk
    for file in tqdm(sensor_files, desc='Downloading and Processing Data Files'):
        file_name = re.findall("deployment.*\.nc$", file)[0]
        response = SESSION.get(file, timeout=(3.05, 120))
        if response.ok:
            # load the data file
            if use_dask:
                ds = xr.open_dataset(io.BytesIO(response.content), decode_cf=False, chunks=10000)
            else:
                ds = xr.load_dataset(io.BytesIO(response.content), decode_cf=False)
            # Preprocess downloaded data
            ds = process_file(ds)
            if 'serial_number' in ds.variables:
                ds = ds.drop_vars('serial_number')
            file_path = os.path.join(folder_path, file_name)
            ds.to_netcdf(file_path)
        else:
            print("Bad request: unable to download file %s" % file_name)
    return sensor_files

def dev1_data_request(site, node, sensor, method, stream, params):
    # Input:
    #   refdes: string containing the site, node, and sensor ID of interest separated by '-'
    #   method: string representing method of data retrieval, either 'recovered_inst', 'recovered_host', or 'telemetered'
    #   stream: string for resulting data product from given method of data retrieval
    #
    # Returns:
    #   data: xarray Dataset containing the concatenated and preprocessed data files
    #
    # To-do: This function will also save the individual data files in the external data folder, organized by site, node, sensor from refdes 

    # Initialize credentials for dev1 server
    # This process is borrowed from ooinet.M2M module
    try:
        nrc = netrc.netrc()
        AUTH = nrc.authenticators('ooinet-dev1-west.intra.oceanobservatories.org')
        devlogin, devpassword = AUTH[0], AUTH[2]
        if AUTH is None:
            raise RuntimeError(
                'No entry found for machine ``ooinet-dev1-west.oceanobservatories.org`` in the .netrc file')
    except FileNotFoundError as e:
        raise OSError(e, os.strerror(e.errno), os.path.expanduser('~'))

    # Sub in ooinet-dev1-west.intra.oceanobservatories.org into the avaialbe API urls
    Dev01_urls = {}
    for key in M2M.URLS:
        url = M2M.URLS.get(key)
        if "opendap" in url:
            dev1_url = re.sub("opendap", "opendap-dev1-west.intra", url)
        else:
            dev1_url = re.sub("ooinet","ooinet-dev1-west.intra", url)
        Dev01_urls[key] = dev1_url
    
    # Use the Dev1 data catalog URL for the request
    api_base_url = Dev01_urls['data']

    # Use the fileServer URL for downloading data files from the thredds server
    tds_url = Dev01_urls['fileServer']

    # Create the request URL
    data_request_url ='/'.join((api_base_url,site,node,sensor,method,stream))

    # Set parameters (optional, but not really)
    # We specify a date range to control the size of the dataset requested 
    params = params

    # Build and send the data request
    r = requests.get(data_request_url, params=params, auth=(devlogin, devpassword))
    dev1_request = r.json()

    # Wait until asynchronous request is completed before attempting to download data
    print('Waiting for Dev 01 to process and prepare data request, this may take up to 20 minutes.')
    url = [url for url in dev1_request['allURLs'] if re.match(r'.*async_results.*', url)][0]
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


    # Download the NetCDF data files from the Dev01 thredds server and load into xarray dataset
    url = dev1_request['outputURL']
    files = common.list_files(url)
    use_dask = False
    frames =[]

    for file in tqdm(files, desc='Downloading and Processing Data Files'):
        file_url = re.sub('catalog.html\?dataset=', tds_url, file)
        r = SESSION.get(file_url, timeout=(3.05, 120))
        if r.ok:
            # load the data file
            if use_dask:
                ds = xr.open_dataset(io.BytesIO(r.content), decode_cf=False, chunks=10000)
            else:
                ds = xr.load_dataset(io.BytesIO(r.content), decode_cf=False)

            # Perform preprocessing on opened dataset
            # In this workflow we are omitting preprocessing to the *_qartod_executed variables as in ooi_data_explorations.common.process_file()
            
            ds = ds.swap_dims({'obs': 'time'})              # Swap the primary dimension
            ds = ds.chunk({'time': 100})                    # Used for optimization
            
            ds = ds.reset_coords()
            keys = ['obs', 'id', 'provenance', 'driver_timestamp', 'ingestion_timestamp',
                'port_timestamp', 'preferred_timestamp']
            for key in keys:
                if key in ds.variables:
                    ds = ds.drop_vars(key)

            # Since the CF decoding of the time is failing, explicitly reset all instances where the units are
            # seconds since 1900-01-01 to the correct CF units and convert the values to datetime64[ns] types

            time_pattern = re.compile(r'^seconds since 1900-01-01.*$')
            ntp_date = np.datetime64('1900-01-01')
            for v in ds.variables:
                if 'units' in ds[v].attrs.keys():
                    if isinstance(ds[v].attrs['units'], str):  # because some units use non-standard characters...
                        if time_pattern.match(ds[v].attrs['units']):
                            del(ds[v].attrs['_FillValue'])  # no fill values for time!
                            ds[v].attrs['units'] = 'seconds since 1900-01-01T00:00:00.000Z'
                            np_time = ntp_date + (ds[v] * 1e9).astype('timedelta64[ns]')
                            ds[v] = np_time

            # Sort by time
            ds = ds.sortby('time') 
            # Clear-up some global attributes we will no longer be using

            keys = ['DODS.strlen', 'DODS.dimName', 'DODS_EXTRA.Unlimited_Dimension', '_NCProperties', 'feature_Type']
            for key in keys:
                if key in ds.attrs:
                    del(ds.attrs[key])

            try: 
                ds.encoding['unlimited_dims']
                del ds.encoding['unlimited_dims']
            except KeyError:
                pass

            # Resetting cdm_data_type from Point to Station and the featureType from point to timeSeries

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
        data = common.merge_frames(frames)
    else:
        print('no data frames merged')
        data = None
    
    return data

def timeseries_dict_to_xarray(dictionary, ds):
    # Input:
    #   dictionary
    #   ds
    #
    # Returns:
    #   dict_as_ds

    # convert dict to data frame
    dict_as_df = pd.DataFrame.from_dict(dictionary)

    # Add time vector to df and set time as index
    dict_as_df = dict_as_df.assign(time=ds.time.values)
    dict_as_df = dict_as_df.set_index("time")
    
    # convert df to xarray
    dict_as_ds = dict_as_df.to_xarray() 

    return dict_as_ds

GITHUB_BASE_URL = "https://raw.githubusercontent.com/oceanobservatories/qc-lookup/master/qartod"

def load_gross_range_qartod_test_values(refdes, stream, ooinet_param):
    """
    Load the gross range QARTOD test from gitHub
    """
    subsite, node, sensor = refdes.split("-", 2)
    sensor_type = sensor[3:8].lower()
    
    # gitHub url to the gross range table
    GROSS_RANGE_URL = f"{GITHUB_BASE_URL}/{sensor_type}/{sensor_type}_qartod_gross_range_test_values.csv"
    
    # Download the results
    download = requests.get(GROSS_RANGE_URL)
    if download.status_code == 200:
        df = pd.read_csv(io.StringIO(download.content.decode('utf-8')))
        df["parameters"] = df["parameters"].apply(ast.literal_eval)
        df["qcConfig"] = df["qcConfig"].apply(ast.literal_eval)
        
    # Next, filter for the desired parameter
    mask = df["parameters"].apply(lambda x: True if x.get("inp") == ooinet_param else False)
    df = df[mask]
    
    # Now filter for the desired stream
    df = df[(df["subsite"] == subsite) & 
            (df["node"] == node) & 
            (df["sensor"] == sensor) &
            (df["stream"] == stream)]
    
    return df

def load_climatology_qartod_test_values(refdes, param):
    """
    Load the OOI climatology qartod test values table from gitHub
    
    Parameters
    ----------
    refdes: str
        The reference designator for the given sensor
    param: str
        The name of the 
    """
    
    site, node, sensor = refdes.split("-", 2)
    sensor_type = sensor[3:8].lower()
    
    # gitHub url to the climatology tables
    CLIMATOLOGY_URL = f"{GITHUB_BASE_URL}/{sensor_type}/climatology_tables/{refdes}-{param}.csv"
    
    # Download the results
    download = requests.get(CLIMATOLOGY_URL)
    if download.status_code == 200:
        df = pd.read_csv(io.StringIO(download.content.decode('utf-8')), index_col=0)
        df = df.applymap(ast.literal_eval)
    else:
        return None
    return df

def qartod_gross_range_test(refdes, stream, test_parameters, ds):
    # Input:
    #   refdes
    #   stream
    #   test_parameters
    #   ds
    #
    # Returns:
    #   gross_range_results

    # Run through all of the parameters which had the QARTOD tests applied by OOINet and
    # run the tests locally, saving the results in a dictionary
    gross_range_results = {}
    for param in test_parameters:
        # Get the ooinet name
        ooinet_name = test_parameters.get(param)
        
        # Load the gross_range_qartod_test_values from gitHub
        gross_range_qartod_test_values = load_gross_range_qartod_test_values(refdes, stream, ooinet_name)
        
        # Get the qcConfig object, the fail_span, and the suspect_span
        qcConfig = gross_range_qartod_test_values["qcConfig"].values[0]
        fail_span = qcConfig.get("qartod").get("gross_range_test").get("fail_span")
        suspect_span = qcConfig.get("qartod").get("gross_range_test").get("suspect_span")
        
        # Run the gross_range_tenst
        param_results = gross_range_test(
            inp = ds[param].values,
            fail_span = fail_span,
            suspect_span = suspect_span)
        
        # Save the results
        gross_range_results.update(
            {param: param_results}
        )
        
    return gross_range_results

def qartod_climatology_test(refdes, stream, test_parameters, ds):
    # Input:
    #   refdes
    #   stream
    #   test_parameters
    #   ds
    #
    # Returns:
    #   climatology_results

    

    # Run through all of the parameters which had the QARTOD tests applied by OOINet and
    # run the tests locally, saving the results in a dictionary
    climatology_results = {}

    for param in test_parameters:
        # Get the ooinet name
        ooinet_name = test_parameters.get(param)

        # This test uses the same fail span from the gross range test config
        gross_range_qartod_test_values = load_gross_range_qartod_test_values(refdes, stream, ooinet_name)
        qcConfig = gross_range_qartod_test_values["qcConfig"].values[0]
        fail_span = qcConfig.get("qartod").get("gross_range_test").get("fail_span")
        
        # Load the gross_range_qartod_test_values from gitHub
        climatology_qartod_test_values = load_climatology_qartod_test_values(refdes, ooinet_name)
        
        if climatology_qartod_test_values is None:
            climatology_results.update({
                param: "Not implemented."
            })
            continue
        
        # Initialize a climatology config object
        c = ClimatologyConfig()
        
        # Iterate through the pressure ranges
        for p_range in climatology_qartod_test_values.index:
            # Get the pressure range
            pmin, pmax = ast.literal_eval(p_range)

            # Convert the pressure range values into a dictionary
            p_values = climatology_qartod_test_values.loc[p_range].to_dict()

            # Check the pressure values. If [0, 0], then set the range [0, 5000]
            if pmax == 0:
                pmax = 5000
                # if "sea_water_pressure" not in ds.variables:
                #     ds["sea_water_pressure"] = [0, 0]

            for tspan in p_values.keys():
                # Get the time span
                tstart, tend = ast.literal_eval(tspan)

                # Get the values associated with the time span
                vmin, vmax = p_values.get(tspan)

                # Add the test to the climatology config object
                c.add(tspan=[tstart, tend],
                    vspan=[vmin, vmax],
                    fspan=[fail_span[0], fail_span[1]],
                    period="month")
                # print([pmin, pmax])
        # Run the climatology test
        time = ds["time"].to_numpy()
        param_results = climatology_test(c,
                                        inp=ds[param],
                                        tinp=time,
                                        zinp=np.full_like(ds[param], np.nan))
        # param_results = climatology_test(c,
        #                                 inp=ds[param],
        #                                 tinp=time,
        #                                 zinp=ds['sea_water_pressure'])
        
        # Append the results
        climatology_results.update({
            param: param_results
        })


    return climatology_results

def parse_qartod_executed(ds, parameters):
    """
    Parses the qartod tests for the given parameter into separate variables.
    
    Parameters
    ----------
    ds: xarray.DataSet
        The dataset downloaded from OOI with the QARTOD flags applied.
    parameters: list[str]
        The name of the parameters in the dataset to parse the QARTOD flags
        
    Returns
    -------
    ds: xarray.DataSet
        The dataset with the QARTOD test for the given parameters split out
        into new seperate data variables using the naming convention:
        {parameter}_qartod_{test_name}
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
    
        # Iterate through the available tests and create separate variables with the results
        for test in test_order:
            test_index = test_order.index(test)
            test_name = f"{param}_qartod_{test.strip()}"
            ds[test_name] = ds[qartod_name].str.get(test_index)

    return ds

def get_test_parameters(ds):
    # Create a dictionary of key-value pairs of dataset variable name:alternate parameter name for parameters that have undergone QARTOD testing.
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

def run_comparison(ds, param, test_results, test):
    """
    Andrew's example:
    Runs a comparison between the qartod gross range results returned as part of the dataset
    and results calculated locally.
    
    edit:
        added parameter - test: string "gross_range" or "climatology" for test name to compare
    """
    # Get the local test results and convert to string type for comparison
    local_results = test_results[param].astype(str)
    
    # Run comparison
    not_equal = np.where(ds[f"{param}_qartod_{test}_test"] != local_results)[0]
    
    if len(not_equal) == 0:
        return None
    else:
        return not_equal
    
def qartod_results_summary(ds, params, test):
    """
    Calculate the statistics for parameter qartod flags.
    
    This function takes in a list of the parameters and
    the associated QARTOD tests to calculate the number
    of each flag and the percent of the flag.
    
    Parameters
    ----------
    ds: xarray.DataSet
        An xarray dataset which contains the data
    params: list[strings]
        A list of the variables/parameters in the given
        dataset that have been tested with QARTOD
    tests: list[strings]
        A list of the QARTOD test names which to parse
        for the given parameters.
        
    Returns
    -------
    results: dict
        A dictionary which contains the number of each
        QARTOD flag and the percent of the total flags
        for each test applied to each parameter in the
        given dataset.
        
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
            
            }
        )
        
        # Save the test results for each parameter
        results.update({
            param: test_results
        })
    
    return results

def qartod_summary_expanded(ds, params, deployment, test):
    """
    Calculate the statistics for parameter qartod flags.
    
    This function takes in a list of the parameters and
    the associated QARTOD tests to calculate the number
    of each flag and the percent of the flag.
    
    Parameters
    ----------
    ds: xarray.DataSet
        An xarray dataset which contains the data
    params: list[strings]
        A list of the variables/parameters in the given
        dataset that have been tested with QARTOD
    tests: list[strings]
        A list of the QARTOD test names which to parse
        for the given parameters.
        
    Returns
    -------
    results: dict
        A dictionary which contains the number of each
        QARTOD flag and the percent of the total flags
        for each test applied to each parameter in the
        given dataset.
        
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
                    f"{param.split('_')[-1]} good": "NaN",
                    f"{param.split('_')[-1]} suspect": "NaN",
                    f"{param.split('_')[-1]} fail": "NaN"
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
                    f"{param.split('_')[-1]} good": (len(good), np.round(len(good)/n*100, 2)),
                    f"{param.split('_')[-1]} suspect": (len(suspect), np.round(len(suspect)/n*100, 2)),
                    f"{param.split('_')[-1]} fail": (len(bad), np.round(len(bad)/n*100, 2))
                    })
            
    return results

def get_mismatched_flags(expected_ds, local_ds, parameters, deployment, test, expected_file):
    """
    Arguments:
        expected_ds: Xarray Dataset with flags for different QARTOD tests parsed into separate variables
        local_ds: Xarray Dataset containing only the resulting flags from running the QARTOD test locally
        mismatch: dictionary that will hold the results of the comparison
        deployment: 2-character string for the deployment number of the subsite
        test: string of the test to use for comparison, either "gross_range" or "climatology"
    ---------------
    Returns:
        mismatch: the updated dictionary with results of the comparison for the current test and deployment added
        
    Version 12 July 2023, Kylene M Cooley
    """
    # initialize dictionary to hold results of comparison and update with current deployment
    mismatch = {}
    mismatch.update({ "deployment" : f"{deployment}" })
    
    # Loop through parameters while updating dictionaries for mismatched flags
    for param in parameters:

        # First, check that the test was applied
        test_name = f"{param}_qartod_{test}_test"
        if test_name not in expected_ds.variables:
            print(f"{test_name} not implemented")
            mismatch.update({ f"{param}" : np.nan })
        
        else:
            # Evaluate comparison of local test and expected test flags to update dictionary of differences in the results
            print("Checking for mismatched QARTOD flags in "f"{param}")
            flag_mismatch = run_comparison(expected_ds, param, local_ds, test)

            if flag_mismatch is None:
                print("No mismatched values found")
                mismatch.update({ f"{param}" : np.nan })

            else:
                mismatch.update({ f"{param}": {
                        'total' : f"{len(flag_mismatch)} ({round(100*len(flag_mismatch)/len(expected_ds['time']))}%)",
                        'datetimes' : expected_ds['time'][flag_mismatch].values,
                        'expected_flags' : expected_ds[f"{param}_qartod_{test}_test"][flag_mismatch].values,
                        'local_flags' : local_ds[param][flag_mismatch].values,
                        'file_name' : f"{expected_file}"
                }})
    return mismatch