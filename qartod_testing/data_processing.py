import os
import re
import xarray as xr
import dask
from dask.diagnostics import ProgressBar
from ooinet import M2M
from ooinet.Instrument.common import process_file
import ooi_data_explorations.common as common
from ooi_data_explorations.common import SESSION
import warnings
from tqdm import tqdm
import numpy as np
import netrc
import requests
import io
import time
import sys

def build_data_path(refdes,method,stream,prefix,folder='interim',suffix='.nc'):
    # Input: 
    #   refdes: string built from OOI site, node, and sensor for chosen dataset
    #   method: 'recovered_inst', 'recovered_host', or 'telemetered'(?) 
    #   stream: name of data stream 
    #   source: 'prod' or 'dev'
    #   folder: 'interim' (default), 'processed', 'raw', or 'external'
    #
    # Returns:
    #   ds_path: relative path to dataset from notebook folder
    
    filename = '-'.join((prefix,refdes,method,stream))+suffix              # build filename from dataset type and source

    data_folder = os.path.relpath('../data')                            # path to data folder from notebook folder

    ds_path=os.path.join(data_folder,folder,filename)                   # build full relative path 
    
    return ds_path

def ooinet_gold_copy_request(refdes, method, stream):
    # Input:
    #   refdes: string containing the site, node, and sensor ID of interest separated by '-'
    #   method: string representing method of data retrieval, either 'recovered_inst', 'recovered_host', or 'telemetered'
    #   stream: string for resulting data product from given method of data retrieval
    #
    # Returns:
    #   data: xarray Dataset containing the concatenated and preprocessed data files
    #
    # To-do: This function will also save the individual data files in the external data folder, organized by site, node, sensor from refdes 

    # Generic preprocessing routine to do some generic dataset cleaning/processing
    @dask.delayed
    def preprocess(ds):
        ds = xr.open_dataset(ds, chunks={})
        ds = process_file(ds)
        return ds

    # Use the gold copy THREDDs datasets
    thredds_url = M2M.get_thredds_url(refdes, method, stream, goldCopy=True)

    # Get the THREDDs catalog
    thredds_catalog = M2M.get_thredds_catalog(thredds_url)
    deployments = M2M.get_deployments(refdes)

    # Clean the THREDDs catalog
    # This step separates entries from thredds_catalog if they do not match the stream. These ancillary files are usually provided 
    # because they are used in calculating a derived variable from the measured variable stream.
    sensor_files = M2M.clean_catalog(thredds_catalog, stream, deployments) 

    # Now build the url to access the data
    sensor_files = [re.sub("catalog.html\?dataset=", M2M.URLS["goldCopy_dodsC"], file) for file in sensor_files]

    # preprocess the data
    zs = [preprocess(file) for file in sensor_files]

    # Build path to folder where data files will be indiviually saved
    # Still not sure how to do the file name
    # filepath = build_data_path(refdes, method, stream, folder='external') 
    # for # I messed up when I tried setting this up wit ds in dask.compute(*zs) because then data was not concatenating all of the data as I thought it would
            # ds.to_netcdf(filepath) # come back to this when I get to the step of trying to save data files individually

    # Load all the datasets
    with ProgressBar():
            data = xr.concat([ds.chunk()for ds in dask.compute(*zs)], dim="time")

    return data

def dev1_request(site, node, sensor, method, stream, params):
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
    # This process is borrowed from ooinet.M2M
    try:
        nrc = netrc.netrc()
        AUTH = nrc.authenticators('ooinet-dev1-west.intra.oceanobservatories.org')
        login, password = AUTH[0], AUTH[2]
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
    r = requests.get(data_request_url, params=params, auth=(login, password))
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