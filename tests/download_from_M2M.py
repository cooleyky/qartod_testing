# Note: this is not a test of any qartod_testing project functions, just a test that the contents of the gold copy request function work together as expected for one deployment of an example instrument

# Import libraries
import xarray as xr
import os
import re
import io
import stat
import requests
from requests.adapters import HTTPAdapter
from urllib3.util import Retry
from ooinet import M2M
from ooinet.Instrument.common import process_file

# Setup session object for data request
SESSION = requests.Session()
retry = Retry(connect=5, backoff_factor=0.5)
adapter = HTTPAdapter(max_retries=retry)
SESSION.mount('https://', adapter)

# Setup parameters needed to request data
refdes = "GA01SUMO-RII11-02-CTDBPP032"
method = "recovered_inst"
stream = "ctdbp_cdef_instrument_recovered"

# Use the gold copy THREDDs datasets
thredds_url = M2M.get_thredds_url(refdes, method, stream, goldCopy=True)

# Get the THREDDs catalog
thredds_catalog = M2M.get_thredds_catalog(thredds_url)
deployments = M2M.get_deployments(refdes)

# Remove ancillary files from the requested the THREDDs catalog
sensor_files = M2M.clean_catalog(thredds_catalog, stream, deployments) 

# Now build the url to access the data
sensor_files = [re.sub("catalog.html\?dataset=", M2M.URLS["goldCopy_fileServer"], file) for file in sensor_files]

# build path to folder where data will be saved
folder_path = os.path.join(os.path.abspath('../data/external'), method, stream, refdes)
# make folder if it does not already exist
if not os.path.exists(folder_path):
    os.makedirs(folder_path)
os.chmod(folder_path, stat.S_IWRITE)

# Try first file in catalog
file = sensor_files[0]

# Grab file name
file_name = re.findall("deployment.*\.nc$", file)[0]

# Request the data from file server then download
data_request = SESSION.get(file, timeout=(3.05, 120))
ds = xr.load_dataset(io.BytesIO(data_request.content), decode_cf=False)

# Preprocess and write to disk
ds = process_file(ds)
file_path = os.path.join(folder_path, file_name)
# file_path = re.sub('c', 'C', file_path, 1) # I tried swapping c for C in "C:\\Users"
ds.to_netcdf(file_path)