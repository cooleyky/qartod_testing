import requests
from requests.adapters import HTTPAdapter
from urllib3.util import Retry

# Setup template for requests to OOINet
SESSION = requests.Session()                    # Session object
retry = Retry(connect=5, backoff_factor=0.5)    # Retry class object
adapter = HTTPAdapter(max_retries=retry)
SESSION.mount('https://', adapter)

format='application/netcdf'

refdes = 'CP01CNSM-MFD37-03-CTDBPD000'
method = 'recovered_inst'
stream = 'ctdbp_cdef_instrument_recovered'