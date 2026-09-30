""" qc_completion.py
A module containing functions for
assessment of automated QC test
completion in comparison to test
parameters documented in the 
ocean-observatories/qc-lookup repo
on GH.

Version: 0.2 (8 Aug 2025)
Previous Versions: 0.1 (31 Oct 2024)

Author: Kylene Cooley (WHOI/OOI-CGSN)
"""
# Import libraries
import io
import requests
import ast
from glob import glob
import numpy as np
import pandas as pd

from ooi_data_explorations.common import m2m_request, m2m_collect
from ooinet import M2M

# Define lists of keywords and function to skip certain data streams 
SKIP_STREAM_KW = ["power", "metadata", "blank", "diagnostic", "dcl_eng",
                  "cpm_eng", "metbk_hourly", "hyd_o", "wavss_a_dcl_fourier",
                  "wavss_a_dcl_motion", "wavss_a_dcl_non_dir", "mopak_o_dcl_rate",
                  "wave_burst", "wfp_eng", "offset", "sio_eng", "glider_eng",
                  "glider_gps", "adcp_config", "imodem_control"]
# Some classes are ignored below while skipping these two classes
# serves a purpose.
SKIP_CLASS_KW = ["FDCHP", "MOPAK"] # , "HYDGN", "DCLENG", "CPMENG"]

def check_skip_kw(value, category):
    """ Check a value against a list of keywords
    that should be skipped for the given
    category. 
    
    Parameters
    -----------
        value: str
            String to search for values in the
            selected list of keywords.
        category: str
            Indicates the keyword list to use
            (e.g., stream, class).
    
    Revision History
    -----------------
        [2025-02-10] K. Cooley, Wrote original function.
    """
    kw_list = category.upper()
    skip_value = [value.find(x) for x in eval(f"SKIP_{kw_list}_KW")]
    return np.mean(skip_value)

# Define functions to load lookup table entries
# Functions from local_qc_test modified to load existing
# tests for multiple parameters
GITHUB_BASE_URL = "https://raw.githubusercontent.com/oceanobservatories/qc-lookup/master/qartod"

def load_gross_range_qartod_test_list(refdes, stream):
    """ Load the list of OOI gross range QARTOD tests developed
    from lookup table on GitHub. 
    
    Parameters
    ----------
    refdes: str
        The reference designator for the sensor to be checked.
    stream: str
        The data stream to be checked.
    
    Revision History
    ----------------
        [YYYY-MM-DD] A. Reed, Original code for qc-lookup table
        import.
        [2024-10-18] K. Cooley, Update table filtering to load
        any test developed for refdes-datastream combination
        and drop unused columns.
    """
    subsite, node, sensor = refdes.split("-", 2)
    sensor_type = sensor[3:8].lower()
    
    # gitHub url to the gross range table
    GROSS_RANGE_URL = f"{GITHUB_BASE_URL}/{sensor_type}/{sensor_type}_qartod_gross_range_test_values.csv"
    
    # Download the results
    download = requests.get(GROSS_RANGE_URL)
    if download.status_code == 200:
        df = pd.read_csv(io.StringIO(download.content.decode('utf-8')))
    else:
        return False
    
    # Now filter for the desired stream
    df = df[(df["subsite"] == subsite) & 
            (df["node"] == node) & 
            (df["sensor"] == sensor) &
            (df["stream"] == stream)]
    if df.empty:
        return False
    else:
        # Next, change parameter field to parameter names
        df["parameters"] = df["parameters"].apply(ast.literal_eval)
        df["parameters"] = df["parameters"].apply(lambda x: x.get("inp"))
        
        # Drop columns for qcConfig, source, notes
        df.drop(columns=["qcConfig", "source", "notes"], inplace=True)
    return df


def load_climatology_qartod_test_list(refdes, stream):
    """ Load the list of OOI climatology QARTOD tests developed
    from lookup table on GitHub. 
    
    Parameters
    ----------
    refdes: str
        The reference designator for the sensor to be checked.
    stream: str
        The data stream to be checked.
        
    Revision History
    ----------------
        [YYYY-MM-DD] A. Reed, Original code for qc-lookup table
        import.
        [2024-10-18] K. Cooley, Update table filtering to load
        any test developed for refdes-datastream combination
        and drop unused columns.
    """
    
    subsite, node, sensor = refdes.split("-", 2)
    sensor_type = sensor[3:8].lower()
    
    # gitHub url to the climatology test tables
    CLIMATOLOGY_URL = f"{GITHUB_BASE_URL}/{sensor_type}/{sensor_type}_qartod_climatology_test_values.csv"

    # Get the correct climatologyTable
    download = requests.get(CLIMATOLOGY_URL)
    
    # Exit function if there is no climatology test table for the instrument class
    if download.status_code == 200:
        df = pd.read_csv(io.StringIO(download.content.decode('utf-8')))
    else:
        return False
    
    # Now filter for the desired stream
    df = df[(df["subsite"] == subsite) & 
            (df["node"] == node) & 
            (df["sensor"] == sensor) &
            (df["stream"] == stream)]
    if df.empty:
        return False
    else:
        # Next, change parameter field to parameter name
        df["parameters"] = df["parameters"].apply(ast.literal_eval)
        df["parameters"] = df["parameters"].apply(lambda x: x.get("inp"))
        
        # Drop columns for climatologyTable, source, notes
        df.drop(columns=["climatologyTable", "source", "notes"], inplace=True)
    return df


# Create a dictionary of key-value pairs of dataset variable
# name:alternate parameter name
def make_test_parameter_dict(data):
    """ Store executed test parameter names and alternate
    ooinet names in a dictionary. Reverse key, value
    order for comparison with developed test list.
    
    Parameters
    ----------
    data: xarray.DataSet
        Dataset containing variables with qartod tests
        executed.
        
    Revision History
    ----------------
        [YYYY-MM-DD] A. Reed, Original code for identifying
        QARTOD parameters.
        [2024-10-18] K. Cooley, Add step for reversing order
        of key, value pairs for comparison with lookup table
        tests.
    """
    test_parameters={}
    for var in data.variables:
        if "qartod_executed" in var:
            # Get the parameter name
            # param = var.split("_qartod")[0]

            # Check if the parameter has an alternative ooinet_name
            if "alternate_parameter_name" in data[var].attrs:
                ooinet_name = data[var].attrs["alternate_parameter_name"]
            else:
                ooinet_name = var

            # Save the results in a dictionary
            test_parameters.update({
                var: ooinet_name
            })
    # Swap order of the results
    test_parameters = dict([(value, key) for key, value in test_parameters.items()])
    return test_parameters


def check_tests_exe(data, test_parameters, grt_table=False, ct_table=False):
    """ Loop through lookup table parameters and compare
    against tests executed in the current dataset. If
    the lookup table parameter is in the list of
    parameters with a test executed, then the list of
    tests executed (currently gross range and/or
    climatology) is saved in a dictionary.
    
    Parameters
    ----------
    data: xarray.Dataset
        Dataset containing variables that may or may not have
        variables ending in "_qartod_executed"
    test_parameters: dict
        Key, value pairs of ooinet test (alternative) names and
        dataset parameter names.
    grt_table: pandas.DataFrame (optional)
        If tests exist, a data frame containing all parameters
        for which a gross range test has been developed.
    ct_table: pandas.DataFrame (optional)
        If tests exist, a data frame containing all parameters
        for which a climatology test has been developed.
    
    Revision History
    ----------------
    [2024-10-18] K. Cooley, Original loop code.
    [2024-10-22] K. Cooley, Revised for use as function.
    [2025-08-08] K. Cooley, Improve error handling for 
        datasets on kdata missing the attributes for 
        QARTOD variables.
    """
    test_exe = {}
    load_m2m = False
    if grt_table is not False:
        for param in grt_table.parameters:
            qartod = param+"_qartod_executed"
            if qartod in test_parameters.keys():
                var = test_parameters[qartod]
                try:
                    test_exe.update({param: data[var].tests_executed})
                except AttributeError:
                    print(f"Dataset is missing {var} QARTOD variable attributes.")
                    test_exe.update({param: "AttributeError"}) # might be unnecessary
                    load_m2m = True
            else:
                test_exe.update({param: "none"})
    if ct_table is not False:
        for param in ct_table.parameters:
            qartod = param+"_qartod_executed"
            if qartod in test_parameters.keys():
                var = test_parameters[qartod]
                try:
                    test_exe.update({param: data[var].tests_executed})
                except AttributeError:
                    print(f"Dataset is missing {var} QARTOD variable attributes.")
                    test_exe.update({param: "AttributeError"})
                    load_m2m = True
            elif param not in grt_table.parameters:
                test_exe.update({param: "none"})
                # print(qartod)
    return test_exe, load_m2m


def make_results_table(grt_table=False, ct_table=False):
    if grt_table is not False:
        table = grt_table.reset_index(drop=True)
        table["GRTtable"] = True
    if ct_table is not False:
        try:
            table["CTtable"] = table["parameters"].isin(list(ct_table["parameters"]))
        finally:
            ct_table["GRTtable"] = False
            ct_table["CTtable"] = True
            ct_table = ct_table[np.bitwise_not(ct_table["parameters"].isin(list(table["parameters"])))]
            table = pd.concat([table, ct_table], ignore_index=True, sort=False)        
    else:
        table["CTtable"] = False
    return table


def add_test_exe(table, test_exe):
    """ Add column to cross-ref results table
    with QARTOD tests executed by parameter.
    """
    table["testsExecuted"] = "none"
    for k in table.index: 
        param = table.at[k, "parameters"]
        table.at[k, "testsExecuted"] = test_exe.get(param)
    return table


def write_results(table, csv_name="test_cross-ref_results.csv", csv_dir="/../data/processed/"):
    """ Write QARTOD test cross-reference results
    table to a CSV. Appends results to existing table if
    found at the resulting csv_path.
    """
    csv_path = csv_dir+csv_name
    if glob(csv_path)==[]:
        file = open(csv_path, mode='w')
        table.to_csv(csv_path, mode='a', index=False)
    else: 
        file = open(csv_path, mode='a')
        table.to_csv(csv_path, mode='a', header=False, index=False)
    # close file 
    file.close()
    print(f"results saved to {csv_path}")
    return


def load_m2m_data(refdes, method, stream, deploy):
    """ Consolidated steps for loading a dataset
    from the OOINet M2M API for the specific use
    case of loading a single deployment. Start
    and end times are included as arguments to
    m2m_result to reduce the length of time spent
    waiting for the request to be completed.
    """
    deploy_info = M2M.get_deployments(refdes, deploy_num=str(deploy))
    dt_start = deploy_info.deployStart[0] + pd.Timedelta(5, 'D')
    dt_end = deploy_info.deployStart[0] + pd.Timedelta(20, 'D')
    dt_start = dt_start.isoformat(timespec='milliseconds')+'Z'
    dt_end = dt_end.isoformat(timespec='milliseconds')+'Z'
    site, node, sensor = refdes.split("-", 2)
    instclass = sensor[3:8]
    m2m_result = m2m_request(site, node, sensor, method, stream,
         start=dt_start, stop=dt_end)
    m2m_data = m2m_collect(m2m_result, tag=(
                            '.*deployment%04d.*%s.*\.nc$' % (deploy, instclass)))
    return m2m_data