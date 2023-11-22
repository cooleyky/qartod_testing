"""
    @author Kylene Cooley
    @brief Functions for running a QC test locally and compare wtih
            expected QC test results. Borrowed functions are attributed
            to the original author in the docstring.
"""

# Import modules used in this notebook
import io
import ast
import os

import numpy as np
import xarray as xr
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from ioos_qc.qartod import gross_range_test, climatology_test, \
    ClimatologyConfig
from urllib3.util import Retry

from qartod_testing.qc_flag_statistics import get_test_parameters, \
    timeseries_dict_to_xarray, get_deployment_ds

# Initialize Session object for M2M data requests
SESSION = requests.Session()
retry = Retry(connect=5, backoff_factor=0.5)
adapter = HTTPAdapter(max_retries=retry)
SESSION.mount('https://', adapter)

# Set base URL for QARTOD test lookup tables
GITHUB_BASE_URL = "https://raw.githubusercontent.com/oceanobservatories/qc-lookup/master/qartod"


def load_gross_range_qartod_test_values(refdes, stream, ooinet_param):
    """Load the gross range QARTOD lookup table from GitHub repo
    oceanobservatories/qc-lookup and return the table as a Pandas
    DataFrame.
    Original function created by Andrew Reed.
    
    Input:
    ------
    
    Returns:
    --------
    
    """
    subsite, node, sensor = refdes.split("-", 2)
    sensor_type = sensor[3:8].lower()
    
    # GitHub url to the gross range table
    GROSS_RANGE_URL = f"{GITHUB_BASE_URL}/{sensor_type}/{sensor_type}_qartod_gross_range_test_values.csv"
    
    # Download the results
    download = requests.get(GROSS_RANGE_URL)
    if download.status_code == 200:
        df = pd.read_csv(io.StringIO(download.content.decode('utf-8')))
        df["parameters"] = df["parameters"].apply(ast.literal_eval)
        df["qcConfig"] = df["qcConfig"].apply(ast.literal_eval)
        
    # Next, filter for the desired parameter
    mask = df["parameters"].apply(lambda x: True if x.get("inp") == 
                                  ooinet_param else False)
    df = df[mask]
    
    # Now filter for the desired stream
    df = df[(df["subsite"] == subsite) & 
            (df["node"] == node) & 
            (df["sensor"] == sensor) &
            (df["stream"] == stream)]   
    return df


def load_climatology_qartod_test_values(refdes, stream, param):
    """Load the OOI climatology QARTOD test values table from GitHub
    Original function created by Andrew Reed.
    
    Parameters
    ----------
    refdes: str
        The reference designator for the given sensor
    param: str
        The name of the 
    """
    
    subsite, node, sensor = refdes.split("-", 2)
    sensor_type = sensor[3:8].lower()
    
    # GitHub url to the climatology test tables
    CLIMATOLOGY_URL = f"{GITHUB_BASE_URL}/{sensor_type}/{sensor_type}_qartod_climatology_test_values.csv"
    
    # Get the correct climatologyTable
    download = requests.get(CLIMATOLOGY_URL)
    df = pd.read_csv(io.StringIO(download.content.decode('utf-8')))
    df["parameters"] = df["parameters"].apply(ast.literal_eval)
    # Next, filter for the desired parameter
    mask = df["parameters"].apply(lambda x: True if x.get("inp") == param else False)
    df = df[mask]
    
    # Now filter for the desired stream
    df = df[(df["subsite"] == subsite) & 
            (df["node"] == node) & 
            (df["sensor"] == sensor) &
            (df["stream"] == stream)]
    
    # Get the "zinp" as a check
    zinp = df["parameters"].values[0].get('zinp')
    
    # Get the correct climatologyTable
    climatologyTable = df["climatologyTable"].values[0]

    # Construct the url to the climatologyTable
    CLIMATOLOGY_TABLE_URL = f"{GITHUB_BASE_URL}/{sensor_type}/{climatologyTable}"

    # Download the results
    download = requests.get(CLIMATOLOGY_TABLE_URL)
    if download.status_code == 200:
        df = pd.read_csv(io.StringIO(download.content.decode('utf-8')), index_col=0)
        df = df.applymap(ast.literal_eval)
    else:
        return None
    return df, zinp


def run_qartod_gross_range(refdes, stream, test_parameters, ds):
    """ Run through all of the parameters which had the QARTOD tests
    applied by OOINet and run the tests locally, saving the results in
    a dictionary.
    Based on workflow developed by Andrew Reed.
    
    Input:
    ------
      refdes
      stream
      test_parameters
      ds
    
    Returns:
    --------
      gross_range_results
    """

    gross_range_results = {}
    for param in test_parameters:
        # Get the ooinet name
        ooinet_name = test_parameters.get(param)
        
        # Load the gross_range_qartod_test_values from gitHub
        gross_range_qartod_test_values = load_gross_range_qartod_test_values(
            refdes, stream, ooinet_name)
        
        # Get the qcConfig object, the fail_span, and the suspect_span
        qcConfig = gross_range_qartod_test_values["qcConfig"].values[0]
        fail_span = qcConfig.get("qartod").get("gross_range_test").get(
            "fail_span")
        suspect_span = qcConfig.get("qartod").get("gross_range_test").get(
            "suspect_span")
        
        # Run the gross_range_test
        param_results = gross_range_test(
            inp = ds[param].fillna(999999).values,
            fail_span = fail_span,
            suspect_span = suspect_span)
        
        # Save the results
        gross_range_results.update(
            {param: param_results}
        )
    return gross_range_results


def run_qartod_climatology(refdes, stream, test_parameters, ds):
    """ Run through all of the parameters which had the QARTOD tests
    applied by OOINet and run the QARTOD climatology tests locally, saving the results in
    a dictionary.
    Based on workflow developed by Andrew Reed.
    
    Input:
    ------
      refdes
      stream
      test_parameters
      ds
    
    Returns:
    --------
      climatology_results
    """
    # Allocate an empty dictionary to store results
    climatology_results = {}

    for param in test_parameters:
        # Get the ooinet name
        ooinet_name = test_parameters.get(param)
        
        # Load the climatology QARTOD test values from GitHub
            # try:
        climatology_qartod_test_values, zinp = \
            load_climatology_qartod_test_values(refdes, stream, ooinet_name)
            # except:
            #     climatology_results.update({
            #         param: "Not implemented."
            #     })
            #     continue
        
        if climatology_qartod_test_values is None:
            climatology_results.update({
                param: "Not implemented."
            })
            continue
        
        # Check that the 'zinp' is in the dataset. If not, need to add a dummy variable
        if zinp not in ds.variables:
            ds['zinp'] = (["time"], np.ones(ds["time"].shape))
            zinp = 'zinp'
        
        # This test uses the same fail span from the gross range test config
        gross_range_qartod_test_values = load_gross_range_qartod_test_values(
            refdes, stream, ooinet_name)
        qcConfig = gross_range_qartod_test_values["qcConfig"].values[0]
        fail_span = qcConfig.get("qartod").get("gross_range_test").get(
            "fail_span")
        
        # Initialize a climatology config object
        c = ClimatologyConfig()
        
        # Iterate through the pressure ranges
        for p_range in climatology_qartod_test_values.index:
            # Get the pressure range
            pmin, pmax = ast.literal_eval(p_range)

            # Convert the pressure range values into a dictionary
            p_values = climatology_qartod_test_values.loc[p_range].to_dict()

            # Check the pressure values. If current range is [0, 0],
            # then set the range to [0, 5000]
            if pmax == 0:
                pmax = 5000
               
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
                               
        # Run the climatology test
        param_results = climatology_test(c,
                                         inp=ds[param].fillna(999999),
                                         tinp=ds["time"].values,
                                         zinp=ds[zinp])
        
        # Append the results to the end of the dictionary
        climatology_results.update({
            param: param_results
        })
    return climatology_results


def save_local_qc_tests(file_paths, refdes, stream, method, test_name):
    """ For a set of file paths and QC test name to run, this function
    loads each dataset individually and runs the requested local QC
    test. The parameters to test are determined from the variable names
    in the dataset containing "qartod_executed". The results of the test
    on all parameters that have relevant lookup tables on the
    oceanobservatories/qc-lookup repo on GitHub are saved in an
    xarray Dataset to a NetCDF file in the
    data/interim/{method}/{stream}/{refdes} directory in this project.
    
    Input:
    ------
    file_paths: list of strings, the paths to the external dataset files
        downloaded from OOINet.
    test_name: string, either "climatology" or "gross_range" for the
        desired QC test to run. These options both currently call the
        QARTOD tests from the US IOOS ioos_qc package.
        
    Returns:
    --------
    None
    """
    # Create copy of list of file paths for loading multi-file, single-deployment datasets
    paths_copy = file_paths.copy()
    
    while len(paths_copy) > 0:
        # Load data from a single deployment
        deploy_ds, deployment, paths_copy = get_deployment_ds(paths_copy) 
    
        # Create a dictionary of key-value pairs of dataset variable name:alternate parameter name
        test_parameters = get_test_parameters(deploy_ds)
        
        # Run local QARTOD test
        print(f'Running local QARTOD {test_name} test for deployment {deployment}.')  
        if len(test_parameters) > 0:
            qc_test_results = eval(f'run_qartod_{test_name}(refdes, stream, test_parameters, deploy_ds)')
            qc_test_results = timeseries_dict_to_xarray(qc_test_results, deploy_ds)
            # Add variable containing deployment number used in sorting and merging overlapping deployments
            qc_test_results['deployment'] = deploy_ds['deployment']
        else:
            pass

        # Build file name and directory for local QARTOD test results
        folder_path = os.path.join(os.path.abspath('../data/interim'), method, stream, refdes)
        os.makedirs(folder_path, exist_ok=True)
        test_results_path = os.path.join(folder_path, f"{test_name}_test-deployment00{deployment}.nc")

        # Save local test results
        qc_test_results.to_netcdf(test_results_path, mode='w')
    return


def run_comparison(ds, param, test_results, test):
    """Runs a comparison between the qartod gross range results returned
    as part of the dataset and results calculated locally.
    Original function created by Andrew Reed.
    
    edit:
        added parameter - test: string "gross_range" or "climatology"
        for test name to compare
    """
    # Get the local test results and convert to string type for comparison
    local_results = test_results[param].astype(str)
    
    # Run comparison
    not_equal = np.where(ds[f"{param}_qartod_{test}_test"] != local_results)[0]
    
    if len(not_equal) == 0:
        return None
    else:
        return not_equal


def get_mismatched_flags(expected_ds, local_ds, parameters, deployment,
                         test, expected_file):
    """ Runs comparison of local and expected QARTOD test flags, then
    exports local flags, expected flags, expected parameter value, and
    the corresponding datetime for any mismatched flags as an Xarray
    Dataset. Also adds percentage of total flags that disagree to a
    summary dictionary.
    
    Input:
    ------
        expected_ds: Xarray Dataset with flags for different QARTOD
            tests parsed into separate variables
        local_ds: Xarray Dataset containing only the resulting flags
            from running the QARTOD test locally
        mismatch: dictionary that will hold the results of the
            comparison
        deployment: 2-character string for the deployment number of the
            subsite
        test: string of the test to use for comparison, either
            "gross_range" or "climatology"
    
    Returns:
    --------
        mismatch: the updated dictionary with percentage of all flags
            in local QARTOD test that don't match expected QARTOD test
            flags for the current test and deployment
        
    Version 12 Sept 2023, Kylene M Cooley
    """
    # initialize dictionary to hold results of comparison and update
    # with the current deployment
    mismatch = {}
    mismatch.update({ "deployment" : f"{deployment}" })
    
    # Loop through parameters while updating dictionaries for
    # mismatched flags
    for param in parameters:

        # First, check that the test was applied
        test_name = f"{param}_qartod_{test}_test"
        if test_name not in expected_ds.variables:
            print(f"{test_name} not implemented")
            mismatch.update({ f"{param}" : np.nan })
        
        else:
            # Evaluate comparison of local test and expected test flags
            # to update dictionary of differences in the results
            print("Checking for mismatched QARTOD flags in "f"{param}")
            flag_mismatch = run_comparison(expected_ds, param, local_ds, test)

            if flag_mismatch is None:
                print("No mismatched values found")
                mismatch.update({ f"{param}" : np.nan })

            else:
                # Check if the parameter has an alternative ooinet_name
                # for later downloading QARTOD lookup tables
                if "alternate_parameter_name" in expected_ds[param].attrs:
                    ooinet_name = expected_ds[param].attrs[
                        "alternate_parameter_name"]
                else:
                    ooinet_name = param
                
                # Save datetime, expected, and local flag values as
                # DataArrays in a Dataset
                compare = xr.Dataset(
                    data_vars=dict(
                        expected_flags=(["time"], expected_ds[
                            f"{param}_qartod_{test}_test"][
                            flag_mismatch].values),
                        local_flags=(["time"], local_ds[param][
                            flag_mismatch].values),
                        expected_values=(["time"], expected_ds[f"{param}"][
                            flag_mismatch].values),
                    ),
                    coords=dict(
                        time=expected_ds['time'][flag_mismatch].values,
                    ),
                    attrs=dict(
                        parameter_name=f"{param}",
                        ooinet_name=f"{ooinet_name}",
                        percent_mismatched=f"{round(100*len(flag_mismatch)/len(expected_ds['time']))}%",
                        file_name=f"{expected_file}",
                    ),
                )
                export_filename="{}_comparison-deployment{}.nc".format(
                    test_name, deployment)
                export_path=os.path.join(os.path.abspath("../data/processed"),
                                         expected_ds.collection_method,
                                         expected_ds.stream,
                                         expected_ds.id[0:27], export_filename)
                compare.to_netcdf(export_path)
            
                mismatch.update({ f"{param}": {
                        'total' : f"{len(flag_mismatch)} ({round(100*len(flag_mismatch)/len(expected_ds['time']))}%)",
                }})
    return mismatch