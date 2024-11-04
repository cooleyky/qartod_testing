#!/usr/bin/env qartod_test
#
# run_qartod_test_cross-ref.py

# Import libraries
import numpy as np
import pandas as pd
import xarray as xr
import requests
import io
import ast
from glob import glob

# Import functions from ooinet and ooi_data_explorations libraries
from ooi_data_explorations.common import load_kdata, get_vocabulary
from ooinet import M2M

# Import functions from ooinet and ooi_data_explorations libraries
from ooi_data_explorations.common import load_kdata, get_vocabulary
from ooinet import M2M

# Import functions from project qc_completion module
from qartod_testing.qc_completion import load_gross_range_qartod_test_list, \
    load_climatology_qartod_test_list, make_test_parameter_dict, \
    check_tests_exe, make_results_table, add_test_exe, write_results

# Define site for refdes search and find datasets available
site = 'GA01SUMO'
datasets = M2M.search_datasets(site)
datasets.reset_index(inplace=True)
datasets.drop(labels="index", axis=1, inplace=True)

# Set csv save directory and file name for results
csv_name = f"{site}_test_cross-ref_results.csv"
csv_dir = "./data/processed/"
# loop through sensors to check and find datastreams available
for k in datasets.index:
    refdes = datasets.refdes[k]
    datastreams = M2M.get_datastreams(refdes)
    # loop through datastreams and first deployment available
    site, node, sensor = refdes.split("-", 2)
    for m in datastreams.index:
        method = datastreams.method[m]
        stream = datastreams.stream[m]
        deploy = datasets.deployments[k][0]
        instclass = sensor[3:8]
        # Load gross range and climatology test tables
        grt_table = load_gross_range_qartod_test_list(refdes, stream)
        ct_table = load_climatology_qartod_test_list(refdes, stream)
        if (grt_table is False) and (ct_table is False):
            print(f"No existing qc-lookup table for {refdes}-{stream}.")
        else:
            # Load data
            get_vocabulary(site, node, sensor)
            data = load_kdata(site, node, sensor, method, stream, ('*deployment%04d*%s*.nc' % (deploy, instclass)))
            try:
                print(data.id) # to check data stream loaded
            except AttributeError:
                print(f"No dataset available for {refdes}-{stream} tests.")
                del [grt_table, ct_table]
            else:
                # Create a dictionary of key-value pairs of dataset variable name:alternate parameter name
                test_parameters = make_test_parameter_dict(data)
                # Use test parameters to check for tests executed in dataset
                test_exe = check_tests_exe(data, test_parameters, grt_table, ct_table)
                # Make table for cross-ref results
                table = make_results_table(grt_table, ct_table)
                # Add column with QARTOD tests executed by parameter
                table = add_test_exe(table, test_exe)
                # Write QARTOD test cross-reference results table to a CSV
                write_results(table, csv_name, csv_dir)
                del [grt_table, ct_table, data, test_parameters,
                     test_exe, table]