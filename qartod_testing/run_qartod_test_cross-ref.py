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
from ooi_data_explorations.common import load_kdata, get_vocabulary, \
    m2m_request, m2m_collect
from ooinet import M2M

# Import functions from project qc_completion module
from qartod_testing.qc_completion import load_gross_range_qartod_test_list, \
    load_climatology_qartod_test_list, make_test_parameter_dict, \
    check_tests_exe, make_results_table, add_test_exe, write_results, \
    check_skip_kw

# Define site for refdes search and find datasets available
site = 'GS02HYPM'
prefix = f"{site[0:2]}"
datasets = M2M.search_datasets(site)
datasets.reset_index(inplace=True)
datasets.drop(labels="index", axis=1, inplace=True)

# Set csv save directory and file name for results
csv_name = f"{site}_test_cross-ref_results.csv"
csv_dir = f"./data/processed/{prefix}_tests_completed3/"
# loop through sensors to check and find datastreams available
for k in datasets.index:
    refdes = datasets.refdes[k]
    # Skip this refdes if it contains a class keyword
    if check_skip_kw(refdes, "class") != -1:
        print("skipped "+refdes)
        continue
    datastreams = M2M.get_datastreams(refdes)
    # loop through datastreams and first deployment available
    site, node, sensor = refdes.split("-", 2)
    for m in datastreams.index:
        method = datastreams.method[m]
        stream = datastreams.stream[m]
        deploy = datasets.deployments[k][0]
        instclass = sensor[3:8]

        # Skip this stream if it contains a stream keyword
        if check_skip_kw(stream, "stream") != -1:
            continue
        
        # Load gross range and climatology test tables
        grt_table = load_gross_range_qartod_test_list(refdes, stream)
        ct_table = load_climatology_qartod_test_list(refdes, stream)
        if (grt_table is False) and (ct_table is False):
            print(f"No existing qc-lookup table for {refdes}-{stream}.")
        else:
            # Load data
            get_vocabulary(site, node, sensor)
            print(f"Loading deployment {deploy}")
            try:
                data = load_kdata(site, node, sensor, method, stream, ('*deployment%04d*%s*.nc' % (deploy, instclass)))
            except:
                print(f"Loading deployment {deploy} from kdata failed")
            while data is None:
                deploy+=1
                get_vocabulary(site, node, sensor)
                if deploy > datasets.deployments[k][-1]:
                    print(f"No dataset available for {refdes}-{stream} tests.")
                    break
                else:
                    print(f"Loading deployment {deploy}")
                    try:
                        data = load_kdata(site, node, sensor, method, stream, ('*deployment%04d*%s*.nc' % (deploy, instclass)))
                    except:
                        print(f"Loading deployment {deploy} from kdata failed")
            try:
                print(data.id) # to check data stream loaded
            except AttributeError:
                del [grt_table, ct_table]
            else:
                # Create a dictionary of key-value pairs of dataset variable name:alternate parameter name
                test_parameters = make_test_parameter_dict(data)
                # If QARTOD test not executed on kdata then load data with M2M request
                if len(test_parameters)==0:
                    deploy_info = M2M.get_deployments(refdes, deploy_num=str(deploy))
                    dt_start = deploy_info.deployStart[0] + pd.Timedelta(5, 'D')
                    dt_end = deploy_info.deployStart[0] + pd.Timedelta(20, 'D')
                    dt_start = dt_start.isoformat(timespec='milliseconds')+'Z'
                    dt_end = dt_end.isoformat(timespec='milliseconds')+'Z'
                    m2m_result = m2m_request(site, node, sensor, method, stream,
                         start=dt_start, stop=dt_end)
                    try:
                        m2m_data = m2m_collect(m2m_result, tag=('.*deployment%04d.*%s.*\.nc$' % (deploy, instclass)))
                    except:
                        print("M2M data request or collection failed")
                    else:
                        data = m2m_data
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