#!/usr/bin/bash
# qartod_crossref - runs the qartod tests cross-reference
# Python script and saves the standard output to a 
# text file named with the OOI array site tested.
#
# 15 Jan 2025
# Kylene Cooley
#
# Updated 30 Jan 2026: Uses absolute path variables
# TO-DO: figure out how to run all commands needed within script in a function.

# inside_script() {
#     conda activate qartod_test
#     python run_qartod_test_cross-ref.py $site_name
# }


qartod_crossref() {
    local site_name=$1
    script ~/code/qartod_testing/data/processed/$site_name-output.txt -c "conda init && conda activate qartod_test && python run_qartod_test_cross-ref.py $site_name"
}