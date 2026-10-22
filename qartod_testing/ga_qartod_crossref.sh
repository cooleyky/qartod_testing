#!/usr/bin/bash
# ga_qartod_crossref.sh - runs the qartod tests 
# cross-reference Python script for all Global
# Argentine Basin Array sites. Adapted from 
# gs_qartod_crossref.sh.
#
# 6 Oct 2026
# Kylene Cooley
#

ga_qartod_crossref() {
    conda activate qartod_test 
    python run_qartod_test_cross-ref.py GA01SUMO
    python run_qartod_test_cross-ref.py GA02HYPM
    python run_qartod_test_cross-ref.py GA03FLMA
    python run_qartod_test_cross-ref.py GA03FLMB
    python run_qartod_test_cross-ref.py GA05MOAS
}
