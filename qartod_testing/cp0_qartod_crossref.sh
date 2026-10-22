#!/usr/bin/bash
# cp0_qartod_crossref.sh - runs the qartod tests 
# cross-reference Python script for all Coastal
# Pioneer NES Array sites. Adapted from 
# gs_qartod_crossref.sh.
#
# 6 Oct 2026
# Kylene Cooley
#

cp0_qartod_crossref() {
    conda activate qartod_test 
    python run_qartod_test_cross-ref.py CP01CNPM
    python run_qartod_test_cross-ref.py CP01CNSM
    python run_qartod_test_cross-ref.py CP01CNSP
    python run_qartod_test_cross-ref.py CP02PMCI
    python run_qartod_test_cross-ref.py CP02PMCO
    python run_qartod_test_cross-ref.py CP02PMUI
    python run_qartod_test_cross-ref.py CP02PMUO
    python run_qartod_test_cross-ref.py CP03ISPM
    python run_qartod_test_cross-ref.py CP03ISSM
    python run_qartod_test_cross-ref.py CP03ISSP
    python run_qartod_test_cross-ref.py CP04OSPM
    python run_qartod_test_cross-ref.py CP04OSSM
    python run_qartod_test_cross-ref.py CP05MOAS
}
