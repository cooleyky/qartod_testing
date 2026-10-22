#!/usr/bin/bash
# gi_qartod_crossref.sh - runs the qartod tests 
# cross-reference Python script for all Global
# Irminger Sea Array sites. Adapted from 
# gs_qartod_crossref.sh.
#
# 6 Oct 2026
# Kylene Cooley
#

gi_qartod_crossref() {
    conda activate qartod_test 
    python run_qartod_test_cross-ref.py GI01SUMO
    python run_qartod_test_cross-ref.py GI02HYPM
    python run_qartod_test_cross-ref.py GI03FLMA
    python run_qartod_test_cross-ref.py GI03FLMB
    python run_qartod_test_cross-ref.py GI05MOAS
}
