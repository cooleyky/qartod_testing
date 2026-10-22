#!/usr/bin/bash
# gp_qartod_crossref.sh - runs the qartod tests 
# cross-reference Python script for all Global
# Station PAPA Array sites. Adapted from 
# gs_qartod_crossref.sh.
#
# 6 Oct 2026
# Kylene Cooley
#

gp_qartod_crossref() {
    conda activate qartod_test 
    python run_qartod_test_cross-ref.py GP01SUMO # FDCHP only
    python run_qartod_test_cross-ref.py GP02HYPM
    python run_qartod_test_cross-ref.py GP03FLMA
    python run_qartod_test_cross-ref.py GP03FLMB
    python run_qartod_test_cross-ref.py GP05MOAS
}
