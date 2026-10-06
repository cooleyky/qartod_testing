#!/usr/bin/bash
# gs_qartod_crossref.sh - runs the qartod tests 
# cross-reference Python script for all Global
# Southern Ocean Array sites.
#
# 15 Jan 2025
# Kylene Cooley
#
# Updated 30 Jan 2026: Uses absolute path variables
# Updated 6 Oct 2026: Run Python script for specifically 
# the CGSN GS array, all sites. Python script now logs
# all messages at level logging.info and above.

gs_qartod_crossref() {
    conda activate qartod_test 
    python run_qartod_test_cross-ref.py GS01SUMO
    python run_qartod_test_cross-ref.py GS02HYPM
    python run_qartod_test_cross-ref.py GS03FLMA
    python run_qartod_test_cross-ref.py GS03FLMB
    python run_qartod_test_cross-ref.py GS05MOAS
}