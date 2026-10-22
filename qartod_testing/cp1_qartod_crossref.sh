#!/usr/bin/bash
# cp1_qartod_crossref.sh - runs the qartod tests 
# cross-reference Python script for all Coastal
# Pioneer MAB Array sites. Adapted from 
# cp0_qartod_crossref.sh.
#
# 6 Oct 2026
# Kylene Cooley
#

cp1_qartod_crossref() {
    conda activate qartod_test 
    python run_qartod_test_cross-ref.py CP10CNSM
    python run_qartod_test_cross-ref.py CP11NOSM
    python run_qartod_test_cross-ref.py CP11SOSM
    python run_qartod_test_cross-ref.py CP12CNSW
    python run_qartod_test_cross-ref.py CP12WESW
    python run_qartod_test_cross-ref.py CP13EAPM
    python run_qartod_test_cross-ref.py CP13NOPM
    python run_qartod_test_cross-ref.py CP13SOPM
    python run_qartod_test_cross-ref.py CP14NEPM
    python run_qartod_test_cross-ref.py CP14SEPM
    python run_qartod_test_cross-ref.py CP15MOAS
    python run_qartod_test_cross-ref.py CP16MOAS
}
