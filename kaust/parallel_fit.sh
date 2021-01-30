#!/bin/bash
parallel Rscript ~/GitHub/RandomFields_cracked/kaust/fit_kaust2.R ::: $(seq 9 16) &