#!/bin/bash
parallel Rscript kaust/fit_2a.R ::: $(seq 1 2) &