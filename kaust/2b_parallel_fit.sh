#!/bin/bash
parallel Rscript kaust/fit_2b.R ::: $(seq 1 2) &