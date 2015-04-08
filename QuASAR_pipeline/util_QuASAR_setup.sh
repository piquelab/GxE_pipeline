#!/bin/bash
set -v
set -e

plate=$1

## results directory & QuASAR scripts
mkdir -p ../../jointGenotyping/QuASAR_results_${plate}/data
cp QuASAR_* ../../jointGenotyping/QuASAR_results_${plate}
## this is redundant. make it better by checking for util_QuASAR_runAll.sh before copying it.
cp util_QuASAR_runAll.sh ../../jointGenotyping/

## link to all files for analysis
cd ../../jointGenotyping/QuASAR_results_${plate}/data
ls /wsu/home/groups/piquelab/charvey/GxE/derived_data/${plate}/pileups/*.clean.bed.gz | while read f; do ln -sv $f; done
