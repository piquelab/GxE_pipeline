#!/bin/bash
set -v
set -e

## Set up joint genotyping directory
mkdir -p ../../jointGenotyping/
cp QuASAR_* ../../jointGenotyping/

cat ../../derived_data/covariates/GxE_DP* | grep -v 'Plate.ID' | cut -f 9 | sort | uniq | \
    while read plate; do

    ## results directory & QuASAR scripts
    mkdir -p ../../jointGenotyping/QuASAR_results_${plate}/data

    done


