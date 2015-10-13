#!/bin/bash
set -v
set -e

## Set up joint genotyping directory
mkdir -p ../../jointGenotyping/
cp QuASAR_prep.* ../../jointGenotyping/

cat ../../derived_data/covariates/GxE_DP* | grep -v 'Plate.ID' | cut -f 1 | sort | uniq | \
    while read plate; do

    ## results directory & QuASAR scripts
    mkdir -p ../../jointGenotyping/QuASAR_results_${plate}/data
    cp QuASAR_pipeline.* ../../jointGenotyping/QuASAR_results_${plate}

    ## link to all files for analysis/
    cd ../../jointGenotyping/QuASAR_results_${plate}/data/
    ls ../../../derived_data/${plate}/pileups/*.pileup.clean.bed.gz | while read f; do ln -svf $f; done
    cd -

    done


