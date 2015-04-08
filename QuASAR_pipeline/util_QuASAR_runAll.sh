#!/bin/bash
set -v
set -e

myArray=( 1 2 4 5 6 7 8 9 11 12 )

for ii in ${myArray[@]}; do 
    cd QuASAR_results_DP${ii}; 
    sh QuASAR_pipeline_all.sh DP${ii}; 
    cd ../; 
done
