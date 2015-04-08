#!/bin/bash
set -e 
set -v

topdir=$1
myArray=( 1 2 4 5 6 7 8 9 11 12 )

for ii in ${myArray[@]}; do ls ../../QuASAR_results_DP${ii}/output/*masterTable_logFC.txt | while read f; do less $f | grep -v "plate"; done; done > MESH_QuASAR_master_logFC.txt

## attach controls to the treatment data
Rscript MESH_QuASAR_master_assignControls.R ${topdir}
