#!/bin/bash
## Run QuASAR on all cell types

less ../../derived_data/covariates/GxE_DP*_covariates.txt \
	| awk '$9!="CellType" {print $9}' \
	| sort \
	| uniq \
	| while read cl; do sh QuASAR_pipeline.sh $cl; done  
