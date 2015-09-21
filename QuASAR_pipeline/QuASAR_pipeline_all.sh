#!/bin/bash
plate=$1
less ../../derived_data/covariates/GxE_${plate}_covariates.txt \
	| awk '$10!="CellLine" {print $10}' \
	| sort \
	| uniq \
	| while read cl; do sh QuASAR_pipeline.sh ${plate} $cl; done  
