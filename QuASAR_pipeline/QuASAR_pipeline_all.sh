#!/bin/bash
plate=$1
less ../../derived_data/covariates/GxE_${plate}_covariates.txt \
	| awk '$13!="CellLine" {print $13}' \
	| sort \
	| uniq \
	| while read cl; do sh QuASAR_pipeline.sh ${plate} $cl; done  
