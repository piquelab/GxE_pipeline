#!/bin/bash
set -e
set -v

plate=$1 # e.g., DP1
bams=$2  # e.g., 'bwa', indicaing aligner used
word=$3  # e.g., 'counts.fpkm.adj', indicating processing steps
type=$4  # for makeMasterTable, 'all' transcripts or one 'perGene'

cores=$5
if [ -z "$cores" ]; then 
    cores=1;
fi

## Convert counts into fpkms
R --vanilla --args ${plate} ${bams} ${cores} < expr_countsToFpkms.R

## Get the mean expression per transcipt (per condition) and the top transcript per gene
R --vanilla --args ${plate} ${bams}.${word} ${cores} < expr_getMeanAndTopExpr.R

## Add logFC data to make a per-plate 'master' expression table
R --vanilla --args ${plate} ${bams}.${word} ${type} < expr_makeExpressionTable.R
