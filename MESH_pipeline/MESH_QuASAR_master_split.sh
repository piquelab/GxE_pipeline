#!/bin/bash
## warning: each parition of the data needs to be run in its own analysis folder
set -v
set -e

data=$1
outdir=../analysis
mkdir -p $outdir

less ${data} | grep -v "snp" | awk '{print $1,$2,$3,$4"\n"$1,$7,$8,$9}' | tr " " "\t" > ${outdir}/allPlates_allSnps.txt
cp ../../../GxE_pipeline/MESH_pipeline/MESH/* ${outdir}
echo Processing complete
