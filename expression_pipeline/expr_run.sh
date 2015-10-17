#!/bin/bash
set -e

ct=$1
al=$2
wd=$3
qq=mmtxq
qp=12

for pt in `cat ../../derived_data/covariates/GxE_DP*.txt | grep $ct | cut -f 1 | sort | uniq`
do
    jid=`echo "cd ${PWD}; Rscript expr_countsToFpkms.R platePrefix=${pt} bamMethod=${al} cores=${qp}" | qsub -q ${qq} -l ncpus=${qp} -N ${pt}_fpkm`
    jobs=("${jobs[@]}" "${jid}")
done

## Submit the next script to run after the previous ones finish
echo "cd ${PWD}; Rscript expr_getMeanAndTopExpr.R celltype=${ct} word=${al}.${wd}" | qsub -q ${qq} -l ncpus=${qp} -N expr_${ct} -W depend=afterok:`echo ${jobs[@]} | sed s,\ ,:,g`
