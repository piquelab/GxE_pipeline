#!/bin/bash
plate=$1
ncpus=4
queue=mmtxq
jobName=QC_${plate}

echo "cd ${PWD}; Rscript counts_QC_logs.R ${plate} ../../bams/" | qsub -q ${queue} -l nodes=1:ppn=${ncpus} -N ${jobName} -o ${jobName}.Qsub -e ${jobName}.Qsub.e 
