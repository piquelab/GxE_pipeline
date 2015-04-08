#!/bin/bash
plate=$1
ncpus=4
queue=mmtxq
jobName=QC_counts

echo "cd ${PWD}; Rscript counts_QC_logs.R ${plate}" | qsub -q ${queue} -l nodes=1:ppn=${ncpus} -N ${jobName} -o ${jobName}.Qsub -e ${jobName}.Qsub.e 
