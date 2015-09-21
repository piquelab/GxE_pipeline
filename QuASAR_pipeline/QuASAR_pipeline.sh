#!/bin/bash
plate=$1
cellLine=$2
ncpus=2
queue=mmtxq
jobName=${plate}-${cellLine}

echo "cd ${PWD}; Rscript QuASAR_pipeline.R ${plate} ${cellLine}" | qsub -q ${queue} -l nodes=1:ppn=${ncpus} -N ${jobName} -o ${jobName}.Qsub -e ${jobName}.Qsub.e 
