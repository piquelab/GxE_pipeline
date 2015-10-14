#!/bin/bash
## Run QuASAR for all individuals of a given cell type

cellType=$1
ncpus=4
queue=mmtxq
jobName=QuASAR-${cellType}

echo "cd ${PWD}; Rscript QuASAR_pipeline.R ${cellType}" | qsub -q ${queue} -l nodes=1:ppn=${ncpus} -N ${jobName} -o ${jobName}.Qsub -e ${jobName}.Qsub.e 
