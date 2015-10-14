#!/bin/bash
## Run QuASAR for all individuals of a given cell type

cellType=$1
ncpus=4
queue=mmtxq
jobName=Qprep-${cellType}

echo "cd ${PWD}; Rscript QuASAR_prep.R ${cellType}" | qsub -q ${queue} -l nodes=1:ppn=${ncpus} -N ${jobName} -o ${jobName}.Qsub -e ${jobName}.Qsub.e 
