#!/bin/bash
dat=$1
tag=$2
ncpus=4
queue=mmtxq
jobName=multlrall


echo "cd ${PWD}; Rscript multiclasslr_batch.R ${dat} ${tag}" | qsub -q ${queue} -l nodes=1:ppn=${ncpus} -N ${jobName} -o ${jobName}.Qsub -e ${jobName}.Qsub.e 
