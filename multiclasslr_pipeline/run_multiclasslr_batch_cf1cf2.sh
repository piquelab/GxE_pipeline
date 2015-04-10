#!/bin/bash
dat=$1
tag=$2
ncpus=4
queue=mmtxq
jobName=multlr


echo "cd ${PWD}; Rscript multiclasslr_batch_cf1cf2.R ${dat} ${tag}" | qsub -q ${queue} -l nodes=1:ppn=${ncpus} -N ${jobName} -o ${jobName}.Qsub -e ${jobName}.Qsub.e 
