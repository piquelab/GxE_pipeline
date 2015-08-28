#!/bin/bash
plate=$1
ncpus=12
queue=mmtxq
jobName=${plate}_DEG

echo "cd ${PWD}; Rscript GxE_DESeq2_automate.R ${plate} ${ncpus}" | qsub -q ${queue} -l nodes=1:ppn=${ncpus} -N ${jobName} -o ${jobName}.Qsub -e ${jobName}.Qsub.e
