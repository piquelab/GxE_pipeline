#!/bin/bash
data=$1
ncpus=4
queue=mmtxq
jobName=MESH_all

echo "cd ${PWD}; python run_analysis.py ${data}; Rscript calc_posteriors.R" | qsub -q ${queue} -l nodes=1:ppn=${ncpus} -N ${jobName} -o ${jobName}.Qsub -e ${jobName}.Qsub.e 
