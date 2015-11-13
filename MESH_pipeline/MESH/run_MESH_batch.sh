#!/bin/bash
data=$1
table=$2
ncpus=4
queue=mmtxq
jobName=MESH_all

echo "cd ${PWD}; python run_analysis.py ${data}; Rscript MESH_calc_posteriors.R; Rscript MESH_plots.R ${table}" | qsub -q ${queue} -l nodes=1:ppn=${ncpus} -N ${jobName} -o ${jobName}.Qsub -e ${jobName}.Qsub.e 
