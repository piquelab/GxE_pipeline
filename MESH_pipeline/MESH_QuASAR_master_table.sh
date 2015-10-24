#!/bin/bash
set -e
set -v

ncpus=12
queue=mmtxq
jobName=MESH-MT

## MESH_QuASAR_makeMasterTable.R is very customizable, to use non-default parameters
## modify call here, or run interactively
echo "cd ${PWD}; Rscript MESH_QuASAR_makeMasterTable.R cores=${ncpus}" | qsub -q ${queue} -l nodes=1:ppn=${ncpus} -N ${jobName} -o ${jobName}.Qsub -e ${jobName}.Qsub.e
