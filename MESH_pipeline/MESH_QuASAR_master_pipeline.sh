#!/bin/bash

## topdir := directory with all plate specific analysis, such as ~/piquelab/charvey/GxE/new_system/jointGenotyping
topdir=$1

ncpus=4
queue=mmtxq
jobName=mstrtbl

echo "cd ${PWD}; sh MESH_QuASAR_master.sh  ${topdir}" | qsub -q ${queue} -l nodes=1:ppn=${ncpus} -N ${jobName} -o ${jobName}.Qsub -e ${jobName}.Qsub.e 
