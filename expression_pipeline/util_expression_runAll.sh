#!/bin/bash
set -e
set -v

bams=$1 # e.g., 'bwa', indicating aligner used
word=$2 # e.g., 'counts.fpkm.adj', indicating processing steps
type=$3 # for makeMasterTable, 'all' transcripts or one 'perGene'
queue=$4

cores=$5
if [ -z "$cores" ]; then
    cores=1;
fi

## Call the expression pipeline once per plate
myArray=( 1 2 4 5 6 7 8 9 10 11 12 )
for ii in ${myArray[@]}
do
    echo "cd ${PWD}; sh expr_run.sh DP${ii} ${bams} ${word} ${type} ${cores}" | qsub -q ${queue} -l ncpus=${cores} -N expr_DP$ii;
done

