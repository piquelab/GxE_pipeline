#!/bin/bash
set -v 
set -e

platenum=$1
root=../../derived_data/$platenum;

mydirs=(fastqs bams counts counts/GC counts/QC pileups)

for ii in ${mydirs[@]}; do
  mkdir -p $root/${ii}
done

cp Makefile $root/bams
cp counts_DEG.* $root/counts/GC
cp counts_QC_logs.* $root/counts/QC
cp Makefile.pileup $root/pileups/Makefile
cp ai_preprocessing.R $root/pileups

cd $root/fastqs
pwd

if [[ $platenum == "D"* ]]; then
    platenum1=${platenum//D/}; 
    findstrprefix=D*${platenum1};
else
    findstrprefix=$platenum
fi

for bc in {1..96}; do 
    find /wsu/home/groups/piquelab/OurData/ -name "${findstrprefix}-HT${bc}_*R1*.fastq.gz" \
    	| awk '{print $0,NR}' \
    	| while read f r; do  
    	    b=${f//R1/R2}; 
    	    echo "${platenum}-HT${bc}_L${r}";
    	    ln -s $f ${platenum}-HT${bc}_L${r}_R1.fastq.gz; 
    	    ln -s $b ${platenum}-HT${bc}_L${r}_R2.fastq.gz; 
    done; 
done;
