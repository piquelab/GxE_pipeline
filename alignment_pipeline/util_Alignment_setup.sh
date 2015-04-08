#!/bin/bash
set -v 
set -e

plate=$1 

mkdir -p ../D$plate/fastqs
mkdir -p ../D$plate/bams
mkdir -p ../D$plate/counts
mkdir -p ../D$plate/counts/GC
mkdir -p ../D$plate/counts/QC
mkdir -p ../D$plate/pileups

cp Makefile ../D$plate/bams
cp counts_DEG.* ../D$plate/counts/GC
cp counts_QC_logs.* ../D$plate/counts/QC
cp Makefile1 ../D$plate/pileups

cd ../D$plate/fastqs
pwd

for bc in {1..96}; do 
	#find /nfs/rprscratch/OurData/ -name "*${plate}-HT${bc}_*R1*.fastq.gz" \
	## DANGER the line below resulted from a mis-labeling of the original fastqs. They correspond to D2P1
	## If you run this script with another plate name they will be named with that prefix, so don't. Please & thank you.
	#find ~/piquelab/OurData/Nextseq/ -name "D*${plate}-HT${bc}_*R1*.fastq.gz" \
	## D2P1 was named with 'Luca-HT . . .'
	#find /wsu/home/groups/piquelab/OurData/140404_SN7001329_0366_AH94AWADXX/FastQ/ -name "Luca-HT${bc}_*R1*.fastq.gz" \	
	#find /wsu/home/groups/piquelab/OurData/140404_SN7001329_0366_AH94AWADXX/FastQ/ -name "Luca-HT${bc}_*R1*.fastq.gz" \
	find /wsu/home/groups/piquelab/OurData/ -name "D*${plate}-HT${bc}_*R1*.fastq.gz" \
		| awk '{print $0,NR}' \
		| while read f r; do  
			b=${f//R1/R2}; 
			echo "D${plate}-HT${bc}_L${r}";
			ln -s $f D${plate}-HT${bc}_L${r}_R1.fastq.gz; 
			ln -s $b D${plate}-HT${bc}_L${r}_R2.fastq.gz; 
		done; 
done;
