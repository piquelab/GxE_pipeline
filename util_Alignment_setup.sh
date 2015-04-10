#!/bin/bash
set -v 
set -e

platenumb=$1 

if [ "$2" = "-w" ] 
then
  plate=D$1-WASP;
else
  plate=D$1;
fi

mkdir -p ../$plate/fastqs
mkdir -p ../$plate/bams
mkdir -p ../$plate/counts
mkdir -p ../$plate/counts/GC
mkdir -p ../$plate/counts/QC
mkdir -p ../$plate/pileups

if [ "$2" = "-w" ] 
then
  cp /WASP/Makefile ../$plate/bams;
else
  cp Makefile ../$plate/bams;
fi


cp counts_DEG.* ../$plate/counts/GC
cp counts_QC_logs.* ../$plate/counts/QC
cp Makefile1 ../$plate/pileups

cd ../$plate/fastqs
pwd

for bc in {1..96}; do 
	#find /nfs/rprscratch/OurData/ -name "*${plate}-HT${bc}_*R1*.fastq.gz" \
	## DANGER the line below resulted from a mis-labeling of the original fastqs. They correspond to D2P1
	## If you run this script with another plate name they will be named with that prefix, so don't. Please & thank you.
	#find ~/piquelab/OurData/Nextseq/ -name "D*${plate}-HT${bc}_*R1*.fastq.gz" \
	## D2P1 was named with 'Luca-HT . . .'
	#find /wsu/home/groups/piquelab/OurData/140404_SN7001329_0366_AH94AWADXX/FastQ/ -name "Luca-HT${bc}_*R1*.fastq.gz" \	
	#find /wsu/home/groups/piquelab/OurData/140404_SN7001329_0366_AH94AWADXX/FastQ/ -name "Luca-HT${bc}_*R1*.fastq.gz" \
	find /wsu/home/groups/piquelab/OurData/ -name "D*${platenum}-HT${bc}_*R1*.fastq.gz" \
		| awk '{print $0,NR}' \
		| while read f r; do  
			b=${f//R1/R2}; 
			echo "D${plate}-HT${bc}_L${r}";
			ln -s $f D${platenum}-HT${bc}_L${r}_R1.fastq.gz; 
			ln -s $b D${platenum}-HT${bc}_L${r}_R2.fastq.gz; 
		done; 
done;
