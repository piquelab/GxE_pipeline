#!/bin/bash

vb=0
while getopts "p:vh" opt; do
    case $opt in
	p) plate=$OPTARG;;
	v) vb=1;;
	h)
	    echo "useage: checkAlignment.sh -p PLATE [-v]"
	    echo ""
	    echo "-p      PLATE to check, e.g., DP1"
	    echo "-v      Verbose mode"
	    echo "-h      Print help and exit"
	    exit
	    ;;
	\?) echo "Invalid option: -$OPTARG" >2&;;
    esac
done

fqDir=../../derived_data/${plate}/fastqs/
bmDir=../../derived_data/${plate}/bams/

## Check that all fastq pairs were aligned, and that none of the aligned 
## bams are truncated (samtools will warn, at least in v1.2), by saying
## "[bam_header_read] EOF marker is absent. The input is probably truncated."
echo "Checking initial alignments..."
find ${fqDir} -name '*_R1.fastq.gz' | while read fq; do
    stem=`echo ${fq%_R1.fastq.gz} | sed s,.*/,,`
    if [ $vb -eq 1 ]; then
    	echo "Checking ${stem}"
    fi
    samtools view ${bmDir}/${stem}.out/${stem}.Aligned.sort.bam | head -0
done

## Check that all bams were merged & quality filtered properly. Merged bams are
## removed (intermediate to quality), but check that the index exists to be sure
## the merge occured without errors. Then check the quality file.
echo "Checking merged/quality bams..."
find ${fqDir} -name '*_L1_R1.fastq.gz' | while read fq; do
    stem=`echo ${fq%_L1_R1.fastq.gz} | sed s,.*/,,`
    if [ $vb -eq 1 ]; then
    	echo "Checking ${stem}"
    fi
    if [ -f "${bmDir}/${stem}_merged.bam.bai" ]; then
	samtools view ${bmDir}/${stem}_quality.bam | head -0
    else
	echo "${stem}_merged.bam.bai not found!"
    fi
done
