
Alignment pipeline
===================
This directory contains scripts related to the alignment and quality filtering of sequencing data for the GxE data.

## Setup
Initialize the alignment directory (one per plate) using the utility script:
```
 $ sh util_Alignment_setup.sh DP1
```
For the indicated plate, the setup script will setup the directory and copy the relevant scripts. Additionally, it will create links to the fastq files in the fastqs/ directory.


## Prerequisites
* A [Covariate file](../misc/) explaining each barcde and linking treatments to controls
* A [gene annotation](../misc/) file with gene and transcript information

### Software packages needed
* [Samtools](http://www.htslib.org/)
* [Bedtools](http://bedtools.readthedocs.org/en/latest/)
* [STAR Aligner](https://github.com/alexdobin/STAR)
* [WASP Suite](https://github.com/bmvdgeijn/WASP)

## Outline
1. Align the raw sequencing reads to the genome
2. Merge and quality filter aligned reads
3. Determine which positions need remapping
4. Filter reads with alleleic bias & remove PCR duplicates
5. Compile read counts across filtering stages (QC)
6. Count reads overlapping genes
7. Count reads and alleles overlapping known variants


## Notes
Each step is run on one plate at a time. Multiple plates are setup through the utility script and placed in separate directories.

By default, the pipeline is built to submit jobs to a PBS scheduling system. This is built into the Makefile, and is achieved elsewhere through shell scripts.

---
## Full Pipeline
##### Step 1: Align the raw sequencing reads to the genome

```
directory: bams/
script: Makefile
dependencies: the folder ../fastqs with symlinked fastqs (done in setup script)
in: -
out: one sorted bam file per *R1*/*R2* fastq pair
```

Example:
```
 $ make all_bams
```

##### Step 2: Merge and quality filter aligned reads
Merge technical replicates corresponding to the same barcode and apply a quality filter

```
directory: bams/
script: Makefile
dependencies: Aligned and sorted bam files
in: -
out: clean bam files & .txt files with read counts from each stage of filtering
```

Example:
```
 $ make all_quality
```

##### Step 3: Determine which positions need remapping
Using WASP, determine which reads overlap SNPs and should be remapped to reduce allele-specific mapping biases. This step may take considerable time, depending on the number of reads and the number of SNPs.

```
directory: bams/
script: Makefile
dependencies: Quality filtered bam files
in: -
out: quality filtered bams, list of SNP positions
```

Example:
```
 $ make all_fq
```

##### Step 4: Filter reads with alleleic bias & remove PCR duplicates
With the results from Step 3, remove reads that map to different positions depending on the allele. Finally, remove PCR duplicates.

```
directory: bams/
script: Makefile
dependencies: WASP output from Step 3
in: -
out: clean bams
```

Example:
```
 $ make all_clean
```


##### Step 5: Compile read counts across filtering stages (QC)
Create a histogram & table from read counts at each stage of .bam processing

```
directory: counts/QC/
script: counts_QC_logs.R (counts_QC_logs.sh to submit)
dependencies: clean bam files
in: plate number
out: <plate>_QC_counts.pdf & <plate>_QC_counts.txt
```

Example:
```
 $ sh counts_QC_logs.sh DP1
```
or,
```
 $ Rscript counts_QC_logs.R DP1
```

##### Step 6: Count reads overlapping genes
Use bedtools coverage to calculate the number of reads overlapping transcripts per sample

```
directory: counts/GC/
script: counts_DEG.R (counts_DEG.sh to submit)
dependencies: cleam bam files
in: plate number, location of bam files, # of cores
out: <plate>_QC_counts.pdf & <plate>_QC_counts.txt
```

Example:
```
 $ sh counts_DEG.sh DP1
```
or,
```
 $ Rscript counts_DEG.R 12 DP1 ../../bams
```

##### Step 7: Count reads and alleles overlapping known variants
Make pileup files from each clean bam file. Requires a list of SNPs to overlap (see [misc/](./misc/)).

``` 
directory: pileups/ 
script: Makefile
dependencies: clean bam files, SNP list
in: - 
out: clean pileups ready for input to QuASAR
```

Example:
```
 $ make all_pileups
```

