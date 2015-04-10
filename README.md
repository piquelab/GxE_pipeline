# GxE-pipeline
This repo contains scripts that comprise the Gene Environment Interaction (GxE) pipeline, for alignment & statistical modeling of allelic imbalances in RNA-Seq data across many cellular environments.

## Outline
### I. Align RNA-seq reads and perform QC
  * 1.) Create directories, symbolic link to fastqs, and distribute scripts for alignment
  * 2.) Align reads, merge across barcodes, & do read counts
  * 3.) Combine read counts from all stages of filtering
  * 4.) Gene counts using DESeq2
  * 5.) Make pileups for QuASAR

### II. QuASAR (EM) pipeline to infer individual ASE
  * 6.) Create directories, symbolic link to pileups, and distribute scripts for QuASAR
  * 7.) Run the QuASAR pipeline for joint genotyping ASE inference   

### II. MESH (MCMC) to infer condition specific ASE
  * 8.) Create directories, setup the master table, and distribute scripts for MESH 
  * 9.) Combine QuASAR data into a master table and process data for MESH
  * 10.) Split master table data on logFC for MESH analysis
  * 11.) Run the MESH pipeline for condition specific ASE

### V. Multiclass logistic regression on MESH output and DESeq2 counts 
  * 12.) Run multiclass logistic regression

## Notes
  * I is performed per plate, so we will use DP1 as an example. All other plates can be processed similarly by replacing DP1 with the desired plate.

## Details
##### 1.) Create directories, symbolic link to fastqs, and distribute scripts for alignment
    description: create directory structure for read alignment, copy alignment/QC scripts
                 into relevant directories, & symlink to relevant .fastqs.
    script: /alignemnt_pipeline/util_Alignment_setup.sh
    dependencies: fastq files, named by plate number such as DP1_<things here>L1_R1.fastq 
    in: A plate number, such as DP1
    out: null

##### 2.) Align reads, merge across barcodes, & do read counts 
    description: align fastq files, merge barcodes, quality filter, then remove 
                 duplicates. Do read counts at each stage of filtering
    script: /alignemnt_pipeline/bams/Makefile 
    dependencies: the folder ../fastqs with symlinked fastqs.
    in: null
    out: clean bam files & .txt files with read counts from each stage of filtering

##### 3.) Combine read counts from all stages of filtering
    description: creates a histogram & table from read counts at each stage of .bam processing
    script: /alignemnt_pipeline/counts_QC_logs.sh 
    scripts: counts_QC_logs.R  
    dependencies: 2.) 
    in: A plate number, such as DP1
    out: DP1_QC_counts.pdf & DP1_QC_counts.txt

##### 4.) Gene counts using DESeq2
    description: use DESeq2 to measure gene expression, at the transcript level, from processed bams
    script: /alignemnt_pipeline/counts_DEG.sh 
    script1: counts_DEG.R
    dependencies: 2.) 
    in: A plate number, such as DP1
    out: DP1.data.gz & DP1.data.Rd

##### 5.) Make pileups for QuASAR
    description: make pileup/bed files from each clean bam file
    script: /alignemnt_pipeline/Makefile1
    dependencies: 2.) 
    in: null
    out: clean DP1*.pileup.clean.bed.gz files 

##### 6.) Create directories, symbolic link to pileups, and distribute scripts for QuASAR
    description: create directory structure for QuASAR pipeline, copy scripts
                 into relevant directories, & symlink to relevant *.clean.bed.gz files.
    script: /QuASAR_pipeline/util_QuASAR_setup.sh
    dependencies: 5.)	
    in: A plate number, such as DP1
    out: null

##### 7.) Run the QuASAR pipeline for joint genotyping ASE inference
    description: a full QuASAR analysis for each plate/cellLine
    script: /QuASAR_pipeline/util_QuASAR_runAll.sh (run all plates)
    script1: /QuASAR_pipeline/QuASAR_pipeline_all.sh (run all cell lines)
    script2: /QuASAR_pipeline/QuASAR_pipeline.sh (run an individual experiment)
    script3: /QuASAR_pipeline/QuASAR_pipeline.R
    script4: /QuASAR_pipeline/QuASAR_pipeline_masterTable.R
    dependencies: 6.) & a properly formatted covariate file  
    in: A plate number, such as DP1, & the name of the analysis script, QuASAR_pipeline.R
    out: QuASAR output for each plate:cellLine:treatment, plate:cellLine, QQ plots, and manhattan plots
    notes: gene annotations fro all SNPs are done at the QuASAR step resulting in *_newLogFC.txt.tid

##### 8.) Create directories, setup the master table, and distribute scripts for MESH
    description: bins data based on logFC and deposits analysis data into correct MESH directory for analysis
    script: /MESH_pipeline/util_MESH_setup.sh
    dependencies: 
    in: master table
    out: data is split and added to the MESH analysis directory


##### 9.) Combine QuASAR data into a master table and process data for MESH
    description: Create a master table with MESH input (betas & SEs) & MESH output (bayes-factors & posteriors) 
    script: /MESH_pipeline/MESH_QuASAR_master_pipeline.sh
    script1: /MESH_pipeline/MESH_QuASAR_master.sh
    script2: /MESH_pipeline/MESH_QuASAR_master_assignControls.R
    dependencies: 7.) 8.) 
    in: null
    out: Master_table_betas_bfs.txt

##### 10.) Split master table data on logFC for MESH analysis
    description: Split master table on the level of logFC output to MESH analysis folder
    script: /MESH_pipeline/MESH_QuASAR_master_split.sh
    dependencies: 9.) 
    in: Master_table_betas_bfs.txt
    out: an array of plots concerning condition specific ASE (biplots, dotplots, etc.)

##### 11.) Run the MESH pipeline for condition specific ASE
    description: Run MESH analysis, calculate BFs/posteriors, & plot results
    script: /MESH_pipeline/MESH/run_MESH_batch.sh
    script1: /MESH_pipeline/MESH/run_analysis.py
    script2: /MESH_pipeline/MESH/analysis.pl
    script3: /MESH_pipeline/MESH/calc_posteriors.R
    dependencies: 10.) 
    in: Master_table_betas_bfs.txt
    out: an array of plots concerning condition specific ASE (biplots, dotplots, etc.)

##### 12. Run multiclass logistic regression
    description: using the posteriors from MESH as labels, and plate, treatment, fpkm, gene expression, 
                 ... etc. as variables, run a multinomial logistic regression and plot the results (forest plots)
    script: multiclasslr_pipeline/run_multiclasslr_batch.sh
    script1: multiclasslr_pipeline/multiclasslr_batch.R
    dependencies: 11.) 
    in: data= masterTable_multinomial.txt (from 12.)) & a tag to lable the run (may include parameters 
        used for instance)
    out: summary data from the model run and plots in the directory ./multinomialPlots
    notes: optional model with combined cf1 & cf2