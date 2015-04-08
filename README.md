# GxE-pipeline
This repo contains scripts that comprise the Gene Environment Interaction (GxE) pipeline, for alignment & statistical modeling of allelic imbalances in RNA-Seq data across many cellular environments.

## Outline
### I. Align RNA-seq reads and perform QC
  * 1.) Create directories, symbolic link to fastqs, and distribute scripts for alignment
  * 2.) Align reads, merge barcodes, and perform QC
  * 3.) Read counts from all stages of filtering
  * 4.) Gene counts using DESeq2
  * 5.) Make pileups for QuASAR

### II. QuASAR (EM) pipeline to infer individual ASE
  * 6.) Create directories, symbolic link to pileups, and distribute scripts for QuASAR
  * 7.) Run the QuASAR pipeline for joint genotyping and inference   
  * 8.) Add gene annotations to QuASAR output

### III. Process QuASAR output and prepare for MESH 
  * 9.) Combine all QuASAR output into a master table
  * 10.) Add logFC annotations (Gregs script)
  * 11.) Split data for MESH and distribute to analysis directory

### IV. MESH (MCMC) to infer condition specific ASE
  * 12.) Run the MESH pipeline
  * 13.) Combine all MESH data
  * 14.) Plot condition specific ASE

### V. Multiclass logistic regression on MESH output and DESeq2 counts 
  * 15.) Run multiclass logistic regression

## Notes
  * Label the top directory as: td=/wsu/home/groups/piquelab/charvey/GxE
  * Scripts that call other scripts are chained together using ->
  * I is performed per plate, so we will use DP1 as an example. All other plates can be processed similarly by replacing DP1 with the desired plate.

## Details
##### 1.) Create directories, symbolic link to fastqs, and distribute scripts for alignment
    description: create directory structure for read alignment, copy alignment/QC scripts
                 into relevant directories, & symlink to relevant .fastqs.
    script: td/derived_data/scripts/Alignment_util_makelinks.sh 
    dependencies: fastq files, named by plate number such as DP1_<things here>L1_R1.fastq 
    in: A plate number, such as DP1
    out: null

##### 2.) Align reads, merge barcodes, and perform QC 
    description: align fastq files, merge barcodes, quality filter, then remove duplicates. Do read counts at each stage of filtering
    script: td/derived_data/DP1/bams/Makefile 
    dependencies: the folder ../fastqs with symlinked fastqs.
    in: null
    out: clean bam files & .txt files with read counts from each stage of filtering

##### 3.) Read counts from all stages of filtering
    description: creates a histogram & table from read counts at each stage of .bam processing
    script: td/derived_data/DP1/counts/QC counts_QC_logs.sh -> counts_QC_logs.R  
    dependencies: 2.) 
    in: A plate number, such as DP1
    out: DP1_QC_counts.pdf & DP1_QC_counts.txt

##### 4.) Gene counts using DESeq2
    description: use DESeq2 to measure gene expression, at the transcript level, from processed bams
    script: td/derived_data/DP1/counts/GC counts_DEG.sh -> counts_DEG.R
    dependencies: 2.) 
    in: A plate number, such as DP1
    out: DP1.data.gz & DP1.data.Rd

##### 5.) Make pileups for QuASAR
    description: make pileup/bed files from each clean bam file
    script: td/derived_data/DP1/pileups/Makefile1
    dependencies: 2.) 
    in: null
    out: clean DP1*.pileup.clean.bed.gz files 

##### 6.) Create directories, symbolic link to pileups, and distribute scripts for QuASAR
    description: create directory structure for QuASAR pipeline, copy QuASAR/MESH scripts
                 into relevant directories, & symlink to relevant *.clean.bed.gz files.
    script: td/jointGenotyping/scripts/QuASAR_util_makeLinks.sh
    dependencies: 5.) 
    in: A plate number, such as DP1
    out: null

##### 7.) Run the QuASAR pipeline for joint genotyping and inference
    description: a full QuASAR analysis for each plate/cellLine
    script: td/jointGenotyping/QuASAR_results_DP1/ QuASAR_pipeline_all.sh ->  QuASAR_pipeline.R
    dependencies: 6.) & a properly formatted covariate file in td/derived_data/covariates 
    in: A plate number, such as DP1, & the name of the analysis script, QuASAR_pipeline.R
    out: QuASAR output for each plate:cellLine:treatment, plate:cellLine, QQ plots, and manhattan plots

##### 8.) Add gene annotations to QuASAR output
    description: add annotations of each gene to every heterozygote SNP identified with QuASAR
    script: td/jointGenotyping/QuASAR_results_DP1/output/add_annotations.sh 
    dependencies: 7.) 
    in: null
    out: QuASAR output tables with gene name annotation on each SNP

##### 9. Combine all QuASAR output into a master table
    description: create a master table from QuASAR output for each plate:cellLine 
    script: td/jointGenotyping/QuASAR_results_DP1/ QuASAR_pipeline_all.sh ->  QuASAR_pipeline_masterTable.R
    dependencies: 7.) 8.) 
    in: A plate number, such as DP1, & the name of the analysis script, QuASAR_pipeline_masterTable.R
    out: an output file, DP1_<cellLine>_masterTable_logFC.txt, ready to be combined into a master masterTable

#### 9.a combine masterTables?
    description:  
    script: 
    dependencies: 
    in:
    out:

##### 10. Add logFC annotations (Greg's script)
    description: annotate logFC and fpkm data to each SNP in the masterTable  
    script: ??
    dependencies: 10.) is succesfully completed
    in: masterTable from 9.a
    out: masterTable annotated with logFC and fpkm data


##### 11. Split data for MESH and distribute to analysis directory
    description: bins data based on logFC and deposits analysis data into correct MESH directory for analysis
    script: QuASAR_assignControls.R (deprecated?)
    script1: QuASAR_MESH_split_logFC.sh 
    dependencies: 10.)
    in: master table
    out: data is split and added to the MESH analysis directory

##### 12. Run the MESH pipeline
    description: run MESH on the processed master table data, calculate posteriors, calculate bayes factors, plot data, & prep for multinomial.
    script: td/jointGenotyping/QuASAR_results_masterTable_4/data_MESH_all/run_MESH_batch.sh 
    script1: run_analysis.py
    script2: calc_posteriors.R
    dependencies: 11.) 
    in: the name of the data set to be anlaysed (should be in the same directory)
    out: MESH output, posteriors, bayes factors, configuration estimates, & input data for a multinomial analysis

##### 13. Combine all MESH data 
    description: Create a master table with MESH input (betas & SEs) & MESH output (bayes-factors & posteriors) 
    script: td/jointGenotyping/QuASAR_results_masterTable_4/data_MESH_all/Master_table_betas_bfs.R
    script1: master_MESH_outTable.sh || concatenate posterior CIs of configs and pi0
    dependencies: 12.) 
    in: null
    out: Master_table_betas_bfs.txt

##### 14. Plot condition specific ASE
    description: Plot condition specific ASE form MESH and QuASAR data
    script: td/jointGenotyping/QuASAR_results_masterTable_4/data_MESH_all/MESH_plots.R
    dependencies: 13.) 
    in: Master_table_betas_bfs.txt
    out: an array of plots concerning condition specific ASE (biplots, dotplots, etc.)

##### 15. Run multiclass logistic regression
    description: using the posteriors from MESH as labels, and plate, treatment, fpkm, gene espression, ... etc. as variables,
                 run a multinomial logistic regression and plot the results (forst plots)
    script: td/jointGenotyping/QuASAR_results_masterTable_4/data_MESH_all/multinomial/run_multiclasslr_batch.sh
    script1: multiclasslr_batch.R
    dependencies: 12.) 
    in: data= masterTable_multinomial.txt (from 12.)) & a tag to lable the run (may include parameters used for instance)
    out: summary data from the model run and plots in the directory ./multinomialPlots
