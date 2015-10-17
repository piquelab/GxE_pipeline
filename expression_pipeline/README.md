
Expression pipeline
===================
This directory contains scripts related to the expression analysis of the GxE pipeline.

## Setup
Initialize the expression analysis directory using the utility script:
```
 $ sh util_expression_setup.sh 
```


## Prerequisites
* The [alignment pipeline](../alignment_pipeline/) must be run first
* A [covariate file](../misc/) explaining each barcde and linking treatments to controls
* A [gene annotation](../misc/) file with transcript information
* The [getArgs](../misc/) helper function to better process command line arguments in R

## Outline
1. Asses differential gene expression with DESeq2
2. Calculate transcript FPKMs
3. Select transcripts representative of a gene, for gene-based analyses
4. Combine expression data with the QuASAR output

## Notes
Each step is run on one plate at a time.

#### Definitions
The following terms are used, often as input to scripts or parts of file names.
<dl>
	<dt>Plate</td>
	<dd>Generally refers to a single sequencing run, comprised of three individual samples of a single cell type</dd>
	<dt>Aligner / 'bams'</dt>
	<dd>A short, file-name friendly string refering to the aligner used, e.g., 'bwa', 'tophat', or 'star'</dd>
	<dt>Word</dt>
	<dd>A string indicating processing steps, used to distinguish output from <a href=expr_countsToFpkms.R>expr_countsToFpkms.R</a></dd>
	<dt>Type</dt>
	<dd>Adjustment type for transcripts with NA logFC values, either remove ("NaRm") or 0'd out ("Na0s")</dd>
</dl>

## Express Pipeline
To quickly apply the pipeline to multiple samples, use the utility scripts provided. Below is an example workflow:
```
 $ sh util_expression_runAll.sh bwa counts.fpkm.adj perGene mmtxq 12
```
This will submit the expression pipeline script, [expr_run.R](./expr_run.R) for each sample.

---
## Full Pipeline
##### Step 1: Asses differential gene expression using DESeq2
Determine genes that are differentially expressed in a response to a treatment.

```
script: DESeq2_automate.R
dependencies: gene counts
in: plate number, aligner used, number of cores [1]
out: out_data_<plate>/, containing DEseq2 output & various plots
```

Required arguments
* Plate ID to process
* Number of cores to use [1]

Example:
```
[you@node DESeq2]$ Rscript DESeq2_automate.R DP1 12
```

Optionally, you can submit the job to the PBS scheduling system:
```
[you@grid DESeq2]$ sh DESeq2_automate.sh DP1
```

###### Output


##### Step 2: Calculate transcript FPKMs
Calculate expression values (FPKMs) from counts data. Saves output at multiple levels of post-processing, including raw FPKMs, normalized FPKMs, and FPKMs with individual effect removed.

```
script: expr_countsToFpkms.R
dependencies: gene counts, gene annotations (see misc)
```

Required arguments
* platePreifx - Plate ID to process, e.g., DP1
* bamMethod - Aligner used, e.g., bwa

Optional arguments
* cores - number of processors to use [1]
* rmLow - Switch to remove lowly expressed transcripts [false]
* bedTranscriptome - location of the bed-formatted transcriptome
* gcContentFile - location of the GC content information of transcriptome

Example:
```
[you@node fpkms]$ Rscript expr_countsToFpkms.R platePreix=DP1 bamMethod=bwa cores=12
```

To remove lowly expressed transcripts,
```
[you@node fpkms]$ Rscript expr_countsToFpkms.R platePreix=DP1 bamMethod=bwa cores=12 rmLow
```

##### Step 3: Select transcripts representative of a gene, for gene-based analyses
Collapse expression data across individuals, and report the higest expressed transcript per gene.

```
script: expr_getMeanAndTopExpr.R
dependencies: step 2
```

Required arguments
* celltype - The cell type to process
* word - the file infix denoting aligner used & processing steps

Optional arguments
* cores - number of processors to use
* bedTranscriptome - location of the bed-formatted transcriptome

Example:
```
 $ Rscript expr_getMeanAndTopExpr.R celltype=LCL word=bwa.counts.fpkm.adj cores=12
```

##### Optional: Combine Step 3 & 4 and submit to scheduler
To simplify steps 3 and 4, use the provided shell script to configure the inputs and schedule the job. The script will lookup plates beloning to an input cell type and run expr_countsToFpkms.R for each before running expr_getMeanAndTopExpr.R

```
script: expr_run.sh
dependencies: step 1
```

Required arguments (in order to be provided)
* Cell type to process
* Aligner used
* Adjustment method, e.g., which output from countToFpkms to use
  * Most commonly 'counts.fpkm.gc' (or 'counts.fpkm.rmlow.gc'), indicating the FPKMs have been normalized and adjusted for length & GC content biases

Example:
```
[you@grid fpkms]$ sh expr_run.sh LCL bwa counts.fpkm.gc 12
```

##### Optional: Generate PCA and Heatmap plots from expression values
Generate PCA and heatmap/clustering plots from expression data
```
script: expr_plots.R
dependences: step 1
in: expression file (.Rd)
out: pdf plots
```
Example:
```
 $ Rscript expr_plots.R DP1.counts.fpkm.adj.Rd
```

