
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
in: plate number, number of cores [1]
out: out_data_<plate>/, containing DEseq2 output & various plots
```

Example:
```
 $ Rscript DESeq2_automate.R DP1 12
```

##### Step 2: Calculate transcript FPKMs
Calculate expression values (FPKMs) from counts data. Saves output at multiple levels of post-processing, including raw FPKMs, normalized FPKMs, and FPKMs with individual effect removed.
```
script: expr_countsToFpkms.R
dependencies: gene counts, gene annotations (see misc)
in: plate number, aligner, cores
out: <plate>.<aligner>.counts.fpkm*.Rd
```

Example:
```
 $ Rscript expr_countsToFpkms.R DP1 bwa 12
```

##### Step 3: Select transcripts representative of a gene, for gene-based analyses
Collapse expression data across individuals, and report the higest expressed transcript per gene.

```
script: expr_getMeanAndTopExpr.R
dependencies: step 2
in: cell type, file indicator, cores
out: <plate>.<indicator>.meanExpr.[all | perGene].perGene and <plate>.<indicator>.topTx.Rd [per plate of given cell type]
```
Example:
```
 $ Rscript expr_getMeanAndTopExpr.R LCL bwa.counts.fpkm.adj 12
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

