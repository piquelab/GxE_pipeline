
Expression pipeline
===================

This directory contains scripts related to the expression analysis of the GxE pipeline.

## Prerequisites
* The [alignment pipeline](../alignment_pipeline/) must be run first
* A [covariate file](../misc/) explaining each barcde and linking treatments to controls
* A [gene annotation](../misc/) file with transcript information

## Outline
1) Asses differential gene expression with DESeq2
2) Calculate transcript FPKMs
3) Select transcripts representative of a gene, for gene-based analyses
4) Combine expression and differential expression into a master table
5) Combine expression data with the QuASAR output

## Notes
Each step is run on one plate at a time.

### Definitions
The following terms are used, often as input to scripts or parts of file names.
<dl>
	<dt>Plate</td>
	<dd>Generally refers to a single sequencing run, comprised of three individual samples of a single cell type</dd>
	<dt>Aligner / 'bams'</dt>
	<dd>A short, file-name friendly string refering to the aligner used, e.g., 'bwa', 'tophat', or 'star'</dd>
	<dt>Word</dt>
	<dd>A string indicating processing steps, used to distinguish output from [expr_countsToFpkms.R](./expr_countsToFpkms.R)</dd>
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
in: plate number, file indicator, cores
out: <plate>.<indicator>.meanExpr.[all | perGene].perGene and <plate>.<indicator>.topTx.Rd
```
Example:
```
 $ Rscript expr_getMeanAndTopExpr.R DP1 bwa.counts.fpkm.adj 12
```

##### Step 4: Combine expression and differential expression into a master table
Merge expression and differential expression data. Outputs two files with different approaches to transcripts with NA logFC values. Either remove them and drop transcripts ("NaRm"), or set them to 0 ("Na0s").

```
script: expr_makeExpressionTable.R
dependencies: step 1 & 3
in: plate number, file indicator, type (from step 3)
out: <plate>.<indicator>.<type>.master.[NaRm | Na0s]
```

Example:
```
 $ Rscript expr_makeExpressionTable.R DP1 bwa.counts.fpkm.adj all
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

