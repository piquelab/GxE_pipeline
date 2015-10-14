
Expression pipeline
===================
This directory contains scripts related to the joint genotyping & Allele-specific expression testing.

## Setup
Initialize the expression analysis directory using the utility script:
```
[you@grid ~] $ cd <base>/GxE_pipeline/QuASAR_pipeline/
[you@grid <base>] $ sh util_QuASAR_setup.sh
```
where <base> represents your project directory, where you cloned the GxE_pipeline.

This sets up the directory <base>/jointGenotyping/, as well as <base>/jointGenotyping/QuASAR_results_<cellType> for each <cellType> found in your covariate files.

## Prerequisites
* The [alignment pipeline](../alignment_pipeline/) must be run first
* A [covariate file](../misc/) explaining each barcde and linking treatments to controls
* A [gene annotation](../misc/) file with transcript information

## Outline
1. Prep the samples for genotyping
2. Run QuASAR

## Notes
Each step is run on one cell type (e.g., HUVECs) at a time, so be sure that all plates have been run through the alignment pipeline.

#### Definitions
The following terms are used, often as input to scripts or parts of file names.
<dl>
	<dt>Plate</td>
	<dd>Generally refers to a single sequencing run or group of shared runs, comprised of three individual samples of a single cell type</dd>
	<dt>Cell line / Individual</dt>
	<dd>An identifer tracing a sample back to the cell line or individual from which the sample was derived</dd>
</dl>

## Express Pipeline
To quickly apply the pipeline to multiple samples, use the utility script provided:
```
[you@grid jointGenotyping]$ sh QuASAR_pipeline_all.sh
```
This will read through the covariate files and automatically submit all cell types through QuASAR.

---
## Full Pipeline
##### Step 1: Prep the samples for genotyping
Process the filtered pileups into a format for QuASAR.

```
script: QuASAR_prep.sh / QuASAR_prep.R
dependencies: pileup files
in: cell type, number of cores [4]
out: QuASAR_results_<celltype>/data/*_quasarIn.Rd files
```

Example:
```
[you@grid jointGenotyping]$ sh QuASAR_prep.sh HUVEC
```

Alternatively, you can call the R script directly:
```
[you@node jointGenotyping]$ Rscript QuASAR_prep.R HUVEC 4
```

##### Step 2: Run QuASAR
Run QuASAR, first jointly genotyping each individual across all data for that individual, then testing for ASE within each individual-treatment experiment.

```
script: QuASAR_pipeline.sh / QuASAR_pipeline.R
dependencies: pileups, *_quasarIn.Rd files from Step 1
in: cell type, number of cores [4]
out: out_data_<cellType>/, containing genotypes and ASE test results
```

Example:
```
[you@grid jointGenotyping]$ sh QuASAR_pipeline.sh HUVEC
```

Alternatively, you can call the R script directly:
```
[you@node jointGenotyping]$ Rscript QuASAR_pipeline.sh HUVEC 4
```
