
#### Covariate files

Covariate files are used to keep track of meta information associated with samples. These files are used in almost every stage of the pipeline, so it is important that they are properly formatted. Covariate files for deep sequencing plates are automatically derived from the shallow ones in the [QC counts](../alignment_pipeline/counts_QC_logs.R) script.

Here is an example file, [GxE_PX_covariates.txt](GxE_PX_covariates.txt), where 'X' would normally denote the plate number. The blanks fields in this example file are not required, but useful for record keeping. 

#### Gene annotation files

Several areas of the pipeline require a gene annotation file, a bed file of all transcripts in the transcriptome. This also may be necessary to build a genome index for alignment (e.g., with STAR).

##### Downloading the annotation data
The transcript annotations are avaialbe in the [UCSC Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables). To download it,
1) Select genome assembly
2) Select group "Genes and Gene Predictions"
3) Select track "Ensembl Genes"
4) Select table "ensGene"
5) Enter an output filename, then click "get output"
6) From the primary table, select the following: name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score, name2
7) From the "ensemblToGeneName" table, select name and value

##### Converting to a bed file
To convert the above download into the bed transcriptome file, use the provided [ensgeneToBed.py](ensgeneToBed.py) script.

Example:
```
[you@node ~]$ python ensgeneToBed.py <input> transcriptome.bed
```

##### Calculating transcript length and GC content
The transcript length and GC content are necessary for FPKM calculations. To obtain this information,
1) use the bed transcriptome file to extract fasta records for the transcripts
2) calculate length and GC content using the faCount tool

Example:
```
[you@node ~]$ bedtools getfasta -fi <genome>.fa -bed transcriptome.bed -fo  stdout -name -s -split | gzip > transcriptome.fa.gz
[you@node ~]$ faCount transcriptome.fa | gzip > transcriptome.faCount.gz
```

#### Checking alignments
As an optional step in the alignment pipeline, you can double check that the fastqs aligned properly. This step is meant to be run after the initial quality filtering, to also check that the merge and quality filtering finished without error.

To do so, run the [checkAlignment.sh](checkAlignment.sh) script as such:
```
[you@node misc]$ sh checkAlignment.sh -p <plate>
```
Optionally, you can use the -v (verbose) switch to print each barcode as it is checked.
