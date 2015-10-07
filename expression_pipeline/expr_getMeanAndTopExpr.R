##################################################################
## 
## GMB 09/09/15
## Takes a matrix of FPKM data and a gene annotation file,
## and returns a matrix condensed to one transcript per gene
##
## Required Arguments
##   celltype - character; cell type to process, e.g., HUVEC
##   word - character; file infix denoting processing steps, etc
##
## Optional Arguments
##   cores - number of cores for parallel sapplys
##   bedTranscriptome - character; transcriptome bed file
##
##################################################################

## Libraries ##
library(parallel)
library(reshape)
source('../../GxE_pipeline/misc/getArgs.R')

## Get command-line arguments
defaultList = list(
  cores=1,
  bedTranscriptome="/wsu/home/groups/piquelab/data/RefTranscriptome/ensGene.hg19.2014.bed.gz"
  )
args <- getArgs(defaults=defaultList)
#platePrefix      <- args$platePrefix # change to cell.type
cell.type        <- args$celltype
word             <- args$word
cores            <- as.numeric(args$cores)
bedTranscriptome <- args$bedTranscriptome

## Check for required arguments
stopifnot(!is.null(cell.type) & !is.null(word))

print(args)

ParallelSapply <- function(...,mc.cores=cores){
  simplify2array(mclapply(...,mc.cores=mc.cores))
}

##################################################################
## specify covariate file and pileup directories
##################################################################
covDir <- '../../derived_data/covariates/'
cmd <- paste0('cat ', covDir, 'GxE_DP*_covariates.txt | grep -w ', cell.type)
cv <- read.table(pipe(cmd), as.is=TRUE, sep='\t')
names(cv) <- c("Plate.ID", "Barcode.ID", "Raw.Reads", "Quality.Reads", "Clean.Reads",
                              "Treatment.ID", "Treatment", "BarcodeSequence", "CellType", "CellLine",
                              "Control.ID", "Control")
cv$Treatment <- gsub(' ',  '_', cv$Treatment)


## Load the gene annontation file
anno <- read.table(as.is=T, sep="\t",
                   file=pipe(paste0("less ", bedTranscriptome, " | cut -f 1-4,7,8,13,14")))
colnames(anno) <-  c("chr", "start", "stop", "t.id", "c.start", "c.stop",
                     "ensg", "g.id")
rownames(anno) <- anno$t.id

## Data comes from the counts_to_rpkms.R script
## Assumes colnames are "<indiv>.<cond>"

##load(paste0(platePrefix, '/', platePrefix, '.', word, '.Rd'))
fpkms <- do.call('cbind', lapply(unique(cv$Plate.ID), function(p) {
  load(paste0(p, '/', p, '.', word, '.Rd'))
  colnames(fpkms) <- paste(p, colnames(fpkms), sep='.')
  fpkms
}))

## Get the mean expression, across transcripts, per condition
## NOTE: Breaks if '.' in cell.line ID (or txID or plateID)
x <- strsplit(colnames(fpkms), '\\.')
plates      <- unique(sapply(x, function(y) { y[1] }))
individuals <- unique(sapply(x, function(y) { y[2] }))
conditions  <- unique(sapply(x, function(y) { y[3] }))
meanExpr <- as.data.frame(ParallelSapply(conditions, function(c) {
  apply(fpkms[ , grep(c, colnames(fpkms))], 1, mean)
}))
colnames(meanExpr) <- conditions

## Add the gene annotation information
m <- match(rownames(meanExpr), anno$t.id)
meanExprAnno <- as.data.frame(cbind(anno[m , c(1,2,3,4,7,8)], meanExpr))

## Get the list of top expressed transcrips
topTx <- ParallelSapply(unique(meanExprAnno$ensg), function(g){
      aux <- meanExprAnno[ meanExprAnno$ensg==g, ]
      m = apply(aux[, 7:dim(aux)[2]], 1, mean)
      n = names(which( m == max(m) ))
      if (length(n)>1) {
        #return (sample(n, 1))
        print(n)
        n = sample(n, 1) ## If multiple have same avg, choose one randomly
      }
      m[n]
    })

## Ouput the mean gene expression for each transcript. Split back
## into separate files per plate
for ( p in plates ) {
  tx <- unique(cv$Treatment.ID[ cv$Plate.ID == p ])  
  outFile = paste0(p, "/", p, '.', word, '.meanExpr.all.gz')
  write.table(meanExprAnno[, c(colnames(meanExprAnno)[1:6], sort(tx))],
              quote=F, sep='\t', row.names=F, file=gzfile(outFile))
  outFile = paste0(p, '/', p, '.', word, '.meanExpr.perGene.gz')
  write.table(meanExprAnno[names(topTx), ], quote=F, sep='\t', row.names=F,
              file=gzfile(outFile))
  save(topTx, file=paste0(p, '/', p, '.', word, '.topTx.Rd'))
}

## the end ##
