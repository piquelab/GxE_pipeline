##################################################################
## 
## GMB 09/09/15
## Takes a matrix of FPKM data and a gene annotation file,
## and returns a matrix condensed to one transcript per gene
##
## Required Arguments
##   platePrefix - character; experiment to run, e.g., DP1
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
platePrefix      <- args$platePrefix # change to cell.type
word             <- args$word
cores            <- as.numeric(args$cores)
bedTranscriptome <- args$bedTranscriptome

print(args)

ParallelSapply <- function(...,mc.cores=cores){
  simplify2array(mclapply(...,mc.cores=mc.cores))
}

## Load the gene annontation file
anno <- read.table(as.is=T, sep="\t",
                   file=pipe(paste0("less ", bedTranscriptome, " | cut -f 1-4,7,8,13,14")))
colnames(anno) <-  c("chr", "start", "stop", "t.id", "c.start", "c.stop",
                     "ensg", "g.id")
rownames(anno) <- anno$t.id

## Data comes from the counts_to_rpkms.R script
## Assumes colnames are "<indiv>.<cond>"
load(paste0(platePrefix, '/', platePrefix, '.', word, '.Rd'))

## Get the mean expression, across transcripts, per condition
individuals = unique(gsub("\\..*", "", colnames(fpkms)))
conditions = unique(gsub(".*\\.", "", colnames(fpkms)))
meanExpr <- as.data.frame(ParallelSapply(conditions, function(c) {
  apply(fpkms[ , grep(c, colnames(fpkms))], 1, mean)
}))
colnames(meanExpr) <- conditions

## Ouput the mean gene expression for each transcript
stopifnot(sum(rownames(meanExpr) == anno$t.id) == dim(meanExpr)[1])
meanExprAnno <- as.data.frame(cbind(anno[ , c(1,2,3,4,7,8)], meanExpr))
outFile = paste0(platePrefix, "/", platePrefix, '.', word, '.meanExpr.all.gz')
write.table(meanExprAnno, quote=F, sep='\t', row.names=F, file=gzfile(outFile))

## Get the list of top expressed transcrips
topTx <- ParallelSapply(unique(meanExprAnno$ensg), function(g){
      aux <- meanExprAnno[ meanExprAnno$ensg==g, ]
      m = apply(aux[, 7:dim(aux)[2]], 1, mean)
      n = names(which( m == max(m) ))
      if (length(n)>1) {
        #return (sample(n, 1)) 
        n = sample(n, 1) ## If multiple have same avg, choose one randomly
      }
      m[n]
    })

## Save the list of max transcripts and their platewide mean fpkm
save(topTx, file=paste0(platePrefix, '/', platePrefix, '.', word, '.topTx.Rd'))

## Output the data twice, once for all transcripts, once for the top set
outFile = paste0(platePrefix, '/', platePrefix, '.', word, '.meanExpr.all.gz')
write.table(meanExprAnno, row.names=F, quote=F, sep='\t',
            file=gzfile(outFile))
outFile = paste0(platePrefix, '/', platePrefix, '.', word, '.meanExpr.perGene.gz')
write.table(meanExprAnno[names(topTx), ], quote=F, sep='\t', row.names=F,
  file=gzfile(outFile))

## the end ##
