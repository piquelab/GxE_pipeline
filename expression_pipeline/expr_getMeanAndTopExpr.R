##################################################################
## 
## GMB 09/09/15
## Takes a matrix of FPKM data and a gene annotation file,
## and returns a matrix condensed to one transcript per gene
##
## INPUT/ARGS
##   platePrefix - experiment to run, e.g., DP1
##   cores - number of cores for parallel sapplys
##
##################################################################

library(parallel)
library(reshape)
LPG <- Sys.getenv("LPG")

cargs<-commandArgs(trail=TRUE);
if (length(cargs)>=1) { platePrefix <- cargs[1] }
if (length(cargs)>=2) { word        <- cargs[2] }
if (length(cargs)>=3) { cores       <- cargs[3] }

if (cores < 1){cores <- 1}
ParallelSapply <- function(...,mc.cores=cores){
  simplify2array(mclapply(...,mc.cores=mc.cores))
}

## Load the gene annontation file
bedtranscript <- paste0("less ", LPG,
  "/data/RefTranscriptome/ensGene.hg19.v2.bed.gz", " | cut -f 1-4,7,8,13,14")
anno <- read.table(file=pipe(bedtranscript),as.is=T,sep="\t")
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
