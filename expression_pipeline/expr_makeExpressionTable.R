#########################################################################
##
## GMB 09/09/15
## Combine gene expression (e.g., FPKM) and differential expression
## into a master table of data per transcript (or per gene)
##
## INPUT/ARGS:
##   platePrefix - plate to analyzie (e.g., DP1)
##   type - 'all' transcripts, or one 'perGene'
##   word - middle of file name, indicates aligner,
##          normalization, adjustments, etc
##
#########################################################################

library(reshape)

cargs<-commandArgs(trail=TRUE);
if (length(cargs)>=1) { platePrefix <- cargs[1] }
if (length(cargs)>=2) { word        <- cargs[2] }
if (length(cargs)>=3) { type        <- cargs[3] }

## Load the treatment and plate keys
treatmentKey <- read.table('../treatmentKey.txt', comment.char="", sep='\t',
                           as.is=TRUE, header=TRUE)
plateKey <- read.table('../plateKey.txt', as.is=T, header=TRUE)

## Load the gene expression data
geneSetFile = paste0(platePrefix, '/', platePrefix, '.', word, '.meanExpr.',
  type, '.gz')
geneSet <- read.table(as.is=T, sep='\t', header=T, file=geneSetFile)
rownames(geneSet) <- geneSet$t.id
treatments <- colnames(geneSet)[grep('T.*C.*', colnames(geneSet))]

## Annotate to the logFC for that Tx in each tx-control pair
lfc_dir = '../DESeq2/'
load(paste0(lfc_dir, '/out_data_', platePrefix, '/data_objects/',
            'DESeq2_', platePrefix, '.RData'))

## Load in the platewide average fpkm values
load(paste0(platePrefix, '/', platePrefix, '.', word, '.topTx.Rd'))

## Initialize the master table
tx = treatments[1]
ctrl <- treatmentKey[tx, ]$Control.ID
mTable <-data.frame(t.id = geneSet$t.id,
                    treatment = rep(tx, dim(geneSet)[1]),
                    control = rep(ctrl, dim(geneSet)[1]),
                    fpkm.t = geneSet[,tx],
                    fpkm.c = geneSet[,ctrl],
                    fpkm.avg=topTx[rownames(geneSet)],
                    stringsAsFactors=F)
lfc <- as.data.frame(res[tx])
lfc <- lfc[mTable$t.id, ]
mTable$lfc <- lfc[ , paste0(tx,'.log2FoldChange')]
mTable$lfcSE <- lfc[ , paste0(tx,'.lfcSE')]
mTable$lfcPval <- lfc[ , paste0(tx,'.pvalue')]
mTable$lfcPadj <- lfc[ , paste0(tx,'.padj')]

## Iteratively add data for each treatment
for (tx in treatments[2:length(treatments)]) {
  ctrl <- treatmentKey[tx, ]$Control.ID
  expr <-data.frame(t.id=geneSet$t.id,
                    treatment=rep(tx, dim(geneSet)[1]),
                    control=rep(ctrl, dim(geneSet)[1]),
                    fpkm.t=geneSet[,tx],
                    fpkm.c=geneSet[,ctrl],
                    fpkm.avg=topTx[rownames(geneSet)], stringsAsFactors=F)
  lfc <- as.data.frame(res[tx])
  lfc <- lfc[expr$t.id, ]
  expr$lfc <- lfc[ , paste0(tx,'.log2FoldChange')]
  expr$lfcSE <- lfc[ , paste0(tx,'.lfcSE')]
  expr$lfcPval <- lfc[ , paste0(tx,'.pvalue')]
  expr$lfcPadj <- lfc[ , paste0(tx,'.padj')]
  mTable <- rbind(mTable, expr)
}

## Sanity checks
summary(factor(mTable$treatment))
summary(factor(mTable$control))

## Annotate plate and cell type
mTable$plate <- platePrefix
mTable$cell.type <- unique(plateKey[ plateKey$plate == gsub('D', '', platePrefix),
                                       "cell.type" ])

## Reorder the fields
mTable <- mTable[ ,c("t.id", "plate", "cell.type", "treatment",
                     "control", "fpkm.t", "fpkm.c", "fpkm.avg", "lfc",
                     "lfcSE", "lfcPval", "lfcPadj")]

## Save two version of the master table file:
## 1) Drop all transcripts w/o logFC data
mTableNaRm <- mTable[!is.na(mTable$lfc), ]
outFile = paste0(platePrefix, '/', platePrefix, '.', word, '.',
  type, '.master.NaRm.txt.gz')
write.table(mTableNaRm, col.names=T, row.names=F, quote=F, sep='\t',
  file=gzfile(outFile))

## 2) Output everything, where lfc of NA -> 0
mTable$lfc[is.na(mTable$lfc)] <- 0
outFile = paste0(platePrefix, '/', platePrefix, '.', word, '.',
  type, '.master.Na0s.txt.gz')
write.table(mTable, col.names=T, row.names=F, quote=F, sep='\t',
	file=gzfile(outFile))
