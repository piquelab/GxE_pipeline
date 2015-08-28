##################################################################
##
## System anlaysis of differential gene expression using DESeq2
##
## Created by CTH, modified by RPR & GMB
## Wayne State University, 2015
##################################################################

library(DESeq2)
library(qvalue)
require(ggplot2)
require(plyr)
library(reshape)
require(parallel)
require(BiocParallel)

## Required argument: plate prefix 
cargs <- commandArgs(trail=TRUE);
stopifnot(length(cargs)>=1)
platePrefix <- cargs[1]

## Optional argument: number of cores
cores <- as.integer(Sys.getenv("NCPUS"))
if(length(cargs)>=2)
  cores <- as.numeric(cargs[2])

timestamp()
cat(platePrefix,cores,"\n")

if(cores<1){cores <- 1}
ParallelSapply <- function(...,mc.cores=cores){
  simplify2array(mclapply(...,mc.cores=mc.cores))
}

## To run DESeq2 in parallel, using the
## BiocParallel library
register(MulticoreParam(cores))

LPG <- Sys.getenv("LPG")

## Gene counts: this is our data for anlaysis
readCounts <- paste('../../derived_data/', platePrefix, '/counts/GC/',
                    platePrefix, '.data.Rd', sep='')
load(readCounts)
n.barcodes <- dim(data)[2]

## master list of treatment names and IDs
treatmentKey <- read.table(paste('../treatmentKey.txt', sep=''),
                           as.is=TRUE, header=TRUE)
row.names(treatmentKey) <- treatmentKey$Treatment_ID

##################################################################
## assign variables, load data, and load experiment information
topDirectory <- 'out_data_'
outDir <- paste(topDirectory, platePrefix, sep='')
system(paste("mkdir -p",outDir))
##
plotsDir <- paste(outDir, '/plots', sep='')
system(paste("mkdir -p", plotsDir))
##
statsDir <- paste(outDir, '/stats', sep='')
system(paste("mkdir -p", statsDir))
##
dataDir <- paste(outDir, '/data_objects', sep='')
system(paste("mkdir -p", dataDir))
##
degDir <- paste(outDir, '/data_DEG', sep='')
system(paste("mkdir -p", degDir))
##
qcDir <- paste(outDir, '/QC', sep='')
system(paste("mkdir -p", qcDir))

## reference transcriptome
bedtranscript <- "/wsu/home/groups/piquelab/data/RefTranscriptome/ensGene.hg19.v2.bed.gz" 

############ QC ############
myNormedData <- ParallelSapply(seq_along(1:n.barcodes), FUN=function(ii){
      qqnorm(rank(data[, ii], ties.method = "random"), plot = F)$x
  })
myCor <- cor(myNormedData, method='pearson')
medianCor <- apply(myCor,1,median)
medmedCor <- median(medianCor)

## median absolute deviation
robust.sd <- median(abs(medianCor-medmedCor))*1.4826
thresh2 <- medmedCor-robust.sd*5

thresh <- .2
idxBlackListed <- medianCor<thresh
cat("##Num blacklist:",sum(idxBlackListed),"\n")

hist.name <- paste(qcDir, '/', platePrefix, "_QC_corr_hist", ".pdf", sep="")
plotCor <- data.frame(correlations=c(myCor))
pdf(hist.name)
m <- ggplot(plotCor, aes(x=correlations))
m + geom_histogram() + ggtitle(paste("Library correlations  ", platePrefix, sep='')) + geom_vline(xintercept = thresh2, colour="red", linetype="longdash")
dev.off()

blackBarcodes <- which(idxBlackListed)
if(length(blackBarcodes > 0)){
  fullBarcodes <- data.frame(blackListedBC=paste('HT', blackBarcodes, sep=''))
  bcName <- paste(qcDir, '/', platePrefix, "_QC_blackList", ".txt", sep="")
  write.table(fullBarcodes, file=bcName, quote=FALSE, row.names=FALSE, col.names=FALSE)
}

heat.name <- paste(qcDir, '/', platePrefix, "_QC_corr_heat", ".pdf", sep="")
pdf(heat.name)
image(x=1:n.barcodes, y=1:n.barcodes, myCor, axes=FALSE,xlab="",ylab="")
axis(1,at=c(1,(1:12)*8),las=1,cex=0.5,lwd.ticks=2)
axis(1,at=1:n.barcodes,label=NA,las=1,cex=0.5)
axis(2,at=c(1,(1:12)*8),las=1,cex=0.5,lwd.ticks=2)
axis(2,at=1:n.barcodes,label=NA,las=1,cex=0.5)
title(main=paste("Quantile Normalized Correlation ", platePrefix, sep=''))
dev.off()
############ /QC ############

## covariates file with experimental information on the samples
cov.file <- paste('../../derived_data/covariates/GxE_', platePrefix, '_covariates.txt', sep='')
cv <- read.table(cov.file , as.is=T, sep="\t", header=T, comment="")
cv <- cv[order(cv$Barcode.ID),]
cv <- cv[grep(platePrefix,cv$Plate.ID),] ## remove information not corresponding to $platePrefix
cv$BlackList=FALSE
cv$BlackList[idxBlackListed] <- TRUE
 
##################################################################
## anno is an object in P*.data.Rd 
## Naming the rows with the transcript ID. 
rownames(anno) <- anno$t.id

## Annotating transcript length in bp and GC content. 
gcContentFile <- "/wsu/home/groups/piquelab/data/RefTranscriptome/ensGene.hg19.faCount.gz"
anno2 <- read.table(gcContentFile,as.is=T,sep="\t",header=T,comment="")
rownames(anno2) <- gsub("hg19_ensGene_", "", anno2$X.seq)
anno2$avg.cg <- (anno2$C+anno2$G)/anno2$len
#coding length 
anno$codLen <- anno2[anno$t.id,"len"]
#transcript length
anno$txLen <- (anno$c.start-anno$start)+(anno$stop-anno$c.stop)+anno$codLen
anno2$avg.cg <- (anno2$C+anno2$G)/anno2$len
#Avg. CG on the coding part. 
anno$avg.cg <- anno2[anno$t.id,"avg.cg"]
rm(anno2)

##################################################################
## Manual conversion to R factor objects:
cv$Treatment.ID <- factor(cv$Treatment.ID)
TreatmentLevels <- levels(cv$Treatment.ID)
cat("#",TreatmentLevels,"\n")
cv$Control.ID <- factor(cv$Control.ID)
ControlLevels <- levels(cv$Control.ID)
cat("#",ControlLevels,"\n")
TreatmentOnlyLevels <- TreatmentLevels[!(TreatmentLevels %in% ControlLevels)]
cat("#",TreatmentOnlyLevels,"\n")
ControlOnlyLevels <- TreatmentLevels[(TreatmentLevels %in% ControlLevels)]
cat("#",ControlOnlyLevels,"\n")

cv$CellType <- factor(cv$CellType)
cv$CellLine <- factor(cv$CellLine)
CellLineLevels <- levels(cv$CellLine)
cat("#",CellLineLevels,"\n")

## Just print some tables to see that cv looks fine:
table(cv$Control.ID,cv$Treatment.ID)
table(cv$Treatment.ID,cv$CellLine)

##################################################################
## Preparing data for DEseq:
## Combine processed data into a DESeqDataSet
## & remove genes with very low coverage/expression 

allColSamples <- paste(cv$CellLine, cv$Treatment.ID, sep=".")
cat("#",allColSamples,"\n")
ddsFull <- DESeqDataSetFromMatrix(
    countData = round(data),
    colData = cv,
    design = ~ CellLine + Treatment.ID)
keep <- rowSums(counts(ddsFull)) > 20
ddsFull <- ddsFull[keep,]

## Fit the model on the whole plate
system.time(ddsFull <- DESeq(ddsFull,parallel=TRUE))

## Contrast each treatment to its respective control
res <- ParallelSapply(TreatmentOnlyLevels,function(t){

  ## Select appropiate control for this treatment
  c <- ControlLevels[cv$Control.ID[cv$Treatment.ID==t][1]]
  cat("-----------------------------------------------------------","\n")
  cat("Processing treatment:",t,treatmentKey[t,"Treatment_Name"],"\n")
  cat("Using control:",c,treatmentKey[c,"Treatment_Name"],"\n")            

  res <- results(ddsFull, contrast=c("Treatment.ID",t,c))
  head(res)

  ## Use a filter to help with power
  fThresh <- attr(res,"filterThreshold")
  use <- res$baseMean > fThresh

  ## Plot filter
  fname=paste(plotsDir, '/', platePrefix, "_", gsub(' ', '_', t),"_Filter", ".pdf", sep="")
  pdf(fname);
  aux = attr(res,"filterNumRej")
  plot(aux,type="b",ylab="number of rejections")
  abline(v=aux$theta[which.max(aux$numRej)],lty=3)  
  dev.off();

  ## QQ-plot
  fname=paste(plotsDir, '/', platePrefix, "_", gsub(' ', '_', t), ".pdf", sep="")
  pdf(fname, height=500, width=500)
  qqplot(-log10(ppoints(length(res$pval))),-log10(res$pval),pch=20,cex=0.5)
  title(main=paste0(t,", FDR10%=",sum(res$padj<0.1,na.rm=T),", ",treatmentKey[t,"Short_Name"],", ",sum(use),"/",length(use)))
  abline(0, 1, col='red')
  dev.off()

  ## Get the top differentially expressed genes
  fname=paste(statsDir, '/', platePrefix,"_","DEG_stats_", gsub(' ', '_', t), ".txt",sep="")
  sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
  names(sub.table) <- c('t.id', 'padj', 'pval', 'logFC')
  sub.table$ensg <- anno[sub.table$t.id, 'ensg']
  sub.table$g.id <- anno[sub.table$t.id, 'g.id']
  sub.table <- sub.table[!is.na(sub.table$padj), ]
  write.table(sub.table, file=fname, quote=FALSE, row.names=FALSE)
  cat("BH diff. expressed trx.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
  cat("BH unique genes diffexpr.", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), length(unique(anno$g.id[keep][(res$padj)<tr])) )),"\n")

  ## Return table with p-values, fold diff, qvalue and so on.  
  res 
}); names(res) <- TreatmentOnlyLevels

save(res,ddsFull, file=paste(dataDir, '/DESeq2_', platePrefix, '.RData', sep=''))

##############################
## Plots and summary tables ##
##############################
fname <- paste(plotsDir, '/', platePrefix, "_plotDispEsts.pdf", sep="")
pdf(fname)
plotDispEsts(ddsFull)
dev.off()

n.treats <- length(TreatmentOnlyLevels)
degs <- t(ParallelSapply(1:n.treats, FUN=function(this){
          ParallelSapply(c(0.01, 0.05, 0.1, 0.2), FUN=function(tr){
            DE_transcripts <- which(na.omit(res[[this]]$padj)<tr)
            DE_genes <- unique(anno[rownames(res[[this]])[DE_transcripts], 'ensg'])
            if(tr==0.1 & length(DE_genes)>0){
              outFile <- paste(degDir, '/', platePrefix, '_DEG_DE10_', TreatmentOnlyLevels[this], '.txt', sep='')
              write.table(DE_genes, file=outFile, quote=FALSE, row.names=FALSE, col.names=FALSE)
            } 
            length(DE_genes)
          })
}))

Treatment <- sapply(1:n.treats, function(ii){
      index <- which(treatmentKey$Treatment_ID==TreatmentOnlyLevels[ii])
      treatmentKey[index, "Short_Name"]
  })

degs <- data.frame(Treatment.ID=TreatmentOnlyLevels, Treatment, degs)
colnames(degs) <- c('Treatment.ID', 'Treatment', '0.01', '0.05', '0.1', '0.2')  
outFile <- paste(degDir, '/', platePrefix, '_DEG.txt', sep='')
write.table(degs, file=outFile, quote=FALSE, sep='\t', row.names=FALSE)

## bubble plot ##
countData <- paste('../../derived_data/', platePrefix, '/counts/QC/', platePrefix, '_QC_counts.txt', sep='')
cleanReads <- read.table(countData, header=TRUE, as.is=TRUE)
cleanReads <- cleanReads[, c('Barcode.ID', 'clean')]
names(cleanReads) <- c('Barcode.ID', 'readCount')
cleanReads$Barcode.ID <- as.numeric(sub('_.*','',sub(".*HT","",cleanReads$Barcode.ID)))

## merge read counts with covariate table so we have the treatment names
tempDat <- merge(cv[, c('Barcode.ID', 'Treatment.ID')], cleanReads, by='Barcode.ID')
counts <- ddply(tempDat, .(Treatment.ID), function(df){sum(df[, 3])})
counts$Treatment.ID <- as.character(counts$Treatment.ID)
counts <- counts[!(counts$Treatment.ID %in% c('CTRL (DMSO)', 'CTRL (ETOH)', 'CTRL (Water)')), ]
names(counts) <- c('Treatment.ID', 'CleanReads')

plotDat <- data.frame(merge(degs[, c('Treatment.ID', '0.1')], counts, by='Treatment.ID'), stringsAsFactors=FALSE)
plotDat$Treatment <- sapply(seq_along(1:dim(plotDat)[1]), FUN=function(ii){
    index <- which(rownames(treatmentKey)==plotDat$Treatment.ID[ii])
    treatmentKey$Short_Name[index]
  })
names(plotDat) <- c('Treatment.ID', 'hits', 'CleanReads', 'Treatment')
plotDat$hits[which(plotDat$hits==0)] <- 1
plotDat$log10.hits <- log10(plotDat$hits)

plotDatOrdered <- plotDat
plotDatOrdered$Treatment <- factor(plotDat$Treatment, levels = plotDat[order(plotDat$log10.hits), "Treatment"])

## standard bubbble plot
fname <- paste(plotsDir, '/', platePrefix, "_depth_power.pdf", sep="")
pdf(fname)
myPlot <- ggplot(plotDat, aes(x=Treatment, y=log10.hits, fill=Treatment), show_guide=FALSE)
myPlot + coord_flip() + geom_bar(stat='identity', width=.2, show_guide=FALSE) + geom_point(aes(size=CleanReads, colour=Treatment), show_guide=TRUE) + ggtitle(paste(platePrefix, ' Depth & Empirical Power (Q<.1)', sep='')) + scale_y_continuous(limits = c(0, 4.5)) + guides(colour="none")
dev.off()

## ordered bubble plot
fname <- paste(plotsDir, '/', platePrefix, "_depth_power_ordered.pdf", sep="")
pdf(fname)
myPlot <- ggplot(plotDatOrdered, aes(x=log10.hits, y=Treatment))
myPlot + geom_point(aes(size=CleanReads, colour=Treatment), show_guide=TRUE) + ggtitle(paste(platePrefix, ' Depth & Empirical Power (Q<.1)', sep=''))  + scale_x_continuous(limits = c(0, 4.5)) + guides(colour="none")
dev.off()

##
## THE END ##
##
