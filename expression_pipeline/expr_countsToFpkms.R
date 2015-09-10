#########################################################################
##
## GMB 09/09/15
## Convert read counts to FPKM measurements. Uses a gene model file
## such as one from gencode or ensembl 
##
## INPUT/ARGS:
##   platePrefix - plate to analyzie (e.g., DP1)
##   bamMethod - aligner used, part of input filename
##   cores - number of cores for parallel sapply
##   
#########################################################################

## Libraries ##
library(parallel)
LPG <- Sys.getenv("LPG")

## Get command-line arguments ##
cargs<-commandArgs(trail=TRUE);
if (length(cargs)>=1) { platePrefix <- cargs[1] }
if (length(cargs)>=2) { bamMethod   <- cargs[2] }
if (length(cargs)>=3) { cores       <- cargs[3] }

## Create the output directory
system(paste0('mkdir -p ', platePrefix))

## Get base plate ID for deep plates
plateBase <- gsub('DP', 'P', platePrefix) 

if (cores < 1){cores <- 1}
ParallelSapply <- function(...,mc.cores=cores){
  simplify2array(mclapply(...,mc.cores=mc.cores))
}

## Get the FPKM data
readCounts = paste0('../../derived_data/', platePrefix, '/counts/GC/',
                    platePrefix, '.data.gz')
countsData <- read.table(readCounts, as.is=T, sep='\t', header=T)
barcodes <- gsub(paste0(platePrefix, '.HT'), '', colnames(countsData))

## Get the covariate file for the indicated plate
## If analyzing a deep plate, remove unused barcodes
cv <- read.table(paste0('../../derived_data/covariates/GxE_', plateBase,
                        '_covariates.txt'), as.is=T, sep='\t', header=T)
cv <- cv[ cv$Barcode.ID %in% barcodes, ]

## Load the reference transcriptome (ensembl)
bedtranscript <- paste0("less ", LPG,
  "/data/RefTranscriptome/ensGene.hg19.v2.bed.gz", " | cut -f 1-4,7,8,13,14")
transcriptAnno <- read.table(file=pipe(bedtranscript),as.is=T,sep="\t")
colnames(transcriptAnno) <-  c("chr", "start", "stop", "t.id", "c.start",
                               "c.stop", "ensg", "g.id")
rownames(transcriptAnno) <- transcriptAnno$t.id

## Annotate transcript length and average GC content
gcContentFile <- paste0(LPG, "/data/RefTranscriptome/ensGene.hg19.faCount.gz")
gcAnno <- read.table(gcContentFile,as.is=T, sep="\t", header=T, comment="")
rownames(gcAnno) <- gsub('hg19_ensGene_', '', gcAnno$X.seq)
gcAnno <- gcAnno[1:dim(gcAnno)[1]-1,] # trim the "totals" column
gcAnno$avg.cg <- (gcAnno$C + gcAnno$G) / gcAnno$len
transcriptAnno$txLen <- gcAnno[transcriptAnno$t.id, "len"]
transcriptAnno$avg.cg <- gcAnno[transcriptAnno$t.id, "avg.cg"]
rm(gcAnno)

## Manual conversion to R factor objects - block from CTH
cv$Treatment.ID <- factor(cv$Treatment.ID)
TreatmentLevels <- levels(cv$Treatment.ID)
cv$Control.ID <- factor(cv$Control.ID)
ControlLevels <- levels(cv$Control.ID)
TreatmentOnlyLevels <- TreatmentLevels[!(TreatmentLevels %in% ControlLevels)]
ControlOnlyLevels <- TreatmentLevels[(TreatmentLevels %in% ControlLevels)]
cv$CellType <- factor(cv$CellType)
cv$CellLine <- factor(cv$CellLine)
CellLineLevels <- levels(cv$CellLine)

##  Collapsing technical replicates
allColSamples <- paste(cv$CellLine, cv$Treatment.ID, sep=".")
sp <- split( seq(along=allColSamples), allColSamples )
cdata <- ParallelSapply(sp, function(columns)
                        round(rowSums( countsData[,columns,drop=FALSE] )))
cv2 <- cv[sapply(sp, `[` , 1),]

## Convert the count data to FPKMs
cs <- colSums(cdata)/1E6
cdata2 <- ParallelSapply(1:ncol(cdata), function(jj){
	cdata[,jj] / cs[jj] * 1000 / transcriptAnno$txLen
})
colnames(cdata2) <-  colnames(cdata)
geneAnno2 <- transcriptAnno[rownames(cdata2), ]

## Remove the lowly expressed transcripts
keep  <- rowMeans(cdata2 > 0.1) > 0.9

## Save the fpkm values
fpkms <- cdata2
save(list=c("cv2", "fpkms"),
     file=paste0(platePrefix, "/", platePrefix, ".", bamMethod, ".counts.fpkm.all.Rd"))

fpkms <- cdata2[ keep, ]
save(list=c("cv2", "fpkms"),
     file=paste0(platePrefix, "/", platePrefix, ".", bamMethod, ".counts.fpkm.Rd"))

txInd <- rownames(cdata2)
fit.lm <- lm(log10(cdata2[,1] + 1e-6) ~
             transcriptAnno[txInd,]$avg.cg + transcriptAnno[txInd,]$txLen +
             transcriptAnno[txInd,]$avg.cg * transcriptAnno[txInd,]$txLen)
summary(fit.lm)

cdata3 <- ParallelSapply(1:ncol(cdata2),function(ii){
  fit.lm <- lm(log10(cdata2[,ii] + 1e-6) ~
               transcriptAnno[txInd,]$avg.cg + transcriptAnno[txInd,]$txLen +
               transcriptAnno[txInd,]$avg.cg * transcriptAnno[txInd,]$txLen)
  x <- residuals(fit.lm)
})
colnames(cdata3) <- colnames(cdata2)

fpkms <- cdata3
save(list=c("cv2", "fpkms"),
     file=paste0(platePrefix, "/", platePrefix, ".", bamMethod, ".counts.fpkm.full.gc.Rd"))

fpkms <- cdata3[ keep, ]
save(list=c("cv2", "fpkms"),
     file=paste0(platePrefix, "/", platePrefix, ".", bamMethod, ".counts.fpkm.gc.Rd"))

## Quantile normalization (by sample)
cdata4 <- apply(cdata3,2,function(x){
  qqnorm(rank(x, ties.method = "random"), plot = F)$x
})
colnames(cdata4) <- colnames(cdata3)
rownames(cdata4) <- rownames(cdata3)

fpkms <- cdata4
save(list=c("cv2", "fpkms"),
     file=paste0(platePrefix, "/", platePrefix, ".", bamMethod, ".counts.fpkm.full.norm.Rd"))

fpkms <- cdata4[ keep, ]
save(list=c("cv2", "fpkms"),
     file=paste0(platePrefix, "/", platePrefix, ".", bamMethod, ".counts.fpkm.norm.Rd"))


## ADJUSTMENTS - two separate methods
## 1) Remove means
fpkms <- cdata4
for (cl in 1:length(unique(as.numeric(cv2$CellLine)))) {
  avg <- apply(fpkms[,as.numeric(cv2$CellLine)==cl], 1, mean)
  fpkms[,as.numeric(cv2$CellLine)==cl] <-
    fpkms[,as.numeric(cv2$CellLine)==cl] - avg
}
save(list=c("cv2", "fpkms"),
     file=paste0(platePrefix, "/", platePrefix, '.', bamMethod, ".counts.fpkm.full.adj.Rd"))

fpkms <- cdata4[ keep, ]
for (cl in 1:length(unique(as.numeric(cv2$CellLine)))) {
  avg <- apply(fpkms[,as.numeric(cv2$CellLine)==cl], 1, mean)
  fpkms[,as.numeric(cv2$CellLine)==cl] <-
    fpkms[,as.numeric(cv2$CellLine)==cl] - avg
}
save(list=c("cv2", "fpkms"),
     file=paste0(platePrefix, "/", platePrefix, '.', bamMethod, ".counts.fpkm.adj.Rd"))

## 2) Remove controls
fpkms <- cdata4
drop = c()
for (cl in levels(cv2$CellLine)) {
  for (cnt in levels(cv2$Control.ID)) {
    fpkms[,cv2$CellLine==cl & cv2$Control.ID==cnt] <-
      (fpkms[,cv2$CellLine==cl & cv2$Control.ID==cnt] -
       fpkms[,cv2$CellLine==cl & cv2$Treatment.ID==cnt])
    drop <- c(drop, which(cv2$CellLine==cl & cv2$Treatment.ID==cnt))
  }
}
fpkms <- fpkms[,-drop] # Drop the "empty" control column
save(list=c("cv2", "fpkms"),
     file=paste0(platePrefix, "/", platePrefix, '.', bamMethod, ".counts.fpkm.full.cntRm.Rd"))

fpkms <- cdata4[ keep, ]
drop = c()
for (cl in levels(cv2$CellLine)) {
  for (cnt in levels(cv2$Control.ID)) {
    fpkms[,cv2$CellLine==cl & cv2$Control.ID==cnt] <-
      (fpkms[,cv2$CellLine==cl & cv2$Control.ID==cnt] -
       fpkms[,cv2$CellLine==cl & cv2$Treatment.ID==cnt])
    drop <- c(drop, which(cv2$CellLine==cl & cv2$Treatment.ID==cnt))
  }
}
fpkms <- fpkms[,-drop] # Drop the "empty" control column
save(list=c("cv2", "fpkms"),
     file=paste0(platePrefix, "/", platePrefix, '.', bamMethod, ".counts.fpkm.cntRm.Rd"))

## the end ##

