###########################################################################################
## System: Counts per transcript 
## CTH 09-Apr-14
##
## Prepatory step for DESeq2 Count reads for each transcript.
## 
###########################################################################################
require(parallel)

#cores <- 2
#platePrefix <- "DP5"
#dataFolder <- "../../bams"

cargs <- commandArgs(trail=TRUE)
if(length(cargs)>=1)
  cores <- cargs[1]
if(length(cargs)>=2)
  platePrefix <- cargs[2]
if(length(cargs)>=3)
  dataFolder <- cargs[3]


ParallelSapply <- function(...,mc.cores=cores){
  simplify2array(mclapply(...,mc.cores=mc.cores))
}

## Genome index was sorted by size, make sure the transcriptome file is too
bedtranscript <- "/wsu/home/groups/piquelab/data/RefTranscriptome/ensGene.hg19.2014.bed.gz";
anno <- read.table(bedtranscript,as.is=T,sep="\t")
anno <- anno[,-c(9:12)]
colnames(anno) <-  c("chr","start","stop","t.id","score","strand","c.start","c.stop","ensg","g.id")

## get the barcodes
barcodes <- as.integer(scan(pipe(paste('ls ../../bams/', platePrefix, '-HT*_clean.bam | sed s/_.*//g | sed s/.*', platePrefix, '-HT//g | sort | uniq', sep='')), character(0)))
barcodes <- sort(barcodes)

## Note: This only works with bedtools 2.24.0+
data <- ParallelSapply(barcodes, function(ii) {
  command <- paste0("coverageBed -counts -sorted -split -s -a ", bedtranscript,
                    " -b ", dataFolder, "/", platePrefix, "-HT", ii, "_clean.bam | cut -f15")
  scan(pipe(command))
})

rownames(data) <- anno$t.id;
colnames(data) <- paste(platePrefix,"-HT",barcodes,sep="")


dataFile=paste(platePrefix,".data.gz",sep="")
write.table(data,gzfile(dataFile),sep="\t",quote=F)

RdataFile=paste(platePrefix,".data.Rd",sep="")
save(anno,data,file=RdataFile)

##
## THE END
##
