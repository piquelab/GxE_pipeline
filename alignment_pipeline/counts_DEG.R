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

bedtranscript <- "/wsu/home/groups/piquelab/data/RefTranscriptome/ensGene.hg19.v2.bed.gz";
anno <- read.table(bedtranscript,as.is=T,sep="\t")
anno <- anno[,-c(9:12)]
colnames(anno) <-  c("chr","start","stop","t.id","score","strand","c.start","c.stop","ensg","g.id")

## get the barcodes
barcodes <- as.integer(scan(pipe(paste('ls ../../fastqs/', platePrefix, '-HT* | sed s/_L.*//g | sed s/.*', platePrefix, '-HT//g | sort | uniq', sep='')), character(0)))

barcodes <- sort(barcodes)

data <- ParallelSapply(barcodes, function(ii){
	command <- paste("samtools bedcov ", bedtranscript, " ", dataFolder, "/", platePrefix, "-HT", ii, "_clean.bam | cut -f15",sep="") 
	cat("#",command,"\n");
	aux <- scan(pipe(command))
	aux / 100
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
