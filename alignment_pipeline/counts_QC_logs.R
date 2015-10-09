###########################################################################################
## System diagnostics
## CTH 4.8.2014
###########################################################################################
## all files made with commands such as
## ls *quality.bam | while read f; do echo $f `samtools view -c $f`; done > counts_quality.txt
require('ggplot2')
#library(ggplot2)

##################################################################
## takes a list of inference results and recursively merges them
## across common variables
##################################################################
MergeList <-  function(index, list, mergeCol){
  stopifnot(index >= 2)
  stopifnot(is.list(list))
  l.len <- length(list)
  stopifnot(l.len >= 2)
  cat(paste("Merging sample ", l.len, " \n",sep=""))
  if(l.len==2){
    return(merge(list[[2]][, c(1, index)], list[[1]][, c(1, index)], by = mergeCol, all=TRUE))
  }else{
    temp <- list[[l.len]]
    list[[l.len]] <- NULL
    merge(temp[, c(1, index)], Recall(index, list, mergeCol), by = mergeCol, all=TRUE)
  }
}

#x11(display="localhost:10.0" ,type="Xlib")

cargs <- commandArgs(trail=TRUE);
if(length(cargs)>=1)
  platePrefix <- cargs[1]
if(length(cargs)>=2)
  bamDir <- cargs[2]

## for debugging
#platePrefix <- "D1P6"
#bamDir <- "../../bams"

filters <- c('raw', 'merged', 'quality', 'clean')

## collect all of the read counts from log files
all_counts <- lapply(filters, function(this_filter){
  #this_filter <- filters[1]
  tmpFile <- scan(pipe("mktemp -t"),character(0))
  system(paste('ls ', bamDir, '/*', this_filter, '_count.txt | while read f; do cat $f | xargs echo $f; done > ', tmpFile, sep=''))
  counts <- read.table(tmpFile, sep=" ", as.is=T, header=FALSE)
  colnames(counts) <- c('sample', 'count')
  counts$sample <- as.integer(gsub(paste('_', this_filter, '_count.txt', sep=''), '', gsub(paste(bamDir, '/', platePrefix, '-HT', sep=''), '', counts$sample)))
  counts[order(counts$sample), ]
})
names(all_counts) <- filters

counts_merged <- MergeList(2, all_counts, mergeCol = 'sample')
names(counts_merged) <- c('sample', rev(filters))
n.barcodes <- dim(counts_merged)[1]

plot_sample <- rep(counts_merged$sample, 3)
plot_count <-c(counts_merged[, 2], counts_merged[, 3], counts_merged[, 4])
plot_filter <- factor(c(rep('clean', n.barcodes), rep('quality', n.barcodes), rep('merged', n.barcodes)), level=c('merged', 'quality', 'clean'))

counts_long <- data.frame(barcode=as.factor(plot_sample), counts=plot_count, filter=plot_filter)

## plot
plotName <- paste(platePrefix, '_QC_counts.pdf', sep='')
pdf(plotName, width=17)
my.plot <- ggplot(counts_long, aes(x=barcode, y=counts, fill=filter))
my.plot + geom_bar(position='dodge', stat='identity') + ggtitle(paste(platePrefix, ' QC counts',sep='')) + theme_bw()
dev.off()

cov.file <- paste('../../../covariates/GxE_', gsub('DP', 'P', platePrefix),
                  '_covariates.txt', sep='')
covariates <- read.table(cov.file, header=TRUE, sep='\t')
cov_short <- covariates[, c('Plate.ID', 'Barcode.ID', 'CellLine', 'Treatment.ID')]
names(counts_merged) <- c('Barcode.ID', filters[c(3, 2, 1)])
out_data <- merge(cov_short, counts_merged, by='Barcode.ID')

## For deep plates, update the Plate.ID. Additionally, write out an abbreviate
## covariate table, using the barcodes read in from the bam files. This is particularily
## useful for updating the covariate file when new treatments are added.
if ( strtrim(platePrefix, 2) == "DP" ) {
  out_data$Plate.ID <- platePrefix
  deepCovariates <- covariates[ covariates$Barcode.ID %in% out_data$Barcode.ID,
                                c("Plate.ID", "Barcode.ID", "RawReads", "QualityReads",
                                  "CleanReads", "Treatment.ID", "Treatment", "BarcodeSequence",
                                  "CellType", "CellLine", "Control.ID", "ControlCategory") ]
  deepCovariates$Plate.ID <- platePrefix
  deepCovariates$RawReads     <- out_data$merged
  deepCovariates$QualityReads <- out_data$quality
  deepCovariates$CleanReads   <- out_data$clean
  write.table(deepCovariates, file=paste0('../../../covariates/GxE_', platePrefix, '_covariates.txt'),
              quote=FALSE, row.names=FALSE, sep='\t')
}
write.table(out_data, file=paste(platePrefix, '_QC_counts.txt', sep=''), quote=FALSE, row.names=FALSE, sep='\t')

## Finally, output the aggregate sums of raw reads for each barcode/individual
## combination to check if we've hit the desired depth
tb <- aggregate(out_data$raw, by=list(out_data$CellLine, out_data$Treatment.ID), sum)
colnames(tb) <- c("Cell.Line", "Treatment.ID", "Raw.Counts")
write.table(tb, file=paste0(platePrefix, '_coverage.txt'), quote=FALSE, row.names=FALSE, sep='\t')

##
# THE END
##
