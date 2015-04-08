###########################################################################################
## System diagnostics
## CTH 4.8.2014
###########################################################################################
## all files made with commands such as
## ls *quality.bam | while read f; do echo $f `samtools view -c $f`; done > counts_quality.txt
require('ggplot2')
#library(ggplot2)

source('~/piquelab/charvey/source/aseR/aseSuite_functions_v0.0.R')

#x11(display="localhost:10.0" ,type="Xlib")

cargs <- commandArgs(trail=TRUE);
if(length(cargs)>=1)
  platePrefix <- cargs[1]
if(length(cargs)>=2)
  bamDir <- cargs[2]

## for debugging
#platePrefix <- "D1P6"
#bamDir <- "../../bams"

filters <- c('merged', 'quality', 'clean')

## collect all of the read counts from log files
all_counts <- lapply(filters, function(this_filter){
  #this_filter <- filters[1]
  tmpFile <- scan(pipe("mktemp -t"),character(0))
  system(paste('ls ../../bams/*', this_filter, '_count.txt | while read f; do cat $f | xargs echo $f; done > ', tmpFile, sep=''))
  counts <- read.table(tmpFile, sep=" ", as.is=T, header=FALSE)
  colnames(counts) <- c('sample', 'count')
  counts$sample <- as.integer(gsub(paste('_', this_filter, '_count.txt', sep=''), '', gsub(paste('../../bams/', platePrefix, '-HT', sep=''), '', counts$sample)))
  counts[order(counts$sample), ]
})
names(all_counts) <- filters

counts_merged <- MergeList(2, all_counts, mergeCol='sample')
names(counts_merged) <- c('sample', filters[c(3, 2, 1)])
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

cov.file <- paste('../../../covariates/GxE_', platePrefix, '_covariates.txt', sep='')
covariates <- read.table(cov.file, header=TRUE, sep='\t')
cov_short <- covariates[, c('Plate.ID', 'Barcode.ID', 'CellLine', 'Treatment.ID')]
names(counts_merged) <- c('Barcode.ID', filters[c(3, 2, 1)])
out_data <- merge(cov_short, counts_merged, by='Barcode.ID')

write.table(out_data, file=paste(platePrefix, '_QC_counts.txt', sep=''), quote=FALSE, row.names=FALSE, sep='\t')

##
# THE END
##
