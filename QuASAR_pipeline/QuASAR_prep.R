##################################################################
## Prepare data for input to QuASAR for joint genotyping and ASE
##
## Arguments: cell.type
## Return Values: 
##################################################################

library('QuASAR')
library('ggplot2')

cargs <- commandArgs(trail=TRUE);
if(length(cargs)>=1)
  cell.type <- cargs[1] 

system(paste0('mkdir -p ', cell.type))
system(paste0('mkdir -p ', cell.type, '/output'))

##################################################################
## specify covariate file and pileup directories
##################################################################    
covDir <- '../derived_data/covariates/'
cmd <- paste0('cat ', covDir, 'GxE_DP*_covariates.txt | grep -w ', cell.type)
cv <- read.table(pipe(cmd), as.is=TRUE, sep='\t')
names(cv) <- c("Plate.ID", "Barcode.ID", "Raw.Reads", "Quality.Reads", "Clean.Reads",
               "Treatment.ID", "Treatment", "BarcodeSequence", "CellType", "CellLine",
               "Control.ID", "Control")
cv$Treatment <- gsub(' ',  '_', cv$Treatment)

## Loop over individuals
sapply(unique(cv$CellLine), function(cell.line) {

  cov.file <- cv[ cv$CellLine == cell.line, ]
  
  files <- apply(cov.file, 1, function(x) {
    f <- paste0('../derived_data/', x["Plate.ID"], '/pileups/', x["Plate.ID"], '-HT',
                x["Barcode.ID"], '.pileup.clean.bed.gz')
    f <- gsub(' ', '', f) # not sure why it inserts space into single-digit BCs
  })

  ## TODO: place in helper script, source in
  UnionExtractFields <- function (fileList, combine = FALSE) 
    {
      tmpFile <- scan(pipe("mktemp -t"), character(0))
      system(
        paste("zcat ",
              paste(fileList, collapse = " "),
              " | grep -v -w '^chr\\|^NA' | cut -f 1-4,6-7 | sortBed -i stdin | uniq | gzip > ",
              tmpFile))
      anno <- read.table(gzfile(tmpFile), sep = "\t", as.is = T)
      aux <- sapply(fileList, function(fn) {
        cat("Processing:", fn, "\n")
        command = paste("intersectBed -a ", tmpFile, " -b ", 
          fn, " -wao | cut -f 1-3,14-16 ", sep = "")
        aa <- read.table(pipe(command), sep = "\t", as.is = T, 
                         na.strings = ".")
        aa[is.na(aa)] <- 0
        stopifnot(identical(aa[, 1:3], anno[, 1:3]))
        aa[, -(1:3)]
      })
      colnames(anno) = c("chr", "pos0", "pos", "ref", "rsID", "af")
      Ref <- as.matrix(do.call(cbind, aux[1, ]))
      Alt <- as.matrix(do.call(cbind, aux[2, ]))
      Err <- as.matrix(do.call(cbind, aux[3, ]))
      return.list <- list(ref = Ref, alt = Alt, err = Err, anno = anno)
      if (combine == TRUE) {
        allRef <- apply(Ref, MARGIN = 1, sum)
        allAlt <- apply(Alt, MARGIN = 1, sum)
        allErr <- apply(Err, MARGIN = 1, sum)
        return.list$all <- as.matrix(cbind(allRef, allAlt, allErr))
      }
      return(return.list)
    }

  ##################################################################    
  ## prepare the relevant samples
  ##################################################################    
  ase.dat <- UnionExtractFields(files, combine=TRUE)
  for( ii in seq_along(1:3) ) {
    colnames(ase.dat[[ii]]) <- gsub('.*pileups/', '', gsub('\\.pileup.*', '', files))
  }
  save(list=c('ase.dat', 'cov.file'),
       file=paste0(cell.type, '/output/', cell.type, '_', cell.line, '_quasarIn.Rd'))
})  

## the end ##
