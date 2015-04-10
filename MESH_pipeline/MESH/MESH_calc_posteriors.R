#####################################################
##
## CTH 012215
## QUASAR -> MESH -> multiclassLR pipeline
## From MESH output, calculate bayes facotrs and posteriors
##
#####################################################
## check default ncpus for myParallel.R
#source('~/piquelab/charvey/source/myRScripts/myParallel.R')
library(reshape2)

parmfiles <- scan(pipe('ls *.hm_config.est'), character(0))
bffiles <- scan(pipe('ls *.hm_config.in'), character(0))
n.samples <- length(parmfiles)

for(ii in 1:n.samples){
  # ii <- 1  
  #parmf <- "P1_18507_T15C1.hm_config.est"
  #bffile <- "P1_18507_T15C1.hm_config.in"

  parmf <- parmfiles[ii]
  bffile <- bffiles[ii]   

  pi0 <- as.numeric(scan(pipe(paste("less ", parmf, " | awk '$1==\"pi0:\" {print $2}'", sep='')),character(0)))
  pi00 <- 1-pi0

  config1 <- as.numeric(scan(pipe(paste("less ", parmf, " | awk '$1==\"config:\" {print $2}'", sep='')),character(0)))
  config2 <- as.numeric(scan(pipe(paste("less ", parmf, " | awk '$1==\"config:\" {print $4}'", sep='')),character(0)))
  config3 <- as.numeric(scan(pipe(paste("less ", parmf, " | awk '$1==\"config:\" {print $6}'", sep='')),character(0)))
  config <- c(config1, config2, config3)

  config <- c(0.1,0.1,0.8)

  priors <- c(pi00 * config, pi0)
  sum(priors)

  grid1 <- as.numeric(scan(pipe(paste("less ", parmf, " | awk '$1==\"grid:\" {print $2}'", sep='')),character(0)))
  grid2 <- as.numeric(scan(pipe(paste("less ", parmf, " | awk '$1==\"grid:\" {print $4}'", sep='')),character(0)))
  grid3 <- as.numeric(scan(pipe(paste("less ", parmf, " | awk '$1==\"grid:\" {print $6}'", sep='')),character(0)))
  grid <- c(grid1, grid2, grid3)

  data <- read.table(bffile, as.is=TRUE)
  names(data) <- c('rsID', 'config', 'bf1', 'bf2', 'bf3')

  data$abfs <- 10^as.matrix(data[, c('bf1', 'bf2', 'bf3')]) %*% grid
  bfs <- cbind(acast(data, rsID ~ config ,value.var="abfs"), 1)
  const <- 1/ (bfs %*% priors)

  post <- sapply(1:4,function(col){
    bfs[,col] * priors[col] * const
  })

  MyPaste <- function(...){paste(..., sep='_')}
  split_name <- function(name)
  {
    names_snps <- strsplit(name, '_')[[1]][1:5]
    name_out <- Reduce(MyPaste, names_snps)
    name_out
  }
  
  names_snps <- sapply(rownames(bfs), split_name, USE.NAMES=FALSE)

  rownames(post) <- names_snps
  rownames(bfs) <- names_snps
  colnames(bfs) <- c('bf_c1', 'bf_c2', 'bf_c3', 'bf_c4')
  colnames(post) <- c('post_c1', 'post_c2', 'post_c3', 'post_c4')

  bfs_out <- data.frame(bfs, stringsAsFactors=FALSE)
  post_out <- data.frame(post, stringsAsFactors=FALSE)

  myDir <- './posteriors_bfs'
  system(paste0('mkdir -p ', myDir))
  root <- gsub('.hm_config.est', '', parmf)

  write.table(post_out, file=paste(myDir, '/', root, '_posteriors.txt', sep=''), quote=FALSE, row.names=TRUE, sep='\t')
  write.table(bfs_out, file=paste(myDir, '/', root, '_bayesFactors.txt', sep=''), quote=FALSE, row.names=TRUE, sep='\t')

}

############################################################################################
## create a master table with MESH input (betas & SEs) & MESH output (bayes-factors & posteriors)
############################################################################################
#install.packages('data.table')
library(data.table)

## capture MESH input
meshinput <- scan(pipe('ls allPlates_*.txt'), character(0))
if(length(meshinput) != 1){stop("Pipeline is set to process a single partition of data.")}
dfbetas <- read.table(meshinput, stringsAsFactors=FALSE, header=FALSE)  
dfbetas_wide <- data.frame(dfbetas[seq(1, dim(dfbetas)[1], 2), c(1, 3, 4)], dfbetas[seq(2, dim(dfbetas)[1], 2), c(3, 4)])
dtbetas <- data.table(dfbetas_wide)
setnames(dtbetas, c('snp', 'beta.c', 'se.c', 'beta.t', 'se.t'))
setkey(dtbetas, snp)

## capture Bayes Factors
dfbfs <- read.table("posteriors_bfs/allPlates_allSNPs_bayesFactors.txt", stringsAsFactors=FALSE, header=TRUE)
dtbfs <- data.table(dfbfs)
setkey(dtbfs, snp)

## capture posteriors
dfpost <- read.table("posteriors_bfs/allPlates_allSNPs_posteriors.txt", stringsAsFactors=FALSE, header=TRUE)
dtpost <- data.table(dfpost)
setkey(dtpost, snp)

## check the size of tables that have been loaded
tables()

dtfull <- dtbetas[dtbfs[dtpost]]
write.table(dtfull, file='./posteriors_bfs/Master_table_betas_bfs.txt', quote=FALSE, col.names=TRUE, row.names=FALSE)

##         ##
## THE END ##
##         ##
