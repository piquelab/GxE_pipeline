##################################################################
## Software suite for joint genotyping and ASE inference on multiple
## experiments and samples
##
## Combine multiple ASE analyses & differential expression data
##
## Arguments: plate, cell.line, covariate.file, pileup.dir, out.dir
## Return Values: Genotypes, model.convergence, inference, metaData
##################################################################

LPG <- '/wsu/home/groups/piquelab'

require(parallel)

cargs <- commandArgs(trail=TRUE);
if(length(cargs)>=1)
  plate <- cargs[1]
if(length(cargs)>=2)
  cell.line <- cargs[2]
## plate <- "DP1"
## cell.line <- "18507" 

cores <- as.integer(Sys.getenv("NCPUS"))
if(cores<1 || is.na(cores)){cores <- 16}

ParallelSapply <- function(...,mc.cores=cores){
    simplify2array(mclapply(...,mc.cores=mc.cores))
  }
ParallelLapply <- function(...,mc.cores=cores){
    mclapply(...,mc.cores=mc.cores)
  }


## extract covariates table
cov.name <- paste('../../derived_data/covariates/GxE_', plate, '_covariates.txt', sep='')
cv <- read.table(file=cov.name, sep="\t", header=TRUE, stringsAsFactors=FALSE)
cv <- cv[cv$Plate.ID==plate & cv$CellLine==cell.line, ]
ids <- unique(cv$Treatment.ID)[!(unique(cv$Treatment.ID) %in% c('CO1', 'CO2', 'CO3'))]
#ids <- unique(cv$Treatment.ID)

cat('Processing Cell line:: ', cell.line, '\n')

out_dat <- Reduce(rbind, ParallelLapply(ids, function(id){
   cat('Processing plate: ', plate, ' id: ', id, '\n')
  #id <- ids[1]

  ## d_dat: DGE data from DESeq2 at the level of transcripts and, same as q_dat, multiple gene_IDs
  d_dat <- read.table(paste0('../../expression/DESeq2/out_data_', plate, '/stats/',
                             plate, '_DEG_stats', '_', id, '.txt'), stringsAsFactors=FALSE, header=TRUE)
  rownames(d_dat) <- d_dat$t.id 
  #head(d_dat)
  
  ## q_dat: complete QuASAR output with trasnscript and gene IDs attached. Note that SNPs will be overlapped by
  ## multiple transcripts, hence repeated gene_IDs
  q_dat <- read.table(paste0('./output/', paste(plate, cell.line, id, 'allOutput.txt.tid', sep='_')), stringsAsFactors=FALSE)
  names(q_dat) <- c('chr', 'pos0', 'pos', 'rsID', 'beta', 'beta.se', 'pval', 'qval', 't.id', 'g.id', 'ensg')
  #head(q_dat)

  ## include data on logFC
  q_dat$padj.deg <- d_dat[q_dat$t.id, c('padj')]
  q_dat$pval.deg <- d_dat[q_dat$t.id, c('pval')]
  q_dat$logFC <- d_dat[q_dat$t.id, c('logFC')]
  
  temp <- Reduce(rbind, lapply(unique(q_dat$rsID), function(this_rsID){
    ##this_rsID <- "rs7523549"
    ##cat('\n', this_rsID, '\n')
    temp <- q_dat[which(q_dat$rsID==this_rsID), ]
    temp <- temp[complete.cases(temp), ] ## remove NAs
    minMax <- min(temp$logFC) 
    maxMax <- max(temp$logFC)
    ##cat(this_rsID, '\n')
    if(abs(minMax)>maxMax){
    	myMax <- minMax 
    } else{
    	myMax <- maxMax
    }	
    keep_ind <- which(temp$logFC==myMax)[1]
    ret_val <- temp[keep_ind,]
    row.names(ret_val) <- NULL
    ret_val
  }))

  cat('Completing plate: ', plate, ' id: ', id, '\n')
  temp$plate <- plate
  temp$cell.line <- cell.line
  temp$treatment <- id
  temp
}));

out_dat <- out_dat[, c('plate', 'cell.line', 'treatment', 'chr', 'pos0', 'pos', 'rsID', 't.id', 'g.id', 'ensg', 'beta', 'beta.se', 'pval', 'qval', 'logFC', 'padj.deg', 'pval.deg')]
fileName <- paste0('./output/', plate, '_', cell.line, '_masterTable_logFC.txt')
write.table(out_dat, fileName, col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')

##
cat("###### THE END ######", "\n")
##
