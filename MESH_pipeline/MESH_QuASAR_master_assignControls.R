##################################################################
## Software suite for joint genotyping and ASE inference on multiple
## experiments and samples
##
## Created 07.14.2014 
## changes from *meshprep: all controls are co
##
## derived from aseSuite_v0.0_D1P6_meshprep2.R
## Author: CTH
##
## Version 0.0: Preliminary
## Arguments: plate, cell.line, covariate.file, pileup.dir, out.dir
## Return Values: Genotypes, model.convergence, inference, metaData
##################################################################
#library('ggplot2', lib.loc=myRlib)
#library('ggplot2')
##x11(display="localhost:12.0" ,type="Xlib")

## pass in the directory contianing analysis for each plate
cargs <- commandArgs(trail=TRUE);
if(length(cargs)>=1){
  top_dir <- cargs[1]
} else {
  top_dir <- '/wsu/home/groups/piquelab/charvey/GxE/new_system/jointGenotyping/'  
}

require(parallel)
cores <- as.integer(Sys.getenv("NCPUS"))
if(cores<1 || is.na(cores)){cores <- 4}
## need to get this working for lapply
ParallelSapply <- function(...,mc.cores=cores){
    simplify2array(mclapply(...,mc.cores=mc.cores))
  }
ParallelLapply <- function(...,mc.cores=cores){
    mclapply(...,mc.cores=mc.cores)
  }

##################################################################
## use the covariate table and use the barcodes to select samples
##################################################################    
cov_data <- function(plate, cell.line, treat, fields='all'){
  cov.file <- paste('~/piquelab/scratch/charvey/GxE/derived_data/covariates/GxE_', plate, '_covariates.txt', sep='')
  cv <- read.table(file=cov.file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  cv$Treatment <- gsub(' ',  '_', cv$Treatment)
  if(fields=='all'){
    cv[cv$Plate.ID==plate & cv$CellLine==cell.line & cv$Treatment.ID==treat, ]
  } else {
    cv[cv$Plate.ID==plate & cv$CellLine==cell.line & cv$Treatment.ID==treat, fields]
  }
}

#treat_dat <- read.table('./data/MASTER_treat_control.txt', stringsAsFactors=FALSE, header=FALSE)
treat_dat <- read.table('./MESH_QuASAR_master_logFC.txt', stringsAsFactors=FALSE, header=FALSE)
treat_dat <- treat_dat[complete.cases(treat_dat), ]
names(treat_dat) <- c('plate', 'cell.line', 'treatment', 'chr', 'pos0', 'pos', 'rsID', 't.id', 'g.id', 'ensg', 'beta', 'beta.se', 'pval', 'qval', 'logFC', 'padj.deg', 'pval.deg')


## output master table for MESH with both treatment and control data
########################################################
########################################################
  data <- treat_dat
  plates <- unique(data$plate)

  ## for each plate
  aux3 <- Reduce(rbind, ParallelLapply(plates, function(plate){  
    #plate <- plates[1]

    lines <- unique(data[which(data$plate==plate), 'cell.line'])

    ## for each cell line
    aux2 <- Reduce(rbind, lapply(lines, function(line){
      #line <- lines[1]
      
      treats <- unique(data[which(data$plate==plate & data$cell.line==line), 'treatment'])

      ## for each treatment
      aux1 <- Reduce(rbind, lapply(treats, function(treat){
        #treat <- treats[1]
        cat("Processing:: plate:", plate, " | cellLine:", line, " | treatment:", treat, "\n")
        
        # get the control ID
        cont <- cov_data(plate=plate, cell.line=line, treat=treat, fields='Control.ID')
        stopifnot(length(cont)==1)

        # get the treatment data
        this_treat <- data[which(data$plate==plate & data$cell.line==line & data$treatment==treat), c('rsID', 'treatment', 'beta', 'beta.se', 'logFC')]
        #str(treat_dat)
        
        # get the control data
        target <- paste0(top_dir, 'QuASAR_results_', plate, '/output/', plate, '_', line, '_cov15_', 'inference.RData')
        load(target)
        cont_dat <- inference.data[[which(names(inference.data)==cont)]][[1]][, c('annotations.rsID', 'betas', 'betas.se')]
        names(cont_dat) <- c('rsID', 'beta', 'beta.se')
        cont_dat$rsID <- as.character(cont_dat$rsID) 
        cont_dat$treatment <- cont
        cont_dat <- cont_dat[, c('rsID', 'treatment', 'beta', 'beta.se')]
        #str(cont_dat)
 
        ## merge 
        totaldata <- merge(cont_dat, this_treat, by='rsID')
        names(totaldata) <- c('rsID', 'cont', 'beta.con', 'se.con', 'treat', 'beta.treat', 'se.treat', 'logFC')
        totaldata
        
      }))

      aux1                         
                               
    }))

    aux2              
                  
  }));

  write.table(aux3, file=paste0('./MESH_QuASAR_master_logFC_controlTreat.txt'), row.names=FALSE, col.names=FALSE, quote=FALSE)

##
cat("###### THE END ######")
##
