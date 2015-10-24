##################################################################
## GMB 09/12/15
##
## Combine ASE analyses & expression data across samples
##
## Arguments: 
##################################################################

library(parallel)
library(reshape)
library(data.table)
source('../../../GxE_pipeline/misc/getArgs.R')

## e.g., DP*.star.counts.fpkm.adj.meanExpr.perGene.gz
defaultList <- list(
  aligner="star",
  word="counts.fpkm.gc.meanExpr",
  type="perGene",
  cores=1,
  exprDir='../../../expression/fpkms/',
  bedTranscriptome="/wsu/home/groups/piquelab/data/RefTranscriptome/ensGene.hg19.2014.bed.gz",
  snpGeneTable="~/piquelab/gmb/SNPs/1KG/20130502/rsID_to_txID.txt"
  )
args <- getArgs(defaults=defaultList)
exprFileSuffix   <- paste(args$aligner, args$word, args$type, "gz", sep='.')
cores            <- as.numeric(args$cores)
bedTranscriptome <- args$bedTranscriptome
snpGeneTable     <- args$snpGeneTable

print(args)

ParallelSapply <- function(...,mc.cores=cores){
    simplify2array(mclapply(...,mc.cores=mc.cores))
  }
ParallelLapply <- function(...,mc.cores=cores){
    mclapply(...,mc.cores=mc.cores)
  }

## Gene annotation file
geneAnno <- read.table(bedTranscriptome, as.is=T,sep="\t")
geneAnno <- geneAnno[,-c(9:12)]
colnames(geneAnno) <-  c("chr","start","stop","t.id","score","strand","c.start",
                     "c.stop","ensg","g.id")
rownames(geneAnno) <- geneAnno$t.id

## Load in SNP-Gene annotation file, made by intersecting 1KG SNPs with bedTranscriptome
snpGT = fread(snpGeneTable, header=F)
setnames(snpGT, colnames(snpGT), c( "rsID", "t.id"))

## Iterate over plates
cov.files <- list.files('../../../derived_data/covariates/', pattern="DP", full.names=TRUE)

out_tab <- lapply(cov.files, function(cov.file) {

  ## Extract information from covariates table
  cv <- read.table(file=cov.file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  cell.type <- unique(cv$CellType)
  cell.lines <- unique(cv$CellLine)
  ids <- unique(cv$Treatment.ID)
  plate <- gsub('.*GxE_', '', gsub('_covariates.txt', '', cov.file))

  ## Look up corresponding control IDs
  c.tab <- unique(cv[, c("Treatment.ID", "Control.ID")])
  rownames(c.tab) <- c.tab$Treatment.ID

  ## Load the mean expression data (per cell type)
  expr_dat <- read.table(as.is=TRUE, sep='\t', header=TRUE,
                         paste0('../../../expression/fpkms/', plate, '/', plate, '.', exprFileSuffix))
  expr_dat$avg <- rowMeans(expr_dat[, 7:ncol(expr_dat)])
  
  ## Load the individual fpkm data
  ## Colnames are <indiv>.<treatment>
  load(paste0('../../../expression/fpkms/', plate, '/', plate, '.',
              gsub('meanExpr.*', '', exprFileSuffix), 'Rd'))
  
  ## Load the DESeq2 data, "ddsFull" and "res"
  ## List of objects, named with treatment ID
  load(paste0('../../../expression/DESeq2/out_data_', plate,
              '/data_objects/DESeq2_', plate, '.RData'))
  lfc_res <- res ## to avoid confusion/overwriting
  
  ## Get QuASAR output for each plate / cell type / treatment
  out_dat <- lapply(cell.lines, function(cell.line) {
    
    tx_dat <- Reduce(rbind, ParallelLapply(ids, function(id) {

      ## Sometimes a combo from the covariate file might not exist (e.g.,
      ## removed for QC reasons). If so, print a warning and skip.
      q_file <- paste0('../../QuASAR_results_', cell.type, "/",
                       paste(cell.type, cell.line, paste(plate, id, sep='.'),
                             'allOutput.txt', sep='_'))
      if ( !file.exists(q_file) ) {
        print(paste0("## Warning: ", paste(plate, cell.line, id, sep='-'), " not found!"))
        NULL
      } else {
        q_dat <- read.table(q_file, stringsAsFactors=FALSE, header=TRUE)      
        cat('Processing: ', plate, '/', cell.line, ' id: ', id, '\n')
        q_dat$plate <- plate
        q_dat$cell.line <- cell.line
        q_dat$treatment <- id
        q_dat$control <- c.tab[ q_dat$treatment, "Control.ID" ]
                
        ## Expand q_dat to include all transcripts overlapping a SNP
        ## May take a while, depending on number of SNPs. This will also
        ## drop SNPs that are not in annotated transcripts
        dd <- merge(q_dat, snpGT, by="rsID")
        dd$ensg <- geneAnno[ dd$t.id, "ensg"]
        
        ## Add the mean expression data
        m <- match(dd$ensg, expr_dat$ensg) # use m on expr_dat
        dd$fpkm.avg <- expr_dat$avg[m]
        dd$fpkm <- expr_dat[m, id ]

        ## Add the individual expression data
        dd$fpkm.indiv <- fpkms[ dd$t.id, paste( cell.line, id, sep='.')]
        
        ## Add the logFC data. Add 0s and NAs for controls
        if (id %in% names(lfc_res)) {
          lfc_dat <- as.data.frame(lfc_res[[id]])
          m <- match(dd$t.id, rownames(lfc_dat)) # will insert NA's where missing
          dd$lfc     <- lfc_dat[ dd$t.id, "log2FoldChange" ]
          dd$lfcSE   <- lfc_dat[ dd$t.id, "lfcSE" ]
          dd$lfcPval <- lfc_dat[ dd$t.id, "pvalue" ]
          dd$lfcPadj <- lfc_dat[ dd$t.id, "padj" ]
          dd$lfc[ is.na(dd$lfc) ] <- 0
          dd$lfcSE[ is.na(dd$lfcSE) ] <- 0
        } else {
          dd$lfc     <- NA
          dd$lfcSE   <- NA
          dd$lfcPval <- NA
          dd$lfcPadj <- NA
        }          

        ## Set a SNP id, e.g., DP1:19239:rs13469:T15C1
        dd$snp <- paste(dd$plate, dd$cell.line, dd$rsID, dd$treatment, sep=';')
        
        ## Need to choose one gene per plate-snp-treatment (snp-experiment)
        ## RPR says to keep snp-gene with max fpkm, so order by snpID then
        ## by fpkm, then remove duplicates
        ddX <- dd[order(dd$snp, dd$fpkm, decreasing=T),]
        dup <- duplicated(ddX$snp)
        dd2 <- ddX[!dup, ]
        dd2
      }
    }))
    
  })
  out_dat <- do.call('rbind', out_dat)
})
full_dat <- do.call('rbind', out_tab)
  
## Output the full master table
system('mkdir -p ../../../analysis/')
outCols <- c("snp", "ref.reads", "alt.reads", "treatment", "control", "beta",
             "beta.se", "pval", "qval", "lfc", "lfcSE", "lfcPval", "lfcPadj",
             "fpkm", "fpkm.indiv", "fpkm.avg", "cell.type", "t.id", "ensg")
write.table(full_dat[order(full_dat$plate, full_dat$chr, full_dat$pos), outCols],
            file='../../../analysis/gxeMasterTable.txt', row.names=FALSE, quote=FALSE, sep='\t')

## Prepare the data for MESH. Reshape to include the controls in the same
## record as the treatment, instead of on their own line (i.e., none of the
## controls should be listed in the "treatment" column)

## Match each SNP to its control & add control data to the table
c.tab <- unique(full_dat[, c("treatment", "control")])
rownames(c.tab) <- c.tab$treatment
ids <- sapply(full_dat$snp, function(snp) {
  tx <- unlist(strsplit(snp, ';'))[4]
  gsub(tx, c.tab[tx, "control"], snp)
})
m <- match(ids, full_dat$snp)

## Some SNPs are covered in the treatment but not the control, remove them
ind <- is.na(m)
meshIn <- full_dat[ !ind, ]
m <- na.omit(m)
cat("Warning:", sum(ind), "records thrown out due to low-coverage controls\n")

## double-check
meshIn$beta.c <- meshIn[m, "beta"]
meshIn$se.c   <- meshIn[m, "beta.se"]
meshIn$pval.c <- meshIn[m, "pval"]
meshIn$qval.c <- meshIn[m, "qval"]
meshIn$fpkm.c <- meshIn[m, "fpkm"]
meshIn$fpkm.indiv.c <- meshIn[m, "fpkm.indiv"]
  
## Update the colnames
colnames(meshIn)[grep("^beta$", colnames(meshIn))] <- "beta.t"
colnames(meshIn)[grep("^beta.se$", colnames(meshIn))] <- "se.t"
colnames(meshIn)[grep("^pval$", colnames(meshIn))] <- "pval.t"
colnames(meshIn)[grep("^qval$", colnames(meshIn))] <- "qval.t"
colnames(meshIn)[grep("^fpkm$", colnames(meshIn))] <- "fpkm.t"
colnames(meshIn)[grep("^fpkm.indiv$", colnames(meshIn))] <- "fpkm.indiv.t"

## Finally, drop the controls from the table
meshIn <- meshIn[ !(meshIn$treatment %in% meshIn$control), ]

## Order & output the completed table
outCols <- c("snp", "control", "beta.c", "se.c", "pval.c", "qval.c", "treatment",
             "beta.t", "se.t", "pval.t", "qval.t", "lfc", "lfcSE", "lfcPval",
             "lfcPadj", "fpkm.c", "fpkm.t", "fpkm.indiv.c", "fpkm.indiv.t",
             "fpkm.avg", "cell.type", "t.id", "ensg")
write.table(meshIn[ order(meshIn$chr, meshIn$pos, meshIn$plate, meshIn$treatment), outCols],
            file="MESH_input.txt", row.names=FALSE, quote=FALSE, sep=" ")

##
cat("###### THE END ######", "\n")
##
