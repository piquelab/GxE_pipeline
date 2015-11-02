##################################################################
## Software suite for joint genotyping and ASE inference on multiple
## experiments and samples
##
## Arguments: plate, cell.line, covariate.file, pileup.dir, out.dir
## Return Values: Genotypes, model.convergence, inference, metaData
##################################################################

library('QuASAR')
library('ggplot2')
library('parallel')
library('qqman')
require('qvalue')

cores = 1 # default
cargs <- commandArgs(trail=TRUE);
if(length(cargs)>=1)
    cell.type <- cargs[1]
if(length(cargs)>=2)
    cores <- cargs[2]

## Helper functions
ParallelSapply <- function(...,mc.cores=cores){
    simplify2array(mclapply(...,mc.cores=mc.cores))
  }

output.folder <- paste0('QuASAR_results_', cell.type, '/')
#system(paste0('mkdir -p ', output.folder)) # already made with prep script
system(paste0('mkdir -p ', output.folder, '/plots/QQ'))
system(paste0('mkdir -p ', output.folder, '/plots/QC/oraf'))
system(paste0('mkdir -p ', output.folder, '/plots/QC/coverage'))

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

ParallelSapply(unique(cv$CellLine), function(cell.line) {

  ## Load in the input data
  ## Loads 'ase.dat' and 'cov.file'
  load(paste0(output.folder, '/data/', cell.type, '_', cell.line, '_quasarIn.Rd'))
  
  barcodes <- paste0(cov.file$Plate.ID, "-HT", cov.file$Barcode.ID) #cov.file$Barcode.ID
  treatments <- cov.file$Treatment
  n.treatments <- length(treatments)
  treatment.IDs <- cov.file$Treatment.ID
  controls <- unique(cov.file$Control.ID)
  
  min.cov <- 15
  ase.dat.gt <- PrepForGenotyping(ase.dat, min.coverage=min.cov)
  str(ase.dat.gt)
  
  ## prep unamed SNPs for MESH
  tempchr <- ase.dat.gt$annotations[which(ase.dat.gt$annotations$rsID=='.'), 'chr']
  temppos <- ase.dat.gt$annotations[which(ase.dat.gt$annotations$rsID=='.'), 'pos']
  ase.dat.gt$annotations[which(ase.dat.gt$annotations$rsID=='.'), 'rsID'] <- paste(tempchr, temppos, sep='-')
  
  ##################################################################
  ## QC
  ##################################################################
  if(TRUE){
    for(ii in (1:n.treatments)){
      ##ii <- 1
      cat(ii)
      
      oraf <- ase.dat$ref[, ii]/(ase.dat$ref[, ii]+ase.dat$alt[, ii])
      ind <- ((ase.dat$ref[, ii] + ase.dat$alt[, ii]) > min.cov) & (oraf>0) & (oraf<1)
      oraf <- oraf[ind]
      aux <- data.frame(oraf=oraf)
      
      pdf.file <- paste0(output.folder, '/plots/QC/oraf/', cell.type, '_', cell.line, '_', treatment.IDs[ii], '_',  ii, '_QC_oraf.pdf', sep='')
      pdf(file=pdf.file)
      ##hist(oraf, breaks=100)
      qc_plot <- ggplot(aux, aes(x=oraf))
      print(qc_plot +
            geom_histogram(binwidth=0.01) +
            ggtitle(paste(cell.type, '_', cell.line, '_', treatments[ii], sep=''))) +
              theme_bw()
      dev.off()
    }}
  
  ########################################################################
  ##chrCol <- rainbow(length(chrList))
  chrList <- paste("chr",sort(as.numeric(gsub("chr","",unique(ase.dat$anno$chr)))),sep="")
  
  for(ii in (1:n.treatments)){
    ##ii <- 3
    depth <- (ase.dat$ref[, ii]+ase.dat$alt[, ii])
    ind <- ((depth>min.cov) & (depth<1000))
    af <- ase.dat$anno$af
    oraf <- ase.dat$ref[, ii]/(ase.dat$ref[, ii]+ase.dat$alt[, ii])
    chr <- ase.dat$anno$chr
    sName <- treatment.IDs[ii]
    chrCol <- rep(c("orange","darkblue"),length(chrList))[1:length(chrList)]
    names(chrCol) <- chrList
    pdf.file <- paste0(output.folder, '/plots/QC/coverage/', cell.type, '_', cell.line, '_', sName, '_',  ii, '_QC_manhattan.pdf')
    pdf(file=pdf.file,width=16,height=8)
    layout(t(c(1,1,1,2)))
    par(cex=1.0)
    oldmar <- par("mar")
    mar <- oldmar
    mar[4] <- 0.0
    par(mar=mar)
    ind2 <- ind & (abs(af-0.5)<0.4)
    ##x <- 1:sum(ind2)
    plot(oraf[ind2],xlab="",ylab="Obs. reference allele Freq.",pch='.',cex=3,col=chrCol[chr[ind2]],axes=F)
    axis(2)
    x.at <- c(which(!duplicated(chr[ind2])),sum(ind2))
    axis(1,at=x.at,labels=FALSE,cex=0.3)
    abline(v=x.at,lty=3)
    text((x.at[-1]+x.at[-length(x.at)])*0.5, -0.15, labels = chrList, srt = 45, pos =3, xpd = TRUE,cex=0.7)
    title(treatments[ii])
    abline(h=0.5,lty=3)
    mar <- oldmar
    mar[2] <- 0.0
    par(mar=mar)
    aux <- hist(jitter(oraf[ind2]),breaks=101,plot=F)
    barplot(aux$counts+1,horiz=TRUE,log="x")
    title(xlab="Freq.")
    par(mar=oldmar)
    dev.off()
  }
  
  ##################################################################
  ## Collapse the controls
  ##################################################################
  plates <- unique(gsub('-HT.*', '', barcodes))
  
  controlref <- do.call('cbind', lapply(plates, function(p) {
    do.call('cbind', lapply(controls, function(this){rowSums(ase.dat.gt$ref[, barcodes[ intersect(grep(this, treatment.IDs), grep(p, barcodes))  ] ])}))
  }))
  colnames(controlref) <- sapply(plates, function(p) { sapply(controls, function(c) { paste0(p, '-', c) }) })
  
  controlalt <- do.call('cbind', lapply(plates, function(p) {
    do.call('cbind', lapply(controls, function(this){rowSums(ase.dat.gt$alt[, barcodes[ intersect(grep(this, treatment.IDs), grep(p, barcodes))  ]  ])}))
  }))
  colnames(controlalt) <- sapply(plates, function(p) { sapply(controls, function(c) { paste0(p, '-', c) }) })

  ## If any are zero, i.e., in one plate & not the other, remove
  torm = which(apply(controlref, 2, sum)==0)
  if ( length(torm) > 0 ) {
    controlref <- controlref[, -which(apply(controlref, 2, sum) == 0) ]
  }
  torm = which(apply(controlalt, 2, sum)==0)
  if ( length(torm) > 0 ) {
    controlalt <- controlalt[, -which(apply(controlalt, 2, sum) == 0) ]
  }
  
  finalref <- cbind(controlref, ase.dat.gt$ref[, barcodes[ -grep("CO", treatment.IDs) ]])
  finalalt <- cbind(controlalt, ase.dat.gt$alt[, barcodes[ -grep("CO", treatment.IDs) ]])

  txTab <- unique(cv[, c("Treatment.ID", "Treatment")])
  rownames(txTab) <- txTab$Treatment.ID

  ## Create a lookup table to get the final treatment names
  ## rbind is fickle: needs same types (no coersion), same names (or none)
  a <- unname(unique(cv[-grep('CTRL', cv$Treatment), c("Plate.ID", "Barcode.ID", "Treatment.ID", "Treatment")]))
  a[, 2] <- as.character(a[, 2])
  b <- unname(unique(cv[grep('CTRL', cv$Treatment), c("Plate.ID", "Treatment.ID", "Treatment.ID", "Treatment")]))
  names(a) <- names(b) <- c("Plate.ID", "ID", "Treatment.ID", "Treatment")
  txTab <- rbind(a, b)
  txTab$collapsed <- paste(txTab$Plate.ID, txTab$Treatment, sep='_')
  txTab$final <- paste(txTab$Plate.ID, txTab$Treatment.ID, sep='.')
  rownames(txTab) <- paste(txTab$Plate.ID, txTab$ID, sep='-')

  treatments_collapsed <- txTab[gsub('HT', '', colnames(finalref)), "collapsed"]
  treatmentIDs_final <- txTab[gsub('HT', '', colnames(finalref)), "final"]
  ##treatments_collapsed <- c(unique(unlist(subset(cov.file, Treatment.ID %in% controls, Treatment))), treatments[ -grep("CO", treatment.IDs) ])
  ##treatmentIDs_final <- c(controls, treatment.IDs[ -grep("CO", treatment.IDs) ])
  
  colnames(finalref) <- treatments_collapsed
  colnames(finalalt) <- treatments_collapsed
  
  ##################################################################
  ## QuASAR Model Fitting
  ## ase.joint ~ object for joint genotyoping across all samples
  ##################################################################
  ## Ensure the input to fitAseNullMulti is identical to ase.dat.final below
  ase.joint <- fitAseNullMulti(finalref, finalalt, log.gmat=log(ase.dat.gt$gmat))

  ## Mandatory objects for inference
  ase.joint.eps <- ase.joint$eps
  n.eps <- length(ase.joint.eps)
  sample.names <- treatments_collapsed
  treatments.final <- treatments_collapsed
  ase.dat.final <- list(ref=finalref, alt=finalalt, gmat=ase.dat.gt$gmat, annotations=ase.dat.gt$annotations)
  
  ##################################################################
  ## Output model data; genotypes, etc.
  ## ase.joint ~ object for joint genotyoping across all samples
  ##################################################################
  out.gts <- data.frame(rsID=ase.dat.final$annotations$rsID, g0=ase.joint$gt[, 'g0'], g1=ase.joint$gt[, 'g1'], g2=ase.joint$gt[, 'g2'])
  dat.name <- paste(output.folder, "/", cell.type, "_",cell.line,'_genotypes.txt', sep='')
  write.table(out.gts, file=dat.name, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  
  ##################################################################
  ## inference
  ##################################################################
  inference.data <- lapply(seq_along(1:n.eps), function(ii){
    ##################################################################
    ## sample ~ current sample to assess ASE
    ## this.sample ~ sample name
    ## coverage ~ coverage for this sample
    ## coverage.floor ~ minimum sample wide coverage
    ## coverage.ind ~ indicator for sufficient coverage of this sample
    ## ref ~ reference count for this sample with sufficient covergae
    ## alt ~ alternate count for this sample with sufficient covergae
    ## phi ~ genotype priors for this sample
    ## eps ~ jointly inferred error rate for this sample
    ## het ~ jointly inferred heterozygote probabilities for sites with
    ##       sufficient coverage
    ## het.ind ~ indicator for heterozygotes with p > .99
    ## annotations ~ annotations filtered by coverage and het probability
    ## q.thresh ~ q-value threshold for declaring signifigance
    ## DEBUGGING
    ##ii <- 1
    sample <- ii
    this.sample <- sample.names[sample]
    coverage <- (ase.dat.final$ref[, sample] + ase.dat.final$alt[, sample])
    coverage.floor <- 5
    coverage.ind <- (coverage>coverage.floor)
    ref <- ase.dat.final$ref[coverage.ind, sample]
    alt <- ase.dat.final$alt[coverage.ind, sample]
    tot <- ref + alt
    phi <- ase.dat.final$gmat[coverage.ind]
    eps <- ase.joint.eps[sample]
    het <- ase.joint$gt[coverage.ind, 2]
    het.ind <- (het > 0.99)
    numb.hets <- sum(het.ind)
    annotations <- ase.dat.final$annotations[coverage.ind, ][het.ind, ]
    q.thresh <- 0.2
    this.treatment <- treatmentIDs_final[sample]
    
    cat("#Hets:", this.sample, this.treatment, numb.hets, "\n")
    
    ##################################################################
    ## M ~ a grid of possible dispersion values for the Beta-binomial model
    ## aux ~ loglikelihood of the beta bionmial model across D values
    ##               with null \rho value
    ## Mmax ~ disperison which maximizes the llk
    cov_breaks <- c((0:6)*10,80,100,250,500,1000,100000)
    bin <- cut(tot,cov_breaks)
    tapply(tot,bin,mean)
    M <- exp((0:500)/50)
    Mvec <- sapply(levels(bin),function(mybin){
      aux <- sapply(M,function(M){
        sum(logLikBetaBinomialRhoEps(0.5,eps,M,ref[het.ind & bin==mybin ],alt[het.ind & bin==mybin]))
      })
      Mmax <- M[which.max(aux)]
      Mmax
    })
    
    cat("#Disp: ", round(Mvec, 2), "\n")
    
    ##################################################################
    ## Find MLE for \rho using the the Beta-Binomial model
    ## Dmax2 ~ the dispersian parameter estimed from the llk in the
    ##         previous step
    ## auxLogis ~ optimization of the Beeta-biomial model in terms of
    ##            logit(\rho)
    ## rho3 ~ vector of rho estimates from expit(logit(\rho))
    ## lrt3 ~ Recaluclate Het LRT using the beta-bionomial
    ## pval3 ~ pval of thr llk ratio test using the chi-square approximation
    ##         (does not include uncertainty in genotyping)
    ## qv3 ~ qvalue object from llkRT p-values
    ## qvals.qv3 ~ qvalues from the above p-values
    ## betas.beta.binom ~ logit(\rho) or beta value for heterozygotes
    ## betas.se ~ standard error of the beta value
    ## betas.z ~ Z scores of the heterozygote beta values
    ## betas.pval ~ pvalues for the above z-scores
    ##auxLogis <- optim(rep(0,sum(het.ind)),fn=logLikBetaBinomial2,
    ##gr=gLogLikBetaBinomial,D=Dmax2,R=ref[het.ind],A=alt[het.ind],
    ##method="L-BFGS-B", hessian=TRUE)

    ##
    ## unbounded optimization
    ##
    aux2 <- t(sapply(1:sum(het.ind),function(ii){
      auxLogis <- optim(0,
                        fn=logLikBetaBinomial2,
                        gr=gLogLikBetaBinomial,
                        D=Mvec[bin[het.ind][ii]],
                        R=ref[het.ind][ii],
                        A=alt[het.ind][ii],
                        method="L-BFGS-B",
                        hessian=TRUE,
                        lower=-5,
                        upper=5)
      c(auxLogis$par,1/(auxLogis$hessian)^.5)
    }))
    
    rho3 <- plogis(aux2[,1])
    betas.beta.binom <- aux2[,1]
    lrt3 <-
      logLikBetaBinomialRhoEps(rho3,eps,Mvec[bin[het.ind]],ref[het.ind],alt[het.ind]) -
        logLikBetaBinomialRhoEps(0.5,eps,Mvec[bin[het.ind]],ref[het.ind],alt[het.ind])
    pval3 <- (1-pchisq(2*lrt3,df=1))
    betas.se <- abs(betas.beta.binom/qnorm(pval3/2))
    betas.se[which(betas.se=='NaN')] <- aux2[, 2][which(betas.se=='NaN')]
    betas.z <- betas.beta.binom/betas.se
    betas.pval <- 2*pnorm(-abs(betas.z))
    
    ## Stratified qvalue. Save for further testing.
    ## qvals.qv3bins <- rep(1,length(pval3));
    ## pi0.vec <- sapply(levels(bin),function(mybin){
    ##   ind<-bin[het.ind]==mybin
    ##   myqv <- qvalue(pval3[ind],pi0.method='bootstrap')
    ##   qvals.qv3bins[ind] <<- myqv$qv
    ##   myqv$pi0
    ## })
    qv3 <- qvalue(pval3, pi0.method='bootstrap')
    qvals.qv3 <- qv3$qv
    this.pi0 <- qv3$pi0
    qv.05 <- sum(qvals.qv3<0.05)
    qv.1 <- sum(qvals.qv3<0.1)
    qv.2 <- sum(qvals.qv3<0.2)
    qv.5 <- sum(qvals.qv3<0.5)
    
    ##cat("#Loci ", annotations$rsID[which(qvals.qv3<q.thresh)], "\n", sep=" ")
    cat("#ASE:", this.treatment, qv.2, paste0("(", q.thresh*100, "%FDR)"), "Pi0:", round(this.pi0, 3), "\n")

    ## output complete information for every heterozygote
    complete.dat <- annotations
    complete.dat$cell.type <- cell.type
    complete.dat$cell.line <- cell.line
    complete.dat$treatment <- this.treatment
    complete.dat$ref.reads <- ref[het.ind]
    complete.dat$alt.reads <- alt[het.ind]
    complete.dat$beta <- betas.beta.binom
    complete.dat$beta.se <- betas.se
    complete.dat$pval <- pval3
    complete.dat$qval <- qvals.qv3
    filename.all <- paste0(output.folder, '/', cell.type, '_', cell.line, '_', this.treatment, '_allOutput.txt')
    write.table(complete.dat, file=filename.all, row.names=FALSE, col.name=TRUE, quote=FALSE)
    
    ## return data frame
    rsID <- annotations
    betas <- betas.beta.binom
    dat <- data.frame(annotations$rsID, annotations$chr, annotations$pos0, rho=rho3, betas, betas.se, qv=qvals.qv3, pval=pval3, refCount=ref[het.ind], altCount=alt[het.ind])
    meta.dat <- data.frame(hets=numb.hets, pi0=this.pi0, qv.05=qv.05, qv.1=qv.1, qv.2=qv.2, qv.5=qv.5, mean.rho=mean(rho3), median.rho=median(rho3)) # ,dispersion=Dmax2)
    temp <- list(dat=dat, meta.dat=meta.dat)
    
  }) ## Returns a list of data & metaData

  names(inference.data) <- treatmentIDs_final
  dat.name <- paste(output.folder, "/", cell.type, "_",cell.line, "_cov", min.cov, '_inference.RData', sep='')
  save(inference.data, file=dat.name)
  str(inference.data)
  ##load(dat.name)
  
  ##########################################
  ## 0.) plots for average rho aross the individual
  ##########################################
  all_rho_hat <- Reduce(c, sapply(seq_along(1:length(inference.data)), FUN=function(ii){inference.data[[ii]]$dat$rho}))
  mean_rho_hat <- round(mean(all_rho_hat), 4)
  median_rho_hat <- round(median(all_rho_hat), 4)
  se_rho_hat <- sd(all_rho_hat)
  
  all_ref <- Reduce(c, sapply(seq_along(1:length(inference.data)), FUN=function(ii){inference.data[[ii]]$dat$refCount}))
  all_alt <- Reduce(c, sapply(seq_along(1:length(inference.data)), FUN=function(ii){inference.data[[ii]]$dat$altCount}))
  all_coverage <- '+'(all_ref, all_alt)
  emp_rho <- '*'(all_ref, all_coverage^(-1))

  rho_title <- paste0(cell.type, '-', cell.line, ' : Rho_hat across all treatements', '\n mean.Rho=', mean_rho_hat, ' | median.Rho=', median_rho_hat)
  pdf.file <- paste0(output.folder, '/plots/QC/', cell.type, '_', cell.line, '_averageRho', '.pdf', sep='')
  pdf(file=pdf.file)
  hist(all_rho_hat[all_coverage>100], main=rho_title, xlim=c(0,1), breaks=seq(0,1,0.01), col='darkgrey', axes=FALSE)
  abline(v=mean_rho_hat, lty=2, col='red')
  axis(1, at=seq(0, 1, .1)); axis(2)
  dev.off()
  
  ##########################################
  ## 1.) QQ-plots for all treatments
  ##########################################
  for(ii in seq_along(1:length(inference.data))){
    ##ii <- 1
    treatment <- names(inference.data)[ii]
    pvals <- inference.data[[ii]]$dat$pval
    pi0 <- round(inference.data[[ii]]$meta.dat$pi0, 2)
    hets <- inference.data[[ii]]$meta.dat$hets
    qv.2 <- inference.data[[ii]]$meta.dat$qv.2
    coverage <- inference.data[[ii]]$dat$refCount + inference.data[[ii]]$dat$altCount
    avg.depth <- floor(mean(coverage))
    ##disp <- round(mean(inference.data[[ii]]$meta.dat$dispersion), 2)
    
    ## extract pvalues from only the high coverage loci
    pval_high <- inference.data[[1]]$dat$pval[which(coverage>100)]
    qqp <- qqplot(-log10(ppoints(length(pval_high))),-log10(pval_high), plot.it=F)
    
    pdf.file <- paste(output.folder, '/plots/QQ/', cell.type, '_', cell.line, '_', treatment, "_cov", min.cov, '_', ii, '_QQ', '.pdf', sep='')
    title <- paste0(cell.line, " | ", treatments_collapsed[ii], " | Pi0=", pi0, " | #hets=", hets, "\n #qv.2=", qv.2, " | avg.depth=", avg.depth) #, " | disp=", disp)
      pdf(file=pdf.file)
    qq(pvals)
    points(qqp,pch='.',cex=5,col='blue')
    title(main=title)
    legend("topleft", c(paste("all hets"), paste("hets with minCov=100")),
           fill=c("black","blue"))
    dev.off()
  }

  ##########################################
  ## 2.) Expression table across all treatments
  ##########################################
  asetable <- t(sapply(seq_along(1:length(inference.data)), FUN=function(ii){
    ##ii <- 1
    sapply(c(.01, .05, .1, .2), FUN=function(jj){sum(inference.data[[ii]]$dat$qv < jj)})
  }))
  
  rownames(asetable) <- names(inference.data)
  colnames(asetable) <- c('Q<.01', 'Q<.05', 'Q<.1', 'Q<.2')
  
  outfile <- paste(output.folder, '/', cell.type, "_",cell.line, '_Qhits.txt', sep='')
  write.table(asetable, file=outfile, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
})

##                          ##
cat("###### THE END ######")
##                          ##
