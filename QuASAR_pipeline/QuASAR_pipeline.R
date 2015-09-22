##################################################################
## Software suite for joint genotyping and ASE inference on multiple
## experiments and samples
##
## Arguments: plate, cell.line, covariate.file, pileup.dir, out.dir
## Return Values: Genotypes, model.convergence, inference, metaData
##################################################################

library('QuASAR')
library('ggplot2')
library(qqman)
require(qvalue)


##################################################################
## use the covariate table and use the barcodes to select samples
##################################################################    
cargs <- commandArgs(trail=TRUE);
if(length(cargs)>=1)
  platePrefix <- cargs[1]
if(length(cargs)>=2)
  cell.line <- cargs[2]

system('mkdir -p output')
system('mkdir -p plots/QQ')
system('mkdir -p plots/QC/oraf')
system('mkdir -p plots/QC/manhattan')

##################################################################
## specify covariate file and pileup directories
##################################################################    
output.folder <- paste('./output',sep='')
## extract covariates table
cov.file <- paste('../../derived_data/covariates/GxE_', platePrefix, '_covariates.txt', sep='')
cv <- read.table(file=cov.file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
cv$Treatment <- gsub(' ',  '_', cv$Treatment)
cov.file <- cv[cv$Plate.ID==platePrefix & cv$CellLine==cell.line, ]

barcodes <- cov.file$Barcode.ID
treatments <- cov.file$Treatment
n.treatments <- length(treatments)
treatment.IDs <- cov.file$Treatment.ID
controls <- unique(cov.file$Control.ID)


UnionExtractFields <- function (fileList, combine = FALSE) 
{
    tmpFile <- scan(pipe("mktemp -t"), character(0))
    system(
      paste("zcat ",
            paste(fileList, collapse = " "),
#            " | grep -v -w '^chr\\|^NA' | cut -f 1-3,6-7 | sortBed -i stdin | uniq | gzip > ",
            " | grep -v -w '^chr\\|^NA' | cut -f 1-4,6-7 | sortBed -i stdin | uniq | gzip > ",
            tmpFile))
    anno <- read.table(gzfile(tmpFile), sep = "\t", as.is = T)
    aux <- sapply(fileList, function(fn) {
        cat("Processing:", fn, "\n")
        command = paste("intersectBed -a ", tmpFile, " -b ", 
#            fn, " -wao | cut -f 1-3,13-15 ", sep = "")
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
files <- paste0('../../derived_data/', platePrefix, '/pileups/', platePrefix,
                "-HT", barcodes, ".pileup.clean.bed.gz")
ase.dat <- UnionExtractFields(files, combine=TRUE)
for(ii in seq_along(1:3)){colnames(ase.dat[[ii]])<-paste("HT", barcodes,sep='')}
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
		#ii <- 1
		cat(ii)

		oraf <- ase.dat$ref[, ii]/(ase.dat$ref[, ii]+ase.dat$alt[, ii])
		ind <- ((ase.dat$ref[, ii] + ase.dat$alt[, ii]) > min.cov) & (oraf>0) & (oraf<1)
		oraf <- oraf[ind]
		aux <- data.frame(oraf=oraf)

		pdf.file <- paste('./plots/QC/oraf/', platePrefix, '_', cell.line, '_', treatment.IDs[ii], '_',  ii, '_QC_oraf.pdf', sep='')	
		pdf(file=pdf.file)
		#hist(oraf, breaks=100)
		  qc_plot <- ggplot(aux, aes(x=oraf))
		  print(qc_plot +
                        geom_histogram(binwidth=0.01) +
                        ggtitle(paste(platePrefix, '_', cell.line, '_', treatments[ii], sep=''))) +
                        theme_bw()
		dev.off()
	}}

########################################################################
##chrCol <- rainbow(length(chrList))
chrList <- paste("chr",sort(as.numeric(gsub("chr","",unique(ase.dat$anno$chr)))),sep="")

for(ii in (1:n.treatments)){
	#ii <- 3
	depth <- (ase.dat$ref[, ii]+ase.dat$alt[, ii])
	ind <- ((depth>min.cov) & (depth<1000))
	af <- ase.dat$anno$af
	oraf <- ase.dat$ref[, ii]/(ase.dat$ref[, ii]+ase.dat$alt[, ii])
	chr <- ase.dat$anno$chr
	sName <- treatment.IDs[ii]
	chrCol <- rep(c("orange","darkblue"),length(chrList))[1:length(chrList)]
	names(chrCol) <- chrList
	pdf.file <- paste('./plots/QC/manhattan/', platePrefix, '_', cell.line, '_', sName, '_',  ii, '_QC_manhattan.pdf', sep='')	
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
controlref <- sapply(controls, function(this){rowSums(ase.dat.gt$ref[, paste0("HT", unlist(subset(cov.file, Treatment.ID==this, Barcode.ID)))])}) 
controlalt <- sapply(controls, function(this){rowSums(ase.dat.gt$alt[, paste0("HT", unlist(subset(cov.file, Treatment.ID==this, Barcode.ID)))])}) 

treat.ind <- !(colnames(ase.dat.gt$ref) %in%  paste0("HT", unlist(subset(cov.file, Treatment.ID %in% controls, Barcode.ID))))
finalref <- cbind(controlref, ase.dat.gt$ref[, treat.ind])
finalalt <- cbind(controlalt, ase.dat.gt$alt[, treat.ind])
treatments_collapsed <- c(unique(unlist(subset(cov.file, Treatment.ID %in% controls, Treatment))), treatments[treat.ind])
treatmentIDs_final <- c(controls, treatment.IDs[treat.ind])
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
dat.name <- paste(output.folder, "/",platePrefix, "_",cell.line,'_genotypes.txt', sep='')
write.table(out.gts, file=dat.name, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

for(ii in 1:length(treatments)){
	this.treat <- treatment.IDs[ii]	
	out.counts <- data.frame(rsID=ase.dat.gt$annotations$rsID, ref=ase.dat.gt$ref[, ii], alt=ase.dat.gt$alt[, ii])
	dat.name <- paste(output.folder, "/",platePrefix, "_",cell.line, '_', this.treat, '_readCounts.txt', sep='')
	write.table(out.counts, file=dat.name, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
}

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
			##		sufficient coverage
			## het.ind ~ indicator for heterozygotes with p > .99
			## annotations ~ annotations filtered by coverage and het probability
			## q.thresh ~ q-value threshold for declaring signifigance
			## DEBUGGING
			#ii <- 1
			sample <- ii					
			this.sample <- sample.names[sample]
			coverage <- (ase.dat.final$ref[, sample] + ase.dat.final$alt[, sample])
			coverage.floor <- 5
			coverage.ind <- (coverage>coverage.floor)
			ref <- ase.dat.final$ref[coverage.ind, sample]
			alt <- ase.dat.final$alt[coverage.ind, sample]
			phi <- ase.dat.final$gmat[coverage.ind]
			eps <- ase.joint.eps[sample]
			het <- ase.joint$gt[coverage.ind, 2]
			het.ind <- (het > 0.99)
			numb.hets <- sum(het.ind)
			annotations <- ase.dat.final$annotations[coverage.ind, ][het.ind, ]
			q.thresh <- 0.2
			this.treatment <- treatmentIDs_final[sample]
			
			cat("============================================================\n", sep="")
			cat("==========Processing Sample: ", this.sample, "==========\n", sep="")		
			cat("==========Treatment: ", this.treatment, " ==========\n", sep="")
			cat("==========", numb.hets, " heterozygotes with [P(het)>.99]==========","\n", sep="")
			
			################################################################## 
			## rho ~ calculate rho using eps
			## rho0 ~ simple rho estimate
			## lrt ~ likelihood ratio test using the simple rho estimate
			## pval ~ pval of the llk ratio test using the chi-square approximation
			## qv ~ qvalues calucalted from the llkRT using the Storey method
			rho <- comp.rho(ref,alt,eps)
			rho0 <- plogis(log(ref)-log(alt))
			lrt <- lrtEpsRhoBinom(ref,alt,eps,rho)
			pval <- (1-pchisq(2*lrt,df=1))
			qv <- qvalue(pval,pi0.method="bootstrap")
			
			cat("==========Simple rho estimate LLRT with [Q<", q.thresh, "]: ", sum(qv$qv<q.thresh), "==========\n", sep="")	
			
			################################################################## 
			## D ~ a grid of possible dispersion values for the Beta-binomial model
			## aux ~ loglikelihood of the beta bionmial model across D values 
			##		 with null \rho value
			## Dmax ~ disperison which maximizes the llk
			D <- exp((0:500)/50)
			aux <- sapply(D,function(D){
						sum(logLikBetaBinomialRhoEps(0.5,eps,D,ref[het.ind],alt[het.ind]))
					})
			Dmax <- D[which.max(aux)]
			
			cat("==========Dispersion estimate: ", round(Dmax, 3), "==========\n", sep="")	
			
			################################################################## 
			## Find MLE for \rho using the the Beta-Binomial model 
			## Dmax2 ~ the dispersian parameter estimed from the llk in the 
			##		   previous step
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
			Dmax2 <- Dmax
			#auxLogis <- optim(rep(0,sum(het.ind)),fn=logLikBetaBinomial2,
			#		gr=gLogLikBetaBinomial,D=Dmax2,R=ref[het.ind],A=alt[het.ind],
			#		method="L-BFGS-B", hessian=TRUE)

                        ##
                        ## unbounded optimization
                        ##
                        aux2 <- t(sapply(1:sum(het.ind),function(ii){
                           auxLogis <- optim(0,
                                             fn=logLikBetaBinomial2,
                                             gr=gLogLikBetaBinomial,
                                             D=Dmax2,
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
                        #betas.se <- aux2[,2]
			#rho3 <- plogis(auxLogis$par)
			lrt3 <- logLikBetaBinomialRhoEps(rho3,eps,Dmax2,ref[het.ind],alt[het.ind]) - 
					logLikBetaBinomialRhoEps(0.5,eps,Dmax,ref[het.ind],alt[het.ind])
			pval3 <- (1-pchisq(2*lrt3,df=1))
			betas.se <- abs(betas.beta.binom/qnorm(pval3/2))
			betas.se[which(betas.se=='NaN')] <- aux2[, 2][which(betas.se=='NaN')]
			qv3 <- qvalue(pval3, pi0.method='bootstrap')
			qvals.qv3 <- qv3$qv
                        ##
                        ##
                        
                        
                        
			##
			# output the genes which are differentially expressed
			##
			if(sum(qvals.qv3<0.1) > 0){
				bed <- data.frame(annotations[which(qvals.qv3<0.1),],strand="+")
				tmpFile <- scan(pipe("mktemp -t"),character(0))
				write.table(bed,file=tmpFile,sep="\t",quote=F,row.names=F,col.names=F)
				# should be using /wsu/home/groups/piquelab/data/SNPs/1KG_SNPs_filt.011614.bed.gz but no genes
				geneAnnoFile <- "/wsu/home/groups/piquelab/data/RefTranscriptome/ensGene.hg19.v2.bed.gz"
				command <- paste("intersectBed -a ",tmpFile," -b ",geneAnnoFile," -split -wao",sep="")
				aa <- read.table(pipe(command),sep="\t",as.is=T,na.strings=".")
				filename <- paste('./output/', platePrefix, '_', cell.line, '_', this.treatment, '_ensg.txt', sep='')
				write.table(data.frame(unique(aa$V19)), file=filename, row.names=FALSE, col.names=FALSE, quote=FALSE)
				filename <- paste('./output/', platePrefix, '_', cell.line, '_', this.treatment, '_gene.txt', sep='')
				write.table(data.frame(unique(aa$V20)), file=filename, row.names=FALSE, col.names=FALSE, quote=FALSE)
			}
                        
			#betas.beta.binom <- auxLogis$par
			#betas.se <- 1/(diag(auxLogis$hessian)^.5)
			betas.z <- betas.beta.binom/betas.se
			betas.pval <- 2*pnorm(-abs(betas.z))
			this.pi0 <- qv3$pi0
			qv.05 <- sum(qvals.qv3<0.05)
			qv.1 <- sum(qvals.qv3<0.1)
			qv.2 <- sum(qvals.qv3<0.2)
			qv.5 <- sum(qvals.qv3<0.5)
			
			cat("Loci ", annotations$rsID[which(qvals.qv3<q.thresh)], "\n", sep=" ")
			cat("==========Beta-Binomial LLRT with [Q<", q.thresh, "]: ", qv.2," Pi0=", round(this.pi0, 3),  "==========\n", sep="")

                        ## output complete information for every heterozygote
                        #complete.dat <- annotations[,-5]  ## remove allele frequency
                        complete.dat <- annotations
                        complete.dat$beta <- betas.beta.binom 
                        complete.dat$beta.se <- betas.se 
                        complete.dat$pval <- pval3
                        complete.dat$qval <- qvals.qv3 
                        filename.all <- paste('./output/', platePrefix, '_', cell.line, '_', this.treatment, '_allOutput.txt', sep='')
                        write.table(complete.dat, file=filename.all, row.names=FALSE, col.name=TRUE, quote=FALSE)
                        system(
                          paste0("less ",
                                 filename.all,
                                 " | tr ' ' '\t' | grep -v 'pos0' | intersectBed -wo -a stdin -b ~/piquelab/data/RefTranscriptome/ensGene.hg19.v2.bed.gz -split | cut -f1-3,5,7-10,14,23-24 - > ",
                                 paste0(filename.all, ".tid")))
                        
			################################################################## 
			## rho2 ~ reassign rho calculated with a simple estimate from epsilon 
			## rho2[het.ind] ~ heterozygotes are re-assigned rho values
			##      		   calucalted from grid optimization above
			## aux ~ choose the null model whith largest probability  
			## lrt2 ~ llkRT with Beta-bionmial and the possibility that 
			##        that the heterozygote is of a different genotype
			## pval2 ~ pval of thr llk ratio test using the chi-square approximation
			##         (includes uncertainty in genotyping)
			## qv2 ~ qvalue object calucalted from the llkRT 
			## qv2.qvals ~ qvalues from the previous calculation
			rho2 <- rho
			rho2[het.ind] <- rho3
			aux <- pmax(logLikBetaBinomialRhoEps(0.0,eps,Dmax,ref,alt),
					logLikBetaBinomialRhoEps(1.0,eps,Dmax,ref,alt),
					logLikBetaBinomialRhoEps(0.5,eps,Dmax,ref,alt))
			lrt2 <- logLikBetaBinomialRhoEps(rho2,eps,Dmax2,ref,alt) - aux;
			pval2 <- (1-pchisq(2*lrt2,df=1))
			qv2 <- qvalue(pval2[het.ind],pi0.method="bootstrap")
			qvals.qv2 <- qv2$qv
			#qv.05 <- sum(qvals.qv2<0.05)
			#qv.1 <- sum(qvals.qv2<0.1)
			#qv.2 <- sum(qvals.qv2<0.2)
			#qv.5 <- sum(qvals.qv2<0.5)
                        #this.pi0 <- qv2$pi0
                        #betas.se <- abs(betas.beta.binom/qnorm(pval2[het.ind]/2))
			#betas.se[which(betas.se=='NaN')] <- aux2[, 2][which(betas.se=='NaN')]
			                        
			#100-qv2$pi0*100
			cat("==========LLRT considering uncertain GT with [Q<", q.thresh, "]: ", sum(qv2$qv<q.thresh), "==========\n\n", sep="")		
			
			##################################################################
			## return data frame
			rsID <- annotations
			betas <- betas.beta.binom
                        dat <- data.frame(annotations$rsID, annotations$chr, annotations$pos0, rho=rho3, betas, betas.se, qv=qvals.qv3, pval=pval3, refCount=ref[het.ind], altCount=alt[het.ind])
                        meta.dat <- data.frame(hets=numb.hets, pi0=this.pi0, qv.05=qv.05, qv.1=qv.1, qv.2=qv.2, qv.5=qv.5, mean.rho=mean(rho3), median.rho=median(rho3), dispersion=Dmax2)
			temp <- list(dat=dat, meta.dat=meta.dat)
                        
}) ## Returns a list of data & metaData

names(inference.data) <- treatmentIDs_final
dat.name <- paste(output.folder, "/",platePrefix, "_",cell.line, "_cov", min.cov, '_inference.RData', sep='')
save(inference.data, file=dat.name)
str(inference.data)
#load(dat.name)

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

rho_title <- paste0(platePrefix, '-', cell.line, ' : Rho_hat across all treatements', '\n mean.Rho=', mean_rho_hat, ' | median.Rho=', median_rho_hat)
pdf.file <- paste('./plots/QC/', platePrefix, '_', cell.line, '_averageRho', '.pdf', sep='')
pdf(file=pdf.file)
hist(all_rho_hat[all_coverage>100], breaks=80, main=rho_title, xlim=c(0,1), breaks=seq(0,1,0.01), color='darkgrey', axes=FALSE)
abline(v=mean_rho_hat, lty=2, col='red')
axis(1, at=seq(0, 1, .1)); axis(2)
dev.off()

##########################################
## 1.) QQ-plots for all treatments
##########################################
for(ii in seq_along(1:length(inference.data))){
        #ii <- 1
        treatment <- names(inference.data)[ii]
	pvals <- inference.data[[ii]]$dat$pval
	pi0 <- round(inference.data[[ii]]$meta.dat$pi0, 2)
	hets <- inference.data[[ii]]$meta.dat$hets
	qv.2 <- inference.data[[ii]]$meta.dat$qv.2
        coverage <- inference.data[[ii]]$dat$refCount + inference.data[[ii]]$dat$altCount
	avg.depth <- floor(mean(coverage))
        disp <- round(mean(inference.data[[ii]]$meta.dat$dispersion), 2)

        ## extract pvalues from only the high coverage loci
        pval_high <- inference.data[[1]]$dat$pval[which(coverage>100)]
        qqp <- qqplot(-log10(ppoints(length(pval_high))),-log10(pval_high), plot.it=F)
         
        pdf.file <- paste('./plots/QQ/', platePrefix, '_', cell.line, '_', treatment, "_cov", min.cov, '_', ii, '_QQ', '.pdf', sep='')
        title <- paste(cell.line, " | ", treatments_collapsed[ii], ' | Pi0=', pi0, ' | #hets=', hets, '\n #qv.2=', qv.2, ' | avg.depth=', avg.depth, ' | disp=', disp, sep='')
       
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
                        #ii <- 1 
			sapply(c(.01, .05, .1, .2), FUN=function(jj){sum(inference.data[[ii]]$dat$qv < jj)})
	}))

rownames(asetable) <- names(inference.data)
colnames(asetable) <- c('Q<.01', 'Q<.05', 'Q<.1', 'Q<.2')

outfile <- paste('./output/', platePrefix, "_",cell.line, '_Qhits.txt', sep='')
write.table(asetable, file=outfile, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t") 

##                          ##
cat("###### THE END ######")
##                          ##
