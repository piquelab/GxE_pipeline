###########################################################################
## multiclasslr_batch.R
## 
## This package implements a multiclass (softmax regression)
## Outcomes are one of four configurations from the MESH model & covariates are indicators for treatments &
## indicators for plates. Recovering the baseline predictor for each set of indicators is not possible from
## a single intercept. Instead we will try the following approach
## 1.) fit the model with an intercept only for each of k-1 classes
## 2.) use this intercept as a parameter in the final model
## 3.) 
##
## CTH 020815
##
##
###########################################################################
library(devtools)
library(MASS)
library(ggplot2)
#install_github('piquelab/multiclassLR')
#install_github('ctharve/multiclassLR')
library('multiclassLR')
#source("./TestFunctions_cf1cf2.R")
source("~/piquelab/charvey/source/multiclassLR/R/TestFunctions.R")

 cargs <- commandArgs(trail=TRUE);
if(length(cargs)>=1)
  dat.file <- cargs[1]
if(length(cargs)>=2)
  tag <- cargs[2]

# dat.file <- "masterTable_multinomial.txt"
# tag <- "2.0"
  
## directory for all plots
plot.dir <- "multinomialPlots"
system(paste0('mkdir -p ', plot.dir))

## at this iteration we will run snps that have all been run through MESH together
dd <- read.table(dat.file, stringsAsFactors=FALSE, header=TRUE)

## extract class probabilities
classlabs <- c('post_c1', 'post_c2', 'post_c3', 'post_c4')
#classlabs <- c('post1', 'post2', 'post3', 'post4')
Y <- dd[, classlabs]
classes <- dim(Y)[2] - 1

## prep data for model
#X <- dd[, !(colnames(dd) %in% classlabs)]
X <- dd[, !(colnames(dd) %in% c('snp', 'intercept', classlabs))]
n.obs <- dim(X)[1]

###############################
## 1.) run an intercept only model (.io) to estimate the global mean of each outcome
###############################
#X.io <- matrix(X[, 'intercept'])
X.io <- matrix(rep(1, n.obs))
Y.io <- Y
n.params.io <- dim(X.io)[2]
#classes.io <- dim(Y.io)[2] - 1
## remove another class to combine cf1 & cf2
classes.io <- dim(Y.io)[2] - 2
initweights.io <- rep(0, n.params.io*classes.io)
var.names.io <- c('inter.cf1cf2', 'inter.cf3')
multifinal.io <- MyMultinomial(initweights.io, X.io, Y.io, var.names.io, io.cf1cf2=TRUE)
head(multifinal.io$results$dat.sig)
io.betas <- multifinal.io$results$dat.all$parms
names(io.betas) <- rownames(multifinal.io$results$dat.sig)
## TODO change the naming within the function ##

###############################
## 2.) design matrix with all categorical predictors & no intercept 
###############################
## extract class probabilities. ensure the baseline class is the last class 
#X.ni <- X[, -which(colnames(X)=='neg')]
X.ni <- X
Y.ni <- Y

## initialize weights & create var names
initweights <- InitialWeights(X.ni, classes.io)
configurations <- c('_cf1cf2', '_cf3') ## should probably use class labs here
var.names.ni <- VarNames(X.ni, configurations)

## run model with initial intercepts
## we need to check all parameters before running
multifinal.ni <- MyMultinomial.ni(intercepts=io.betas, parms=initweights, X.ni, Y.ni, var.names.ni, cf1cf2=TRUE, maxit=1E5)

save(multifinal.ni, file=paste0("multifinal.ni_", tag, "_cf1cf2.RData"))

plot_dat <- CreatePlotDat(multifinal.ni)
plot.title <- paste0("multinomial_parms_", tag, "_cf1cf2.pdf")
PlotMultiForest(plot_dat, plot.title, plot.dir)

plot_dat <- CreatePlotDat(multifinal.ni, sig.parms=TRUE)
plot.title <- paste0("multinomial_parms_sig_", tag, "_cf1cf2.pdf")
PlotMultiForest(plot_dat, plot.title, plot.dir)

## output table
sigout <- data.frame(format(multifinal.ni$results$dat.sig, digits=5, scientific=TRUE))
write.table(sigout, file=paste0("multinomial_parms_sig_", tag, "_cf1cf2.txt"), row.names=TRUE, quote=FALSE, sep='\t')

##         ##
## THE END ##
##         ##
