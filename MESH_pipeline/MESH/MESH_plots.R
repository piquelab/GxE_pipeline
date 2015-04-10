###########################################################################
##
## CTH 111114
## bi-plots for z-scores coloured by Bayes Factor
##
###########################################################################
library(ggplot2)
library(reshape2)
source('~/piquelab/charvey/GxE/jointGenotyping/scripts/misc/mesh_analysis_functions.R')
source('/wsu/home/groups/piquelab/charvey/tools/R_tools/myParallel.R')

## x11(display="localhost:11.0" ,type="Xlib")

###########################################################################
system('mkdir -p plots')

###########################################################################
## lets grab gene id information so we can annotate significant SNPs
###########################################################################
ddgene <- read.table("MESH_QuASAR_master_logFC_controlTreat_correctlfc.txt", header=TRUE, stringsAsFactors=FALSE)[, c('snp', 'g.id')]

###########################################################################
## snp | beta.con | se.con | beta.treat | se.treat | bf1 | bf2 | bf3 | post4
dd <- read.table('./posteriors_bfs/MESH_masterTable_betas_bfs.txt', header=TRUE, stringsAsFactors=FALSE)

###########################################################################
## forest plot of betas for treatment only ASE
###########################################################################
betaTreat_dat <- read.table(pipe("find ./posteriors_bfs -name \"*bayesFactors*.txt\" | while read bf; do less $bf | grep -v \"snp\" | awk \'($3>$4*2) && ($3>20) && ($2<0.5){print $0\"\t\"$3/$4}\' | sort -k6 -n | cut -f1 | tr \"_\" \"\t\" | cut -f1,2,3 | tr \"\t\" \"_\"| while read f; do grep $f allPlates*.txt | tr \" \" \"\t\"; done; done"), stringsAsFactors=FALSE)
names(betaTreat_dat) <- c('snp', 'treat', 'beta', 'beta.se')
ddtreat <- merge(betaTreat_dat, dd, by='snp')
ddtreat$lb <- ddtreat$beta - 1.96*ddtreat$beta.se
ddtreat$ub <- ddtreat$beta + 1.96*ddtreat$beta.se
ddtreat$treat[(seq(1, dim(ddtreat)[1], 2) + 1)] <- "treatment"
ddtreat$treat[(seq(1, dim(ddtreat)[1], 2))] <- "control"

snps_treat <- do.call(rbind, lapply(ddtreat$snp, function(snp){
  # snp <- ddtreat$snp[1]
  strsplit(snp, ":")[[1]]
})); colnames(snps_treat) <- c('plate', 'line', 'rsID', 'treatment')

plotdat_treat <- data.frame(snp=ddtreat[, 'snp'], snps_treat, ddtreat[, -1], stringsAsFactors=FALSE)
plotdat_treat <- merge(plotdat_treat, ddgene, by='snp')

## forest plot for treatment configurations
pd <- position_dodge(width=0.5, height=NULL)
p <- ggplot(plotdat_treat, aes(x=snp, y=beta, colour=treat))
pdf('./plots/Mesh_allBetas_treatment_config.pdf', width=9, height=12)
p + geom_abline(intercept=0, slope=0, colour='red',linetype=4) +
    layer(geom="point", position=pd, geom_params=list(size=3.2)) +
    geom_errorbar(aes(ymin=lb , ymax=ub), width=0.5, size=1, position=pd) +
    theme_bw() +
    coord_flip() +
    scale_color_manual(values=c("black", "red")) +
    theme(legend.position="bottom",
                  legend.direction="horizontal",
                  legend.key=element_blank(),
                  plot.title=element_text(angle=0,size=16, face="bold"),
                  axis.title.x=element_text(size=14, face="bold"),
                  axis.text.x=element_text(size=12),
                  axis.text.y=element_text(size=12)) +
    labs(title="Estimates of Treatment Only ASE", y=expression(hat(beta)), x="")
dev.off()
write.table(plotdat_treat, file="./plots/Mesh_allBetas_treatment_config.txt", col.names=TRUE, row.names=FALSE, quote=FALSE)

###########################################################################
## forest plot of betas for control only ASE
###########################################################################
betaControl_dat <- read.table(pipe("find ./posteriors_bfs -name \"*bayesFactors*.txt\" | while read bf; do less $bf | grep -v \"snp\" | awk \'($2>$4*2) && ($2>20) && ($3<0.5){print $0\"\t\"$2/$4}\' | sort -k6 -n | cut -f1 | tr \"_\" \"\t\" | cut -f1,2,3 | tr \"\t\" \"_\"| while read f; do grep $f allPlates*.txt | tr \" \" \"\t\"; done; done"), stringsAsFactors=FALSE)
names(betaControl_dat) <- c('snp', 'treat', 'beta', 'beta.se')
ddcontrol <- merge(betaControl_dat, dd, by='snp')
ddcontrol$lb <- ddcontrol$beta - 1.96*ddcontrol$beta.se
ddcontrol$ub <- ddcontrol$beta + 1.96*ddcontrol$beta.se
ddcontrol$treat[(seq(1, dim(ddcontrol)[1], 2) + 1)] <- "treatment"
ddcontrol$treat[(seq(1, dim(ddcontrol)[1], 2))] <- "control"

snps_control <- do.call(rbind, lapply(ddcontrol$snp, function(snp){
  # snp <- ddtreat$snp[1]
  strsplit(snp, ":")[[1]]
})); colnames(snps_control) <- c('plate', 'line', 'rsID', 'treatment')

plotdat_control <- data.frame(snp=ddcontrol[, 'snp'], snps_control, ddcontrol[, -1], stringsAsFactors=FALSE)
plotdat_control <- merge(plotdat_control, ddgene, by='snp')

## forest plot for treatment configurations
pd <- position_dodge(width=0.5, height=NULL)
p <- ggplot(plotdat_control, aes(x=snp, y=beta, colour=treat))
pdf('./plots/Mesh_allBetas_control_config.pdf', width=9, height=12)
p + geom_abline(intercept=0, slope=0, colour='red',linetype=4) +
    layer(geom="point", position=pd, geom_params=list(size=3.2)) +
    geom_errorbar(aes(ymin=lb , ymax=ub), width=0.5, size=1, position=pd) +
    theme_bw() +
    coord_flip() +
    scale_color_manual(values=c("black", "red")) +
    theme(legend.position="bottom",
                  legend.direction="horizontal",
                  legend.key=element_blank(),
                  plot.title=element_text(angle=0,size=16, face="bold"),
                  axis.title.x=element_text(size=14, face="bold"),
                  axis.text.x=element_text(size=12),
                  axis.text.y=element_text(size=12)) +
    labs(title="Estimates of Control Only ASE", y=expression(hat(beta)), x="")
dev.off()
write.table(plotdat_control, file="./plots/Mesh_allBetas_control_config.txt", col.names=TRUE, row.names=FALSE, quote=FALSE)

##         ##
## THE END ##
##         ##
