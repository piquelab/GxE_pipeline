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

cargs <- commandArgs(trail=TRUE)
if( length(cargs) >= 1) {
    meshinput <- cargs[1]
}

###########################################################################
system('mkdir -p plots')

###########################################################################
## lets grab gene id information so we can annotate significant SNPs
###########################################################################
ddgene <- read.table(meshinput, header=TRUE, stringsAsFactors=FALSE)[, c('snp', 'ensg')]

###########################################################################
## snp | beta.con | se.con | beta.treat | se.treat | bf1 | bf2 | bf3 | post4
dd <- read.table('./posteriors_bfs/Master_table_betas_bfs.txt', header=TRUE, stringsAsFactors=FALSE)

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
  strsplit(snp, ";")[[1]]
})); colnames(snps_treat) <- c('plate', 'line', 'rsID', 'treatment')

plotdat_treat <- data.frame(snp=ddtreat[, 'snp'], snps_treat, ddtreat[, -1], stringsAsFactors=FALSE)
plotdat_treat <- merge(plotdat_treat, ddgene, by='snp')

## forest plot for treatment configurations
pd <- position_dodge(width=0.5) #, height=NULL) # 'height' not in 1.0.1
p <- ggplot(plotdat_treat, aes(x=snp, y=beta, colour=treat))
pdf('./plots/Mesh_allBetas_treatment_config.pdf', width=9, height=12)
p + geom_abline(intercept=0, slope=0, colour='red',linetype=4) +
    ##layer(geom="point", position=pd, geom_params=list(size=3.2)) +
    geom_point(size=3.2, position=pd) +
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
  strsplit(snp, ";")[[1]]
})); colnames(snps_control) <- c('plate', 'line', 'rsID', 'treatment')

plotdat_control <- data.frame(snp=ddcontrol[, 'snp'], snps_control, ddcontrol[, -1], stringsAsFactors=FALSE)
plotdat_control <- merge(plotdat_control, ddgene, by='snp')

## forest plot for treatment configurations
pd <- position_dodge(width=0.5) #, height=NULL)
p <- ggplot(plotdat_control, aes(x=snp, y=beta, colour=treat))
pdf('./plots/Mesh_allBetas_control_config.pdf', width=9, height=12)
p + geom_abline(intercept=0, slope=0, colour='red',linetype=4) +
    ##layer(geom="point", position=pd, geom_params=list(size=3.2)) +
    geom_point(position=pd, size=3.2) +
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

#############################################################################################################################
dd$z.con <- dd$beta.c/dd$se.c
dd$z.treat <- dd$beta.t/dd$se.t
dd$z.ind <- ((dd$z.treat^2 + dd$z.con^2)>4)

dd_filt <- dd[dd$post_c4<.5, ]
controlAse   <- (dd_filt$bf_c1 - dd_filt$bf_c3) > 1 & dd_filt$bf_c1 > 2
treatmentAse <- (dd_filt$bf_c2 - dd_filt$bf_c3) > 1 & dd_filt$bf_c2 > 2
dd_filt$ASE <- '~'
dd_filt$ASE[ controlAse ] <- 'control'
dd_filt$ASE[ treatmentAse ] <- 'treatment'

bf.cutoff <- 2

system('mkdir -p plots')

##
pdf('./plots/Mesh_bf4_posteriorProb_hist.pdf')
ggplot(dd, aes(x=post_c4)) + geom_histogram(binwidth=0.01) +
  labs(title="Posterior of BF4 per SNP") +
  geom_vline(xintercept=median(dd$post_c4), colour="red", size=0.7, lty=2) +
  annotate("text", x=.75, y=150000, label="median post(BF4)=.96", colour='red')
dev.off()

## all ASE z-scores coloured by BF
pdf('./plots/Mesh_zscores_allASE_stringent_cutOff.5.pdf')
#png('./plots/Mesh_zscores_allASE_stringent_cutOff.5.png')
p <- ggplot(dd_filt, aes(z.con, z.treat, color=ASE))
p + #layer(geom="point", geom_params=list(size=2)) +
  geom_point(size=2) + 
  scale_color_manual(values=c("transparent", "blue", "red")) + 
  theme_bw() +
  theme(legend.position="bottom",
        legend.direction="horizontal",
        legend.key=element_blank(),
        plot.title=element_text(angle=0,size=16,face="bold"),
        axis.title.x=element_text(color="black", size=14),
        axis.title.y=element_text(color="black", size=14)) +
  labs(title="Z-scores with Posterior(BF4)<0.5", x="Control Z score", y="Treatment Z score")
dev.off()


## bi-plot with z-scores coloured by BF4
pdf('./plots/Mesh_zscores_BF4.cutOff.5.pdf')
#png('./plots/Mesh_zscores_BF4.cutOff.5.png')
p <- ggplot(dd, aes(z.con, z.treat))
p + geom_point(size=2) + 
  theme_bw() +
  #scale_colour_gradient2(midpoint=0.5, high='white', low='red') +
  theme(legend.position="bottom",
        legend.direction="horizontal",
        legend.key=element_blank(),
        plot.title=element_text(angle=0,size=16,face="bold"),
        axis.title.x=element_text(color="black", size=14),
        axis.title.y=element_text(color="black", size=14)) +
  labs(title="Z-scores with Posterior(BF4)<0.5", x="Control Z score", y="Treatment Z score")
dev.off()

## control ASE
pdf('./plots/Mesh_zscores_controlASE.pdf')
p <- ggplot(dd[dd$z.ind, ], aes(z.con, z.treat, color=(bf_c1>bf.cutoff)))
p + #layer(geom="point", geom_params=list(size=2)) +
  geom_point(size=2) +
  scale_color_manual(values=c("black", "blue")) + 
  theme_bw() +
  theme(legend.position="bottom",
        legend.direction="horizontal",
        legend.key=element_blank(),
        plot.title=element_text(angle=0,size=16,face="bold"),
        axis.title.x=element_text(color="black", size=14),
        axis.title.y=element_text(color="black", size=14)) +
  labs(title="Control ASE", x="Control Z score", y="Treatment Z score", color="log10(BF.1)>2")
dev.off()

## treatment ASE
pdf('./plots/Mesh_zscores_treatmentASE.pdf')
#png('./plots/Mesh_zscores_treatmentASE.png')
p <- ggplot(dd[dd$z.ind, ], aes(z.con, z.treat, color=(bf_c2>bf.cutoff)))
p + #layer(geom="point", geom_params=list(size=2)) +
  geom_point(size=2) +
  scale_color_manual(values=c("black", "red")) +
  theme_bw() +
  theme(legend.position="bottom",
        legend.direction="horizontal",
        legend.key=element_blank(),
        plot.title=element_text(angle=0,size=16,face="bold"),
        axis.title.x=element_text(color="black", size=14),
        axis.title.y=element_text(color="black", size=14)) +
  labs(title="Treatment ASE", x="Control Z score", y="Treatment Z score", color="log10(BF.2)>2")
dev.off()

## more stringent cutoffs ##
## control ASE
pdf('./plots/Mesh_zscores_controlASE_stringent1.pdf')
#png('./plots/Mesh_zscores_controlASE_stringent1.png')
p <- ggplot(dd[dd$z.ind, ], aes(z.con, z.treat, color=((bf_c1-bf_c3)>1 & bf_c1>2)))
p + #layer(geom="point", geom_params=list(size=2)) +
  geom_point(size=2) +
  scale_color_manual(values=c("black", "blue")) + 
  theme_bw() +
  theme(legend.position="bottom",
        legend.direction="horizontal",
        legend.key=element_blank(),
        plot.title=element_text(angle=0,size=16,face="bold"),
        axis.title.x=element_text(color="black", size=14),
        axis.title.y=element_text(color="black", size=14)) +
  labs(title="Control ASE", x="Control Z score", y="Treatment Z score", color="(log10(bf1)-log10(bf3))>1 & log10(bf1)>2)")
dev.off()

## more stringent cutoffs ##
## treatment ASE
pdf('./plots/Mesh_zscores_treatmentASE_stringent1.pdf')
#png('./plots/Mesh_zscores_treatmentASE_stringent1.png')
p <- ggplot(dd[dd$z.ind, ], aes(z.con, z.treat, color=((bf_c2-bf_c3)>1 & bf_c2>2)))
p + #layer(geom="point", geom_params=list(size=2)) +
  geom_point(size=2) + 
  scale_color_manual(values=c("black", "red")) + 
  theme_bw() +
  theme(legend.position="bottom",
        legend.direction="horizontal",
        legend.key=element_blank(),
        plot.title=element_text(angle=0,size=16,face="bold"),
        axis.title.x=element_text(color="black", size=14),
        axis.title.y=element_text(color="black", size=14)) +
  labs(title="Treatment ASE", x="Control Z score", y="Treatment Z score", color="(log10(bf2)-log10(bf3))>1 & log10(bf2)>2)")
dev.off()

## control only ASE
pdf('./plots/Mesh_zscores_controlOnlyASE.pdf')
#png('./plots/Mesh_zscores_controlOnlyASE.png')
p <- ggplot(dd[dd$z.ind, ], aes(z.con, z.treat, color=(bf_c1>bf.cutoff & bf_c3<bf.cutoff)))
p + #layer(geom="point", geom_params=list(size=2)) +
  geom_point(size=2) + 
  scale_color_manual(values=c("black", "blue")) +
  theme_bw() +
  theme(legend.position="bottom",
        legend.direction="horizontal",
        legend.key=element_blank(),
        plot.title=element_text(angle=0,size=16,face="bold"),
        axis.title.x=element_text(color="black", size=14),
        axis.title.y=element_text(color="black", size=14)) +
  labs(title="Control Only ASE", x="Control Z score", y="Treatment Z score", color="log10(BF.1)>2 & log10(BF.3)<2")
dev.off()

## treatment only ASE
pdf('./plots/Mesh_zscores_treatmentOnlyASE.pdf')
#png('./plots/Mesh_zscores_treatmentOnlyASE.png')
p <- ggplot(dd[dd$z.ind, ], aes(z.con, z.treat, color=(bf_c2>bf.cutoff & bf_c3<bf.cutoff)))
p + #layer(geom="point", geom_params=list(size=2)) +
  geom_point(size=2) +
  scale_color_manual(values=c("black", "red")) +
  theme_bw() +
  theme(legend.position="bottom",
        legend.direction="horizontal",
        legend.key=element_blank(),
        plot.title=element_text(angle=0,size=16,face="bold"),
        axis.title.x=element_text(color="black", size=14),
        axis.title.y=element_text(color="black", size=14)) +
  labs(title="Treatment Only ASE", x="Control Z score", y="Treatment Z score",  color="log10(BF.2)>2 & log10(BF.3)<2")
dev.off()


##         ##
## THE END ##
##         ##
