#########################################################################
##
## GMB 08/19/15
## Plot expression data
## 
## INPUT/ARGS:
##   exprFile - FPKM data to use for plots
##   
#########################################################################

## Libraries ##
library(ggplot2)
library(gplots)

## Get command-line arguments ##
cargs<-commandArgs(trail=T)
if (length(cargs)>0){
  exprFile=cargs[1]}

platePrefix = gsub('.*/', '', gsub('\\..*', '', exprFile)) 

## Get data
load(exprFile) # loads 'fpkms' and 'cv2'

## Load the treatment key
treatmentKey <- read.table("../treatmentKey.txt", as.is=T, sep='\t',
                           header=T, comment.char="")

## PCA
corMtx <- cor(fpkms) #cor(cdataAdj)
pca <- princomp(corMtx, cor=T)
var <- pca$sdev^2
var.p <- round(var/sum(var) * 100, 2)
PC1 <- pca$scores[,1]
PC2 <- pca$scores[,2]
PC3 <- pca$scores[,3]
PC4 <- pca$scores[,4]

cell.lines <- gsub('\\..*', '', colnames(fpkms))
treatments <- gsub('.*\\.', '', colnames(fpkms))
cols = treatmentKey[unique(treatments)[order(unique(treatments))],]$DeepColor

## PC1 v PC2
fname=paste0(platePrefix, '/', platePrefix, '_pc1pc2.pdf')
pdf(fname)
dat <- data.frame(x=PC1, y=PC2, tx=treatments, cl=cell.lines)
p1 <- ggplot(dat, aes(x=x, y=y, col=tx, shape=cl))
p1 + 
  geom_point(size=4) +
  xlab(paste0("PC1 - ", (var.p[1]), "%")) + 
  ylab(paste0("PC2 - ", (var.p[2]), "%")) + 
  labs(shape="Cell Line", color="Treatment") + 
  scale_color_manual(values=cols) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    )
dev.off()

## PC2 v PC3
fname=paste0(platePrefix, '/', platePrefix, '_pc2pc3.pdf')
pdf(fname)
dat <- data.frame(x=PC2, y=PC3, tx=treatments, cl=cell.lines)
p1 <- ggplot(dat, aes(x=x, y=y, col=tx, shape=cl))
p1 + 
	geom_point(size=4) +
	xlab(paste0("PC2 - ", (var.p[2]), "%")) + 
	ylab(paste0("PC3 - ", (var.p[3]), "%")) + 
	labs(shape="Cell Line", color="Treatment") + 
	scale_color_manual(values=cols) +
	theme_bw() +
	theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
          )
dev.off()

## PC3 v PC4
fname=paste0(platePrefix, '/', platePrefix, '_pc3pc4.pdf')
pdf(fname)
dat <- data.frame(x=PC3, y=PC4, tx=treatments, cl=cell.lines)
p1 <- ggplot(dat, aes(x=x, y=y, col=tx, shape=cl))
p1 + 
	geom_point(size=4) +
	xlab(paste0("PC3 - ", (var.p[3]), "%")) + 
	ylab(paste0("PC4 - ", (var.p[4]), "%")) + 
	labs(shape="Cell Line", color="Treatment") + 
	scale_color_manual(values=cols) +
	theme_bw() +
	theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
          )
dev.off()


##############
## Heatmaps
## Todo - implement heatmap.3

## Color the columns by treatment and the rows by cell type
fname=paste0(platePrefix, '/', platePrefix, '.pdf')
pdf(fname)
#hc <- hclust(dist(corMtx),method="complete")
hc <- hclust(dist(t(fpkms)),method="complete")
n.cols <- dim(fpkms)[2]
heatmap.2(corMtx,
          Rowv = as.dendrogram(hc),
          Colv = as.dendrogram(hc), 
          trace = "none",
          col = redblue(n.cols), 
          ColSideColors = treatmentKey[treatments,]$DeepColor,
          labRow = cell.lines, 
          labCol = treatmentKey[ treatments, "Short_Name" ]
          )
dev.off()

## the end ##
