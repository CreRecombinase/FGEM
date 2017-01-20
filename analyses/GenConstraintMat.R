##6/12/15
##NWK
#Generate annotation matrix for use with FGEM

library(reshape2)
library(topGO)
library(Matrix)
setwd("~/Dropbox/BayesianDA/FGEM/Data/")
BFile <- "TADA_ASC_SSC_results_Dec23.csv"
ConstraintFile <- "forweb_cleaned_exac_r03_2015_03_16_z_data.txt"
outfile <- "ConstraintMat.RDS"
args <- commandArgs(trailingOnly=T)

BFile <- args[1]
ConstraintFile <- args[2]
outfile <- args[3]


BF.G <- read.table(BFile,header=T,sep=",")
CF.df <- read.table(ConstraintFile,header=T,sep="\t",stringsAsFactors = F)
rownames(CF.df) <- CF.df$gene
CF.mat <- data.matrix(CF.df[,-c(1:3)])

BF.G$sig <- ifelse(BF.G$BF>10,1,0)


allgenes <- data.frame(Gene=BF.G$Gene,B=BF.G$BF,stringsAsFactors=F)
rownames(allgenes) <- allgenes$Gene
allgenes$Gi <- as.integer(factor(allgenes$Gene))
CF.mat <- CF.mat[rownames(CF.mat) %in% rownames(allgenes),]



saveRDS(CF.mat,file=outfile)
