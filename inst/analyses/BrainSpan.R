##Code to usefully integrate gene expression measurements for FGEM

library(edgeR)
setwd("~/Dropbox/BayesianDA/FGEM/Data/")
colmetaf <- "columns_metadata.csv"
rowmetaf <- "rows_metadata.csv"
dataf <- "expression_matrix.csv"
BayesFactorFile <- "TADA_ASC_SSC_results_Dec23.csv"

BF.G <- read.table(BayesFactorFile,header=T,sep=",")
colmeta <- read.table(colmetaf,sep=",",header=T,stringsAsFactors=F)
rowmeta <- read.table(rowmetaf,sep=",",header=T,stringsAsFactors=F)
expdata <- read.table(dataf,sep=",")

expdata <- data.matrix(expdata[,-1])

#Remove regions where we don't have at least 20 samples, and where the gene is not present in the BF file
multsamp <- names(table(colmeta$structure_name))[table(colmeta$structure_name)>20]

mexpdata <- expdata[,colmeta$structure_name %in% multsamp]
mexpdata <- mexpdata[rowmeta$gene_symbol %in% BF.G$Gene,]
rowmeta <- rowmeta[rowmeta$gene_symbol %in% BF.G$Gene,]

mcolmeta <- colmeta[ colmeta$structure_name %in% multsamp,]
colnames(mexpdata) <- mcolmeta$column_num
rownames(mexpdata) <- rowmeta$ensembl_gene_id

maleexp <- mexpdata[,mcolmeta$gender=="M"]
malecolmeta <- mcolmeta[mcolmeta$gender=="M",]

fmaleexp <- mexpdata[,mcolmeta$gender!="M"]
fmalecolmeta <- mcolmeta[mcolmeta$gender!="M",]

mcolmetasub <- mcolmeta[,c("structure_name"),drop=F]
texp <- t(mexpdata)
expdf <- data.frame(mcolmetasub,texp)

expv <- daply(expdf,c("structure_name"),colwise(sd),.progress = "text")



mexpv <- matrix(unlist(expv,recursive = T),dim(expv)[1],dim(expv)[2],dimnames = dimnames(expv))
dexpv <- data.frame(t(mexpv))
dexpv$gid <- rowmeta$gene_symbol
dexpv <- daply(dexpv,"gid",colwise(max))
mexpv <- matrix(unlist(dexpv,recursive = T),nrow(dexpv),ncol(dexpv),dimnames=dimnames(dexpv))
mexpv <- mexpv[rownames(mexpv)%in% BF.G$Gene,]
saveRDS(mexpv,"BrainSpanVar.RDS")


of <- lm(Y~X+0)
B <- coef(of)
sum((Y-B*X)^2)
sum(Y^2-(mean(X*Y)^2/mean(X^2)))
sum(Y^2)-mean(X*Y)^2/mean(X^2)
sum((X-mean(X))*Y)/sum((X-mean(X))^2)
sum((Y^2-mean(X*Y)^2/mean(X^2)))
