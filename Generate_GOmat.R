##6/12/15
##NWK
                                        #Generate annotation matrix for use with FGEM

library(reshape2)
library(topGO)
library(Matrix)

# args <- commandArgs(trailingOnly=T)

# BFile <- args[1]
# outfile <- args[2]
setwd("~/Dropbox/BayesianDA/FGEM/Data")
BFile <- "TADA_ASC_SSC_results_Dec23.csv"
outfile <- "GOmat.RDS"
BF.G <- read.table(BFile,header=T,sep=",")
BF.G$sig <- ifelse(BF.G$BF>10,1,0)


allgenes <- data.frame(Gene=BF.G$Gene,B=BF.G$BF,stringsAsFactors=F)
rownames(allgenes) <- allgenes$Gene
allgenes$Gi <- as.integer(factor(allgenes$Gene))


GOMF.G <- melt(annFUN.org("MF",feasibleGenes=BF.G$Gene,mapping="org.Hs.eg.db",ID="symbol"))
colnames(GOMF.G) <- c("Gene","GOT")
GOBP.G <- melt(annFUN.org("BP",feasibleGenes=BF.G$Gene,mapping="org.Hs.eg.db",ID="symbol"))
colnames(GOBP.G) <- c("Gene","GOT")

GO.G <- rbind(GOMF.G,GOBP.G)
GO.G$Gi <- allgenes[GO.G$Gene,"Gi"]
GO.G$GOi <- as.integer(factor(GO.G$GOT))


GOmat <- sparseMatrix(i=GO.G$Gi,j = GO.G$GOi,x = rep(1,nrow(GO.G)),dimnames=list(as.character(allgenes$Gene),as.character(levels(factor(GO.G$GOT)))))
GOmat <- GOmat[,colSums(GOmat)>3]

B <- allgenes[rownames(GOmat),"B"]



APBeta <- numeric(ncol(GOmat))
Bmat <- matrix(B,nrow=length(B),ncol=1)
GOb <- Matrix(GOmat,forceCheck = T,sparse = F)
GOb <- matrix(GOb,nrow(GOb),ncol(GOb))
cB <- cor(GOb,B)^2
 
GOmat <- GOmat[,order(cB[,1],decreasing = T)]

saveRDS(GOmat,file=outfile)
