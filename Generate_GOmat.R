##6/12/15
##NWK
                                        #Generate annotation matrix for use with FGEM

library(reshape2)
library(topGO)
library(Matrix)

# args <- commandArgs(trailingOnly=T)

# BFile <- args[1]
# outfile <- args[2]
setwd("./Data/")
BFile <- "TADA_ASC_SSC_results_Dec23.csv"
outfile <- "GOmat.RDS"
BF.G <- read.table(BFile,header=T,sep=",")
BF.G$sig <- ifelse(BF.G$BF>10,1,0)


allgenes <- data.frame(Gene=BF.G$Gene,B=BF.G$BF,stringsAsFactors=F)
rownames(allgenes) <- allgenes$Gene
allgenes$Gi <- 1:nrow(allgenes)
genelist <- allgenes$Gene

genGOdf <- function(genelist,GOC){
  GO.G <- melt(annFUN.org(GOC,feasibleGenes=genelist,mapping="org.Hs.eg.db",ID="symbol"))
  colnames(GO.G) <- c("Gene","GOT")
  return(GO.G)
}
#  GOBP.G <- melt(annFUN.org("BP",feasibleGenes=BF.G$Gene,mapping="org.Hs.eg.db",ID="symbol"))
#colnames(GOBP.G) <- c("Gene","GOT")


GO.G <- rbind(genGOdf(as.character(genelist),"BP"),genGOdf(as.character(genelist),"MF"))
GOT <- data.frame(Gi=1:length(unique(GO.G$GOT)),GOT=unique(GO.G$GOT))
rownames(GOT) <- GOT$GOT
GO.G$Gi <- allgenes[as.character(GO.G$Gene),"Gi"]
GO.G$GOi <- GOT[GO.G$GOT,"Gi"]




GOmat <- sparseMatrix(i=GO.G$Gi,j = GO.G$GOi,x=1,dims=c(nrow(allgenes),max(GO.G$GOi)),
                       dimnames=list(as.character(allgenes$Gene),rownames(GOT)))


saveRDS(GOmat,file=outfile)
