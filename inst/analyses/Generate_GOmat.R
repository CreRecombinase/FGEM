##6/12/15
##NWK
                                        #Generate annotation matrix for use with FGEM
library(rhdf5)
library(reshape2)
library(topGO)
library(Matrix)

# args <- commandArgs(trailingOnly=T)

# BFile <- args[1]
# outfile <- args[2]
setwd("~/Dropbox/BayesianDA/FGEM_Data/")
BFile <- "TADA_ASC_SSC_results_Dec23.csv"
outfile <- "GOmat.RDS"
h5file <- "/home/nwknoblauch/Dropbox/BayesianDA/FGEM/Data/Annomat.h5"
h5createFile(h5file)
h5createGroup(h5file,"Annotations")
h5createGroup(h5file,"Annotations/GO")
h5createGroup(h5file,"BF")
h5createGroup(h5file,"GeneList")


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
GL <- rownames(GOmat)
h5write(GL,h5file,"GeneList/GL")
fid <- H5Fopen(h5file)
did <- H5Dopen(fid,"Annotations/GO/GOmat")
h5write(1:nrow(GOmat),h5file,"GeneList/Gi")
h5write(colnames(GOmat),h5file,"Annotations/GO/GOterms")
h5write(1:ncol(GOmat),h5file,"Annotations/GO/Goi")
GOsm <- data.matrix(GO.G[,c("Gi","GOi")])
h5write(GOsm,h5file,"Annotations/GO/GOmat")
h5readAttributes(h5file,"Annotations/GO/GOmat")
h5writeAttribute.integer(1,h5obj=did,name="sparseMat")
h5write(BF.G$BF,h5file,"BF/BayesFactors")
H5Dclose(did)
H5Fclose(fid)

saveRDS(GOmat,file=outfile)
