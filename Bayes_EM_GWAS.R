library(topGO)
library(reshape2)
library(Matrix)

args <- commandArgs(trailingOnly=T)
GOmatfile <- args[1]
Startterm <- as.integer(args[2])
Numterms <- as.integer(args[3])

Pi <- function(x)Reduce("*",x)
GOmat <- readRDS(GOmatfile)
BFile <- "~/Downloads/TADA_ASC_SSC_results_Dec23.csv"
BF.G <- read.table(BFile,header=T,sep=",")
BF.G$sig <- ifelse(BF.G$BF>10,1,0)
Z <- BF.G$sig

allgenes <- data.frame(Gene=BF.G$Gene,Z,B=BF.G$BF,stringsAsFactors=F)
rownames(allgenes) <- allgenes$Gene
allgenes$Gi <- as.integer(factor(allgenes$Gene))
# GOMF.G <- melt(annFUN.org("MF",feasibleGenes=BF.G$Gene,mapping="org.Hs.eg.db",ID="symbol"))
# colnames(GOMF.G) <- c("Gene","GOT")
# GOBP.G <- melt(annFUN.org("BP",feasibleGenes=BF.G$Gene,mapping="org.Hs.eg.db",ID="symbol"))
# colnames(GOBP.G) <- c("Gene","GOT")
# 
# GO.G <- rbind(GOMF.G,GOBP.G)
# GO.G$Gi <- allgenes[GO.G$Gene,"Gi"]
# GO.G$GOi <- as.integer(factor(GO.G$GOT))


# GOmat <- sparseMatrix(i=GO.G$Gi,j = GO.G$GOi,x = rep(1,nrow(GO.G)),dimnames=list(as.character(allgenes$Gene),as.character(levels(factor(GO.G$GOT)))))
# GOmat <- GOmat[,colSums(GOmat)>3]
# 
# 
# 

Z <- allgenes[rownames(GOmat),"Z"]
B <- allgenes[rownames(GOmat),"B"]



# APBeta <- numeric(ncol(GOmat))
# for(i in 1:ncol(GOmat)){
#   if(i%%100==0){
#     print(i)
#   }
#   APBeta[i] <- summary(lm(B~GOmat[,i]))[["coefficients"]][2,4]
# }
# 
# 
# GOmat <- GOmat[,order(APBeta)]
####Starting EM
Betadf <- data.frame(Intercept=numeric(0),Beta=numeric(0),GO=character(0),liknull=numeric(0),likalt=numeric(0),stringsAsFactors=F)

for(j in 1:50){
  print(j)
  x <- as.matrix(GOmat[,j,drop=F])
  xm <- cbind(rep(1,length(Z)),GOmat[,j])
  colnames(xm) <- c("Intercept",colnames(GOmat)[j])
  Bofit <- glm(Z~x,family=binomial())
  Beta <- coefficients(Bofit)
  for(i in 1:40){
    if(i%%10==0){
      print(i)
    }
    pvec <- 1/(1+exp(-(xm%*%Beta)))
    uvec <- (pvec*B)/((pvec*B)+(1-pvec))
    Beta <- coefficients(glm(uvec~xm-1,family=quasibinomial(link="logit")))
  }
  nullBeta <- c(Beta[1],0)
  altpvec <- 1/(1+exp(-(xm%*%Beta)))
  altlikvec <- altpvec*B+(1-altpvec)
  likalt <- Pi(altlikvec)
  nullpvec <-  1/(1+exp(-(xm%*%nullBeta)))
  nulllvec <- nullpvec*B+(1-nullpvec)
  liknull <- Pi(nulllvec)
  Betadf[j,c("Intercept","Beta")]=Beta
  Betadf[j,"GO"] <- colnames(GOmat)[j]
  Betadf[j,"liknull"] <- liknull
  Betadf[j,"likalt"]<- likalt
}

Betadf$chisq<- -2*log(Betadf$liknull/Betadf$likalt)

Betadf$p <- dchisq(Betadf$chisq,df = 1)



