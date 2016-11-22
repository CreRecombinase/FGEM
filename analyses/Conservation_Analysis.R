##Analysis of Conservation Features (And GO features)
library(plyr)
library(Matrix)

setwd("~/Dropbox/BayesianDA/FGEM/Data/")
BFile <- "TADA_ASC_SSC_results_Dec23.csv"
fmat  <- readRDS("ConstraintMat.RDS")
gmat <- readRDS("GOmat.RDS")
allg <- intersect(rownames(gmat),rownames(fmat))
bfmat <- fmat[allg,]
bgmat <- gmat[allg,]

BF.G <- read.table(BFile,header=T,sep=",")
BF.G <- BF.G[BF.G$Gene %in% rownames(bfmat),]
B <- BF.G$BF




setwd("~/Dropbox/BayesianDA/FGEM/Data/Conservation/")
fgem.res.files <- dir()
fgem.res <- ldply(fgem.res.files,readRDS)
fgem.res$pval <- pchisq(fgem.res$Chisq,df=1,lower.tail=F)
fgem.res <- fgem.res[order(fgem.res$pval),]
fgem.res$adjp <- p.adjust(fgem.res$pval,method = "fdr")

sig.fgem.res.cons <- fgem.res[fgem.res$pval<0.05,]
bfmat.sig <- bfmat[,colnames(bfmat)%in% sig.fgem.res.cons$X1]

setwd("~/Dropbox/BayesianDA/FGEM/Data/FGEM_GO_RDS/")
fgem.res.files <- dir()
fgem.res <- ldply(fgem.res.files,readRDS)
fgem.res$pval <- pchisq(fgem.res$Chisq,df=1,lower.tail=F)
fgem.res <- fgem.res[order(fgem.res$pval),]
fgem.res$adjp <- p.adjust(fgem.res$pval,method = "fdr")
sig.fgem.res.go <- fgem.res[fgem.res$pval<0.05,]
bgmat.sig <- bgmat[,colnames(bgmat) %in% sig.fgem.res.go$X1]
bsigmat <- cbind(as.matrix(bgmat.sig),bfmat.sig)
scalebsigm <- scale(bsigmat)




dfrep <- ldply(repfgem,data.frame)
hist(-2*(dfrep$NullLogLikelihood-dfrep$LogLikelihood))
hist(rchisq(n=60,df = 1))

BFCons2 <- FGEM(scalebsigm,B,iters=30)
BFCons <- BFCons2
Beta <- t(t(BFCons$Beta))
mx <- cbind(1,scalebsigm[,rownames(BFCons[-1,])])
BFCons$adjp <- p.adjust(BFCons$p,method="fdr")
sigBFCons <- BFCons[ BFCons$adjp<0.05,]
nsigBFCons <- sigBFCons[order(abs(sigBFCons$Beta),decreasing = T),]

upp <- 1/(1+exp(-(mx%*%Beta)))
nupp <- upp*B/(upp*B+1-upp)
par(bg="#D6D6CE")
plot(nupp~log(B),main="Overall model",xlab="log(B)",ylab="posterior",cex.lab=2,cex.axis=2,cex.main=2)
