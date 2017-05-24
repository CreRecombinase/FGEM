library(plyr)
library(reshape2)
library(topGO)
setwd("~/Dropbox/BayesianDA/FGEM/Data/")
  analysesdir <- "~/Dropbox/BayesianDA/FGEM/Data/Conservation/"
featmatf <- "~/Dropbox/BayesianDA/FGEM/Data/ConstraintMat.RDS"
BayesFactorFile <- "TADA_ASC_SSC_results_Dec23.csv"

featmat <- readRDS(featmatf)
resfiles <- dir(analysesdir,full.names=T)
genelist <- rownames(GOmat)
BF.df <- read.table(BayesFactorFile,header=T,sep=",")
rownames(BF.df) <- BF.df$Gene
BF.df <- BF.df[rownames(featmat),]
# GOMF.G <- melt(annFUN.org("MF",feasibleGenes=genelist,mapping="org.Hs.eg.db",ID="symbol"))
# colnames(GOMF.G) <- c("Gene","GOT")
# GOMF.G$Source <- "MF"
# GOBP.G <- melt(annFUN.org("BP",feasibleGenes=genelist,mapping="org.Hs.eg.db",ID="symbol"))
# colnames(GOBP.G) <- c("Gene","GOT")
# GOBP.G$Source <- "BP"
# GO.G <- rbind(GOMF.G,GOBP.G)

resdf <- ldply(resfiles,readRDS)
resdf <- rename(resdf,replace=c(X1="Feature"))
resdf <- resdf[!duplicated(resdf$Feature),]
rownames(resdf) <- resdf$Feature
resdf$pval <- pchisq(resdf$Chisq,df = 1,lower.tail = F)
resdf$adjp <- p.adjust(resdf$pval,method="fdr")
resdf <- resdf[order(resdf$pval),]
Beta <- c(resdf[1,"Intercept"],resdf[1,"Beta"])

B <- BF.df$BF
x <- featmat[,"GO:0042733"]
pvec <- 1/(1+exp(-(Beta[1]+Beta[2]*x)))
post <- (pvec*B)/(pvec*B+(1-pvec))
Bpdiff <- rank(B)-rank(pvec)
resdf <- resdf[order(resdf$pval),]
summary(abs(resdf$Beta))
which.max(abs(resdf$Beta))


rcq <- rchisq(length(resdf$pval),df = 1)
qqplot(qchisq(ppoints(length(resdf$pval)),df=1),resdf$Chisq,main="QQ plot",xlab="Theoretical X^2",ylab="Observed Chisq")
qqline(resdf$Chisq,distribution=function(p)qchisq(p,df=1))
summary(resdf$adjp)

GO
