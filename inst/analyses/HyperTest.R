library(plyr)
library(dplyr)
library(ggplot2)
library(GOstats)
GOmf <- "~/Dropbox/BayesianDA/FGEM/Data/GOmat.RDS"
BayesFactorFile <- "~/Dropbox/BayesianDA/FGEM/Data/TADA_ASC_SSC_results_Dec23.csv"
RDSres <- dir("~/Dropbox/BayesianDA/FGEM/Data/FGEM_LOG_RDS_B/",full.names = T)
rdsdf <- bind_rows(lapply(RDSres,function(x){
  tdf <- readRDS(x)
  tdf$got <- rownames(tdf)
  return(tdf)
}))
rdsdf <- mutate(rdsdf,chisq=-2*(Null-Likelihood),
                pch=pchisq(chisq,df=1,lower.tail=F),
                pa=p.adjust(pch,method = "fdr")) %>% arrange(pa)

BF.df <- read.table(BayesFactorFile,header=T,sep=",")
sigB <- BF.df$BF>5
GOmat <- readRDS(GOmf)
siggo <- GOmat[,head(rdsdf$got)]


nafe <- adply(GOmat,2,function(x,B)fisher.test(table(x,B),alternative = "greater")$p.value,B=sigB,.progress="text")
afe <- adply(siggo,2,function(x,B)fisher.test(table(x,B))$p.value,B=sigB,.progress="text")
afe <- rename(nafe,got=X1,fishp=V1)

rdsdf <- left_join(rdsdf,afe)

ggplot(rdsdf)+geom_point(aes(x=fishp,y=pch))

allgenes <- data.frame(Gene=BF.df$Gene,B=BF.df$BF,stringsAsFactors=F)
rownames(allgenes) <- allgenes$Gene
allgenes <- allgenes[rownames(GOmat),]

library("org.Hs.eg.db")
 frame = toTable(org.Hs.egGO)
 goframeData = data.frame(frame$go_id, frame$Evidence, frame$gene_id)
 goFrame=GOFrame(goframeData,organism="Homo sapiens")
 goAllFrame=GOAllFrame(goFrame)
 gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
GOparams <- new("GOHyperGParams",
                 geneIds=rownames(GOmat)[sigB],
                 universeGeneIds=rownames(GOmat),
                 annotation="org.Hs.eg.db",
                 ontology="BP",
                 pvalueCutoff=0.001,
                 conditional=TRUE,
                 testDirection="over")

hgOver <- hyperGTest(GOparams)

