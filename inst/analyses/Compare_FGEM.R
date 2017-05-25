#Compare methods
library(dplyr)
library(SQUAREM)
library(FGEM)
GOmatfile <- "~/Dropbox/BayesianDA/FGEM_Data/GOmat.RDS"
GOmat <- readRDS(GOmatfile)
Consfile <- "~/Dropbox/BayesianDA/FGEM_Data/ConstraintMat.RDS"

Feat_file <- "~/Dropbox/BayesianDA/FGEM/SFARIgen_Exp_Features.RDS"
Bfile <- "~/Dropbox/BayesianDA/FGEM_Data/TADA_ASC_SSC_results_Dec23.csv"
BFile <- "TADA_ASC_SSC_results_Dec23.csv"
resf <- dir("~/Dropbox/BayesianDA/FGEM_Data/FGEM_GO_RDS/",full.names = T)
EXPfile <- "~/Dropbox/BayesianDA/FGEM_Data/BrainSpanVar.RDS"
resf <- c("~/Dropbox/BayesianDA/FGEM_Data/FGEM_Gene_Constraint.RDS",
          resf,
          "~/Dropbox/BayesianDA/FGEM_Data/Expression_Res.RDS")
# outfile <- "GOmat.RDS"
# trf <- readRDS(resf[1])
BF.G <- read.table(Bfile,header=T,sep=",")
mres <- bind_rows(lapply(resf,readRDS))
mres <- mutate(mres,feature=ifelse(is.na(feature),X1,feature)) %>% select(-X1,-method)

mres <- mutate(mres,p=pchisq(Chisq,df=1,lower.tail = F),isGO=grepl("GO",feature)) %>% group_by(isGO) %>% mutate(p.adj=p.adjust(p,method="fdr")) %>% ungroup() %>% arrange(p.adj)

X1 <- filter(mres,p.adj<0.05)[["feature"]]
fml <- gen_final_mat(X1,GOmatfile,Consfile,Bfile,EXPfile)
fmlm <- fml[[1]]
fmlm <- scale(fmlm,center=T,scale=F)

par <- fullsem(fmlm,fml[[2]])
tbeta <- gen_p(par,fml[[1]],fml[[2]])
B <- fml$B
tgp <- GenPrior(tbeta,fml[[1]])
plot(tgp,log(B))
tpost <- (tgp*B)/(tgp*B+(1-tgp))
# plot(log(B),tpost)
fmldf <- as_data_frame(data.frame(fml[[1]],check.names=F))
fmldf <- mutate(fmldf,Gene=rownames(fml[[1]]),BF=B,prior=tgp,posterior=tpost)
ggplot(fmldf) + geom_point(aes(x=prior,y=posterior))
ggplot(fmldf) + geom_point(aes(x=log(BF),y=posterior))
fmldf <- mutate(fmldf,rank_change = rank(BF)-rank(posterior),prior=c(prior))
# ggplot(fmldf)+geom_histogram(aes(x=posterior),binwidth=0.01)
# write.table(ofmldf,"~/Dropbox/BayesianDA/FGEM/ASD_Priors_scaled_061.txt",sep="\t",col.names=T,row.names=F,quote=F)
# oofmldf <- read.table("~/Dropbox/BayesianDA/FGEM/ASD_Priors_061.txt",sep="\t",header=T) %>% select(Gene,oldprior=prior)
# nfmldf <- inner_join(oofmldf,fmldf) %>% mutate(prior=c(prior))
# nfmldf <- as_data_frame(nfmldf)
# arrange(nfmldf,desc(prior)) %>% mutate(prior_rank=rank(prior),B_rank=rank(B),prior_diff=rank(oldprior)-prior_rank) %>% arrange(desc(prior_diff)) %>% select(-contains("."),-contains("GO"),-striatum) %>% head
# select(nfmldf,prior) %>% head
# plot(nfmldf$oldprior,nfmldf$prior)
# filter(fmldf,post<0.1,B>5) %>% arrange(desc(abs(rchange))) %>% head(10)
# sigfml <- filter(fmldf,post>0.9) %>% arrange(B)
# ggplot(fmldf)+geom_point(aes(x=log(B),y=post,colour=mis_z))+theme(axis.text=element_text(size=12,face="bold"),
#                                                                   axis.title=element_text(size=12,face="bold"))
