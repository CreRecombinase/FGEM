library(dplyr)
library(ggplot2)
source("~/Dropbox/BayesianDA/FGEM/SFARIFunc.R")
expf <- "~/Dropbox/BayesianDA/GSE/expd.RDS"
consf <- "~/Dropbox/BayesianDA/GSE/consd.RDS"
outdir <- "~/Dropbox/BayesianDA/GSE/"
featdir <- "~/Desktop/SFARI/"


#exp <- genexp(fl)
#saveRDS(exp,expf)
consdf <- gen_annotations(filelist,featdir)
consdf <- mutate(consdf,TADAq=ifelse(is.na(TADAq),1,TADAq))
consdf <- mutate(consdf,TADAi=ifelse(TADAq<0.15,1,0))
consdf <- mutate(consdf,UASD=pmax(TADAi,ASDc))
#consdf <- mutate(consdf,TADAi=ifelse(TADAq<0.1,1,0))
#saveRDS(consdf,consf)
#exp <- readRDS(expf)
#consdf <- readRDS(consf)

sfarif <- dir("~/Dropbox/BayesianDA/GSE/",pattern="*3.RDS",full.names=T)

sfarip <- dir("~/Dropbox/BayesianDA/GSE/",pattern="*3-PRED.RDS",full.names=T)

result_df <-readRDS("~/Dropbox/BayesianDA/GSE/UASD-3.RDS")
pred_df <- readRDS("~/Dropbox/BayesianDA/GSE/UASD-3-PRED.RDS")
err_df <- group_by(result_df,iter) %>% summarise(err=err[1])
ggplot(err_df)+geom_histogram(aes(x=err),binwidth=.005)+ggtitle(label="misclassification error across iterations")


betasum <- result_df %>% group_by(feat,coef) %>% summarise(medB=median(B),pb=sum(B!=0)/n()) %>% ungroup() %>% arrange(desc(pb))
mpred <- group_by(pred_df,gene) %>% summarise(varp=var(pred),meanp=mean(pred),medp=median(pred),minp=min(pred),maxp=max(pred)) %>% ungroup() %>% arrange(desc(minp))
predf <- inner_join(consdf,mpred) %>% arrange(desc(minp))
predf <- arrange(predf,desc(medp))
predf <- mutate(predf,meanpredz=(meanp-mean(meanp))/sqrt(varp),
                medpredz=(medp-mean(medp))/sqrt(varp),
                minpredz=(minp-mean(minp))/sqrt(varp),
                maxpredz=(maxp-mean(maxp))/sqrt(varp)) %>%
  mutate(meanpredz=rank(meanpredz),medpredz=rank(medpredz),minpredz=rank(minpredz),maxpredz=rank(maxpredz))

predf <- arrange(predf,desc(maxpredz))
filter(predf,UASD==1) %>%ggplot()+geom_point(aes(x=meanpredz,y=medpredz,colour=factor(UASD)))


medp <- predf[["minp"]]
y=predf[["UASD"]]
lowerb <- as.numeric(filter(predf,UASD==1) %>% summarise(min(minp)))+0.01
upperb <- max(medp)-.04
optimize(matthewc,lower=lowerb,upper=upperb,medp=medp,y=y)
head(filter(predf,minp<0.5305693))                 
