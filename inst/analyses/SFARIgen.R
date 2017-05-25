#options(java.parameters = "-Xmx25g")
#library(bartMachine)
#library(randomForest)

outdir <- "~/Dropbox/BayesianDA/GSE/"
featdir <- "~/Desktop/SFARI/"

source("~/Dropbox/BayesianDA/FGEM/SFARIFunc.R")
cl <- makeCluster(detectCores())
registerDoParallel(cl)

# exp <- genexp(filelist)
# print("exp read!")
consdf <- gen_annotations(filelist,featdir)
print("annotations read!")

consdf <- mutate(consdf,TADAq=ifelse(is.na(TADAq),1,TADAq))
consdf <- mutate(consdf,TADAi=ifelse(TADAq<0.15,1,0))
consdf <- mutate(consdf,UASD=pmax(TADAi,ASDc))


nconsdf <-  select(consdf,-ASDc,-TADAq,-TADAi) 
nconsdf <- rename(nconsdf,y=UASD)
nconsdf <- combofeat(nconsdf,sep="5")

nconsdf <- nconsdf[complete.cases(nconsdf),]
predgene <- nconsdf$gene
yl<-filter(nconsdf,y>0)
rl <- filter(nconsdf,y==0)
bl <- list()
gl <- list()
featn <- "COMBUASD"
nx <- data.matrix(select(nconsdf,-y,-gene))
print(paste0("Analyzing ",featn))
for(i in 1:300){
  print(i)
  tdf <- bind_rows(yl,sample_n(rl,size = nrow(yl),replace=F))
  r <- data.matrix(select(tdf,-y,-gene))
  y <- tdf[["y"]]
  tmod <- cv.glmnet(r,y,family="binomial",alpha=0.95,maxit=750000,parallel = T,type.measure="class")
  tml <- tmod$lambda.min
  err <- min(tmod$cvm)
  tpred  <- predict.lognet(tmod$glmnet.fit,newx = nx,s=tml,type="response")
  bb <- tmod$glmnet.fit$beta[,which.min(abs((tmod$glmnet.fit$lambda-tmod$lambda.min)))]
  gl[[i]] <- data_frame(pred=tpred[,1],gene=predgene,iter=i,coef=featn)
  bl[[i]] <- data_frame(coef=names(bb),B=as.numeric(bb),iter=i,feat=featn,err=err)
}
print("Done!")
bdf <- bind_rows(bl)
gdf <- bind_rows(gl)
err_df <- group_by(bdf,iter) %>% summarise(err=err[1])
group_by(bdf,coef) %>% summarise(pB=sum(B!=0)/n(),medB=median(B),meanB=mean(B)) %>% arrange(desc(pB)) %>% data.frame

ggplot(bdf)+geom_histogram(aes(x=err))
#bdfs <- bdf %>% separate(coef,c("class","rest","summ"),sep="(_{3})|(-{3})",fill="right") %>% mutate(rest=ifelse(is.na(rest),class,rest),summ=ifelse(is.na(summ),rest,summ))
foutfile <- file.path(outdir,paste0(featn,"-3.RDS"))
soutfile <- file.path(outdir,paste0(featn,"-3-PRED.RDS"))
saveRDS(bdf,foutfile)
saveRDS(gdf,soutfile)

