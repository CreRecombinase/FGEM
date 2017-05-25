library(dplyr)
library(glmnet)
library(ggplot2)
library(tidyr)
library(doParallel)
library(methods)


featdir <- "~/Desktop/SFARI/"
filelist=list(
  gsf = "master_score.csv",
  gscf = "gene-score.csv",
  sczf ="Purcell_NG_2014_composite_set_genelist.txt",
  fmrpf = "Darnell_Cell_2011_FMRP_targets_genelist.txt",
  hapf = "Huang_Plosgen_haploinsufficiency_score_by_gene.txt",
  IDf = "Pinto_AJHG_ID_genelist.txt",
  PSdf = "Bayes_postsynaptic_genelist.txt",
  gomf = "GOmat.RDS",
  dbf = "CCTD_VarDB.db",
  NE_Exacf = "NP_Exac_Exon_Cons.txt",
  expdir = featdir,
  NG_Exacf = "NP_Exac_Gene_Cons.txt"
)




tauf <- function(x){
  xl <- x/max(x)
  return(sum(1-xl)/(length(x)-1))
}

summ_tau_g <- function(df,gene_rows,field){
  df<- group_by(df,gene) %>% summarise(tau_mean=tauf(mean)) 
  tdf <-   mutate(df,gene=as.integer(gene)) %>% 
    inner_join(gene_rows,by=c("gene"="row_num")) %>% 
    group_by(gene_symbol) %>% 
    summarise(min_t=min(tau_mean,na.rm=T),max_t=max(tau_mean,na.rm=T),mean_t=mean(tau_mean,na.rm=T)) %>% rename(gene=gene_symbol)
  ct <- colnames(tdf) 
  colnames(tdf)<- c(ct[1],paste0(field,"___",ct[-1]))
  return(tdf)
}

broad_exp <- function(df,field,sep="---"){
  ncn <- colnames(df)[!colnames(df) %in% c("gene",field)]
  tdf <- gather_(df,"feat","val",ncn)
  nname <- paste0(field,"_feat")
  tdf <- unite_(tdf,nname,c(field,"feat"),sep=sep)
  tdfm <- spread_(tdf,nname,"val")
  ct <- colnames(tdfm)
  colnames(tdfm) <- c(ct[1],paste0(field,"___",ct[-1]))
  return(tdfm)
}

broad_gene <- function(df,gene_rows,field,ge){
  df <- mutate(df,gene=as.integer(gene)) %>% inner_join(gene_rows,by=c("gene"="row_num")) %>% 
    select(-entrez_id,-ensembl_gene_id,-gene_id,-gene) %>% rename(gene=gene_symbol) %>% group_by_(field,"gene") %>% 
    summarise(min=min(min,na.rm=T),max=max(max,na.rm=T),mean=mean(mean,na.rm=T),median=mean(median,na.rm=T),mad=mean(mad,na.rm=T),var=mean(var,na.rm=T)) %>% ungroup()
  df <- group_by_(df,field) %>% mutate(perc_mean=percent_rank(mean)) %>% ungroup()
  ct <- colnames(df)
#  colnames(df) <- c(ct[1],ct[2],paste0(ge,"_",ct[-c(1,2)]))
  dfm <- broad_exp(df,field)
  return(dfm)
}

ID_feat <- function(exp,asddf){
  
  min_feat <- select(exp,gene,ends_with("mad"),ends_with("min_t"),ends_with("gene_perc_mean"))
  mindf <- left_join(asddf,exp)
  mindf <- mindf[,apply(mindf,2,function(x)sum(is.na(x))/length(x))!=1]
  mindf <- mindf[complete.cases(mindf),]
  minX <- data.matrix(select(mindf,-gene,-ASDc))
  minY <- mindf$ASDc
  dfl <- list()
  cvl <- list()
  for(i in 1:10){
    print(i)
    cvl[[i]] <- cv.glmnet(minX,minY,family="binomial",alpha=.95,maxit=1000000)
    mbl <- cvl[[i]]$lambda.min
    minbetas <- cvl[[i]]$glmnet.fit$beta[,which.min(abs(cvl[[i]]$glmnet.fit$lambda-mbl))]
    dfl[[i]] <- data_frame(betas=minbetas,names=names(minbetas)) 
    dfl[[i]] <- dfl[[i]]%>% separate(names,c("class","rest","feat"),sep="(_{3})|(-{3})",fill="right")
  }
  for( i in 1:10){
    dfl[[i]] <- mutate(dfl[[i]],iter=i)
  }
  fdfl <- bind_rows(dfl)
  fdfl <- mutate(fdfl,feat=ifelse(is.na(feat),rest,feat))
  fdfl %>% group_by(class) %>% summarise(mean(betas))
  fdfl %>% filter(class=="structure_acronym") %>% group_by(rest) %>% summarise(mb=mean(betas)) %>% arrange(desc(abs(mb))) %>% data.frame
  fdc <- fdfl %>% group_by(rest,feat) %>% summarise(class=class[1],mb=abs(mean(betas)),fnz=sum(betas!=0)/n()) %>% ungroup() %>%  arrange(desc(fnz)) %>% filter(fnz>.1)
  return(fdc)
}

filt_feat <- function(exp,fdc){
  fdc <- mutate(fdc,col=ifelse(grepl("*_t$",feat),paste0(class,"?",feat),paste0(class,"?",rest,":",feat)))
  nexp <- select(exp,gene,one_of(fdc$col))
  return(nexp)
}




genexp <- function(filelist,structnum=15){
  
  
  expdir <- filelist[["expdir"]]
  gene_stagef <- file.path(expdir,"gene_stage_summary.RDS")
  gene_stagetf <- file.path(expdir,"gene_stage_tau.RDS")
  gene_structf <- file.path(expdir,"gene_struct_summary.RDS")
  gene_structtf <- file.path(expdir,"gene_struct_tau.RDS")
  gene_rowsf <- file.path(expdir,"gene_row_metadata.RDS")
  gene_cf <- file.path(expdir,"both_column_metadata.RDS")
  
  gene_cols <- readRDS(gene_cf)
  goodstruct <- group_by(gene_cols,structure_acronym) %>% summarise(nsamp=n()) %>% filter(nsamp>=structnum)
  
  gene_rows <- readRDS(gene_rowsf)
  gene_stage <- readRDS(gene_stagef)
  gene_struct <- readRDS(gene_structf)
  gene_staget <- readRDS(gene_stagetf)
  gene_structt <- readRDS(gene_structtf)
  
  gene_struct <- semi_join(gene_struct,goodstruct)
  
  ngene_staget <- summ_tau_g(gene_stage,gene_rows,"stage")
  ngene_structt <- summ_tau_g(gene_struct,gene_rows,"structure_acronym")
  tau <- inner_join(ngene_staget,ngene_structt)
  stage <- broad_gene(gene_stage,gene_rows,"stage","gene")
  struct <- broad_gene(gene_struct,gene_rows,"structure_acronym","gene")
  exp <- inner_join(stage,struct)
  exp <- left_join(exp,tau)
  return(exp)
}




gen_annotations <- function(filelist,featdir){
    gsf <- file.path(featdir,filelist[["gsf"]])
    gscf <- file.path(featdir,filelist[["gscf"]])
    sczf <- file.path(featdir,filelist[["sczf"]])
    hapf <- file.path(featdir,filelist[["hapf"]])
    IDf <- file.path(featdir,filelist[["IDf"]])
    PSdf <- file.path(featdir,filelist[["PSdf"]])
    fmrpf <- file.path(featdir,filelist[["fmrpf"]])
    gomf <- file.path(featdir,filelist[["gomf"]])
    dbf <- file.path(featdir,filelist[["dbf"]])
    NE_Exacf <- file.path(featdir,filelist[["NE_Exacf"]])
    NG_Exacf <- file.path(featdir,filelist[["NG_Exacf"]])
    
    
    
    NEE <- read.table(NE_Exacf,header=T,sep="\t",stringsAsFactors = F)
    NEE <- group_by(NEE,gene) %>% summarise(NP_exon_max_z=max(z_sign_mis))
    gez <- NEE


   
    mdb <- src_sqlite(dbf,create=F)
    
    TADA <- data.frame(collect(tbl(mdb,"TADABF"))) %>% select(gene,TADAq=qvalue)
    
    bg  <- data.frame(collect(tbl(mdb,"BrainGO")))
    fmrpd <- scan(fmrpf,what=character())
    hapd <- read.table(hapf,stringsAsFactors = F)
    hapd <- rename(hapd,gene=V1,hs=V2)
    idd <- scan(IDf,what=character())
    psdd <- scan(PSdf,what=character())
    sczd <- scan(sczf,what=character())
    
    
    GOmat <- readRDS(gomf)
    bgm <- GOmat[,unique(bg$got)]
    bgm <- as.matrix(bgm)
    gsd <- read.table(gsf,header = F,sep=",",fill=T,stringsAsFactors = F,comment.char = '',quote = "\"")
    gsc <- trimws(scan(gscf,what=character(),sep=",",nlines = 1))
    gsc <- gsub(pattern = " ",replacement = ".",gsc)
    colnames(gsd) <- gsc
    gsd <- gsd[!duplicated(gsd$Gene.Symbol),]
    rownames(gsd) <- gsd$Gene.Symbol
    
    sczdf <- data_frame(gene=sczd,SCZ=1)
    fmrpdf <- data_frame(gene=fmrpd,FMRP=1)
    psdf <- data_frame(gene=psdd,PSD=1)
    iddf <- data_frame(gene=idd,IDG=1)
    bgmdf <- as_data_frame(data.frame(bgm)) %>% mutate(gene=rownames(bgm))
    gsddf <- select(gsd,gene=Gene.Symbol,Score) %>% mutate(is1=ifelse(Score=="1",1,0),
                                                         is1S=ifelse(Score=="1S",1,0),
                                                         is2=ifelse(Score=="2",1,0),
                                                         is2S=ifelse(Score=="2S",1,0),
                                                         is3=ifelse(Score=="3",1,0),
                                                         is3S=ifelse(Score=="3S",1,0),
                                                         is4=ifelse(Score=="4",1,0),
                                                         is4S=ifelse(Score=="4s",1,0),
                                                         is5=ifelse(Score=="5",1,0),
                                                         is6=ifelse(Score=="6",1,0),
                                                         isS=as.integer((is1S|is2S|is3S)))
    gsddf <- mutate(gsddf,ASDc=as.integer(is1|is1S|is2|is2S))%>% select(-starts_with("is"),-Score)
    
    
    any_bg <- apply(select(bgmdf,-gene),1,function(x)max(x))
    any_bgdf <- select(bgmdf,gene) %>% mutate(any_Brain_GO=any_bg)
    
    consdf <- gsddf%>% full_join(iddf) %>% full_join(hapd) %>% full_join(sczdf) %>% full_join(fmrpdf) %>%
      full_join(psdf) %>%    full_join(gez) %>%
      left_join(any_bgdf) %>% left_join(TADA) %>%
      replace_na(replace=list(IDG=0,PSD=0,any_Brain_GO=0,ASDc=0,hs=0,SCZ=0,FMRP=0)) %>%  filter(!gene %in% c("AMBIGUOUS","NOT_FOUND"))
    return(consdf)
}

multinom.glm <- function(consdf,exp,outfile){
  
  ndat <- select(consdf,gene,ASDc,IDG,PSD,any_Brain_GO,TADAi,NP_exon_max_z) %>% inner_join(exp)
  ndat <- ndat[complete.cases(ndat),]
  nx <- select(ndat,-ASDc,-IDG,-PSD,-any_Brain_GO,-TADAi,-gene) %>% data.matrix()
  ny <- select(ndat,ASDc,IDG,PSD,any_Brain_GO,TADAi) %>% data.matrix()
  nym <- t(t(as.integer(rowSums(ny)==0)))
  colnames(nym) <- c("None")
  nym <- cbind(ny,nym)
  nmod <- cv.glmnet(nx,nym,alpha=0.95,family="multinomial",parallel=T,maxit=750000)
  tmod <- cv.glmnet(nx,ny,alpha=0.95,family="multinomial",parallel=T)
  fmod <- cv.glmnet(nx,ny,alpha=0.95,family="multinomial",type.multinomial="grouped",parallel=T,maxit=750000)
  nmlam <- nmod$lambda.min
  tmlam <- tmod$lambda.min
  bflam <- fmod$lambda.min
  npredt <- predict.lognet(nmod$glmnet.fit,newx=nx,s=nmlam,type="response")
  npredc <- predict.lognet(nmod$glmnet.fit,newx=nx,s=nmlam,type="class")
  str(nmod$glmnet.fit$beta)
  abeta <- lapply(nmod$glmnet.fit$beta,"[",1:203,which(nmod$lambda==nmlam))
  nbeta <- do.call("cbind",abeta)
  colnames(nbeta) <- c(colnames(nbeta)[-ncol(nbeta)],"None")
  tnbeta <- 
  
  nbd <- as_data_frame(as.data.frame(nbeta))
  nbd <- mutate(nbd,feat=rownames(nbeta))
  filter(nbd,abs(IDG)>0)
  
  predt <- predict.lognet(tmod$glmnet.fit,newx=nx,s=tmlam,type="response")
  predtc <- predict.lognet(tmod$glmnet.fit,newx=nx,s=tmlam,type="class")
  predf <- predict.lognet(fmod$glmnet.fit,newx=nx,s=bflam,type="response")
  predfc <- predict.lognet(fmod$glmnet.fit,newx=nx,s=bflam,type="class")
  predtfg <- predict.glmnet(fmod$glmnet.fit,newx=nx,s=bflam,type="response")
  predtcg<-  predict.glmnet(tmod$glmnet.fit,newx=nx,s=tmlam,type="response")
  predtcg <- predict(tmod,nx,s="lambda.min",type="response")
}

matthewc <- function(x,medp,y){
  g <- as.integer(medp>x)
  tp <- sum((y==1)&(g==1))
  tn <- sum((y==0)&(g==0))
  fp <- sum((y==0)&(g==1))
  fn <- sum((y==1)&(g==0))
  mcc <- ((tp*tn)-(fp*fn))/exp(0.5*log(tp+fp)+log(tp+fn)+log(tn+fp)+log(tn+fn))
  return(-mcc*10^6)
}

combofeat <- function(nconsdf,sep="5"){
  ncombm <- t(apply(select(nconsdf,IDG,SCZ,hs,FMRP,PSD,NP_exon_max_z,any_Brain_GO),1,function(x){
    combn(x,2,function(y){y[1]*y[2]})
  }))
  colnames(ncombm) <- combn(c("IDG","SCZ","hs","FMRP","PSD","NP_exon_max_z","any_Brain_GO"),2,function(y)paste(y[1],y[2],sep=sep))
  nconsdf <- bind_cols(nconsdf,as_data_frame(data.frame(ncombm)))
  nconsdf <- nconsdf[complete.cases(nconsdf),]
  return(nconsdf)
}


percf <- function(x,pred_df,ydf){
  cat(x)
  cat("\n")
  tdf <- group_by(pred_df,gene) %>% summarise(medp=quantile(pred,x[1])) %>% inner_join(ydf)
  mc <- matthewc(x[2],tdf[["medp"]],tdf[["y"]])
  cat(mc)
  cat("\n")
  return(mc)
}



