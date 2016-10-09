


#Given a sparse matrix, write it to hdf5 format
# sparseAnnowrite <- function(h5f,)

#Code to read sparse matrix annotations from an HDF5 file ()
sparseAnnoread <- function(h5f,datapath,rowpath="GeneList/",colpath){
  require(rhdf5)
  require(Matrix)
  rowname <- h5read(h5f,rowpath)
  colname <- h5read(h5f,colpath)
  dat <-h5read(h5f,datapath)
  datmat <- sparseMatrix(i=dat[,1],j=dat[,2],x=1,
                         dims=c(length(rowname),length(colname)),
                         dimnames=list(rowname,colname))
  return(datmat)
}

#Return a matrix with two columns,z_score, and Bayes Factor (genes are rownames)
cons_to_bf <- function(consmat,allgenes){
  require(dplyr)
  consdf <- data_frame(Gene=rownames(consmat),mis_z=consmat[,"mis_z"])
  acons <- inner_join(consdf,allgenes)
  acm <- as.matrix(select(acons,-Gene))
  rownames(acm) <- acons$Gene
  return(acm)
}



GenPrior <- function(Beta,x){
  mx <- cbind(1,x)
  pvec <- 1/(1+exp(-(mx%*%Beta)))
  return(pvec)
}


Logit_FGEM <- function(u,x,B){
  mx <- cbind(1,x)
  Beta <- coefficients(glm(u~x,family=quasibinomial(link="logit")))
  pvec  <- 1/(1+exp(-(mx%*%Beta)))
  uvec <- (pvec*B)/((pvec*B)+(1-pvec))
  return(uvec)
}

B_Logit_FGEM <- function(Beta,x,B,isIntercept=F){
  if(!isIntercept){
    mx <- cbind(1,x)
  }else{
    mx <- cbind(x)
  }
  pvec  <- 1/(1+exp(-(mx%*%Beta)))
  uvec <- (pvec*B)/((pvec*B)+(1-pvec))
  if(!isIntercept){
    Beta <- coefficients(glm(uvec~x,family=quasibinomial(link="logit")))
  }else{
    Beta <- coefficients(glm(uvec~x+0,family=quasibinomial(link="logit")))
  }
  return(Beta)
}

B_Logit_log_lik <- function(Beta,x,B,isIntercept=F){
  if(!isIntercept){
    mx <- cbind(1,x)
  }else{
    mx <- cbind(x)
  }
  pvec  <- 1/(1+exp(-(mx%*%Beta)))
  # opvec  <- 1/(1+exp((mx%*%Beta)))
  uvec <- (pvec*B)/((pvec*B)+(1-pvec))
  return(-sum(uvec*log(pvec+B)+(1-uvec)*log(1-pvec)))
}



Logit_log_lik <- function(u,x,B){
  mx <- cbind(1,x)
  Beta <- coefficients(glm(u~x,family=quasibinomial(link="logit")))
  pvec  <- 1/(1+exp(-(mx%*%Beta)))
  # opvec  <- 1/(1+exp((mx%*%Beta)))
  uvec <- (pvec*B)/((pvec*B)+(1-pvec))
  return(-sum(uvec*log(pvec+B)+(1-uvec)*log(1-pvec)))
}


Beta_log_lik <- function(Beta,x,B,isIntercept=F){
  if(!isIntercept){
    mx <- cbind(1,x)
  }else{
    mx <- cbind(x)
  }
  pvec  <- 1/(1+exp(-(mx%*%Beta)))
  uvec <- (pvec*B)/((pvec*B)+(1-pvec))
  return(-sum(uvec*log(pvec+B)+(1-uvec)*log(1-pvec)))
}

gen_p <- function(u,x,B){
  mx <- cbind(1,x)
  pvec  <- 1/(1+exp(-(mx%*%Beta)))
  Beta <- coefficients(glm(u~x,family=quasibinomial(link="logit")))
  return(Beta)
}

Bayesfile <- function(BayesFactorFile){
  BF.df <- read.table(BayesFactorFile,header=T,sep=",")
  allgenes <- data.frame(Gene=BF.df$Gene,B=BF.df$BF,stringsAsFactors=F)
  allgenes$Gene <- as.character(allgenes$Gene)
  rownames(allgenes) <- allgenes$Gene
  return(allgenes)
}

gen_final_mat <- function(X1,GOmatfile,Consfile,BayesFactorFile,expfeatfile=NULL){
  gox <- X1[grepl("GO",X1)]
  consx <- X1[!grepl("GO",X1)]
  GOmat <- readRDS(GOmatfile)
  consdat <- readRDS(Consfile)
  consg <- rownames(consdat)
  consdat <- as_data_frame(consdat) %>% select(one_of(X1)) %>% mutate(Gene=consg)
allgenes <- Bayesfile(BayesFactorFile)
  consdat <- inner_join(consdat,allgenes)
  GOdf <- as_data_frame(data.frame(as.matrix(GOmat[,gox,drop=F]),check.names = F,stringsAsFactors = F)) %>% mutate(Gene=rownames(GOmat))
  if(!is.null(expfeatfile)){
    expdat <- readRDS(expfeatfile)
    expg <- rownames(expdat)
    expdf <- as_data_frame(expdat) %>% select(one_of(X1)) %>% mutate(Gene=expg)
    consdf <- inner_join(expdf,consdat)
  }
  findf <- inner_join(consdf,GOdf)
  fmat <- data.matrix(select(findf,-Gene,-B))
  rownames(fmat) <- findf$Gene
  return(list(fmat,B=findf$B))
}

anno2df <- function(annodf,feat.name="ExAC"){
  require(dplyr)
  require(tidyr)
  stopifnot(any(colnames(annodf)=="Gene"))
  retdf <- gather(annodf,feature,value,-Gene) %>% mutate(class=feat.name)
}


exp2df<- function(expmat,feat.name="Expression"){
  require(dplyr)
  require(tidyr)
  Genes <-rownames(expmat)
  expdf <- as_data_frame(expmat) %>% mutate(gene=Genes) %>% gather(feature,value,-gene) %>% mutate(class=feat.name)
  return(expdf)
}

glist2df <- function(Genes,feat.name){
  require(dplyr)
  retdf <- data_frame(gene=Genes) %>% mutate(feature=feat.name,value=1,class=feat.name)
  return(retdf)
}

GO2df <-function(Genes,terms="BP"){
  require(topGO)
  require(org.Hs.eg.db)
  require(dplyr)
  require(reshape2)
  mdf <- melt(annFUN.org(terms,feasibleGenes=Genes,mapping="org.Hs.eg.db",ID="symbol"))
  mdf <- as_data_frame(mdf)
  mdf <- dplyr::select(mdf,gene=value,feature=L1) %>% mutate(value=1,class=paste0("GO_",terms))
  return(mdf)
}


cfeat_df <- function(annodf,datadf,impute=F){
  require(dplyr)
  require(SQUAREM)
  stopifnot(length(unique(annodf$feature))==1)
   isbin <- all(annodf$value==1)
   if(isbin){
     full_feat <- left_join(datadf,annodf)
     full_feat <- mutate(full_feat,feature=feature[!is.na(feature)][1],
                         class=class[!is.na(class)][1],
                         value=ifelse(is.na(value),0,1))
   }else{
     if(impute){
       full_feat <- left_join(datadf,annodf)
       full_feat <- mutate(full_feat,value=ifelse(is.na(value),mean(value,na.rm=T),value))
     }else{
       full_feat <- inner_join(datadf,annodf)
     }
   }
   return(full_feat)
}

pmean <-function(Beta,feat_mat,isIntercept=F){
  if(!isIntercept){
    mx <- cbind(1,feat_mat)
  }else{
    mx <- cbind(feat_mat)
  }
  pvec  <- 1/(1+exp(-(mx%*%Beta)))
  return(mean(pvec))
}


pem_df <-function(full_feat){
  cat(length(unique(full_feat$feature)))
  stopifnot(length(unique(full_feat$feature))!=1)
  feat_mat <-select(full_feat,Gene,feature,value,BF) %>% spread(feature,value) %>% select(-Gene) %>% filter(complete.cases(.))

  BF <-feat_mat$BF
  tmu <-BF/(BF+exp(-BF))
  feat_mat <- data.matrix(select(feat_mat,-BF))
  fBeta <- coefficients(glm(tmu~feat_mat,family=quasibinomial(link="logit")))
  opf <- squarem(par=fBeta,fixptfn=B_Logit_FGEM,
                 objfn=B_Logit_log_lik,x=feat_mat,B=BF)
  pBeta <- opf$par
#  opf$par <- NULL
  cn <- c(colnames(feat_mat))
  nret <- data.frame(opf) %>% mutate(feature=c("Intercept",cn),value.objfn=-value.objfn) %>%
    rename(LogLik=value.objfn,Beta=par)
  nNullLogLik <- -Beta_log_lik(c(pBeta[1],rep(0,ncol(feat_mat))),feat_mat,BF)
  prior_mean <- pmean(pBeta,feat_mat,isIntercept = T)
  nret <- mutate(nret,NullLogLik=nNullLogLik,Chisq=2*(NullLogLik-LogLik),
                 pval=pchisq(Chisq,df=1,lower.tail = F),prior_mean=prior_mean)
  return(nret)
}


Null_Intercept <-function(full_feat){
  stopifnot(length(unique(full_feat$feature))==1)
  feat_mat <-select(full_feat,Gene,feature,value,BF) %>% mutate(value=1,feature="Intercept") %>% spread(feature,value) %>% select(-Gene)
  BF <-feat_mat$BF
  tmu <-BF/(BF+exp(-BF))
  feat_mat <- data.matrix(select(feat_mat,-BF))
  fBeta <- coefficients(glm(tmu~feat_mat+0,family=quasibinomial(link="logit")))
  opf <- squarem(par=fBeta,fixptfn=B_Logit_FGEM,
                 objfn=B_Logit_log_lik,x=feat_mat,B=BF,isIntercept=T)
  pBeta <- opf$par
  opf$par <- NULL
  cn <- c(colnames(feat_mat))
  nret <- data.frame(opf) %>% mutate(feature=cn,value.objfn=-value.objfn,Intercept=pBeta[1]) %>%
    rename(LogLik=value.objfn)
  nNullLogLik <- -Beta_log_lik(c(pBeta[1]),feat_mat,BF,isIntercept = T)
  prior_mean <- pmean(pBeta,feat_mat,isIntercept = T)
  nret <- mutate(nret,NullLogLik=nNullLogLik,prior_mean=prior_mean)
}


sem_df <-function(full_feat,NullIntercept=NULL){
  require(dplyr)
  require(SQUAREM)
  require(tidyr)
  cat(length(unique(full_feat$feature)))
  stopifnot(length(unique(full_feat$feature))==1)
  feat_mat <-select(full_feat,Gene,feature,value,BF) %>% spread(feature,value) %>% select(-Gene)
  BF <-feat_mat$BF
  tmu <-BF/(BF+exp(-BF))
  feat_mat <- data.matrix(select(feat_mat,-BF))
  fBeta <- coefficients(glm(tmu~feat_mat,family=quasibinomial(link="logit")))
  opf <- squarem(par=fBeta,fixptfn=B_Logit_FGEM,
                 objfn=B_Logit_log_lik,x=feat_mat,B=BF)
  pBeta <- opf$par
  opf$par <- NULL
  cn <- c(colnames(feat_mat))
  nret <- data.frame(opf) %>% mutate(feature=cn,value.objfn=-value.objfn,Beta=pBeta[2],Intercept=pBeta[1]) %>%
    rename(LogLik=value.objfn)
  nNullLogLik <- -Beta_log_lik(c(pBeta[1],rep(0,ncol(feat_mat))),feat_mat,BF)
  prior_mean <- pmean(pBeta,feat_mat)
  nret <- mutate(nret,NullLogLik=nNullLogLik,Chisq=2*(NullLogLik-LogLik),
                 pval=pchisq(Chisq,df=1,lower.tail = F),prior_mean=prior_mean)
  return(nret)
}

fisher_comp <- function(full_feat,prior=0.02){

    feat_res <- group_by(full_feat,feature) %>% mutate(isbin=(length(unique(value))==2),
                                           feat_bin=factor(ifelse(isbin,factor(value),value>quantile(value,prior))),
                                           BF_bin=factor(BF>quantile(BF,prior))) %>%
do(par=chisq.test(.$feat_bin,.$BF_bin,simulate.p.value = T)) %>% summarise(feature=feature[1],prior=prior,statistic=par$statistic,pval=par$p.value)
    return(feat_res)
}



fullsem <- function(fmat,B){
  pf <- squarem(par=B/(B+exp(-B)),fixptfn=Logit_FGEM,
                objfn=Logit_log_lik,x=fmat,B=B)
  return(pf$par)
}

compute_full_FGEM <-function(GOmatfile,Consfile,EXPfile,Bfile,GOresfiles){

  resf <- c(Consfile,
            GOresfiles,
            "~/Dropbox/BayesianDA/FGEM/Data/Expression_Res.RDS")

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
  fmldf <- mutate(fmldf,rank_change = rank(BF)-rank(posterior),prior=c(prior))
  return(fmldf)
}




