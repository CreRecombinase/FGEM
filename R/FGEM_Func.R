


FGEM_Logit <- function(Beta,x,B){
  pvec  <- 1/(1+exp(-(x%*%Beta)))
  uvec <- (pvec*B)/((pvec*B)+(1-pvec))
  Beta <- coefficients(glm.fit(x = x,y = uvec,family=quasibinomial(link="logit"),control=list(maxit=50)))
  return(Beta)
}

FGEM_Logit_fit <- function(Beta,x,B){
  pvec  <- 1/(1+exp(-(x%*%Beta)))
  uvec <- (pvec*B)/((pvec*B)+(1-pvec))

  return(Beta)
}



FGEM_Logit_log_lik <- function(Beta,x,B){
  pvec  <- 1/(1+exp(-(x%*%Beta)))
  # opvec  <- 1/(1+exp((mx%*%Beta)))
  uvec <- (pvec*B)/((pvec*B)+(1-pvec))
  return(-sum(uvec*(log(pvec)+log(B))+(1-uvec)*log(1-pvec)))
}



gen_p <- function(Beta,x){
  pvec  <- 1/(1+exp(-(x%*%Beta)))
  return(pvec)
}
gen_u <- function(Beta,x,B){
  p <- gen_p(Beta,x)
  return((p*B)/((p*B)+(1-p)))
}

anno2df <- function(annodf,feat.name="ExAC"){
  require(dplyr)
  require(tidyr)
  stopifnot(any(colnames(annodf)=="Gene"))
  retdf <- gather(annodf,feature,value,-Gene)  %>% mutate(class=feat.name)
}

exp2df<- function(expmat,feat.name="Expression"){
  require(dplyr)
  require(tidyr)
  Genes <-rownames(expmat)
  expdf <- as_data_frame(expmat) %>% mutate(Gene=Genes) %>% gather(feature,value,-Gene) %>% mutate(class=feat.name)
  return(expdf)
}

glist2df <- function(Genes,feat.name){
  require(dplyr)
  retdf <- data_frame(Gene=Genes) %>% mutate(feature=feat.name,value=1,class=feat.name)
  return(retdf)
}




cfeat_mat <- function(annodf,datadf){
  require(tidyr)
  require(dplyr)
  tannodf <- select(annodf,-one_of("feat_ind")) %>% spread(key = feature,value = value)
  return(inner_join(tannodf,datadf))
}

cfeat_df <- function(annotation_df,datadf){
  require(dplyr)
  require(tidyr)
  # stopifnot(length(unique(annodf$feature))==1)
  annotation_df <- group_by(annotation_df,feature) %>% mutate(isBin=(n_distinct(value)==1),ngenes=n()) %>%ungroup()
  mingenes <-filter(annotation_df,isBin==FALSE) %>% group_by(Gene) %>% summarise(nfeat=n_distinct(feature)) %>% ungroup() %>% filter(nfeat==max(nfeat))
  if(all(annotation_df$isBin)){
    mingenes <-distinct(annotation_df,Gene)
  }
  annotation_df <- annotation_df %>% inner_join(mingenes,by="Gene") %>% select(-one_of(c("feat_ind","nfeat","isBin","ngenes"))) %>% spread(feature,value,fill = 0)
  datadf <- select(datadf,Gene,BF)
  full_feat <- inner_join(datadf,annotation_df,by="Gene")
  return(full_feat)
}

gen_prior_df <- function(annodf,datadf,scale=F){
  tfeatdf <- distinct(annodf,feature,Beta,Intercept)
  tannodf <- select(annodf,-Beta,-Intercept) %>% group_by(feature) %>% do(cfeat_df(.,datadf))  %>%  ungroup() %>%
    inner_join(tfeatdf)
  am <-anno2mat(tannodf)
  bv <- result2beta(annodf)
  mprior <- 1/(1+exp(-am%*%bv[colnames(am)]))

  retdf <- group_by(tannodf,Gene) %>% summarise(prior=1/(1+exp(-(mean(Intercept)+sum(value*Beta))))) %>% filter(Gene %in% rownames(am))
  return(retdf)
}


pmean <-function(Beta,feat_mat){
  pvec  <- 1/(1+exp(-(feat_mat%*%Beta)))
  return(mean(pvec))
}



EM_mat <-function(fBeta,feat_matrix,BF){
  require(SQUAREM)
  require(dplyr)
  require(tidyr)
  opf <- squarem(par=fBeta,fixptfn=FGEM_Logit,x=feat_matrix,B=BF)
  cn <- colnames(feat_matrix)
  opf$LogLik <- -FGEM_Logit_log_lik(opf$par,feat_matrix,BF)
  opf$feature_name <- names(opf$par)
  opdf <- as_data_frame(data.frame(opf)) %>% select(-value.objfn,-iter,-fpevals,-objfevals)
  return(nest(opdf,par,feature_name))
}


result2beta <- function(result_df){
  aIntercept <- data_frame(feature="Intercept",Beta=mean(result_df$Intercept))
  rdf <- select(result_df,Beta,feature) %>% distinct() %>% bind_rows(aIntercept)
  Betavec <- rdf$Beta
  names(Betavec) <- rdf$feature
  return(Betavec)
}


anno2mat <- function(full_feat){
  feat_mat <-select(full_feat,Gene,feature,value,BF) %>% spread(feature,value)  %>% filter(complete.cases(.))
  Genes <- select(feat_mat,Gene)
  feat_mat <- mutate(Genes,Intercept=1) %>% inner_join(feat_mat) %>% select(-Gene)
  BF <-feat_mat$BF
  tmu <-BF/(BF+exp(-(-5+1.5*log(BF))))
  feat_mat <- data.matrix(select(feat_mat,-BF))
  rownames(feat_mat) <- Genes$Gene
  return(feat_mat)
}

FGEM <-function(fBeta,feat_mat,BF,null_features="Intercept"){
  retdf <- EM_mat(fBeta,feat_mat,BF) %>% mutate(nac=NA)
  retdf <- EM_mat(fBeta[paste0("feat_mat",null_features)],feat_mat[,null_features,drop=F],BF = BF) %>%
    select(NullLogLik=LogLik) %>% mutate(nac=NA)  %>% inner_join(retdf,by="nac") %>% select(-nac)
  retdf <- mutate(retdf,Chisq=-2*(NullLogLik-LogLik),
                  pval=pchisq(Chisq,df=ncol(feat_mat)-1,lower.tail=F))
  return(retdf)
}


 FGEM_bootstrap <- function(full_feat,scale=F,prior_mean=0.02,iterations=100){
   require(dplyr)
   require(tidyr)
   require(foreach)
   feat_mat <-select(full_feat,Gene,feature,value,BF) %>% spread(feature,value)  %>% filter(complete.cases(.))
   Genes <- select(feat_mat,Gene)
   feat_mat <- mutate(Genes,Intercept=1) %>% inner_join(feat_mat) %>% select(-Gene)

   BF <-feat_mat$BF
   tmu <-(prior_mean*BF)/((prior_mean*BF)+(1-prior_mean))
   feat_mat <- data.matrix(select(feat_mat,-BF))
   feat_mat <- scale(feat_mat,scale,scale)
   if(sum(is.na(feat_mat))==nrow(feat_mat)){
     feat_mat[is.na(feat_mat)] <- 1
   }
   fBeta <- coefficients(glm(tmu~feat_mat+0,family=quasibinomial(link="logit")))

   while(any(is.na(fBeta))){
     feat_mat <- feat_mat[,!is.na(fBeta)[-1]]
     fBeta <- coefficients(glm(tmu~feat_mat+0,family=quasibinomial(link="logit")))
   }
  foreach(i=1:iterations,.combine = "bind_rows") %do%{
    subsamp <- sample(1:nrow(feat_mat),nrow(feat_mat),replace = T)
    bretdf <- FGEM(fBeta=fBeta,feat_mat=feat_mat[subsamp,],BF=BF[subsamp]) %>% mutate(iter=i)
  }
}

sem_df <-function(full_feat,scale=F,prior_mean=0.02){
  require(dplyr)
  require(tidyr)

  feat_mat <-select(full_feat,Gene,feature,value,BF) %>% spread(feature,value)  %>% filter(complete.cases(.))
  Genes <- select(feat_mat,Gene)
  feat_mat <- mutate(Genes,Intercept=1) %>% inner_join(feat_mat) %>% select(-Gene)

  BF <-feat_mat$BF
  tmu <-(prior_mean*BF)/((prior_mean*BF)+(1-prior_mean))
  feat_mat <- data.matrix(select(feat_mat,-BF))
  feat_mat <- scale(feat_mat,scale,scale)
  if(sum(is.na(feat_mat))==nrow(feat_mat)){
    feat_mat[is.na(feat_mat)] <- 1
  }
  fBeta <- coefficients(glm(tmu~feat_mat+0,family=quasibinomial(link="logit")))

  while(any(is.na(fBeta))){
    feat_mat <- feat_mat[,!is.na(fBeta)[-1]]
    fBeta <- coefficients(glm(tmu~feat_mat+0,family=quasibinomial(link="logit")))
  }
  retdf <- FGEM(fBeta = fBeta,feat_mat=feat_mat,BF=BF)


  return(retdf)
}

em_df <-function(data_df,prior_mean=0.02,scale=F,null_feature){

  stopifnot(nrow(datadf)>1,
            all(null_feature %in% colnames(data_df)),
            all(!is.na(data_df)))
  data_df <-mutate(data_df,Intercept=1)
  null_feature <- unique(c(null_feature,"Intercept"))
  BF <- data_df$BF
  tmu <-(prior_mean*BF)/((prior_mean*BF)+(1-prior_mean))

  feat_mat <- data.matrix(select(data_df,-one_of("BF","Gene","class")))
  feat_mat <- scale(feat_mat,scale,scale)
  if(sum(is.na(feat_mat))==nrow(feat_mat)){
    feat_mat[is.na(feat_mat)] <- 1
  }
  fBeta <- coefficients(glm(tmu~feat_mat+0,family=quasibinomial(link="logit")))

  while(any(is.na(fBeta))){
    feat_mat <- feat_mat[,!is.na(fBeta)[-1]]
    fBeta <- coefficients(glm(tmu~feat_mat+0,family=quasibinomial(link="logit")))
  }
  fretdf <- FGEM(fBeta = fBeta,feat_mat=feat_mat,BF=BF,null_features = null_feature)
  return(fretdf)
}





chisq_comp <- function(full_feat,prior=0.02){
  isbin <- length(unique(full_feat$value))==2
  if(isbin){
    return(chisq.test(full_feat$value,factor(quantile(full_feat$BF,prior))))
  }
}

gen_model <- function(feat_list,annodf,datadf,scale=F){
  retdf <- filter(annodf,feature %in% feat_list) %>%
    group_by(feature) %>%
    do(cfeat_df(.,datadf)) %>%
    ungroup() %>%
    do(sem_df(.,scale = scale))
  return(retdf)
}


posterior_results <- function(full_model,datadf,anno_df,scale=F){
    retdf <- select(full_model,Beta,Intercept,feature) %>% inner_join(anno_df) %>% do(gen_prior_df(.,datadf,scale=scale)) %>%
    mutate(old_prior=mean(prior))  %>% rename(new_prior=prior) %>% inner_join(datadf) %>%
    mutate(new_posterior=new_prior*BF/((new_prior*BF)+(1-new_prior)),
           old_posterior=old_prior*BF/((old_prior*BF)+(1-old_prior)),
           post_improvement=new_posterior-old_posterior)
  return(retdf)
}


#     feat_res <- group_by(full_feat,feature) %>% mutate(isbin=(length(unique(value))==2),
#                                            feat_bin=factor(ifelse(isbin,factor(value),value>quantile(value,prior))),
#                                            BF_bin=factor(BF>quantile(BF,prior))) %>%
# do(par=chisq.test(.$feat_bin,.$BF_bin,simulate.p.value = T)) %>% summarise(feature=feature[1],prior=prior,statistic=par$statistic,pval=par$p.value)
#     return(feat_res)
# }



