


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
  # uvec <- (pvec*B)/((pvec*B)+(1-pvec))
  return(-sum(log(pvec*B+(1-pvec))))
  # return(-sum(uvec*(log(pvec)+log(B))+(1-uvec)*log(1-pvec)))
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
  retdf <- gather(annodf,feature,value,-Gene)  %>% dplyr::mutate(class=feat.name)
}

exp2df<- function(expmat,feat.name="Expression"){
  require(dplyr)
  require(tidyr)
  Genes <-rownames(expmat)
  expdf <- as_data_frame(expmat) %>% dplyr::mutate(Gene=Genes) %>% gather(feature,value,-Gene) %>% dplyr::mutate(class=feat.name)
  return(expdf)
}

glist2df <- function(Genes,feat.name){
  require(dplyr)
  retdf <- data_frame(Gene=Genes) %>% dplyr::mutate(feature=feat.name,value=1,class=feat.name)
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
  annotation_df <- group_by(annotation_df,feature) %>% dplyr::mutate(isBin=(n_distinct(value)==1),ngenes=n()) %>%ungroup()
  mingenes <-filter(annotation_df,isBin==FALSE) %>% group_by(Gene) %>% summarise(nfeat=n_distinct(feature)) %>% ungroup() %>% filter(nfeat==max(nfeat))
  if(all(annotation_df$isBin)){
    mingenes <-distinct(annotation_df,Gene)
  }
  annotation_df <- annotation_df %>% inner_join(mingenes,by="Gene") %>% select(-one_of(c("isBin","ngenes"))) %>% spread(feature,value,fill = 0)
  datadf <- select(datadf,Gene,BF)
  full_feat <- inner_join(datadf,annotation_df,by="Gene")
  return(full_feat)
}

gen_prior_df <- function(annodf,datadf,scale=F){
  tfeatdf <- dplyr::distinct(annodf,feature,Beta,Intercept)
  tannodf <- dplyr::select(annodf,-Beta,-Intercept) %>% dplyr::group_by(feature) %>% dplyr::do(cfeat_df(.,datadf))  %>%  dplyr::ungroup() %>%
    dplyr::inner_join(tfeatdf)
  am <-anno2mat(tannodf)
  bv <- result2beta(annodf)
  mprior <- 1/(1+exp(-am%*%bv[colnames(am)]))

  retdf <- dplyr::group_by(tannodf,Gene) %>% dplyr::summarise(prior=1/(1+exp(-(mean(Intercept)+sum(value*Beta))))) %>% dplyr::filter(Gene %in% rownames(am))
  return(retdf)
}


pmean <-function(Beta,feat_mat){
  pvec  <- 1/(1+exp(-(feat_mat%*%Beta)))
  return(mean(pvec))
}



EM_mat <-function(Beta0,feat_mat,BF){
  opf <- SQUAREM::squarem(par=Beta0,fixptfn=FGEM_Logit,x=feat_mat,B=BF)
  cn <- colnames(feat_mat)
  opf$LogLik <- -FGEM_Logit_log_lik(opf$par,feat_mat,BF)
  opf$feature_name <- names(opf$par)
  opdf <- dplyr::as_data_frame(data.frame(opf))%>% dplyr::select(-value.objfn,-iter,-fpevals,-objfevals) %>% rename(Beta=par)
  return(tidyr::nest(opdf,Beta,feature_name))
}


result2beta <- function(result_df){
  aIntercept <- data_frame(feature="Intercept",Beta=mean(result_df$Intercept))
  rdf <- dplyr::select(result_df,Beta,feature) %>% distinct() %>% bind_rows(aIntercept)
  Betavec <- rdf$Beta
  names(Betavec) <- rdf$feature
  return(Betavec)
}


anno2mat <- function(full_feat){
  feat_mat <-dplyr::select(full_feat,Gene,feature,value,BF) %>% spread(feature,value)  %>% filter(complete.cases(.))
  Genes <- dplyr::select(feat_mat,Gene)
  feat_mat <- dplyr::mutate(Genes,Intercept=1) %>% inner_join(feat_mat) %>% dplyr::select(-Gene)
  BF <-feat_mat$BF
  tmu <-BF/(BF+exp(-(-5+1.5*log(BF))))
  feat_mat <- data.matrix(dplyr::select(feat_mat,-BF))
  rownames(feat_mat) <- Genes$Gene
  return(feat_mat)
}


#' Run FGEM with a dataframe as input
#' This version of FGEM takes a single dataframe with 
#' both gene-level bayes factors and gene-level features 
#' and returns a dataframe with effect size estimates 
#' for every feature.
#' @param feat_df dataframe with data.  `feat_df` must have a column called `BF` specifying bayes factors.
#'   Any additional columns (besides an optional gene-name column that must be titled `Gene`) are treated as annotations.
#' @param prior_mean scalar between 0 and 1 specifying the starting 
#' prior probability (this is for initializing the EM algorithm). If unspecified, it defaults to 0.2
#' @param null_features this is a character vector specifying which features should be included in the null model.
#' The default is to only include an intercept term.
#' @param Beta0 vector specifying starting guesses for effect sizes of annotations. 
#' In almost all cases the default (NULL) should be kept



FGEM_df <- function(feat_df,prior_mean=0.02,null_features="Intercept",Beta0=NULL){
  
  BF <- unlist(dplyr::select(feat_df,BF))
  
  tmu <-(prior_mean*BF)/((prior_mean*BF)+(1-prior_mean))

  if(any(!null_features %in% colnames(feat_df))){
    warning("Some (or all) of null_features not found in feat_df, creating Intercept column\n ")
    feat_df <- dplyr::mutate(feat_df,Intercept=1)
    null_features <- "Intercept"
  }
  data_mat_df <- dplyr::select(feat_df,-Gene,-BF) %>% mutate(tmu=tmu)
  data_mat_df <- mutate_all(data_mat_df,function(x){
    if(all(is.na(x))){
    return(rep(1,length(x)))
           }else{
             return(x)
           }
  })

  if(is.null(Beta0)){
    Beta0 <- coefficients(glm(tmu~.+0,data=data_mat_df,family=quasibinomial(link="logit")))    
  }
  
  while(any(is.na(Beta0))){
    bad_betas <- names(Beta0)[is.na(Beta0)]
    data_mat_df <- dplyr::select(data_mat_df,-one_of(bad_betas[1]))
    # cat("Removing NA Feature,",sum(is.na(Beta0))," still NA\n")
    # data_mat <- data_mat[,!is.na(Beta0)[-1]]
    Beta0 <- coefficients(glm(tmu~.+0,data=data_mat_df,family=quasibinomial(link="logit")))
  }
  data_mat <- data.matrix(dplyr::select(data_mat_df,-tmu))
  ret <- FGEM(Beta0 = Beta0,feat_mat=data_mat,BF=BF,null_features = null_features)
  return(ret)
}


#' Run FGEM with a matrix as input
#' This version of FGEM takes a single dataframe with 
#' both gene-level bayes factors and gene-level features 
#' and returns a dataframe with effect size estimates 
#' for every feature.
#' @param feat_df dataframe with data.  `feat_df` must have a column called `BF` specifying bayes factors.
#'   Any additional columns are treated as annotations.
#' @param prior_mean scalar between 0 and 1 specifying the starting 
#' prior probability (this is for initializing the EM algorithm). If unspecified, it defaults to 0.2
#' @param null_features this is a character vector specifying which features should be included in the null model.
#' The default is to only include an intercept term.
#' @param Beta0 vector specifying starting guesses for effect sizes of annotations. 
#' In almost all cases the default (NULL) should be kept
#'
#'
#' 
FGEM <-function(Beta0,feat_mat,BF,null_features="Intercept"){
  retdf <- EM_mat(Beta0,feat_mat,BF) %>% dplyr::mutate(nac=NA)
  retdf <- EM_mat(Beta0[null_features],feat_mat[,null_features,drop=F],BF = BF)%>%
    dplyr::select(NullLogLik=LogLik) %>% dplyr::mutate(nac=NA)  %>% dplyr::inner_join(retdf,by="nac") %>% dplyr::select(-nac)
  retdf <- dplyr::mutate(retdf,Chisq=-2*(NullLogLik-LogLik),
                  pval=pchisq(Chisq,df=ncol(feat_mat)-1,lower.tail=F))
  return(retdf)
}


 FGEM_bootstrap <- function(full_feat,scale=F,prior_mean=0.02,iterations=100){
   require(dplyr)
   require(tidyr)
   require(foreach)
   feat_mat <-select(full_feat,Gene,feature,value,BF) %>% spread(feature,value)  %>% filter(complete.cases(.))
   Genes <- select(feat_mat,Gene)
   feat_mat <- dplyr::mutate(Genes,Intercept=1) %>% inner_join(feat_mat) %>% select(-Gene)

   BF <-feat_mat$BF
   tmu <-(prior_mean*BF)/((prior_mean*BF)+(1-prior_mean))
   feat_mat <- data.matrix(select(feat_mat,-BF))
   feat_mat <- scale(feat_mat,scale,scale)
   if(sum(is.na(feat_mat))==nrow(feat_mat)){
     feat_mat[is.na(feat_mat)] <- 1
   }
   Beta0 <- coefficients(glm(tmu~feat_mat+0,family=quasibinomial(link="logit")))

   while(any(is.na(Beta0))){
     feat_mat <- feat_mat[,!is.na(Beta0)[-1]]
     Beta0 <- coefficients(glm(tmu~feat_mat+0,family=quasibinomial(link="logit")))
   }
  foreach(i=1:iterations,.combine = "bind_rows") %do%{
    subsamp <- sample(1:nrow(feat_mat),nrow(feat_mat),replace = T)
    bretdf <- FGEM(Beta0=Beta0,feat_mat=feat_mat[subsamp,],BF=BF[subsamp]) %>% dplyr::mutate(iter=i)
  }
}

 

 
 
sem_df <-function(full_feat,scale=F,prior_mean=0.02){
  require(dplyr)
  require(tidyr)

  feat_mat <-select(full_feat,Gene,feature,value,BF) %>% spread(feature,value) # %>% filter(complete.cases(.))
  Genes <- select(feat_mat,Gene)
  feat_mat <- dplyr::mutate(Genes,Intercept=1) %>% inner_join(feat_mat) %>% select(-Gene)

  BF <-feat_mat$BF
  tmu <-(prior_mean*BF)/((prior_mean*BF)+(1-prior_mean))
  feat_mat <- data.matrix(select(feat_mat,-BF))
  feat_mat <- scale(feat_mat,scale,scale)
  if(sum(is.na(feat_mat))==nrow(feat_mat)){
    feat_mat[is.na(feat_mat)] <- 1
  }
  Beta0 <- coefficients(glm(tmu~feat_mat+0,family=quasibinomial(link="logit")))
  
  while(any(is.na(Beta0))){
    feat_mat <- feat_mat[,!is.na(Beta0)[-1]]
    Beta0 <- coefficients(glm(tmu~feat_mat+0,family=quasibinomial(link="logit")))
  }
  retdf <- FGEM(Beta0 = Beta0,feat_mat=feat_mat,BF=BF)


  return(retdf)
}

em_df <-function(data_df,prior_mean=0.02,scale=F,null_feature){

  stopifnot(nrow(datadf)>1,
            all(null_feature %in% colnames(data_df)),
            all(!is.na(data_df)))
  data_df <-dplyr::mutate(data_df,Intercept=1)
  null_feature <- unique(c(null_feature,"Intercept"))
  BF <- data_df$BF
  tmu <-(prior_mean*BF)/((prior_mean*BF)+(1-prior_mean))

  feat_mat <- data.matrix(select(data_df,-one_of("BF","Gene","class")))
  feat_mat <- scale(feat_mat,scale,scale)
  if(sum(is.na(feat_mat))==nrow(feat_mat)){
    feat_mat[is.na(feat_mat)] <- 1
  }
  Beta0 <- coefficients(glm(tmu~feat_mat+0,family=quasibinomial(link="logit")))

  while(any(is.na(Beta0))){
    feat_mat <- feat_mat[,!is.na(Beta0)[-1]]
    Beta0 <- coefficients(glm(tmu~feat_mat+0,family=quasibinomial(link="logit")))
  }
  fretdf <- FGEM(Beta0 = Beta0,feat_mat=feat_mat,BF=BF,null_features = null_feature)
  return(fretdf)
}


chisq_comp <- function(full_feat,prior=0.02){
  isbin <- length(unique(full_feat$value))==2
  if(isbin){
    return(chisq.test(full_feat$value,factor(quantile(full_feat$BF,prior))))
  }
}
#sfinal_model <- gen_model(feat_list = unique(sig_anno$feature),annodf = sig_anno,datadf = datadf,scale = F)
gen_model <- function(feat_list,annodf,datadf,scale=F){
  anno_mat <- filter(annodf,feature %in% feat_list)
  n_annomat <- spread(anno_mat,feature,value,fill=NA)
  inp_df <- anno_mat%>% group_by(feature) %>% do(cfeat_df(.,datadf)) %>%
    
  ungroup() 
    retdf <- filter(annodf,feature %in% feat_list) %>%
  group_by(feature) %>%
  do(cfeat_df(.,datadf)) %>%
  ungroup() %>%
    do(sem_df(.,scale = scale))
  return(retdf)
}


posterior_results <- function(full_model,datadf,anno_df,scale=F){
    retdf <- select(full_model,Beta,Intercept,feature) %>% inner_join(anno_df) %>% do(gen_prior_df(.,datadf,scale=scale)) %>%
    dplyr::mutate(old_prior=mean(prior))  %>% rename(new_prior=prior) %>% inner_join(datadf) %>%
    dplyr::mutate(new_posterior=new_prior*BF/((new_prior*BF)+(1-new_prior)),
           old_posterior=old_prior*BF/((old_prior*BF)+(1-old_prior)),
           post_improvement=new_posterior-old_posterior)
  return(retdf)
}






