


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

EM_mat <-function(Beta0,feat_mat,BF){
  opf <- SQUAREM::squarem(par=Beta0,fixptfn=FGEM_Logit,x=feat_mat,B=BF)
  cn <- colnames(feat_mat)
  opf$LogLik <- -FGEM_Logit_log_lik(opf$par,feat_mat,BF)
  opf$feature_name <- names(opf$par)
  opdf <- dplyr::as_data_frame(data.frame(opf))%>% dplyr::select(-value.objfn,-iter,-fpevals,-objfevals) %>% rename(Beta=par)
  return(tidyr::nest(opdf,Beta,feature_name))
}


FGEM_Logit_log_lik <- function(Beta,x,B){
  pvec  <- 1/(1+exp(-(x%*%Beta)))
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





#' Generate prior from fitted FGEM_df result and annotation dataframe
#' @param feat_df annotation dataframe (as described inf `FGEM_df`)
#' @param fgem_result_df result of running FGEM_df
gen_prior <- function(feat_df,fgem_result_df){
  
  
  Beta0 <- t(as.matrix(tidyr::unnest(fgem_result_df) %>% select(Beta,feature_name) %>% spread(feature_name,Beta)))
  # BF <- feat_df[["BF"]]
  feat_df <- mutate(feat_df,Intercept=1)
  data_mat <- data.matrix(feat_df[,rownames(Beta0)])
  prior_n <- c(gen_p(Beta = Beta0,x = data_mat))
  return(mutate(feat_df,prior=prior_n) %>% select(-Intercept))
}


#' Generate posterior from fitted FGEM_df result and annotation dataframe
#' @param feat_df annotation dataframe (as described inf `FGEM_df`)
#' @param fgem_result_df result of running FGEM_df
gen_posterior <- function(feat_df,fgem_result_df){
  Beta0 <- t(as.matrix(tidyr::unnest(fgem_result_df) %>% select(Beta,feature_name) %>% spread(feature_name,Beta)))
  BF <- feat_df[["BF"]]
  feat_df <- mutate(feat_df,Intercept=1)
  data_mat <- data.matrix(feat_df[,rownames(Beta0)])
  post_n <- c(gen_u(Beta = Beta0,x = data_mat,B = BF))
  return(mutate(feat_df,posterior=post_n) %>% select(-Intercept))
}

pmean <-function(Beta,feat_mat){
  pvec  <- 1/(1+exp(-(feat_mat%*%Beta)))
  return(mean(pvec))
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


 







