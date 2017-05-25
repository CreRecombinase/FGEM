
True_Logit <- function(Beta,x,XH1,XH0){
  pvec  <- 1/(1+exp(-(x%*%Beta)))
  uvec <- (pvec*XH1)/((pvec*XH1)+(1-pvec)*XH0)
  Beta <- coefficients(glm(uvec~x+0,family=quasibinomial(link="logit")))
  return(Beta)
}

True_Logit_log_lik <- function(Beta,x,XH1,XH0){
  pvec  <- 1/(1+exp(-(x%*%Beta)))
  # opvec  <- 1/(1+exp((mx%*%Beta)))
  uvec <- (pvec*XH1)/((pvec*XH1)+(1-pvec)*XH0)
  return(-sum(uvec*(log(pvec)+log(XH1))+(1-uvec)*(log(1-pvec)+log(XH0))))
}


True_EM_mat <-function(fBeta,feat_matrix,XH1,XH0){
  require(SQUAREM)
  require(dplyr)
  opf <- squarem(par=fBeta,fixptfn=True_Logit,x=feat_matrix,XH1=XH1,XH0=XH0)
  cn <- colnames(feat_matrix)
  opf$LogLik <- -True_Logit_log_lik(opf$par,feat_matrix,XH1,XH0)
  retdf <-data.frame(opf) %>% select(-value.objfn) %>% mutate(Intercept=opf$par[1],feature=cn) %>% rename(Beta=par)
  if(nrow(retdf)>1){
    retdf <- filter(retdf,feature!="Intercept")
  }
  return(retdf)
}

true_est <- function(fBeta,A,XH1,XH0){
  retdf <- True_EM_mat(fBeta,A,XH1,XH0) %>% mutate(nac=NA)
  retdf <- True_EM_mat(fBeta[1],A[,"Intercept",drop=F],XH1,XH0) %>%
    select(NullLogLik=LogLik) %>% mutate(nac=NA) %>%
    inner_join(retdf) %>% select(-nac)
  retdf <- mutate(retdf,Chisq=-2*(NullLogLik-LogLik),
                  pval=pchisq(Chisq,df=ncol(A)-1,lower.tail=F))
}


sim_dat <- function(tprior,ngenes,samplesize,betap){
  Z <- rbinom(ngenes,size=1,prob=tprior)
  x <-  rbinom(n = ngenes,size = samplesize,prob = ifelse(Z==1,rbeta(n = ngenes,shape1 = betap[1],shape2 = betap[2]),0.5))
  XH1 <-choose(samplesize,x)*((beta(a = x+betap[1],b = samplesize-x+betap[2])/beta(betap[1],betap[2])))
  XH0 <- dbinom(x,size=samplesize,prob=0.5)
  BF <-XH1/XH0

  return(data_frame(Z=Z,x=x,XH1=XH1,XH0=XH0,BF=BF,tprior=c(tprior),
                    mu=(tprior*XH1)/((tprior*XH1)+(1-tprior)*XH0)))
}

t_loglik <- function(sample_data){
  summarise(sample_data,lik=sum(mu*(log(tprior)+log(XH1))+(1-mu)*(log(1-tprior)+log(XH0))))
}


testb <- function(mu,A){
  return(coefficients(glm(mu~A+0,family=quasibinomial(link="logit"))))
}



comp_model <- function(tprior,betavec,betap,ngenes,samplesize,A,i){
  sample_data <- sim_dat(tprior,ngenes,samplesize,betap)
  sample_data <- mutate(sample_data,tmu= (gprior*BF)/((gprior*BF)+(1-gprior)),
                        bmu=(tprior*BF)/((tprior*BF)+(1-gprior)))
  tfvec <-coefficients(glm(sample_data$mu~A+0,family=quasibinomial(link="logit")))
  bfvec <-coefficients(glm(sample_data$bmu~A+0,family=quasibinomial(link="logit")))
  return(data_frame(Intercept=betavec[1],Beta=betavec[2],Intercept_t=tfvec[1],
                    Beta_t=tfvec[2],Intercept_b=bfvec[1],Beta_b=bfvec[2],i=i))
}





gen_est <- function(tprior,betap,ngenes,samplesize,A,i,marg=FALSE){
  gprior <- mean(tprior)
  sample_data <- sim_dat(tprior,ngenes,samplesize,betap)
  sample_data <- mutate(sample_data,tmu= (gprior*BF)/((gprior*BF)+(1-gprior)),
                        bmu=(tprior*BF)/((tprior*BF)+(1-gprior)))
  if(marg){
    retdfl <- list()
    for(j in 2:(ncol(A))){
      fbeta <- coefficients(glm(sample_data$tmu~A[,c(1,j)]+0,family=quasibinomial(link="logit")))
      retdfl[[j]] <- EM_mat(fBeta = fbeta,feat_matrix = A[,c(1,j)],BF=BF)
    }
    retdf <- bind_rows(retdfl)
  }else{
    fbeta <- coefficients(glm(sample_data$tmu~A+0,family=quasibinomial(link="logit")))
    retdf <- EM_mat(fBeta = fbeta,feat_matrix = A,BF = BF)
  }
  retdf <- mutate(retdf,iter=i,gprior=gprior,marg=marg)
  return( retdf)
}

gen_anno <- function(annofunc,annofeat,ngenes){
  A <- mapply(function(x,y){
    do.call(x,c(list(n=ngenes),y))
  },annofunc,annofeat)
  return(A)
}


simfunc <- function(ngenes=20000,samplesize=1000,iterations=100,betavec,
                    annofunc=list("Intercept"=rbinom,"a"=runif),
                    annofeat=list("Intercept"=list(size=1,prob=1),"a"=list(min=0,max=1)),betap=c(0,1,0.1),marg=FALSE){
  require(dplyr)
  require(ggplot2)
  require(foreach)
  require(doParallel)
  library(parallel)
  betap <- c(0.1,0.1)
  nfeat <- length(betavec)
  A <- gen_anno(annofunc,annofeat,ngenes)
  tprior <- gen_p(Beta = betavec,x=A)
  gprior <- mean(tprior)
  cl <- makeCluster(8)
  registerDoParallel(cl)
  pdf <- foreach(i= 1:iterations,.combine = 'bind_rows',.packages = c("FGEM","dplyr")) %dopar% gen_est(tprior = tprior,betap = betap,ngenes = ngenes,samplesize=samplesize,A=A,i=i,marg=marg)
  stopCluster(cl)
  return(pdf)
}
