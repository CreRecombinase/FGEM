

sim_dat <-
gen_est <- function(tprior,betap,ngenes,samplesize,A,i,marg=FALSE){
  gprior <- mean(tprior)
  Z <- rbinom(ngenes,size=1,prob=tprior)
  x <-  rbinom(n = ngenes,size = samplesize,prob = ifelse(Z==1,rbeta(n = ngenes,shape1 = betap[1],shape2 = betap[2]),0.5))
  BF <- (beta(a = x+betap[1],b = samplesize-x+betap[2])/beta(betap[1],betap[2]))/dbinom(x,size=samplesize,prob=0.5)
  tmu <- (gprior*BF)/((gprior*BF)+(1-gprior))

  if(marg){
    retdfl <- list()
    for(j in 1:(ncol(A)-1)){
      fbeta <- coefficients(glm(tmu~A[,c(1,j)]+0,family=quasibinomial(link="logit")))
      retdfl[[j]] <- EM_mat(fBeta = fbeta,feat_matrix = A[,c(1,j)],BF=BF)
    }
    retdf <- bind_rows(retdfl)
  }else{
    fbeta <- coefficients(glm(tmu~A+0,family=quasibinomial(link="logit")))
    retdf <- EM_mat(fBeta = fbeta,feat_matrix = A,BF = BF)
  }
  retdf <- mutate(retdf,iter=i,gprior=gprior,marg=marg)
  return( retdf)
}


simfunc <- function(ngenes=20000,samplesize=1000,iterations=100,betavec,
                    annofunc=list("Intercept"=rbinom,"a"=runif),
                    annofeat=list("Intercept"=list(size=1,prob=1),"a"=list(min=-1,max=1)),marg=FALSE){
  require(dplyr)
  require(ggplot2)
  require(foreach)
  require(doParallel)
  library(parallel)
  betap <- c(0.5,0.5)
  nfeat <- length(betavec)
  A <- mapply(function(x,y){
    do.call(x,c(list(n=ngenes),y))
  },annofunc,annofeat)
  tprior <- gen_p(Beta = betavec,x=A)
  pdfl <- list()
  gprior <- mean(tprior)
  cl <- makeCluster(8)
  registerDoParallel(cl)
  pdf <- foreach(i= 1:iterations,.combine = 'bind_rows',.packages = c("FGEM","dplyr")) %dopar% gen_est(tprior = tprior,betap = betap,ngenes = ngenes,samplesize=samplesize,A=A,i=i,marg=marg)
  stopCluster(cl)
  return(pdf)
}
