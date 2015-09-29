library(methods)
library(Matrix)
library(plyr)
# #nodes <- detectCores()
# 
# registerDoMC(nodes)
                                        
# GOmatfile <- GOmatfile
# BayesFactorFile <- Bayes
# Starterm <- 1001
# Numterms <- 100
# EMiter <- 50
# Outdir <- "~/temp"
#Read from the command line to get which iteration to use
args <- commandArgs(trailingOnly=T)
GOmatfile <- args[1]
BayesFactorFile <- args[2]
Startterm <- as.integer(args[3])
Numterms <- as.integer(args[4])
EMiter <- as.integer(args[5])
Outdir <- args[6]

outfile <- file.path(Outdir,paste0("FGEM_DF_",Startterm,"_",Startterm+Numterms,".RDS"))
if(!all(file.exists(c(GOmatfile,BayesFactorFile)))){
    stop(paste0("File not found: ",c(GOmatfile,BayesFactorFile)[!file.exists(c(GOmatfile,BayesFactorFile))]))
}
GOmat <- readRDS(GOmatfile)
BF.df <- read.table(BayesFactorFile,header=T,sep=",")


#To start the EM algorithm, we need a starting guess for Beta. To do this we will use 
#logistic regression on thresholded Bayes Factors (Bayes Factor>3) 

allgenes <- data.frame(Gene=BF.df$Gene,B=BF.df$BF,stringsAsFactors=F)
rownames(allgenes) <- allgenes$Gene
allgenes <- allgenes[rownames(GOmat),]

#Pull the relevant rows from the dataframe

B <- allgenes[rownames(GOmat),"B"]


####Starting EM
#Computing the null model
likfun <- function(x,B){
  sum(log((1/(1+exp(-x)))*B+(1-(1/(1+exp(-x))))))
}
NullLogLikelihood <- optimize(likfun,interval = c(-10,10),B=B,maximum=T)[["objective"]]

#If the number of terms to be tested is greater than the number of remaining terms
# we change the number of terms to be tested
if(Startterm+Numterms>ncol(GOmat)){
  Numterms <- ncol(GOmat)-Startterm
}

#Subset the GO matrix to only the terms we're going to test
GOmat <- GOmat[,Startterm:(Startterm+Numterms)]


#Function that performs EM
FGEM <- function(x,B,Z=NULL,iters=NULL,tol=NULL,keep=F,NullLogLikelihood=NULL){
  #x is a vector of length I with relevant annotations
  #B is a vector of length I with Bayes Factors
  #Z is a vector of length I with the initial guess for membership (1's or 0's)
  ## Z can be NULL, in which case Z will be randomly sampled from a binomial with probabilities according to the Bayes Factor.
  
  #iters is the number of EM iterations.  If not specified, EM will be run untill convergence specified by tol.
  
  if(is.null(NullLogLikelihood)){
    likfun <- function(x,B){
      sum(log((1/(1+exp(-x)))*B+(1-(1/(1+exp(-x))))))
    }
    NullLogLikelihood  <- optimize(likfun,interval = c(-10,10),B=B,maximum=T)$objective
  }
  
  if(is.null(tol)&is.null(iters)){
    stop("Must specify tol or iters")
  }
  if(ifelse(is.null(Z),
            length(x)!=length(B),
            (length(x)!=length(B))|(length(x)!=length(Z))|(length(B)!=length(Z)))){
    stop("Length of x and B (and Z, if specified) must be equal")
  }
  
  if(is.null(Z)){
    pB <- ecdf(B)(B)
    sigcount <- sum(pB)/sum(B>10)
    Z <- rbinom(length(B),size = 1,prob = pB/sigcount)
  }
  if(!is.null(iters)&is.null(tol)){
    tol <- 0
  }
  if(!is.null(tol)&is.null(iters)){
    iters <- 1000
  }
  if(keep){
    matdim <- iters
    Beta <- matrix(0,2,matdim)
    Beta[,1] <- coefficients(glm(Z~x,family=binomial))
    pvec <- 1/(1+exp(-(Beta[1,1]+Beta[2,1]*x)))
    LogLikelihood <- numeric(matdim)
    LogLikelihood[1] <- sum(log(pvec*B+(1-pvec)))
    i <- 1
    while ((i < iters)&(i<1000)){
     #This is the mixture proportion
      uvec  <- (pvec*B)/((pvec*B)+(1-pvec))
      #This is where we optimize the Q function (in this case, using logistic regression)
      Beta[,i] <- coefficients(glm(uvec~x,family=quasibinomial(link="logit")))
      pvec  <- 1/(1+exp(-(Beta[1,i]+Beta[2,i]*x)))
      #If the difference between the current log likelihood and the previous log likelihood is below the tolerance, we quit
      rettol <- sum(log(pvec*B+(1-pvec)))-LogLikelihood[i]
      if(rettol<tol){
        LogLikelihood[i] <- sum(log(pvec*B+(1-pvec)))
        break
      }
      LogLikelihood[i] <- sum(log(pvec*B+(1-pvec)))
      #If we're about to stop and we still aren't better than the null model, we keep going (up to 1k iterations)
      if(i==(iters-1)){
        if(LogLikelihood[i]<NullLogLikelihood){
          iters <- iters+10
        }
      }
      i <- i+1
    }
    Chisq <- -2*(NullLogLikelihood-LogLikelihood[matdim])
    retlist <- list(Intercept=Beta[1,matdim],Beta=Beta[2,matdim],
                    LogLikelihood=LogLikelihood[matdim],
                    NullLogLikelihood=NullLogLikelihood,
                    Chisq=Chisq,tol=rettol)
  }
  else{
    Beta <- coefficients(glm(Z~x,family=binomial))
    pvec <- 1/(1+exp(-(Beta[1]+Beta[2]*x)))
    LogLikelihood <- sum(log(pvec*B+(1-pvec)))
    i <- 1
    while ((i < iters)&(i<1000)){
      #This is the mixture proportion
      uvec  <- (pvec*B)/((pvec*B)+(1-pvec))
      #This is where we optimize the Q function (in this case, using logistic regression)
      Beta <- coefficients(glm(uvec~x,family=quasibinomial(link="logit")))
      pvec  <- 1/(1+exp(-(Beta[1]+Beta[2]*x)))
      #If the difference between the current log likelihood and the previous log likelihood is below the tolerance, we quit
      rettol <- sum(log(pvec*B+(1-pvec)))-LogLikelihood
      if(rettol<tol){
        LogLikelihood <- sum(log(pvec*B+(1-pvec)))
        break
      }
      LogLikelihood <- sum(log(pvec*B+(1-pvec)))
      #If we're about to stop and we still aren't better than the null model, we keep going (up to 1k iterations)
      if(i==(iters-1)&LogLikelihood<NullLogLikelihood){
        iters <- iters+10
      }
      i <- i+1
    }
    Chisq <- -2*(NullLogLikelihood-LogLikelihood)
    retlist <- list(Intercept=Beta[1],Beta=Beta[2],
                    LogLikelihood=LogLikelihood,
                    NullLogLikelihood=NullLogLikelihood,
                    Chisq=Chisq,tol=rettol)
  } 
  if(keep){
    mats <- list(Betas=Beta,LogLikelihood=LogLikelihood)
    return(list(retlist,mats))
  }else{
    return(data.frame(retlist))
  }
}

FGEM.df <- adply(GOmat,2,FGEM,B=B,iters=50,keep=F,NullLogLikelihood=NullLogLikelihood,.parallel=F)
print("Done!")
saveRDS(FGEM.df,outfile)



