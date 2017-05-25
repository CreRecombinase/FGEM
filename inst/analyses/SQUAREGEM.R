#FGEM using SQUAREM
library(SQUAREM)
library(methods)
library(Matrix)
library(plyr)

setwd("~/Dropbox/BayesianDA/FGEM/Data/")
GOmatfile <- "GOmat.RDS"
BayesFactorFile <- "TADA_ASC_SSC_results_Dec23.csv"
Startterm <- 1
Numterms <- 15
EMiter <- 50
featureClassName <- "GOterms"
Outdir <- "./"
outfile <- file.path(Outdir,paste0("FGEM_DF_",featureClassName,"_",Startterm,"_",Startterm+Numterms,".RDS"))
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
print("Subsetting GO matrix")
GOmat <- GOmat[,Startterm:(Startterm+Numterms)]
print("GO matrix is now of size ")
print(dim(GOmat))


# Do a single step of EM, take in p and B, give back new pvec
fgem.em <- function(pvec,B,annodat){
  uvec  <- (pvec*B)/((pvec*B)+(1-pvec))
  #This is where we optimize the Q function (in this case, using logistic regression)
  Beta <- coefficients(bayesglm(uvec~annodat+0,family=quasibinomial(link="logit")))
  pvec  <- 1/(1+exp(-(mx%*%Beta)))
  return(pvec)
}
fgem.obj <- function(pvec,B,annodat){
  LogLikelihood <- sum(log(pvec*B+(1-pvec)))
  return(-LogLikelihood)
}

Beta <- matrix(0,ncol(mx),1)
Beta[] <- coefficients(bayesglm(Z~mx+0,family=binomial))
pvec <- 1/(1+exp(-(mx%*%Beta)))




#Fixed point EM
pf1 <- fpiter(pvec, B=B,annodat=mx, fixptfn=fgem.em, objfn=fgem.obj)
pf2 <- squarem(pvec,B=B,annodat=mx,fixptfn = fgem.em,objfn=fgem.obj)

finp <- pf1$par
sinp <- pf2$par
plot(finp,sinp)

head(finp)
plot(finp[,1],pvec[,1])

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



