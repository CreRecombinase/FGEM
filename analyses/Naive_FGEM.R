library(SQUAREM)
library(plyr)

Logit_FGEM <- function(u,x,B){
  
  mx <- cbind(1,x)
  Beta <- coefficients(glm(u~x,family=quasibinomial(link="logit")))
  pvec  <- 1/(1+exp(-(mx%*%Beta)))
  uvec <- (pvec*B)/((pvec*B)+(1-pvec))
  return(uvec)
}

Naive_FGEM <- function(u,x,B){
  #x is a vector of length I with relevant annotations
  #B is a vector of length I with Bayes Factors
  #Z is a vector of length I with the initial guess for membership (1's or 0's)
  ## Z can be NULL, in which case Z will be randomly sampled from a binomial with probabilities according to the Bayes Factor.
  
  ##If x is a vector, we'll add an intercept and make it a matrix. If it's a matrix, we'll just add an an intercept column
  n <- length(x)
  pa <- mean(x)
  pz <- mean(u)
  paz <- sum(u[x==1])/sum(u)
  prz=(pa*pz)/pa
  
  pvec <- prz^x*(1-prz)^(1-x)
  uvec <- (pvec*B)/((pvec*B)+(1-pvec))
  return(uvec)
}

log_lik <-function(u,x,B){
  n <- length(x)
  pa <- mean(x)
  pz <- mean(u)
  paz <- sum(u[x==1])/sum(u)
  prz=(pa*pz)/pa
  pvec <- prz^x*(1-prz)^(1-x)
  uvec <- (pvec*B)/((pvec*B)+(1-pvec))
  lik <- -sum(log(uvec*(pvec+B)+(1-uvec)*(1-pvec)))
  return(lik)
}


returnp <- function(u,x){
  n <- length(x)
  pa <- mean(x)
  pz <- mean(u)
  paz <- sum(u[x==1])/sum(u)
  prz=(pa*pz)/pa
  pvec <- prz^x*(1-prz)^(1-x)
  return(pvec)
}



lik <- numeric(ncol(GOmat))
for ( i in 1:ncol(GOmat)){
  if(i%%10==0){
    print(i)
  }
 lik[i]  <- log_lik(GOMr[i,],GOmat[,i],B)
}
prmat <- matrix(0,nrow(GOmat),ncol(GOmat))
for( i in 1:ncol(GOmat)){
  if( i %%100==0){
    print(i)
  }
  prmat[,i] <-returnp(GOMr[i,],GOmat[,i])
}
fpr <- rowSums(log(prmat))
bm <- apply(prmat,2,max)
bmin <- apply(prmat,2,min)
brat <- bm/bmin


hist(lik)

N <- 18735
pz <- 0.053
P <- 200
#
# paz <- rbeta(n = P,shape1 = 1/10,shape2 = 10)
# pa <- mean(colMeans(GOmat))
# pz <- 0.053
# annomat <- matrix(0,N,P)
# for(i in 1:P){
#   annomat[,i] <- rbinom(N,size=1,prob = pa)
# }
# Z <- numeric(N)
# pmat <- matrix(0,N,P)
# for(i in 1:P){
#   pmat[,i] <- paz[i]^annomat[,i]*(1-paz[i])^(1-annomat[,i])
# }
# amat <- matrix(0,N,P)
# for(i in 1:P){
#   amat[,i] <- pa^annomat[,i]*(1-pa)^(1-annomat[,i])
# }
# postZ <- log(pz)+rowSums(log(pmat)-log(amat))
# sum(exp(postZ))
# sum(exp(postZ))/N
#
# mB <- sample(B,size = N,replace=T,prob = exp(postZ)/sum(exp(postZ)))
#
#
# testZ <- runif(N)
# pf <- squarem(par=testZ,fixptfn = Naive_FGEM,objfn=log_lik,x=annomat[,1],B=mB)

GOmatfile <- "~/Dropbox/BayesianDA/FGEM/Data/GOmat.RDS"
GOMff <- "~/Dropbox/BayesianDA/FGEM/Data/oGOmat.RDS"
GOmat <- readRDS(GOmatfile)
Bfile <- "~/Dropbox/BayesianDA/FGEM/Data/TADA_ASC_SSC_results_Dec23.csv"
BF <- read.table(Bfile,header=T,sep=",")
B <- BF$BF
testZ <- 1*(B>5)
GOMr <- aaply(GOmat,2,function(x){
  pf <- squarem(par=B/sum(B),fixptfn = Naive_FGEM,objfn=log_lik,x=x,B=B)
  return(pf$par)},.progress = "text")

nullpost <- squarem(par=B/sum(B),fixptfn = Naive_FGEM,objfn=log_lik,x=rep(1,length(B)),B=B)


GOlik <- numeric(nrow(GOMr)) 
for(i in 1:length(GOlik)){
  if(i%%100==0){
    print(i)
  }
  GOlik[i] <- log_lik(GOMr[i,],GOmat[,i],B)
}

fprior <- numeric(nrow(GOMr))
for(i in 1:nrow(GOMr)){
  
  fprior <- 
}
trim.GOm <- GOMr[colSums(GOmat)>3,]
saveRDS(trim.GOm,GOMff)
GOMff <- "~/Dropbox/CompBio/GOMrf.RDS"


tpzs <- rowMeans(GOMr)
pzs <- colMeans(GOMr)
pas <- rowMeans(GOmat)

tGO <- readRDS("~/Dropbox/BayesianDA/FGEM/Data/FGEM_GO_RDS/FGEM_DF_1_151.RDS")

tgB <- tGO$Beta
plot(brat[1:nrow(tGO)],exp(tGO$Beta))
postlog <- matrix(0,nrow(GOmat),nrow(tGO))
1/(1+exp(-(mx%*%Beta)))
for(i in 1:ncol(postlog)){
  mx <- cbind(1,GOmat[,i])
  Beta <- c(tGO$Intercept[i],tGO$Beta[i])
  pvec <- 1/(1+exp(-(mx%*%Beta)))
  postlog[,i] <-(pvec*B)/((pvec*B)+(1-pvec))
}

