library(rhdf5)
library(Matrix)
library(methods)
library(Matrix)
library(dplyr)
library(SQUAREM)


#Given a sparse matrix, write it to hdf5 format
# sparseAnnowrite <- function(h5f,)

#Code to read sparse matrix annotations from an HDF5 file ()
sparseAnnoread <- function(h5f,datapath,rowpath="GeneList/",colpath){
  rowname <- h5read(h5f,rowpath)
  colname <- h5read(h5f,colpath)
  dat <-h5read(annof,datapath)
  datmat <- sparseMatrix(i=dat[,1],j=dat[,2],x=1,
                         dims=c(length(rowname),length(colname)),
                         dimnames=list(rowname,colname))
  return(datmat)
}
#Return a matrix with two columns,z_score, and Bayes Factor (genes are rownames)
cons_to_bf <- function(consmat,allgenes){
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


Logit_log_lik <- function(u,x,B){
  mx <- cbind(1,x)
  Beta <- coefficients(glm(u~x,family=quasibinomial(link="logit")))
  pvec  <- 1/(1+exp(-(mx%*%Beta)))
  uvec <- (pvec*B)/((pvec*B)+(1-pvec))
  return(-sum(uvec*log(pvec+B)+(1-uvec)*log(1-pvec)))
}


Beta_log_lik <- function(Beta,x,B){
  mx <- cbind(1,x)
  pvec  <- 1/(1+exp(-(mx%*%Beta)))
  uvec <- (pvec*B)/((pvec*B)+(1-pvec))
  return(-sum(uvec*log(pvec+B)+(1-uvec)*log(1-pvec)))
}

gen_p <- function(u,x,B){
  mx <- cbind(1,x)
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

gen_final_mat <- function(X1,GOmatfile,Consfile,BayesFactorFile){
  gox <- X1[grepl("GO",X1)]
  consx <- X1[!grepl("GO",X1)]
  GOmat <- readRDS(GOmatfile)
  consmat <- readRDS(Consfile)
  allgenes <- Bayesfile(BayesFactorFile)
  cbf <- cons_to_bf(consmat = consmat,allgenes = allgenes)
  consdf <- data_frame(Gene=rownames(cbf),mis_z=cbf[,"mis_z"],B=cbf[,"B"])
  tGOdf <- as_data_frame(data.frame(as.matrix(GOmat[,gox,drop=F]),check.names = F,stringsAsFactors = F))
  GOdf <- mutate(tGOdf,Gene=rownames(GOmat))
  findf <- inner_join(consdf,GOdf)
  fmat <- data.matrix(select(findf,-Gene,-B))
  rownames(fmat) <- findf$Gene
  return(list(fmat,B=findf$B))
}


fullsem <- function(fmat,B){
  pf <- squarem(par=B/(B+exp(-B)),fixptfn=Logit_FGEM,
                objfn=Logit_log_lik,x=fmat,B=B)
  return(pf$par)
}
 





