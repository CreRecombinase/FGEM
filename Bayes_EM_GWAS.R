library(methods)
library(Matrix)
library(plyr)
library(SQUAREM)
library(arm)
# #nodes <- detectCores()
# 
# registerDoMC(nodes)
Exacff <- "~/VariantPrioritization/gene_exac.RDS"
Ex <- readRDS(Exacff)
rownames(Ex) <- Ex$gene
tallgenes <- allgenes[rownames(allgenes) %in% rownames(Ex),]
tEx <- Ex[rownames(tallgenes),]
tB <- tallgenes$B
         
 #GOmatfile <- "/mnt/gluster/home/nwknoblauch/Gene_Annotation/GOmat.RDS"
 #BayesFactorFile <- "/mnt/gluster/home/nwknoblauch/Gene_Annotation/TADA_ASC_SSC_results_Dec23.csv"
 #Startterm <- 1
 #Numterms <- 3

# Outdir <- "/mnt/gluster/home/nwknoblauch/FGEM_GO_RDS/"
#Read from the command line to get which iteration to use
args <- commandArgs(trailingOnly=T)
GOmatfile <- args[1]
BayesFactorFile <- args[2]
Startterm <- as.integer(args[3])
Numterms <- as.integer(args[4])
Outdir <- args[5]

outfile <- file.path(Outdir,paste0("FGEM_DF_",Startterm,"_",Startterm+Numterms,".RDS"))
if(!all(file.exists(c(GOmatfile,BayesFactorFile)))){
    stop(paste0("File not found: ",c(GOmatfile,BayesFactorFile)[!file.exists(c(GOmatfile,BayesFactorFile))]))
}
GOmat <- readRDS(GOmatfile)
print("Size of Annotation Matrix is ")
print(dim(GOmat))
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
## likfun <- function(x,B){
#  sum(log((1/(1+exp(-x)))*B+(1-(1/(1+exp(-x))))))
#}
#NullLogLikelihood <- optimize(likfun,interval = c(-10,10),B=B,maximum=T)[["objective"]]

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

Logit_FGEM <- function(u,x,B,useBayes=T){
  
  mx <- cbind(1,x)
  if(useBayes){
    Beta <- coefficients(bayesglm(u~x,family=quasibinomial(link="logit")))
  }else{
    Beta <- coefficients(glm(u~x,family=quasibinomial(link="logit")))
  }
  pvec  <- 1/(1+exp(-(mx%*%Beta)))
  uvec <- (pvec*B)/((pvec*B)+(1-pvec))
  return(uvec)
}


Logit_log_lik <- function(u,x,B,useBayes=T){
    mx <- cbind(1,x)
    if(useBayes){
      Beta <- coefficients(bayesglm(u~x,family=quasibinomial(link="logit")))
    }else{
      Beta <- coefficients(glm(u~x,family=quasibinomial(link="logit")))
    }
    pvec  <- 1/(1+exp(-(mx%*%Beta)))
    uvec <- (pvec*B)/((pvec*B)+(1-pvec))
    return(-sum(log(uvec*(pvec+B)+(1-uvec)*(1-pvec))))
}


Beta_log_lik <- function(Beta,x,B){
    mx <- cbind(1,x)
#    Beta <- coefficients(glm(u~x,family=quasibinomial(link="logit")))
    pvec  <- 1/(1+exp(-(mx%*%Beta)))
    uvec <- (pvec*B)/((pvec*B)+(1-pvec))
    return(-sum(log(uvec*(pvec+B)+(1-uvec)*(1-pvec))))
}

Logit.uvec <- apply(GOmat,2,function(x,B){
    pf <- squarem(par=B/(B+exp(-B)),fixptfn=Logit_FGEM,
                  objfn=Logit_log_lik,x=x,B=B)
    return(pf$par)},B=B)
exac.B <- squarem(par=tB/(tB+exp(-tB)),fixptfn=Logit_FGEM,
                  objfn=Logit_log_lik,x=tEx$mis_z,B=tB)

eB <- bayesglm(exac.B$par~tEx$mis_z,family=quasibinomial(link="logit"))
Bmat <- matrix(0,ncol(Logit.uvec),4)
colnames(Bmat) <- c("Intercept","Beta","Likelihood","Null")
rownames(Bmat) <- colnames(Logit.uvec)

for(i in 1:nrow(Bmat)){
    tbeta <- gen_p(Logit.uvec[,i],GOmat[,i],B)
    Bmat[i,"Intercept"] <- tbeta[1]
    Bmat[i,"Beta"] <- tbeta[2]
    Bmat[i,"Likelihood"] <- Beta_log_lik(tbeta,GOmat[,i],B)
    Bmat[i,"Null"] <- Beta_log_lik(c(tbeta[1],0),GOmat[,i],B)
}
Bmat <- data.frame(Bmat)
    

print("Done!")
#rownames(Logit.uvec) <- rownames(GOmat)
saveRDS(Bmat,outfile)



