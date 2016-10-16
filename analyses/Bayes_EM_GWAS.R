library(methods)
library(Matrix)
library(dplyr)
library(SQUAREM)
library(FGEM)

# #nodes <- detectCores()
#
# registerDoMC(nodes)
# Exacff <- "~/VariantPrioritization/gene_exac.RDS"
# Ex <- readRDS(Exacff)
# rownames(Ex) <- Ex$gene
#
# tallgenes <- allgenes[rownames(allgenes) %in% rownames(Ex),]
# tEx <- Ex[rownames(tallgenes),]
# tB <- tallgenes$B
#
#funcfile <- "~/Dropbox/BayesianDA/FGEM/FGEM_Func.R"
#GOmatfile <- "~/Dropbox/BayesianDA/FGEM/Data/ConstraintMat.RDS"
GOmatfile <- "~/Dropbox/BayesianDA/FGEM/Data/BrainSpanVar.RDS"
BayesFactorFile <- "~/Dropbox/BayesianDA/FGEM/Data/TADA_ASC_SSC_results_Dec23.csv"
outfile <- "~/Dropbox/BayesianDA/FGEM/Data/Expression_Res.RDS"
 #GOmatfile <- "/mnt/gluster/home/nwknoblauch/Gene_Annotation/GOmat.RDS"
 #BayesFactorFile <- "/mnt/gluster/home/nwknoblauch/Gene_Annotation/TADA_ASC_SSC_results_Dec23.csv"
 #Startterm <- 1
 #Numterms <- 3

# Outdir <- "/mnt/gluster/home/nwknoblauch/FGEM_GO_RDS/"
Outdir <- "~/Dropbox/BayesianDA/FGEM/Data/FGEM_GO_RDS/"

#Read from the command line to get which iteration to use
args <- commandArgs(trailingOnly=T)
GOmatfile <- args[1]
BayesFactorFile <- args[2]
Startterm <- as.integer(args[3])
Numterms <- as.integer(args[4])
Outdir <- args[5]
funcfile <- args[6]

source(funcfile)
outfile <- file.path(Outdir,paste0("FGEM_DF_",Startterm,"_",Startterm+Numterms,".RDS"))
if(!all(file.exists(c(GOmatfile,BayesFactorFile)))){
    stop(paste0("File not found: ",c(GOmatfile,BayesFactorFile)[!file.exists(c(GOmatfile,BayesFactorFile))]))
}
GOmat <- readRDS(GOmatfile)
print("Size of Annotation Matrix is ")
print(dim(GOmat))
allgenes <- Bayesfile(BayesFactorFile)
intergenes <- intersect(rownames(GOmat),rownames(allgenes))
GOmat <- GOmat[intergenes,]
allgenes <- allgenes[intergenes,]



#To start the EM algorithm, we need a starting guess for Beta. To do this we will use
#logistic regression on thresholded Bayes Factors (Bayes Factor>3)


# allgenes <- allgenes[rownames(GOmat),]

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

Logit.uvec <- apply(GOmat,2,function(x,B){
  pf <- squarem(par=B/(B+exp(-B)),fixptfn=Logit_FGEM,
                x=x,B=B)
  return(pf$par)},B=B)


Bmat <- matrix(0,ncol(Logit.uvec),4)
colnames(Bmat) <- c("Intercept","Beta","LogLikelihood","NullLogLikelihood")
rownames(Bmat) <- colnames(Logit.uvec)


for(i in 1:nrow(Bmat)){
  tbeta <- gen_p(Logit.uvec[,i],GOmat[,i],B)
  Bmat[i,"Intercept"] <- tbeta[1]
  Bmat[i,"Beta"] <- tbeta[2]
  Bmat[i,"LogLikelihood"] <- Beta_log_lik(tbeta,GOmat[,i],B)
  Bmat[i,"NullLogLikelihood"] <- Beta_log_lik(c(tbeta[1],0),GOmat[,i],B)
}
Bmat <- as_data_frame(Bmat)
Bmat$Chisq <- -2*(Bmat$NullLogLikelihood-Bmat$LogLikelihood)
Bmat <- mutate(Bmat,feature=colnames(Logit.uvec),method="objfn")
print("Done!")
#rownames(Logit.uvec) <- rownames(GOmat)
saveRDS(Bmat,outfile)



