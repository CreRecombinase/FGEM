#Compare methods
GOmatfile <- "~/Dropbox/BayesianDA/FGEM/Data/GOmat.RDS"
Bfile <- "~/Dropbox/BayesianDA/FGEM/Data/TADA_ASC_SSC_results_Dec23.csv"
BFile <- "TADA_ASC_SSC_results_Dec23.csv"
outfile <- "GOmat.RDS"
BF.G <- read.table(Bfile,header=T,sep=",")


allgenes <- data.frame(Gene=BF.G$Gene,B=BF.G$BF,stringsAsFactors=F)
rownames(allgenes) <- allgenes$Gene

Exacff <- "~/Dropbox/BayesianDA/RRDS/gene_exac.RDS"
Ex <- readRDS(Exacff)
rownames(Ex) <- Ex$gene
tallgenes <- allgenes[rownames(allgenes) %in% rownames(Ex),]
tEx <- Ex[rownames(tallgenes),]
tB <- tallgenes$B

adat <- readRDS("~/Dropbox/CompBio/badd.RDS")
bdat <- readRDS("~/Dropbox/CompBio/bdd.RDS")
GOmat <- readRDS(GOmatfile)
bdat <- data.frame(collect(bdat))
adat <- data.frame(collect(adat))
rownames(bdat) <- bdat$id
rownames(adat) <- adat$id
bdat <- bdat[rownames(adat),]

bdat$chisq <- -bdat$chisq
bdat$pch <- pchisq(bdat$chisq,df=1,lower.tail = F)
bdat$pad <- p.adjust(bdat$pch,method="fdr")
plot(bdat$Beta,adat$Beta)
plot(bdat$Likelihood,adat$Likelihood)
plot(bdat$chisq,adat$chisq)

bdat$name <- Gon[bdat$id,"name"]


fbdat <- head(bdat[order(bdat$pch),],40)

fGOm <- cbind(as.matrix(GOmat[rownames(tEx),rownames(fbdat)]),Exac=tEx$mis_z)



Logit_FGEM <- function(u,x,B){
  mx <- cbind(1,x)
  Beta <- coefficients(glm(u~x,family=quasibinomial(link="logit")))
  pvec  <- 1/(1+exp(-(mx%*%Beta)))
  uvec <- (pvec*B)/((pvec*B)+(1-pvec))
  return(uvec)
}

Logit_log_lik <- function(uvec,x,B){
  mx <- cbind(1,x)
  Beta <- coefficients(glm(uvec~x,family=quasibinomial(link="logit")))
  pvec  <- 1/(1+exp(-(mx%*%Beta)))
 # uvec <- (pvec*B)/(A(pvec*B)+(1-pvec))

  return(-sum(log(uvec*(pvec+B)+(1-uvec)*(1-pvec))))
}

lik <- numeric(10)
mx <- cbind(1,x)
uvec <- B/(B+exp(-B))
for(i in 1:10){
  Beta <- coefficients(glm(uvec~x,family=quasibinomial(link="logit")))
  pvec  <- 1/(1+exp(-(mx%*%Beta)))
  uvec <- (pvec*B)/((pvec*B)+(1-pvec))
  lik[i] <- -sum(log(uvec*(pvec+B)+(1-uvec)*(1-pvec)))
  print(lik[i])
}


tx <- crossprod(fGOm)
ee <- eigen(tx)
evals <- zapsmall(ee$values)
evecs <- zapsmall(ee$vectors)
nx <- fGOm[,-which(evecs[,which(evals==0)]!=0)]
bgfg <- squarem(par=B/(B+exp(-B)),fixptfn=Logit_FGEM,objfn=Logit_log_lik,x=fGOm,B=B,control=list(trace=TRUE))

gen_p <- function(u,x,B){
  Beta <-bayesglm(u~x,family=quasibinomial(link="logit"))
  return(Beta)
}

fbeta <- gen_p(bgfg$par,fGOm,B)
nB <- coef(fbeta)
nX <- cbind(1,fGOm)
pvec  <- 1/(1+exp(-(mx%*%Beta)))
plot(B,bgfg$par)
predict(fbeta)
coef(summary(fbeta))




