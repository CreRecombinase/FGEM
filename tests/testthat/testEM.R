context("simulation")
set.seed(123)
test_that("FGEM works on simulated data",{
  g <- 1000
  sz <- 200
  f <- 15
  ps <- 0.5
  X <- do.call(cbind,purrr::map(runif(f,min = 0.3,max=0.7),~as.numeric(rbinom(n = g,size=1,prob=.x))))
  
  eff <- c(-5,abs(rnorm(f,mean=0.7))*rbinom(f,size=1,prob=ps))
  fn <- c("Intercept",paste0("V",seq_len(f)))
  true_value <- tibble::tibble(feature_name=fn,true_Beta=eff)  
  pr <- fgem:::gen_p(eff,X)
  y <- rbinom(g,size=1,prob=pr)
  q <- dplyr::if_else(y==1,0.5,runif(g))
  xg <- rbinom(g,size=sz,prob = q)
  pxz1 <- dbinom(x = xg,size = sz,prob=0.5)
  pxz2 <- beta((xg+1),sz-xg+1)*choose(sz,xg)
  bf <- log(pxz1)-log(pxz2)
  ebf <- exp(bf)

#   new_grad <- function(par,X,eBF){
#       myenv <- new.env()
#       assign("Beta",par, envir = myenv)
# #      assign("Beta0",par[1] , envir = myenv)
#       assign("eBF", eBF, envir = myenv)
#       assign("X", X, envir = myenv)
#       rd <- numericDeriv(quote(sum(log(eBF + (1 - eBF)/t(1 + exp(X%*%Beta[-1]+Beta[1]))))),c("Beta"),myenv)
#       return(c(attr(rd,"gradient")))
#   }
  plan(multisession)
  tr <- fgem:::fgem_bfgs(X,bf,log_BF=TRUE ,alpha=.5,lambda=lambda,hess=TRUE,grad=TRUE,verbose=FALSE)
  ctr <- fgem::cv_fgem(X,ebf,log_BF=FALSE ,alpha=.5,hess=TRUE,grad=TRUE,verbose=FALSE,nlambda=4,v=4)
  rlct <- cv_relax_fgem(X,ebf,log_BF=FALSE ,alpha=.5,hess=TRUE,grad=TRUE,verbose=FALSE,nlambda=4,v=4)
  sctr <- fgem:::summarise_cv_lik(ctr)
  check_lik <- purrr::pmap_dbl(tr,function(Beta,l1,l2,...){
    B <- Beta$Beta
    l <- -fgem::fgem_lik(B,X,bf,log_BF=TRUE)
    l1p <- fgem:::l1_norm_penalty(B,l1 = l1)
    l2p <- fgem:::l2_norm_penalty(B,l2 = l2)
    return(-(l+l1p+l2p))
  })
  
  tB0 <- tr$Beta[[1]]$Beta
  check_grad <- purrr::pmap(tr,function(l1,l2,...){
    g <- fgem::fgem_grad(tB0,X,bf,log_BF=TRUE,l1=l1,l2=0)
    return(g)
  })
  gmat0 <- do.call(cbind,check_grad)
  
  rgrd <- map(fgem::fgem_grad(par = tr$Beta[[1]]$Beta,X,bf,l2=0,l1=tr$l1[1],log_BF=TRUE)[-1])
  
  
  
  mutate(tr,check_lik=check_lik) %>% filter(convergence==0,lik<max(lik)) %>% ggplot(aes(x=lik,y=check_lik))+geom_point()
  fgem:::fgem_lik(tr$Beta[[3]]$Beta,X,bf,0,0,log_BF=TRUE)
  new_lik(tr$Beta[[3]]$Beta,X,exp(bf))
  tg <- fgem:::fgem_grad(tr$Beta[[3]]$Beta,X,bf,0,0,log_BF=TRUE)
   ebf <- exp(bf)
  cg <- new_grad(par,X,eBF)

  mb <- microbenchmark::microbenchmark(  nl = fgem:::fgem_lik(tr$Beta[[3]]$Beta,X,ebf,0,0,log_BF=FALSE),
                                         ll = fgem:::fgem_lik(tr$Beta[[3]]$Beta,X,bf,0,0,log_BF=TRUE))
    
  
  wl <- which(tr$l0n>1)[1]
  l <- lambda[wl]
  prec <- tr$l2[wl]
  l1 <- tr$l1[wl]
  l <- lambda[wl]
  fgl <- fgem:::fgem_elasticnet(X,bf,alpha=0.5,lambda=l,log_BF=TRUE,grad=TRUE,verbose=TRUE)
  
  


  gr <- fgem::fgem_grad(fgl$Beta[[1]]$Beta,X,BF=bf,log_BF=TRUE,prec=fgl$l2,l1=fgl$l1)
  sctg <- sign(gr)
  i <- 4
  par <- fgl$Beta[[1]]$Beta
  prec <- fgl$l2
  l1 <- fgl$l1
  BF <- bf
  log_BF <- TRUE
  lik_idx <- function(i,par,X,BF,log_BF,prec,l1,sdd=2){
    par[i] <- rnorm(1,mean=par[i],sd=sdd)
    nlik <- fgem_lik(par,X,BF,prec=prec,log_BF=log_BF,l1=l1)
    uplik <-fgem_lik(par,X,BF,prec=0.0,log_BF=log_BF,l1=0.0)
    gr <- fgem_grad(par,X,BF,prec=prec,log_BF=log_BF,l1=l1)
    upgr <- fgem_grad(par,X,BF,prec=0.0,log_BF=log_BF,l1=0.0)
    return(tibble::tibble(par=par[i],lik=nlik,uplik=uplik,grad=gr[i],upgrad=upgr[i]))
  }
  table(y[X[,1]==1])
  tdf <- purrr::map_dfr(rep(2,1000),lik_idx,par=par,X=X,BF=BF,log_BF=log_BF,prec=prec,l1=l1,sdd=.5)
  ggplot2::ggplot(tdf,ggplot2::aes(x=par,y=lik,col=grad))+ggplot2::geom_point()+ggplot2::geom_vline(xintercept=fgl$Beta[[1]]$Beta[2])+ggplot2::scale_color_gradient2()
  ggplot(tdf,aes(x=par,y=uplik,col=upgrad))+geom_point()+geom_vline(xintercept=fgl$Beta[[1]]$Beta[4])+scale_color_gradient2()
  
 
})




