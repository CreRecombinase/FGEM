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

  rlct <- fgem::cv_relax_fgem(X,ebf,log_BF=FALSE ,alpha=.5,hess=TRUE,grad=TRUE,verbose=FALSE,nlambda=4,v=4)
  sctr <- fgem:::summarise_cv_lik(rlct)
  
})


testthat::test_that("uninformative BF + uninformative prior -> uninformative posterior",{
  
  BF <- 1
  hBF <- BF + 0.1
  lBF <- BF - 0.1
  
  prior <- 0.5  
  hprior <- prior + 0.1
  lprior <- prior - 0.1
  
  
  balanced_post <- fgem:::BF2posterior(BF,prior=prior,log_BF=FALSE)
  expect_equal(balanced_post,0.5)
  
  hbhp <- fgem:::BF2posterior(hBF,prior=hprior,log_BF=FALSE)
  expect_gt(hbhp,balanced_post)
  
  lblp <- fgem:::BF2posterior(lBF,prior=lprior,log_BF=FALSE)
  expect_lt(lblp,balanced_post)
  
  
  lbalanced_post <- fgem:::BF2posterior(log(BF),prior=log(prior),log_BF=TRUE)
  expect_equal(exp(lbalanced_post),0.5)
  
  
  lhbhp <- exp(fgem:::BF2posterior(log(hBF),prior=log(hprior),log_BF=TRUE))
  expect_equal(hbhp,lhbhp)
  
  
  llblp <- exp(fgem:::BF2posterior(log(lBF),prior=log(lprior),log_BF=TRUE))
  expect_equal(lblp,llblp)
  
})



