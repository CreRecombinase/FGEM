library(testthat)
library(FGEM)
x <-cbind(rep(1,3))
B <- rep(2,3)
Beta <- c(2)
FGEM_Logit(Beta = Beta,x = x,B = B)


test_that("FGEM estimation machinery works with uninformative BF and a single (constant) annotation",{
  expect_equivalent(FGEM_Logit(Beta = 1,x = cbind(rep(1,3)),B = rep(1,3)),1)
  expect_equivalent(FGEM_Logit(Beta = 2,x = cbind(rep(1,3)),B = rep(1,3)),2)
  expect_equivalent(FGEM_Logit(Beta = -3.5,x = cbind(rep(1,3)),B = rep(1,3)),-3.5)
  expect_equivalent(FGEM_Logit(Beta = 0,x = cbind(rep(1,3)),B = rep(1,3)),0)
})





test_that("FGEM estimation machinery works with uninformative BF and two (binary, continuous, or constant) annotations",{
  expect_equivalent(FGEM_Logit(Beta=c(-2,1),
                               x=cbind(rep(c(1),6),rep(c(0,1),3)),
                               B=rep(1,6)),c(-2,1))
  expect_equivalent(FGEM_Logit(Beta=c(1,-2),
                               x=cbind(rep(c(1),6),rep(c(0,1),3)),
                               B=rep(1,6)),c(1,-2))
  expect_equivalent(FGEM_Logit(Beta=c(-2,1),
                               x=cbind(runif(100),runif(100)),
                               B=rep(1,100)),c(-2,1))
  expect_equivalent(FGEM_Logit(Beta=c(-2,1),
                               x=cbind(rep(c(0,1),50),runif(100)),
                               B=rep(1,100)),c(-2,1))
})








