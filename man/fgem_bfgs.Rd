% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FGEM_Func.R
\name{fgem_bfgs}
\alias{fgem_bfgs}
\title{fit fgem using lbfgs}
\usage{
fgem_bfgs(
  X,
  BF,
  Beta0 = c(0, rep(0, ncol(X))),
  verbose = FALSE,
  alpha = 1,
  lambda = c(2, 1, 5 * 10^(-(seq(1, 6, length.out = 75))), 0),
  add_lambda = FALSE,
  progress = FALSE,
  log_BF = FALSE,
  ...
)
}
\arguments{
\item{X}{feature matrix}

\item{BF}{vector of bayes factors}

\item{Beta0}{initial guess (defaults to 0)}

\item{verbose}{(whether to give verbose output)}

\item{alpha}{alpha (as in glmnet). Alpha  is the (scalar) proportion of
`lambda` applied to l1 penalization, while `1-alpha` is applied to l2}

\item{lambda}{vector of shrinkage parameters}

\item{...}{currently unused}
}
\value{
tibble with results of fit
}
\description{
fit fgem using lbfgs
}
