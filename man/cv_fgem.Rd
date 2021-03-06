% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FGEM_Func.R
\name{cv_fgem}
\alias{cv_fgem}
\title{fit cross-validated FGEM model}
\usage{
cv_fgem(
  X,
  BF,
  log_BF = FALSE,
  alpha = 1,
  nlambda = 100,
  lambda = NULL,
  stratify_BF = TRUE,
  v = 10,
  ...
)
}
\arguments{
\item{X}{Feature matrix}

\item{BF}{vector of bayes factors}

\item{alpha}{alpha (as in glmnet). Alpha  is the (scalar) proportion of
`lambda` applied to l1 penalization, while `1-alpha` is applied to l2}

\item{lambda}{vector of shrinkage parameters}
}
\value{

}
\description{
fit cross-validated FGEM model
}
