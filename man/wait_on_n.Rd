% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FGEM_Func.R
\name{wait_on_n}
\alias{wait_on_n}
\title{take a list of futures and block until at least `n` future is resolved}
\usage{
wait_on_n(x, n = 1, pause = 0.02)
}
\arguments{
\item{x}{container of futures}
}
\value{
x
}
\description{
take a list of futures and block until at least `n` future is resolved
}
