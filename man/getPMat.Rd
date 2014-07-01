\name{getPMat}
\alias{getPMat}
\title{Get p-value matrix from fit}
\usage{
getPMat(spdata, fit, setOne = FALSE, alpha = 0.05)
}
\arguments{
  \item{spdata}{The SPdata object the fit was made from}

  \item{fit}{The fit returned by calling nnReg}

  \item{setZero}{Sets all entries to zero that don't pass
  the threshold}

  \item{alpha}{The significance level for the tests.}
}
\description{
This function constructs a nChannel by nChannel directed
matrix of p values
}

