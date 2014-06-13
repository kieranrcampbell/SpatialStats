\name{getAdj}
\alias{getAdj}
\title{Get adjacency matrix from fit}
\usage{
getAdj(spdata, fit, alpha = 0.05)
}
\arguments{
  \item{spdata}{The SPdata object the fit was made from}

  \item{fit}{The fit returned by calling nnReg}

  \item{alpha}{The significance level for the tests.}
}
\description{
This function constructs a nProt by nProt directed
adjacency matrix from the fit
}

