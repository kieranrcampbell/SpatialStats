\name{weightedSubsetBoundaryRegression}
\alias{weightedSubsetBoundaryRegression}
\title{Performs weighted (cell edge weights) regression along the tumour-stromal
boundary for a given set of response channels and a class (tumour or stromal)
to focus on in the reponses}
\usage{
weightedSubsetBoundaryRegression(sp, cellClasses, responseClass, boundary,
  weights, responseSubset, dependentSubset)
}
\arguments{
  \item{sp}{The SPData class for the experiment}

  \item{cellClasses}{A vector of length nCells(sp) where
  each entry is either 1 or 2}

  \item{responseClass}{A numeric 1 or 2 indicating which
  cells should be included in the response variable}

  \item{boundary}{A vector with the indices of cells that
  lie along the beoundary}

  \item{weights}{A list of length nCells where each entry
  gives the length of the boundary between it and all
  nearest neighbour cells}

  \item{responseSubset}{A vector containing the indices of
  interesting proteins to use in the response variables}

  \item{dependentSubset}{A vector containing the indices of
  proteins to include in regression coefficients}
}
\description{
Performs weighted (cell edge weights) regression along the
tumour-stromal boundary for a given set of response
channels and a class (tumour or stromal) to focus on in the
reponses
}

