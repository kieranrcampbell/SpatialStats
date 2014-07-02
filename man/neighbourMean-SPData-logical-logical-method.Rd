\docType{methods}
\name{neighbourMean,SPData,logical,logical-method}
\alias{neighbourMean,SPData,logical,logical-method}
\title{Averages over the nearest neighbour cells, going from a
list of length cell to a cell by channel matrix}
\usage{
\S4method{neighbourMean}{SPData,logical,logical}(object, useWeights, normalise)
}
\arguments{
  \item{useWeights}{If TRUE then nearest neighbours are
  weighted by cell boundary size}

  \item{normalise}{If TRUE then then each channel is
  normalised to mean 0 sd 1}
}
\description{
Averages over the nearest neighbour cells, going from a
list of length cell to a cell by channel matrix
}

