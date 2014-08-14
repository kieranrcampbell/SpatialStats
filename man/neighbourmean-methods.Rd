\docType{methods}
\name{neighbourMean}
\alias{neighbourMean}
\alias{neighbourMean,SPData,logical,logical-method}
\alias{neighbourMean,SPData-methods}
\title{Average over nearest neighbours}
\usage{
neighbourMean(object, useWeights, normalise)

\S4method{neighbourMean}{SPData,logical,logical}(object, useWeights, normalise)
}
\arguments{
  \item{object}{The SPData object to use}

  \item{useWeights}{If TRUE then the mean is weighted using
  the relative boundary weights. If FALSE then the normal
  mean is taken over the nearest neighbours.}

  \item{normalise}{If TRUE each channel is centre-scaled to
  mean 0 and standard deviation 1.}
}
\value{
A matrix of dimension \emph{n} by \emph{p} for \emph{n}
cells and \emph{p} nearest neighbours.
}
\description{
Given the nearest neighbour measurements, average over the
nearest neighbours of each cell, with optional boundary
weights.
}

