\docType{methods}
\name{neighbours}
\alias{neighbours}
\alias{neighbours,SPData-method}
\alias{neighbours,SPData-methods}
\alias{neighbours<-}
\alias{neighbours<-,SPData-methods}
\title{List of nearest neighbour readouts}
\usage{
neighbours(object)

neighbours(object) <- value

\S4method{neighbours}{SPData}(object)
}
\arguments{
  \item{object}{The SPData object to use}

  \item{value}{The new list of nearest neighbour readouts}
}
\description{
Returns a list of length nCells(object). The ith entry is
an n by m matrix, for cell i having n neighbours each of
which have nChannel(object) channels.

Set the list of nearest neighbour readouts.
}

