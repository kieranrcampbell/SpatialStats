\docType{methods}
\name{[}
\alias{[}
\alias{[,SPData-methods}
\title{Subset an SPData set}
\usage{
\S4method{[}{SPData,ANY,ANY,ANY}(x, i, j)
}
\arguments{
  \item{i}{Cells to subset}

  \item{j}{Channels to subset}

  \item{x}{ object from which to extract element(s) or in
  which to replace element(s).  }
}
\value{
An SPData object reduced to cells \code{i} and channels
\code{j}
}
\description{
Select SPData[i,j] for cells \code{i} and channels
\code{j}. Note that this does not subset out nearest
neighbours also.
}

