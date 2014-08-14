\docType{methods}
\name{cellClass}
\alias{cellClass}
\alias{cellClass,SPData-method}
\alias{cellClass,SPData-methods}
\alias{cellClass<-}
\alias{cellClass<-,SPData-methods}
\title{Classes of each cell}
\usage{
cellClass(object)

cellClass(object) <- value

\S4method{cellClass}{SPData}(object)
}
\arguments{
  \item{object}{The SPData object to use}
}
\value{
A numeric vector of cell classes
}
\description{
If cells within a given tissue have been classified into
different types then get/set the classes.

Classes of each cell
}

