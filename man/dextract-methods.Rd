\docType{methods}
\name{[[}
\alias{[[}
\alias{[[,SPExp-methods}
\alias{[[<-}
\title{Access a sample}
\usage{
\S4method{[[}{SPExp}(x, i)
}
\arguments{
  \item{x}{The SPExp instance to use}

  \item{i}{The index of the \code{SPData} object within the
  \code{SPExp} to extract}

  \item{x}{The SPExp instance to use}

  \item{i}{The index of the \code{SPData} object within the
  \code{SPExp} to set}

  \item{value}{An object of class \code{SPData} to replace
  slot \code{i}}
}
\value{
The \code{SPData} object in slot \code{i} of the
\code{SPExp} instance
}
\description{
Extracts a single sample from a \code{SPExp} object

Sets a single sample in a \code{SPExp} object
}

