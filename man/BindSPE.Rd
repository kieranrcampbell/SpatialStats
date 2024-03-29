% Generated by roxygen2 (4.0.1): do not edit by hand
\name{BindSPE}
\alias{BindSPE}
\title{Bind multiple sample data together}
\usage{
BindSPE(SPE, choose.class = NULL, use.weights = FALSE, normalize = TRUE)
}
\arguments{
\item{SPE}{The SPExp object to use}

\item{choose.class}{Choose a particular class of cells. Should be a vector of length
\code{length{SPE}} where the ith entry indicates the class of cells in sample i to pick.}

\item{use.weights}{Passed to neighbourMeans to optionally weight the regression by relative boundary size}

\item{normalize}{If true predictor columns are centre-scaled}
}
\value{
A list with three components:
\itemize{
\item{\code{X} }{The bound predictor matrix}
\item{\code{Y} }{The bound response matrix}
\item{\code{sizes} }{A vector with the number of cells selected from each \code{SPE}}
}
}
\description{
The low power and high noise in some analyses makes it necessary to bind multiple
samples together. Given an SPExp object, this method separates out each cell type,
calculate neighbour means and binds the samples together, i.e. the individual response
matrices for different samples \eqn{Y1}, \eqn{Y2}, \eqn{Y3} are row-bound together to
form one large response matrix \eqn{Y = [Y1' Y2' Y3']'}.
}
\details{
If \code{normalize = TRUE} it is necessary to introduce sample factors as constructed
by \code{\link{ConstructSampleFactors}}, otherwise the regression estimates will be affected
by the relative means of different samples.
}

