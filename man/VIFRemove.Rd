\name{VIFRemove}
\alias{VIFRemove}
\title{Takes a matrix of predictors and sequentially removes them until the largest
VIF is below a threshold (default=10)}
\usage{
VIFRemove(X, vif.threshold = 10)
}
\arguments{
  \item{X}{Matrix of predictors}

  \item{vif.threshold}{Maximum VIF threshold allowed}
}
\description{
Takes a matrix of predictors and sequentially removes them
until the largest VIF is below a threshold (default=10)
}

