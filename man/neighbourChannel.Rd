\name{neighbourChannel}
\alias{neighbourChannel}
\title{Filter out nearest neighbours by cell class. If cell i has no
nearest neighbours of class cell.class then numeric(0) is returned.}
\usage{
neighbourChannel(NN, channel.ids)
}
\arguments{
  \item{NN}{Neighbour list (e.g. via neighbours(sp))}

  \item{channel.ids}{A channel list to select out}
}
\description{
Filter out nearest neighbours by cell class. If cell i has
no nearest neighbours of class cell.class then numeric(0)
is returned.
}

