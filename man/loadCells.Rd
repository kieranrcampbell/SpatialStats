\name{loadCells}
\alias{loadCells}
\title{Loads an Xell matlab file into the SPData format}
\usage{
loadCells(filename, id = -1, control.isotopes = c("Xe131", "Cs133",
  "Ir193"), rescale.data = FALSE, scale.factor = 10000)
}
\arguments{
  \item{filename}{The matlab file}

  \item{control.isotopes}{The isotopes used for control to
  exclude from analysis}

  \item{log.data}{Boolean of whether to log the data}
}
\description{
This function parses the matlab files, pulling out relevant
proteins in the 'D' channel and validates that the correct
proteins are present. Matlab files can be found at
https://s3.amazonaws.com/supplemental.cytobank.org/report_data/report_113/Figure_5/Figure_5_raw_image_files.zip
}

