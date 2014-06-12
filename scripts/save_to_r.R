
library(R.utils)

sourceDirectory("SpatialPRo-master/R")

edir <- "data/"

spe <- SPExperiment(edir)

spe <- loadExp(spe)


save(spe, file="spe.RData")
