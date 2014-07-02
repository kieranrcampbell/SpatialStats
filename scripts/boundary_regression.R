library(EMCluster)
library(devtools)
library(rgl)
library(R.matlab)
library(gplots)

load_all("..")

#load("~/ebi/sp/data/sp5.RData")


clust <- cellClass(sp)

cl1 <- which(clust == 1)
cl2 <- which(clust == 2)

Y <- cells(sp)

## find out which is tumour and which is stromal
kerindex <- grep("Keratin",channels(sp))
kerReads <- Y[,kerindex]
tumourID <- which.max( c(mean(kerReads[cl1]), mean(kerReads[cl2])))

boundary <- NULL


doLMTest <- function(sp, classWithin,clust, alpha=0.01) {
    source("parse-nn.R")

###############################################################################
    ## We want to regress the response variables (responseSubset) in the tumour  ##
    ## (group 2) with EMT interesting values on the phospho signalling molecules ##
    ## in the stromal region (group 1)                                           ##
###############################################################################

    ## regressing only on phosphates, so
    phInd <- c(3,4,13,14,15,17,18,27)
    phNames <- channels(sp)[phInd]
    print(phNames)




    protein.names <- channels(sp)
    getProteinIds <- function(proteinNames,proteinList) {
        p.indices <- sapply(proteinList, grep, proteinNames)
        p.indices
    }

    responseNames <- c("Cadherin",
                       "bcat",
                       "Vimentin",
                       "CD44","Twist",
                       "Slug", "NFkB",
                       "EGFR")

    responseSubset <- getProteinIds(protein.names, responseNames)


    ## select out response channels
    Y <- Y[,responseSubset]

    ## select out response cells
    Y <- Y[responseCells,]
    NN <- NN[responseCells]

    ## average over nearest neighbours
    X <- lapply(NN, function(nn) {
        if(is.matrix(nn)) colMeans(nn) else nn
    })

    X <- matrix(unlist(X), ncol=length(phInd), byrow=TRUE)

    ## finally add colnames and construct linear model

    fit <- lm(Y ~ X)

    ntrial <- 100

    ## let's sample 80% of cells

    nCellSample <- round(0.8 * dim(Y)[1])

    res <- matrix(0, nrow=8, ncol=8)
    colnames(res) <- responseNames
    rownames(res) <- phNames

    for(n in 1:ntrial) {
        s <- sample(1:(dim(Y)[1]), nCellSample)
        y <- Y[s,] ; x <- X[s,]

        fit <- lm(y ~ x)
        suma <- summary(fit)

        for(i in 1:8) {
            co <- suma[[i]]$coefficients
            co <- co[-1,4]
            co <- as.numeric(co < alpha)
            res[,i] <- res[,i] + co
        }

    }
    return(res)
}
