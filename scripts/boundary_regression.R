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


doLMTest <- function(sp,tumourID=NULL,alpha=0.05, useWeights=TRUE) {
    #source("parse-nn.R")
    set.seed(31)

###############################################################################
    ## We want to regress the response variables (responseSubset) in the tumour  ##
    ## (group 2) with EMT interesting values on the phospho signalling molecules ##
    ## in the stromal region (group 1)                                           ##
###############################################################################

    ## regressing only on phosphates, so
    phInd <- c(3,4,13,14,15,17,18,27)
    phNames <- channels(sp)[phInd]
    print(phNames)

    Y <- cells(sp)
    X <- neighbourMean(sp, useWeights, TRUE)
    X <- X[,phInd]
    colnames(X) <- phNames

    ### want only cells that lie well within the tumour

    responseCells <- NULL
    if(!is.null(tumourID)) {
        cl <- which(cellClass(sp) == tumourID)
        responseCells <- setdiff(cl, findBoundary(sp))
    } else {
        responseCells <- 1:nCells(sp)
    }

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

    X <- X[responseCells,]

    ## finally add colnames and construct linear model

    fit <- lm(Y ~ X)

    ntrial <- 1000

    ## let's sample 80% of cells

    nCellSample <- round(0.7 * dim(Y)[1])

    res <- matrix(0, nrow=dim(X)[2], ncol=dim(Y)[2])
    colnames(res) <- responseNames

    rownames(res) <- colnames(X)

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


plsr <- function(sp, tumourID=NULL, useWeights=TRUE) {
    library(plsdepot)
    phInd <- c(3,4,13,14,15,17,18,27)
    phNames <- channels(sp)[phInd]
    print(phNames)

    Y <- cells(sp)
    X <- neighbourMean(sp, useWeights, TRUE)
    X <- X[,phInd]
    colnames(X) <- phNames

    ### want only cells that lie well within the tumour

    if(!is.null(tumourID)) {
        cl <- which(cellClass(sp) == tumourID)
        responseCells <- setdiff(cl, findBoundary(sp))
    } else {
        responseCells <- 1:nCells(sp)
    }

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

    X <- X[responseCells,]


    pls <- plsreg2(X, Y)
    return(pls)
}
