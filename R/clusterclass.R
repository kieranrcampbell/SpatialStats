
########################################
## Tools for clustering by cell class ##
########################################

#' Clusters cells into different classes using
#' Expectation Maximisation (EM) clustering
#'
#' @param Y A cell by channel matrix of readouts
#' @param doPCA If true, first performs PCA and uses the top
#' 3 principle components for clustering
#'
#' @export
clusterClass <- function(Y, doPCA = TRUE, nclass=2) {
    require(EMCluster)
    d <- NULL
    if(doPCA) {
        ypca <- princomp(Y)
        d <- ypca$scores[,1:3]
    } else {
        d <- Y
    }

    emobj <- simple.init(d, nclass=nclass)
    emobj <- shortemcluster(d, emobj)
    ret <- emcluster(d, emobj, assign.class = TRUE)
    ret$class
}

#' Given a 2 class sample, findTumourClass returns one of the indices
#' (1 or 2) depending on which has the higher overall keratin concentration
#'
#' @param sp An SPData object
#'
#' @export
findTumourID <- function(sp) {
    Y <- rawData(sp)

    keratin.index <- grep("Keratin", channels(sp))
    Yk <- Y[,keratin.index]

    if(length(keratin.index) == 0) stop("Keratin not found in channels!")
    if(length(keratin.index) > 1) stop("More than 1 keratin sample found!")

    cl1 <- which(cellClass(sp) == 1)
    cl2 <- which(cellClass(sp) == 2)

    meanCl1 <- mean(Yk[cl1])
    meanCl2 <- mean(Yk[cl2])


    which.max(c(meanCl1, meanCl2))
}
