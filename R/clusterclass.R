
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
