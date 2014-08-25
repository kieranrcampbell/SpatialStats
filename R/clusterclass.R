########################################
## Tools for clustering by cell class ##
########################################


#' Clusters cells into different classes using PCA-Gaussian Mixture model clustering
#' 
#' A method of using high-dimensional proteomics data to classify cellular phenotypes (e.g. 
#' epithelial and stromal cells may be present within a tumour biopsy). First principle component
#' analysis is performed on the raw cell-by-channel matrix \code{Y}. Then the top \code{n.pc}
#' (default 3) principle components are used in Gaussian mixture modelling through 
#' expectation maximisation, which clusters the cellular readouts into \code{nclass} (default 2)
#' different classes. The underlying assumption is that in the dimensionality-reduced PC-space
#' the cells will be phenotypically distinct enough to appear drawn from two separate 
#' classes. 
#' 
#' This is an unsupervised method - other methods may be used (e.g. training a SVM
#' on known cell phenotypes). The cell classes can always be set by alternative methods using
#' \code{cellClass(sp) <- classes}
#' for an \code{SPData} \code{sp} and a vector of cell classes \code{classes}.
#'
#' @param Y A cell by channel matrix of readouts (e.g. from \code{cells(sp)})
#' @param doPCA If true, first performs PCA and uses the top
#' principle components for clustering
#' @param nclass The number of classes into which to cluster the cells (default = 2)
#' @param n.pc The number of principle components to use in Gaussian mixture modelling
#' (default = 3)
#' 
#' @references
#' Chen, W.-C., Maitra, R., Melnykov, V. (2012) EMCluster: EM Algorithm for Model-Based
#' Clustering of Finite Mixture Gaussian Distribution. R Package, URL
#' http://cran.r-project.org/package=EMCluster
#' 
#' @return A vector of length number of cells with numeric values corresponding
#' to distinct classes for each cell.
#'
#' @export
clusterClass <- function(Y, doPCA = TRUE, nclass=2, n.pc =3) {
    d <- NULL
    if(doPCA) {
        ypca <- princomp(Y)
        d <- ypca$scores[,1:n.pc]
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
