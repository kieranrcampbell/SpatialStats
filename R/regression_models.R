## regression models for nearest neighbour data
## kieranrcampbell@gmail.com

#' Nearest neighbour regression
#' @export
setMethod("nnReg", "SPData", function(object) {
    Y <- cells(object)

    X <- t(sapply(object@X, function(x) {
        if(!is.matrix(x)) x
        else colMeans(x)
    }))

    fit <- lm(Y ~ X)
    fit
})
