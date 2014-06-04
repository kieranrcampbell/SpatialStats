## creates small mock dataset


test_that("SPData is initialised properly", {
    ncells <- 5
    nprot <- 10

    Y <- t(replicate(5, rpois(10,5), simplify=T))

    protein.names <- paste("prot", 1:nprot, sep="")

    X <- lapply(1:ncells, function(i) {
        Y[sample(1:ncells, sample(2:4,1)),]

    })

    sp <- SPData(protein.names=protein.names,
                 Y=Y,X=X)

    expect_that(nCells(sp), equals(ncells))
    expect_that(nProt(sp), equals(nprot))
    expect_that(length(pNames(sp)), equals(nprot))
})

