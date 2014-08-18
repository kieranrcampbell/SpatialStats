
data(SPE_report)


ncells <- 5
nprot <- 10
sID <- 1

Y <- t(replicate(5, rpois(10,5), simplify=T))

protein.names <- paste("prot", 1:nprot, sep="")

X <- lapply(1:ncells, function(i) {
    Y[sample(1:ncells, sample(2:4,1)),]

})

pos <- matrix(sample(1:10, ncells * 2, TRUE), ncol=2)

sp <- SPData(channelNames=protein.names,
             readouts=Y,raw=Y, cellNeighbours=X,
             size=rep(1, ncells), id=sID, weights=list(0),
             pos=pos, cellClass=-1)

spe <- SPExperiment("dir","files",list(sp,sp,sp))

test_that("SPExp is initialised properly", {
    expect_equal(files(spe), "files")
    expect_equal(getDir(spe), "dir")
    expect_equal(IDs(spe), rep(1,3))

    expect_true(all.equal(list(sp,sp,sp), SPlist(spe)))

    expect_true(all.equal(cells(sp), cells(spe[[1]])))

    expect_equal(length(spe), 3)

    expect_equal(3, getIDfromTMAname("TMA_ID3_HI"))
})

test_that("Manipulation functions used in regression work", {
    XY <- BindSPE(spe, FALSE, FALSE, FALSE)

    cell <- cells(spe[[1]])
    y <- rbind(cell, cell, cell)
    expect_equal(y, XY$Y)

    nm <- neighbourMean(spe[[1]], FALSE, FALSE)
    x <- rbind(nm, nm, nm)
    expect_equal(x, XY$X)

    factors <- ConstructSampleFactors(XY, 1:3)

    expect_equal(ncol(factors), 2)
    expect_equal(sum(factors), 2*nCells(spe[[1]]))

})
