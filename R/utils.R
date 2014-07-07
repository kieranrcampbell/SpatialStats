
######################################
## Utility functions for SpatialPRo ##
######################################

#' This takes the entire set of cells and returns sets whose
#' members aren't neighbours of each other
findDisjointCellSet <- function(sp) {
    nc <- nCells(sp)
    nn <- neighbourIDs(sp)

    x <- allCells <- 1:nc
    listofSets <- list()
    counter <- 1

    while(length(x) > 0) {
        disSet <- removeNeighbours(nn, x)
        listofSets[[ counter ]] <- disSet
        counter <- counter + 1 ;
        x <- setdiff(allCells, unlist(listofSets))
    }

    return(listofSets)
}

removeNeighbours <- function(nn, x) {
    xcopy <- x
    for(i in xcopy) {
        if(i %in% x) {
            x <- setdiff(x, nn[[i]])
        }
    }
    x
}
