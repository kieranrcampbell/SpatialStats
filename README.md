SpatialPRo - spatial proteomics analysis

### Usage

For an example dataset call `load('data/sp5.RData')`. This loads an example object of class `SPData` called `sp`. The following can then be used to access data elements:

* `channels(sp)` returns the channel names
* `nChannel(sp)` returns the number of channels
* `nCells(sp)` returns the number of cells
* `cells(sp)` returns a cell by channel matrix with the readouts for each cell
* `neighbours(sp)` returns a list of length nCells, where entry *i*  is an *n* by channel matrix for cell *i* having *n* neighbours 
* `neighbourIDs(sp)` returns a list of the nearest neighbour cell identifiers
* `size(sp)` returns a vector containing the cell sizes (in pixels)
* `weights(sp)` returns the boundary sizes for each cell to its nearest neighbours
* `id(sp)` returns the sample ID