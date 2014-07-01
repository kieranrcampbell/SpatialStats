SpatialPRo - spatial proteomics analysis

### Installation

Using devtools call `install_github("kieranrcampbell/SpatialPRo")`

### Usage

There is an example dataset included - to load it first load the library via `library(SpatialPRo)` then call `data(sp5)`. This loads an example object of class `SPData` called `sp`. The following can then be used to access data elements:

* `channels(sp)` returns the channel names
* `nChannel(sp)` returns the number of channels
* `nCells(sp)` returns the number of cells
* `cells(sp)` returns a cell by channel matrix with the readouts for each cell
* `neighbours(sp)` returns a list of length nCells, where entry *i*  is an *n* by channel matrix for cell *i* having *n* neighbours
* `neighbourIDs(sp)` returns a list of the nearest neighbour cell identifiers
* `size(sp)` returns a vector containing the cell sizes (in pixels)
* `weights(sp)` returns the boundary sizes for each cell to its nearest neighbours
* `id(sp)` returns the sample ID
* `cellClass(sp)` returns a vector of length nCells where each cell is assigned a class. Here there are 2 classes, with 1 corresponding to stromal and 2 corresponding to tumour.