SpatialPRo - spatial proteomics analysis

### Installation

Using devtools call `install_github("kieranrcampbell/SpatialPRo")`

### Usage

`SpatialPRo` contains two different classes - `SPData` that represents a single (tissue) sample and `SPExp` that represents the entire experiment, storing a list of `SPData`s, the sample ids and the directory and filenames of the original files.

### SPData

There is an example dataset included - to load it first load the library via `library(SpatialPRo)` then call `data(sp5)`. This loads an example object of class `SPData` called `sp`. The following can then be used to access data elements:

#### Generic SPData methods
* `channels(sp)` returns the channel names
* `nChannel(sp)` returns the number of channels
* `nCells(sp)` returns the number of cells
* `rawData(sp)` returns a cell by channel matrix of the raw data
* `cells(sp)` returns a cell by channel matrix with the readouts for each cell, normalised wrt cell size & concentration
* `neighbours(sp)` returns a list of length nCells, where entry *i*  is an *n* by channel matrix for cell *i* having *n* neighbours
* `neighbourIDs(sp)` returns a list of the nearest neighbour cell identifiers
* `size(sp)` returns a vector containing the cell sizes (in pixels)
* `weight(sp)` returns the boundary sizes for each cell to its nearest neighbours
* `ID(sp)` returns the sample ID
* `cellClass(sp)` returns a vector of length nCells where each cell is assigned a class. Here there are 2 classes, with 1 corresponding to stromal and 2 corresponding to tumour.
* `xy(sp)` returns a nCell by 2 matrix of centre-of-mass locations for each cell
* An sp object can be partitioned using the [ operator, so `sp[1:2, c(3,5,7)]` returns a new SPData object
using cells 1 & 2 and channels 3, 5 and 7

#### Methods for dealing with different cell classes
* `cellClass(sp)` returns a vector of length nCells(sp) with a numeric classification for each cell (currently only 2 types are supported)
* `neighbourClass(sp, cell.class)` functions similarly to cellClass but returns a list similar to neighbours(sp) but with only neighbours of class cell.class remaining
* `neighbourChannel(sp, c(1,3,5)` returns the neighbour list as in neighbours(sp) but including only channels 1,3 and 5
* `findBoundary(sp)` returns the (numeric identifiers of) cells that lie on the boundary between 2 different classes


The class `SPData` also contains several plotting routines:

* `boxplots(sp,4,8)` plots the ranges of each of the 32 channels (in 4 rows by 8 columns)
* `channelPlot(sp, 1:4)` provides boxplots of the first 4 channels but highlighting differences between different classes of data (class1 vs class2)

### SPExp

An SPExp has four slots:

1. `dir` Directory location of experiment files
2. `files` List of files used in the experiment
3. `spdata` List containing objects of class `SPData`
4. `ids` Vector of sample ids

#### Generic SPExp Methods
* `files(spe)` shows the experiment files
* `getDir(spe)` experiment directory
* `IDs(spe)` sample IDs
* `spe[2:3]` subset samples 2-3
* `spe[[4]]` return an `SPData` object of the fourth sample
* `SPlist(spe)` the underlying list holding the `SPData` objects
* `length(spe)` number of samples

#### Example workflows
Download all the original matlab files from
https://s3.amazonaws.com/supplemental.cytobank.org/report_data/report_113/Figure_5/Figure_5_raw_image_files.zip
and extract to
~/myfolder

We can then initialise an empty `SPE` with
`spe <- SPExperimentfromDir("~/myfolder")`
If the `files` option is left blank then all files are used.

We can then load all the files into the object using
`spe <- loadExp(spe)`
(this may take some time for matlab -> R conversion)

