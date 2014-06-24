library(EMCluster)
library(devtools)
library(rgl)
library(R.matlab)

load_all("..")

load("~/ebi/data/spe.RData")

sp <- spe[[5]]

Y <- cells(sp)

Y <- Y - min(Y) + 1

ycor <- cor(t(Y))

y.pc <- princomp(Y)

y.scores <- y.pc$scores


d <- y.scores[,1:3]

emobj <- simple.init(d, nclass=2)
emobj <- shortemcluster(d, emobj)
ret <- emcluster(d, emobj, assign.class=TRUE)

nnids <- nnID(sp)

clust <- ret$class

cl1 <- which(clust == 1)
cl2 <- which(clust == 2)


## find out which is tumour and which is stromal
kerindex <- grep("Keratin",pNames(sp))
kerReads <- Y[,kerindex]
tumourID <- which.max( c(mean(kerReads[cl1]), mean(kerReads[cl2])))

boundary <- NULL

for(i in cl1) {
    nn <- nnids[[i]]
    if(length(intersect(nn,cl2)) > 0) boundary <- c(boundary, i)
}

for(i in cl2) {
    nn <- nnids[[i]]
    if(length(intersect(nn, cl1)) > 0) boundary <- c(boundary, i)
}

protein.names <- pNames(sp)
getProteinIds <- function(proteinNames,proteinList) {
    p.indices <- sapply(proteinList, grep, proteinNames)
    p.indices
}


responseSubset <- getProteinIds(protein.names,
                                c("Cadherin",
                                  "bcat",
                                  "Vimentin",
                                  "CD44","Twist",
                                  "Slug", "NFkB",
                                  "EGFR"))


dependentRemoveSubset <- getProteinIds(protein.names,
                                       c("Keratin",
                                         "ER", "PR",
                                         "Ki67","CD"))

dependentRemoveSubset <- unlist(dependentRemoveSubset)

pset <- 1:length(protein.names)

pset <- setdiff(pset, dependentRemoveSubset)

cell.class <- ret$class

weights <- readMat("../matlab/boundary_sizes.mat")
weights <- weights[[1]]
weights <- lapply(weights, as.vector)

fit <- weightedSubsetBoundaryRegression(sp, cell.class, 2, boundary,
                                        weights, responseSubset, pset)

print("Reponse variables")
print(protein.names[responseSubset])
print("Dependent variables")
print(protein.names[pset])
