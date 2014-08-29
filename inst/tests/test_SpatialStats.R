## Test functions in spatialstats

library(SpatialPRo)

test_that("Sample binding is working properly", {

  ## tumour classes from looking are 2 1 2 1 2
  visual.class <- c(2,1,2,1,2)
  tumour.classes <- FindIDs(SPE, "Vimentin", "lower", "median") 
  
  expect_equal(visual.class, tumour.classes)
  
  bound <- BindSPE(SPE, choose.class = tumour.classes)
  
  expect_equal(ncol(bound$Y), nChannel(SPE[[1]]))
  expect_equal(sum(bound$sizes), nrow(bound$Y))
    
  sample.factors <- ConstructSampleFactors(bound, IDs(SPE),intercept = TRUE)
  
  f <- sapply(1:length(SPE), function(i) rep(ID(SPE[[i]]), bound$sizes[i]))
  f <- as.factor(unlist(f))
  m <- model.matrix(~ f)
  m <- m[,-1]
  colnames(m) <- IDs(SPE)[2:5]  

  expect_equal(sample.factors, m)
})

test_that("FindOverlap", {
  m1 <- matrix(c(1,1,2,2,2,3,3,4),ncol=2,byrow=TRUE)
  m2 <- matrix(c(1,1,3,4), ncol=2, byrow=TRUE)
  results <- list(m1, m2)
  olap.diff <- FindOverlap(results, "different")  
  olap.same <- FindOverlap(results, "same")
  
  expect_equal(olap.diff, c(3,4))
  expect_equal(olap.same, c(1,1))
})
