# Utility function unit tests

test_that("loading raw data works", {
  # Load the example mzML file and count the number of spectra
  experiment <- msFeatureCmp:::loadMSFile(system.file("extdata",
    "20190122_HeLa_QC_Slot1-47_1_3228_800-810.mzML", package = "msFeatureCmp"))
  numSpectra <- reticulate::py_to_r(experiment$getNrSpectra())
  expect_equal(numSpectra, 8)
})

test_that("loading feature data works", {
  # Load the first example featureXML file and count the number of features
  featureSet <- msFeatureCmp:::loadFeatureFile(
    system.file("extdata", "featureSetA.featureXML", package = "msFeatureCmp"))
  numFeatures <- reticulate::py_to_r(featureSet$size())
  expect_equal(numFeatures, 720)
})

test_that("threshold checking works", {
  # Simple distance checking with multiple thresholds
  x <- 5.0
  y <- 5.049
  t1 <- 0.05
  t2 <- 0.04

  expect_true(msFeatureCmp:::withinThreshold(x, y, t1))
  expect_false(msFeatureCmp:::withinThreshold(x, y, t2))

  # Try with random thresholds to check determinism
  x <- runif(1)
  y <- runif(1)
  thresholds <- runif(5)

  for (t in thresholds) {
    dist <- abs(x - y)
    expect_equal(msFeatureCmp:::withinThreshold(x, y, t), dist < t)
  }
})

test_that("similarity checking works", {
  ropenms <- reticulate::import("pyopenms", convert = FALSE)

  # Try checking Feature objects
  f1 = ropenms$Feature()
  f1$setRT(800)
  f1$setMZ(300)
  f2 = ropenms$Feature()
  f2$setRT(801)
  f2$setMZ(300.001)
  f3 = ropenms$Feature()
  f3$setRT(802)
  f3$setMZ(300.011)

  expect_true(msFeatureCmp:::similarFeatures(f1, f2))
  expect_false(msFeatureCmp:::similarFeatures(f1, f3))

  # Try checking vectors
  f1 <- c(800, 300)
  f2 <- c(801, 300.001)
  f3 <- c(802, 300.011)

  expect_true(msFeatureCmp:::similarFeatures(f1, f2))
  expect_false(msFeatureCmp:::similarFeatures(f1, f3))
})

test_that("matrix sorting works", {
  # Sorting a small matrix by column number
  m <- matrix(1:9, 3, 3)
  m[1, 1] <- 2
  m[2, 1] <- 1
  m[2, 2] <- 6
  m[3, 2] <- 5
  # [[ 2 4 7 ]
  #  [ 1 6 8 ]
  #  [ 3 5 9 ]]

  # First column as key, ascending
  sortedMatrix <- msFeatureCmp:::sortMatrixByColumn(m, 1L, FALSE)
  expect_equal(sortedMatrix[1, ], c(1, 6, 8))
  expect_equal(sortedMatrix[2, ], c(2, 4, 7))
  expect_equal(sortedMatrix[3, ], c(3, 5, 9))

  # Third column as key, descending
  sortedMatrix <- msFeatureCmp:::sortMatrixByColumn(m, 3L, TRUE)
  expect_equal(sortedMatrix[1, ], c(3, 5, 9))
  expect_equal(sortedMatrix[2, ], c(1, 6, 8))
  expect_equal(sortedMatrix[3, ], c(2, 4, 7))
})

test_that("feature conversion works", {
  # Convert the first example featureXML file into a matrix
  featureSet <- msFeatureCmp:::loadFeatureFile(
    system.file("extdata", "featureSetA.featureXML", package = "msFeatureCmp"))
  featureMatrix <- msFeatureCmp:::convertFeaturesToSortedMatrix(featureSet)

  expect_equal(reticulate::py_to_r(featureSet$size()), nrow(featureMatrix))
  expect_equal(featureMatrix[1, ],
               c(800.002794376452, 614.596298389786, 10900.7001953125))
  expect_equal(featureMatrix[2, ],
               c(800.008514955597, 501.302382748776, 19653.0))
  expect_equal(featureMatrix[nrow(featureMatrix), ],
               c(809.9793354586872, 564.20839028551, 5604.33984375))
})

test_that("binary search works", {
  # Binary searching a small matrix
  m <- matrix(1:10, 5, 2)
  m[2, 1] <- 3
  m[4, 1] <- 3
  # [[1, 6], [3, 7], [3, 8], [3, 9], [5, 10]]

  # First occurrence of 3 is at index 2
  expect_equal(msFeatureCmp:::findFirstFeature(m, 1, 3), 2)
  # First occurrence of 1 is at index 1
  expect_equal(msFeatureCmp:::findFirstFeature(m, 1, 1), 1)
  # First occurrence of 5 is at index 5
  expect_equal(msFeatureCmp:::findFirstFeature(m, 1, 5), 5)
  # No occurrences of 0, so returns the rank of 0 in the matrix, which is 0
  expect_equal(msFeatureCmp:::findFirstFeature(m, 1, 0), 0)
  # No occurrences of 6, so returns the rank of 6 in the matrix, which is 5
  expect_equal(msFeatureCmp:::findFirstFeature(m, 1, 6), 5)
})
