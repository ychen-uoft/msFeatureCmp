# This file contains the private helper functions for the msFeatureCmp package.
# Examples in this file cannot be run by roxygen because all the functions are
# private.

#' Loads a raw mass spectrometry run into memory as an MSExperiment object. The
#' input MS file must be in mzML (OpenMS) format.
#'
#' @param filePath The location of the mzML file, as a string.
#'
#' @return The in-memory representation of the MS run, as an MSExperiment.
#'
#' @examples
#' \dontrun{
#' experiment <- loadMSFile("inst/extdata/20190122_HeLa_QC_Slot1-47_1_3228_800-860.mzML")
#' experiment$getNrSpectra()  # Returns 47
#' }
loadMSFile <- function(filePath) {
  experiment <- ropenms$MSExperiment()
  ropenms$MzMLFile()$load(filePath, experiment)
  return(experiment)
}

#' Loads a mass spectrometry feature set (containing found features) into
#' memory as a FeatureMap object. The feature file must be in featureXML
#' (OpenMS) format.
#'
#' @param filePath The location of the featureXML file, as a string.
#'
#' @return The in-memory representation of the feature set, as a FeatureMap.
#'
#' @examples
#' \dontrun{
#' featureSet <- loadFeatureFile("inst/extdata/featureSetA.featureXML")
#' featureSet$size()  # Returns 3770
#' }
loadFeatureFile <- function(filePath) {
  featureSet <- ropenms$FeatureMap()
  ropenms$FeatureXMLFile()$load(filePath, featureSet)
  return(featureSet)
}

#' Checks if two values are within a certain distance of each other.
#'
#' @param value1 The first value, as a number.
#' @param value2 The second value, as a number.
#' @param threshold The maximum distance between the values, as a number.
#'
#' @return TRUE if the values are within the distance; FALSE otherwise.
#'
#' @examples
#' \dontrun{
#' x1 <- 5.1
#' x2 <- 5.05
#' maxDist <- 0.01
#' withinThreshold(x1, x2, maxDist)  # Returns FALSE
#' }
withinThreshold <- function(value1, value2, threshold) {
  isWithinThreshold <- FALSE
  if (abs(value1 - value2) <= threshold) {
    isWithinThreshold <- TRUE
  }
  return(isWithinThreshold)
}

#' Compares two features to determine if they are similar.
#'
#' We consider features to be similar if their retention times and mass-to-
#' charge values are within some fixed thresholds of each other (which can be
#' provided by the user).
#'
#' The provided features can be ropenms Feature objects or vectors (e.g. from
#' a feature matrix).
#'
#' @param feature1 The first feature to compare.
#' @param feature2 The second feature to compare.
#' @param rtThreshold The RT threshold to use, as a number.
#' @param mzThreshold The m/z threshold to use, as a number.
#'
#' @return TRUE if the two features are similar; FALSE otherwise.
#'
#' @examples
#' \dontrun{
#' featureA <- ropenms$Feature()
#' featureA$setRT(5)
#' featureA$setMZ(300)
#' featureB <- ropenms$Feature()
#' featureB$setRT(6)
#' featureB$setMZ(300.005)
#' similarFeatures(featureA, featureB)  # Returns TRUE
#' }
similarFeatures <- function(feature1, feature2, rtThreshold = RT_THRESHOLD,
                            mzThreshold = MZ_THRESHOLD) {
  isSimilar <- FALSE
  # "Function overloading" - if the features are rows in a feature matrix
  if (is.vector(feature1) & is.vector(feature2) &
      withinThreshold(feature1[RT_IDX], feature2[RT_IDX], rtThreshold) &
      withinThreshold(feature1[MZ_IDX], feature2[MZ_IDX], mzThreshold)) {
    isSimilar <- TRUE
  }
  # If the features are Features in a feature set
  else if (withinThreshold(feature1$getRT(), feature2$getRT(), rtThreshold) &
           withinThreshold(feature1$getMZ(), feature2$getMZ(), mzThreshold)) {
    isSimilar <- TRUE
  }
  # TODO: maybe need to enforce a type check?
  return(isSimilar)
}

#' Sorts a matrix by a given column number.
#'
#' @param matrix The matrix to sort.
#' @param column The column number to use as key.
#' @param descending If TRUE, sorts the matrix by descending values.
#'
#' @return The sorted matrix.
#'
#' @examples
#' \dontrun{
#' m <- matrix(1:9, nrow = 3, ncol = 3)
#' m[1, 1] <- 2
#' m[2, 1] <- 1
#' m[2, 2] <- 6
#' m[3, 2] <- 5
#' sortedMatrix <- sortMatrixByColumn(m, descending = TRUE)
#' sortedMatrix
#' }
sortMatrixByColumn <- function(matrix, column = 1, descending = FALSE) {
  # Need to enforce that the column is an integer
  if (typeof(column) != "integer") {
    cat("msFeatureCmp::sortMatrixByColumn must take an int as column number\n")
  }
  key <- paste("V", column, sep = "")
  sortedMatrix <- matrix[order(matrix[ , key, descending]), ]
  return(sortedMatrix)
}

#' Converts a loaded set of features into a matrix of features.
#'
#' Each row of the matrix represents a feature, and contains its retention
#' time, mass-to-charge, and signal intensity values, in that order. The matrix
#' will be also sorted by ascending RT.
#'
#' @param featureSet The set of features to convert, as a FeatureMap.
#'
#' @return The corresponding sorted matrix for the feature set.
#'
#' @examples
#' \dontrun{
#' featureSet <- loadFeatureFile("inst/extdata/featureSetA.featureXML")
#' featureMatrix <- convertFeaturesToSortedMatrix(featureSet)
#' featureMatrix[1, ]  # Returns [802.979781716543, 861.395085396801, 280332.0]
#' }
convertFeaturesToSortedMatrix <- function(featureSet) {
  # Start with an empty matrix of the correct size, and fill it one feature
  # at a time
  featureMatrix <- matrix(nrow = featureSet$size(), ncol = 3)
  rowNum <- 1
  for (feature in featureSet) {
    # Global constants are defined in comparator.R
    featureMatrix[rowNum, RT_IDX] <- feature$getRT()
    featureMatrix[rowNum, MZ_IDX] <- feature$getMZ()
    featureMatrix[rowNum, IT_IDX] <- feature$getIntensity()
    rowNum <- rowNum + 1
  }

  # Sort by the first column (RT), ascending
  sortedFeatureMatrix <- sortMatrixByColumn(featureMatrix, RT_IDX, FALSE)
  return(sortedFeatureMatrix)
}

#' Finds the first feature in a feature matrix that has its element at a given
#' key equal to a given target. If no feature satisfies the requirement, counts
#' how many features have its element at that same key less than the target
#' (its rank in the matrix).
#'
#' This function uses binary search to find the feature, so the feature matrix
#' must be sorted (ascending, in the key column).
#'
#' @param sortedFeatureMatrix The sorted feature matrix to search.
#' @param key The column in the matrix to search in.
#' @param target The target value to search for, as a number.
#'
#' @return The index of the first feature that has its element at the given key
#' equal to the given target, if it exists; the feature's rank otherwise.
#'
#' @examples
#' \dontrun{
#' m <- matrix(1:10, nrow = 5, ncol = 2)
#' m[2, 1] <- 3
#' m[4, 1] <- 3
#' rank <- findFirstFeature(m, key = 1, target = 3)
#' rank  # Returns 2
#' }
findFirstFeature <- function(sortedFeatureMatrix, key = RT_IDX, target = 0) {
  firstIdx <- 1
  lastIdx <- nrow(sortedFeatureMatrix)

  # Repeatedly halve the search space
  while (firstIdx < lastIdx) {
    middleIdx <- trunc((firstIdx + lastIdx) / 2)
    if (sortedFeatureMatrix[middleIdx, key] < target) {
      # Repeat the search on the "lower" half of the matrix
      firstIdx <- middleIdx + 1
    }
    else {
      # Repeat the search on the "upper" half of the matrix
      lastIdx <- middleIdx
    }
  }

  # Error checking is left to the caller
  resultIdx <- 0
  if (firstIdx > 1) {
    resultIdx <- firstIdx - 1
  }
  return(resultIdx)
}
