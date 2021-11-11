# This file contains the helper functions for the msFeatureCmp package.

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
#' }
loadMSFile <- function(filePath) {
  experiment <- ropenms$MSExperiment()
  ropenms$MzMLFile()$load(filePath, experiment)
  return(experiment)
}

#' Loads a mass spectrometry feature set (containing found features) into
#' memory as a FeatureMap. The feature file must be in featureXML (OpenMS)
#' format.
#'
#' @param filePath The location of the featureXML file, as a string.
#'
#' @return The in-memory representation of the feature set, as a FeatureMap.
#'
#' @examples
#' \dontrun{
#' featureSet <- loadFeatureFile("inst/extdata/featureSetA.featureXML")
#' }
loadFeatureFile <- function(filePath) {
  featureSet <- ropenms$FeatureMap()
  ropenms$FeatureXMLFile()$load(filePath, featureSet)
  return(featureSet)
}

#' Compares two features to determine if they are similar.
#'
#' We consider features to be similar if their retention times and mass-to
#' charge values are within some fixed thresholds of each other (which can be
#' provided by the user).
#'
#' The default thresholds have been experimentally determined to provide the
#' "best" results.
#'
#' @param feature1 The first feature to compare, as a Feature.
#' @param feature2 The second feature to compare, as a Feature.
#' @param rtThreshold The retention time threshold to use, as a float.
#' @param mzThreshold The mass-to-charge threshold to use, as a float.
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
similarFeatures <- function(feature1, feature2, rtThreshold = 5,
                              mzThreshold = 0.01) {
  isSimilar <- FALSE
  if (abs(feature1$getRT() - feature2$getRT()) <= rtThreshold &
      abs(feature1$getMZ() - feature2$getMZ()) <= mzThreshold) {
    isSimilar <- TRUE
  }
  return(isSimilar)
}

#' Converts a set of features, as a FeatureMap, into a matrix of features.
#'
#' Each row of the matrix represents a feature, and contains its retention
#' time, mass-to-charge, and signal intensity values, in that order. The
#' resulting matrix will be also sorted by descending signal intensity.
#'
#' @param featureSet The set of features to convert, as a FeatureMap.
#'
#' @return The corresponding sorted matrix for the features.
#'
#' @examples
#' \dontrun{
#' featureSet <- loadFeatureFile("inst/extdata/featureSetA.featureXML")
#' featureMatrix <- convertFeaturesToSortedMatrix(featureSet)
#' featureMatrix[1]  # TODO: what does this return?
#' }
convertFeaturesToSortedMatrix <- function(featureSet) {
  # Start with an empty matrix of the correct size, and fill it one feature
  # at a time.
  featureMatrix <- matrix(nrow = 3, ncol = featureSet$size())
  rowIter <- 1
  for (feature in featureSet) {
    # Global constants are defined in comparator.R
    featureMatrix[rowIter, RT_IDX] <- feature$getRT()
    featureMatrix[rowIter, MZ_IDX] <- feature$getMZ()
    featureMatrix[rowIter, IT_IDX] <- feature$getIntensity()
    rowIter <- rowIter + 1
  }

  # Sort by the third column (signal intensity), descending
  sortedFeatureMatrix <- featureMatrix[order(featureMatrix[ , "V3",
                                             decreasing = TRUE]), ]
  return(sortedFeatureMatrix)
}

sortMatrixByColumn <- function(matrix, column = 1, descending = FALSE) {
  # Need to enforce that the column is an integer
  if (typeof(column) != "numeric") {
    cat("msFeatureCmp::sortMatrixByColumn must take a numeric column number\n")
  }
  key <- cat("V", column)
  sortedMatrix <- matrix[order(matrix[ , key, descending]), ]
  return(sortedMatrix)
}

findFirstFeature <- function(sortedFeatureMatrix, key = RT_IDX, target = 0) {
  firstIdx <- 1
  lastIdx <- nrow(sortedFeatureMatrix)

  # Continuously halve the search space
  while (firstIdx < lastIdx) {
    middleIdx <- trunc((firstIdx + lastIdx) / 2)
    if (sortedFeatureMatrix[middleIdx, key] < target) {
      # Repeat the search on the lower half of the matrix
      firstIdx <- middleIdx + 1
    }
    else {
      # Repeat the search on the upper half of the matrix
      lastIdx <- middleIdx
    }
  }

  resultIdx <- 0
  if (firstIdx >= 0) {
    # The target was found, so update the result's index
    resultIdx <- firstIdx - 1
  }
  return(resultIdx)
}

findFirstFeatureByRT <- function(featureMatrix, target) {
  sortedFeatureMatrix <- sortMatrixByColumn(featureMatrix, RT_IDX)
  targetIdx <- findFirstFeature(sortedFeatureMatrix, RT_IDX, target)
  return(sortedFeatureMatrix[targetIdx])
}

findFirstFeatureByMZ <- function(featureMatrix, target) {
  sortedFeatureMatrix <- sortMatrixByColumn(featureMatrix, MZ_IDX)
  targetIdx <- findFirstFeature(sortedFeatureMatrix, MZ_IDX, target)
  return(sortedFeatureMatrix[targetIdx])
}
