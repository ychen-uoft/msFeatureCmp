#' Loads a raw mass spectrometry run into memory. The MS file must be mzML
#' format.
#'
#' @param filePath The location of the mzML file.
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

#' Loads a mass spectrometry feature set into memory. The feature file must be
#' featureXML format.
#'
#' @param filePath The location of the featureXML file.
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

#' Compares two features to determine if they are identical.
#'
#' Features can be considered to be identical if their retention times and
#' mass-to-charges are within some thresholds of each other (which can be
#' provided by the user).
#'
#' @param feature1 The first feature to compare
#' @param feature2 The second feature to compare
#' @param rtThreshold The retention time threshold to use
#' @param mzThreshold The mass-to-charge threshold to use
#'
#' @return TRUE if the features are identical; FALSE otherwise.
#'
#' @examples
#' \dontrun{
#' featureA <- ropenms$Feature()
#' featureA$setRT(5)
#' featureA$setMZ(300)
#' featureB <- ropenms$Feature()
#' featureB$setRT(6)
#' featureB$setMZ(300.005)
#' identicalFeatures(featureA, featureB)  # TRUE
#' }
identicalFeatures <- function(feature1, feature2, rtThreshold = 5, mzThreshold = 0.01) {
  isIdentical <- FALSE
  if (abs(feature1$getRT() - feature2$getRT()) <= rtThreshold &
      abs(feature1$getMZ() - feature2$getMZ()) <= mzThreshold) {
    isIdentical <- TRUE
  }
  return(isIdentical)
}

#' Converts a set of features, as a FeatureMap, to a matrix of features.
#'
#' Each row of the matrix represents a feature, and contains its retention
#' time (RT), mass-to-charge (m/z), and signal intensity in that order. The
#' resulting matrix will be sorted by descending signal intensity.
#'
#' @param features The set of features to convert, as a FeatureMap.
#'
#' @return The corresponding sorted matrix for the input features.
#'
#' @examples
#' \dontrun{
#' read feature map
#' convert feature map
#' }
convertFeaturesToSortedMatrix <- function(features) {
  RT_IDX <- 1
  MZ_IDX <- 2
  IT_IDX <- 3

  featureMatrix <- matrix(nrow = 3, ncol = features$size())
  rowIter <- 1
  for (feature in features) {
    featureMatrix[rowIter, RT_IDX] <- feature$getRT()
    featureMatrix[rowIter, MZ_IDX] <- feature$getMZ()
    featureMatrix[rowIter, IT_IDX] <- feature$getIntensity()
    rowIter <- rowIter + 1
  }

  sortedFeatureMatrix <- featureMatrix[order(featureMatrix[,"V3", decreasing = TRUE]),]
  return(sortedFeatureMatrix)
}
