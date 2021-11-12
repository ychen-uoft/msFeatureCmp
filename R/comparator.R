# This file contains the public APIs for the msFeatureCmp package.

# The global pyOpenMS python import
ropenms <- reticulate::import("pyopenms", convert = FALSE)

# Other global constants
RT_IDX <- 1  # Retention time index
MZ_IDX <- 2  # Mass-to-charge index
IT_IDX <- 3  # Signal intensity index

RT_THRESHOLD <- 5  # For feature similarity
MZ_THRESHOLD <- 0.01  # For feature similarity

#' Compares two sets of mass spectrometry features corresponding to a single MS
#' run. For example, to compare the results of different feature finding
#' algorithms run on the same raw dataset.
#'
#' TODO: it should print a comparison, including statistical scores. What else?
#'
#' @param featureFilePath1 The location of the first featureXML file, as a string.
#' @param featureFilePath2 The location of the second featureXML file, as a string.
#' @param rawDataFilePath The location of the mzML file, as a string.
#'
#' @examples
#' \dontrun{
#' compareFeatures("inst/extdata/featureSetA.featureXML",
#'                 "inst/extdata/featureSetB.featureXML",
#'                 "inst/extdata/20190122_HeLa_QC_Slot1-47_1_3228_800-860.mzML")
#' }
#'
#' @export
compareFeatures <- function(featureFilePath1, featureFilePath2,
                            rawDataFilePath) {
  # Load the raw MS data first, followed by both feature sets
  experiment <- ropenms$MSExperiment()
  ropenms$MzMLFile()$load(rawDataFilePath, experiment)
  featureSetA <- ropenms$FeatureMap()
  ropenms$FeatureXMLFile()$load(featureFilePath1, featureSetA)
  featureSetB <- ropenms$FeatureMap()
  ropenms$FeatureXMLFile()$load(featureFilePath2, featureSetB)

  # Convert both feature sets to feature matrices for easier computation
  featureMatrixA <- convertFeaturesToSortedMatrix(featureSetA)
  featureMatrixB <- convertFeaturesToSortedMatrix(featureSetB)

  # Variables for statistical analysis
  numCommonFeatures <- 0  # A AND B size
  numUnmatchedFeatures <- 0  # A NAND B size
  numSinglyMatchedFeatures <- 0  # One-to-one matches
  numMultiplyMatchedFeaturesAB <- 0  # One-to-many matches (A to B)
  numMultiplyMatchedFeaturesBA <- 0  # One-to-many matches (B to A)

  # 1 - Find the intersection of the two sets
  # a) Use the first feature set as the reference set
  for (i in 1:featureSetA$size()) {
    featureA <- featureMatrixA[i, ]
    similarIdx <- findFirstFeature(featureMatrixB, target = featureA[RT_IDX])
    numSimilarFeatures <- 0

    featureB <- featureMatrixB[similarIdx, ]
    while (withinThreshold(featureA[RT_IDX], featureB[RT_IDX], RT_THRESHOLD)) {
      if (similarFeatures(featureA, featureB)) {
        numSimilarFeatures <- numSimilarFeatures + 1
      }
      similarIdx <- similarIdx + 1
    }

    if (numSimilarFeatures == 0) {
      numUnmatchedFeatures <- numUnmatchedFeatures + 1
    }
    else if (numSimilarFeatures == 1) {
      numSinglyMatchedFeatures <- numSinglyMatchedFeatures + 1
      numCommonFeatures <- numCommonFeatures + 1
    }
    else {
      numMultiplyMatchedFeaturesAB <- numMultiplyMatchedFeaturesAB + 1
      numCommonFeatures <- numCommonFeatures + 1
    }
  }

  # b) Use the second feature set as the reference set
  for (i in 1:featureSetB$size()) {
    featureB <- featureMatrixB[i, ]
    similarIdx <- findFirstFeature(featureMatrixA, target = featureB[RT_IDX])
    numSimilarFeatures <- 0

    featureA <- featureMatrixA[similarIdx, ]
    while (withinThreshold(featureB[RT_IDX], featureA[RT_IDX], RT_THRESHOLD)) {
      if (similarFeatures(featureB, featureA)) {
        numSimilarFeatures <- numSimilarFeatures + 1
      }
      similarIdx <- similarIdx + 1
    }

    # No adding to numCommonFeatures, no else statement to avoid overcounting
    if (numSimilarFeatures == 0) {
      numUnmatchedFeatures <- numUnmatchedFeatures + 1
    }
    else if (numSimilarFeatures > 1) {
      numMultiplyMatchedFeaturesBA <- numMultiplyMatchedFeaturesBA
    }
  }
}

# TODO: add more functions for comparisons and plotting
