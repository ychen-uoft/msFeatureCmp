# This file contains the public APIs for the msFeatureCmp package.

# The global pyOpenMS Python import. This allows the package to interface with
# the raw mass spectrometry data.
ropenms <- reticulate::import("pyopenms", convert = FALSE)

# Global constants for feature matrices
RT_IDX <- 1  # Retention time index
MZ_IDX <- 2  # Mass-to-charge index
IT_IDX <- 3  # Signal intensity index
# IM_IDX <- 4  # Ion mobility index

# Global constants for feature similarity
RT_THRESHOLD <- 5
MZ_THRESHOLD <- 0.01
# IM_THRESHOLD <- 0.031

#' Compares two sets of mass spectrometry features corresponding to a single MS
#' run.
#'
#' This can be used to compare the results of different feature finding
#' algorithms that have been run on the same raw dataset.
#'
#' TODO: what does this function do?
#'
#' @param featureFilePath1 The location of the first featureXML file, as a
#' string.
#' @param featureFilePath2 The location of the second featureXML file, as a
#' string.
#' @param rawDataFilePath The location of the mzML file, as a string.
#'
#' @examples
#' \dontrun{  # TODO: should this example be run?
#' compareFeatures("inst/extdata/featureSetA.featureXML",
#'                 "inst/extdata/featureSetB.featureXML",
#'                 "inst/extdata/20190122_HeLa_QC_Slot1-47_1_3228_800-860.mzML")
#' }
#'
#' @export
compareFeatures <- function(featureFilePath1, featureFilePath2,
                            rawDataFilePath) {
  # TODO: remove rawDataFilePath, unless there's some analysis that needs it
  # experiment <- ropenms$MSExperiment()
  # ropenms$MzMLFile()$load(rawDataFilePath, experiment)

  # Load both feature sets from their files
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
  numMultiplyMatchedFeaturesAB <- 0  # One-to-many matches (from A to B)
  numMultiplyMatchedFeaturesBA <- 0  # One-to-many matches (from B to A)

  # 1 - Find the intersection of the two feature sets
  # a) Use the first feature set as the reference set
  for (i in 1:featureSetA$size()) {
    featureA <- featureMatrixA[i, ]  # A vector of [RT, m/z, intensity]
    similarIdx <- findFirstFeature(featureMatrixB,
                                   target = featureA[RT_IDX] - RT_THRESHOLD)
    numSimilarFeatures <- 0

    # Prevent index out of range errors
    if (similarIdx == 0) {
      similarIdx <- 1
    }
    featureB <- featureMatrixB[similarIdx, ]

    while (withinThreshold(featureA[RT_IDX], featureB[RT_IDX], RT_THRESHOLD)) {
      if (similarFeatures(featureA, featureB)) {
        numSimilarFeatures <- numSimilarFeatures + 1
      }

      # Advance to the next feature in the second set while still in range
      similarIdx <- similarIdx + 1
      if (similarIdx > featureSetB$size()) {
        break
      }
    }

    # featureA was not similar to any feature in the second set
    if (numSimilarFeatures == 0) {
      numUnmatchedFeatures <- numUnmatchedFeatures + 1
    }
    # featureA was similar to a single feature in the second set
    else if (numSimilarFeatures == 1) {
      numSinglyMatchedFeatures <- numSinglyMatchedFeatures + 1
      numCommonFeatures <- numCommonFeatures + 1
    }
    # featureA was similar to multiple features in the second set
    else {
      numMultiplyMatchedFeaturesAB <- numMultiplyMatchedFeaturesAB + 1
      numCommonFeatures <- numCommonFeatures + 1
    }
  }

  # b) Use the second feature set as the reference set
  for (i in 1:featureSetB$size()) {
    featureB <- featureMatrixB[i, ]
    similarIdx <- findFirstFeature(featureMatrixA,
                                   target = featureB[RT_IDX] - RT_THRESHOLD)
    numSimilarFeatures <- 0

    if (similarIdx == 0) {
      similarIdx <- 1
    }
    featureA <- featureMatrixA[similarIdx, ]

    while (withinThreshold(featureB[RT_IDX], featureA[RT_IDX], RT_THRESHOLD)) {
      if (similarFeatures(featureB, featureA)) {
        numSimilarFeatures <- numSimilarFeatures + 1
      }

      similarIdx <- similarIdx + 1
      if (similarIdx > featureSetA$size()) {
        break
      }
    }

    # No adding to numCommonFeatures; no else statement to avoid overcounting
    if (numSimilarFeatures == 0) {
      numUnmatchedFeatures <- numUnmatchedFeatures + 1
    }
    else if (numSimilarFeatures > 1) {
      numMultiplyMatchedFeaturesBA <- numMultiplyMatchedFeaturesBA
    }
  }

  # Numbers have been collected; now do analysis with them
}

# TODO: add more functions for comparisons and plotting
