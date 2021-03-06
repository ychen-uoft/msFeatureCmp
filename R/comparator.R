# comparator.R
# Package: msFeatureCmp
# Author: Yijia Chen
# Date: 2021-12-08
# Version: 0.1.0

# This file contains all the public APIs for the msFeatureCmp package.

# Global constants for feature matrices
RT_IDX <- 1  # Retention time index
MZ_IDX <- 2  # Mass-to-charge index
IT_IDX <- 3  # Signal intensity index
# IM_IDX <- 4  # Ion mobility index

# Global constants for feature similarity
RT_THRESHOLD <- 5
MZ_THRESHOLD <- 0.01
# IM_THRESHOLD <- 0.031

#' Feature set comparator.
#'
#' Compares two sets of mass spectrometry features corresponding to a single MS
#' run. This can be used to compare the results of different feature finding
#' algorithms that have been run on the same raw dataset.
#'
#' A short summary, including the number of features that were found to be
#' common to both sets of features, as well as statistical values such as
#' recall, precision, and F_1 score, will be printed to the terminal.
#'
#' This function may take a while to run.
#'
#' @param rawDataFilePath The location of the mzML raw data file, as a string.
#' @param featureFilePath1 The location of the first featureXML file, as a
#' string.
#' @param featureFilePath2 The location of the second featureXML file, as a
#' string.
#'
#' @return Nothing. A basic analysis is printed to the screen instead.
#'
#' @examples
#' \dontrun{
#' compareFeatures("inst/extdata/20190122_HeLa_QC_Slot1-47_1_3228_800-810.mzML",
#'                 "inst/extdata/featureSetA.featureXML",
#'                 "inst/extdata/featureSetB.featureXML")
#' }
#'
#' @export
#' @import reticulate
compareFeatures <- function(rawDataFilePath, featureFilePath1,
                            featureFilePath2) {
  ropenms <- reticulate::import("pyopenms", convert = FALSE)
  # browser()

  # Load the raw data and both feature sets
  experiment <- loadMSFile(rawDataFilePath)
  featureSetA <- loadFeatureFile(featureFilePath1)
  featureSetB <- loadFeatureFile(featureFilePath2)

  # Convert both feature sets to feature matrices for easier computation
  featureMatrixA <- convertFeaturesToSortedMatrix(featureSetA)
  featureMatrixB <- convertFeaturesToSortedMatrix(featureSetB)

  # Variables for statistical analysis
  numCommonFeatures <- 0  # In A and B
  numUnmatchedFeaturesA <- 0  # In A but not in B
  numUnmatchedFeaturesB <- 0  # In B but not in A
  numSinglyMatchedFeatures <- 0  # One-to-one matches
  numMultiplyMatchedFeaturesAB <- 0  # One-to-many matches (from A to B)
  numMultiplyMatchedFeaturesBA <- 0  # One-to-many matches (from B to A)

  # 1 - Find the intersection of the two feature sets
  # a) Use the first feature set as the reference set
  for (i in 1:reticulate::py_to_r(featureSetA$size())) {
    featureA <- featureMatrixA[i, ]  # A vector of [RT, m/z, intensity]
    similarIdx <- findFirstFeature(featureMatrixB,
                                   target = featureA[RT_IDX] - RT_THRESHOLD)
    numSimilarFeatures <- 0

    # Prevent index out of range errors
    if (similarIdx == 0) {
      similarIdx <- 1
    }
    featureB <- featureMatrixB[similarIdx, ]

    # Iterate through the features in the second set from the lower bound
    while (TRUE) {
      if (similarFeatures(featureA, featureB)) {
        numSimilarFeatures <- numSimilarFeatures + 1
      }

      similarIdx <- similarIdx + 1
      if (similarIdx > reticulate::py_to_r(featureSetB$size())) {
        # Reached the end of the second set
        break
      }
      featureB <- featureMatrixB[similarIdx, ]
      if (featureB[RT_IDX] > featureA[RT_IDX] + RT_THRESHOLD) {
        # Reached the end of the RT threshold window
        break;
      }
    }

    # featureA was not similar to any feature in the second set
    if (numSimilarFeatures == 0) {
      numUnmatchedFeaturesA <- numUnmatchedFeaturesA + 1
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

  # b) Use the second feature set as the reference set. The process is similar
  # to using the first feature set as the reference set.
  for (i in 1:reticulate::py_to_r(featureSetB$size())) {
    featureB <- featureMatrixB[i, ]
    similarIdx <- findFirstFeature(featureMatrixA,
                                   target = featureB[RT_IDX] - RT_THRESHOLD)
    numSimilarFeatures <- 0

    if (similarIdx == 0) {
      similarIdx <- 1
    }
    featureA <- featureMatrixA[similarIdx, ]

    # It is provable that this loop (and the one above) will never result in an
    # infinite loop
    while (TRUE) {
      if (similarFeatures(featureB, featureA)) {
        numSimilarFeatures <- numSimilarFeatures + 1
      }

      similarIdx <- similarIdx + 1
      if (similarIdx > reticulate::py_to_r(featureSetA$size())) {
        break
      }
      featureA <- featureMatrixA[similarIdx, ]
      if (featureA[RT_IDX] > featureB[RT_IDX] + RT_THRESHOLD) {
        break
      }
    }

    # No adding to numCommonFeatures; no else statement to avoid overcounting
    if (numSimilarFeatures == 0) {
      numUnmatchedFeaturesB <- numUnmatchedFeaturesB + 1
    }
    else if (numSimilarFeatures > 1) {
      numMultiplyMatchedFeaturesBA <- numMultiplyMatchedFeaturesBA
    }
  }

  # 2 - Number crunching
  # This output can be captured into a character vector by wrapping a call with
  # capture.output
  cat("\nmsFeatureCmp results\n")
  cat("Raw data file :", rawDataFilePath, "\n")
  cat("Feature file 1:", featureFilePath1, "\n")
  cat("Feature file 2:", featureFilePath2, "\n\n")

  numSpectra <- reticulate::py_to_r(experiment$getNrSpectra())
  numPeaks <- 0

  for (i in 0:(numSpectra - 1)) {
    spectrum <- experiment$getSpectrum(i)
    numPeaks <- numPeaks + reticulate::py_to_r(spectrum$size())
  }

  cat("Number of spectra:", numSpectra, "\n")
  cat("Number of peaks:  ", numPeaks, "\n")
  cat("Features in set 1:", reticulate::py_to_r(featureSetA$size()), "\n")
  cat("Features in set 2:", reticulate::py_to_r(featureSetB$size()), "\n\n")

  # a) Direct comparison between both feature sets
  totalNumFeatures <- reticulate::py_to_r(featureSetA$size()) +
    reticulate::py_to_r(featureSetB$size())
  totalUnmatchedFeatures <- numUnmatchedFeaturesA + numUnmatchedFeaturesB
  numMultiplyMatchedFeatures <- numMultiplyMatchedFeaturesAB +
    numMultiplyMatchedFeaturesBA

  percentCommon <- (numCommonFeatures / totalNumFeatures) * 100
  percentUnmatched <- (totalUnmatchedFeatures / totalNumFeatures) * 100
  percentSinglyMatched <- (numSinglyMatchedFeatures / totalNumFeatures) * 100
  percentMultiplyMatched <- (numMultiplyMatchedFeatures /
                               totalNumFeatures) * 100

  cat("Total features (both sets):", totalNumFeatures, "\n")
  cat("Number of common features: ", numCommonFeatures, "\n")
  cat("Percent common:             ", percentCommon, "%\n", sep = "")
  cat("Number of zero matches:    ", totalUnmatchedFeatures, "\n")
  cat("Percent unmatched:          ", percentUnmatched, "%\n", sep = "")
  cat("Number of single matches:  ", numSinglyMatchedFeatures, "\n")
  cat("Percent single:             ", percentSinglyMatched, "%\n", sep = "")
  cat("Number of multiple matches:", numMultiplyMatchedFeatures, "\n")
  cat("Percent multiple:           ",
      percentMultiplyMatched, "%\n\n", sep = "")

  # b) Use the first feature set as the ground truth
  numCommonFeatures <- numSinglyMatchedFeatures + numMultiplyMatchedFeaturesBA
  recall <- numCommonFeatures / reticulate::py_to_r(featureSetA$size())
  precision <- numCommonFeatures / reticulate::py_to_r(featureSetB$size())
  f1Score <- (2 * precision * recall) / (precision + recall)

  percentUnmatched <- (numUnmatchedFeaturesA /
                         reticulate::py_to_r(featureSetA$size())) * 100
  percentSinglyMatched <- (numSinglyMatchedFeatures /
                             reticulate::py_to_r(featureSetA$size())) * 100
  percentMultiplyMatched <- (numMultiplyMatchedFeaturesAB /
                               reticulate::py_to_r(featureSetA$size())) * 100

  # Prints some common statistics. I chose to not make this a global helper
  # function because it is more like a macro or an inline function that will
  # only ever be used in this function.
  # zeroMatches: The number of unmatched features, as an integer
  # multipleMatches: The number of multiply matched features, as an integer
  printSomeStats <- function(zeroMatches, multipleMatches) {
    cat("Number of common features: ", numCommonFeatures, "\n")
    cat("Recall:                    ", recall, "\n")
    cat("Precision:                 ", precision, "\n")
    cat("F_1 score:                 ", f1Score, "\n")
    cat("Number of zero matches:    ", zeroMatches, "\n")
    cat("Percent unmatched:          ", percentUnmatched, "%\n", sep = "")
    cat("Number of single matches:  ", numSinglyMatchedFeatures, "\n")
    cat("Percent single:             ", percentSinglyMatched, "%\n", sep = "")
    cat("Number of multiple matches:", multipleMatches, "\n")
    cat("Percent multiple:           ",
        percentMultiplyMatched, "%\n\n", sep = "")
  }

  cat("Using set 1 as the ground truth\n")
  printSomeStats(numUnmatchedFeaturesA, numMultiplyMatchedFeaturesAB)

  # c) Use the second feature set as the ground truth
  numCommonFeatures <- numSinglyMatchedFeatures + numMultiplyMatchedFeaturesAB
  recall <- numCommonFeatures / reticulate::py_to_r(featureSetB$size())
  precision <- numCommonFeatures / reticulate::py_to_r(featureSetA$size())
  f1Score <- (2 * precision * recall) / (precision + recall)

  percentUnmatched <- (numUnmatchedFeaturesB /
                         reticulate::py_to_r(featureSetB$size())) * 100
  percentSinglyMatched <- (numSinglyMatchedFeatures /
                             reticulate::py_to_r(featureSetB$size())) * 100
  percentMultiplyMatched <- (numMultiplyMatchedFeaturesBA /
                               reticulate::py_to_r(featureSetB$size())) * 100

  cat("Using set 2 as the ground truth\n")
  printSomeStats(numUnmatchedFeaturesB, numMultiplyMatchedFeaturesBA)
  return(invisible(NULL))
}

#' Specific feature retrieval.
#'
#' Gets information (in particular, retention time, mass-to-charge, and signal
#' intensity) for the feature at the given index in the given feature file. If
#' the feature does not exist, an error is returned (in the string).
#'
#' This function may not necessarily be directly helpful in comparing multiple
#' feature sets, but if a specific feature is known, or a rough idea of a
#' range of features is needed, this function serves as a quick lookup for
#' that information.
#'
#' @param featureFilePath The location of the featureXML file, as a string.
#' @param idx The index of the required feature in the set, as an integer.
#'
#' @return The specified feature's RT, m/z, and intensity, as a string. An
#' error message if the provided index is out of range.
#'
#' @examples
#' \dontrun{
#' getFeatureByIdx("inst/extdata/featureSetA.featureXML", 250)
#' }
#'
#' @export
#' @import reticulate
getFeatureByIdx <- function(featureFilePath, idx) {
  ropenms <- reticulate::import("pyopenms", convert = FALSE)
  retInfo <- NULL

  featureSet <- loadFeatureFile(featureFilePath)
  featureSet$sortByRT()
  if (idx < 1 | idx > reticulate::py_to_r(featureSet$size())) {
    retInfo <- paste("Error: cannot get feature from", featureFilePath,
                     "with index", idx, " (out of range)")
  }
  else
  {
    feature <- featureSet[idx]
    retInfo <- c(reticulate::py_to_r(feature$getRT()),
                 reticulate::py_to_r(feature$getMZ()),
                 reticulate::py_to_r(feature$getIntensity()))
  }

  return(retInfo)
}

# [END]
