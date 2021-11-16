# This file contains the public visualization functions for the msFeatureCmp
# package.

#' Generates a data frame from a raw mass spectrometry dataset.
#'
#' Every data point (peak) is represented as a row containing its retention
#' time, mass-to-charge, and signal intensity.
#'
#' @param rawDataFilePath The location of the mzML file, as a string.
#'
#' @return The corresponding data frame for the dataset.
#'
#' @examples
#' \dontrun{
#' msFrame <- generateMSDataFrame(
#' "inst/extdata/20190122_HeLa_QC_Slot1-47_1_3228_800-810.mzML")
#' }
generateMSDataFrame <- function(rawDataFilePath) {
  ropenms <- reticulate::import("pyopenms", convert = FALSE)

  experiment <- loadMSFile(rawDataFilePath)
  numSpectra <- reticulate::py_to_r(experiment$getNrSpectra())
  numPeaks <- 0

  # Need to preallocate memory for the raw data matrix
  for (i in 0:(numSpectra - 1)) {
    spectrum <- experiment$getSpectrum(i)
    numPeaks <- numPeaks + reticulate::py_to_r(spectrum$size())
  }

  dataMatrix <- matrix(nrow = numPeaks, ncol = 3)
  currIdx <- 1

  # Add a new row in the matrix for each new data point (peak)
  for (i in 0:(numSpectra - 1)) {
    spectrum <- experiment$getSpectrum(i)
    peaks <- spectrum$get_peaks()

    # Each spectrum is a 2D slice of the 3D data
    pyMtZ <- peaks[0]
    pyInt <- peaks[1]
    massToCharges <- reticulate::py_to_r(pyMtZ)
    intensities <- reticulate::py_to_r(pyInt)

    for (j in 1:length(massToCharges)) {
      dataMatrix[currIdx, RT_IDX] <- reticulate::py_to_r(spectrum$getRT())
      dataMatrix[currIdx, MZ_IDX] <- massToCharges[j]
      dataMatrix[currIdx, IT_IDX] <- intensities[j]
      currIdx <- currIdx + 1
    }
  }

  # Convert the matrix into a data frame for plotting
  msFrame <- as.data.frame(dataMatrix)
  colnames(msFrame) <- c(c("RetentionTime", "MassToCharge", "Intensity"))
  return(msFrame)
}

#' Plots a raw mass spectrometry dataset onto a scatter plot with colour.
#'
#' The x-axis represents retention time, the y-axis represents mass-to-charge,
#' and the z-axis (the colour) represents intensity.
#'
#' This function may take a while to run.
#'
#' @param rawDataFilePath The location of the mzML file, as a string.
#'
#' @examples
#' \dontrun{
#' plotRawData("inst/extdata/20190122_HeLa_QC_Slot1-47_1_3228_800-810.mzML")
#' }
#'
#' @export
plotRawData <- function(rawDataFilePath) {
  msDataFrame <- generateMSDataFrame(rawDataFilePath)
  ggplot2::ggplot(msDataFrame,
                  ggplot2::aes(x = RetentionTime,
                               y = MassToCharge,
                               colour = Intensity)) + ggplot2::geom_point()
}

#' Generates a data frame from a feature set.
#'
#' Every feature is represented as a row containing its retention time, mass-
#' to-charge, and signal intensity.
#'
#' @param featureFilePath The location of the featureXML file, as a string.
#'
#' @return The corresponding data frame for the feature set.
#'
#' @examples
#' \dontrun{
#' featureFrame <- generateFeatureDataFrame(
#' "inst/extdata/featureSetA.featureXML")
#' }
generateFeatureDataFrame <- function(featureFilePath) {
  ropenms <- reticulate::import("pyopenms", convert = FALSE)

  featureSet <- loadFeatureFile(featureFilePath)
  featureMatrix <- convertFeaturesToSortedMatrix(featureSet)
  featureFrame <- as.data.frame(featureMatrix)
  colnames(featureFrame) <- c(c("RetentionTime", "MassToCharge", "Intensity"))
  return(featureFrame)
}

#' Plots a feature set onto a blank scatter plot with colour.
#'
#' The x-axis represents retention time, the y-axis represents mass-to-charge,
#' and the z-axis (the colour) represents intensity.
#'
#' This function may take a while to run.
#'
#' @param featureFilePath The location of the featureXML file, as a string.
#'
#' @examples
#' \dontrun{
#' plotSingleFeatureSet("inst/extdata/featureSetA.featureXML")
#' }
#'
#' @export
plotSingleFeatureSet <- function(featureFilePath) {
  featureDataFrame <- generateFeatureDataFrame(featureFilePath)
  ggplot2::ggplot(featureDataFrame,
                  ggplot2::aes(x = RetentionTime,
                               y = MassToCharge,
                               colour = Intensity)) + ggplot2::geom_point()
}

#' Plots two feature sets (overlapping) onto a blank scatter plot,
#' distinguished by colour.
#'
#' The x-axis represents retention time and the y-axis represents mass-to-
#' charge.
#'
#' This function may take a while to run.
#'
#' @param featureFilePath1 The location of the first featureXML file, as a
#' string.
#' @param featureFilePath2 The location of the second featureXML file, as a
#' string.
#'
#' @examples
#' \dontrun{
#' plotTwoFeatureSets("inst/extdata/featureSetA.featureXML",
#' "inst/extdata/featureSetB.featureXML")
#' }
#'
#' @export
plotTwoFeatureSets <- function(featureFilePath1, featureFilePath2) {
  # First convert both feature sets into data frames
  featureDataFrameA <- generateFeatureDataFrame(featureFilePath1)
  # Add an additional column to distinguish between the two sets
  typeColumn <- rep("A", nrow(featureDataFrameA))
  cbind(featureDataFrameA, typeColumn)
  # Make sure the column names are consistent
  colnames(featureDataFrameA) <- c(c("RetentionTime", "MassToCharge",
                                     "Intensity", "Type"))

  featureDataFrameB <- generateFeatureDataFrame(featureFilePath2)
  typeColumn <- rep("B", nrow(featureDataFrameB))
  cbind(featureDataFrameB, typeColumn)
  colnames(featureDataFrameB) <- c(c("RetentionTime", "MassToCharge",
                                     "Intensity", "Type"))

  # Add the second data frame to the first and plot the result
  featureDataFrame <- rbind(featureDataFrameA, featureDataFrameB)
  ggplot2::ggplot(featureDataFrame,
                  ggplot2::aes(x = RetentionTime,
                               y = MassToCharge,
                               colour = Type)) + ggplot2::geom_point()
}

#' Plots a feature set layered on top of its raw mass spectrometry dataset,
#' distinguished by colour.
#'
#' The x-axis represents retention time and the y-axis represents mass-to-
#' charge.
#'
#' This function may take a while to run.
#'
#' @param rawDataFilePath The location of the mzML file, as a string.
#' @param featureFilePath The location of the featureXML file, as a string.
#'
#' @examples
#' \dontrun{
#' plotFeatureSetOnRawData(
#' "inst/extdata/20190122_HeLa_QC_Slot1-47_1_3228_800-810.mzML",
#' "inst/extdata/featureSetA.featureXML")
#' }
#'
#' @export
plotFeatureSetOnRawData <- function(rawDataFilePath, featureFilePath) {
  # First convert the raw data and feature set into data frames
  msDataFrame <- generateMSDataFrame(rawDataFilePath)
  # R = raw data, F = feature
  typeColumn <- rep("R", nrow(msDataFrame))
  cbind(msDataFrame, typeColumn)
  colnames(msDataFrame) <- c(c("RetentionTime", "MassToCharge",
                               "Intensity", "Type"))

  featureDataFrame <- generateFeatureDataFrame(featureFilePath)
  typeColumn <- rep("F", nrow(featureDataFrame))
  cbind(featureDataFrame, typeColumn)
  colnames(featureDataFrame) <- c(c("RetentionTime", "MassToCharge",
                                    "Intensity", "Type"))

  totalDataFrame <- rbind(msDataFrame, featureDataFrame)
  ggplot2::ggplot(totalDataFrame,
                  ggplot2::aes(x = RetentionTime,
                               y = MassToCharge,
                               colour = Type)) + ggplot2::geom_point()
}
