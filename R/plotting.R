# plotting.R
# Package: msFeatureCmp
# Author: Yijia Chen
# Date: 2021-12-04
# Version: 0.1.0

# This file contains all the public visualization functions for the
# msFeatureCmp package.

#' Mass spectrometry data frame generator.
#'
#' Generates a data frame from a raw mass spectrometry dataset. Every data
#' point/peak is represented as a row containing its retention time, mass-to-
#' charge, and signal intensity.
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
#'
#' @import reticulate
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

#' Raw data plotter.
#'
#' Plots a raw mass spectrometry dataset onto a scatter plot.
#'
#' The x-axis represents retention time (in sec), the y-axis represents mass-
#' to-charge (in Th), and the z-axis (with darker colours for greater values)
#' represents signal intensity (unitless).
#'
#' This function may take a while to run.
#'
#' @param rawDataFilePath The location of the mzML file, as a string.
#'
#' @return The raw data scatter plot.
#'
#' @examples
#' \dontrun{
#' plotRawData("inst/extdata/20190122_HeLa_QC_Slot1-47_1_3228_800-810.mzML")
#' }
#'
#' @export
#' @import ggplot2
plotRawData <- function(rawDataFilePath) {
  msDataFrame <- generateMSDataFrame(rawDataFilePath)
  rawDataPlot <- ggplot2::ggplot(msDataFrame,
                                 ggplot2::aes(x = RetentionTime,
                                              y = MassToCharge,
                                              colour = Intensity)) +
    ggplot2::geom_point()

  return(rawDataPlot)
}

#' Data frame generator.
#'
#' Generates a data frame from a feature set. Every feature is represented as a
#' row, containing its retention time, mass-to-charge, and signal intensity, in
#' that order.
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
#'
#' @import reticulate
generateFeatureDataFrame <- function(featureFilePath) {
  ropenms <- reticulate::import("pyopenms", convert = FALSE)

  featureSet <- loadFeatureFile(featureFilePath)
  featureMatrix <- convertFeaturesToSortedMatrix(featureSet)
  featureFrame <- as.data.frame(featureMatrix)
  colnames(featureFrame) <- c(c("RetentionTime", "MassToCharge", "Intensity"))
  return(featureFrame)
}

#' Single feature set plotter.
#'
#' Plots a single feature set onto a blank scatter plot.
#'
#' The x-axis represents retention time (in sec), the y-axis represents mass-
#' to-charge (in Th), and the z-axis (with darker colours for greater values)
#' represents signal intensity (unitless).
#'
#' This function may take a while to run.
#'
#' @param featureFilePath The location of the featureXML file, as a string.
#'
#' @return The feature set scatter plot.
#'
#' @examples
#' \dontrun{
#' plotSingleFeatureSet("inst/extdata/featureSetA.featureXML")
#' }
#'
#' @export
#' @import ggplot2
plotSingleFeatureSet <- function(featureFilePath) {
  featureDataFrame <- generateFeatureDataFrame(featureFilePath)
  featureDataPlot <- ggplot2::ggplot(featureDataFrame,
                                     ggplot2::aes(x = RetentionTime,
                                                  y = MassToCharge,
                                                  colour = Intensity)) +
    ggplot2::geom_point()

  return(featureDataPlot)
}

#' Overlapping feature set plotter.
#'
#' Plots two feature sets, at the same time, onto a blank scatter plot. The two
#' feature sets are distinguished by colour.
#'
#' The x-axis represents retention time (in sec) and the y-axis represents
#' mass-to-charge (in Th).
#'
#' This function may take a while to run.
#'
#' @param featureFilePath1 The location of the first featureXML file, as a
#' string.
#' @param featureFilePath2 The location of the second featureXML file, as a
#' string.
#'
#' @return The overlapping feature scatter plot.
#'
#' @examples
#' \dontrun{
#' plotTwoFeatureSets("inst/extdata/featureSetA.featureXML",
#' "inst/extdata/featureSetB.featureXML")
#' }
#'
#' @export
#' @import ggplot2
plotTwoFeatureSets <- function(featureFilePath1, featureFilePath2) {
  # First convert both feature sets into data frames
  featureDataFrameA <- generateFeatureDataFrame(featureFilePath1)
  # Add an additional column to distinguish between the two sets
  Type <- rep("A", nrow(featureDataFrameA))
  featureDataFrameA <- cbind(featureDataFrameA, Type)

  featureDataFrameB <- generateFeatureDataFrame(featureFilePath2)
  Type <- rep("B", nrow(featureDataFrameB))
  featureDataFrameB <- cbind(featureDataFrameB, Type)

  # Add the second data frame to the first and plot the result
  featureDataFrame <- rbind(featureDataFrameA, featureDataFrameB)
  featureDataPlot <- ggplot2::ggplot(featureDataFrame,
                                     ggplot2::aes(x = RetentionTime,
                                                  y = MassToCharge,
                                                  colour = Type)) +
    ggplot2::geom_point()

  return(featureDataPlot)
}

#' Layered feature set plotter.
#'
#' Creates a scatter plot, where the given feature set is layered on top of its
#' raw mass spectrometry dataset. The data is distinguished from the features
#' by colour.
#'
#' The x-axis represents retention time (in sec) and the y-axis represents
#' mass-to-charge (in Th).
#'
#' This function may take a while to run.
#'
#' @param rawDataFilePath The location of the mzML file, as a string.
#' @param featureFilePath The location of the featureXML file, as a string.
#'
#' @return The overlayed feature/raw data scatter plot.
#'
#' @examples
#' \dontrun{
#' plotFeatureSetOnRawData(
#' "inst/extdata/20190122_HeLa_QC_Slot1-47_1_3228_800-810.mzML",
#' "inst/extdata/featureSetA.featureXML")
#' }
#'
#' @export
#' @import ggplot2
plotFeatureSetOnRawData <- function(rawDataFilePath, featureFilePath) {
  # First convert the raw data and feature set into data frames
  msDataFrame <- generateMSDataFrame(rawDataFilePath)
  # R = raw data, F = feature
  Type <- rep("R", nrow(msDataFrame))
  msDataFrame <- cbind(msDataFrame, Type)

  featureDataFrame <- generateFeatureDataFrame(featureFilePath)
  Type <- rep("F", nrow(featureDataFrame))
  featureDataFrame <- cbind(featureDataFrame, Type)

  totalDataFrame <- rbind(msDataFrame, featureDataFrame)
  totalDataPlot <- ggplot2::ggplot(totalDataFrame,
                                   ggplot2::aes(x = RetentionTime,
                                                y = MassToCharge,
                                                colour = Type)) +
    ggplot2::geom_point()

  return(invisible(NULL))
}

# [END]
