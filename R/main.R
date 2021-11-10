# TODO: separate code into files

ropenms <- reticulate::import("pyopenms", convert = FALSE)


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

#' Compares two sets of mass spectrometry features on a raw mass spectrometry
#' data set.
#'
#' @param featureFilePath1 The location of the first featureXML file.
#' @param featureFilePath2 The location of the second featureXML file.
#' @param rawDataFilePath The location of the mzML file.
#'
#' @examples
#' \dontrun{
#' compareFeatures("inst/extdata/featureSetA.featureXML",
#'                 "inst/extdata/featureSetB.featureXML",
#'                 "inst/extdata/20190122_HeLa_QC_Slot1-47_1_3228_800-860.mzML")
#' }
#'
#' @export
compareFeatures <- function(featureFilePath1, featureFilePath2, rawDataFilePath) {
  experiment <- ropenms$MSExperiment()
  ropenms$MzMLFile()$load(rawDataFilePath, experiment)
  featureSetA <- ropenms$FeatureMap()
  ropenms$FeatureXMLFile()$load(featureFilePath1, featureSetA)
  featureSetB <- ropenms$FeatureMap()
  ropenms$FeatureXMLFile()$load(featureFilePath2, featureSetB)
}
