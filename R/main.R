ropenms <- reticulate::import("pyopenms", convert = FALSE)


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
