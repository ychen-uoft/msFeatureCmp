# This file contains the primary points of access for the msFeatureCmp package.

# The global pyOpenMS python import
ropenms <- reticulate::import("pyopenms", convert = FALSE)

# Other global constants
RT_IDX <- 1  # Retention time index
MZ_IDX <- 2  # Mass-to-charge index
IT_IDX <- 3  # Signal intensity index

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
  # Load the raw MS run first, followed by both feature sets
  experiment <- ropenms$MSExperiment()
  ropenms$MzMLFile()$load(rawDataFilePath, experiment)
  featureSetA <- ropenms$FeatureMap()
  ropenms$FeatureXMLFile()$load(featureFilePath1, featureSetA)
  featureSetB <- ropenms$FeatureMap()
  ropenms$FeatureXMLFile()$load(featureFilePath2, featureSetB)
}

# TODO: add functions to access the plotting APIs
