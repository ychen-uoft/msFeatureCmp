% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{loadFeatureFile}
\alias{loadFeatureFile}
\title{Feature set loader.}
\usage{
loadFeatureFile(filePath)
}
\arguments{
\item{filePath}{The location of the featureXML file, as a string.}
}
\value{
The in-memory representation of the feature set, as a FeatureMap.
}
\description{
Loads a mass spectrometry feature set, containing found features, into
memory as a FeatureMap object. The feature file must be in featureXML
(OpenMS) format.
}
\examples{
\dontrun{
featureSet <- loadFeatureFile("inst/extdata/featureSetA.featureXML")
featureSet$size()  # Returns 720
}

}
