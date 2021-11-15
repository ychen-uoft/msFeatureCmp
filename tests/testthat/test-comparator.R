# Comparator integration tests

test_that("feature comparator pipeline works", {
  skip("Too long")

  cmpOutput <- capture.output(
    msFeatureCmp::compareFeatures(
      system.file("extdata",
                  "20190122_HeLa_QC_Slot1-47_1_3228_800-810.mzML",
                  package = "msFeatureCmp"),
      system.file("extdata",
                  "featureSetA.featureXML",
                  package = "msFeatureCmp"),
      system.file("extdata",
                  "featureSetB.featureXML",
                  package = "msFeatureCmp")))
  outputStr <- paste(cmpOutput, collapse = "\n")

  # refInput <- readLines(system.file("extdata", "pipeline-output.txt",
  #                                   package = "msFeatureCmp"))
  # expect_equal(outputStr, refInput)
})

test_that("feature retrieval works", {
  featureInfo <- msFeatureCmp::getFeatureByIdx(
    system.file("extdata", "featureSetA.featureXML", package = "msFeatureCmp"),
    250)

  expect_equal(featureInfo[RT_IDX], 803.795599611922)
  expect_equal(featureInfo[MZ_IDX], 814.443167554736)
  expect_equal(featureInfo[IT_IDX], 4485.009765625)
})
