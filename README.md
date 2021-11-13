
<!-- README.md is generated from README.Rmd. Please edit that file -->

# msFeatureCmp

<!-- badges: start -->
<!-- badges: end -->

## Description

The goal of msFeatureCmp is to provide a simple and fast way to evaluate
the performance of different mass spectrometry feature finding
algorithms, relative to one another. This package provides a way to
compare features (peptide signals) found by different algorithms,
perform some basic statistical analysis on the results, and visualize
the different sets of features and their locations on the raw data.
Existing MS tools (like OpenMS and MaxQuant) do not provide methods that
allow for cross-tool comparisons–likely because they are competitors–so
this package aims to make the process more accessible.

msFeatureCmp 0.1.0 was developed with R 4.1.1 in addition to Python
3.7.9 and the Python package pyOpenMS 2.7.0 on Windows 10 19043.1348.

## Installation

Before you can use this package, you need to ensure that you have a
version of Python from 3.5 to 3.7 installed, as well as the Python
package `pyOpenMS` (at least version 2.4.0). Python installers can be
downloaded from the Python downloads page
[here](https://www.python.org/downloads/), and `pyOpenMS` can be
installed from PyPI by running `pip install pyopenms` in a
Python-enabled shell (e.g. bash, cmd, etc.).

You can install the development version of msFeatureCmp from
[GitHub](https://github.com/ychen-uoft/msFeatureCmp) with:

``` r
install.packages("devtools")
devtools::install_github("ychen-uoft/msFeatureCmp", build_vignettes = TRUE)
library("msFeatureCmp")
```

To run the Shiny app:

``` r
Under construction
```

Raw data (used for examples) is located in the “inst/extdata” directory.
The three files (one mzML and two featureXMLs) need to be unzipped
before they can be used.

If R has trouble finding the correct Python environment (and you have
checked that all the Python prerequisites outlined above have been
installed properly), then you can set it by hand with
`Sys.setenv(RETICULATE_PYTHON = "<path to the correct Python binary>")`
before installing the package.

## Overview

``` r
ls("package:msFeatureCmp")
```

-   Describe each of the public APIs here

For more details, refer to the package vignettes.

``` r
browseVignettes("msFeatureCmp")
```

An overview of the package is illustrated below.

-   Provide an overview image (e.g. start with an mzML file and two
    featureXML files, and then show the possible options and outputs for
    each in a flowchart)
-   Images should be stored in “inst/extdata” and accessed with
    <!-- ![](./inst/extdata/image.png) -->

The package tree structure is provided below.

-   Provide the package tree (just use bash tree when the package is
    finalized)

## Contributions

The author of this package is Yijia Chen. All of the functions in this
package rely on the R package `reticulate` (which imports the Python
package `pyOpenMS`) in order to access the OpenMS mass spectrometry APIs
(e.g. to load the raw MS data into memory), but no other R packages
provide additional functionality (i.e. all algorithms were designed and
implemented by hand).

The example mzML file (raw MS data, in “inst/extdata”) was downloaded
from EMBL-EBI’s PRIDE archive as a Bruker (MaxQuant) d file, converted
to the open-source OpenMS mzML format using the Python package
`diapysef`, and cut to an arbitrary 60-second slice (in this case,
800s-860s) from the entire 2-hour run using OpenMS’s FileFilter tool.
The example featureXML file (containing found features, in
“inst/extdata”) “featureSetA” was obtained by running my own feature
finding algorithm (developed for my BCB330 project with the Röst Lab at
UofT) on the data, and the second example featureXML file “featureSetB”
was obtained by converting and filtering a tsv file, containing results,
that was provided with the raw data (which was obtained by running one
of MaxQuant’s feature finding algorithms on it).

-   Don’t forget to add whatever package was used for plotting

## References

Chen, Y., Xu, L., & Röst, H. (2020). *Improving feature finding in
LC-MS/MS runs by using additional ion mobility data*. \[Unpublished
BCB330Y1 final report\]. Department of Cell & Systems Biology,
University of Toronto, Toronto, Canada.

Meier, F., Brunner, A.-D., Frank, M., Ha, A., Bludau, I., Voytik, E.,
Kaspar-Schoenefeld, S., Lubeck, M., Raether, O., Bache, N., Aebersold,
R., Collins, B., Röst, H., & Mann, M. (2020). diaPASEF: parallel
accumulation-serial fragmentation combined with data-independent
acquisition. *Nat Methods, 17*, 1229-1236.
<https://doi.org/10.1038/s41592-020-00998-0>

Prianichnikov, N., Koch, H., Koch, S., Lubeck, M., Heilig, R., Brehmer,
S., Fischer, R., & Cox, J. (2020). MaxQuant software for ion mobility
enhanced shotgun proteomics. *Mol Cell Proteomics, 19*, 1058-1069.
<https://doi.org/10.1074/mcp.TIR119.001720>

Röst, H., Schmitt, U., Aebersold, R., & Malmström, M. (2014). pyOpenMS:
a Python-based interface to the OpenMS mass-spectrometry algorithm
library. *Proteomics, 14(1)*, 74-77.
<https://doi.org/10.1002/pmic.201300246>

Röst, H., Sachsenberg, T., Aiche, S., Bielow, C., Weisser, H., Aicheler,
F., Andreotti, S., Ehrlich, H.-C., Gutenbrunner, P., Kenar, E., Liang,
X., Nahnsen, S., Nilse, L., Pfeuffer, J., Rosenberger, G., Rurik, M.,
Schmitt, U., Veit, J., Walzer, M., … Kohlbacher, O. (2016). OpenMS: a
flexible open-source software platform for mass spectrometry data
analysis. *Nat Methods, 13*, 741-748.
<https://doi.org/10.1038/nmeth.3959>

Ushey, K., Allaire, J., & Tang, Y. (2021). *reticulate: interface to
‘Python’*. R package version 1.22.
<https://CRAN.R-project.org/package=reticulate>

## Acknowledgements

This package was developed as part of an assessment for 2021 BCB410H1:
Applied Bioinformatics, University of Toronto, Toronto, CANADA.

## Notes, to be removed before final check in

``` r
devtools::build_readme()
```

<img src="man/figures/README-pressure-1.png" width="100%" />

-   Try msFeatureCmp:::function() if not exported

-   usethis::use\_test(“test-testName”)

-   devtools::test() to test

-   usethis::use\_vignette(“vignette\_name”)

-   stop()

-   warning()

-   message()

-   ggplot2
