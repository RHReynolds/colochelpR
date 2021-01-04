<!-- badges: start -->

[![Lifecycle:
dormant](https://img.shields.io/badge/lifecycle-dormant-blue.svg)](https://www.tidyverse.org/lifecycle/#dormant)

<!-- badges: end -->

# Coloc helper functions

Helper/wrapper functions for executing [coloc](https://github.com/chr1swallace/coloc) & associated analyses.

- Author: [David Zhang](https://github.com/dzhang32)
- Maintainer: [Regina H. Reynolds](https://github.com/RHReynolds) 

## Installation instructions

There is no plan to ever submit this code to `CRAN` or `Bioconductor`. This code was developed for use by the [Ryten Lab](https://github.com/rytenlab) and collaborators thereof. While most of the code has been separately tested by David and Regina, it is highly likely bugs exist. 

If you would like to install the development version from [GitHub](https://github.com/) you can use the following command:

``` r
if(!require("remotes"))
   install.packages("remotes") # if necessary
library(remotes
install_github("RHReynolds/colochelpR")
```

## Background

The functions herein are simply wrapper functions for executing [coloc](https://github.com/chr1swallace/coloc). Please refer to the  [coloc vignette](https://chr1swallace.github.io/coloc/) for more background on use of coloc.

## Example

For an example of how the wrapper functions in this package have been used in analyses, please refer to the following repository: *To be added*
