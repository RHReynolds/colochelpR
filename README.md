<!-- badges: start -->

[![Lifecycle:
dormant](https://img.shields.io/badge/lifecycle-dormant-blue.svg)](https://www.tidyverse.org/lifecycle/#dormant)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5011869.svg)](https://doi.org/10.5281/zenodo.5011869)

<!-- badges: end -->

# Coloc helper functions

Helper/wrapper functions for executing [coloc](https://github.com/chr1swallace/coloc) & associated analyses.

- Author: [David Zhang](https://github.com/dzhang32)
- Maintainer: [Regina H. Reynolds](https://github.com/RHReynolds) 
- Contributor(s): [Alejandro Martinez Carrasco](https://github.com/AMCalejandro)

## Citation

If you use any of the code from this repository, please cite the original [coloc R package](https://github.com/chr1swallace/coloc) and our [doi](https://doi.org/10.5281/zenodo.5011869).

## License

The code in this repository is released under an MIT license. This repository is distributed in the hope that it will be useful to the wider community, but without any warranty of any kind. Please see the [LICENSE](LICENSE.md) file for more details.

## Installation instructions

There is no plan to ever submit this code to `CRAN` or `Bioconductor`. This code was developed for use by the [Ryten Lab](https://github.com/rytenlab) and collaborators thereof. While most of the code has been separately tested by David and Regina, it is highly likely bugs exist. 

If you would like to install the development version from [GitHub](https://github.com/) you can use the following command:

``` r
if(!require("remotes"))
   install.packages("remotes") # if necessary
library(remotes)
install_github("RHReynolds/colochelpR")
```

## Background

The functions herein are simply wrapper functions for executing [coloc](https://github.com/chr1swallace/coloc). Please refer to the  [coloc vignette](https://chr1swallace.github.io/coloc/) for more background on use of coloc. Please remember to cite [coloc](https://github.com/chr1swallace/coloc), if you use this helper package.

## Usage example

For an example of how the wrapper functions in this package have been used in analyses, please refer to the following [repository](https://github.com/RHReynolds/RBD-GWAS-analysis).

### Note on column names

For these helper functions to work, it is important that GWASs/QTL datasets have been "tidied" and the necessary columns are present. So far, we have only used coloc for GWAS and eQTL datasets, the necessary columns for which are presented in tables below.

#### GWAS

**Column name** | **Description**
--------------- | -------------------------------------------------------------------- 
GWAS | GWAS name.
SNP | SNP identifier. The format of these will need to match the second GWAS/dataset you intend to run in coloc i.e. both location-based (CHR:BP) or RS-ID-based. **If location-based, ensure that you are mapping to the same genome build in both datasets.** 
beta | Regression coefficient.
se | Standard error of regression coefficient.
varbeta | This is generated through use of the function `colochelpR::get_varbeta()`, provided either `se` is available or `beta` and `t.stat`.
p.value | P-value of association.
Al1 | Effect allele.
Al2 | Alternate allele.
maf | Minor allele frequency.

- According to the original [coloc paper](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383), it is possible to compute posterior probabilities from single-variant association `p.value` and `maf`, but estimated single SNP regression co-efficients (`beta`) and their variances (`varbeta`) or standard errors (`se`) are preferred.
- If you do not have the columns `beta` and `se`, but you do have the odds ratio (often OR) and the Z-score, you can generate `beta` and `se` via the following formulas:

``` r
GWAS %>% 
  dplyr::mutate(beta = log(OR, base = exp(1)),
                se = log(OR, base = exp(1))/Z_SCORE)
```

- You will also need the number of individuals within your GWAS. 
    - This does not need to be a column in your GWAS dataframe. A vector with the number will suffice.
    - If the GWAS is a case-control dataset, you will additionally need the number of cases and controls, in order to calculate the proportion of individuals that are cases (i.e. n_cases/n_total).

#### eQTL

**Column name** | **Description**
--------------- | -------------------------------------------------------------------- 
eQTL_dataset | Name of eQTL dataset.
gene | Ensembl gene ID for gene regulated by SNP. Often the original eQTL dataset will provide HGNC symbols, or if you are using a microarray-based dataset, it might provide probe IDs. For HGNC symbols, you can either use `biomaRt`, or if you will be converting many gene ids, you can always use a `.GTF` for GRCh37/38 (depending on what your eQTL dataset was mapped to). For probe IDs, these will need to be mapped back to genes. This is may be provided by the eQTL generator, although for common microarray plates (e.g. Affymetrix arrays) (i) `biomaRt` can be used and (ii) Bioconductor offers offers a range of [annotation packages](http://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData) that can be used to convert probe IDs.
SNP | SNP identifier. The format of these will need to match the first GWAS/dataset you intend to run in coloc i.e. both location-based (CHR:BP) or RS-ID-based. **If location-based, ensure that you are mapping to the same genome build in both datasets.**
beta | Regression coefficient.
se | Standard error of regression coefficient.
varbeta | This is generated through use of the function `colochelpR::get_varbeta()`, provided either `se` is available or `beta` and `t.stat`.
p.value | P-value of association.
Al1 | Effect allele.
Al2 | Alternate allele.
maf | Minor allele frequency.
N | Number of individuals.

- As with the GWAS dataset, if columns `beta` and `se` are not available, `p.values` and `maf` can be used.
- Many eQTL datasets do not provide MAFs. It is, however, possible to use a reference database to look up MAFs for SNPs within the eQTL datasets. For example, the package [`MafDb.1Kgenomes.phase3.hs37d5`](https://bioconductor.org/packages/release/data/annotation/html/MafDb.1Kgenomes.phase3.hs37d5.html) contains MAFs for a number of populations (**remember to select the population that matches your eQTL dataset**). See example code below:

``` r
mafdb <- MafDb.1Kgenomes.phase3.hs37d5

example_df <- data.frame(rs_id = c("rs12921634", "rs1476958", "rs56189750"))

mafs <- GenomicScores::gscores(x = mafdb, ranges = example_df$rs_id %>% as.character(), pop = "EUR_AF")
mafs <- mafs %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "SNP") %>%
  dplyr::rename(maf = EUR_AF) %>% 
  dplyr::select(SNP, maf)

example_df <- example_df %>% 
  inner_join(mafs, by = c("rs_id" = "SNP"))
```

### Note on SNP notations

- The format that SNPs are provided in often varies between location-based (CHR:BP) or RS-ID-based formats. 
- As a result, there may sometimes be a need to convert between formats. 
- The most common conversion we have run into is conversion of RS IDs to chromosome locations, and have therefore written a function to do this (`colochelpR::convert_rs_to_loc()`). 
    - This functions uses the `BSgenome::snpsById()` function, and therefore requires a SNPlocs object, which is a container for storing known SNP locations. 
    - Use `BSgenome::available.SNPs()` to return a vector of the available SNPlocs packages. 
    - **Remember to use a SNPlocs package with the same genome build as the genome build you are intending to match these RS IDs to.**
- We have also written a wrapper function for conversion of chromosome locations to RS ids (`colochelpR::convert_loc_to_rs()`).

    - This functions uses the `BSgenome::snpsByOverlaps()` function, and therefore requires a SNPlocs object, which is a container for storing known SNP locations (as described above).
    - Users should note that some chromosome-bp locations may have more than one associated RS-ID, thus the user may have to check for duplicates after conversion. We have chosen not to do this, as some users may wish to simply remove duplicated SNPs, or keep one of the duplicates.
    - Users should also note that some chromosome-bp locations may not have an associated RS ID. These will be returned as NA in the RS ID column (`SNP`) after the conversion.
