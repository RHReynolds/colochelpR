# colochelpR 0.99.2

User visible changes:

* `get_genes_within_1Mb_of_signif_SNPs()` - added arguement, `pvalue_threshold` to permit users to enter their own significance threshold should they wish to. Default set to 5e-8, such that the function is backwards compatible with previous versions.

Bug fixes:

* [[#2](https://github.com/RHReynolds/colochelpR/issues/2)] - prevent error in `convert_loc_to_rs()` upon join if `strand` column supplied by user
* [[#4]](https://github.com/RHReynolds/colochelpR/issues/4) - prevent error in coloc's `check_dataset()` caused by default `NA` supplied by `colochelpR::get_coloc_results()` when `N` is not supplied by user for a case control trait 

# colochelpR 0.99.0

Initial release.
