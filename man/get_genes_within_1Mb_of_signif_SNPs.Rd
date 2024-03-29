% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_genes_within_1Mb_signif_SNPs.R
\name{get_genes_within_1Mb_of_signif_SNPs}
\alias{get_genes_within_1Mb_of_signif_SNPs}
\title{Get genes within +/- 1Mb of significant SNPs}
\usage{
get_genes_within_1Mb_of_signif_SNPs(
  GWAS,
  pvalue_column,
  CHR_column,
  BP_column,
  pvalue_threshold = 5e-08,
  mart = 38
)
}
\arguments{
\item{GWAS}{dataframe. Dataframe with GWAS results, containing, as a minimum,
the following columns: (1) SNP identifier, (2) chromosome, in integer
format, (3) basepair position, in integer format and (4) p-value column.}

\item{pvalue_column}{chr. Name of column (in quotation marks) containing
p-values in GWAS dataframe.}

\item{CHR_column}{chr. Name of column (in quotation marks) containing
chromosome name in GWAS dataframe.}

\item{BP_column}{chr. Name of column (in quotation marks) containing basepair
position in GWAS dataframe.}

\item{pvalue_threshold}{num. P-value threshold for significance. Default
value is 5e-8.}

\item{mart}{int. Specify genome build.}
}
\value{
All genes within +/- 1Mb of significant SNPs.
}
\description{
As colochelpR has primarily been used together with cis-eQTL datasets (where
testing of SNP-gene associations is limited to genes within 1 Mb of the SNP),
we define a region of +/- 1 Mb.
}
