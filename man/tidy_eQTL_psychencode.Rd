% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coloc_eQTL_tidying.R
\name{tidy_eQTL_psychencode}
\alias{tidy_eQTL_psychencode}
\title{Tidy psychencode eQTLs for use with coloc.}
\usage{
tidy_eQTL_psychencode(psychencode_QTLs, psychencode_SNP_info, add_colnames = T)
}
\arguments{
\item{psychencode_QTLs}{df. Full set of cis-eQTLs with no p-value or FDR
filtering as downloaded from \url{http://resource.psychencode.org/}.}

\item{psychencode_SNP_info}{df. SNP information for all QTLs considered,
including rsIDs (if available), location, and reference and alternate
alleles, as downloaded from \url{http://resource.psychencode.org/}.}

\item{add_colnames}{logical. Should columne names be added to the file?
Default is TRUE i.e. assumes user has not derived column names for full
summary statistic from another psychencode-derived file.}
}
\value{
Tidy dataframe with psychencode eQTLs. Note that the dataframe
  includes an additional column (SNP_rs). This is necessary for MAF lookup,
  but should be removed prior to submitting to
  \code{colochelpR::get_coloc_results()}.
}
\description{
Tidy psychencode eQTLs for use with coloc.
}
