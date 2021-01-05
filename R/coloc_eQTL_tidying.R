#' Tidy psychencode eQTLs for use with coloc.
#'
#' @param psychencode_QTLs df. Full set of cis-eQTLs with no p-value or FDR
#'   filtering as downloaded from \url{http://resource.psychencode.org/}.
#' @param psychencode_SNP_info df. SNP information for all QTLs considered,
#'   including rsIDs (if available), location, and reference and alternate
#'   alleles, as downloaded from \url{http://resource.psychencode.org/}.
#' @param add_colnames logical. Should columne names be added to the file?
#'   Default is TRUE i.e. assumes user has not derived column names for full
#'   summary statistic from another psychencode-derived file.
#'
#' @return Tidy dataframe with psychencode eQTLs. Note that the dataframe
#'   includes an additional column (SNP_rs). This is necessary for MAF lookup,
#'   but should be removed prior to submitting to
#'   \code{colochelpR::get_coloc_results()}.
#' @export
#' 

tidy_eQTL_psychencode <- function(psychencode_QTLs, psychencode_SNP_info, add_colnames = T){
  
  if(add_colnames == TRUE){
    
    # While column names are simply written in a vector below, these were orginally deduced from another psychencode-derived file (DER-08a_hg19_eQTL.significant.txt)
    # From the file description provided by psychencode, we know that no p-value/FDR filtering has occurred in the full summary statistics.
    # Therefore would not expect an FDR column. 
    # This reduces the column names derived from DER-08a_hg19_eQTL.significant.txt 14. 
    # This matches the number found in the full summary statistics.
    colnames(psychencode_QTLs) <- c("gene_id", "gene_chr", "gene_start", "gene_end", "strand", "number_of_SNPs_tested",
                                    "SNP_distance_to_TSS", "SNP_id", "SNP_chr", "SNP_start", "SNP_end", "nominal_pval",
                                    "regression_slope", "top_SNP")
  }
  
  tidy <- 
    psychencode_QTLs %>% 
    dplyr::mutate(gene_id = str_remove(gene_id, "\\..*"), 
                  eQTL_dataset = "psychencode",
                  # Number is retrieved from FAQ section of psychencode resource: https://faq.gersteinlab.org/category/capstone4/page/1/
                  # See post entitled: "Inquiry regarding PsychENCODE eQTL resource" from May 3, 2019.
                  N = 1387) %>% 
    dplyr::select(eQTL_dataset, gene = gene_id, SNP = SNP_id, beta = regression_slope, p.value = nominal_pval, N) %>% 
    dplyr::inner_join(psychencode_SNP_info %>% 
                        # As Psychencode is based on GTEx pipeline, and GTEx betas refer to ALT allele, A1 = ALT here
                        dplyr::select(SNP = PEC_id, SNP_rs = Rsid, Al1 = ALT, Al2 = REF))
  
  return(tidy)
  
}

#' Tidy eQTLGen eQTLs for use with coloc.
#'
#' @param eQTLGen_QTLs df. Full cis-eQTL summary statistics, as downloaded from
#'   \url{https://www.eqtlgen.org/cis-eqtls.html}. Should be post 2019-12-23
#'   changelog update.
#'
#' @return Tidy dataframe with eQTLGen eQTLs. Note that the dataframe includes
#'   an additional column (SNP_rs). This is necessary for MAF lookup, but should
#'   be removed prior to submitting to \code{colochelpR::get_coloc_results()}.
#' @export
#' 

tidy_eQTL_eQTLGen <- function(eQTLGen_QTLs){
  
  tidy <- 
    eQTLGen_QTLs %>% 
    dplyr::mutate(SNPloc = str_c(SNPChr, ":", SNPPos), 
                  eQTL_dataset = "eQTLGen") %>% 
    dplyr::select(eQTL_dataset, gene = Gene, SNP = SNPloc, SNP_rs = SNP, p.value = Pvalue, Al1 = AssessedAllele, Al2 = OtherAllele, N = NrSamples)
  
  return(tidy)
  
}