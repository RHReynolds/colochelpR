#' Get genes within +/- 1Mb of significant SNPs
#'
#' @param GWAS dataframe. Dataframe with GWAS results, containing, as a minimum,
#'   the following columns: (1) SNP identifier, (2) chromosome, in integer
#'   format, (3) basepair position, in integer format and (4) p-value column.
#' @param pvalue_column chr. Name of column (in quotation marks) containing
#'   p-values in GWAS dataframe.
#' @param CHR_column chr. Name of column (in quotation marks) containing
#'   chromosome name in GWAS dataframe.
#' @param BP_column chr. Name of column (in quotation marks) containing basepair
#'   position in GWAS dataframe.
#' @param mart int. Specify genome build.
#'
#' @return All genes within +/- 1Mb of significant SNPs.
#'

get_genes_within_1Mb_of_signif_SNPs <- function(GWAS,
                                                pvalue_column,
                                                CHR_column,
                                                BP_column,
                                                mart = 38){
  # tidy evaluation
  pvalue_column_var <- rlang::sym(pvalue_column)
  CHR_column_var <- rlang::sym(CHR_column)
  BP_column_var <-  rlang::sym(BP_column)

  GWAS_signif_SNPs <-
    GWAS %>%
    dplyr::filter(!!pvalue_column_var <= 5e-8)

  GWAS_signif_SNPs_max_min_bp <-
    GWAS_signif_SNPs %>%
    dplyr::group_by(!!CHR_column_var) %>%
    dplyr::filter(!!BP_column_var == max(!!BP_column_var) | !!BP_column_var == min(!!BP_column_var)) %>%
    dplyr::mutate(signif_1Mb_window_min = !!BP_column_var - 1000000,
                  signif_1Mb_window_max = !!BP_column_var + 1000000) %>%
    dplyr::group_by(!!CHR_column_var) %>%
    dplyr::summarise(seqnames = !!CHR_column_var %>% unique(),
                     start = signif_1Mb_window_min %>% min(),
                     end = signif_1Mb_window_max %>% max())

  ensembl_all_genes_start_stop <-
    .query_biomart(mart = mart, attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"), filter = "ensembl_gene_id", values = "")

  ensembl_all_genes_start_stop_gr <-
    GenomicRanges::GRanges(
      seqnames = ensembl_all_genes_start_stop[["chromosome_name"]],
      ranges =IRanges::IRanges(start =  ensembl_all_genes_start_stop[["start_position"]], end = ensembl_all_genes_start_stop[["end_position"]]),
      strand = "*",
      ensembl_gene_id = ensembl_all_genes_start_stop[["ensembl_gene_id"]])

  GWAS_signif_SNPs_max_min_bp_gr <-
    GenomicRanges::GRanges(
      seqnames = GWAS_signif_SNPs_max_min_bp$seqnames,
      ranges = IRanges::IRanges(start = GWAS_signif_SNPs_max_min_bp$start,
                       end = GWAS_signif_SNPs_max_min_bp$end),
      strand = "*"
    )

  overlaps_hits <-
    GenomicRanges::findOverlaps(GWAS_signif_SNPs_max_min_bp_gr, ensembl_all_genes_start_stop_gr, minoverlap = 1, type = "any")

  ensembl_all_genes_start_stop_gr_overlapping <- ensembl_all_genes_start_stop_gr[overlaps_hits %>% subjectHits()]

  ensembl_gene_ids_overlapping_1Mb_window_hit <- ensembl_all_genes_start_stop_gr_overlapping$ensembl_gene_id %>%
    unique()

  return(ensembl_gene_ids_overlapping_1Mb_window_hit)

}


#' Vector-based biomart query.
#'
#' @param mart Specify genome build.
#' @param attributes Vector of attributes to extract from BioMart.
#' @param filter Vector of filter to be used for BioMart query.
#' @param values Values of the filter, e.g. vector of ensembl gene IDs.
#'
#' @return Original dataframe together with extracted biomart query.
#'

.query_biomart <- function(mart = 38, attributes, filter, values){

  if(mart != 38 && mart != 37) stop("Mart must be 38 or 37...")

  if(mart == 38){

    ensembl_mart <-
      biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

  }else if(mart == 37){

    ensembl_mart <-
      biomaRt::useMart(host = "grch37.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

  }

  # BioMart search
  biomart_query <- biomaRt::getBM(attributes = attributes, filters = filter, values = values , mart = ensembl_mart)

  return(biomart_query)

}

#' Dataframe-based biomart query.
#'
#' @param dataframe Dataframe with gene names.
#' @param columnToFilter Name of column in dataframe, which contains gene names.
#' @param mart Specify genome build.
#' @param attributes Vector of attributes to extract from BioMart.
#' @param filter Vector of filter to be used for BioMart query.
#'
#' @return Original dataframe together with extracted biomart query.
#' @export
#'

biomart_df <- function(dataframe, columnToFilter, mart = 38, attributes, filter){

  if(mart != 38 && mart != 37) stop("Mart must be 38 or 37...")

  if(mart == 38){

    ensembl_mart <-
      biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

  }else if(mart == 37){

    ensembl_mart <-
      biomaRt::useMart(host = "grch37.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

  }

  # Query genes as a vector
  genes <- dataframe %>% .[[columnToFilter]] %>% unique()
  print(str_c("Number of unique genes to search: ", length(genes)))

  # BioMart search
  biomart_query <- biomaRt::getBM(attributes = attributes, filters = filter, values = genes , mart = ensembl_mart)
  print(str_c("Number of matches found:", nrow(biomart_query)))

  # Create new data frame with positional information + remainder of the original dataframe
  # First requires creating join vector for the by argument in inner_join
  join_vector <- filter
  names(join_vector) <- columnToFilter
  merged <- dplyr::inner_join(dataframe, biomart_query, by = join_vector)

  return(merged)

}

