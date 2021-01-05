#' Convert rs ids to CHR:BP locations.
#'
#' @param df dataframe. Dataframe containing SNPs as rs ids.
#' @param SNP_column chr. Name of column (in quotation marks) containing rs ids
#'   in dataframe.
#' @param dbSNP BS genome reference snps (choose appropriate dbSNP build
#'   dependent on genome build)
#'
#' @return Dataframe with rs ids and CHR:BP locations.
#' @export
#'

convert_rs_to_loc <- function(df, SNP_column, dbSNP){

  rs <- BSgenome::snpsById(dbSNP, df[[SNP_column]], ifnotfound = "drop") %>%
    as.data.frame() %>%
    tidyr::unite(col = "loc", seqnames, pos, sep = ":", remove = T) %>%
    dplyr::rename(rs = RefSNP_id) %>%
    dplyr::select(rs, loc)

  filter_vector <- c("rs")
  names(filter_vector) <- SNP_column

  df <- df %>%
    dplyr::inner_join(rs, by = filter_vector)

  return(df)

}
