#' Convert rs ids to CHR:BP locations.
#'
#' @param df dataframe. Dataframe containing SNPs as rs ids.
#' @param SNP_column chr. Name of column (in quotation marks) containing rs ids
#'   in dataframe.
#' @param dbSNP BS genome reference snps (choose appropriate dbSNP build
#'   dependent on genome build).
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

#' Convert CHR:BP locations to rs ids.
#'
#' @description Function will convert genomic co-ordinates to rs ids.
#'
#' @section Warning:
#' \itemize{
#'   \item Some CHR:BP locations have more than one
#'   associated rs id, thus some filtering for duplicates may have to occur
#'   after conversion. We leave this to the user to decide how they wish to
#'   filter.
#'   \item Some CHR:BP locations may not have an associated rs id --
#'   these will be represented by NA in the SNP column after conversion.
#'   }
#'
#' @param df dataframe. Dataframe containing SNPs as genomic locations. Must
#'   contain 2 columns labelled \code{CHR} and \code{BP}, with chromosome and
#'   base pair positions, respectively. If the dataframe contains additional
#'   columns included these will be preserved.
#' @param dbSNP BS genome reference snps (choose appropriate dbSNP build
#'   dependent on genome build).
#'
#' @return Dataframe with rs ids and CHR:BP locations.
#' @export
#'

convert_loc_to_rs <- function(df, dbSNP){

  # If df CHR column has "chr" in name, remove
  if(stringr::str_detect(df$CHR[1], "chr")){

    df <-
      df %>%
      dplyr::mutate(CHR = stringr::str_replace(CHR, "chr", ""))

  }

  # If columns CHR are not correct format, this can cause problems with later join
  df <-
    df %>%
    dplyr::mutate(CHR = as.factor(CHR),
                  BP = as.integer(BP))

  # Convert df to GRanges object
  df_gr <-
    GenomicRanges::makeGRangesFromDataFrame(df,
                                            keep.extra.columns = FALSE,
                                            ignore.strand = TRUE,
                                            seqinfo = NULL,
                                            seqnames.field = "CHR",
                                            start.field = "BP",
                                            end.field = "BP",
                                            starts.in.df.are.0based = FALSE)

  # Genomic position object as dataframe with SNP locations converted to RS id.
  df_gr <-
    BSgenome::snpsByOverlaps(dbSNP, df_gr, minoverlap = 1L) %>%
    # Note that the default value for minoverlap is 0 which means that, by default, in addition to the SNPs that are
    # located within the genomic regions specified thru the ranges argument, snpsByOverlaps also returns SNPs that are
    # adjacent to these regions. Use minoverlap=1L to omit these SNPs.
    as.data.frame()

  combined <-
    df_gr %>%
    dplyr::rename(
      SNP = RefSNP_id,
      CHR = seqnames,
      BP = pos,
      # Rename strand column just in case this is a column in the inputted df
      gr_strand = strand) %>%
    dplyr::right_join(df, by = c("CHR", "BP")) %>%
    dplyr::select(-gr_strand, -alleles_as_ambig)

  return(combined)

}


