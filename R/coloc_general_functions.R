# Author and source: David Zhang (GitHub: https://github.com/dzhang32)
# Minor tweaks by: RHReynolds

# Functions -------------------------------------------------------------------------------------------

##### Top level #####

#' Make results directory
#'
#' Generates output folder to store coloc results, does not overwrite the folder
#' if it already exists.
#'
#' @param results_path chr. Path to folder.
#' @param folder_name chr. Name of folder to make.
#'
#' @return The full path to the new folder.
#' @export
#'

make_results_dir <- function(results_path, folder_name){

  results_dir_path <- stringr::str_c(results_path, "/", folder_name)

  if(!dir.exists(results_dir_path)){

    dir.create(results_dir_path)

  }else {

    print(stringr::str_c(results_dir_path, " directory already exists.."))

  }

  return(results_dir_path)

}

#' Get varbeta
#'
#' Generates varbeta from the beta and the t.stat/se
#'
#' @param df Dataframe with the columns beta and t.stat or se.
#'
#' @return Dataframe with the column varbeta added as a calculation from the
#'   beta and se/t.stat.
#' @export
#'

get_varbeta <- function(df){

  if("varbeta" %in% colnames(df)){

    stop("varbeta column already present")

  }else if(!any(c("t.stat", "se") %in% colnames(df))){

    stop("No t.stat or se column found to generate varbeta")

  }

  if("se" %in% colnames(df)){

    df[["varbeta"]] <- as.numeric(df[["se"]]) ^ 2

  }else if("t.stat" %in% colnames(df)){

    if(!"beta" %in% colnames(df)) stop("No beta column found to calculate se using t.stat")

    df[["varbeta"]] <- (as.numeric(df[["beta"]])/as.numeric(df[["t.stat"]])) ^ 2

  }

  return(df)

}

#' Check coloc data format
#'
#' Checks that the necessary columns for coloc are present and are the correct
#' type (num, character) currently for the betas requires: (i) SNP - chr/pos or
#' rsids used to join the two datasets; (ii) beta; (iii) varbeta; (iv) maf
#' (option depending on the check_maf argument).
#'
#' @param df Dataframe to check the columns of
#' @param beta_or_pval Either "beta" or "pval" depending on whether you want to
#'   use betas or pvals in the coloc analysis.
#' @param check_maf lgl. Should function check if the maf is present?
#'
#' @return Inputted dataframe with columns the correct type.
#' @export
#'

check_coloc_data_format <- function(df, beta_or_pval, check_maf){

  # always check if SNP present
  if(!"SNP" %in% colnames(df)) stop("Warning: SNP column not found")
  if(!is.character(df[["SNP"]])){

    print("SNP column is not a character (probably a factor), converting to character...")

    df[["SNP"]] <- as.character(df[["SNP"]])

  }

  if(beta_or_pval == "beta"){

    if(!"beta" %in% colnames(df)) stop("Warning: beta column not found")
    if(!"varbeta" %in% colnames(df)) stop("Warning: varbeta column not found")

    if(!is.numeric(df[["beta"]])){

      print("Converting beta column to numeric...")

      df[["beta"]] <- as.numeric(df[["beta"]])

    }

    if(!is.numeric(df[["varbeta"]])){

      print("Converting varbeta column to numeric...")

      df[["varbeta"]] <- as.numeric(df[["varbeta"]])

    }

  }else if(beta_or_pval == "pval"){

    if(!"p.value" %in% colnames(df)) stop("Warning: p.value column not found")

    if(!is.numeric(df[["p.value"]])){

      print("Converting p.value column to numeric...")

      df[["p.value"]] <- as.numeric(df[["p.value"]])

    }

  }

  if(check_maf == T){

    if(!"maf" %in% colnames(df)) stop("Warning: maf column not found")

    if(any(is.na(df[["maf"]]))) stop("Warning: maf column contains missing values -- these must be removed")

    if(!is.numeric(df[["maf"]])){

      print("Converting maf column to numeric...")

      df[["maf"]] <- as.numeric(df[["maf"]])

    }

    if(any(df[["maf"]] > 0.5)) stop("some mafs > 0.5")

  }

  return(df)

}

#' Run coloc
#'
#' Wrapper function that joins the two formatted datasets and performs coloc
#' analysis then annotates coloc results currently only allows for df1_type to
#' be "quant" or "cc" and df2_type to be "quant" so if running GWAS vs eQTL,
#' eQTL must be df2. Note: "quant" refers to a quantitative trait, while "cc"
#' refers to "case-control".
#'
#' @param df1 df. First formatted dataset for coloc analysis.
#' @param df2 df. Second formatted dataset for coloc analysis.
#' @param harmonise logical. Whether you would like to modify the sign of the
#'   beta so both datasets are with respect to the same allele
#' @param df1_type Either "quant" or "cc".
#' @param df2_type Only "quant" currently.
#' @param df1_beta_or_pval Either "beta" or "pval" depending on whether you want
#'   to use betas or pvals in the coloc analysis
#' @param df2_beta_or_pval Either "beta" or "pval" depending on whether you want
#'   to use betas or pvals in the coloc analysis
#' @param df1_N num. Number of samples in df1 (only required for quant)
#' @param df2_N num. Number of samples in df2 (only required for quant)
#' @param df_1_propor_cases num. Option to put in the proportion of cases for a
#'   GWAS in \code{coloc::coloc.abf()}.
#' @param annotate_signif_SNP_df1_df2 lgl. Whether to annotate results with the
#'   signif SNPs from df1/df2 and from coloc
#' @param key_cols all columns needed to uniquely identify analysis - requires
#'   appending "_1" or "_2" depending on whether col from df1 or df2
#'   respectively
#' @param df_1_name chr. df1 name that will be added to results list (e.g.
#'   "GWAS")
#' @param df_2_name chr. df2 name that will be added to results list (e.g.
#'   "eQTL")
#' @param df1_path chr. Path to df1
#' @param df2_path chr. Path to df2
#' @param p1 num. Prior probability a SNP is associated with trait 1, set to
#'   coloc.abf default of 1e-4
#' @param p2 num. Prior probability a SNP is associated with trait 2, set to
#'   coloc.abf default of 1e-4
#' @param p12 num. Prior probability a SNP is associated with both traits, set
#'   to coloc.abf default of 1e-5
#'
#' @return List containing coloc results annotated or NULL if there are no
#'   overlapping SNPs or all SNPs are removed through harmonisation
#' @export
#'

get_coloc_results <- function(df1, df2, harmonise = F, df1_type, df2_type, df1_beta_or_pval, df2_beta_or_pval, df1_N = NA, df2_N = NA, df_1_propor_cases,
                              annotate_signif_SNP_df1_df2 = F, key_cols, df_1_name, df_2_name, df1_path, df2_path,
                              p1 = 1e-04, p2 = 1e-04, p12 = 1e-05){

  df1_df2_joined <-
    join_coloc_datasets(df1 = df1, df2 = df2, harmonise = harmonise)

  if(is.null(df1_df2_joined)){

    print("No overlapping SNPs between df1 and df2")

    return(NULL)

  }

  coloc_results <-
    run_coloc_abf(df1_df2_joined = df1_df2_joined,
                  df1_type = df1_type,
                  df2_type = df2_type,
                  df1_beta_or_pval = df1_beta_or_pval,
                  df2_beta_or_pval = df2_beta_or_pval,
                  df1_N = df1_N,
                  df2_N = df2_N,
                  df_1_propor_cases = df_1_propor_cases,
                  p1 = p1, p2 = p2, p12 = p12)

  coloc_results_annotated <-
    annotate_coloc_results(coloc_results = coloc_results,
                           df1_df2_joined = df1_df2_joined,
                           annotate_signif_SNP_df1_df2 = annotate_signif_SNP_df1_df2,
                           key_cols = key_cols,
                           df_1_name = df_1_name,
                           df_2_name = df_2_name,
                           df1_path = df1_path,
                           df2_path = df2_path)

  return(coloc_results_annotated)

}

#' Save coloc results
#'
#' @param coloc_results_annotated Results, as returned by
#'   \code{colochelpR::get_coloc_results}.
#' @param results_dir_path chr. The directory in which to save the merged
#'   results.
#'
#' @return Saves coloc results in provided directory.
#' @export
#'

save_coloc_results <- function(coloc_results_annotated, results_dir_path){

  if(!is.null(coloc_results_annotated)){

    keys <-
      coloc_results_annotated[["keys"]] %>%
      unlist() %>%
      stringr::str_c(collapse = "_")

    save(coloc_results_annotated,
         file = stringr::str_c(results_dir_path, "/", keys, ".rda"))

  }else{

    print("No coloc results found, likely due to no overlapping SNPs... ")

  }


}

#' Merge coloc summaries
#'
#' Merges the raw annotated coloc results into a summary dataframe.
#'
#' @param results_folder_path chr. Path to the folder where the results are
#'   stored
#' @param add_signif_SNP lgl. Determines whether to add the details of the most
#'   signif SNPs to the summary. Default is FALSE.
#' @param recursive lgl. Argument used in \code{list.files}, which determines
#'   whether to recursively look for results in every other folder within the
#'   results_folder_path. Default is FALSE.
#' @param pattern lgl. Argument used in \code{list.files}, which specifies a
#'   pattern to match for the results (e.g. "\\.rda"). Default is NULL.
#'
#' @return The dataframe with the summaries merged.
#' @export
#'

merge_coloc_summaries <- function(results_folder_path, add_signif_SNP = F, recursive = F, pattern = NULL){

  coloc_result_paths <-
    list.files(path = results_folder_path, recursive = recursive, pattern = pattern, full.names = T)

  coloc_results_summaries_w_keys <-
    coloc_result_paths %>%
    purrr::map(~extract_coloc_summary(coloc_result_path = ., add_signif_SNP))

  coloc_results_summaries_w_keys_all <-
    do.call(dplyr::bind_rows, coloc_results_summaries_w_keys)

  coloc_results_summaries_w_keys_all_w_sum_ratio_PPH4_PPH3 <-
    coloc_results_summaries_w_keys_all %>%
    dplyr::mutate(sum_PPH3_PPH4 = PP.H3.abf + PP.H4.abf,
                  ratio_PPH4_PPH3 = PP.H4.abf/PP.H3.abf)

  return(coloc_results_summaries_w_keys_all_w_sum_ratio_PPH4_PPH3)

}

# work in progress...
summarise_coloc_results <- function(coloc_results_merged){

  num_tests_total <- coloc_results_merged %>% nrow()
  PPH4_ab_0.75 <- coloc_results_merged %>% dplyr::filter(PP.H4.abf >= 0.75) %>% nrow()
  PPH4_ab_0.9 <- coloc_results_merged %>% dplyr::filter(PP.H4.abf >= 0.9) %>% nrow()
  sum_PPH3_PPH4_ab_0.8_ratio_PPH4_PPH3_ab_2 <- coloc_results_merged %>% dplyr::filter(sum_PPH3_PPH4 >= 0.8, ratio_PPH4_PPH3 >= 2) %>% nrow()
  sum_PPH3_PPH4_ab_0.8_ratio_PPH4_PPH3_ab_5 <- coloc_results_merged %>% dplyr::filter(sum_PPH3_PPH4 >= 0.8, ratio_PPH4_PPH3 >= 5) %>% nrow()

  summary_results <-
    data_frame(`Cut off` = c("Number tests total", "PPH4 above 0.75", "PPH4 above 0.90", "PPH4 + PPH3 above 0.8 & ratio PPH4/PPH3 above 2", "PPH4 + PPH3 above 0.8 & ratio PPH4/PPH3 above 5"),
             `Number of results` = c(num_tests_total, PPH4_ab_0.75, PPH4_ab_0.9, sum_PPH3_PPH4_ab_0.8_ratio_PPH4_PPH3_ab_2, sum_PPH3_PPH4_ab_0.8_ratio_PPH4_PPH3_ab_5)) %>%
    ggpubr::ggtexttable(rows = NULL)

  return(summary_results)

}

add_most_signif_SNP_to_signif_results <- function(coloc_results_merged, results_folder_path){


  return(summary_results)

}



##### Second level #####



#' Join coloc datasets
#'
#' Wrapper function that joins the two formatted datasets.
#'
#' @param df1 df. First formatted dataset for coloc analysis.
#' @param df2 df. Second formatted dataset for coloc analysis.
#' @param harmonise logical. Whether you would like to modify the sign of the
#'   beta so both datasets are with respect to the same allele
#'
#' @return Joined dataframe with data from both datasets.
#' @export
#'

join_coloc_datasets <- function(df1, df2, harmonise = F){

  df1_overlap <-
    df1 %>%
    dplyr::semi_join(df2, by = "SNP") %>%
    dplyr::arrange(SNP)

  df2_overlap <-
    df2 %>%
    dplyr::semi_join(df1, by = "SNP") %>%
    dplyr::arrange(SNP)

  # incase of no SNPs overlapping
  if(nrow(df2_overlap) == 0){

    return(NULL)

  }

  if(!identical(df1_overlap[["SNP"]], df2_overlap[["SNP"]])) stop("SNPs do not match")
  if(any(duplicated(df1_overlap[["SNP"]]))) stop("Some SNPs are duplicated")

  if(harmonise == T){

    df1_overlap[["REF_Al_1"]] <- df2_overlap[["Al1"]]
    df1_overlap[["REF_Al_2"]] <- df2_overlap[["Al2"]]

    df1_overlap_harmonised_list <-
      match_alleles(df1_overlap,
                    A1.ref= "REF_Al_1", A2.ref= "REF_Al_2",
                    A1.data = "Al1", A2.data = "Al2", BETA.data= "beta", flip = TRUE)

    df1_overlap_harmonised <-
    df1_overlap_harmonised_list %>%
    .[[1]] %>%
      dplyr::mutate(REF_Al_1 = NULL,
             REF_Al_2 = NULL) %>%
      dplyr::arrange(SNP)

    # incase of all overlapping SNPs removed through harmonisation
    if(nrow(df1_overlap_harmonised) == 0){

      return(NULL)

    }

    df2_overlap_harmonised <-
      df2_overlap %>%
      dplyr::semi_join(df1_overlap_harmonised, by = "SNP") %>%
      dplyr::arrange(SNP)

    colnames(df1_overlap_harmonised) <- stringr::str_c(colnames(df1_overlap_harmonised), "_1")
    colnames(df2_overlap_harmonised) <- stringr::str_c(colnames(df2_overlap_harmonised), "_2")

    df1_df2_joined <-
      inner_join(df1_overlap_harmonised, df2_overlap_harmonised, by = c("SNP_1" = "SNP_2")) %>%
      dplyr::rename(SNP = SNP_1)

    return(df1_df2_joined)

  }

  colnames(df1_overlap) <- stringr::str_c(colnames(df1_overlap), "_1")
  colnames(df2_overlap) <- stringr::str_c(colnames(df2_overlap), "_2")

  df1_df2_joined <-
    inner_join(df1_overlap, df2_overlap, by = c("SNP_1" = "SNP_2")) %>%
    dplyr::rename(SNP = SNP_1)

  return(df1_df2_joined)

}

run_coloc_abf <- function(df1_df2_joined, df1_type, df1_beta_or_pval, df2_type, df2_beta_or_pval, df1_N = NA, df2_N = NA, df_1_propor_cases,
                          p1 = p1, p2 = p2, p12 = p12){

  if(df1_type == "cc" && df2_type == "quant"){

    if(is.na(df2_N) ) stop("No N provided for df2 - needed for quant trait")

    # here if the GWAS MAF is NULL, it doesnt seem to change the result, tested changing the MAF of the GWAS and makes no difference

    if(df1_beta_or_pval == "beta" & df2_beta_or_pval == "beta"){
      coloc_results <-
        coloc.abf(dataset1 = list( type = "cc",
                                   snp = df1_df2_joined[["SNP"]],
                                   beta = df1_df2_joined[["beta_1"]],
                                   varbeta = df1_df2_joined[["varbeta_1"]],
                                   MAF = df1_df2_joined[["maf_1"]],
                                   s = df_1_propor_cases,
                                   N = df1_N),
                  dataset2 = list( type = "quant",
                                   snp = df1_df2_joined[["SNP"]],
                                   beta = df1_df2_joined[["beta_2"]],
                                   varbeta = df1_df2_joined[["varbeta_2"]],
                                   MAF = df1_df2_joined[["maf_2"]],
                                   N = df2_N),
                  p1 = p1, p2 = p2, p12 = p12
        )

    }else if(df1_beta_or_pval == "beta" & df2_beta_or_pval == "pval"){
      coloc_results <-
        coloc.abf(dataset1 = list( type = "cc",
                                   snp = df1_df2_joined[["SNP"]],
                                   beta = df1_df2_joined[["beta_1"]],
                                   varbeta = df1_df2_joined[["varbeta_1"]],
                                   MAF = df1_df2_joined[["maf_1"]],
                                   s = df_1_propor_cases,
                                   N = df1_N),
                  dataset2 = list( type = "quant",
                                   snp = df1_df2_joined[["SNP"]],
                                   pvalues = df1_df2_joined[["p.value_2"]],
                                   MAF = df1_df2_joined[["maf_2"]],
                                   N = df2_N),
                  p1 = p1, p2 = p2, p12 = p12
        )

    }else if(df1_beta_or_pval == "pval" & df2_beta_or_pval == "beta"){
      coloc_results <-
        coloc.abf(dataset1 = list( type = "cc",
                                   snp = df1_df2_joined[["SNP"]],
                                   pvalues = df1_df2_joined[["p.value_1"]],
                                   MAF = df1_df2_joined[["maf_1"]],
                                   s = df_1_propor_cases,
                                   N = df1_N),
                  dataset2 = list( type = "quant",
                                   snp = df1_df2_joined[["SNP"]],
                                   beta = df1_df2_joined[["beta_2"]],
                                   varbeta = df1_df2_joined[["varbeta_2"]],
                                   MAF = df1_df2_joined[["maf_2"]],
                                   N = df2_N),
                  p1 = p1, p2 = p2, p12 = p12
        )

    }else{
      coloc_results <-
        coloc.abf(dataset1 = list( type = "cc",
                                   snp = df1_df2_joined[["SNP"]],
                                   pvalues = df1_df2_joined[["p.value_1"]],
                                   MAF = df1_df2_joined[["maf_1"]],
                                   s = df_1_propor_cases,
                                   N = df1_N),
                  dataset2 = list( type = "quant",
                                   snp = df1_df2_joined[["SNP"]],
                                   pvalues = df1_df2_joined[["p.value_2"]],
                                   MAF = df1_df2_joined[["maf_2"]],
                                   N = df2_N),
                  p1 = p1, p2 = p2, p12 = p12
        )

    }

  } else if(df1_type == "quant" && df2_type == "quant"){

    if( is.na(df1_N) | is.na(df2_N) ) stop("No N provided for df1 or df2 - needed for quant trait")

    if(df1_beta_or_pval == "beta" & df2_beta_or_pval == "beta"){
      coloc_results <-
        coloc.abf(dataset1 = list( type = "quant",
                                   snp = df1_df2_joined[["SNP"]],
                                   beta = df1_df2_joined[["beta_1"]],
                                   varbeta = df1_df2_joined[["varbeta_1"]],
                                   MAF = df1_df2_joined[["maf_1"]],
                                   N = df1_N),
                  dataset2 = list( type = "quant",
                                   snp = df1_df2_joined[["SNP"]],
                                   beta = df1_df2_joined[["beta_2"]],
                                   varbeta = df1_df2_joined[["varbeta_2"]],
                                   MAF = df1_df2_joined[["maf_2"]],
                                   N = df2_N),
                  p1 = p1, p2 = p2, p12 = p12
        )

    }else if(df1_beta_or_pval == "beta" & df2_beta_or_pval == "pval"){
      coloc_results <-
        coloc.abf(dataset1 = list( type = "quant",
                                   snp = df1_df2_joined[["SNP"]],
                                   beta = df1_df2_joined[["beta_1"]],
                                   varbeta = df1_df2_joined[["varbeta_1"]],
                                   MAF = df1_df2_joined[["maf_1"]],
                                   N = df1_N),
                  dataset2 = list( type = "quant",
                                   snp = df1_df2_joined[["SNP"]],
                                   pvalues = df1_df2_joined[["p.value_2"]],
                                   MAF = df1_df2_joined[["maf_2"]],
                                   N = df2_N),
                  p1 = p1, p2 = p2, p12 = p12
        )

    }else if(df1_beta_or_pval == "pval" & df2_beta_or_pval == "beta"){
      coloc_results <-
        coloc.abf(dataset1 = list( type = "quant",
                                   snp = df1_df2_joined[["SNP"]],
                                   pvalues = df1_df2_joined[["p.value_1"]],
                                   MAF = df1_df2_joined[["maf_1"]],
                                   N = df1_N),
                  dataset2 = list( type = "quant",
                                   snp = df1_df2_joined[["SNP"]],
                                   beta = df1_df2_joined[["beta_2"]],
                                   varbeta = df1_df2_joined[["varbeta_2"]],
                                   MAF = df1_df2_joined[["maf_2"]],
                                   N = df2_N),
                  p1 = p1, p2 = p2, p12 = p12
        )

    }else{
      coloc_results <-
        coloc.abf(dataset1 = list( type = "quant",
                                   snp = df1_df2_joined[["SNP"]],
                                   pvalues = df1_df2_joined[["p.value_1"]],
                                   MAF = df1_df2_joined[["maf_1"]],
                                   N = df1_N),
                  dataset2 = list( type = "quant",
                                   snp = df1_df2_joined[["SNP"]],
                                   pvalues = df1_df2_joined[["p.value_2"]],
                                   MAF = df1_df2_joined[["maf_2"]],
                                   N = df2_N),
                  p1 = p1, p2 = p2, p12 = p12
        )

    }

  }

  return(coloc_results)

}

annotate_coloc_results <- function(coloc_results, df1_df2_joined, annotate_signif_SNP_df1_df2 = F,
                                   key_cols, df_1_name, df_2_name, df1_path, df2_path){

  if(!any(key_cols %in% colnames(df1_df2_joined))) print("key columns not found")

  coloc_results[["df1_name"]] <- df_1_name

  coloc_results[["df2_name"]] <- df_2_name

  coloc_results[["df1_path"]] <- df1_path

  coloc_results[["df2_path"]] <- df2_path

  coloc_results[["keys"]] <- df1_df2_joined[, key_cols][1,]

  coloc_results[["signif_coloc_SNP"]] <-
    coloc_results[["results"]] %>%
    dplyr::filter(SNP.PP.H4 == max(SNP.PP.H4)) %>%
    dplyr::select(snp, SNP.PP.H4) %>%
    left_join(df1_df2_joined, by = c("snp" = "SNP"))

  if(annotate_signif_SNP_df1_df2 == T){

    if(!all(c("p.value_1", "p.value_2") %in% colnames(df1_df2_joined))) stop("No pvalue found to get signif SNP details for df1/df2")

    coloc_results[["signif_df1_SNP"]] <-
      df1_df2_joined %>%
      dplyr::filter(as.numeric(p.value_1) == min(as.numeric(p.value_1))) %>%
      dplyr::select(SNP, contains("_1"))

    coloc_results[["signif_df2_SNP"]] <-
      df1_df2_joined %>%
      dplyr::filter(as.numeric(p.value_2) == min(as.numeric(p.value_2))) %>%
      dplyr::select(SNP, contains("_2"))

  }

  return(coloc_results)

}

extract_coloc_summary <- function(coloc_result_path, add_signif_SNP = F){

  load(coloc_result_path)

  coloc_results_summary_df <-
    coloc_results_annotated[["summary"]] %>%
    as.matrix() %>%
    t() %>%
    as.data.frame()

  coloc_results_summary_df_w_keys <-
    bind_cols(coloc_results_annotated[["keys"]], coloc_results_summary_df)

  df1_name <- coloc_results_annotated[["df1_name"]]
  df2_name <- coloc_results_annotated[["df2_name"]]

  if(add_signif_SNP == T){

    signif_coloc_snp <-
      coloc_results_annotated[["signif_coloc_SNP"]] %>%
      dplyr::filter(!duplicated(p.value_1), !duplicated(p.value_2)) %>%
      dplyr::select(snp, SNP.PP.H4,
                    beta_1, p.value_1,
                    beta_2, p.value_2)

    signif_df1_snp <-
      coloc_results_annotated[["signif_df1_SNP"]] %>%
      dplyr::filter(!duplicated(as.numeric(p.value_1))) %>%
      dplyr::select(SNP, beta_1, p.value_1)

    signif_df2_snp <-
      coloc_results_annotated[["signif_df2_SNP"]] %>%
      dplyr::filter(!duplicated(as.numeric(p.value_2))) %>%
      dplyr::select(SNP, beta_2, p.value_2)

    colnames(signif_coloc_snp) <-
      colnames(signif_coloc_snp) %>%
      stringr::str_c("signif_coloc_snp_", .)

    colnames(signif_df1_snp) <-
      colnames(signif_df1_snp) %>%
      stringr::str_c("signif_", df1_name, "_", .)

    colnames(signif_df2_snp) <-
      colnames(signif_df2_snp) %>%
      stringr::str_c("signif_", df2_name, "_", .)

    coloc_results_summary_df_w_signif_snps <-
      coloc_results_summary_df_w_keys %>%
      bind_cols(signif_coloc_snp, signif_df1_snp, signif_df2_snp)

    return(coloc_results_summary_df_w_signif_snps)

  }

  return(coloc_results_summary_df_w_keys)

}

##### Third level (courtesy of moloc package) #####

match_alleles <- function(data, A1.ref="A1.ref", A2.ref="A2.ref",  A1.data = "A1", A2.data = "A2", BETA.data="BETA", flip = TRUE) {

  match_correct = data[,A1.ref] == data[,A1.data] & data[,A2.ref]== data[,A2.data]
  match_flip = data[,A1.ref] == data[,A2.data] & data[,A2.ref] == data[,A1.data]
  match_comp_one = data[,A1.ref] == complement_snp(data[,A1.data]) & data[,A2.ref]== complement_snp(data[,A2.data])
  match_comp_two = data[,A1.ref] == complement_snp(data[,A2.data]) & data[,A2.ref] == complement_snp(data[,A2.data])
  snp_allele_match = match_flip | match_correct | match_comp_one | match_comp_two
  message(sum(snp_allele_match), " SNPs out of ", length(snp_allele_match), " had the correct alleles, discarding SNPs without the correct alleles")

  if (flip) {
    if (any(which(match_flip)>0)) {
      data[match_flip, A1.data]=data[match_flip, A1.ref]
      data[match_flip, A2.data]=data[match_flip, A2.ref]
      data[match_flip, BETA.data]=-data[match_flip, BETA.data]
    }
  }

  removed = data[!snp_allele_match,]
  data = data[snp_allele_match,]

  return(list(data, removed))

}

complement_snp <- function(x){

  as = x =="A"
  ts = x == "T"
  gs = x == "G"
  cs = x == "C"
  ins = x == "I"
  dels = x == "D"
  x[as] = "T"
  x[ts] = "A"
  x[gs] = "C"
  x[cs] = "G"
  x[ins] = "NA"
  x[dels] = "NA"

  return(x)

}

change_indels <- function(data) {

  data$A2[nchar(data$A1)>1] = "D"
  data$A1[nchar(data$A1)>1] = "I"
  data$A1[nchar(data$A2)>1] = "D"
  data$A2[nchar(data$A2)>1] = "I"

  return(data)

}
