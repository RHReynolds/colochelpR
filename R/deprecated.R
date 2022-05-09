
# Functions -------------------------------------------------------------------------------------------

#' Deprecated function: tidy GWASs
#'
#' This is only really meant to be used on the Ryten Lab server. If using
#' different GWASs/different servers, data tidying will have to be performed by
#' the user.
#'
#' @keywords internal
#' @rdname deprecated

get_GWAS_details <- function(){

  GWAS_details_df <-
    data_frame(
      GWAS_disease = c("PD", "AD", "SCZ", "MS", "HD_progression", "ALS", "HD_GEM", "SCZ_Euro_49", "AD_phase_3", "SCZ_2018_all", "PD_chang_ex_23andme", "intelligence_savage_2018"),
      GWAS_path = c("/data/LDScore/GWAS/PD2017_meta5/resultsForCojo_april17th2017.tab.gz",
                    "/data/LDScore/GWAS/AD/AD_all_sum_stats/IGAP_stage_1.txt",
                    "/data/LDScore/GWAS/SCZ/SCZ_all_sum_stats/ckqny.scz2snpres",
                    "/data/LDScore/GWAS/MS/MS_all_sum_stats/hg19_gwas_ms_imsgc_4_19_2.tab.txt",
                    "/data/LDScore/GWAS/HD_Progression_GWAS/TRACK_REGISTRY_final_RS.txt",
                    "/data/LDScore/GWAS/ALS_2016/Summary_Statistics_GWAS_2016_merged/ALS_2016_meta_sum_stats_all_chrs_merged.txt",
                    "/data/LDScore/GWAS/HD_GEM2015/gem.moa.snp.table.160907_tab_delim.txt",
                    "/data/LDScore/GWAS/SCZ/SCZ_EUR49/daner_PGC_SCZ49.sh2_mds10_1000G-frq_2.gz",
                    "/data/LDScore/GWAS/AD_Phase3_meta_sumstats_12_03_2018/Phase3_meta_sumstats_beta.gz",
                    "/data/LDScore/GWAS/SCZ2018/RS_clozuk_pgc2.meta.sumstats.txt",
                    "/data/LDScore/GWAS/PD2017_Chang_ex23andMe/META_no231_formatted_mappedtoRS.txt",
                    "/data/LDScore/GWAS/intelligence_savage_2018/SavageJansen_2018_intelligence_metaanalysis.txt"),
      n_cases = c(25987, 17008, 36989, 14498, 218+1773, 12577, 4082, 33640, 71880, 40675, 13708, 269867),
      n_controls = c(403177, 37154, 113075, 24091, 0, 23475, 0, 43456, 383378, 1474036, 95282, 0),
      cc_or_quant = c("cc", "cc", "cc", "cc", "quant", "cc", "quant", "cc", "cc", "cc", "cc", "quant"),
      rsid_or_chr_pos = c("chr_pos", "both", "both", "both", "both", "both", "both", "both", "both", "both", "both", "both"),
      freq_present = c(T, F, F, F, T, T, T, T, T, T, T, T),
      build = c("hg19", "hg19", "hg19", "hg19", "hg19", NA, NA, "hg19", "hg19", "hg19", NA, "hg19")) %>%
    dplyr::mutate(n_total = n_cases + n_controls,
           ratio_cases = n_cases/n_total)

  return(GWAS_details_df)

}

tidy_GWAS <- function(GWAS, GWAS_disease){

  GWAS_tidy <-
    switch (GWAS_disease,

            PD = GWAS %>%
              dplyr::mutate(freq = as.numeric(freq),
                     maf = ifelse(freq > 0.5, 1-freq, freq),
                     GWAS = GWAS_disease) %>%
              dplyr::select(GWAS, SNP, beta = b, se, p.value = p, Al1 = A1, Al2 = A2, maf),

            AD = GWAS %>%
              dplyr::mutate(GWAS = GWAS_disease) %>%
              dplyr::select(GWAS, SNP = MarkerName, beta = Beta, se = SE, p.value = Pvalue, Al1 = Effect_allele, Al2 = Non_Effect_allele),

            SCZ = GWAS %>%
              dplyr::mutate(or = as.numeric(or),
                     beta = log(or, base = exp(1)),
                     GWAS = GWAS_disease) %>%
              dplyr::select(GWAS, SNP = snpid, beta, se, p.value = p, Al1 = a1, Al2 = a2),

            MS = GWAS %>%
              tidyr::separate(col = "Alleles(Maj>Min)", into =  c("Al2", "Al1"), sep = ">") %>%
              dplyr::mutate(`OR(MinAllele)` = as.numeric(`OR(MinAllele)`),
                     LowerOR = as.numeric(LowerOR),
                     beta = log(`OR(MinAllele)`, base = exp(1)),
                     lower_beta = log(LowerOR, base = exp(1)),
                     se = abs((beta - lower_beta)/1.96),
                     GWAS = GWAS_disease) %>%
              dplyr::filter(!is.nan(se)) %>%
              dplyr::select(GWAS, SNP = Marker, beta, se, p.value = PValue, Al1, Al2),

            HD_progression = GWAS %>%
              dplyr::mutate(FRQ = as.numeric(FRQ),
                     maf = ifelse(FRQ > 0.5, 1-FRQ, FRQ),
                     A1 = stringr::str_to_upper(A1),
                     A2 = stringr::str_to_upper(A2),
                     GWAS = GWAS_disease) %>%
              dplyr::select(GWAS, SNP, beta = Beta, se = SE, p.value = P, Al1 = A1, Al2 = A2, maf),

            ALS = GWAS %>%
              dplyr::mutate(freq = as.numeric(freq),
                     maf = ifelse(freq > 0.5, 1-freq, freq),
                     a1 = stringr::str_to_upper(a1),
                     a2 = stringr::str_to_upper(a2),
                     GWAS = GWAS_disease) %>%
              dplyr::select(GWAS, SNP = snp, beta = b, se, p.value = p, Al1 = a1, Al2 = a2, maf),

            HD_GEM = GWAS %>%
              dplyr::mutate(freq = as.numeric(effect.allele.frequency),
                     maf = ifelse(freq > 0.5, 1-freq, freq),
                     GWAS = GWAS_disease) %>%
              dplyr::select(GWAS, SNP = snp, beta = effect.size, se, p.value = pval, Al1 = effect.allele, Al2 = other.allele, maf),

            SCZ_Euro_49 = GWAS %>%
              dplyr::mutate(freq = as.numeric(FRQ_A_33640),
                     maf = ifelse(freq > 0.5, 1-freq, freq),
                     OR = as.numeric(OR),
                     beta = log(OR, base = exp(1)),
                     GWAS = GWAS_disease) %>%
              dplyr::filter(INFO > 0.6) %>%
              dplyr::select(GWAS, SNP, beta = beta, se = SE, p.value = P, Al1 = A1, Al2 = A2, maf),

            AD_phase_3 = GWAS %>%
              dplyr::mutate(GWAS = GWAS_disease) %>%
              dplyr::select(GWAS, SNP, beta = BETA, se = SE, p.value = P, Al1 = A1, Al2 = A2, maf = MAF_HRC) %>%
              dplyr::filter(!is.na(maf)),

            SCZ_2018_all = GWAS %>%
              dplyr::mutate(or = as.numeric(OR),
                     beta = log(or, base = exp(1)),
                     freq = as.numeric(Freq.A1),
                     maf = ifelse(freq > 0.5, 1-freq, freq),
                     GWAS = GWAS_disease) %>%
              dplyr::filter(!is.nan(beta)) %>%
              dplyr::select(GWAS, CHR, BP, SNP, beta, se = SE, p.value = P, Al1 = A1, Al2 = A2, maf),

            PD_chang_ex_23andme = GWAS %>%
              dplyr::mutate(beta = as.numeric(b),
                     se = as.numeric(se),
                     freq = as.numeric(freq),
                     maf = ifelse(freq > 0.5, 1-freq, freq),
                     GWAS = GWAS_disease) %>%
              dplyr::filter(!is.nan(beta)) %>%
              dplyr::select(GWAS, CHR, BP, SNP = RefSNP_id, beta, se, p.value = p, Al1 = A1, Al2 = A2, maf),

            intelligence_savage_2018 = GWAS %>%
              dplyr::mutate(beta = as.numeric(stdBeta),
                     se = as.numeric(SE),
                     freq = as.numeric(EAF_HRC),
                     maf = ifelse(freq > 0.5, 1-freq, freq),
                     GWAS = GWAS_disease,
                     A1 = stringr::str_to_upper(A1),
                     A2 = stringr::str_to_upper(A2)) %>%
              dplyr::filter(!is.nan(beta)) %>%
              dplyr::select(GWAS, CHR, POS, SNP, beta, se, p.value = P, Al1 = A1, Al2 = A2, maf)


    )

  return(GWAS_tidy)

}
