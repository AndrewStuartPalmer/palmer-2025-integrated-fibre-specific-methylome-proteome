## Author: Andrew Palmer
## Title: proteomics functions
## 2Date: 2024-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
## See the 'LICENSE' file in the root directory for full license details.

####################################################################
## This code contains custom functions for use in proteomic analysis.
####################################################################

####################################################################
## Filter missing values and non homo sapien proteins ##
####################################################################

## This function can be used to filter proteomic samples prior to
## normalisation and missing value imputation on log 2 Max lfq intensity
## values. User provides a combined protein tsv and experimental design
## data frame output from msfragger.

## NOTE: This function relies on specific column names and filtering criteria
## unique to the current proteomic dataset. Review internal parameters before reuse.

filter.nonHomo.missVal <- function(df_x, df_y, number1 = 1, number2 =1) {

  filtered <- df_x[grep("Homo", df_x$Organism), ]

  intensity_names <- grep(".Intensity|!MaxLFQ", colnames(df_x), value = TRUE)

  intensity_names <- grep("MaxLFQ", (intensity_names), value = TRUE,
                          invert = FALSE)
    holding <- filtered[, intensity_names]
    holding <- t(holding)
    holding <- as.data.frame(holding)
    holding <- holding > 0
    holding <- apply(holding, 2, as.numeric)
    holding <- as.data.frame(holding)
    df_y <- df_y[order(df_y$label),]
    holding$condition <- df_y$description

  missingness <- aggregate(holding[, colnames(holding)
                                   [colnames(holding) != "condition"]],
                           list(holding$condition), FUN = sum)

  missval <- t(missingness)
  colnames(missval) <- missval[1, ]
  missval <- missval[-1, ]
  missval <- as.data.frame(missval)
  missval[] <- lapply(missval, as.numeric)

     filter_subset <- missval$Myh1_pooled_fibres_post > number1 | missval$Myh1_pooled_fibres_pre > number1 |  missval$Myh2_pooled_fibres_post > number2 | missval$Myh2_pooled_fibres_pre > number2 | missval$Myh7_pooled_fibres_post > number2 | missval$Myh7_pooled_fibres_pre > number2


    data_raw_filtered_subset <- subset(filtered, filter_subset)

}



####################################################################
## Data normalisation function ##
####################################################################

## This function can be used to normalise proteomic samples prior to
## missing value imputation on log 2 Max lfq intensity values 

centre.median <- function(df_x) {
  ## df_x = data frame containing log2 columns for normalisation
  log2_names <- grep("log2", names(df_x), value = TRUE)

  df_x[, log2_names] <- lapply(log2_names,
                               function(x) {
                                 log2 <- df_x[[x]]
                                 log2[!is.finite(log2)] <- NA
                                 get_median <- median(log2, na.rm = TRUE)
                                 log2 - get_median
                               }
                             )

  return(df_x)
}

####################################################################
## Data imputation function ##
####################################################################

## This function can be used to impute missing values in proteomic
## samples. The imputation function draws random numbers from a normal
## distribution that is left shifted. This assumes that missing values
## are missing not at random (MNAR), but rather are missing due to low
## abundance. This is a global approach which is slightly more conservative 
## than groupwise imputation.

impute.normal.global <- function(df_x, width = 0.3, downshift = 1.8, seed = 100) {

  if (!is.matrix(df_x)) df_x <- as.matrix(df_x)

  set.seed(seed)
  df_x <- apply(df_x, 2, function(temp) {
    temp[!is.finite(temp)] <- NA
    temp_sd <- stats::sd(temp, na.rm = TRUE)
    temp_mean <- mean(temp, na.rm = TRUE)
    shrunk_sd <- width * temp_sd   # shrink sd width
    downshifted_mean <- temp_mean - downshift * temp_sd  #shift mean of imputed values
    n_missing <- sum(is.na(temp))
    temp[is.na(temp)] <- stats::rnorm(n_missing,
                                      mean = downshifted_mean,
                                      sd = shrunk_sd)
    temp
  })
  return(df_x)

}

####################################################################
## Function to extract coefficient of choice and annotate them ##
####################################################################

## function to extract all test coefficients in limma 
## adds annotations for uniprot, ensembl and symbol

extract.fit.proteome <- function(coef = coef, fit, df.m) {

    diff <- topTable(fit, coef = coef,  number = nrow(df.m), adjust.method = "BH")

    diff$symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, rownames(diff), keytype =  "SYMBOL", column = "SYMBOL")

    diff$uniprot <- AnnotationDbi::mapIds(org.Hs.eg.db, rownames(diff), keytype =  "SYMBOL", column = "UNIPROT")

    diff$ensembl <- AnnotationDbi::mapIds(org.Hs.eg.db, rownames(diff), keytype =  "SYMBOL", column = "ENSEMBL")
   
    return(diff)
}


## This file is saved as a proteomics-function script
## can be loaded during analysis.
