## Author: Andrew Palmer
## Title: proteomics load and transform data
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

####################################################################
## This code loads and transforms filtered combined protein data
####################################################################

# This script is the second step in the proteomic analysis workflow.
# It takes a filtered_raw_data and log2 transforms it.
# It then performs centre median normalisation.
# followed by left downshifted global imputation
# It outputs a data_normalised_imputed tsv file for subsequent analysis steps.

# Source custom helper functions from 'proteomics-functions.R'

source("./proteome/scripts/proteomics-functions.R")

####################################################################
## Load Data ##
####################################################################

## data import of the filtered combined protein file output from 01- protoemics
## script

data_raw <- read.delim("./proteome/data/processed-data/filtered_combined_protein.tsv",
                       header = TRUE,
                       sep = "\t")

data_raw_filtered_subset <- data_raw

####################################################################
## Select MaxLFQ columns and Log transform them ##
####################################################################

## Extract names of intensity columns
intensity_names <- grep(".Intensity|!MaxLFQ",
                        colnames(data_raw_filtered_subset), value = TRUE)

## Extract only MAXlfq columns
intensity_names <- grep("MaxLFQ", (intensity_names),
                        value = TRUE,
                        invert = FALSE)

## Assign column names for log2-transformed data
log_names <- paste(intensity_names, "log2", sep = ".")

##subset data frame to include only needed columns
needed_columns <- c("Protein.ID", "Gene", intensity_names)

data_raw_filtered_subset <- data_raw_filtered_subset[, needed_columns]


data_raw_filtered_subset[intensity_names] <-
  sapply(data_raw_filtered_subset[intensity_names], as.numeric)

####################################################################
## Transform Data ##
####################################################################

## Log 2 Transformation of MaxLFQ data
data_raw_filtered_subset[log_names] <-
  log2(data_raw_filtered_subset[intensity_names])

####################################################################
## Data Normalisation ##
####################################################################

data_normalised <- centre.median(data_raw_filtered_subset)

####################################################################
## Missing Value Imputation ##
####################################################################

log2_names <- grep("log2", names(data_normalised), value = TRUE)

data_normalised_imputed <- impute.normal.global(data_normalised[, log2_names])

####################################################################
## Clean data and Save ##
####################################################################

data_normalised_imputed <- as.data.frame(data_normalised_imputed)

gene <- data_raw$Gene

rownames(data_normalised_imputed) <- gene

data_normalised_imputed$proteinID <- data_raw$Protein.ID

write.table(data_normalised_imputed, "./proteome/data/processed-data/data_normalised_imputed.tsv",
            row.names = TRUE, sep = "\t",
            col.names = NA)


