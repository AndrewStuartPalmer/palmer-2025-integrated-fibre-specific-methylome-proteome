## Author: Andrew Palmer
## Title: proteomics load and filter data
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

####################################################################
## Proteomics load and filter 
####################################################################

# This script is the first step in the proteomic analysis workflow.
# It takes a raw, combined-protein TSV file as input (generated with msfragger).
# It filters proteins based on defined criteria (e.g., missing values, non-Homo sapiens).
# It outputs a filtered-combined-protein TSV file for subsequent analysis steps.


####################################################################
## Load required functions script
####################################################################

source("./proteome/scripts/proteomics-functions.R")

####################################################################
## Load Data ##
####################################################################

## data import of the combined protein file output from fragpipe
## The required file is a combined protein tsv file

data_raw <- read.delim("./proteome/data/processed-data/P22_0579_exp2_combined_protein.tsv", header = TRUE, sep = "\t")

data_raw_copy <- data_raw ## make a copy of the data

experiment_design <- read.delim("./proteome/experiment-design/lfq_analyst_experimental_design.txt")

####################################################################
## Data Filtering ##
####################################################################

## Remove contaminating proteins
## by removing any proteins not belonging to homo-sapiens
## and remove proteins that are not identified in two or
## more condition replicates

data_raw_filtered_subset <- filter.nonHomo.missVal(data_raw_copy,
                                                   experiment_design,
                                                   number1 = 1,
                                                   number2 = 4)

dim(data_raw_filtered_subset)

####################################################################
## Make duplicated gene names unique ##
####################################################################

data_raw_filtered_subset <-
  transform(data_raw_filtered_subset,
            protein.name = paste(data_raw_filtered_subset$Gene,
                                 data_raw_filtered_subset$Protein.ID,
                                 sep = "-"))

data_raw_filtered_subset$Gene <- replace(data_raw_filtered_subset$Gene,
                                         data_raw_filtered_subset$Gene == 0, NA)

data_raw_filtered_subset$Gene[is.na(data_raw_filtered_subset$Gene)] <-
  data_raw_filtered_subset$protein.name[is.na(data_raw_filtered_subset$Gene)]

####################################################################
## Save filtered data for loading in next r script ##
####################################################################

dir.create("./proteome/data/processed-data", recursive = TRUE)
write.table(data_raw_filtered_subset, "./proteome/data/processed-data/filtered_combined_protein.tsv",
            row.names = FALSE, sep = "\t")


