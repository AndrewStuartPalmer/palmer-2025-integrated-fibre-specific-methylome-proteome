## Author: Andrew Palmer
## Title: Methylation sample filtering
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

######################################################################
## Script for filtering 4wk samples and averaging replicate samples ##
## Preparation of CpG annotations                                   ##
######################################################################


####################################################################
## Load required packges ##
####################################################################

library("data.table")

####################################################################
## Load Data ##
####################################################################

## using fread from the data table package is quicker
## load m_vals
m_val <-  data.frame(data.table::fread("./methylome/data/processed-data/full-dasen-normalised-filtered-m.csv"), row.names = 1)

beta <-  data.frame(fread("./methylome/data/processed-data/full-dasen-normalised-filtered-beta.csv"), row.names =1)

## load phenotype data

pheno <- read.delim("./methylome/experiment-design/pheno-table.txt")

## combine cell_type and timepoint for limma analysis

pheno$group <- paste(pheno$cell_type, pheno$timepoint, sep = "_")

####################################################################
## Subset data to exclude 4wWM and pooled fibre samples ##
####################################################################

## take average m values for replicate samples
m_val_rep <- m_val[,c(17,18,19,90,91)]
m_val$WM.f28.preWM.1 <- rowMeans(m_val_rep[,4:5])
m_val$WM.s201.preWM.1 <- rowMeans(m_val_rep[,1:3])

## remove replicate samples and keep only average m values
pheno_filt <- pheno[-c(1,2,17,18,30,57,58,90),] 

m_val_filt <- m_val[,- c(1,2,17,18,30,57,58,90)]
row.names(pheno_filt) <- 1:nrow(pheno_filt)

## same for beta values

beta_rep <- beta[,c(17,18,19,90,91)]
beta$WM.f28.preWM.1 <- rowMeans(beta_rep[,4:5])
beta$WM.s201.preWM.1 <- rowMeans(beta_rep[,1:3])

beta_filt <- beta[,- c(1,2,17,18,30,57,58,90)]

## remove 4wWM samples

subset_position <- !endsWith(colnames(m_val_filt), "4wWM")
m_val_filt_subset <- subset(m_val_filt, select = subset_position)
pheno_filt_subset <- subset(pheno_filt, subset = subset_position)
beta_filt_subset <- subset(beta_filt, select = subset_position)

pheno_filt_subset <- pheno_filt_subset[-c(57),] 

m_val_filt_subset  <- m_val_filt_subset[,- c(57)]
beta_filt_subset <- beta_filt_subset[,- c(57)]


####################################################################
## Save average replicate and filtered m and beta vals ##
####################################################################


write.csv(beta_filt_subset, "./methylome/data/processed-data/beta-filt-subset.csv")
write.csv(m_val_filt_subset, "./methylome/data/processed-data/mval-filt-subset.csv")
write.csv(pheno_filt_subset, "./methylome/experiment-design/pheno-table-subset.csv")



