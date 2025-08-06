## Author: Andrew Palmer
## Title: Methylation Preprocessing and QC 
## Date: 14-06-2025

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

####################################################################
## Methylation preprocessing and QC
####################################################################

# This script is the first step in the methylome analysis workflow.
# It takes raw IDAT files and an meth array sheet.
# loads and performs methylation preprocessing and QC steps.
# saves normalised and filtered Beta and M-values

####################################################################
## Load required functions and packages
####################################################################

source("./methylome/scripts/util-functions/filter-functions.R")
source("./methylome/scripts/util-functions/remap-functions.R")

####################################################################
## Load data and IDATs into rg_set object
####################################################################

## data dir needs to contain IDATs and sample sheet csv file from illumina
data_dir <- "./methylome/data/processed-data/MET2022-354-014/"

## Sample Loading

## Load targets file to enable idat loading
targets <- minfi::read.metharray.sheet(data_dir)

## create red greeen set for filtering and normalisation
## extended needs to be set to TRUE to filter beads 

rg_set <- minfi::read.metharray.exp(targets = targets, extended = TRUE)

minfi::annotation(rg_set)["array"] <- "IlluminaHumanMethylationEPICv2"

minfi::annotation(rg_set)["annotation"] <- "20a1.hg38"

## set meaningful sample names for downstream plots
targets$ID <- paste(targets$Sample_Group, targets$Sample_Name, sep = ".")
minfi::sampleNames(rg_set) <- targets$ID


####################################################################
## Preprocessing
####################################################################

m_set <- minfi::preprocessRaw(rg_set)

meth_inten <- getMeth(m_set)
unmeth_inten <- getUnmeth(m_set)
detPValues <- detectionP(m_set)                          

unmeth_df <- as.data.frame(unmeth_inten)
meth_df <- as.data.frame(meth_inten)
detP_df <- as.data.frame(detPValues)
combined_df <- merge(unmeth_df, meth_df, by = "row.names", suffixes = c("_Unmeth", "_Meth"))
rownames(combined_df) <- combined_df$Row.names
combined_df <- combined_df[,-1]
colnames(detP_df)[-1] <- paste0(colnames(detP_df)[-1], "_detP")
combined_df <- merge(combined_df, detP_df, by = "row.names")
rownames(combined_df) <- combined_df$Row.names
combined_df <- combined_df[,-1]
df_reordered <- combined_df[, order(names(combined_df), decreasing = TRUE)]

write.table(df_reordered, file = "./methylome/data/processed-data/full_meth_unmeth_detp.tsv", sep = "\t", col.names = NA)
## filter out probes

m_set <- filter.EPICv2(m_set, rg_set, Sex = TRUE, DetP = TRUE, lBead = TRUE, bead.number = 3, NonCpg = TRUE, flagged = TRUE)

## normalised probes using dasen normalistaion in the wateRmelon package

m_set_dasen <- wateRmelon::dasen(m_set)
g_set_dasen <- mapToGenome(m_set_dasen, mergeManifest = TRUE)

## remove SNP probes
g_set_dasen <- dropLociWithSnps(g_set_dasen)

## extract m values
m_norm <- getM(g_set_dasen)

## remap off target probes and average replicated probes based off
## Peters 2024 https://link.springer.com/article/10.1186/s12864-024-10027-5
## code modified from the DMRcate package to enable compatability with DMRcate 2.
## this enables a custom cpg.annotate object.

m_final <- remap.Offtarget(m_norm)

m_collapsed <- collapse.Probes(m_final)

## convert M values to beta values

beta <- 2 ^ m_collapsed/(2 ^ m_collapsed + 1)


####################################################################
## Save preprocessed and filtered Beta and M values
####################################################################

write.table(beta, file = "./methylome/data/processed-data/full-dasen-normalised-filtered-beta.csv", sep = ",", col.names = NA)

write.table(m_collapsed, file = "./methylome/data/processed-data/full-dasen-normalised-filtered-m.csv", sep = ",",
            col.names = NA)


writeLines(capture.output(sessionInfo()),
           "./methylome/results/session-info/methylation-pre-processing-sessionInfo.txt")
