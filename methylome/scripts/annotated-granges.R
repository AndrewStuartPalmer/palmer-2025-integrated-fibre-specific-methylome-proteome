## Author: Andrew Palmer
## Title: Creating an annotated Granges for EPICv2
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Script for generating an annotated Granges Object ##
#####################################################################

## Generates an annotated Granges object for the CpGs in the EPICv2

####################################################################
## Load required packages ##
####################################################################

load_packages <- c("data.table",
                   "DMRcate")

lapply(load_packages, library, character.only = TRUE)

####################################################################
## Load required data ##
####################################################################

dmps <- data.frame(fread("./methylome/results/final-dasen-dmps-T1preVsT2pre.csv"), row.names = 1)
beta_filt_subset <-  data.frame(fread("./methylome/data/processed-data/beta-filt-subset.csv"), row.names = 1)
pheno_filt_subset <-  data.frame(fread("./methylome/experiment-design/pheno-table-subset.csv"))

DMPs <- dmps
cpgs <- rownames(dmps)
selections <- c("probeID","CpG_chrm", "CpG_beg", "CpG_end", "genesUniq")

## load accompanying annotations file

anno_epic_v2_sub <-  data.frame(fread("./methylome/annotations/my-annotation-subset.csv"), row.names = 1)
my_annotation <- anno_epic_v2_sub[,selections]
my_annotation_subset <- my_annotation[match(rownames(DMPs), my_annotation$probeID), ]

dmps$CpG_chrm <- my_annotation_subset$CpG_chrm
dmps$CpG_beg<- my_annotation_subset$CpG_beg
dmps$CpG_end <- my_annotation_subset$CpG_end

my_annotation_subset <- subset(my_annotation_subset, !is.na(CpG_chrm))
my_annotation_subset <- subset(my_annotation_subset, !is.na(CpG_end))
my_annotation_subset <- subset(my_annotation_subset, !is.na(CpG_beg))
my_annotation_subset <- subset(my_annotation_subset, !is.na(probeID))

rownames(my_annotation_subset) <- my_annotation_subset$probeID

keep <- rownames(DMPs) %in% (my_annotation_subset$probeID)
DMPs <-  DMPs[keep, ]
cpgs <- rownames(my_annotation_subset)

## set stat and diff to both logFC beta and t_beta for beta files

my_annotated_basic <- GRanges(as.character(my_annotation_subset[, "CpG_chrm"], genome = "hg38"),
                              IRanges(my_annotation_subset[, "CpG_beg"],
                                      my_annotation_subset[, "CpG_end"]),
                              stat = DMPs[, "t"], #t-statistic
                              diff = DMPs[, "logFC_beta"],
                              ind.fdr = DMPs$adj.p.val,
                              is.sig = DMPs$adj.p.val < 0.001,
                              genome = "hg38") #p-value threshold

names(my_annotated_basic) <- cpgs
##my_annotated_basic <- sort(my_annotated_basic)
my_annotated_basic$gene <- DMPs$GeneName
my_annotated_basic$logFC_beta <- DMPs$logFC_beta
my_annotated_basic$probeID <- cpgs

## required for newer versions of DMRcate
## my_annotated_basic$rawpval <- DMPs$p.value

## save custom generated annotation file

write.table(my_annotated_basic, file = "./methylome/results/final-dasen-annotated-g-ranges-t1-t2-pre-m-full.csv", sep = ",",
            col.names = NA)