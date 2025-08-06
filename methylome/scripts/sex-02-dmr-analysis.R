## Author: Andrew Palmer
## Title: DMR Sex differences analysis 
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Script for identifying DMRs 
#####################################################################

####################################################################
## Load custom methylation analysis functions ##
####################################################################

source("./methylome/scripts/util-functions/dmr-analysis-functions-final.R")

load_packages <- c("data.table", "RColorBrewer",
                   "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
                   "DMRcate", "ENmix", "ExperimentHub", "pbapply", "AnnotationHub")

lapply(load_packages, library, character.only = TRUE)

####################################################################
## Load required data ##
####################################################################

beta_filt_subset <-  data.frame(fread("./methylome/data/processed-data/beta-filt-subset.csv"), row.names = 1)

pheno_filt_subset <-  data.frame(fread("./methylome/experiment-design/pheno-table-subset.csv"))

####################################################################
## Load required dmps ##
####################################################################

tI_FvM <- data.frame(fread("./methylome/results/final-dasen-dmps-T1preFVsT1preM.csv"), row.names = 1)

tII_FvM <- data.frame(fread("./methylome/results/final-dasen-dmps-T2preFVsT2preM.csv"), row.names = 1)

wm_FvM <- data.frame(fread("./methylome/results/final-dasen-dmps-WMpreFVsWMpreM.csv"), row.names = 1)

interaction <- data.frame(fread("./methylome/results/final-dasen-dmps-Interaction1.csv"), row.names = 1)

dmp_list <- list(tI_FvM = tI_FvM, tII_FvM = tII_FvM, wm_FvM = wm_FvM, interaction = interaction)

####################################################################
## Run DMRcate ##
####################################################################

dmrcate_dmrs <- pbapply::pblapply(dmp_list, run.DMRcate, cl = 4)

####################################################################
## Add annotations using encode 47##
####################################################################

CpGs <- data.frame(fread("./methylome/results/final-dasen-dmps-T1preVsT2pre.csv"), row.names = 1)
CpGs <- CpGs[,c("Seqname","Start","End")]
CpGs$probeID <- rownames(CpGs)
colnames(CpGs) <- c("seqnames", "start", "end", "probeID")
CpGs <- CpGs[!is.na(CpGs$seqnames), ]
CpGs <- makeGRangesFromDataFrame(CpGs, keep.extra.columns = TRUE, na.rm=TRUE)

#####################################################################
## Load Humand Gencode 47 for annotations 
#####################################################################

gtf <- rtracklayer::import("./methylome/annotations/human_gencode_47.gtf")

dmrs_tI_final <- add.Annotations(dmrcate_dmrs$tI_FvM, gtf, ensemble = TRUE, symbols = TRUE, geneType = TRUE, nearestGene = TRUE)
dmrcate_tI_final <- add.Probes(dmrs_tI_final, CpGs)

dmrs_tII_final <- add.Annotations(dmrcate_dmrs$tII_FvM, gtf, ensemble = TRUE, symbols = TRUE, geneType = TRUE, nearestGene = TRUE)
dmrcate_tII_final <- add.Probes(dmrs_tII_final, CpGs)

dmrs_wm_final <- add.Annotations(dmrcate_dmrs$wm_FvM, gtf, ensemble = TRUE, symbols = TRUE, geneType = TRUE, nearestGene = TRUE)
dmrcate_wm_final <- add.Probes(dmrs_wm_final, CpGs)

dmrs_interaction_final <- add.Annotations(dmrcate_dmrs$interaction, gtf, ensemble = TRUE, symbols = TRUE, geneType = TRUE, nearestGene = TRUE)
dmrcate_interaction_final <- add.Probes(dmrs_interaction_final, CpGs)

#####################################################################
##  Save results
#####################################################################

write.table(dmrcate_wm_final, file = "./methylome/results/final-dmrcate-dmrs-wm_FvM.csv", sep = ",", row.names = FALSE)


write.table(dmrcate_interaction_final, file = "./methylome/results/final-dmrcate-dmrs-interaction.csv", sep = ",", row.names = FALSE)


write.table(dmrcate_tI_final, file = "./methylome/results/final-dmrcate-dmrs-tI_FvM.csv", sep = ",", row.names = FALSE)


write.table(dmrcate_tII_final, file = "./methylome/results/final-dmrcate-dmrs-tII_FvM.csv", sep = ",", row.names = FALSE)
