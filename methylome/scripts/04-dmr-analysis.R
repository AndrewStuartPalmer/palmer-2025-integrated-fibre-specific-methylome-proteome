## Author: Andrew Palmer
## Title: Methylation DMR analysis 
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

####################################################################
## DMR analysis with DMRcate
####################################################################

# This script is the fourth step in the methylome analysis workflow.
# It takes DMPs values.
# calculates differentially methylated regions using DMRcate.
# annotates DMRs and saves the results

####################################################################
## Load required functions and packages
####################################################################

source("./methylome/scripts/util-functions/dmr-analysis-functions-final.R")

load_packages <- c("data.table", 
                   "RColorBrewer",
                   "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
                   "DMRcate", 
                   "ENmix", 
                   "ExperimentHub", 
                   "pbapply",
                   "AnnotationHub")

lapply(load_packages, library, character.only = TRUE)

####################################################################
## Load required data ##
####################################################################

beta_filt_subset <-  data.frame(fread("./methylome/data/processed-data/beta-filt-subset.csv"), 
                                row.names = 1)

pheno_filt_subset <-  data.frame(fread("./methylome/experiment-design/pheno-table-subset.csv"))

####################################################################
## Load required dmps ##
####################################################################

TIprevTIIpre <- data.frame(fread("./methylome/results/final-dasen-dmps-T1preVsT2pre.csv"), 
                                row.names = 1)

dmp_list <- list(TIprevTIIpre = TIprevTIIpre)

####################################################################
## Run DMRcate ##
####################################################################

dmrcate_dmrs <- pbapply::pblapply(dmp_list, run.DMRcate, cl = 4)

dmrcate_dmrs_df <- as.data.frame(dmrcate_dmrs)

dmrcate_dmrs_3 <- dmrcate_dmrs$TIprevTIIpre[dmrcate_dmrs$TIprevTIIpre$no.cpgs >=3 ,]

dmrcate_dmrs_4 <- dmrcate_dmrs$TIprevTIIpre[dmrcate_dmrs$TIprevTIIpre$no.cpgs >=4 ,]

pbapply::pbmapply(function (x,y) write.csv(x, 
                                           file = paste0('./methylome/results/final-dmrcate-dmrs-',y, '.csv'), 
                                           row.names = T), 
                                           dmrcate_dmrs, names(dmrcate_dmrs))

####################################################################
## Add annotations using encode 47##
####################################################################

gtf <- rtracklayer::import("./methylome/annotations/human_gencode_47.gtf")

dmrcate_dmrs <- add.Annotations(dmrcate_dmrs$TIprevTIIpre, 
                                gtf, 
                                ensemble = TRUE, 
                                symbols = TRUE, 
                                geneType = TRUE, 
                                nearestGene = TRUE)


####################################################################
## Obtain a list of CpGs within each DMR ##
####################################################################

CpGs <- TIprevTIIpre[,c("Seqname","Start","End")]

CpGs$probeID <- rownames(TIprevTIIpre)

colnames(CpGs) <- c("seqnames", "start", "end", "probeID")

CpGs <- CpGs[!is.na(CpGs$seqnames), ]

CpGs <- makeGRangesFromDataFrame(CpGs, 
                                 keep.extra.columns = TRUE, 
                                 na.rm=TRUE)

dmrcate_dmrs_final <- add.Probes(dmrcate_dmrs, 
                                 CpGs)


####################################################################
## Save results ##
####################################################################

write.table(dmrcate_dmrs_final, 
            file = "./methylome/results/dmrcate-DMRs-tI-pre-vs-tII-pre.csv", 
            sep = ",", 
            row.names = FALSE)




