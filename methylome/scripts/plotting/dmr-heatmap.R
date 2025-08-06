## Author: Andrew Palmer
## Title: DMR heatmap
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Script for the generation of DMR heatmaps
#####################################################################

#####################################################################
## Load Packages
#####################################################################

load_packages <- c("data.table", "RColorBrewer",
                   "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
                   "DMRcate", "GenomicRanges", "pbapply", "ComplexHeatmap")

lapply(load_packages, library, character.only = TRUE)


#####################################################################
## Load required data
#####################################################################

beta_filt_subset <-  data.frame(fread("./methylome/data/processed-data/beta-filt-subset.csv"), row.names = 1)
pheno_filt_subset <-  data.frame(fread("./methylome/experiment-design/pheno-table-subset.csv"))
my_annotated_basic <-  data.frame(fread("./methylome/results/final-dasen-dmps-T1preVsT2pre.csv"))
my_annotated_basic <- my_annotated_basic[,c(1,12,13,14,2)]
colnames(my_annotated_basic) <- c("probeID","seqname","start","end", "genes")


probes <- my_annotated_basic$probeID
results.ranges <-  data.frame(fread("./methylome/results/dmrcate-DMRs-tI-pre-vs-tII-pre.csv"))
results.ranges <- results.ranges[results.ranges$no.cpgs >3,]
my_annotated_basic <- makeGRangesFromDataFrame(my_annotated_basic)
names(my_annotated_basic) <- probes
results.ranges <- results.ranges[results.ranges$meandiff > 0.1 | results.ranges$meandiff < -0.1,]

genes <- results.ranges$symbol
mean_diff <- results.ranges$meandiff
results.ranges <- makeGRangesFromDataFrame(results.ranges)


#####################################################################
## Generate Mean beta values of all DMRs
#####################################################################

## code to generate beta_value means of all DMRs

allDMRcpgs <- pbsapply(1:length(results.ranges), function (x) names(subsetByOverlaps(my_annotated_basic, results.ranges[x])), cl = 5)
names(allDMRcpgs) <- genes

testf <- function(x){
        DMR <- beta_filt_subset[rownames(beta_filt_subset) %in% x,]
        DMRi <- colMeans(DMR)
}

dmr_data <- pbsapply(allDMRcpgs, testf, simplify = TRUE, USE.NAMES = TRUE, cl = 5)

dmr_data <- t(dmr_data)
subset_position <- !endsWith(colnames(dmr_data), "4wWM")
dmr_data <- subset(dmr_data, select = subset_position)
subset_position <- !endsWith(colnames(dmr_data), "12wWM")
dmr_data <- subset(dmr_data, select = subset_position)
subset_position <- !endsWith(colnames(dmr_data), "12wT1")
dmr_data <- subset(dmr_data, select = subset_position)
subset_position <- !endsWith(colnames(dmr_data), "12wT2")
dmr_data <- subset(dmr_data, select = subset_position)

z_scores <- scale(dmr_data, scale = TRUE, center = FALSE)

pheno_small <- pheno_filt_subset[,c("samplename", "cell_type")]
subset_position <- !endsWith(rownames(pheno_small), "4wWM")
pheno_small_subset <- subset(pheno_small, subset = subset_position)
subset_position <- !endsWith(pheno_small_subset$samplename, "12wWM")
pheno_small_subset <- subset(pheno_small_subset, subset = subset_position)
subset_position <- !endsWith(pheno_small_subset$samplename, "12wT1")
pheno_small_subset <- subset(pheno_small_subset, subset = subset_position)
subset_position <- !endsWith(pheno_small_subset$samplename, "12wT2")
pheno_small_subset <- subset(pheno_small_subset, subset = subset_position)

####################################################################
## Generate heatmap of DMRs ##
####################################################################

my_sample_col <- as.data.frame(pheno_small_subset$cell_type)
my_sample_col$Sample <- pheno_small_subset$samplename
colnames(my_sample_col) <- c("Group")
blues <-  brewer.pal(8, name = "Blues")
blues_ramp <- colorRampPalette(blues)
reds <- brewer.pal(8, name = "Reds")
reds_three <- reds[c(8,2,4)]
mat_colors <- list(Group = brewer.pal(3, "Set1"))
ann_cols <- list(Group = c('T1' =  "#99000D",
                           'T2' =  "#FEE0D2",
                           'WM' =  "#FC9272"))

colnames(z_scores) <- pheno_small_subset$samplename

ha = HeatmapAnnotation(df = data.frame(SampleType = pheno_small_subset$cell_type),
                      col = list(SampleType = c('T1' =  "#99000D",
                                 'T2' =  "#FEE0D2",
                                 'WM' =  "#FC9272")),
                      annotation_legend_param = list(labels_gp = gpar(fontsize = 16),
                                                     title_gp = gpar(fontsize=16, fontface = 'bold'),
                                                     title = "Sample Type",
                                                     grid_height = unit(4, "mm"), grid_width = unit(4, "mm")),
                      simple_anno_size = unit(4, "mm"),
                        annotation_name_gp= gpar(fontsize = 0))

ht_opt(legend_gap = unit(40, "mm"))

htmap <- ComplexHeatmap::Heatmap(dmr_data,
                                 name = "Mean methylation",
                                 show_column_names = FALSE,
                                 row_names_side = "left",
                                 col = blues_ramp(20),
                                 cluster_rows = TRUE,
                                 cluster_columns = TRUE,
                                 show_column_dend = TRUE,
                                 show_row_dend = FALSE,
                                 top_annotation = c(ha),
                                 height = 4,
                                 width = 4,
                                 row_km = 2, column_km = 2,
                                 row_title = "Cluster %s",
                                 row_title_gp = gpar(fontsize = 16),
                                 column_title = "Cluster %s",
                                 column_title_gp = gpar(fontsize = 16),
                                 row_km_repeats = 100,
                                 column_km_repeats = 100,
                                 row_names_max_width = unit(8, "cm"),
                                 row_names_gp = grid::gpar(fontsize = 5),
 heatmap_legend_param =
                            list(labels_gp = gpar(fontsize = 16),
                                 title_gp = gpar(fontsize=16, fontface = 'bold')
                                ,
                                 direction = "vertical",
                                 grid_height = unit(0.4, "mm"), grid_width = unit(4, "mm"))                             
                                          )


pdf("./methylome/figures/robust-dmr-heatmap-final.pdf", height = 18, width = 12)
draw(htmap, heatmap_legend_side = "right", annotation_legend_side = "bottom", merge_legends =TRUE)
dev.off()


