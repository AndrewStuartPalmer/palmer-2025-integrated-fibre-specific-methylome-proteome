## Author: Andrew Palmer
## Title: Sex DMR plots
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Sex DMR plots
#####################################################################

## Script to generate DMR plots for the sex comparison

#####################################################################
## Load Packages
#####################################################################

load_packages <- c("data.table", "RColorBrewer",
                   "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
                   "DMRcate", "GenomicRanges", "pbapply", "ComplexHeatmap")

lapply(load_packages, library, character.only = TRUE)

blues <-  brewer.pal(8, name = "Blues")
blues_ramp <- colorRampPalette(blues)


beta_filt_subset <-  data.frame(fread("./methylome/data/processed-data/beta-filt-subset.csv"), row.names = 1)
pheno_filt_subset <-  data.frame(fread("./methylome/experiment-design/pheno-table-subset.csv"))
my_annotated_basic <-  data.frame(fread("./methylome/results/final-dasen-dmps-T1preVsT2pre.csv"))
my_annotated_basic <- my_annotated_basic[,c(1,12,13,14,2)]
colnames(my_annotated_basic) <- c("probeID","seqname","start","end", "genes")

my_annotated_basic <- my_annotated_basic[complete.cases(my_annotated_basic$start),]

probes <- my_annotated_basic$probeID

results_ranges_tI <-  data.frame(fread("./methylome/results/final-dmrcate-dmrs-tI_FvM.csv"))
results_ranges_tI <- results_ranges_tI[results_ranges_tI$no.cpgs > 3 & results_ranges_tI$Fisher < 0.001,]

results_ranges_tII <-  data.frame(fread("./methylome/results/final-dmrcate-dmrs-tII_FvM.csv"))
results_ranges_tII <- results_ranges_tII[results_ranges_tII$no.cpgs > 3 & results_ranges_tII$Fisher < 0.001,]

results_ranges_wm <-  data.frame(fread("./methylome/results/final-dmrcate-dmrs-wm_FvM.csv"))
results_ranges_wm <- results_ranges_wm[results_ranges_wm$no.cpgs > 3 & results_ranges_wm$Fisher < 0.001,]

my_annotated_basic <- makeGRangesFromDataFrame(my_annotated_basic, na.rm = TRUE)
names(my_annotated_basic) <- probes

prepareBetaVals <- function(ranges, annotation, beta_vals) {
## code to generate beta_value means of all DMRs

    genes <- ranges$symbol
    ranges <- makeGRangesFromDataFrame(ranges)
    
    allDMRcpgs <- pbsapply(1:length(ranges), function (x) names(subsetByOverlaps(annotation, ranges[x])), cl = 5)
    
names(allDMRcpgs) <- genes

testf <- function(x){
        DMR <- beta_filt_subset[rownames(beta_filt_subset) %in% x,]
        DMRi <- colMeans(DMR)
}

    dmr_data <- pbsapply(allDMRcpgs, testf, simplify = TRUE, USE.NAMES = TRUE, cl = 5)

    dmr_data <- t(dmr_data)
    subset_position <- !endsWith(colnames(dmr_data), "12wWM")
    dmr_data <- subset(dmr_data, select = subset_position)
    subset_position <- !endsWith(colnames(dmr_data), "12wT1")
    dmr_data <- subset(dmr_data, select = subset_position)
    subset_position <- !endsWith(colnames(dmr_data), "12wT2")
    dmr_data <- subset(dmr_data, select = subset_position)

    z_scores <- scale(dmr_data, scale = TRUE, center = FALSE)

    z_scores

    }


####################################################################
## Generate heatmap of DMRs TI ##
####################################################################

z_scores <- prepareBetaVals(results_ranges_tI, my_annotated_basic, beta_filt_subset)

## prepare pheno table for annotations

    pheno_small <- pheno_filt_subset[,c("samplename", "cell_type", "sex")]
    subset_position <- !endsWith(pheno_small$samplename, "12wWM")
    pheno_small_subset <- subset(pheno_small, subset = subset_position)
    subset_position <- !endsWith(pheno_small_subset$samplename, "12wT1")
    pheno_small_subset <- subset(pheno_small_subset, subset = subset_position)
    subset_position <- !endsWith(pheno_small_subset$samplename, "12wT2")
    pheno_small_subset <- subset(pheno_small_subset, subset = subset_position)
    colnames(z_scores) <- pheno_small_subset$samplename
    
safe_colorblind_palette <- hcl.colors(13, alpha = 0.8)

ha = HeatmapAnnotation(df = data.frame(SampleType = pheno_small_subset$cell_type),
                      col = list(SampleType = c('T1' =  "#99000D",
                                 'T2'=  "#FEE0D2",
                                 'WM' =  "#FC9272")),
                      annotation_legend_param = list(labels_gp = gpar(fontsize = 20),
                                                     title_gp = gpar(fontsize=20, fontface = 'bold'),
                                                     title = "Sample Type",
                                                     grid_height = unit(4, "mm"), grid_width = unit(4, "mm")),
                      simple_anno_size = unit(4, "mm"),
                      annotation_name_gp= gpar(fontsize = 0))

ha2 = HeatmapAnnotation(df = data.frame(SexType = pheno_small_subset$sex),
                        col = list(SexType = c('M' = safe_colorblind_palette[10] ,
                                               'F' =  safe_colorblind_palette[3])),
                        annotation_legend_param = list(labels_gp = gpar(fontsize = 20),
                                                       title_gp = gpar(fontsize=20, fontface = 'bold'),
                                                       title = "Sex",
                                                       grid_height = unit(4, "mm"), grid_width = unit(4, "mm")),
                        simple_anno_size = unit(4, "mm"),
                        annotation_name_gp= gpar(fontsize = 0))

ht_opt(legend_gap = unit(40, "mm"))

htmap <- ComplexHeatmap::Heatmap(z_scores,
                                 name = "Z Scores",
                                 show_column_names = FALSE,
                                 row_names_side = "left",
                                 col = blues_ramp(20),
                                 cluster_rows = TRUE,
                                 cluster_columns = TRUE,
                                 show_column_dend = TRUE,
                                 show_row_dend = FALSE,
                                 top_annotation = c(ha, ha2),
                                 height = 4,
                                 width = 4,
                                 row_km = 2, column_km = 4,
                                 row_title = "Cluster %s",
                                 row_title_gp = gpar(fontsize = 16),
                                 column_title = "Cluster %s",
                                 column_title_gp = gpar(fontsize = 16),
                                 row_km_repeats = 100,
                                 column_km_repeats = 100,                                                          row_names_max_width = unit(8, "cm"),
                                 row_names_gp = grid::gpar(fontsize = 8),
 heatmap_legend_param =
                            list(labels_gp = gpar(fontsize = 20),
                                 title_gp = gpar(fontsize=20, fontface = 'bold')
                                ,
                                 direction = "vertical",
                                 grid_height = unit(0.4, "mm"), grid_width = unit(4, "mm"))                             
                                          )

pdf("./methylome/figures/robust-dmrs-tI-sex-heatmap.pdf", height = 12, width = 12)
draw(htmap, heatmap_legend_side = "right", annotation_legend_side = "bottom", merge_legends =TRUE)
dev.off()



####################################################################
## Generate heatmap of DMRs TII ##
####################################################################

z_scores <- prepareBetaVals(results_ranges_tII, my_annotated_basic, beta_filt_subset)

## prepare pheno table for annotations

    pheno_small <- pheno_filt_subset[,c("samplename", "cell_type", "sex")]
    subset_position <- !endsWith(pheno_small$samplename, "12wWM")
    pheno_small_subset <- subset(pheno_small, subset = subset_position)
    subset_position <- !endsWith(pheno_small_subset$samplename, "12wT1")
    pheno_small_subset <- subset(pheno_small_subset, subset = subset_position)
    subset_position <- !endsWith(pheno_small_subset$samplename, "12wT2")
    pheno_small_subset <- subset(pheno_small_subset, subset = subset_position)
    colnames(z_scores) <- pheno_small_subset$samplename
    
safe_colorblind_palette <- hcl.colors(13, alpha = 0.8)

ha = HeatmapAnnotation(df = data.frame(SampleType = pheno_small_subset$cell_type),
                      col = list(SampleType = c('T1' =  "#99000D",
                                 'T2'=  "#FEE0D2",
                                 'WM' =  "#FC9272")),
                      annotation_legend_param = list(labels_gp = gpar(fontsize = 20),
                                                     title_gp = gpar(fontsize=20, fontface = 'bold'),
                                                     title = "Sample Type",
                                                     grid_height = unit(4, "mm"), grid_width = unit(4, "mm")),
                      simple_anno_size = unit(4, "mm"),
                      annotation_name_gp= gpar(fontsize = 0))

ha2 = HeatmapAnnotation(df = data.frame(SexType = pheno_small_subset$sex),
                        col = list(SexType = c('M' = safe_colorblind_palette[10] ,
                                               'F' =  safe_colorblind_palette[3])),
                        annotation_legend_param = list(labels_gp = gpar(fontsize = 20),
                                                       title_gp = gpar(fontsize=20, fontface = 'bold'),
                                                       title = "Sex",
                                                       grid_height = unit(4, "mm"), grid_width = unit(4, "mm")),
                        simple_anno_size = unit(4, "mm"),
                        annotation_name_gp= gpar(fontsize = 0))

ht_opt(legend_gap = unit(40, "mm"))

htmap <- ComplexHeatmap::Heatmap(z_scores,
                                 name = "Z Scores",
                                 show_column_names = FALSE,
                                 row_names_side = "left",
                                 col = blues_ramp(20),
                                 cluster_rows = TRUE,
                                 cluster_columns = TRUE,
                                 show_column_dend = TRUE,
                                 show_row_dend = FALSE,
                                 top_annotation = c(ha, ha2),
                                 height = 4,
                                 width = 4,
                                 row_km = 2, column_km = 4,
                                 row_title = "Cluster %s",
                                 row_title_gp = gpar(fontsize = 16),
                                 column_title = "Cluster %s",
                                 column_title_gp = gpar(fontsize = 16),
                                 row_km_repeats = 100,
                                 column_km_repeats = 100,                                                          row_names_max_width = unit(8, "cm"),
                                 row_names_gp = grid::gpar(fontsize = 6),
 heatmap_legend_param =
                            list(labels_gp = gpar(fontsize = 20),
                                 title_gp = gpar(fontsize=20, fontface = 'bold')
                                ,
                                 direction = "vertical",
                                 grid_height = unit(0.4, "mm"), grid_width = unit(4, "mm"))                             
                                          )

pdf("./methylome/figures/robust-dmrs-tII-sex-heatmap.pdf", height = 12, width = 12)
draw(htmap, heatmap_legend_side = "right", annotation_legend_side = "bottom", merge_legends =TRUE)
dev.off()

####################################################################
## Generate heatmap of DMRs WM ##
####################################################################

z_scores <- prepareBetaVals(results_ranges_wm, my_annotated_basic, beta_filt_subset)

## prepare pheno table for annotations

    pheno_small <- pheno_filt_subset[,c("samplename", "cell_type", "sex")]
    subset_position <- !endsWith(pheno_small$samplename, "12wWM")
    pheno_small_subset <- subset(pheno_small, subset = subset_position)
    subset_position <- !endsWith(pheno_small_subset$samplename, "12wT1")
    pheno_small_subset <- subset(pheno_small_subset, subset = subset_position)
    subset_position <- !endsWith(pheno_small_subset$samplename, "12wT2")
    pheno_small_subset <- subset(pheno_small_subset, subset = subset_position)
    colnames(z_scores) <- pheno_small_subset$samplename

safe_colorblind_palette <- hcl.colors(13, alpha = 0.8)

ha = HeatmapAnnotation(df = data.frame(SampleType = pheno_small_subset$cell_type),
                      col = list(SampleType = c('T1' =  "#99000D",
                                 'T2'=  "#FEE0D2",
                                 'WM' =  "#FC9272")),
                      annotation_legend_param = list(labels_gp = gpar(fontsize = 20),
                                                     title_gp = gpar(fontsize=20, fontface = 'bold'),
                                                     title = "Sample Type",
                                                     grid_height = unit(4, "mm"), grid_width = unit(4, "mm")),
                      simple_anno_size = unit(4, "mm"),
                      annotation_name_gp= gpar(fontsize = 0))

ha2 = HeatmapAnnotation(df = data.frame(SexType = pheno_small_subset$sex),
                        col = list(SexType = c('M' = safe_colorblind_palette[10] ,
                                               'F' =  safe_colorblind_palette[3])),
                        annotation_legend_param = list(labels_gp = gpar(fontsize = 20),
                                                       title_gp = gpar(fontsize=20, fontface = 'bold'),
                                                       title = "Sex",
                                                       grid_height = unit(4, "mm"), grid_width = unit(4, "mm")),
                        simple_anno_size = unit(4, "mm"),
                        annotation_name_gp= gpar(fontsize = 0))

ht_opt(legend_gap = unit(40, "mm"))

htmap <- ComplexHeatmap::Heatmap(z_scores,
                                 name = "Z Scores",
                                 show_column_names = FALSE,
                                 row_names_side = "left",
                                 col = blues_ramp(20),
                                 cluster_rows = TRUE,
                                 cluster_columns = TRUE,
                                 show_column_dend = TRUE,
                                 show_row_dend = FALSE,
                                 top_annotation = c(ha, ha2),
                                 height = 4,
                                 width = 4,
                                 row_km = 2, column_km = 4,
                                 row_title = "Cluster %s",
                                 row_title_gp = gpar(fontsize = 16),
                                 column_title = "Cluster %s",
                                 column_title_gp = gpar(fontsize = 16),
                                 row_km_repeats = 100,
                                 column_km_repeats = 100,                                                          row_names_max_width = unit(8, "cm"),
                                 row_names_gp = grid::gpar(fontsize = 6),
 heatmap_legend_param =
                            list(labels_gp = gpar(fontsize = 20),
                                 title_gp = gpar(fontsize=20, fontface = 'bold')
                                ,
                                 direction = "vertical",
                                 grid_height = unit(0.4, "mm"), grid_width = unit(4, "mm"))                             
                                          )

pdf("./methylome/figures/robust-dmrs-wm-sex-heatmap.pdf", height = 12, width = 12)
draw(htmap, heatmap_legend_side = "right", annotation_legend_side = "bottom", merge_legends =TRUE)
dev.off()

