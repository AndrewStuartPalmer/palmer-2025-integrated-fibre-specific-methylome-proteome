## Author: Andrew Palmer
## Title: Heatmap Plots
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Script to Generate Heatmap plot
#####################################################################

library(data.table)

## Heatmaps
dmps <- data.frame(fread("./results/dasen-dmps-T1preVsT2pre.csv"), row.names = 1)

betas <- data.frame(fread("../methylome/data/processed-data/dasen-normalised-filtered-beta.csv"),                    row.names=1)

beta_rep <- betas[,c(17,18,19,90,91)]
betas$WM.f28.preWM.1 <- rowMeans(beta_rep[,4:5])
betas$WM.s201.preWM.1 <- rowMeans(beta_rep[,1:3])

beta_filt <- betas[,- c(1,2,17,18,30,57,58,90)]

pheno <- read.delim("../methylome/experiment-design/pheno-table.txt")
pheno_filt <- pheno[-c(1,2,17,18,30,57,58,90),] 


subset_position <- !endsWith(colnames(beta_filt), "4wWM")

pheno_filt_subset <- subset(pheno_filt, subset = subset_position)
beta_filt_subset <- subset(beta_filt, select = subset_position)

## take out post samples

subset_position <- grepl("pre", colnames(beta_filt))

pheno_filt_subset <- subset(pheno_filt, subset = subset_position)
beta_filt_subset <- subset(beta_filt, select = subset_position)

## keep only the top 20 DMPs
dmps$rank <- dmps$adj.p.val * dmps$logFC_beta
dmps_up <- dmps[dmps$logFC_beta > 0 & dmps$adj.p.val < 0.001,]
dmps_dn <- dmps[dmps$logFC_beta < 0 & dmps$adj.p.val < 0.001,]

dmps_up_highestlogFC <- dmps_up[order(dmps_up$rank, decreasing = FALSE),]
dmps_dn_highestlogFC <- dmps_dn[order(dmps_dn$rank, decreasing = TRUE),]

dmps_top50_hyper <- dmps_up_highestlogFC[1:50,]
dmps_top50_hypo <- dmps_dn_highestlogFC[1:50,]

heatmap_dmps <- c(rownames(dmps_top50_hypo), rownames(dmps_top50_hyper))
keep <- rownames(beta_filt_subset) %in% heatmap_dmps

heatmap_subset <- beta_filt_subset[keep,]

min.max.normalize <- function(x, na.rm = TRUE) {
    return((x- min(x)) /(max(x)-min(x)))
}

top100 <- rbind(dmps_top50_hyper, dmps_top50_hypo)

top_100_subset <- beta_filt_subset[rownames(top100),]
top_100_subset_scaled=t(scale(t(top_100_subset)))
top_100_normalised <- min.max.normalize(beta_filt_subset)
heatmap_subset_scaled=t(scale(t(heatmap_subset)))
heatmap_subset_normalised <- min.max.normalize(heatmap_subset)
genenames <- top100$GeneName
cpg <- rownames(top100)
genenames <- as.data.frame(cbind(genenames, cpg))
genenames <- genenames[match(rownames(heatmap_subset_normalised), genenames$cpg),]

my_sample_col <- as.data.frame(pheno_filt_subset$cell_type)
rownames(my_sample_col) <- colnames(heatmap_subset_normalised)
colnames(my_sample_col) <- c("Group")

blues <-  brewer.pal(8, name = "Blues")
blues_ramp <- colorRampPalette(blues)
reds <- brewer.pal(8, name = "Reds")
reds_three <- reds[c(8,2,4)]
mat_colors <- list(group = brewer.pal(3, "Set1"))
ann_cols <- list(Group = c('T1' =  "#99000D",
                           'T2' =  "#FEE0D2",
                           'WM' =  "#FC9272"))

                            

pdf("./figures/heatmap_top_dmps_test.pdf", height = 25, width = 20)
pheatmap(heatmap_subset_normalised, col = blues_ramp(30),
         main = "Top 100 DMPs TI vs TII",
         scale = "none",
         annotation_col = my_sample_col,
         annotation_colors = ann_cols,
         cluster_rows = TRUE,
         treeheight_row = 0,
         border_color = NA,
         show_colnames = FALSE,
         show_rownames = FALSE,
         labels_row = genenames$genenames,
         fontsize=18)
dev.off()
