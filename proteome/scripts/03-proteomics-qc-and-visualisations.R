## Author: Andrew Palmer
## Title: QC and Visualisations
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

####################################################################
## This code generates QC graphs for the proteomic analysis
####################################################################

####################################################################
## Load plot functions and external packages
####################################################################

source("./proteome/scripts/proteomics-plot-functions.R")

## Load ggplot2 ggrepel for generating plots
load_packages <- c("ggplot2", "ggrepel")

lapply(load_packages, library, character.only = TRUE)

####################################################################
## Load data and experimental design file
## ##################################################################

data_raw <- read.delim("./proteome/data/processed-data/P22_0579_exp2_combined_protein.tsv",
                        header = TRUE, 
                        sep = "\t")

data_filtered <- read.delim("./proteome/data/processed-data/filtered_combined_protein.tsv",                           header = TRUE,
                        sep = "\t",
                        row.names = 1)

experiment_design <- read.delim("./proteome/experiment-design/lfq_analyst_experimental_design.txt")

data_normalised_imputed <- read.delim("./proteome/data/processed-data/data_normalised_imputed.tsv",
                                      header = TRUE, 
                                      sep = "\t", 
                                      row.names = 1)

data_normalised_imputed <- data_normalised_imputed[,-51]

#####################################################################
## Proteins per Sample
## ##################################################################

pdf("./proteome/figures/protein_count_barplot_all.pdf", height = 14,width = 12)
par(oma=c(0,2,0,4), cex= 1.75, cex.main=1.75, cex.lab=1.75, cex.axis=1.75)
plot.protein.per.sample.grouped(data_filtered, 
                                experiment_design, 
                                cols = c(reds[c(8,2)],"grey"), "pre")
dev.off()

#####################################################################
## Missing value heatmap
## ##################################################################

pdf("./proteome/figures/missing_value_heatmap_prot.pdf", height = 11, width = 12)
par(oma = c(2,2,2,1), mar = c(2,2,6,2),  cex.main = 2)
heatmap.missingVal(data_filtered, 
                   experiment_design)
dev.off()

#####################################################################
## Sample Correlation Plot
## ##################################################################

data_cor <- as.data.frame(data_normalised_imputed)

cor_dat <- cor(data_cor)

experiment_design <- experiment_design[order(experiment_design$label),]

colnames(cor_dat) <- experiment_design$sample

rownames(cor_dat) <- experiment_design$sample

my_sample_col <- as.data.frame(experiment_design$fibre_type)

rownames(my_sample_col) <- experiment_design$sample

colnames(my_sample_col) <- c("condition")

pdf("./proteome/figures/correlation_heatmap_prot.pdf", height = 16, width = 16)
heatmap(cor_dat, col = (blues_ramp(30)),
         main = "Pearsons Correlation Heatmap",
         scale = "none",
         cellwidth = 15,
         cellheight = 15,
         annotation_col = my_sample_col,
         annotation_colors = list(condition =
                                    c(T2 = reds[2],
                                      T1 = reds[8],
                                      T2x = "grey")),
         cluster_rows = TRUE,
         treeheight_row = 0,
         border_color = NA,
         fontsize=18)
dev.off()

####################################################################
## Sample CVs
## ##################################################################

cv_data_combined <- sample.coefficient.variation(data_filtered,
                                                 experiment_design,
                                                 plot = FALSE)
medians <- sapply(cv_data_combined, median, na.rm = TRUE)

pdf("./proteome/figures/coefficient_of_variation_prot.pdf",
    width = 10,
    height = 8,
    pointsize = 12)
par(mfrow = c(1, 4),
    mar = c(3, 3, 3, 3),
    oma = c(2, 2, 2, 2))

hist(cv_data_combined$combined,
     breaks = seq(0, 250, length.out = 20),
     freq = TRUE,
     col = c("#44015480"),
     cex = 1.5,
     ylim = c(0, 700),
     xlim = c(0,350),
     cex.main = 1.5,
     cex.sub = 1.25,
     cex.lab = 1.25,
     cex.axis = 1.25,
     main = "Combined",
     xlab = NULL,
     ylab = NULL)
grid(nx = 4, ny = 4,
     col = "lightgray",
     lty = "dotted",
     lwd = 2,
     equilogs = TRUE)
abline(v = 20,
       col = "red",
       lty = "dashed",
       lwd = 2)
text(x = 200, y = 370,
     "Median = 27%",
     col =  "black",
     cex = 1.4)


hist(cv_data_combined$T1,
     breaks = seq(0, 275, length.out = 20),
     freq = TRUE,
     col = c("#3B528B80"),
     ylim = c(0, 700),
     xlim = c(0,350),
     cex.main = 1.5,
     cex.sub = 1.25,
     cex.lab = 1.25,
     cex.axis = 1.25,
     main = "T1",
     xlab = NULL,
     ylab = NULL)
grid(nx = 4, ny = 4, col = "lightgray", lty = "dotted",
     lwd = 2, equilogs = TRUE)
abline(v = 20, col = "red", lty = "dashed", lwd = 2)
text(x = 200, y = 370, "Median = 23%", col = "black", cex = 1.4)

hist(cv_data_combined$T2,
     breaks = seq(0, 260, length.out = 20),
     freq = TRUE,
     ylim = c(0, 700),
     xlim = c(0,350),
     col = c("#21908C80"),
     cex.main = 1.5,
     cex.sub = 1.25,
     cex.lab = 1.25,
     cex.axis = 1.25,
     main = "T2",
     xlab = NULL,
     ylab = NULL)
grid(nx = 4, ny = 4, col = "lightgray", lty = "dotted",
     lwd = 2, equilogs = TRUE)
abline(v = 20, col = "red", lty = "dashed", lwd = 2)
text(x = 200, y = 370, "Median = 24%", col =  "black", cex = 1.4)


hist(cv_data_combined$T2x,
     breaks = seq(0, 250, length.out = 20),
     freq = TRUE,
     col = c("#FDE72580"),
     cex = 1.5,
     ylim = c(0, 700),
     xlim = c(0,350),
     cex.main = 1.5,
     cex.sub = 1.25,
     cex.lab = 1.25,
     cex.axis = 1.25,
     main = "T2X",
     xlab = NULL,
     ylab = NULL)
grid(nx = 4, ny = 4, col = "lightgray", lty = "dotted",
     lwd = 2, equilogs = TRUE)
abline(v = 20, col = "red", lty = "dashed", lwd = 2)
text(x = 200, y = 370, "Median = 21%", col =  "black", cex = 1.4)

title(main = "Sample Coefficient of Variation",
      sub = "Coefficient of Variation (%)",
      outer = TRUE,
      cex.main = 2,
      cex.sub = 1.5,
      line = 0.15)

dev.off()

#####################################################################
## Grouped CV plot
## ##################################################################

cv_data_combined_grouped <- sample.coefficient.variation.grouped(data_filtered,
                                                 experiment_design,
                                                 plot = FALSE)
medians_grouped <- sapply(cv_data_combined_grouped, median, na.rm = TRUE)

pdf("./proteome/figures/coefficient_of_variation_groups.pdf",
    width = 10,
    height = 16,
    pointsize = 12)

par(mfrow = c(4, 2))
par(cex = 0.6)
par(mar = c(1, 1, 1, 1), oma = c(6, 6, 6, 6))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

hist(cv_data_combined_grouped$Myh7_pooled_fibres_pre,
     breaks = seq(0, 350, length.out = 30),
     freq = TRUE,
     col = c("#44015480"),
     cex = 1.5,
     ylim = c(0, 700),
     xlim = c(0, 350),
     cex.main = 1.5,
     cex.sub = 1.25,
     cex.lab = 1.25,
     cex.axis = 1.25,
     main = "",
     xaxt = 'n',
     xlab = NULL,
     ylab = NULL)
grid(nx = 4, ny = 4,
     col = "lightgray",
     lty = "dotted",
     lwd = 2,
     equilogs = TRUE)
abline(v = 20,
       col = "red",
       lty = "dashed",
       lwd = 2)
text(x = 200, y = 400,
     "Median = 22%",
     col =  "black",
     cex = 1.4)
text(x = 315, y = 700,
     substitute(paste(bold('TI-Pre'))),
     col =  "black",
     cex = 1.4)

hist(cv_data_combined_grouped$Myh7_pooled_fibres_post,
      breaks = seq(0, 350, length.out = 30),
     freq = TRUE,
     col = c("#44015480"),
     cex = 1.5,
     ylim = c(0, 700),
     xlim = c(0, 350),
     cex.main = 1.5,
      main = "",
     cex.sub = 1.25,
     cex.lab = 1.25,
     cex.axis = 1.25,
     yaxt = 'n',
     xaxt = 'n',
        xlab = NULL,
     ylab = NULL)

grid(nx = 4, ny = 4,
     col = "lightgray",
     lty = "dotted",
     lwd = 2,
     equilogs = TRUE)
abline(v = 20,
       col = "red",
       lty = "dashed",
       lwd = 2)
text(x = 200, y = 400,
     "Median = 21%",
     col =  "black",
     cex = 1.4)
text(x = 315, y = 700,
     substitute(paste(bold('TI-Post'))),
     col =  "black",
     cex = 1.4)


hist(cv_data_combined_grouped$Myh2_pooled_fibres_pre,
    breaks = seq(0, 350, length.out = 30),
     freq = TRUE,
     col = c("#44015480"),
     cex = 1.5,
     ylim = c(0, 700),
     xlim = c(0, 350),
     cex.main = 1.5,
      main = "",
     xaxt = 'n',
     cex.sub = 1.25,
     cex.lab = 1.25,
     cex.axis = 1.25,
     xlab = NULL,
     ylab = NULL)

grid(nx = 4, ny = 4,
     col = "lightgray",
     lty = "dotted",
     lwd = 2,
     equilogs = TRUE)
abline(v = 20,
       col = "red",
       lty = "dashed",
       lwd = 2)

text(x = 200, y = 400,
     "Median = 21%",
     col =  "black",
     cex = 1.4)
text(x = 315, y = 700,
     substitute(paste(bold('TII-Pre'))),
     col =  "black",
     cex = 1.4)


hist(cv_data_combined_grouped$Myh2_pooled_fibres_post,
     breaks = seq(0, 350, length.out = 30),
     freq = TRUE,
     col = c("#FDE72580"),
     cex = 1.5,
     ylim = c(0, 700),
     xlim = c(0, 350),
     xaxt = 'n',
      main = "",
     yaxt = "n",
     cex.main = 1.5,
     cex.sub = 1.25,
     cex.lab = 1.25,
     cex.axis = 1.25,
     xlab = NULL,
     ylab = NULL)

grid(nx = 4, ny = 4,
     col = "lightgray",
     lty = "dotted",
     lwd = 2,
     equilogs = TRUE)
abline(v = 20,
       col = "red",
       lty = "dashed",
       lwd = 2)

text(x = 200, y = 400, "Median = 22%", col =  "black", cex = 1.4)
text(x = 315, y = 700,
     substitute(paste(bold('TII-Post'))),
     col =  "black",
     cex = 1.4)

hist(cv_data_combined_grouped$Myh1_pooled_fibres_pre,
       breaks = seq(0, 350, length.out = 30),
     freq = TRUE,
     ylim = c(0, 700),
     xlim = c(0,350),
     col = c("#21908C80"),
     cex.main = 1.5,
     xaxt = "n",
      main = "",
     cex.sub = 1.25,
     cex.lab = 1.25,
     cex.axis = 1.25,
     xlab = NULL,
     ylab = NULL)

grid(nx = 4, ny = 4,
     col = "lightgray",
     lty = "dotted",
     lwd = 2,
     equilogs = TRUE)
abline(v = 20,
       col = "red",
       lty = "dashed",
       lwd = 2)
text(x = 200, y = 400, "Median = 13%", col =  "black", cex = 1.4)
text(x = 315, y = 700,
     substitute(paste(bold('TIIx-Pre'))),
     col =  "black",
     cex = 1.4)

hist(cv_data_combined_grouped$Myh1_pooled_fibres_post,
     breaks = seq(0, 350, length.out = 30),
     freq = TRUE,
     col = c("#3B528B80"),
     ylim = c(0, 700),
     xlim = c(0,350),
     yaxt = "n",
     cex.main = 1.5,
      main = "",
     cex.sub = 1.25,
     cex.lab = 1.25,
     cex.axis = 1.25,
     xlab = NULL,
     ylab = NULL)

grid(nx = 4, ny = 4,
     col = "lightgray",
     lty = "dotted",
     lwd = 2,
     equilogs = TRUE)
abline(v = 20,
       col = "red",
       lty = "dashed",
       lwd = 2)

text(x = 200, y = 400, "Median = 17%", col = "black", cex = 1.4)
text(x = 315, y = 700,
     substitute(paste(bold('TIIX-Post'))),
     col =  "black",
     cex = 1.4)


hist(cv_data_combined_grouped$combined,
     breaks = seq(0, 350, length.out = 30),
     freq = TRUE,
     col = c("#44015480"),
     cex = 1.5,
     ylim = c(0, 700),
     xlim = c(0,350),
      main = "",
     cex.main = 1.5,
     cex.sub = 1.25,
     cex.lab = 1.25,
     cex.axis = 1.25,
      xlab = NULL,
     ylab = NULL)

grid(nx = 4, ny = 4,
     col = "lightgray",
     lty = "dotted",
     lwd = 2,
     equilogs = TRUE)
abline(v = 20,
       col = "red",
       lty = "dashed",
       lwd = 2)
text(x = 200, y = 400,
     "Median = 27%",
     col =  "black",
     cex = 1.4)
text(x = 315, y = 700,
     substitute(paste(bold('Combined'))),
     col =  "black",
     cex = 1.4)

title(main = "Sample Coefficient of Variation",
      sub = "Coefficient of Variation (%)",
      ylab = "Frequency",
      outer = TRUE,
      cex.main = 2.25,
      cex.sub = 2,
      cex.lab = 2,
      line = 3)
dev.off()

####################################################################
## Ranked Intensities
####################################################################

## needs log 2 intensities

new_data <- data_normalised_imputed

new_data$mean <- apply(new_data, 1, mean, na.rm = TRUE)

new_data_ranked <- new_data[order(new_data$mean, decreasing = TRUE), ]

new_data_ranked$rank <- order(new_data_ranked$mean, decreasing = TRUE)
myosin <- grep("MYH", rownames(new_data_ranked), value = TRUE)
myosin <- myosin[-c(15, 21 ,36)]
mean <- new_data_ranked[myosin, "mean"]
rank <- new_data_ranked[myosin, "rank"]

new_data_ranked$colour[!rownames(new_data_ranked) %in% myosin] <- "lightgray"
new_data_ranked$colour[rownames(new_data_ranked) %in% myosin] <- "red"

myosin_df <- data.frame(mean,rank,myosin)

## plot ranked intensities

vol_plot <- ggplot(new_data_ranked, aes(x = rank ,y = mean)) +
      geom_point(data = new_data_ranked,
             shape = 21,
             size = 3, 
             fill = "grey", 
             colour = "grey")+
    geom_point(data = myosin_df ,
             shape = 21,
             size = 3, 
             fill = "black", 
             colour = "black")+
 
 
    labs(x = ("Protein Rank"), y = "Mean Intensity") +
    
    ylim(c(-5, 10))+
     
    xlim(c(0, 1300)) +
    
    theme_bw() + # Select theme with a white background  
    theme(panel.border = element_blank(),    
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey",
                                        linewidth = 1,
                                        linetype = 2),
          plot.title = element_text(size = 40,
                                      color = "black",
                                    family = "Helvetica",
                                    hjust = 0.5,
                                    face = "bold"),
          axis.title.x = element_text(size = 30,
                                      color = "black",
                                      family = "Helvetica"),
          axis.title.y = element_text(size = 30,
                                      color = "black",
                                      family = "Helvetica"),
          axis.text.x = element_text(color = "black",
                                     family = "Helvetica",
                                     size = 30),
          axis.text.y = element_text(color = "black",
                                     family = "Helvetica",
                                     size = 30),
          legend.title=element_blank(),
          legend.text = element_text(size=30),
          legend.key.size = unit(1, 'cm'),
          legend.margin = margin(0, 0, 0, 0), # turned off for alignment
          legend.justification.bottom = "right",
          legend.justification.inside = c(1,0),
          legend.position = "inside",
             plot.margin = margin(t = 25,  # Top margin
                             r = 25,  # Right margin
                             b = 25,  # Bottom margin
                             l = 25))+
  
  geom_label_repel(data = myosin_df,
                     aes(label = myosin),
                     size = 7,
                     min.segment.length = 0, 
                     seed = 42, 
                     box.padding = 0.5,
                     max.overlaps = Inf,
                     color = "grey50"
                     )
ggsave(
  "./proteome/figures/ranked-intensity-myosin.pdf",
  plot = vol_plot,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 12,
  height = 12,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)



## calculate percentage of signal from top 10 proteins

intensity_names <- grep(".MaxLFQ", colnames(data_filtered), value = TRUE)

holding <- data_filtered[, intensity_names]
rownames(holding) <- data_filtered$Gene
holding$mean <- apply(holding, 1, mean, na.rm = TRUE)
holding_ranked <- holding[order(holding$mean, decreasing = TRUE), ]

holding_ranked$rank <- order(holding_ranked$mean, decreasing = TRUE)
top10 <- head(holding_ranked,10)
bottom10 <- tail(holding_ranked,10)

## plot ranked intensities

pdf("./proteome/figures/ranked_intensities_prot-maxlfq.pdf", width = 10, height = 10)
par(mar = c(4, 4, 4, 4))
plot(holding_ranked$rank, holding_ranked$mean,
     xlim = c(-100, 1250),
     pch = 16,
     cex.main = 1.5,
     col = "white",
     cex.sub = 1.25,
     cex.lab = 1.25,
     cex.axis = 1.25,
     main = "Ranked Intensities",
     xlab = "",
     ylab = "")
grid(col = "lightgray", lty = "dotted",
     lwd = 3, equilogs = TRUE)

points(holding_ranked$rank, holding_ranked$mean,
       col = "lightgray", pch = 16, cex = 1.5)
points(top10$rank, top10$mean, col = "black", pch = 16, cex = 1.5)
text(c(top10$rank), c(top10$mean), labels = c(rownames(top10)),
     col = "black", cex = 1.25, pos = c(4, 3, 2, 1))
mtext("Protein Rank", side = 1, cex = 1.25, line = 2.5)
mtext("Mean Intensity", side = 2, cex = 1.25, line = 2.5)
dev.off()

pdf("./proteome/figures/ranked_intensities_prot.pdf", width = 10, height = 10)
par(mar = c(4, 4, 4, 4))
plot(new_data_ranked$rank, new_data_ranked$mean,
     ylim = c(-5, 12),
     xlim = c(-100, 1500),
     pch = 16,
     cex.main = 1.75,
     col = "white",
     cex.sub = 1.75,
     cex.lab = 1.75,
     cex.axis = 1.75,
     main = "Ranked Intensities",
     xlab = "",
     ylab = "",
     frame.plot = F,
     xaxt = "n",
     yaxt = "n")
grid(col = "darkgray", lty = "dotted",
     lwd = 3.5, equilogs = TRUE)

axis(1, c(0,500,1000,1500), cex.axis = 1.75, col = NA, line = 0)
axis(2, c(-5, 0, 5, 10), cex.axis = 1.75, col = NA, line = 0)

points(new_data_ranked$rank, new_data_ranked$mean,
       col = "lightgray", pch = 16, cex = 1.5)
points(rank, mean, col = "black", pch = 16, cex = 1.5)
text(c(rank, mean), labels = c(rownames(top10)),
     col = "black", cex = 1.25, pos = c(4, 3, 2, 1))
mtext("Protein Rank", side = 1, cex = 1.75, line = 2.5)
mtext("Mean Intensity", side = 2, cex = 1.75, line = 2.5)
dev.off()


labels <- rbind(top10, bottom_10)
labels <- labels[,1:2]
labels$labels <- rownames(labels)

vol_plot <- ggplot(new_data_ranked, aes(x = rank ,y = mean)) +
      geom_point(data = new_data_ranked,
             shape = 21,
             size = 3, 
             fill = "grey", 
             colour = "grey")+
    geom_point(data = top10 ,
             shape = 21,
             size = 3, 
             fill = "black", 
             colour = "black")+
 
 
    labs(x = ("Protein Rank"), y = "Mean Intensity") +
    
    ylim(c(-5, 10))+
     
    xlim(c(0, 1300)) +
    
    theme_bw() + # Select theme with a white background  
    theme(panel.border = element_blank(),    
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey",
                                        linewidth = 1,
                                        linetype = 2),
          plot.title = element_text(size = 40,
                                      color = "black",
                                    family = "Helvetica",
                                    hjust = 0.5,
                                    face = "bold"),
          axis.title.x = element_text(size = 30,
                                      color = "black",
                                      family = "Helvetica"),
          axis.title.y = element_text(size = 30,
                                      color = "black",
                                      family = "Helvetica"),
          axis.text.x = element_text(color = "black",
                                     family = "Helvetica",
                                     size = 30),
          axis.text.y = element_text(color = "black",
                                     family = "Helvetica",
                                     size = 30),
          legend.title=element_blank(),
          legend.text = element_text(size=30),
          legend.key.size = unit(1, 'cm'),
          legend.margin = margin(0, 0, 0, 0), # turned off for alignment
          legend.justification.bottom = "right",
          legend.justification.inside = c(1,0),
          legend.position = "inside",
             plot.margin = margin(t = 25,  # Top margin
                             r = 25,  # Right margin
                             b = 25,  # Bottom margin
                             l = 25))+
  
  geom_label_repel(data = top10,
                     aes(label = rownames(top10)),
                     size = 7,
                     min.segment.length = 0, 
                     seed = 42, 
                     box.padding = 0.5,
                     max.overlaps = Inf,
                     color = "grey50"
                     )
ggsave(
  "./proteome/figures/ranked-intensity.pdf",
  plot = vol_plot,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 12,
  height = 12,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)
