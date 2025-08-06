## Author: Andrew Palmer
## Title: Generate PCA plots of proteomic data set
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025
# See the 'LICENSE' file in the root directory for full license details.

####################################################################
## PCA analysis
####################################################################

## This script is the fourth step in the proteomic analysis workflow.
## This code loads normalised and imputed protein data
## It then performs PCA analysis using prcomp.
## Producing a PCA and PCA rotation plot.

source("./proteome/scripts/proteomics-plot-functions.R")

load_packages <- c("viridis", 
                   "RColorBrewer",
                   "ggplot2",
                   "ggrepel")

lapply(load_packages, library, character.only = TRUE)

####################################################################
## Load and filter data and experimental design file
## ##################################################################

data_raw <- read.delim("./proteome/data/processed-data/P22_0579_exp2_combined_protein.tsv", header = TRUE, sep = "\t")

data_normalised_imputed <- read.delim("./proteome/data/processed-data/data_normalised_imputed.tsv", header = TRUE, sep = "\t", row.names = 1)

experiment_design <- read.delim("./proteome/experiment-design/lfq_analyst_experimental_design.txt")

##reorder experiment_design

experiment_design <- experiment_design[order(experiment_design$label), ]
colnames(data_normalised_imputed) <- experiment_design$sample
data_normalised_imputed <- data_normalised_imputed[,-51]

subset_position <- !endsWith(colnames(data_normalised_imputed), "12wk-T1")
data_normalised_imputed <- subset(data_normalised_imputed,
                                  select = subset_position)
subset_position <- !endsWith(colnames(data_normalised_imputed), "12wk-T2")
data_normalised_imputed <- subset(data_normalised_imputed,
                                  select = subset_position)
subset_position <- !endsWith(colnames(data_normalised_imputed), "12wk")
data_normalised_imputed <- subset(data_normalised_imputed,
                                  select = subset_position)

subset_position <- !endsWith(experiment_design$sample, "12wk-T1")
experiment_design <- subset(experiment_design,
                            subset = subset_position)
subset_position <- !endsWith(experiment_design$sample, "12wk-T2")
experiment_design <- subset(experiment_design,
                            subset = subset_position)
subset_position <- !endsWith(experiment_design$sample, "12wk")
experiment_design <- subset(experiment_design,
                                  subset = subset_position)

data_normalised_imputed <- data_normalised_imputed[,-2]
experiment_design<-experiment_design[-2,]


####################################################################
## Perform PCA analysis using prcomp
####################################################################

my_pca <- prcomp(t(data_normalised_imputed), scale = FALSE, center = TRUE)

## bind PC 1 and PC 2 to tranposed data frame, each row is a sample
data_normalised_imputed_pca <- cbind(t(data_normalised_imputed), my_pca$x[, 1:2])

## coerce to a data frame data_normalised_imputed_pca
data_normalised_imputed_pca <- as.data.frame(data_normalised_imputed_pca)

## add condition and sample from the experimental design table
data_normalised_imputed_pca$condition <- experiment_design$fibre_type
data_normalised_imputed_pca$sample <- experiment_design$sample
data_normalised_imputed_pca$sex <- experiment_design$sex

## generate variance explained by each principal component percent <-
percent <- round(100 * my_pca$sdev^2 / sum(my_pca$sdev^2), 1)

####################################################################
## Calculate PCA rotations
####################################################################

rotations <- as.data.frame(my_pca$rotation)
rotations <- rotations * 10
rotations_sorted <- rotations[order(rotations$PC1, rotations$PC2), ]
top_10 <- head(rotations_sorted, 10)
bottom_10 <- tail(rotations_sorted,10)


####################################################################
## Set plot colours
####################################################################

cols <- viridis_colours[c(1,5,3)]
shapes <- c(17, 16)
cols <- c("#440154FF",  "#7AD151FF",  "#2A788EFF",  "#FE9F6DFF")
reds <- brewer.pal(9, name = "Reds")
blues <- brewer.pal(9, name = "Blues")
cols <- c(reds[c(8,2)], "#808080")

####################################################################
## Generate PCA Plot
####################################################################

pca_plot <- ggplot(data_normalised_imputed_pca, aes(x = PC1,
                                                    y = PC2)) +
    geom_point(aes(color = condition, shape = sex), size = 7) +
    scale_color_manual(values = c("T1" = reds[8], "T2" = reds[2], "T2x" = "grey"),labels=c('TI', 'TII', 'TIIx'))+
    scale_shape_discrete(labels = c("F", "M"))+
    labs(x = paste0("PC", 1, ": ", percent[1], "%"), y = paste0("PC", 2, ": ", percent[2], "%"), title="Proteomics PCA plot") +
    
geom_hline(yintercept = 0,
             linetype = "dashed", col = "grey", linewidth = 1.5) + 
    geom_vline(xintercept = 0,
               linetype = "dashed", col = "grey", linewidth = 1.5) +
     
    ylim(c(-30, 30))+
     
    xlim(c(-30, 30)) +
    
    theme_bw() + # Select theme with a white background  
    theme(panel.border = element_blank(),    
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
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
          legend.text = element_text(size=40),
          legend.key.size = unit(1, 'cm'),
          legend.margin = margin(0, 0, 0, 0), # turned off for alignment
          legend.justification.bottom = "right",
          legend.justification.inside = c(1,0),
          legend.position = "inside",
             plot.margin = margin(t = 25,  # Top margin
                             r = 25,  # Right margin
                             b = 25,  # Bottom margin
                             l = 25))
  

ggsave(
  "./proteome/figures/pca-prots.pdf",
  plot = pca_plot,
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

####################################################################
## Generate PCA Rotations Plot
####################################################################

labels <- rbind(top_10, bottom_10)
labels <- labels[,1:2]
labels$labels <- rownames(labels)

pca_rotation_plot <- ggplot(rotations, aes(x = PC1,y = PC2)) +
      geom_point(data = rotations,
             shape = 21,
             size = 3, 
             fill = "grey", 
             colour = "grey")+
    geom_point(data = top_10,
             shape = 21,
             size = 3, 
             fill = blues[8], 
             colour = blues[8])+
    geom_point(data = bottom_10,
             shape = 21,
             size = 3, 
             fill = blues[6], 
             colour = blues[6]) +
 
 
    labs(x = ("PC1 (1e-1)"), y = "PC2 (1e-1)", title="Proteomics Rotation plot") +
    
geom_hline(yintercept = 0,
             linetype = "dashed", col = "grey", linewidth = 1.5) + 
    geom_vline(xintercept = 0,
               linetype = "dashed", col = "grey", linewidth = 1.5) +
     
    ylim(c(-2, 2))+
     
    xlim(c(-2, 2)) +
    
    theme_bw() + # Select theme with a white background  
    theme(panel.border = element_blank(),    
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
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
  
  geom_label_repel(data = labels,
                     aes(label = labels),
                     size = 7,
                     min.segment.length = 0, 
                     seed = 42, 
                     box.padding = 0.5,
                     max.overlaps = Inf,
                     color = "grey50"
                     )
ggsave(
  "./proteome/figures/pca-rotation-prots.pdf",
  plot = pca_rotation_plot,
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
