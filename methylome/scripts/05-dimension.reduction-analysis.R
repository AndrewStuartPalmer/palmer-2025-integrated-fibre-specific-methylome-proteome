## Author: Andrew Palmer
## Title: Dimension Reduction Methylation Analysis
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Script for dimension reduction methylation analysis ##
#####################################################################

## this code performs dimension reduction analysis
## it requires normalised M values 
## and performs PCA analysis

####################################################################
## Load required packages ##
####################################################################

load_packages <- c("data.table", "RColorBrewer", "ggplot2", "ggrepel")

lapply(load_packages, library, character.only = TRUE)


## using fread from the data table package is quicker
## load m_vals
m_val <-  data.frame(fread("./methylome/data/processed-data/full-dasen-normalised-filtered-m.csv"), row.names = 1)
   
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
pheno_filt <- pheno[-c(17,18,57,58,90),] 

m_val_filt <- m_val[,- c(17,18,57,58,90)]
row.names(pheno_filt) <- 1:nrow(pheno_filt)

## remove 4wWM samples

subset_position <- !endsWith(colnames(m_val_filt), "4wWM")
m_val_filt_subset <- subset(m_val_filt, select = subset_position)
pheno_filt_subset <- subset(pheno_filt, subset = subset_position)


subset_position <- !endsWith(colnames(m_val_filt_subset), "12wWM")
m_val_filt_subset <- subset(m_val_filt_subset, select = subset_position)
pheno_filt_subset <- subset(pheno_filt_subset, subset = subset_position)


subset_position <- !endsWith(colnames(m_val_filt_subset), "12wT1")
m_val_filt_subset <- subset(m_val_filt_subset, select = subset_position)
pheno_filt_subset <- subset(pheno_filt_subset, subset = subset_position)


subset_position <- !endsWith(colnames(m_val_filt_subset), "12wT2")
m_val_filt_subset <- subset(m_val_filt_subset, select = subset_position)
pheno_filt_subset <- subset(pheno_filt_subset, subset = subset_position)

reds <- brewer.pal(8, name = "Reds")
my_cols <- reds[c(8,6,2)]
shapes <- c(16,17)

####################################################################
## Run prcomp ##
####################################################################

pca_mvals <- prcomp(t(m_val_filt_subset), center = TRUE, scale = FALSE)

PC1_and_PC2 <- data.frame(PC1=pca_mvals$x[,1], PC2= pca_mvals$x[,2], type = pheno_filt_subset$cell_type)

PC1_and_PC2$sex <- pheno_filt_subset$sex

PC1_and_PC2$position <-  pheno_filt_subset$array

PC1_and_PC2$slide <- pheno_filt_subset$slide
PC1_and_PC2$slide <- as.integer(factor((PC1_and_PC2$slide)))


PC1_and_PC2$sex[is.na(PC1_and_PC2$sex)] <- "Mix"

percent <- round(100* pca_mvals$sdev^2 / sum(pca_mvals$sdev^2), 1)
shapes <- c(16,17,15)

rotations <- as.data.frame(pca_mvals$rotation)
rotations <- rotations * 100
rotations_sorted <- rotations[order(rotations$PC1, rotations$PC2), ]
top_10 <- head(rotations_sorted, 15)
bottom_10 <- tail(rotations_sorted,15)

annotation_epicv2 <- read.delim("./methylome/annotations/annotation-epic-v2.txt")

gene_names_top_10 <- annotation_epicv2[match(rownames(top_10), annotation_epicv2$probeID),]

gene_names_bottom_10 <- annotation_epicv2[match(rownames(bottom_10), annotation_epicv2$probeID),]

top_10$gene <- gene_names_top_10$genesUniq
bottom_10$gene <- gene_names_bottom_10$genesUniq
bottom_10 <- bottom_10[complete.cases(bottom_10), ]
top_10 <- top_10[complete.cases(top_10), ]

blues <- brewer.pal(9, name = "Blues")

####################################################################
## Produce PCA plot ##
####################################################################

vol_plot <- ggplot(PC1_and_PC2, aes(x = PC1,
                                                    y = PC2)) +
    geom_point(aes(color = type, shape = sex), size = 7) +
    scale_color_manual(values = c("T1" = reds[8], "T2" = reds[2], "WM" = reds[4],"T2X" = "grey"),labels=c('TI', 'TII', 'TIIx', 'WM'))+
    scale_shape_discrete(labels = c("Female", "Male", "Mix"))+
    labs(x = paste0("PC", 1, ": ", percent[1], "%"), , y = paste0("PC", 2, ": ", percent[2], "%"), title="DNAm PCA plot") +
    
geom_hline(yintercept = 0,
             linetype = "dashed", col = "grey", size = 1.5) + 
    geom_vline(xintercept = 0,
               linetype = "dashed", col = "grey", size = 1.5) +
     
    ylim(c(-200, 200))+
     
    xlim(c(-200, 200)) +
    
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
                               l = 25))
  

ggsave(
  "./methylome/figures/pca-dnam.pdf",
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


####################################################################
## Produce PCA rotation plot ##
####################################################################

labels <- rbind(top_10, bottom_10)
labels$gene[c(4,7)] <- "ATP2A2"
labels$gene[17] <- "MYH1;MYHAS" 

vol_plot <- ggplot(rotations, aes(x = PC1,y = PC2)) +
      geom_point(data = head(rotations,10000),
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
 
 
    labs(x = ("PC1 (1e-1)"), , y = "PC2 (1e-1)", title="DNAm Rotation plot") +
    
geom_hline(yintercept = 0,
             linetype = "dashed", col = "grey", size = 1.5) + 
    geom_vline(xintercept = 0,
               linetype = "dashed", col = "grey", size = 1.5) +
     
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
                     aes(label = gene),
                     size = 7,
                     min.segment.length = 0, 
                     seed = 42, 
                     box.padding = 0.5,
                     max.overlaps = Inf,
                     color = "grey50"
                     )
ggsave(
  "./methylome/figures/pca-rotation-dnam.pdf",
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
