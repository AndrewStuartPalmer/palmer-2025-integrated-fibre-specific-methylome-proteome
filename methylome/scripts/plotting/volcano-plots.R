## Author: Andrew Palmer
## Title: Volcano Plot Script
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Volcano Plot Script
#####################################################################

## takes a list of DMPs and generates a volcano plot

#####################################################################
## Load Packages
#####################################################################

load_packages <- c("data.table",
                   "ggplot2",
                   "RColorBrewer",
                   "ggplotify",
                   "ggrepel")

lapply(load_packages, library, character.only = TRUE)


#####################################################################
## Load DMPs
#####################################################################

DMPs <- data.frame(fread("./methylome/results/final-dasen-dmps-T1preVsT2pre.csv"),
                    row.names=1)

#####################################################################
## Prepare DMP data for plotting
#####################################################################

## extract logfc and pvalues for plotting

logfc <- DMPs$logFC_beta *100
pval <- DMPs$p.value
pval <- -log10(pval)
vol_df_2 <- cbind(logfc, pval)
vol_df_2 <- as.data.frame(vol_df_2)
vol_df_3 <- vol_df_2
vol_df_3$names <- DMPs$GeneName
down_reg <- subset(vol_df_2, logfc < 0 & pval > 4.1)
down_reg_names <- subset(vol_df_3, logfc < -10 & pval > 4.1)
up_reg <- subset(vol_df_2, logfc > 0 & pval > 6)
up_reg_names <- subset(vol_df_3, logfc > 10 & pval > 4.1)
down_reg_lab <- (down_reg_names[1:10,])
up_reg_lab <- (up_reg_names[1:10,])
blues <- brewer.pal(9, name = "Blues")
labels <- rbind(up_reg_lab, down_reg_lab)

vol_df_2 <- head(vol_df_2,70000)

blues <- brewer.pal(9, name = "Blues")

labels[9,3] <- "TLN2"
labels[13,3] <- "TNNI1"
labels[14,3] <- "NA"

#####################################################################
## Generate Volcano Plot
#####################################################################

df <- data.frame(
  x = c(0),
  y = c(45),
  text = c("TI vs TII"))

vol_plot <- ggplot(vol_df_2, aes(x = logfc,
                                 y = pval)) +
    geom_point(data = vol_df_2,
             shape = 21,
             size = 2, 
             fill = blues[4], 
             colour = blues[4])+
    geom_point(data = up_reg_names,
             shape = 21,
             size = 2, 
             fill = blues[8], 
             colour = blues[8])+
    geom_point(data = down_reg_names,
             shape = 21,
             size = 2, 
             fill = blues[6], 
             colour = blues[6]) +
    
    geom_hline(yintercept = -log10(8.16290791407026E-05),
             linetype = "dashed", col = "red") + 
    geom_vline(xintercept = c(-10, 10),
               linetype = "dashed", col = "red") +
    labs(x = "Beta Difference %", y = "-log10(Pval)") +
    scale_x_continuous(breaks=c(-50,-40,-30,-20,-10, 0, 10,20,30,40,50)) +
    

    theme_bw() + # Select theme with a white background  
    theme(panel.border = element_blank(),    
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey",
                                          size = 0.75,
                                          linetype = 2),
          axis.title.x = element_text(size = 20,
                                      color = "black",
                                      family = "Helvetica"),
          axis.title.y = element_text(size = 20,
                                      color = "black",
                                      family = "Helvetica"),
          axis.text.x = element_text(color = "black",
                                     family = "Helvetica",
                                     size = 20),
          axis.text.y = element_text(color = "black",
                                     family = "Helvetica",
                                     size = 20))+
    
    geom_label_repel(data = labels,
                     aes(label = names),
                     size = 8,
                     min.segment.length = 0, 
                     seed = 42, 
                     box.padding = 0.5,
                     max.overlaps = Inf,
                     color = "grey50"
                     ) +
    geom_label(data=df, aes(x=x, y=y, label=text),                 , 
              color="black",
              family = "Helvetica",
           size=12)

#####################################################################
## Save Plot
#####################################################################

ggsave(
  "./methylome/figures/volcano-t1-pre-t2pre-dmps.pdf",
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
