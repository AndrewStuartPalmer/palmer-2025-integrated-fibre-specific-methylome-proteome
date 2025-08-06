## Author: Andrew Palmer
## Title: Sex Volcano Plots
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Script to generate Volcano plots for Sex DMPs
#####################################################################

#####################################################################
## Load Packages
#####################################################################

load_packages <- c("data.table",
                   "ggplot2",
                   "RColorBrewer",
                   "ggplotify",
                   "ggrepel")

lapply(load_packages, library, character.only = TRUE)

## load DMPs

DMPs <- data.frame(fread("./methylome/results/final-dasen-dmps-T1preFVsT1preM.csv"),
                    row.names=1)

## volcano pot p val TI-pre TI-post

logfc <- DMPs$logFC_beta * 100
pval <- DMPs$p.value
pval <- -log10(pval)
vol_df_2 <- cbind(logfc, pval)
vol_df_2 <- as.data.frame(vol_df_2)
vol_df_3 <- vol_df_2
vol_df_3$names <- DMPs$GeneName
vol_df_2 <- head(vol_df_2,10000)
non_sig <- subset(vol_df_2,  pval < 6.21)
sig <- subset(vol_df_2, logfc < 10 & pval > 6.21 | logfc > -10 & pval > 6.21)
down_reg <- subset(vol_df_2, logfc < 0 & pval > 6.21)
down_reg_names <- subset(vol_df_3, logfc < -10 & pval > 6.21)
up_reg <- subset(vol_df_2, logfc > 0 & pval > 6.21)
up_reg_names <- subset(vol_df_3, logfc > 10 & pval > 6.21)
down_reg_lab <- (down_reg_names[1:10,])
up_reg_lab <- (up_reg_names[1:10,])
blues <- brewer.pal(9, name = "Blues")
labels <- rbind(up_reg_lab, down_reg_lab)



blues <- brewer.pal(9, name = "Blues")
labels[is.na(labels)] <- "NA"

df <- data.frame(
  x = c(0),
  y = c(33),
  text = c("TI-F vs TI-M"))

vol_plot <- ggplot(vol_df_2, aes(x = logfc,
                                 y = pval)) +
    geom_point(data = non_sig,
             shape = 21,
             size = 2, 
             fill = "grey", 
             colour = "grey")+
    geom_point(data = sig,
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
       
    geom_hline(yintercept = -log10(6.03034071895515E-07),
             linetype = "dashed", col = "red") + 
    geom_vline(xintercept = c(-10, 10),
               linetype = "dashed", col = "red") +
    labs(x = "Beta Difference (%)", y = "-log10(Pval)") +
    ylim(c(0,35))+
 scale_x_continuous(breaks=c(-50,-40,-30,-20,-10, 0, 10,20,30,40,50), limits = c(-50,50)) +
    
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
           size=10)

ggsave(
  "./methylome/figures/volcano-t1f-pre-t1-m-pre-dmps.pdf",
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


## load DMPs

DMPs <- data.frame(fread("./methylome/results/final-dasen-dmps-T2preFVsT2preM.csv"),
                    row.names=1)

## volcano pot p val TII-pre TII-post

cut_off <- -log10(3.17648749037999E-06)

logfc <- DMPs$logFC_beta * 100
pval <- DMPs$p.value
pval <- -log10(pval)
vol_df_2 <- cbind(logfc, pval)
vol_df_2 <- as.data.frame(vol_df_2)
vol_df_3 <- vol_df_2
vol_df_3$names <- DMPs$GeneName
vol_df_2 <- head(vol_df_2,10000)
non_sig <- subset(vol_df_2,  pval < 5.50)
sig <- subset(vol_df_2, logfc < 10 & pval > 5.50 | logfc > -10 & pval > 5.50)
down_reg <- subset(vol_df_2, logfc < 0 & pval > 5.50)
down_reg_names <- subset(vol_df_3, logfc < -10 & pval > 5.50)
up_reg <- subset(vol_df_2, logfc > 0 & pval > 5.50)
up_reg_names <- subset(vol_df_3, logfc > 10 & pval > 5.50)
down_reg_lab <- (down_reg_names[1:10,])
up_reg_lab <- (up_reg_names[1:10,])
blues <- brewer.pal(9, name = "Blues")
labels <- rbind(up_reg_lab, down_reg_lab)


blues <- brewer.pal(9, name = "Blues")
labels[is.na(labels)] <- "NA"

df <- data.frame(
  x = c(0),
  y = c(33),
  text = c("TII-F vs TII-M"))

vol_plot <- ggplot(vol_df_2, aes(x = logfc,
                                 y = pval)) +
    geom_point(data = non_sig,
             shape = 21,
             size = 2, 
             fill = "grey", 
             colour = "grey")+
    geom_point(data = sig,
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
       
    geom_hline(yintercept = 5,50,
             linetype = "dashed", col = "red") + 
    geom_vline(xintercept = c(-10, 10),
               linetype = "dashed", col = "red") +
    labs(x = "Beta Difference (%)", y = "-log10(Pval)") +
    ylim(c(0,35))+
    scale_x_continuous(breaks=c(-50,-40,-30,-20,-10, 0, 10,20,30,40,50), limits = c(-50,50)) +

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
              size=10)

ggsave(
  "./methylome/figures/volcano-t2f-pre-t2-m-pre-dmps.pdf",
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


## load DMPs

DMPs <- data.frame(fread("./methylome/results/final-dasen-dmps-WMpreFVsWMpreM.csv"),
                    row.names=1)

## volcano pot p val TII-pre TII-post

cut_off <- -log10(2.70274518163236E-06)

logfc <- DMPs$logFC_beta * 100
pval <- DMPs$p.value
pval <- -log10(pval)
vol_df_2 <- cbind(logfc, pval)
vol_df_2 <- as.data.frame(vol_df_2)
vol_df_3 <- vol_df_2
vol_df_3$names <- DMPs$GeneName
vol_df_2 <- head(vol_df_2,10000)
non_sig <- subset(vol_df_2,  pval < 5.56)
sig <- subset(vol_df_2, logfc < 10 & pval > 5.56 | logfc > -10 & pval > 5.56)
down_reg <- subset(vol_df_2, logfc < 0 & pval > 5.56)
down_reg_names <- subset(vol_df_3, logfc < -10 & pval > 5.56)
up_reg <- subset(vol_df_2, logfc > 0 & pval > 5.56)
up_reg_names <- subset(vol_df_3, logfc > 10 & pval > 5.56)
down_reg_lab <- (down_reg_names[1:10,])
up_reg_lab <- (up_reg_names[1:10,])
blues <- brewer.pal(9, name = "Blues")
labels <- rbind(up_reg_lab, down_reg_lab)



blues <- brewer.pal(9, name = "Blues")
labels[is.na(labels)] <- "NA"

df <- data.frame(
  x = c(0),
  y = c(33),
  text = c("WM-F vs WM-M"))

vol_plot <- ggplot(vol_df_2, aes(x = logfc,
                                 y = pval)) +
    geom_point(data = non_sig,
             shape = 21,
             size = 2, 
             fill = "grey", 
             colour = "grey")+
    geom_point(data = sig,
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
       
    geom_hline(yintercept = 5,56,
             linetype = "dashed", col = "red") + 
    geom_vline(xintercept = c(-10, 10),
               linetype = "dashed", col = "red") +
    labs(x = "Beta Difference (%)", y = "-log10(Pval)") +
    ylim(c(0,35))+
      scale_x_continuous(breaks=c(-50,-40,-30,-20,-10, 0, 10,20,30,40,50), limits = c(-50,50)) +
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
              size=10)

ggsave(
  "./methylome/figures/volcano-wmf-pre-wm-m-pre-dmps.pdf",
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



DMPs <- data.frame(fread("./methylome/results/final-dasen-dmps-Interaction1.csv"),
                    row.names=1)

## volcano pot p val TII-pre TII-post

cut_off <- -log10(2.4318948413519E-08)

logfc <- DMPs$logFC_beta * 100
pval <- DMPs$p.value
pval <- -log10(pval)
vol_df_2 <- cbind(logfc, pval)
vol_df_2 <- as.data.frame(vol_df_2)
vol_df_3 <- vol_df_2
vol_df_3$names <- DMPs$GeneName
vol_df_2 <- head(vol_df_2,10000)
non_sig <- subset(vol_df_2,  pval < 7.6)
sig <- subset(vol_df_2, logfc < 10 & pval > 7.6 | logfc > -10 & pval > 7.6)
down_reg <- subset(vol_df_2, logfc < 0 & pval > 7.6)
down_reg_names <- subset(vol_df_3, logfc < 0 & pval > 7.6)
up_reg <- subset(vol_df_2, logfc > 0 & pval > 7.6)
up_reg_names <- subset(vol_df_3, logfc > 0 & pval > 7.6)
down_reg_lab <- (down_reg_names[1:10,])
up_reg_lab <- (up_reg_names[1,])
blues <- brewer.pal(9, name = "Blues")
labels <- rbind(up_reg_lab, down_reg_lab)

blues <- brewer.pal(9, name = "Blues")
labels[is.na(labels)] <- "NA"

df <- data.frame(
  x = c(0),
  y = c(33),
  text = c("Interaction DMPs"))

vol_plot <- ggplot(vol_df_2, aes(x = logfc,
                                 y = pval)) +
    geom_point(data = non_sig,
             shape = 21,
             size = 2, 
             fill = "grey", 
             colour = "grey")+
    geom_point(data = sig,
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
       
    geom_hline(yintercept = 7.6,
             linetype = "dashed", col = "red") + 
    geom_vline(xintercept = c(-10, 10),
               linetype = "dashed", col = "red") +
    labs(x = "Beta Difference (%)", y = "-log10(Pval)") +
    ylim(c(0,35))+
  scale_x_continuous(breaks=c(-50,-40,-30,-20,-10, 0, 10,20,30,40,50), limits = c(-50,50)) +
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
              size=10)

ggsave(
  "./methylome/figures/volcano-f-m-ft-interactions-dmps.pdf",
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
