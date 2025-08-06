## Author: Andrew Palmer
## Title: Fibre-type WM overlap plots
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Script to generate fibre-type WM overlap plots
#####################################################################

#####################################################################
## Load Packages
#####################################################################

load_packages <- c("data.table",
                   "ggplot2",
                   "RColorBrewer",
                   "ggplotify",
                   "ggrepel",
                   "ggvenn")

lapply(load_packages, library, character.only = TRUE)

#####################################################################
## Load data
#####################################################################

tI_tII_dmps <- data.frame(fread("./methylome/results/final-dasen-dmps-T1preVsT2pre.csv"), row.names = 1)

tI_wm_dmps <- data.frame(fread("./methylome/results/final-dasen-dmps-T1preVsWMpre.csv"), row.names = 1)

tII_wm_dmps <- data.frame(fread("./methylome/results/final-dasen-dmps-T2preVsWMpre.csv"), row.names = 1)

sig_DMPs_hyper <- tI_tII_dmps[tI_tII_dmps$logFC_beta > 0.1 & tI_tII_dmps$adj.p.val <0.001,]
sig_DMPs_hypo <- tI_tII_dmps[tI_tII_dmps$logFC_beta < -0.1 & tI_tII_dmps$adj.p.val <0.001,]
sig_T1_WM_hypo <- tI_wm_dmps[tI_wm_dmps$logFC_beta < -0.1 & tI_wm_dmps$adj.p.val <0.001,]
sig_T1_WM_hyper <-  tI_wm_dmps[tI_wm_dmps$logFC_beta > 0.1 & tI_wm_dmps$adj.p.val <0.001,]
sig_T2_WM_hypo <-  tII_wm_dmps[tII_wm_dmps$logFC_beta < -0.1 & tII_wm_dmps$adj.p.val <0.001,]
sig_T2_WM_hyper <- tII_wm_dmps[tII_wm_dmps$logFC_beta > 0.1 & tII_wm_dmps$adj.p.val <0.001,]

#####################################################################
## Generate Overlap Lists
#####################################################################

v1 <- list("Hyper-T1vT2" = rownames(sig_DMPs_hyper), "Hyper-T1vWM" = rownames(sig_T1_WM_hyper), "Hypo-T1vWM" = rownames(sig_T1_WM_hypo),  "Hypo-T1vT2" = rownames(sig_DMPs_hypo))

v2 <- list("Hyper-T1vT2" = rownames(sig_DMPs_hyper), "Hyper-T2vWM" = rownames(sig_T2_WM_hyper), "Hypo-T2vWM" = rownames(sig_T2_WM_hypo) , "Hypo-T1vT2" = rownames(sig_DMPs_hypo))

brewer.pal(4, name = "Blues")

#####################################################################
## Generate venn plot
#####################################################################

venn <- ggvenn(
  v1, 
  fill_color = brewer.pal(4, name = "Blues"),
  stroke_size = 1.5, set_name_size = 7, text_size = 7
)

ggsave(
    "./methylome/figures/TIvTIIvWMvTI-final.pdf",
    plot = venn,
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


venn_2 <- ggvenn(
  v2, 
   fill_color = brewer.pal(4, name = "Blues"),
  stroke_size = 1.5, set_name_size = 7, text_size = 7
  )

ggsave(
    "./methylome/figures/TIvTIIvWMvTII-final.pdf",
    plot = venn_2,
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

#####################################################################
## Generate Volcano plots for FT WM comparisons
#####################################################################

## volcano plots

TI_p <- -log10(1.53917521841974E-05)
TII_p <- -log10(5.80063430296825E-06)


logfc <- tI_wm_dmps$logFC_beta * 100
pval <- tI_wm_dmps$p.value
pval <- -log10(pval)
vol_df_2 <- cbind(logfc, pval)
vol_df_2 <- as.data.frame(vol_df_2)
vol_df_3 <- vol_df_2
vol_df_3$names <-  tI_wm_dmps$GeneName
vol_df_2 <- head(vol_df_2,20000)
non_sig <- subset(vol_df_2,  pval < 4.8)
sig <- subset(vol_df_2, logfc < 10 & pval > 4.8 | logfc > -10 & pval > 4.8)
down_reg <- subset(vol_df_2, logfc < 0 & pval > 4.8)
down_reg_names <- subset(vol_df_3, logfc < -10 & pval > 4.8)
up_reg <- subset(vol_df_2, logfc > 0 & pval > 4.8)
up_reg_names <- subset(vol_df_3, logfc > 10 & pval > 4.8)
down_reg_lab <- (down_reg_names[1:6,])
up_reg_lab <- (up_reg_names[1:6,])
blues <- brewer.pal(9, name = "Blues")
labels <- rbind(up_reg_lab, down_reg_lab)

labels[10,3] <- "NA"
blues <- brewer.pal(9, name = "Blues")

df <- data.frame(
  x = c(0),
  y = c(32),
  text = c("TI vs WM"))

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
    
    geom_hline(yintercept =4.8 ,
             linetype = "dashed", col = "red") + 
    geom_vline(xintercept = c(-10, 10),
               linetype = "dashed", col = "red") +
    labs(x = "Beta Difference (%)", y = "-log10(Pval)") +
    scale_x_continuous(breaks=c(-30,-20,-10, 0, 10,20,30), limits = c(-30,30)) +


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
  "./methylome/figures/volcano-t1-pre-wm-pre-dmps.pdf",
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


logfc <- tII_wm_dmps$logFC_beta * 100
pval <- tII_wm_dmps$p.value
pval <- -log10(pval)
vol_df_2 <- cbind(logfc, pval)
vol_df_2 <- as.data.frame(vol_df_2)
vol_df_3 <- vol_df_2
vol_df_3$names <- tII_wm_dmps$GeneName
vol_df_2 <- head(vol_df_2,20000)
non_sig <- subset(vol_df_2,  pval < 5.23)
sig <- subset(vol_df_2, logfc < 10 & pval > 5.23 | logfc > -10 & pval > 5.23)
down_reg <- subset(vol_df_2, logfc < 0 & pval > 5.23)
down_reg_names <- subset(vol_df_3, logfc < -10 & pval > 5.23)
up_reg <- subset(vol_df_2, logfc > 0 & pval > 5.23)
up_reg_names <- subset(vol_df_3, logfc > 10 & pval > 5.23)
down_reg_lab <- (down_reg_names[1:6,])
up_reg_lab <- (up_reg_names[1:6,])
blues <- brewer.pal(9, name = "Blues")
labels <- rbind(up_reg_lab, down_reg_lab)

labels[6,3] <- "NA"
blues <- brewer.pal(9, name = "Blues")

df <- data.frame(
  x = c(0),
  y = c(32),
  text = c("TII vs WM"))

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
    
    geom_hline(yintercept = 5.23,
             linetype = "dashed", col = "red") + 
    geom_vline(xintercept = c(-10, 10),
               linetype = "dashed", col = "red") +
    labs(x = "Beta Difference (%)", y = "-log10(Pval)") +
    scale_x_continuous(breaks=c(-30,-20,-10, 0, 10,20,30), limits = c(-30,30))+

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
  "./methylome/figures/volcano-t2-pre-wm-pre-dmps.pdf",
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
