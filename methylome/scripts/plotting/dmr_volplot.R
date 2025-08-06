## Author: Andrew Palmer
## Title: Volcano plots for DMRs
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Script to generate Volcano plots of DMRs
#####################################################################

#####################################################################
## Load Packages
#####################################################################

load_packages <- c("data.table", "RColorBrewer", "ggplot2", "pbapply")

lapply(load_packages, library, character.only = TRUE)

#####################################################################
## Load DMR ranges
#####################################################################

results.ranges <-  data.frame(fread("./methylome/results/dmrcate-DMRs-tI-pre-vs-tII-pre.csv"))
results.ranges <- results.ranges[results.ranges$no.cpgs >3,]

#####################################################################
## Prepare data
#####################################################################
logfc <- results.ranges$meandiff * 100
pval <- results.ranges$Fisher
pval <- abs(-log10(pval))
cpgs <- results.ranges$no.cpgs
vol_df_2 <- cbind(logfc, pval, cpgs)
vol_df_2 <- as.data.frame(vol_df_2)
vol_df_3 <- vol_df_2
vol_df_3$names <- results.ranges$symbol
down_reg <- subset(vol_df_2, logfc < 0 & pval > 0)
down_reg_names <- subset(vol_df_3, logfc < -10 & pval > 3.03)
up_reg <- subset(vol_df_2, logfc > 0 & pval > 0)
up_reg_names <- subset(vol_df_3, logfc > 10 & pval > 3.03)
down_reg_lab <- (down_reg_names[1:10,])
up_reg_lab <- (up_reg_names[1:10,])
blues <- brewer.pal(9, name = "Blues")
labs <- subset(vol_df_3, logfc > 20 & pval > 3.03 | logfc < -20 & pval >3.03)
labs$labels <- paste0(labs$names, "(", labs$cpgs, ")")
na_points <- subset(vol_df_3, is.na(names))

labs[1,5] <- "SYT8 TNNI2(20)"

cols <- c("hypo" = blues[8], "hyper" = blues[4])

#####################################################################
## Generate Plot
#####################################################################

df <- data.frame(
  x = c(0),
  y = c(280),
  text = c("TI vs TII"))

vol_plot <- ggplot(vol_df_3, aes(x = logfc,
                                 y = pval)) +
    geom_point(data = vol_df_3,
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
     geom_point(data = na_points,
             shape = 21,
             size = 2, 
             fill = "grey",
             colour = "grey")+
    
    geom_hline(yintercept = 3.03,
             linetype = "dashed", col = "red") + 
    geom_vline(xintercept = c(-10, 10),
               linetype = "dashed", col = "red") +
    
    
    ylim(0,320)+
    labs(x = "Mean Beta Difference (%)", y = "-log10(Pval)") +
    scale_x_continuous(breaks=c(-40,-30,-20,-10, 0, 10,20,30,40), limits = c(-40,40) ) +
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
    
    geom_label_repel(data = labs,
                     aes(label = labels),
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
## Save plot
#####################################################################

ggsave(
  "./methylome/figures/volcano-t1-pre-t2pre-DMR-robust-final.pdf",
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
