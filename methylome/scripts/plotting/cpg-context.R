## Author: Andrew Palmer
## Title: CpG context plot
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Script to generate CpG context plot
#####################################################################

#####################################################################
## Load Packages
#####################################################################

library(data.table)
library(viridis)

#####################################################################
## Load data
#####################################################################

dmps <- data.frame(fread("./methylome/results/final-dasen-dmps-T1preVsT2pre.csv"), row.names = 1)

#####################################################################
## Prepare data
#####################################################################

sig_dmps <- subset(dmps, adj.p.val < 0.001 & logFC_beta > 0.1 | logFC_beta < -0.1)

sig_dmps$CGI <- gsub('CGI;', '', sig_dmps$CGI)
CGI <- gsub('CGI;','', sig_dmps$CGI)
cgi_table <- table(sig_dmps$CGI)
cgi_prop_table <- prop.table(cgi_table)*100
Island <- 2589/length(sig_dmps$CGI)
OpenSea <- 48786/length(sig_dmps$CGI)
Shelf <- 5407/length(sig_dmps$CGI)
Shore <- 10395/length(sig_dmps$CGI)

cgi_precent <- as.matrix(cbind(Island,OpenSea,Shelf,Shore))
my_cols <- viridis(n=4, alpha=0.7)

#####################################################################
## Generate Plot
#####################################################################

pdf("./methylome/figures/cg-context-sig-dmps.pdf", width = 12, height = 8)
par(cex= 1.15, cex.main=1.75, cex.lab=1.5, cex.axis=1.5, family = "Helvetica")
barplot(cgi_prop_table,
        ylim = c(0,100),
        col = my_cols,
        axisnames = TRUE,
        names.arg = names(cgi_prop_table),
        yaxt = "n",
        frame.plot= FALSE)
title("")
axis(2, at = seq(0, 100, 20), col = NA, line = 0)
mtext("Percentage of DMPs (%)", line = 2.5, side = 2, cex = 1.75)
grid(lty =1, nx = NA, ny = NULL)
dev.off()
