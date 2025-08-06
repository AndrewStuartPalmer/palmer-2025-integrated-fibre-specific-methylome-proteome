## Author: Andrew Palmer
## Title: Manhatten Plot of DMPs
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Manhatten Plot of DMPs
#####################################################################

#####################################################################
## Load Packages
#####################################################################


my_annotated_basic_subset <- my_annotated_basic[1:1000]
viridis_ramp <- colorRampPalette(viridis_colours[5:1])
transf.pval <- -log10(my_annotated_basic$ind.fdr)
 
points.col <- colByValue(transf.pval, colors=c("#BBBBBB00","#21908CFF"))

pp <- getDefaultPlotParams(plot.type=4)
pp$bottommargin <- 50
pp$rightmargin <- 0.02
png("./methylome/figures/manhatten-plot-cpgs.png", height = 400, width = 1800)
kp <- plotKaryotype(plot.type = 4, chromosomes = "autosomal", cex = 1.75,  srt=90, plot.params = pp)
kp <- kpPlotManhattan(kp, suggestive.col = "red", suggestiveline = (-log10(0.001)), data=my_annotated_basic, pval = my_annotated_basic$ind.fdr, points.col = points.col, genomewideline = (0), suggestive.lwd = 1)
kpAddLabels(kp, labels = "-log10 adj-p-val", srt=90, pos=3, side = "left", label.margin = 0.04, cex = 1.75, offset = -1)
kpAxis(kp, ymin=0, ymax=kp$latest.plot$computed.values$ymax, cex = 1.5, label.margin = 0.1)
kpAbline(kp, h = 0.091, col="orange", lwd=1.25)
dev.off()



## manhatten plot

viridis_ramp <- colorRampPalette(viridis_colours[5:1])
transf.pval <- -log10(my_annotated_basic_subset$ind.fdr)
 
points.col <- colByValue(transf.pval, colors=c("#BBBBBB00","#21908CFF"))

png("./methylome/figures/manhatten-plot-cpgs.png", height = 400, width = 1800)
kp <- plotKaryotype(plot.type = 4, chromosomes = "autosomal", cex = 1,  srt=90)

kp <- kpPlotManhattan(kp, suggestive.col = "red", suggestiveline = (-log10(0.001)), data=my_annotated_basic_subset, pval = my_annotated_basic_subset$ind.fdr, points.col = points.col, genomewideline = (0), suggestive.lwd = 1)
kpAddLabels(kp, labels = "-log10 adj-p-val", srt=90, pos=3, side = "left", label.margin = 0.04, cex.lab = 1.25)
kpAxis(kp, ymin=0, ymax=kp$latest.plot$computed.values$ymax, cex.axis = 1)
kpAbline(kp, h = 0.091, col="orange", lwd=1.25)
dev.off()
