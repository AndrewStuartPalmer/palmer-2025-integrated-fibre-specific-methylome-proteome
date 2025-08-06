## Author: Andrew Palmer
## Title: Sex Karyoploter Plot
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Generation of Sex DMRs distribution across the chromosomes
#####################################################################

####################################################################
## Load packages ##
####################################################################

library(karyoploteR)
library(GenomicRanges)
library(RColorBrewer)
library(data.table)


####################################################################
## Prpare plot area ##
####################################################################

pp <- getDefaultPlotParams(plot.type=4)
pp$bottommargin <- 40
pp$topmargin <- 50
pp$leftmargin <- 0.06
pp$rightmargin <- 0.01


####################################################################
## prepare TI ranges ##
####################################################################

results_ranges_tI <-  data.frame(fread("./methylome/results/final-dmrcate-dmrs-tI_FvM.csv"))
results_ranges_tI <- results_ranges_tI[results_ranges_tI$no.cpgs > 3 & results_ranges_tI$Fisher < 0.001,]

Granges1 <- makeGRangesFromDataFrame(results_ranges_tI)
Granges1$genes <- results_ranges_tI$symbol
Granges1$mean_diff <- results_ranges_tI$meandiff 
Granges1$fdr <- results_ranges_tI$Fisher
Granges1$nProbes <- results_ranges_tI$no.cpgs

log_fdr1 <- abs(log(Granges1$fdr))
cex.val1 <- sqrt(log_fdr1)/3
top.genes1 <- Granges1[c(17,31,34),]


####################################################################
## prepare TII ranges ##
####################################################################

results_ranges_tII <-  data.frame(fread("./methylome/results/final-dmrcate-dmrs-tII_FvM.csv"))
results_ranges_tII <- results_ranges_tII[results_ranges_tII$no.cpgs > 3 & results_ranges_tII$Fisher < 0.001,]

Granges2 <- makeGRangesFromDataFrame(results_ranges_tII)
Granges2$genes <- results_ranges_tII$symbol
Granges2$mean_diff <- results_ranges_tII$meandiff 
Granges2$fdr <- results_ranges_tII$Fisher
Granges2$nProbes <- results_ranges_tII$no.cpgs

log_fdr2 <- abs(log(Granges2$fdr))
cex.val2 <- sqrt(log_fdr2)/3
top.genes2 <- Granges2[c(99,125,106,83,54),]


####################################################################
## prepare WM ranges  ##
####################################################################
results_ranges_wm <-  data.frame(fread("./methylome/results/final-dmrcate-dmrs-wm_FvM.csv"))
results_ranges_wm <- results_ranges_wm[results_ranges_wm$no.cpgs > 3 & results_ranges_wm$Fisher < 0.001,]

Granges3 <- makeGRangesFromDataFrame(results_ranges_wm)
Granges3$genes <- results_ranges_wm$symbol
Granges3$mean_diff <- results_ranges_wm$meandiff 
Granges3$fdr <- results_ranges_wm$Fisher
Granges3$nProbes <- results_ranges_wm$no.cpgs


log_fdr3 <- abs(log(Granges3$fdr))
cex.val3 <- sqrt(log_fdr3)/3



####################################################################
## Generate heatmap of DMRs TI ##
####################################################################


reds <-  brewer.pal(8, name = "Reds")
col.1 <- reds[8]
col.2 <- reds[2]
col.3 <- reds[4]

points.top <- 0.7


pdf("./methylome/figures/chromosome-sex-plot.pdf", height = 12, width =18)
kp <- plotKaryotype(plot.type = 4, genome = "hg38", chromosomes="autosomal", cex = 1.75,  srt=90,  plot.params = pp)

kpAxis(kp, ymax=0.3, ymin=-0.3, r1=points.top,  label.margin = 0.04, cex=2)
kpAddLabels(kp, labels = "Mean Methylation Difference", srt=90, pos=1, label.margin = 0.05, r1=points.top, cex=2)

kpPoints(kp, data=Granges1, y=Granges1$mean_diff, cex=cex.val1, ymax=0.3, ymin=-0.3, col=col.1, r1=points.top)

gene.mean1 <- start(top.genes1) + (end(top.genes1) - start(top.genes1))/2

kpSegments(kp, chr=as.character(seqnames(top.genes1)), x0=gene.mean1, x1=gene.mean1, y0=top.genes1$mean_diff, y1=0.3, ymax=0.3, ymin=-0.3, r1=points.top)

kpPlotMarkers(kp, top.genes1, labels = top.genes1$genes, text.orientation = "vertical", r0=0.5, label.dist = 0.05, label.color="#444444", adjust.label.position=TRUE, max.iter = 5000, cex = 1)

kpPoints(kp, data=Granges2, y=Granges2$mean_diff, cex=cex.val2, ymax=0.3, ymin=-0.3, col=col.2, r1=points.top)

gene.mean2 <- start(top.genes2) + (end(top.genes2) - start(top.genes2))/2

kpSegments(kp, chr=as.character(seqnames(top.genes2)), x0=gene.mean2, x1=gene.mean2, y0=top.genes2$mean_diff, y1=0.3, ymax=0.3, ymin=-0.3, r1=points.top)

kpPlotMarkers(kp, top.genes2, labels = top.genes2$genes, text.orientation = "vertical", r0=0.5, label.dist = 0.05, label.color="#444444", adjust.label.position=TRUE, max.iter = 5000, cex = 1)

kpPoints(kp, data=Granges3, y=Granges3$mean_diff, cex=cex.val3, ymax=0.3, ymin=-0.3, col=col.3, r1=points.top)

dev.off()


