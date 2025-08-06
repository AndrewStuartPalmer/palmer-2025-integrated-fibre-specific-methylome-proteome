## Author: Andrew Palmer
## Title: DMR chromosome plot 
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## DMR chromosome plot using karyoploter
#####################################################################

#####################################################################
## Load Packages
#####################################################################

load_packages <- c("data.table", 
                   "RColorBrewer", 
                   "karyoploteR")

lapply(load_packages, library, character.only = TRUE)

#####################################################################
## Load Granges object
#####################################################################

results.ranges <-  data.frame(fread("./methylome/results/dmrcate-DMRs-tI-pre-vs-tII-pre.csv"))
results.ranges <- results.ranges[results.ranges$no.cpgs >3,]

Granges <- makeGRangesFromDataFrame(results.ranges)
Granges$genes <- results.ranges$symbol
Granges$mean_diff <- results.ranges$meandiff
Granges$fdr <- results.ranges$Fisher
Granges$nProbes <- results.ranges$no.cpgs

Granges <- Granges[Granges$fdr <0.001,]

Granges_hypo <- Granges[Granges$mean_diff < 0,]
Granges_hyper <- Granges[Granges$mean_diff > 0,]

#####################################################################
## Prepare plotting paramaters
#####################################################################

fc.ymax <- ceiling(max(abs(range(Granges$mean_diff))))
fc.ymin <- -fc.ymax
log_fdr <- abs(log(Granges$fdr))
cex.val <- sqrt(log_fdr)/8

blues <-  brewer.pal(8, name = "Blues")
col.over <- blues[8]
col.under <- blues[6]
sign.col <- rep(col.over, length(Granges))
sign.col[Granges$mean_diff <0] <- col.under
top.genes <- Granges[1:30,]
points.top <- 0.7

top.genes$genes[1] <- "SYT8, TNNI2"
top.genes$genes[2] <- "SLC16A3, CSNK1D"
top.genes$genes[7] <- "SPON2, CTBP1"
top.genes$genes[11] <-"MIRLET7, MIR4763"
top.genes$genes[16] <- "LAD1"
top.genes$genes[18] <- "RASA3, RASA3-IT1"
top.genes$genes[27] <- "PDLIM1"
top.genes$genes[29] <- "B3GNT7"
top.genes$genes[30] <- "ATP2A2"

####################################################################
## Generate and save plot
#####################################################################

pp <- getDefaultPlotParams(plot.type=4)
pp$bottommargin <- 40
pp$topmargin <- 50
pp$leftmargin <- 0.06
pp$rightmargin <- 0.01

pdf("./methylome/figures/chromosome-plot-final.pdf", height = 10, width = 24)
kp <- plotKaryotype(plot.type = 4, genome = "hg38", chromosomes="autosomal", cex = 1.75,  srt=90,  plot.params = pp)

kpPoints(kp, data=Granges, y=Granges$mean_diff, cex=cex.val, ymax=0.4, ymin=-0.4, col=sign.col, r1=points.top)

kpAxis(kp, ymax=0.4, ymin=-0.4, r1=points.top,  label.margin = 0.04, cex=1.75)

kpAddLabels(kp, labels = "Mean Methylation Difference", srt=90, pos=1, label.margin = 0.05, r1=points.top, cex=1.75)

gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2

kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$mean_diff, y1=0.4, ymax=0.4, ymin=-0.4, r1=points.top, line.color = "#777777")

kpPlotMarkers(kp, top.genes, labels = top.genes$genes, text.orientation = "vertical", r0=points.top, label.dist = 0.0000009, label.color="#444444", adjust.label.position=TRUE, max.iter = 5000, cex = 1)

dev.off()

