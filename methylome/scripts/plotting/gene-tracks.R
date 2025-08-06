## Author: Andrew Palmer
## Title: Gene Tracks plots
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Script for generating Gene PLot Tracks using GViz
#####################################################################

#####################################################################
## Load Packages
#####################################################################
load_packages <- c("Gviz",
                   "data.table",
                   "viridis",
                   "rtracklayer",
                   "TxDb.Hsapiens.UCSC.hg38.knownGene",
                   "org.Hs.eg.db",
                   "Homo.sapiens",
                   "RColorBrewer")

lapply(load_packages, library, character.only = TRUE)

#####################################################################
## Laod TxDB database
#####################################################################

txDB <- TxDb.Hsapiens.UCSC.hg38.knownGene
reds <- brewer.pal(8, name = "Reds")

#####################################################################
## Laod Granges object
#####################################################################

granges <-  data.frame(data.table::fread("./methylome/results/final-dasen-annotated-g-ranges-t1-t2-pre-m-full.csv"))

granges <- granges[,-5]


granges <- bsseq::data.frame2GRanges(granges, keepColumns = TRUE)

beta <- data.frame(data.table::fread("./methylome/data/processed-data/full-dasen-normalised-filtered-beta.csv"), row.names =1)

dmrs <-  data.frame(fread("./methylome/results/dmrcate-DMRs-tI-pre-vs-tII-pre.csv"))
dmrs <- makeGRangesFromDataFrame(dmrs)

beta_rep <- beta[,c(17,18,19,90,91)]
beta$WM.f28.preWM.1 <- rowMeans(beta_rep[,4:5])
beta$WM.s201.preWM.1 <- rowMeans(beta_rep[,1:3])

beta_filt <- beta[,- c(1,2,17,18,30,57,58,90)]

pheno <- read.delim("./methylome/experiment-design/pheno-table.txt")

pheno$group <- paste(pheno$cell_type, pheno$timepoint, sep = "_")
pheno_filt <- pheno[-c(1,2,17,18,30,57,58,90),] 

pheno_filt <- pheno_filt[-c(57),] 
beta_filt <- beta_filt[,- c(57)]


#####################################################################
## Remove 4wWM samples
#####################################################################

subset_position <- endsWith(colnames(beta_filt), "preT1")
pheno_filt_subset_t1 <- subset(pheno_filt, subset = subset_position)
beta_filt_subset_t1 <- subset(beta_filt, select = subset_position)

subset_position <- endsWith(colnames(beta_filt), "preT2")
pheno_filt_subset_t2 <- subset(pheno_filt, subset = subset_position)
beta_filt_subset_t2 <- subset(beta_filt, select = subset_position)

subset_position <- grepl("preWM", colnames(beta_filt))
pheno_filt_subset_wm <- subset(pheno_filt, subset = subset_position)
beta_filt_subset_wm <- subset(beta_filt, select = subset_position)


#####################################################################
## Prepare beta values
#####################################################################

beta_val_ready <- cbind(beta_filt_subset_t1, beta_filt_subset_t2, beta_filt_subset_wm)

pheno_filt_ready <- rbind(pheno_filt_subset_t1,pheno_filt_subset_t2,pheno_filt_subset_wm)

tI_vals <- rowMeans(beta_filt_subset_t1)
tII_vals <- rowMeans(beta_filt_subset_t2)
wm_vals <- rowMeans(beta_filt_subset_wm)

beta_vals <- as.data.frame(cbind(tI_vals, tII_vals, wm_vals))

#####################################################################
## Prepare Gene regions of interest
#####################################################################

myh_locis_region <- toGRanges("chr14:23,385,520-23,444,000")
myh7_locus <- subsetByOverlaps (granges, myh_locis_region)


myh7_keep <- rownames(beta_vals) %in% myh7_locus$probeID

myh7_keep <- rownames(beta_vals) %in% granges$probeID

beta_val_myh7 <- beta_vals[myh7_keep,]

beta_val_myh7 <- beta_val_myh7[match(myh7_locus$probeID, rownames(beta_val_myh7)),]

myh7_locus_df <- data.frame(myh7_locus[,])

granges_myh7_full <- cbind(myh7_locus_df[,1:3], beta_val_myh7)
granges_myh7_full$seqnames <- as.factor(granges_myh7_full$seqnames)
granges_myh7_full$genome <- as.factor("hg38")
myh7_dmrs <- subsetByOverlaps(dmrs, myh_locis_region)


##chr14:23,385,520-23,444,000
cell_types <- c(pheno_filt_subset_t1$cell_type, pheno_filt_subset_t2$cell_type, pheno_filt_subset_wm$cell_type)

gen <- "hg38"
chr <- "chr14"

myh7_dmrs@elementMetadata@listData$dmr <- "DMR"
groups_dmr <- myh7_dmrs@elementMetadata@listData$dmr

myh7_start <-  23370000
myh7_end <-  23454000

dmr_track <- AnnotationTrack(myh7_dmrs, name = "", chromosome = "chr14", background.title = "white", cex.group = 1.5, feature = groups_dmr, DMR = "darkblue", groupAnnotation = "feature")

## strack <- SequenceTrack(Hsapiens, chromosome = chr)

itrack <- IdeogramTrack(genome = "hg38", chromosome = "chr14")

atrack <- AnnotationTrack(myh7_locus, name = "CpG", chromosome = "chr14", isPaired=FALSE, background.title = "grey", lex = 1)

gtrack <- GenomeAxisTrack()

rtrack <- GeneRegionTrack(txDB, genome = gen,
                          start = myh7_start,
                          end =  myh7_end,
                          chromosome = "chr14",
                          name = "",
                          showId = TRUE,
                          geneSymbol = TRUE,
                          symbol = symbol,
                          background.title = "white",
                          col = "white", fill = "black",
                          collapseTranscripts = "meta",
                          cex = 1.5,
                          cex.group = 1.5,
                          cex.names = 1.5,
                          cex.lab = 1.5)

SYMBOL <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                       gene(rtrack),
                                       keytype =  "ENTREZID",
                                column = "SYMBOL")

symbol(rtrack) <- SYMBOL

beta_data <- DataTrack(data = beta_val_myh7,
                       start = granges_myh7_full$start,
                       end = granges_myh7_full$end,
                       chromosome = "chr14",
                       genome = gen,
                       name = "DNAm Beta Vals",
                       ylim = c(0,1),
                       type = c("g", "p", "a"),
                       groups = c("TI", "TII", "WM"),
                       lty = 2,
                       lwd = 3,
                       cex.axis = 1.75,
                       cex.legend = 1.75,
                       col.line = c(reds[8], reds[2], reds[4]),
                       col.symbol = c(reds[8], reds[2], reds[4]),
                       background.title = "grey",
                       cex = 2)


#####################################################################
## Generate and Save plot
#####################################################################

pdf("./methylome/figures/myh7-myh6-tracks.pdf", height = 12, width = 24)
plotTracks(list(beta_data, dmr_track, rtrack, gtrack, itrack),
           col = NULL,
           from = myh7_start, to = myh7_end, cex.title = 1.5,
            sizes = c(10,1,1,2,2))
dev.off()

#####################################################################
## Complete the same for MFAP4
#####################################################################

reds <- brewer.pal(8, name = "Reds")

mfap4_start <- 19382700
mfap4_end <- 19388127
mfap4_locis_region <- toGRanges("chr17:19,382,509-19,388,127")
mfap4_locus <- subsetByOverlaps (granges, mfap4_locis_region)

mfap4_keep <- rownames(beta_vals) %in% mfap4_locus$probeID

myh7_keep <- rownames(beta_vals) %in% granges$probeID

beta_val_mfap4 <- beta_vals[mfap4_keep,]

beta_val_mfap4 <- beta_val_mfap4[match(mfap4_locus$probeID, rownames(beta_val_mfap4)),]

mfap4_locus_df <- data.frame(mfap4_locus[,])

granges_mfap4_full <- cbind(mfap4_locus_df[,1:3], beta_val_mfap4)
granges_mfap4_full$seqnames <- as.factor(granges_mfap4_full$seqnames)
granges_mfap4_full$genome <- as.factor("hg38")


mfap4_dmrs <- subsetByOverlaps(dmrs, mfap4_locis_region)
mfap4_dmrs@elementMetadata@listData$dmr <- "DMR"
groups_dmr <- mfap4_dmrs@elementMetadata@listData$dmr


itrack <- IdeogramTrack(genome = "hg38", chromosome = "chr17")

atrack <- AnnotationTrack(mfap4_locus, name = "CpG", chromosome = "chr17", isPaired=FALSE, background.title = "grey", lex = 1)

gtrack <- GenomeAxisTrack()


dmr_track <- AnnotationTrack(mfap4_dmrs, name = "", chromosome = "chr17", background.title = "white", cex.group = 1.5, feature = groups_dmr, DMR = "darkblue", groupAnnotation = "feature")

rtrack <- GeneRegionTrack(txDB, genome = gen,
                          start = mfap4_start,
                          end =  mfap4_end,
                          chromosome = "chr17",
                          name = "",
                          showId = TRUE,
                          geneSymbol = TRUE,
                          symbol = symbol,
                          background.title = "white",
                          col = "white", fill = "black",
                          collapseTranscripts = "meta",
                          cex = 1.5,
                          cex.group = 1.5,
                          cex.names = 1.5,
                          cex.lab = 1.5)

SYMBOL <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                       gene(rtrack),
                                       keytype =  "ENTREZID",
                                column = "SYMBOL")

symbol(rtrack) <- SYMBOL

beta_data <- DataTrack(data = beta_val_mfap4,
                       start = granges_mfap4_full$start,
                       end = granges_mfap4_full$end,
                       chromosome = "chr17",
                       genome = gen,
                       name = "DNAm Beta Vals",
                       ylim = c(0,1),
                       type = c("g", "p", "a"),
                       groups = c("TI", "TII", "WM"),
                       lty = 2,
                       lwd = 3,
                       cex.axis = 1.75,
                       cex.legend = 1.75,
                       col.line = c(reds[8], reds[2], reds[4]),
                       col.symbol = c(reds[8], reds[2], reds[4]),
                       background.title = "grey",
                       cex = 2)

pdf("./methylome/figures/mfap4-gene-track.pdf", height = 12, width = 24)
plotTracks(list(beta_data, dmr_track, rtrack, gtrack, itrack),
           col = NULL,
           from = mfap4_start, to = mfap4_end, cex.title = 1.5,
            sizes = c(10,1,1,2,2))
dev.off()


