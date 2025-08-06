## Author: Andrew Palmer
## Title: miss methyl ORA
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

####################################################################
## Hypergeometric overrepresentation testing.
####################################################################

## this code uses modified functions taken from the miss methyl packages
## required
## list of gene sets for testing
## significant CpGs as probe IDs
## all cpgs measured as probe IDs

####################################################################
##  Load required libraries and miss methyl updated functions
####################################################################

load_packages <- c("data.table",
                   "org.Hs.eg.db",
                   "GO.db",
                   "DMRcate",
                   "ExperimentHub",
                   "clusterProfiler",
                   "GenomicRanges",
                   "ggplot2")

lapply(load_packages, library, character.only = TRUE)

source("./methylome/scripts/util-functions/miss-methyl-functions.R")
anno <- read.delim("./methylome/annotations/epicv2.hg38.manifest.gencode.v41.tsv")

####################################################################
## Get GO Terms
####################################################################

go <- getGO()
collection <- go$idList

## remove go terms with more than 500 and less than 5 terms
collection <- collection[sapply(collection, length) > 5]
collection <- collection[sapply(collection, length) < 500]

## remove any missing terms from collection
collection <- lapply(collection, function(x) x[!is.na(x)])

####################################################################
## prepare DMPs (sig.cpgs) and all probes
####################################################################


dmps <- data.frame(fread("./methylome/results/final-dasen-dmps-T1preVsT2pre.csv"), row.names = 1)

## get all.cpgs run with limma
all.cpg <- rownames(dmps)
all.cpg <- as.character(all.cpg)
all.cpg <- all.cpg[!is.na(all.cpg)]
all.cpg <- unique(all.cpg)

## get sig.cpgs 

anno_sub <- anno[match(all.cpg, anno$probeID),]

cpgs <- GenomicRanges::GRanges(seqnames = anno_sub$CpG_chrm, 
                                   ranges = IRanges::IRanges(start = anno_sub$CpG_beg,                                                              end = anno_sub$CpG_end),
                                   strand = anno_sub$probe_strand,
                               name = anno_sub$probeID)

regions <- read.csv("./methylome/results/dmrcate-DMRs-tI-pre-vs-tII-pre.csv")

regions <- makeGRangesFromDataFrame(regions)

overlaps <- GenomicRanges::findOverlaps(cpgs,regions)

sig.cpg <- cpgs$name[from(overlaps)]

out <- getMappedEntrezIDs(sig.cpg, all.cpg)

sorted.eg.sig <- out$sig.eg
eg.universe <- out$universe
freq_genes <- out$freq
test.de <- out$de
frac <- out$fract.counts
equiv <- out$equiv

collection <- lapply(collection, function(x) x[x %in% eg.universe])
## Remove collections with no genes left after universe filter
collection <- collection[sapply(collection, length) > 0]

####################################################################
## estimate PWF and run hypergeometric tests 
####################################################################

pwf <- estimatePWF(D=test.de,bias=as.vector(equiv))


results <- hypergeometricTest(pwf = pwf,
                                  collection = collection,
                                  sorted.eg.sig =  sorted.eg.sig,
                                  eg.universe = eg.universe,
                                  freq_genes = freq_genes,
                                  test.de = test.de,
                                  frac = frac,
                                  equiv = equiv)


terms <- AnnotationDbi::select(GO.db,
                               keys = rownames(results), 
                               columns=c("ONTOLOGY","TERM","DEFINITION"), 
                               keytype= "GOID")


results <- as.data.frame(results)
results$GeneRatio <- results$DE / length(sig.cpg)
results$BgRatio <- results$DE / length(all.cpg)

ora_results <- cbind(results, terms)

results_df_sig <- ora_results[ora_results$FDR < 0.05,]
results_df_sig_0.001 <- ora_results[ora_results$FDR < 0.001,]

results_df_sig <- results_df_sig[order(results_df_sig$FDR),]

results_df_sig_bp_0.001 <- results_df_sig_0.001[results_df_sig_0.001$ONTOLOGY == "BP",]
results_df_sig_cc_0.001 <- results_df_sig_0.001[results_df_sig_0.001$ONTOLOGY == "CC",]
results_df_sig_mf_0.001 <- results_df_sig_0.001[results_df_sig_0.001$ONTOLOGY == "MF",]

ora_results <- ora_results[order(ora_results$FDR), ]

write.csv(ora_results, "./methylome/results/missmethyl_go_dmrs.csv", row.names = TRUE) 

theme_clean <- function() {
  theme_minimal() +
      theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          strip.text = element_text(size = rel(1), hjust = 0),
          strip.background = element_rect(fill = "grey80", color = NA),
          axis.text=element_text(size=14),
          axis.title=element_text(size=14),
          axis.text.y = element_text(size =20),
          panel.spacing = unit(1, "lines"))
}

pdf("./methylome/figures/miss-methyl-go-t1vst2-robust-dmrs-final.pdf", height = 10, width = 12)
ggplot(data = results_df_sig_0.001, aes(x = GeneRatio , y = reorder(TERM, -FDR), color = FDR, size = DE)) +
    geom_point() +
    scale_size_continuous(range = c(4,10), name = "Gene Count")+
    theme_bw() +
    scale_color_viridis_b(name = "Significance \n P.adjust")+
        ylab("") +
    xlab("GeneRatio") +
    ggtitle("GO Pathways TI vs TII")+
    theme_clean() +
    facet_grid(ONTOLOGY ~ ., scales = "free", space = "free")

dev.off()

####################################################################
## Get Reactome Terms Terms
####################################################################
reactome_url <- "https://reactome.org/download/current/ReactomePathways.gmt.zip"
reactome_filename <- "ReactomePathways.gmt.zip"
reactome_path <- "./methylome/genesets/"

download.file(url = reactome_url, destfile = paste0(reactome_path,
                                                    reactome_filename,
                                                    sep = ""))


reactome_set <- read.gmt(unzip("./methylome/genesets/ReactomePathways.gmt.zip"))

ENTRZIID <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                       reactome_set$gene,
                                       keytype =  "SYMBOL",
                                       column = "ENTREZID")

reactome_set$gene <- ENTRZIID

reactome_set <- split(reactome_set$gene, reactome_set$term)

collection <- reactome_set

## remove go terms with more than 500 and less than 5 terms
collection <- collection[sapply(collection, length) > 5]
collection <- collection[sapply(collection, length) < 500]

collection <- lapply(collection, function(x) x[!is.na(x)])


collection <- lapply(collection, function(x) x[x %in% eg.universe])

## Remove collections with no genes left after universe filter
collection <- collection[sapply(collection, length) > 0]

results <- hypergeometricTest(pwf = pwf,
                                  collection = collection,
                                  sorted.eg.sig =  sorted.eg.sig,
                                  eg.universe = eg.universe,
                                  freq_genes = freq_genes,
                                  test.de = test.de,
                                  frac = frac,
                              equiv = equiv)

results <- as.data.frame(results)
results$GeneRatio <- results$DE / length(sig.cpg)
results$BgRatio <- results$DE / length(all.cpg)

reactome_results <- results

reactome_df_sig <- reactome_results[reactome_results$FDR < 0.05,]

reactome_df_sig <- reactome_df_sig[order(reactome_df_sig$FDR),]
reactome_df_sig$Description <- rownames(reactome_df_sig)

pdf("./methylome/figures/miss-methyl-reactome-t1vst2-dmr.pdf", height = 4, width = 12)
ggplot(data = reactome_df_sig, aes(x = GeneRatio , y = reorder(Description, -FDR), color = FDR, size = DE)) +
  geom_point() +
    scale_size_continuous(range = c(4,10), name = "Gene Count")+
    theme_bw() +
    scale_color_viridis_b(name = "Significance \n P.adjust")+    
    ylab("") +
    xlim(c(0,0.0035))+
    xlab("GeneRatio") +
    ggtitle("Reactome Pathways TI vs TII")+
    theme_clean()
  
dev.off()

write.csv(reactome_results, "./methylome/results/missmethyl_reactome_dmr.csv", row.names = TRUE) 
