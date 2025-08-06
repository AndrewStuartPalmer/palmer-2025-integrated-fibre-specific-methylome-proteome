## Author: Andrew Palmer
## Title: Miss Methyl Region Analysis on Sex data 
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Script for identifying GO and reactome terms
#####################################################################
load_packages <- c("data.table",
              "org.Hs.eg.db",
              "GO.db",
              "DMRcate",
              "ExperimentHub",
              "clusterProfiler")

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

dmps <- data.frame(fread("./methylome/results/final-dasen-dmps-T1preFVsT1preM.csv"), row.names = 1)

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

regions_1 <- read.csv("./methylome/results/final-dmrcate-dmrs-tI_FvM.csv")
regions_1 <- makeGRangesFromDataFrame(regions_1, na.rm = TRUE)

regions_2 <- read.csv("./methylome/results/final-dmrcate-dmrs-tII_FvM.csv")
regions_2 <- makeGRangesFromDataFrame(regions_2, na.rm = TRUE)

regions_3 <- read.csv("./methylome/results/final-dmrcate-dmrs-wm_FvM.csv")
regions_3 <- makeGRangesFromDataFrame(regions_3, na.rm = TRUE)

overlaps_1 <- GenomicRanges::findOverlaps(cpgs,regions_1)
overlaps_2 <-  GenomicRanges::findOverlaps(cpgs,regions_2)
overlaps_3 <-  GenomicRanges::findOverlaps(cpgs,regions_3)

sig.cpg1 <- cpgs$name[from(overlaps_1)]
sig.cpg2 <- cpgs$name[from(overlaps_2)]
sig.cpg3 <- cpgs$name[from(overlaps_3)]
sig.cpg <- c(sig.cpg1, sig.cpg2, sig.cpg3)
sig.cpg <- unique(sig.cpg)

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

write.csv(ora_results, "./methylome/results/sex-missmethyl_go_dmrs.csv", row.names = TRUE) 


####################################################################
## Get Reactome Terms Terms
####################################################################

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
