## Author: Andrew Palmer
## Title: Protein Annotations
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Script to add additional annotations to the Protein Results ##
#####################################################################

load_packages <- c("org.Hs.eg.db")

lapply(load_packages, library, character.only = TRUE)

dep <- read.csv("./proteome/results/t1vst2-pre-proteome.csv",
                header = TRUE,
                sep = ",",
                row.names = 1)

annotation_protein <- dep[,c("symbol","uniprot")]


####################################################################
## Append additional Gene IDs ##
####################################################################

annotation_protein$ENTREZIID <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                       rownames(annotation_protein),
                                       keytype =  "SYMBOL",
                                       column = "ENTREZID")

annotation_protein$UNIPROT <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                     rownames(annotation_protein),
                                     keytype =  "SYMBOL",
                                     column = "UNIPROT")

annotation_protein$ALIAS <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                   rownames(annotation_protein),
                                   keytype =  "SYMBOL",
                                   column = "ALIAS")


annotation_protein$SYMBOL <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                    rownames(annotation_protein),
                                    keytype =  "SYMBOL",
                                    column = "SYMBOL")

annotation_protein$ENSEMBLEPROT <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                    rownames(annotation_protein),
                                    keytype =  "SYMBOL",
                                    column = "ENSEMBLPROT")

annotation_protein$ENSEMBLE <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                    rownames(annotation_protein),
                                    keytype =  "SYMBOL",
                                    column = "ENSEMBL")

####################################################################
## Save annotated Proteins ##
####################################################################

write.csv(annotation_protein, "./proteome/results/annotated-proteins.csv")
