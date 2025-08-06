## Author: Andrew Palmer
## Title: Over Representation Analysis: Reactome
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Script for over-representation analysis of proteins ##
#####################################################################

## This script is the sixth step in the proteomic analysis workflow.
## This code uses ClusterProfiler to assess over represented Reactome
## terms amoung the differentially regulated proteins.
## Producing a reactome plot.

####################################################################
## Load required packages ##
####################################################################

load_packages <- c("clusterProfiler", 
                   "limma", 
                   "viridis", 
                   "org.Hs.eg.db", 
                   "ggplot2", 
                   "ReactomePA")

lapply(load_packages, library, character.only = TRUE)

####################################################################
## Load Data ##
####################################################################

dep <- read.csv("./proteome/results/t1vst2-pre-proteome.csv",
                header = TRUE,
                sep = ",",
                row.names = 1)

experiment_design <- read.delim("./proteome/experiment-design/lfq_analyst_experimental_design.txt")


####################################################################
## Download Required Pathways ##
####################################################################

dir.create("./proteome/data/genesets", recursive = TRUE)

## Reactome

reactome_url <- "https://reactome.org/download/current/ReactomePathways.gmt.zip"
reactome_filename <- "ReactomePathways.gmt.zip"
reactome_path <- "./proteome/data/genesets/"

download.file(url = reactome_url, destfile = paste0(reactome_path,
                                                    reactome_filename,
                                                    sep = ""))

reactome_set <- read.gmt(unzip("./proteome/data/genesets/ReactomePathways.gmt.zip"))
reactome_set_fgsea <- gmtPathways(unzip("./proteome/data/genesets/ReactomePathways.gmt.zip"))

####################################################################
## Append additional Gene IDs ##
####################################################################

dep$ENTREZIID <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                       rownames(dep),
                                       keytype =  "SYMBOL",
                                       column = "ENTREZID")

dep$UNIPROT <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                     rownames(dep),
                                     keytype =  "SYMBOL",
                                     column = "UNIPROT")

dep$ALIAS <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                   rownames(dep),
                                   keytype =  "SYMBOL",
                                   column = "ALIAS")


dep$SYMBOL <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                    rownames(dep),
                                    keytype =  "SYMBOL",
                                    column = "SYMBOL")

####################################################################
## set up ggplot2 custom theme ##
####################################################################

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
####################################################################
## Prepare Genes and Genesets for ORA ##
####################################################################

organism <- "org.Hs.eg.db"

defup <- subset(dep, adj.P.Val <= 0.001 & logFC > 0)

defdn <- subset(dep, adj.P.Val <= 0.001 & logFC < 0)

bg <- rownames(dep)
bg_uniprot <- dep$UNIPROT

message("number of genes in each group")
lapply(list("background" = bg,
            "up-regulated" = defup[, 2],
            "down-regulated" = defdn[, 2]),
       length)

bgdf <- data.frame("background", bg)
colnames(bgdf) <- c("term", "gene")
reactome_set <- rbind(reactome_set, bgdf)

####################################################################
## REACTOME Overrepresentation Analysis ##
####################################################################

ora_up_reactome <- enricher(gene = rownames(defup),
                            universe = bg,
                            maxGSSize = 50000,
                            TERM2GENE = reactome_set,
                            pAdjustMethod = "fdr",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2)

ora_up_reactome_005 <- enricher(gene = rownames(defup),
                            universe = bg,
                            maxGSSize = 50000,
                            TERM2GENE = reactome_set,
                            pAdjustMethod = "fdr",
                            pvalueCutoff = 0.005,
                            qvalueCutoff = 0.2)

ora_up_reactome_001 <- enricher(gene = rownames(defup),
                            universe = bg,
                            maxGSSize = 50000,
                            TERM2GENE = reactome_set,
                            pAdjustMethod = "fdr",
                            pvalueCutoff = 0.001,
                            qvalueCutoff = 0.2)

ora_up_reactome_df <- as.data.frame(ora_up_reactome)

ora_dn_reactome <- enricher(gene = rownames(defdn),
                            universe = bg,
                            maxGSSize = 50000,
                            TERM2GENE = reactome_set,
                            pAdjustMethod = "fdr",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2)

ora_dn_reactome_005 <- enricher(gene = rownames(defdn),
                            universe = bg,
                            maxGSSize = 50000,
                            TERM2GENE = reactome_set,
                            pAdjustMethod = "fdr",
                            pvalueCutoff = 0.005,
                            qvalueCutoff = 0.2)

ora_dn_reactome_001 <- enricher(gene = rownames(defdn),
                            universe = bg,
                            maxGSSize = 50000,
                            TERM2GENE = reactome_set,
                            pAdjustMethod = "fdr",
                            pvalueCutoff = 0.001,
                            qvalueCutoff = 0.2)

 ora_dn_reactome_df <- as.data.frame(ora_dn_reactome)   

## up-reg-reactome
n_pw=10

ora_up_reactome_df$geneID <- NULL
ora_up_reactome_df <- subset(ora_up_reactome_df, p.adjust <0.05)
ora_up_reactome_dfs <- rownames(ora_up_reactome_df)

gr <- as.numeric(sapply(strsplit(ora_up_reactome_df$GeneRatio,"/"),"[[",1)) /
  as.numeric(sapply(strsplit(ora_up_reactome_df$GeneRatio,"/"),"[[",2))

br <- as.numeric(sapply(strsplit(ora_up_reactome_df$BgRatio,"/"),"[[",1)) /
  as.numeric(sapply(strsplit(ora_up_reactome_df$BgRatio,"/"),"[[",2))

ora_up_reactome_df$ES <- gr/br
ora_up_reactome_df <- ora_up_reactome_df[order(ora_up_reactome_df$p.adjust),]
ora_up_reactome_df$Description = NULL
ora_up_reactome_10 <- head(ora_up_reactome_df,10)

topup2 <- rev(head(ora_up_reactome_df$ES, 10))
ora_up_reactome_df[3, 1] <- "Repiratory electron transport..."
names(topup2) <- rev(head(ora_up_reactome_df$ID, 10))

## dn-reg-reactome

ora_dn_reactome_df$geneID <- NULL
ora_dn_reactome_df <- subset(ora_dn_reactome_df, p.adjust <0.05)
ora_dn_reactome_dfs <- rownames(ora_dn_reactome_df)

gr <- as.numeric(sapply(strsplit(ora_dn_reactome_df$GeneRatio,"/"),"[[",1)) /
  as.numeric(sapply(strsplit(ora_dn_reactome_df$GeneRatio,"/"),"[[",2))

br <- as.numeric(sapply(strsplit(ora_dn_reactome_df$BgRatio,"/"),"[[",1)) /
  as.numeric(sapply(strsplit(ora_dn_reactome_df$BgRatio,"/"),"[[",2))

ora_dn_reactome_df$ES <- gr / br
ora_dn_reactome_df <- ora_dn_reactome_df[order(ora_dn_reactome_df$p.adjust),]
ora_dn_reactome_df$Description = NULL
ora_dn_reactome_10 <- head(ora_dn_reactome_df,10)

ora_dn_reactome_10$regulation <- "DownRegulated"
ora_up_reactome_10$regulation <- "UpRegulated"
ora_reactome <- rbind(ora_dn_reactome_10, ora_up_reactome_10)

####################################################################
## Produce Reactome ORA plot ##
####################################################################

pdf("./proteome/figures/ora-reactome-t1prevst2pre-prot.pdf", height = 8, width = 14)
ggplot(data = ora_reactome, aes(x = ES ,y = reorder(ID, -p.adjust), color = p.adjust, size = Count)) +
    geom_point() +
    scale_size_continuous(range = c(4,8), name = "Gene Count")+
    theme_bw() +
    scale_color_viridis(name = "Significance \n P.adjust")+
    ylab("") + 
    xlab("ES") +
    xlim(0,10) +
    ggtitle("Reactome Pathways TI vs TII")+
    theme_clean()+
    facet_grid(. ~ regulation)

dev.off()

####################################################################
## Save Results ##
####################################################################

write.table(ora_reactome,
            file = "./proteome/results/t1vst2-pre-proteome-reactome.csv",
            sep = ",",
            col.names = NA)

