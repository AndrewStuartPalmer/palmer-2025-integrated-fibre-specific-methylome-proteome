## Author: Andrew Palmer
## Title: CpG Plots
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Script to generate CpG plots
#####################################################################

#####################################################################
## Load Packages
#####################################################################

load_packages <- c("ggplot2", "RColorBrewer", "data.table",
                   "gghalves","ggbeeswarm", "ggpubr")

lapply(load_packages, library, character.only = TRUE)



#####################################################################
## Laod Data
#####################################################################
DMPs <- data.frame(fread("./methylome/results/final-dasen-dmps-T1preVsT2pre.csv"),
                    row.names=1)

betas <- data.frame(fread("./methylome/data/processed-data/full-dasen-normalised-filtered-beta.csv"), row.names=1)

beta_rep <- betas[,c(17,18,19,90,91)]
betas$WM.f28.preWM.1 <- rowMeans(beta_rep[,4:5])
betas$WM.s201.preWM.1 <- rowMeans(beta_rep[,1:3])

beta_filt <- betas[,- c(1,2,17,18,30,57,58,90)]

pheno <- read.delim("./methylome/experiment-design/pheno-table.txt")
pheno_filt <- pheno[-c(1,2,17,18,30,57,58,90),] 


subset_position <- !endsWith(colnames(beta_filt), "4wWM")

pheno_filt_subset <- subset(pheno_filt, subset = subset_position)
beta_filt_subset <- subset(beta_filt, select = subset_position)

## optional take out post samples

subset_position <- grepl("pre", colnames(beta_filt))

pheno_filt_subset <- subset(pheno_filt, subset = subset_position)
beta_filt_subset <- subset(beta_filt, select = subset_position)


#####################################################################
## Generate data for MYH 1,2,4,7
#####################################################################

## create a loop to define CpGs for plotting and a function to create CpG plots

## loop to prepare a table of CpGs to plot

gene_list <- c("MYH7", "MYH2", "MYH4", "MYH1")
datalist <- list()


for (i in gene_list) {
    # ... make some data
    subset_cpg <-  grep(paste("\\b",i,"\\b", sep = ""), DMPs$GeneName, ignore.case = FALSE, value = FALSE)
    DMPs_subset <- DMPs[subset_cpg,]
    names_gene <- gene_list
    top_cpg <- DMPs_subset[1,]
    cpg <- rownames(top_cpg)
    
    datalist[[cpg]] <- cpg # add it to your list
}

cpgs_to_plot <- as.data.frame(t(Reduce(rbind, datalist)))
names(cpgs_to_plot) <- gene_list

data_subset <- subset(beta_filt_subset, rownames(beta_filt_subset) %in% cpgs_to_plot)
data_subset <- data_subset[match(cpgs_to_plot, rownames(data_subset)),]
rownames(data_subset) <- paste0(rownames(data_subset),sep = "_", colnames(cpgs_to_plot))

t_data_subset <- as.data.frame(t(data_subset))
t_data_subset$cell_type <- pheno_filt_subset$cell_type
t_data_subset$cell_type <- as.factor(t_data_subset$cell_type)


#####################################################################
## Prepare plots
#####################################################################

cpg.plot <- function(r.i){

    reds <- brewer.pal(8, name = "Reds")

theme_clean <- function() {
  theme_minimal() +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_line(size = 1),
          plot.title = element_text(size = 20, hjust = 0.5, face = "bold", color = "black"),
          strip.text = element_text(size = rel(1), hjust = 0),
          strip.background = element_rect(fill = "grey80", color = NA),
          axis.text=element_text(size=20, color = "black"),
          axis.title=element_text(size=20, color = "black"),
          legend.position="bottom",
          legend.key.size = unit(1.5, 'cm'), #change legend key size
          legend.key.height = unit(1.5, 'cm'), #change legend key height
          legend.key.width = unit(1.5, 'cm'), #change legend key width
          legend.title = element_blank(), #change legend title font size
          legend.text = element_text(size=20))
}

    ggplot() +

    geom_half_boxplot(
        data =  t_data_subset,
        aes(x = cell_type, y = !!sym((r.i)) , fill = cell_type), outlier.color = NA, side = "l") +
    geom_half_point(
        data =  t_data_subset,
        aes(x = cell_type, y = !!sym((r.i)) , color = cell_type),transformation = position_quasirandom(width = 0.05), side = "r", size = 2) +
    
    scale_fill_manual(values = c(reds[8],reds[2],reds[4])) +
    scale_color_manual(values = c(reds[8],reds[2],reds[4])) +
    labs(x = "", y = "", title = r.i) +
    scale_x_discrete(limits = c("T1","WM","T2"), labels=c("TI-Pre", "WM","TII-Pre")) +
    scale_y_continuous(name="Beta Values", breaks = c(0,0.2,0.4,0.6,0.8,1), limits=c(0, 1)) +
    theme_clean()
}

## generate all CpG plots for desired probes

p <- lapply(names(t_data_subset[,1:4]), cpg.plot)


#####################################################################
## Generate plot panel
#####################################################################

lapply(seq_along(p), function(i) ggsave(
      plot = p[[i]], 
      filename = paste0(p[[i]]$labels$title,".pdf"),
      path = paste0("./methylome/figures")
    ))

pdf("./methylome/figures/combined-cpg-plots.pdf", height = 10, width = 10)
ggarrange(plotlist = p,
          common.legend = TRUE,
          legend = "bottom"
          ) 
dev.off()


