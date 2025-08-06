## Author: Andrew Palmer
## Title: Proteomics differential protein analysis
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025
# See the 'LICENSE' file in the root directory for full license details.

####################################################################
## Differentially expressed proteins analysis
####################################################################

# This script is the fifth step in the proteomic analysis workflow.
# It takes a normalised and imputed data frame of LFQ intensities .
# It performs regression analysis using limma.
# It outputs a results file with log2FC and adjusted p value for 
# each protein tested.

# Source custom plot functions from 'proteomics-plot-functions.R'

source("./proteome/scripts/proteomics-plot-functions.R")

## Protein-wise linear models combined with empirical Bayes statistics
## were used for the differential expression analyses

####################################################################
## Load Limma Package
####################################################################

## Load ggplot2 ggrepel for generating plots
load_packages <- c("limma", "org.Hs.eg.db")

lapply(load_packages, library, character.only = TRUE)

####################################################################
## Load Data
####################################################################

## We need to load the imputed data frame and the
## experimental design file to run limma.

data_normalised_imputed <-
  read.delim(
    "./proteome/data/processed-data/data_normalised_imputed.tsv",
    header = TRUE,
    sep = "\t",
    row.names = 1
  )

experiment_design <-
  read.delim("./proteome/experiment-design/lfq_analyst_experimental_design.txt")

experiment_design <-
  experiment_design[order(experiment_design$label),]
colnames(data_normalised_imputed) <- experiment_design$sample

experiment_design$group <-
  paste(experiment_design$fibre_type,
        experiment_design$condition,
        sep = "_")

####################################################################
## Remove T2x samples
####################################################################
data_normalised_imputed <-
  data_normalised_imputed[, -c(47, 48, 49, 50, 51)]
experiment_design <- experiment_design[-c(47, 48, 49, 50), ]
data_normalised_imputed <- data_normalised_imputed[, -c(40)]
experiment_design <- experiment_design[-c(40), ]

####################################################################
## Make model matrix
####################################################################

design <- model.matrix( ~ 0 + group + sex, experiment_design)

colnames(design) <-
  c("TI_post", "TI_pre", "TII_post", "TII_pre", "sexm")

cor <- duplicateCorrelation(data_normalised_imputed, design, block = experiment_design$sample_id)


fit <- lmFit(
  object = data_normalised_imputed,
  design = design,
  block = experiment_design$sample_id,
  correlation = cor$consensus.correlation
)

####################################################################
## Make contrasts
####################################################################

contrasts <- makeContrasts(
  TIpreVsTIIpre = TI_pre - TII_pre,
  TIpostVsTIIpost = TI_post - TII_post,
  TIpostVsTIpre = TI_post - TI_pre,
  TIIpostVsTIIpre = TII_post - TII_pre,
  interaction = (TII_post - TII_pre) -  (TI_post - TI_pre),
  levels = colnames(design)
)
####################################################################
## Run limma
####################################################################

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)


decide_tests_summary_0.001 <-
  summary(decideTests(fit2, adjust.method = "BH", p.value = 0.001))

decide_tests_summary_0.005 <-
  summary(decideTests(fit2, adjust.method = "BH", p.value = 0.005))

decide_tests_summary_0.05 <-
  summary(decideTests(fit2, adjust.method = "BH", p.value = 0.05))


####################################################################
## Save T1 vs T2 comparisons
###################################################################

## comparison T1 vs T2

diff_1 <-
  topTable(fit2,
           coef = 1,
           number = 1226,
           adjust.method = "BH")

diff_1$symbol <-
  AnnotationDbi::mapIds(org.Hs.eg.db,
                        rownames(diff_1),
                        keytype =  "SYMBOL",
                        column = "SYMBOL")

diff_1$uniprot <-
  AnnotationDbi::mapIds(org.Hs.eg.db,
                        rownames(diff_1),
                        keytype =  "SYMBOL",
                        column = "UNIPROT")

write.table(diff_1,
            file = "./proteome/results/t1vst2-pre-proteome.csv",
            sep = ",",
            col.names = NA)


###################################################################
## P value histogram T1 vs T2
###################################################################

## T1pre vs T2pre pval hist

pdf("./proteome/figures/p-value-hist-t1vst2-prot.pdf",
    width = 8,
    height = 8)
par(
  cex = 1.15,
  cex.main = 1.75,
  cex.lab = 1.5,
  cex.axis = 1.5,
  mar = c(4.5, 4.5, 1.5, 1),
  family = "Helvetica"
)
hist(
  fit2$p.value[, 1],
  ylim = c(0, 1000),
  xlab = "P-values",
  main = ""
)
dev.off()

###################################################################
## Volcano Plot
###################################################################

## prepare data for volcano plot

logfc <- diff_1$logFC
pval <- diff_1$P.Value
pval <- -log10(pval)
vol_df_2 <- cbind(logfc, pval)
vol_df_2 <- as.data.frame(vol_df_2)
rownames(vol_df_2) <- rownames(diff_1)
down_reg <- subset(vol_df_2, logfc < -1 & pval > -log10(0.001))
up_reg <- subset(vol_df_2, logfc > 1 & pval > -log10(0.001))


## generate groups for labels

unique <- read.csv("./proteome/unique.csv", row.names = 1)
intersection <-
  read.csv("./proteome/intersection.csv", row.names = 1)
unique_index <- rownames(vol_df_2) %in% unique$x
intersection_index <- rownames(vol_df_2) %in% intersection$x
unique_df <- vol_df_2[unique_index, ]
unique_df_up <- unique_df[unique_df$logfc > 0, ]
unique_df_down <- unique_df[unique_df$logfc < 0, ]
intersection_df <- vol_df_2[intersection_index, ]
intersection_df_up <- intersection_df[intersection_df$logfc > 0, ]
intersection_df_down <- intersection_df[intersection_df$logfc < 0, ]

cols <- viridis_colours
cols2 <- hcl.colors(5, palette = "cold")

## we highlight the proteins that were also reported in the following papers
## we also highlight the different proteins we identified
## ref 1 https://doi.org/10.1038/s41467-020-20556-8
## ref 2 https://doi.org/10.1038/s41467-024-50632-2

sig <-
  vol_df_2[vol_df_2$pval > -log10(0.001) &
             vol_df_2$logfc > 0 |
             vol_df_2$pval > -log10(0.001) & vol_df_2$logfc < 0, ]

## generate volcano plot
df <- data.frame(x = c(0),
                 y = c(35),
                 text = c("TI vs TII"))

vol_plot <- ggplot(vol_df_2, aes(x = logfc,
                                 y = pval)) +
  geom_point(
    data = vol_df_2,
    shape = 21,
    size = 2,
    fill = "grey",
    colour = "grey"
  ) +
  geom_point(
    data = sig,
    shape = 21,
    size = 2,
    fill = blues[2],
    colour = blues[2]
  ) +
  geom_point(
    data = up_reg,
    shape = 21,
    size = 2,
    fill = blues[8],
    colour = blues[8]
  ) +
  geom_point(
    data = down_reg,
    shape = 21,
    size = 2,
    fill = blues[6],
    colour = blues[6]
  ) +
  
  geom_hline(
    yintercept = -log10(0.001),
    linetype = "dashed",
    col = "red"
  ) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed",
             col = "red") +
  labs(x = "LogFC", y = "-log10(Pval)") +
  ylim(c(0, 35)) +
  
  xlim(c(-7, 7)) +
  
  theme_bw() + # Select theme with a white background
  theme(
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      color = "grey",
      size = 0.75,
      linetype = 2
    ),
    axis.title.x = element_text(
      size = 20,
      color = "black",
      family = "Helvetica"
    ),
    axis.title.y = element_text(
      size = 20,
      color = "black",
      family = "Helvetica"
    ),
    axis.text.x = element_text(
      color = "black",
      family = "Helvetica",
      size = 20
    ),
    axis.text.y = element_text(
      color = "black",
      family = "Helvetica",
      size = 20
    )
  ) +
  
  geom_label_repel(
    data = unique_df_up[1:10, ],
    aes(label = rownames(unique_df_up[1:10, ])),
    size = 5,
    fill = cols2[1],
    min.segment.length = 0,
    seed = 42,
    box.padding = 0.5,
    max.overlaps = Inf,
    inherit.aes = TRUE
  ) +
  
  geom_label_repel(
    data = unique_df_down[1:6, ],
    aes(label = rownames(unique_df_down[1:6, ])),
    size = 5,
    fill = cols2[1],
    min.segment.length = 0,
    seed = 42,
    box.padding = 0.5,
    max.overlaps = Inf,
    inherit.aes = TRUE
  ) +
  
  geom_label_repel(
    data = intersection_df_up[1:11, ],
    aes(label = rownames(intersection_df_up[1:11, ])),
    size = 5,
    fill = cols2[4],
    min.segment.length = 0,
    seed = 42,
    box.padding = 0.5,
    max.overlaps = Inf,
    inherit.aes = TRUE
  ) +
  
  geom_label_repel(
    data = intersection_df_down[1:11, ],
    aes(label = rownames(intersection_df_down[1:11, ])),
    size = 5,
    fill = cols2[4],
    min.segment.length = 0,
    seed = 42,
    box.padding = 0.5,
    max.overlaps = Inf,
    inherit.aes = TRUE
  ) +
  
  geom_label(
    data = df,
    aes(x = x, y = y, label = text),
    color = "black",
    family = "Helvetica",
    size = 10
  )

# save volcano plot
ggsave(
  "./proteome/figures/volcano-t1-pre-t2pre-deps.pdf",
  plot = vol_plot,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 12,
  height = 12,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)


