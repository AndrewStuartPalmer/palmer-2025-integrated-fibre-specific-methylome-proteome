## Author: Andrew Palmer
## Title: Proteomics differential protein analysis Sex
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025
# See the 'LICENSE' file in the root directory for full license details.

####################################################################
## Differentially expressed proteins analysis between the sexes
####################################################################

# It takes a normalised and imputed data frame of LFQ intensities .
# It performs regression analysis using limma comparing sex.
# It outputs a results file with log2FC and adjusted p value for 
# each protein tested.

source("./proteome/scripts/proteomics-functions.R")

load_packages <- c("limma", "org.Hs.eg.db", "pbapply")

lapply(load_packages, library, character.only = TRUE)

####################################################################
## Load Data
####################################################################

## We need to load the imputed data frame and the
## experimental design file to run limma.

data_normalised_imputed <- read.delim("./proteome/data/processed-data/data_normalised_imputed.tsv",
                                      header = TRUE,
                                      sep = "\t",
                                      row.names = 1)

experiment_design <- read.delim("./proteome/experiment-design/lfq_analyst_experimental_design.txt")



experiment_design <- experiment_design[order(experiment_design$label), ]
colnames(data_normalised_imputed) <- experiment_design$sample

experiment_design$group <- paste(experiment_design$fibre_type, experiment_design$condition, sep = "_")

experiment_design$sex_groups <-  paste(experiment_design$group, experiment_design$sex, sep = "_")

data_normalised_imputed <- data_normalised_imputed[,-c(47,48,49,50,51)]
experiment_design <- experiment_design[-c(47,48,49,50),]
data_normalised_imputed <- data_normalised_imputed[,-c(40)]
experiment_design <- experiment_design[-c(40),]

design <- model.matrix(~ 0 + sex_groups, experiment_design)

colnames(design) <- c("TI_post_F", "TI_post_M", "TI_pre_F", "TI_pre_M", "TII_post_F", "TII_post_M", "TII_pre_F", "TII_pre_M")

cor <- duplicateCorrelation(data_normalised_imputed, design, block = experiment_design$sample_id)

fit <- lmFit(object = data_normalised_imputed,
             design = design,
             block = experiment_design$sample_id,
             correlation = cor$consensus.correlation)

contrasts <- makeContrasts(TIpreFVsTIpreM = TI_pre_F - TI_pre_M,
                           TIIpreFVsTIIpreM = TII_pre_F - TII_pre_M,
                           TIFvsTIM = (TI_pre_F + TI_post_F) -  (TI_pre_M + TI_post_M),
                             TIIFvsTIIM = (TII_pre_F + TII_post_F) -  (TII_pre_M + TII_post_M),
                         interaction = (TI_pre_F - TI_pre_M) -  ( TII_pre_F - TII_pre_M),                           levels = colnames(design))

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

decide_tests_summary_0.001 <- summary(decideTests(fit2, adjust.method = "BH", p.value = 0.001))

decide_tests_summary_0.005 <- summary(decideTests(fit2, adjust.method = "BH", p.value = 0.005))

decide_tests_summary_0.05 <- summary(decideTests(fit2, adjust.method = "BH", p.value = 0.05))


contrasts <- colnames(decide_tests_summary_0.001)
list_contrasts <- as.list(1:length(contrasts))
names(list_contrasts) <- colnames(decide_tests_summary_0.001)

deps <- pbapply::pblapply(list_contrasts, extract.fit.proteome, fit2, data_normalised_imputed, cl = 4)

pbapply::pbmapply(function (x,y) write.csv(x, file = paste0('./proteome/results/final-sex-deps-',y, '.csv'), row.names = T), deps, names(deps))

pdf("./proteome/figures/p-val-hist-sex-prot.pdf", width = 12, height = 8)
par(cex= 1.15, cex.main=1.75, cex.lab=1.5, cex.axis=1.5, mar = c(4.5,4.5,1.5,1), family = "Helvetica")
hist(fit2$p.value[,2],
     ylim = c(0,250),
     xlab = "P-values",
     main = "")
dev.off()

writeLines(capture.output(sessionInfo()),
           "./proteome/scripts/session-info/final-protein-sex-analysis-sessionInfo.txt")
