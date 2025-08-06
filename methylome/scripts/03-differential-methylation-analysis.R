## Author: Andrew Palmer
## Title: Differential methylation analysis 
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Script for differential methylation analysis ##
#####################################################################

####################################################################
## Load custom methylation analysis functions ##
####################################################################

source("./methylome/scripts/util-functions/methyl-analysis-functions.R")

####################################################################
## Load required packages ##
####################################################################

load_packages <- c("data.table", "limma", "RColorBrewer",
                   "IlluminaHumanMethylationEPICv2anno.20a1.hg38", "pbapply")

lapply(load_packages, library, character.only = TRUE)

####################################################################
## Load Data ##
####################################################################
m_val_filt_subset <-  data.frame(data.table::fread("./methylome/data/processed-data/mval-filt-subset.csv"), row.names = 1)
beta_filt_subset  <- data.frame(data.table::fread("./methylome/data/processed-data/beta-filt-subset.csv"),row.names = 1)
pheno_filt_subset <-  data.frame(data.table::fread("./methylome/experiment-design/pheno-table-subset.csv"))

####################################################################
## Set up linear model M-Vals ##
####################################################################

## using contrasts matrix and setting up a group
## set up model matrix on group without an intercept term
design <- model.matrix(~ 0 + group + sex, pheno_filt_subset)

## each participants contributes 7 samples to the study
## we account for this with duplicate correlation on participant ID

cor <- duplicateCorrelation(m_val_filt_subset, design, block = pheno_filt_subset$participant)

## fit the model blocking on participant

fit <- lmFit(object = m_val_filt_subset,
             design = design,
             block = pheno_filt_subset$participant,
             correlation = cor$consensus.correlation)

## contrasts of interest are set to enable comparisons

contrasts <- makeContrasts(T1preVsT2pre = groupT1_Pre - groupT2_Pre,
                           T1preVsWMpre = groupT1_Pre - groupWM_Pre,
                           T2preVsWMpre = groupT2_Pre - groupWM_Pre,
                           FTVsWM = ((groupT1_Pre + groupT2_Pre) / 2) - groupWM_Pre,
                           levels = colnames(design))

fit_2 <- contrasts.fit(fit, contrasts)

fit_3 <- eBayes(fit_2)

####################################################################
## how to caluculate cohensd note this is a ebayes equivalent
####################################################################

cohensd <- fit_3$coefficients / sqrt(fit_3$s2.post)

####################################################################
## Save decideTests summaries ##
####################################################################

decide_tests_summary <- summary(decideTests(fit_3, adjust.method = "BH", p.value = 0.001))

decide_tests_summary_0.05 <- summary(decideTests(fit_3, adjust.method = "BH", p.value = 0.05))

decide_tests_summary_0.005 <- summary(decideTests(fit_3, adjust.method = "BH", p.value = 0.005))

####################################################################
## Set up linear model beta-vals ##
####################################################################

## run same thing with beta vals
## using contrasts matrix and setting up a group
## set up model matrix on group without an intercept term
## each participants contributes 7 samples to the study
## we account for this with duplicate correlation on participant ID

cor_beta <- duplicateCorrelation(beta_filt_subset, design, block = pheno_filt_subset$participant)

## fit the model blocking on participant

fit_beta <- lmFit(object = beta_filt_subset,
             design = design,
             block = pheno_filt_subset$participant,
             correlation = cor_beta$consensus.correlation)

fit_2_beta <- contrasts.fit(fit_beta, contrasts)

fit_3_beta <- eBayes(fit_2_beta)


####################################################################
## Extract and Save all Coefficients of interest ##
####################################################################

annotation_epicv2 <- read.csv("./methylome/annotations/annotation-epic-v2.txt", sep="\t")

anno_epic_v2_sub <- annotation_epicv2[match(rownames(m_val_filt_subset), annotation_epicv2$probeID),]

contrasts <- colnames(decide_tests_summary)
list_contrasts <- as.list(1:length(contrasts))
names(list_contrasts) <- colnames(decide_tests_summary)

DMPs <- pbapply::pblapply(list_contrasts, extract.fit, fit_3, fit_3_beta, m_val_filt_subset, beta_filt_subset, anno_epic_v2_sub, cl = 4)

pbapply::pbmapply(function (x,y) write.csv(x, file = paste0('./methylome/results/final-dasen-dmps-',y, '.csv'), row.names = T), DMPs, names(DMPs))

####################################################################
## P value histogram  ##
####################################################################

## overall fit-p-values of Coefficient 1

## T1pre vs T2pre pval hist
pdf("./methylome/figures/final-p-value-hist-t1vst2-meth.pdf", width = 12, height = 8)
par(cex= 1.15, cex.main=1.75, cex.lab=1.5, cex.axis=1.5, mar = c(4.5,4.5,1.5,1), family = "Helvetica")
hist(fit_3$p.value[,1],
     ylim = c(0,300000),
     xlab = "P-values",
     main = "")
dev.off()

####################################################################
## Save session info ##
####################################################################

writeLines(capture.output(sessionInfo()),
           "./methylome/results/session-info/methylation-analysis-sessionInfo.txt")



