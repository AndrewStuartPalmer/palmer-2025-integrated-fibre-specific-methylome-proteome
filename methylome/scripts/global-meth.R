## Author: Andrew Palmer
## Title: Global Methylation Analysis
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

####################################################################
## Global Methylation Analysis.
####################################################################

####################################################################
## Load packages
####################################################################

load_packages <- c("data.table", "viridis")

lapply(load_packages, library, character.only = TRUE)

m_val <-  data.frame(fread("./methylome/data/processed-data/full-dasen-normalised-filtered-m.csv"), row.names = 1)

####################################################################
## Filter and prepare samples
####################################################################

pheno <- read.delim("./methylome/experiment-design/pheno-table.txt")

## combine cell_type and timepoint for limma analysis

pheno$group <- paste(pheno$cell_type, pheno$timepoint, sep = "_")

## take average m values for replicate samples
m_val_rep <- m_val[,c(17,18,19,90,91)]
m_val$WM.f28.preWM.1 <- rowMeans(m_val_rep[,4:5])
m_val$WM.s201.preWM.1 <- rowMeans(m_val_rep[,1:3])


## remove replicate samples and keep only average m values
pheno_filt <- pheno[-c(1,2,17,18,30,57,58,90),] 

m_val_filt <- m_val[,- c(1,2,17,18,30,57,58,90)]
row.names(pheno_filt) <- 1:nrow(pheno_filt)

## same fore beta values
beta <-  data.frame(fread("./methylome/data/processed-data/full-dasen-normalised-filtered-beta.csv"), row.names =1)
beta_rep <- beta[,c(17,18,19,90,91)]
beta$WM.f28.preWM.1 <- rowMeans(beta_rep[,4:5])
beta$WM.s201.preWM.1 <- rowMeans(beta_rep[,1:3])

beta_filt <- beta[,- c(1,2,17,18,30,57,58,90)]


## remove 4wWM samples

subset_position <- !endsWith(colnames(m_val_filt), "4wWM")
m_val_filt_subset <- subset(m_val_filt, select = subset_position)
pheno_filt_subset <- subset(pheno_filt, subset = subset_position)
beta_filt_subset <- subset(beta_filt, select = subset_position)

pheno_filt_subset <- pheno_filt_subset[-c(57),] 

m_val_filt_subset  <- m_val_filt_subset[,- c(57)]
beta_filt_subset <- beta_filt_subset[,- c(57)]

####################################################################
## Calculate overal methylation values.
####################################################################


mean_beta <- as.data.frame(colMeans(beta_filt_subset))
mean_beta$cell_type <- pheno_filt_subset$cell_type
colnames(mean_beta) <- c("beta", "cell_type")
mean_beta$beta_percent <- mean_beta$beta * 100 
globalmeth <- aggregate(mean_beta$beta, list(mean_beta$cell_type), FUN = mean)
colnames(globalmeth) <- c("CellType","globalmeth")


with(mean_beta,
     colSums(table(beta, cell_type)))

celltype.lm <- lm(formula = beta ~ cell_type,
                  data = mean_beta)
cell_type.I.aov <- aov(celltype.lm)

aov.model <- aov(beta ~ cell_type, data = mean_beta)
TukeyHSD(aov.model)

boxplot(mean_beta$beta_percent ~ mean_beta$cell_type,
        ylim = c(54,55))

aov_globalmeth <- aov(mean_beta$beta ~ mean_beta$cell_type, data = mean_beta)

mean_beta_t2 <- mean_beta[mean_beta$cell_type == "T2",]
mean_beta_t1 <- mean_beta[mean_beta$cell_type == "T1",]
mean_beta_wm <- mean_beta[mean_beta$cell_type == "WM",]


t1vst2_global <- t.test((mean_beta_t1$beta), (mean_beta_t2$beta))
t1vswm_global <- t.test(mean_beta_t1$beta, mean_beta_wm$beta)
