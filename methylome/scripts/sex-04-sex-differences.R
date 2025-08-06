## Author: Andrew Palmer
## Title: Sex differences overlaps 
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Script for comparing overlaps of sex differences in DNAm
#####################################################################

#####################################################################
## Load Packages
#####################################################################
load_packages <- c("data.table",
                   "ggplot2",
                   "ggvenn",
                   "RColorBrewer",
                   "ggplotify",
                   "ggrepel")

lapply(load_packages, library, character.only = TRUE)

#####################################################################
## Load Data
#####################################################################

tI_f_mdmps <- data.frame(fread("./results/dasen-dmps-T1preFVsTpreM.csv"), 
                         row.names = 1)

tII_f_m_dmps <- data.frame(fread("./results/dasen-dmps-T2preFVsT2preM.csv"), 
                           row.names = 1)

t_wm_f_m_dmps <- data.frame(fread("./results/dasen-dmps-WMpreFVsWMpreM.csv"), 
                            row.names = 1)

#####################################################################
## 
#####################################################################

                   
sig_DMPs_hyper_TI <- tI_f_mdmps[tI_f_mdmps$logFC_beta > 0.1 & tI_f_mdmps$adj.p.val <0.001,]
sig_DMPs_hypo_TI <- tI_f_mdmps[tI_f_mdmps$logFC_beta < -0.1 & tI_f_mdmps$adj.p.val <0.001,]

sig_DMPs_hyper_TII <- tII_f_m_dmps[tII_f_m_dmps$logFC_beta > 0.1 & tII_f_m_dmps$adj.p.val <0.001,]
sig_DMPs_hypo_TII <- tII_f_m_dmps[tII_f_m_dmps$logFC_beta < -0.1 & tII_f_m_dmps$adj.p.val <0.001,]

sig_DMPs_hyper_wm <- t_wm_f_m_dmps[t_wm_f_m_dmps$logFC_beta > 0.1 & t_wm_f_m_dmps$adj.p.val <0.001,]
sig_DMPs_hypo_wm <- t_wm_f_m_dmps[t_wm_f_m_dmps$logFC_beta < -0.1 & t_wm_f_m_dmps$adj.p.val <0.001,]


#####################################################################
## 
#####################################################################


v1 <- list("Hyper-T1" = rownames(sig_DMPs_hyper_TI), "Hyper-T2" = rownames(sig_DMPs_hyper_TII), "Hypo-T1" = rownames(sig_DMPs_hypo_TI),  "Hypo-T2" = rownames(sig_DMPs_hypo_TII))

v2 <- list("Hyper-T1" = rownames(sig_DMPs_hyper_TI), "Hyper-T2" = rownames(sig_DMPs_hyper_TII), "Hyper-WM" = rownames(sig_DMPs_hyper_wm))

v3 <- list("Hypo-T1" = rownames(sig_DMPs_hypo_TI),  "Hypo-T2" = rownames(sig_DMPs_hypo_TII), "Hypo-WM" = rownames(sig_DMPs_hypo_wm))

#####################################################################
## 
#####################################################################

venn <- ggvenn(
  v1, 
  fill_color = brewer.pal(5, name = "Blues"),
  stroke_size = 1.5, set_name_size = 7, text_size = 7
)

ggvenn(
  v2, 
  fill_color = brewer.pal(5, name = "Blues"),
  stroke_size = 1.5, set_name_size = 7, text_size = 7
)

ggvenn(
  v3, 
  fill_color = brewer.pal(5, name = "Blues"),
  stroke_size = 1.5, set_name_size = 7, text_size = 7
)

hyper_tII_only <- sig_DMPs_hyper_TII[!(rownames(sig_DMPs_hyper_TII) %in% rownames(sig_DMPs_hyper_TI)),]

hyper_tII_only <- hyper_tII_only[!(rownames(hyper_tII_only) %in% rownames(sig_DMPs_hyper_wm)),]


