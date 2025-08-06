## Author: Andrew Palmer
## Title: Sex Differences DMR plot
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

#####################################################################
## Script to Generate Sex Differences DMR plots
#####################################################################

#####################################################################
## Load Packages
#####################################################################

library(data.table)
library(ggvenn)
library(RColorBrewer)
library(jsonlite)

dmrs_tI <-  data.frame(fread("./results/final-dmrcate-dmrs-tI_FvM.csv"))

dmrs_tII <-  data.frame(fread("./results/final-dmrcate-dmrs-tII_FvM.csv"))

dmrs_wm <-  data.frame(fread("./results/final-dmrcate-dmrs-wm_FvM.csv"))

dmrs_tI <- dmrs_tI[dmrs_tI$no.cpgs > 3 & dmrs_tI$Fisher < 0.001,]
dmrs_tII <- dmrs_tII[dmrs_tII$no.cpgs > 3  & dmrs_tII$Fisher < 0.001,]
dmrs_wm <- dmrs_wm[dmrs_wm$no.cpgs > 3  & dmrs_wm$Fisher < 0.001,]

dmrs_list <- list("TI FvM DMRs" = unique(dmrs_tI$symbol),"TII FvM DMRs" = unique(dmrs_tII$symbol), "WM FvM DMRs" = unique(dmrs_wm$symbol))

venn <- ggvenn(
  dmrs_list, 
  fill_color = brewer.pal(3, name = "Blues"),
  stroke_size = 1.5, set_name_size = 6, text_size = 6
)

ggsave(
    "./figures/sex_dmrs_overlaps.pdf",
    plot = venn,
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


overlapping_dmrs <- calculate.overlap(dmrs_list)
names(overlapping_dmrs) <- c("TI&TII&WM", "TI&TII","TI&WM", "TII&WM", "TI", "TII","WM")


Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

Union <- function (x) {  
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's. 
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}


dmrs_tI$region <- paste0(dmrs_tI$seq,":", dmrs_tI$start,dmrs_tI$end)
dmrs_tII$region <- paste0(dmrs_tII$seq,":", dmrs_tII$start,dmrs_tII$end)
dmrs_wm$region <- paste0(dmrs_wm$seq, ":",dmrs_wm$start,dmrs_wm$end)

dmrs_list2 <- list("TI FvM DMRs" = unique(dmrs_tI$region),"TII FvM DMRs" = unique(dmrs_tII$region), "WM FvM DMRs" = unique(dmrs_wm$region))

json_data <- toJSON(overlapping_dmrs)
write(json_data, "./results/sex-overlaps.json")
