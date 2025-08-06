## Author: Andrew Palmer
## Title: proteomics plot functions
## Date: 2024-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
## See the 'LICENSE' file in the root directory for full license details.

####################################################################
## A script of proteomic plot functions ##
####################################################################

## This script contains proteomic plot functions that can be used
## for proteomic QC and data visualisations

####################################################################
## Colour-scheme Viridis ##
####################################################################

viridis_colours <- c("#440154B3", "#3C508BB3",  "#25838EB3",
                       "#31B57BB3", "#FDE725B3")

viridis_ramp <- colorRampPalette(viridis_colours[5:1])

viridis_6 <- c("#440154B3", "#414487B3", "#2A788EB3", "#22A884B3", "#7AD151B3", "#FDE725B3")
library(viridis)
viridis_colours <- viridis(n=6, alpha = 0.75, end =0.8)
library(RColorBrewer)
reds <- brewer.pal(9, name = "Reds")
blues <- brewer.pal(9, name = "Blues")
blues_ramp <- colorRampPalette(blues)
####################################################################
## Plot number of proteins in each sample ##
####################################################################

##protein.per.sample function, This plot can confirm if some samples
## have less proteins and should be exculded.

 cols <- viridis_colours[1:3]
plot.protein.per.sample <- function(df_x, df_y, cols = viridis_6) {

    ## Extract only MAXlfq columns
    intensity_names <- grep(".MaxLFQ", colnames(df_x), value = TRUE)
    holding <- df_x[, intensity_names]
    holding <- holding > 0
    data_to_hist <- apply(holding, 2, sum)
    data_to_hist <- as.data.frame(data_to_hist)
    df_y <- df_y[order(df_y$label),]
    rownames(data_to_hist) <- df_y$sample

    data_to_hist <- cbind(data_to_hist, df_y$sample_id)
    colnames(data_to_hist) <- c("protein_count", "sample_id")
    data_to_hist$group <- paste(df_y$fibre_type, df_y$condition, sep = "_")
    prot_count <- dim(df_x)
    prot_count <- prot_count[1]
         
    data_to_hist$group <- as.factor(data_to_hist$group)
    data_to_hist$sample_id <- as.factor(data_to_hist$sample_id)
    data_to_hist <- data_to_hist[order(data_to_hist$group),]
        
    barplot(data_to_hist$protein_count,
            beside = TRUE,
            col = cols[data_to_hist$group],
            las = 1,
            pch = 16,
            las = 2,
            main = "Proteins in Samples",
            ylab = "",
            xlab = "",
            ylim = c(0, 1400))

    legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA, title = "Sample Group", legend = levels(data_to_hist$group), pch = 15, col = cols, cex = 1)
    
    abline(h = prot_count, col = "#2db27d", lwd = 3)
    
    text(x = 18, y = 1300, paste("Total Overall Proteins", prot_count),
         col = "#2db27d", cex = 1.5)
    
    mtext("Sample ID", side = 1, line = 2, cex = 1.5)
    mtext("Protein Count", side = 2, line = 4.5, cex = 1.5)

}

cols <- viridis_colours[1:3]
cols <-  c(reds[8],reds[2],"grey")
plot.protein.per.sample.grouped <- function(df_x, df_y, cols = reds[8,6,2], grep.val = "pre") {

    ## Extract only MAXlfq columns
    intensity_names <- grep(".MaxLFQ", colnames(df_x), value = TRUE)
    holding <- df_x[, intensity_names]
    holding <- holding >  0
    data_to_hist <- apply(holding, 2, sum)
    data_to_hist <- as.data.frame(data_to_hist)
    df_y <- df_y[order(df_y$label),]
    rownames(data_to_hist) <- df_y$sample

    data_to_hist <- cbind(data_to_hist, df_y$sample_id)
    colnames(data_to_hist) <- c("protein_count", "sample_id")
        data_to_hist$group <- paste(df_y$fibre_type, df_y$condition, sep = "_")
    prot_count <- dim(df_x)
    prot_count <- prot_count[1]
    data_to_hist_subset <- data_to_hist[grepl(grep.val, data_to_hist$group),] 
       
    data_to_hist_subset$group <- as.factor(data_to_hist_subset$group)
    data_to_hist_subset$sample_id <- as.factor(data_to_hist_subset$sample_id)
    data_to_hist_subset <- data_to_hist_subset[order(data_to_hist_subset$group),]
    reds <- brewer.pal(8, name = "Reds")
  
    barplot(data_to_hist_subset$protein_count,
            beside = TRUE,
            col = cols[data_to_hist_subset$group],
            las = 1,
            pch = 16,
            las = 2,
            main = "",
            cex.main = 1.75,
            names.arg = data_to_hist_subset$sample_id,
            ylab = "",
            xlab = "",
            xaxt = "n",
            ylim = c(0, 1600))

       legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA, title = "Type", legend = c("TI", "TII", "TIIx"), pch = 15, col = cols, cex = 1.5)
    
    abline(h = prot_count, col = "#2db27d", lwd = 3)
    
    text(x = 18, y = 1300, paste("Total Identified Proteins", prot_count),
         col = "#2db27d", cex = 1.5)

    mtext("Protein Count", side = 2, line = 4.5, cex = 2.75)
}


####################################################################
## Heatmap of Missing Values ##
####################################################################

##heatmap.missingVal function
heatmap.missingVal <- function(df_x, df_y, cols = viridis_colours[c(5,1)]) {
  data_filtered <- df_x
  intensity_names <- grep(".MaxLFQ", colnames(df_x), value = TRUE)

  holding <- data_filtered[, intensity_names]
  holding <- holding > 0
  holding <- apply(holding, 2, as.integer)
  df_y <- df_y[order(df_y$label),]
  colnames(holding) <- df_y$sample


  heatmap(holding,
          col = cols,
          scale = "column",
          labRow = NA,
          cexRow = 0.2 + 1/log10(50),
          Colv = TRUE,
          Rowv = NA,
          xpd = TRUE,
          main = "Missing Value Pattern")

  legend("topright", inset = c(0, 0),
         legend = c("Missing Value", "Valid Value"),
         col = cols,
         pch = 16,
         cex = 1.5)
}


####################################################################
## Coeffiecient of Variation calculations ##
####################################################################

sample.coefficient.variation <- function(df_x, df_y, plot = FALSE) {

  if (plot == FALSE) {
    intensity_names <- grep(".MaxLFQ", colnames(df_x), value = TRUE)

    holding <- df_x[, intensity_names]
    holding[intensity_names] <- sapply(holding[intensity_names], as.numeric)
    holding[holding == 0] <- NA
    df_y <- df_y[order(df_y$label),]
    colnames(holding) <- df_y$sample
    holding <- t(holding)
    holding <- cbind.data.frame(holding, condition = df_y$condition)

    ## calculate combined means and CV
    sd_data_combined <- sapply(holding[, colnames(holding)
                                       [colnames(holding) != "condition"]],
                               sd, na.rm = TRUE)
    mean_data_combined <- colMeans(holding[sapply(holding, is.numeric)],
                                   na.rm = TRUE)
    cv_data_combined <- (sd_data_combined / mean_data_combined) * 100

    ## calculate SD
    sd_data <- aggregate(holding[, colnames(holding)
                                 [colnames(holding) != "condition"]],
                         list(holding$condition),
                         FUN = function(i) sd(i, na.rm = TRUE))

    row.names(sd_data) <- sd_data[, 1]
    sd_data <- sd_data[, -1]

    ## calculate mean
    mean_data <- aggregate(holding[, colnames(holding)
                                   [colnames(holding) != "condition"]],
                           list(holding$condition),
                           FUN = function(i) mean(i, na.rm = TRUE))

    row.names(mean_data) <- mean_data[, 1]
    mean_data <- mean_data[, -1]

    ## CV
    cv_data <- (sd_data / mean_data) * 100

    cv_data_t <- t(cv_data)
    cv_data_t <- as.data.frame(cv_data_t)
    cv_data_t <- data.frame(combined = cv_data_combined, cv_data_t)


    return(cv_data_t)
  }
}


sample.coefficient.variation.grouped <- function(df_x, df_y, plot = FALSE) {

  if (plot == FALSE) {
    intensity_names <- grep(".MaxLFQ", colnames(df_x), value = TRUE)

    holding <- df_x[, intensity_names]
    holding[intensity_names] <- sapply(holding[intensity_names], as.numeric)
    holding[holding == 0] <- NA
     df_y <- df_y[order(df_y$label),]
    colnames(holding) <- df_y$sample
    holding <- t(holding)
    holding <- cbind.data.frame(holding, condition = df_y$description)

    ## calculate combined means and CV
    sd_data_combined <- sapply(holding[, colnames(holding)
                                       [colnames(holding) != "condition"]],
                               sd, na.rm = TRUE)
    mean_data_combined <- colMeans(holding[sapply(holding, is.numeric)],
                                   na.rm = TRUE)
    cv_data_combined <- (sd_data_combined / mean_data_combined) * 100

    ## calculate SD
    sd_data <- aggregate(holding[, colnames(holding)
                                 [colnames(holding) != "condition"]],
                         list(holding$condition),
                         FUN = function(i) sd(i, na.rm = TRUE))

    row.names(sd_data) <- sd_data[, 1]
    sd_data <- sd_data[, -1]

    ## calculate mean
    mean_data <- aggregate(holding[, colnames(holding)
                                   [colnames(holding) != "condition"]],
                           list(holding$condition),
                           FUN = function(i) mean(i, na.rm = TRUE))

    row.names(mean_data) <- mean_data[, 1]
    mean_data <- mean_data[, -1]

    ## CV
    cv_data <- (sd_data / mean_data) * 100

    cv_data_t <- t(cv_data)
    cv_data_t <- as.data.frame(cv_data_t)
    cv_data_t <- data.frame(combined = cv_data_combined, cv_data_t)


    return(cv_data_t)
  }
}
