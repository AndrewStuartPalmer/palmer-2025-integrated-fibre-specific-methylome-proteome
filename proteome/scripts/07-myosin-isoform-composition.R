## Author: Andrew Palmer
## Title: Myosin isoform composition
## Date: 2024-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
## See the 'LICENSE' file in the root directory for full license details.

####################################################################
## Myosin isoform content of proteomic samples.
####################################################################

## This code graphs the myosin isoform content of proteomic samples.
## It takes raw Unique intensity files (generated with msfragger).
## It outputs percentages of MYH7, MYH2, MYH1 and MYH4 of samples.

####################################################################
## Load Data
####################################################################
library(RColorBrewer)
reds <- brewer.pal(9, name = "Reds")
blues <- brewer.pal(9, name = "Blues")

folder_path <- getwd()

files <- list.files(folder_path, pattern = "protein.tsv", full.names = TRUE, recursive = TRUE)

experiment_design <- read.csv("~/Documents/pooled-fibre-type-project/lfq_analyst_experimental_design.csv", sep="\t")

####################################################################
## Generate Myosin plot data
####################################################################

plot_data <- data.frame()

for (file in files) {

data <- read.delim(file)

data_filt <- data[,c("Gene","Unique.Intensity")]

myh <- c("MYH1","MYH7","MYH2","MYH4")

data_subset <- subset(data_filt, Gene %in% myh)

data_subset$Unique_percent <- (data_subset$Unique.Intensity / sum(data_subset$Unique.Intensity)) * 100

base_name <- basename(dirname(file))

    data_subset$sample<- unlist(strsplit(base_name, "_"))[length(unlist(strsplit(base_name, "_")))]

    
    plot_data <- rbind(plot_data, data_subset)
    
}

####################################################################
## Prepare data for calculating averages for each fibre type
####################################################################

## reorder factor levels for plotting
experiment_design <- experiment_design[,c(5,8:9)]

merged_df <- merge(plot_data, experiment_design, by.x = "sample", by.y = "sample_no")

merged_df <-  subset(merged_df, condition == "pre")

merged_df$Gene <- factor(merged_df$Gene, levels = c("MYH7", "MYH2", "MYH1", "MYH4"))

## pull out only pre samples

averages <- aggregate(cbind(Unique_percent) ~ Gene + fibre_type, 
                      data = merged_df, 
                      FUN = mean, 
                      na.rm = TRUE)

####################################################################
## Convert data to wide format for plotting
####################################################################

merged_df_small <- merged_df[,c(1,2,4)]

wide_data <- reshape(merged_df_small, 
                     idvar = "Gene", 
                     timevar = "sample", 
                     direction = "wide")

colnames(wide_data) <- gsub("Unique_percent.", "", colnames(wide_data))

rownames(wide_data) <- wide_data$Gene
wide_data <- wide_data[,-1]
wide_data <- as.matrix(wide_data)
H <- apply(wide_data, 2L, cumsum)
H <- H - wide_data / 2

wide_data <- wide_data[,c(1,3,5,7,9,11,13,15,17,19,21,23,2,4,6,8,10,12,14,16,18,20,22,24,25,26)]

wide_data <- wide_data[c(2,3,1,4),]

####################################################################
## Generate stacked bar plot
####################################################################
    
pdf("myosin_isoforms_barplot.pdf", height = 20, width = 20)
par(mar = c(15, 6, 6, 6), cex.main = 3.5, family = "Helvetica")
x <- barplot(wide_data,
        col = blues[2:5],
        border = "white",
        ylab = "",
        yaxt = "n",
        xaxt = "n",
        xlim = c(0,100),
        las = 2,
        cex.names = 3,
        horiz = TRUE,
        main = "Percentage Myosin Isoforms")

axis(side = 1, col = NA, cex.axis = 3, las = 1)
abline(v = c(0,20,40,60,80,100), lty = 3, lwd = 4, col = "lightgrey")
text(round(H, digits = 1), rep(x, each = nrow(H)), labels = ifelse(wide_data > 5, paste0(round(wide_data, digits = 1), "%"), ""), cex =2.75, col = "black")

rect(-3, 0.2, -1, 14.3, col = reds[8], border = reds[8], xpd = TRUE)
rect(-3, 14.6, -1, 28.8, col = reds[2], border = reds[2], xpd = TRUE)
rect(-3, 29.1, -1, 31.0, col = "grey", border = "grey", xpd = TRUE)
text(-5, 7, label = "TI", xpd = TRUE, srt = 90, cex = 4)
text(-5, 22, label = "TII", xpd = TRUE, srt = 90, cex = 4)
text(-5, 30, label = "TIIx", xpd = TRUE, srt = 90, cex = 4)

legend(x = "bottom",
       inset = c(0,-0.15),
       legend = rownames(wide_data), 
       lty = 1,
       col = blues[2:5],
       lwd = 9,
       cex = 3.5,
       xpd = TRUE,
       horiz = TRUE)

dev.off()

write.csv(wide_data, "myosin_unique_content.csv")
write.csv(merged_df, "myosin_unique_content_2.csv")
write.csv(averages, "myosin_unique_fibre_averages.csv")
