## Author: Andrew Palmer
## Title: Methylation Analysis Functions
## Date: 14-06-2025

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

####################################################################
## Prepare annotation function to add annotations to to DMP data ##
####################################################################

add.annotation.data <- function(df.x, df.y){
    
    anno_holding <- df.y[match(rownames(df.x), df.y$probeID),]

    new_data_frame <- cbind(df.x, anno_holding$CpG_chrm, anno_holding$CpG_beg, anno_holding$CpG_end, anno_holding$chromstat, anno_holding$cgi, anno_holding$epic_v1_loci, anno_holding$Methyl_450_loci, anno_holding$Methyl_27_loci)
    colnames(new_data_frame) <- c("GeneName","logFC","AveExp","t","p.value","adj.p.val","B", "logFC_beta", "t_beta", "AveExpr_beta", "Seqname","Start", "End", "ChromStat", "CGI", "EPICv1Loci", "Methyl450KLoci", "Methyl27KLoci")

    return(new_data_frame)

}

####################################################################
## Function to extract coefficient of choice and annotate them ##
####################################################################

extract.fit <- function(coef = coef, fit_m, fit_b, df.m, df.b, anno) {
    ## pull out DMPs identified in TI pre vs TII pre m_val
    DMPs_mval <- topTable(fit_m, num = nrow(df.m), coef = coef,
                          adjust.method = "BH", p.value = 1,
                          genelist = anno$genesUniq)

    ## pull out DMPs identified in TI pre vs TII pre beta_val
    DMPs_betaval <- topTable(fit_b, num = nrow(df.b), coef = coef ,
                             adjust.method = "BH", p.value = 1,
                             genelist = anno$genesUniq)

    ## match mval and beta val rows and add beta logFC, t and AveExpr results to m_val## table

    DMPs_betaval <- DMPs_betaval[match(rownames(DMPs_mval),
                                       rownames(DMPs_betaval)),]

    DMPs_mval$logFC_beta <- DMPs_betaval$logFC
    DMPs_mval$t_beta <- DMPs_betaval$t
    DMPs_mval$AveExpr_beta <- DMPs_betaval$AveExpr

    ## add annotation information

    DMPs_mval <- add.annotation.data(DMPs_mval, anno)

    DMPs_mval$ChromStat <- gsub("ChromHMM;","", DMPs_mval$ChromStat)

    DMPs_mval$CGI <- gsub("CGI;","", DMPs_mval$CGI)

    return(DMPs_mval)
}


