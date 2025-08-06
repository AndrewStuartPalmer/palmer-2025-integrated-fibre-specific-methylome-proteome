## Author: Andrew Palmer
## Title: Methylation EPICv2 filter functions 
## Date: 14-06-2025

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

####################################################################
## Load required pacakges ##
####################################################################

load_packages <- c("data.table",
                   "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
                   "minfi", "wateRmelon")

lapply(load_packages, library, character.only = TRUE)

## remove sex chromosomes before normalisations

filter.Sex <- function(m_set){
    anno_epic_v2 <- minfi::getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
    chromosomes_keep <- !(featureNames(m_set) %in%
                          anno_epic_v2$Name[anno_epic_v2$chr %in%
                                            c("chrX", "chrY")])
    table(chromosomes_keep)
    chromosome_probes <- table(chromosomes_keep)
    m_set <- m_set[chromosomes_keep, ]
    
    message(paste0(sum(chromosomes_keep == FALSE)," sex chromosome probes removed :",   sum(chromosomes_keep == TRUE), " probes remaining..."))

    m_set
}

filter.DetP <- function(m_set, rg_set) {
## make sure samples in same order in det_p object as g_set_quant
    ## before filtering
    message("Calculating detection P Values...")
    det_p <- minfi::detectionP(rg_set)
     message("Finished calculating detection P Values...")
    det_p <- det_p[match(featureNames(m_set), rownames(det_p)), ]
    
## remove probes with low detection p values
    
    good_probes <- rowSums(det_p < 0.01) == ncol(m_set)

    m_set <- m_set[good_probes, ]
    
    message(paste0(sum(good_probes == FALSE)," probes with < 0.01 DetP removed: ",
                     sum(good_probes == TRUE), " probes remaining..."))

    m_set
}

filter.LowBead <- function(m_set, rg_set, bead.number = 3) {

    message("obtain low bead count < ", bead.number," ...")

    low_b <- getNBeads(rg_set) < bead.number

    pi1 <- getProbeInfo(rg_set, type = "I")
    pi2 <- getProbeInfo(rg_set, type = "II")
    ex1 <- pi1$Name[rowMeans(low_b[pi1$AddressA, ] | low_b[pi1$AddressB, ]) > 0.05]
    ex2 <- pi2$Name[rowMeans(low_b[pi2$AddressA, ]) > 0.05]
    exclude_bds <- unique(c(ex1, ex2))
    length(exclude_bds)

## remove probes with low bead count

    probes_keep <- !(featureNames(m_set) %in% exclude_bds)

    m_set <- m_set[probes_keep, ]
    
    message(paste0(sum(probes_keep == FALSE)," probes with < 3 beads removed: ",
                   sum(probes_keep == TRUE), " probes remaining..."))
m_set
}

filter.NonCpg <- function(m_set) {
    anno_epic_v2 <- minfi::getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
    ## filter non cpg probes
    anno_epic_v2_noncpg <- anno_epic_v2[grepl("ch", anno_epic_v2$Name), ]

    exclude_cpgs <- anno_epic_v2_noncpg$Name

    cpgs_keep <- !(featureNames(m_set) %in% exclude_cpgs)
  
    m_set <- m_set[cpgs_keep, ]
    
    message(paste0(sum(cpgs_keep == FALSE)," non CpG probes removed: ",
                   sum(cpgs_keep == TRUE), " probes remaining..."))
    m_set

}

filter.FlaggedProbes <- function(m_set) {
    ## filter inaccurate/underperforming probes as provided by illumina
    ## load flagged EPICv2 underperforming probes provided by Illumina
    message("remove Illumina flagged underperforming probes...")
    
    underperf_probes <- read.csv("./methylome/illumina-probe-files/EPIC-8v2-0_A1-FlaggedProbes.csv",
                                 header = TRUE)

    underperf_probes <- underperf_probes$IlmnID
    perf_probes <- !(featureNames(m_set) %in% underperf_probes)

    table(perf_probes)
 
    m_set <- m_set[perf_probes, ]

    message(paste0(sum(perf_probes == FALSE)," underperforming probes removed: ",
                   sum(perf_probes == TRUE), " probes remaining..."))

    message("remove Illumina flagged inaccurate probes...")
    ## load flagged EPICv2 inaccurate probes provided by Illumina
    inaccurate_probes <- read.csv("./methylome/illumina-probe-files/EPIC-8v2-0_A1-190MappingInaccuracies.csv",
                                  header = TRUE)

    inaccurate_probes <- inaccurate_probes$IlmnID
    accu_probes <- !(featureNames(m_set) %in% inaccurate_probes)
  
    m_set <- m_set[accu_probes, ]

    message(paste0(sum(accu_probes == FALSE)," inaccurate probes removed: ",
                   sum(accu_probes == TRUE), " probes remaining..."))

## filter masked probes

    message("remove masked probes...")
    cpgs_epicv2 <- featureNames(m_set)
    cpgs_split <- strsplit2(cpgs_epicv2, "_")

    masked_probes <- read.delim("./methylome/illumina-probe-files/AppendixD_Zhou_et_al_MASKgeneral_list.txt",
                            header = FALSE)

    mask_keep <- !(cpgs_split[, 1] %in% masked_probes$V1)

    m_set <- m_set[mask_keep, ]

    message(paste0(sum(mask_keep == FALSE)," masked probes removed: ",
                   sum(mask_keep == TRUE), " probes remaining..."))


    ## filter cross reactive probes that remain from EPICv1
    cpgs_epicv2 <- featureNames(m_set)
    cpgs_split_2 <- strsplit2(cpgs_epicv2, "_")
    xr_epicv1 <- read.delim("./methylome/illumina-probe-files/AppendixE_CrossReactiveProbes_EPICv1.txt",
                            header = TRUE)

    non_xr_keep <- !(cpgs_split_2[, 1] %in% xr_epicv1)

    message(paste0(sum(non_xr_keep == FALSE)," cross reactive probes removed: ",
                   sum(non_xr_keep == TRUE), " probes remaining..."))
    m_set
}


filter.EPICv2 <- function(m_set, rg_set, Sex = FALSE, DetP = TRUE, lBead = TRUE, bead.number = 3, NonCpg = TRUE, flagged = FALSE) {
    
    if(Sex == TRUE) { m_set <- filter.Sex(m_set)
    } else{
    m_set
}

     if(DetP == TRUE) { m_set <- filter.DetP(m_set, rg_set)
    } else{
    m_set
}
          if(lBead == TRUE) { m_set <- filter.LowBead(m_set, rg_set, bead.number)
    } else{
    m_set
    }

        if(NonCpg == TRUE) { m_set <- filter.NonCpg(m_set)
    } else{
    m_set
    }

        if(flagged == TRUE) { m_set <- filter.FlaggedProbes(m_set)
    } else{
    m_set
}
    m_set    
  }
