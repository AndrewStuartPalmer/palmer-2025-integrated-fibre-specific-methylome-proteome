## Author: Andrew Palmer
## Title: EPIC probe annotations
## Date: 2025-06-14

## License: GPLv3 or later
## Copyright (c) 2025 
# See the 'LICENSE' file in the root directory for full license details.

####################################################################
## Annotate EPIC probes
####################################################################

## This script annotates EPIC probes with CpG Island and ChromHMMM states.
load_packages <- c("GEOquery", "data.table")
lapply(load_packages, library, character.only = TRUE)



####################################################################
## Download Annotation Files
####################################################################

download.file(url="https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPICv2/EPICv2.hg38.manifest.gencode.v41.tsv.gz",
destfile=("./methylome/annotations/epicv2.hg38.manifest.gencode.v41.tsv.gz"))
              
R.utils::gunzip(filename="./methylome/annotations/epicv2.hg38.manifest.gencode.v41.tsv.gz",
       destname ="./methylome/annotations/epicv2.hg38.manifest.gencode.v41.tsv",
       remove=T)

## download CpGi and ChromHMM states from the following source
## https://github.com/zhou-lab/KYCG_knowledgebase_EPICv2/tree/main
## load all files

download.file(url="https://github.com/zhou-lab/KYCG_knowledgebase_EPICv2/raw/refs/heads/main/hg38/ChromHMM.20220303.gz",
              destfile=("./methylome/annotations/ChromHMM.20220303.gz"))

GEOquery::gunzip(filename="./methylome/annotations/ChromHMM.20220303.gz",
       destname ="./methylome/annotations/ChromHMM.20220303",
       remove=T)

download.file(url="https://github.com/zhou-lab/KYCG_knowledgebase_EPICv2/raw/refs/heads/main/hg38/CGI.20220904.gz",
              destfile=("./methylome/annotations/CGI.20220904.gz"))

gunzip(filename="./methylome/annotations/CGI.20220904.gz",
       destname ="./methylome/annotations/CGI.20220904",
       remove=T)



####################################################################
## Read in annotation files
####################################################################

chrom_hmm <- read.delim("./methylome/annotations/ChromHMM.20220303", header = TRUE)
colnames(chrom_hmm) <- c("probeID", "chromstat")
cgi <- read.delim("./methylome/annotations/CGI.20220904", header = TRUE)
colnames(cgi) <- c("probeID", "cgi")


## read in gnecode manifest file
epic_v2_file <- read.delim("./methylome/annotations/epicv2.hg38.manifest.gencode.v41.tsv")
epic_v2_file <- subset(epic_v2_file, !is.na(CpG_chrm))
epic_v2_file <- subset(epic_v2_file, !is.na(CpG_beg))
epic_v2_file <- subset(epic_v2_file, !is.na(CpG_end))
epic_v2_file <- subset(epic_v2_file, !is.na(probeID))


####################################################################
## Annotate EPIC probes
####################################################################

selections <- c("probeID","CpG_chrm", "CpG_beg", "CpG_end", "genesUniq")

my_annotation <- epic_v2_file[,selections]

my_annotation <- merge(my_annotation, chrom_hmm, by = "probeID", all.x =TRUE)

cgi_unique <- cgi[!duplicated(cgi$probeID), ]
my_annotation <- merge(my_annotation, cgi_unique, by = "probeID", all.x = TRUE)

anno_2 <- data.frame(data.table::fread("./methylome/annotations/epic-v2-illumina-annotation.csv"))

anno_2_subset <- anno_2[match(my_annotation$probeID, anno_2$IlmnID),]

my_annotation$epic_v1_loci <- anno_2_subset$EPICv1_Loci
my_annotation$Methyl_450_loci <- anno_2_subset$Methyl450_Loci
my_annotation$Methyl_27_loci <- anno_2_subset$Methyl27_Loci

 write.table(my_annotation,
                file=("./methylome/annotations/annotation-epic-v2.txt"),
                quote=FALSE,
                row.names=FALSE,
                col.names=TRUE,
                sep='\t')




