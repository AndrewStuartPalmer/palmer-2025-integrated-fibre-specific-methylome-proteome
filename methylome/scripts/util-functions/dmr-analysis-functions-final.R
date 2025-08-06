## Author: Andrew Palmer
## Title: Methylation DMR Functions
## Date: 14-06-2025

## Copyright (c) [2025] [Andrew Palmer]
##
## Portions of this code are derived from:
## DMRcate_3.0.0
## Copyright (c) 2024, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
## The following additional terms apply under clause 7 of that license:
## EXCEPT AS EXPRESSLY STATED IN THIS AGREEMENT AND TO THE FULL EXTENT PERMITTED BY APPLICABLE LAW, THE SOFTWARE IS PROVIDED "AS-IS". CSIRO MAKES NO REPRESENTATIONS, WARRANTIES OR CONDITIONS OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO ANY REPRESENTATIONS, WARRANTIES OR CONDITIONS REGARDING THE CONTENTS OR ACCURACY OF THE SOFTWARE, OR OF TITLE, MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, THE ABSENCE OF LATENT OR OTHER DEFECTS, OR THE PRESENCE OR ABSENCE OF ERRORS, WHETHER OR NOT DISCOVERABLE.
## TO THE FULL EXTENT PERMITTED BY APPLICABLE LAW, IN NO EVENT SHALL CSIRO BE LIABLE ON ANY LEGAL THEORY (INCLUDING, WITHOUT LIMITATION, IN AN ACTION FOR BREACH OF CONTRACT, NEGLIGENCE OR OTHERWISE) FOR ANY CLAIM, LOSS, DAMAGES OR OTHER LIABILITY HOWSOEVER INCURRED. WITHOUT LIMITING THE SCOPE OF THE PREVIOUS SENTENCE THE EXCLUSION OF LIABILITY SHALL INCLUDE: LOSS OF PRODUCTION OR OPERATION TIME, LOSS, DAMAGE OR CORRUPTION OF DATA OR RECORDS; OR LOSS OF ANTICIPATED SAVINGS, OPPORTUNITY, REVENUE, PROFIT OR GOODWILL, OR OTHER ECONOMIC LOSS; OR ANY SPECIAL, INCIDENTAL, INDIRECT, CONSEQUENTIAL, PUNITIVE OR EXEMPLARY DAMAGES, ARISING OUT OF OR IN CONNECTION WITH THIS AGREEMENT, ACCESS OF THE SOFTWARE OR ANY OTHER DEALINGS WITH THE SOFTWARE, EVEN IF CSIRO HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH CLAIM, LOSS, DAMAGES OR OTHER LIABILITY.
## APPLICABLE LEGISLATION SUCH AS THE AUSTRALIAN CONSUMER LAW MAY APPLY REPRESENTATIONS, WARRANTIES, OR CONDITIONS, OR IMPOSES OBLIGATIONS OR LIABILITY ON CSIRO THAT CANNOT BE EXCLUDED, RESTRICTED OR MODIFIED TO THE FULL EXTENT SET OUT IN THE EXPRESS TERMS OF THIS CLAUSE ABOVE "CONSUMER GUARANTEES". TO THE EXTENT THAT SUCH CONSUMER GUARANTEES CONTINUE TO APPLY, THEN TO THE FULL EXTENT PERMITTED BY THE APPLICABLE LEGISLATION, THE LIABILITY OF CSIRO UNDER THE RELEVANT CONSUMER GUARANTEE IS LIMITED (WHERE PERMITTED AT CSIRO'S OPTION) TO ONE OF FOLLOWING REMEDIES OR SUBSTANTIALLY EQUIVALENT REMEDIES:
## (a) THE REPLACEMENT OF THE SOFTWARE, THE SUPPLY OF EQUIVALENT SOFTWARE, OR SUPPLYING RELEVANT SERVICES AGAIN;
## (b) THE REPAIR OF THE SOFTWARE;
## (c) THE PAYMENT OF THE COST OF REPLACING THE SOFTWARE, OF ACQUIRING EQUIVALENT SOFTWARE, HAVING THE RELEVANT SERVICES SUPPLIED AGAIN, OR HAVING THE SOFTWARE REPAIRED.
## IN THIS CLAUSE, CSIRO INCLUDES ANY THIRD PARTY AUTHOR OR OWNER OF ANY PART OF THE SOFTWARE OR MATERIAL DISTRIBUTED WITH IT. CSIRO MAY ENFORCE ANY RIGHTS ON BEHALF OF THE RELEVANT THIRD PARTY.
##
## The following lines were used and modified from [extractRanges.R]:
## Source: https://github.com/timpeters82/DMRcate-devel/blob/master/R/extractRanges.R
## - Lines 19-29
## - Modifications: 
##

####################################################################
## Load required functions and packages
####################################################################

load_packages <- c("data.table",
                   "DMRcate",
                   "ENmix",
                   "pbapply")

lapply(load_packages, library, character.only = TRUE)

####################################################################
## function to run the DMR analysisfrom DMRcate
####################################################################

run.DMRcate <- function(dmps) {
  
## load accompanying annotations file
cpgs <- rownames(dmps)

selections <- c("probeID","CpG_chrm", "CpG_beg", "CpG_end", "genesUniq")
    
anno_epic_v2_sub <-  data.frame(fread("./methylome/annotations/annotation-epic-v2.txt"))

my_annotation <- anno_epic_v2_sub[,selections]

my_annotation_subset <- my_annotation[match(rownames(dmps), my_annotation$probeID), ]

dmps$CpG_chrm <- my_annotation_subset$CpG_chrm
dmps$CpG_beg<- my_annotation_subset$CpG_beg
dmps$CpG_end <- my_annotation_subset$CpG_end

my_annotation_subset <- subset(my_annotation_subset, !is.na(CpG_chrm))
my_annotation_subset <- subset(my_annotation_subset, !is.na(CpG_end))
my_annotation_subset <- subset(my_annotation_subset, !is.na(CpG_beg))
my_annotation_subset <- subset(my_annotation_subset, !is.na(probeID))

rownames(my_annotation_subset) <- my_annotation_subset$probeID

keep <- rownames(dmps) %in% (my_annotation_subset$probeID)
dmps <-  dmps[keep, ]
cpgs <- rownames(dmps)

## set stat and diff to both logFC beta and t_beta for beta files

my_annotated_basic <- GRanges(as.character(my_annotation_subset[, "CpG_chrm"], genome = "hg38"),
                        IRanges(my_annotation_subset[, "CpG_beg"],
                                my_annotation_subset[, "CpG_end"]),
                        stat = dmps[, "t"], #t-statistic
                        diff = dmps[, "logFC_beta"],
                        ind.fdr = dmps$adj.p.val,
                        is.sig = dmps$adj.p.val < 0.001,
                        genome = "hg38") #p-value threshold

names(my_annotated_basic) <- cpgs
my_annotated_basic <- sort(my_annotated_basic)
my_annotated_basic$gene <- dmps$GeneName
my_annotated_basic$logFC_beta <- dmps$logFC_beta
my_annotated_basic$probeID <- cpgs

## convert custom annotation to object of class CpGannotated

my_annotated <- new("CpGannotated",
                    ranges = my_annotated_basic)

## run DMRcate and extract ranges

DMR_T1vsT2_m_val <- dmrcate(my_annotated,
                            C=2,
                            min.cpgs = 3)

    results_ranges <- extractRanges(DMR_T1vsT2_m_val, 
                                    genome = "hg38")

    results_ranges
}

####################################################################
## function to add ensemble id
####################################################################

add.Ensemble <- function(gr, gtf) {

## new test code for annotation of ranges
 ## overlap ensemble ids

    genesidx <- as.data.frame(findOverlaps(gr, gtf))
    message("Annotating Ensemble IDs...")   
       message("Obtaining overlaps...")
    genesover <- pbapply::pbtapply(genesidx$subjectHits, genesidx$queryHits, function(x) gtf$gene_id[x], cl = 5)
       message("Collapsing Ensemble IDs...")
    op.A <- pbapply::pbsapply(genesover, function(l) paste(l, collapse= ", "), cl = 5)
    genesover_unique <- lapply(genesover, function(l) unique(l))
       message("Adding Ensemble IDs...")
    op.A <- pbapply::pbsapply(genesover_unique, function(l) paste(l, collapse= ", "), cl = 5)

    name.A <- names(genesover_unique)

    m.A <- as.numeric(name.A)
    M <- length(gr)
    genes.ensemble <- rep(NA_character_, M)
    genes.ensemble[m.A] <- op.A
    gr$ensemble_genes <- genes.ensemble
    gr

}

####################################################################
##  function to add symbols
####################################################################

add.Symbols <- function(gr, gtf) {
    message("Annotating Gene Symbols...")   
    message("Obtaining overlaps...")
    genesidx <- as.data.frame(findOverlaps(gr, gtf))
    genesover <- pbapply::pbtapply(genesidx$subjectHits, genesidx$queryHits, function(x) gtf$gene_name[x], cl = 5)
         message("Collapsing Gene Symbols...")
op.A <- pbapply::pbsapply(genesover, function(l) paste(l, collapse= ", "), cl = 5)
    genesover_unique <- lapply(genesover, function(l) unique(l))

    message("Adding Gene Symbols...")
op.A <- pbapply::pbsapply(genesover_unique, function(l) paste(l, collapse= ", "), cl = 5)

name.A <- names(genesover_unique)

m.A <- as.numeric(name.A)
M <- length(gr)
overlapping.genes <- rep(NA_character_, M)
overlapping.genes[m.A] <- op.A
gr$symbol <- overlapping.genes
gr
}

####################################################################
##  function to add gene type
####################################################################

add.GeneType <- function(gr, gtf) {
    message("Annotating Gene Types...")
    message("Obtaining overlaps...")

## gene type
genesidx <- as.data.frame(findOverlaps(gr, gtf))
    genesover <- pbapply::pbtapply(genesidx$subjectHits, genesidx$queryHits, function(x) gtf$gene_type[x], cl = 5)
      message("Collapsing Gene Types...")
    op.A <- pbapply::pbsapply(genesover, function(l) paste(l, collapse= ", "), cl = 5)
    genesover_unique <- lapply(genesover, function(l) unique(l))

    message("Adding GeneType...")
op.A <- pbapply::pbsapply(genesover_unique, function(l) paste(l, collapse= ", "), cl = 5)

name.A <- names(genesover_unique)

m.A <- as.numeric(name.A)
M <- length(gr)
overlapping.genes <- rep(NA_character_, M)
overlapping.genes[m.A] <- op.A
gr$biotype <- overlapping.genes
gr
}

####################################################################
##  function to add distance to nearest gene
####################################################################

add.NearestGene <- function(gr, gtf) {
message("Calculating distance to nearest gene...")
    nearest_gene <- as.data.frame(distanceToNearest(gr, gtf))
    gr$nearest_gene <- nearest_gene$distance
    gr

}

####################################################################
##  function to add probes in DMR
####################################################################

add.Probes <- function(DMRs_gr, CpGs_gr) {
  message("Obtaining probes in DMRs...")

  overlaps <- findOverlaps(DMRs_gr, CpGs_gr)

  DMRs_gr$overlapping_probes <- character(length(DMRs_gr)) 

  for (i in seq_along(DMRs_gr)) {
    cpg_indices <- subjectHits(overlaps)[queryHits(overlaps) == i]
    if (length(cpg_indices) > 0) {
      overlapping_probes <- CpGs_gr$probeID[cpg_indices] # Use probeID instead of pos
      DMRs_gr$overlapping_probes[i] <- paste(overlapping_probes, collapse = ";")
    }
  }
  return(DMRs_gr)
}

####################################################################
## add annotations function that strings together the additional
## annotations
####################################################################

add.Annotations <- function(gr, gtf, ensemble = TRUE, symbols = TRUE, geneType = TRUE, nearestGene = FALSE) {


    if(ensemble == TRUE) { gr <- add.Ensemble(gr, gtf)
     } else{
        gr
    }
    if(symbols == TRUE) { gr <- add.Symbols(gr, gtf)
        } else{
        gr
    }
    if(geneType == TRUE) { gr <- add.GeneType(gr, gtf)
     } else{
        gr
    }
    if(nearestGene == TRUE) { gr <- add.NearestGene(gr, gtf)
          } else{
        gr
    }

    gr
    }
