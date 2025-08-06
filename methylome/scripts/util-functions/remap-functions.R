## Author: Andrew Palmer
## Title: Repurposed Remap functions from DMRcate
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
## The following lines were used and modified from [cpg.annotate]:
## Source: https://github.com/timpeters82/DMRcate-devel/blob/master/R/cpg.annotate.R
## - Lines 24-80
## - Modifications: 
##


library(AnnotationHub)
library(sesame)
library(minfi)

remap.Offtarget <- function(M){
grset <- makeGenomicRatioSetFromMatrix(mat = as.matrix(M), 
                                       array = "IlluminaHumanMethylationEPICv2",
                                       annotation = "20a1.hg38",
                                       what = "M")

anno <- getAnnotation(grset)
object <- getM(grset)

ah <- AnnotationHub()
EPICv2manifest <- ah[["AH116484"]]
EPICv2manifest <- EPICv2manifest[rownames(EPICv2manifest) %in% rownames(anno),]
EPICv2manifest <- EPICv2manifest[match(rownames(anno),rownames(EPICv2manifest)),]
                        

anno <- cbind(anno, EPICv2manifest[rownames(anno), 73:80])

    if(any(anno$CH_WGBS_evidence=="Y")){
      
        torm <- sum(anno$CH_WGBS_evidence=="Y")
        message(paste0("Remapping ", torm, " cross-hybridising probes to their more likely offtarget..."))
        anno$chr[anno$CH_WGBS_evidence=="Y"] <- gsub(":.*", "", anno$Suggested_offtarget[anno$CH_WGBS_evidence=="Y"])
        anno$pos[anno$CH_WGBS_evidence=="Y"] <- as.integer(gsub(".*:", "", anno$Suggested_offtarget[anno$CH_WGBS_evidence=="Y"]))
        anno <- anno[!(anno$CH_BLAT=="Y" & anno$CH_WGBS_evidence==""),]
        m_vals <- object[rownames(anno),]
        m_vals
    }
}


collapse.Probes <- function(new_m) {
grset <- makeGenomicRatioSetFromMatrix(mat = as.matrix(new_m), 
                                       array = "IlluminaHumanMethylationEPICv2",
                                       annotation = "20a1.hg38",
                                       what = "M")

anno <- getAnnotation(grset)
object <- getM(grset)

ah <- AnnotationHub()
EPICv2manifest <- ah[["AH116484"]]
EPICv2manifest <- EPICv2manifest[rownames(EPICv2manifest) %in% rownames(anno),]
EPICv2manifest <- EPICv2manifest[match(rownames(anno),rownames(EPICv2manifest)),]
                        

anno <- cbind(anno, EPICv2manifest[rownames(anno), 73:80])

coords <- paste(anno$chr, anno$pos, sep=":")
posreps <- table(coords)

if (any(posreps > 1)){
    message(paste("Replicate probes that map to the same CpG site found. Filtering these by caluculating replicate means"))}
    posreps <- names(posreps)[posreps > 1]
    
    message("Averaging probes that map to the same CpG site...")
    outs <- lapply(posreps, function (x){
        ids <- coords==x
        means <- colMeans(object[ids,])
        retain <- rownames(anno)[ids][1]
        dups <- rownames(anno)[ids][-1]
        list(means, retain, dups)
    })
    means <- do.call("rbind", lapply(outs, function (x) x[[1]]))
    rownames(means) <- unlist(lapply(outs, function (x) x[[2]]))
    object[rownames(means),] <- means
    dups <- unlist(lapply(outs, function (x) x[[3]]))
    anno <- anno[!rownames(anno) %in% dups,]
    object <- object[!rownames(object) %in% dups,]
    object
}

