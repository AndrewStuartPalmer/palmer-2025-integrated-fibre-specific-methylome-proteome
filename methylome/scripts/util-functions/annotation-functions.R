## Author: Andrew Palmer
## Title: Methylation Results Annotation functions
## Date: 14-06-2025

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
