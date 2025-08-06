## Author: Andrew Palmer
## Title: MissMethyl functions altered for use with EPICv2
## Date: 14-06-2025

## Copyright (c) [2025] [Andrew Palmer]
##
## Portions of this code are derived from:
## missMethyl package
## Author: Belinda Phipson and Jovana Maksimovic
## Maintainer: Belinda Phipson <phipson.b@wehi.edu.au>, Jovana Maksimovic <jovana.maksimovic@petermac.org>, Andrew Lonsdale <andrew.lonsdale@petermac.org>
## Source: [https://www.bioconductor.org/packages/release/bioc/html/missMethyl.html]
##         [https://github.com/Oshlack/missMethyl]    
##
## This code is distributed under the terms of the GNU General Public License (GPL) version 2.
## See the LICENSE file for the full text of the license.
##
## The following lines were used and modified from missMethyl to work with the EPICv2 array:

####################################################################
## This code is a mix of modified miss methyl code for hypergeometric
## overrepresentation testing.
####################################################################

## Load required libraries

library(data.table)
library(org.Hs.eg.db)
library(limma)
library(GO.db)

####################################################################
## flattens EPICv2 array ##
####################################################################

## code obtained and modified from https://github.com/Oshlack/missMethyl/blob/master/R/gometh.R
## 
## - Modifications: [Describe the changes, e.g., "Modified probe annotation to support EPICv2", "Adjusted normalization function for EPICv2 probes", etc.]
   

flattenEPICv2 <- function(array.type=c("EPIC2"),anno=NULL)
    ## flatten EPICv2 array annotation
    ## Andrew Palmer
    ## 10 March 2024
    ## Updated 20 MArch 2024
    ## Modified version of Belinda Phipson's .flattenAnn code

    ## outputs flattened annotation with the following
    ## cpg = probe id
    ## symbol = gene symbol annotated to CpG probe
    ## alias = gene alias annotated to CpG probe
    ## entrezid = gene entrezid annotated to CpG probe
    
{
  if(array.type=="EPIC2"){
    anno <- read.delim("./methylome/annotations/epicv2.hg38.manifest.gencode.v41.tsv")

  rownames(anno) <- anno$probeID
 # get rid of the non-CpG sites
  ann.keep<-anno[grepl("^cg",anno$probeID),]
  
  # get rid of CpGs that are not annotated
  missing<-ann.keep$geneNames==""
  ann.keep<-ann.keep[!missing,]
  
  # get individual gene names for each CpG
  geneslist<-ann.keep$geneNames
  names(geneslist)<-rownames(ann.keep)
    
  flat<-data.frame(symbol=unlist(geneslist))

   splits <- limma::strsplit2(flat$symbol, split=";")
   splits <- as.data.frame(splits)
   splits <- splits[,1:2]
   splits$V1 <- ifelse(grepl("ENSG000", splits$V1), splits$V2, splits$V1)
   gene_replace <- splits$V1
   flat$symbol<-as.character(gene_replace)

    
  flat$cpg<- rownames(flat)
  
  #flat$cpg <- rownames(flat)
  flat$alias <- suppressWarnings(limma::alias2SymbolTable(flat$symbol))
  
  #eg <- toTable(org.Hs.egSYMBOL2EG)
  eg <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, 
                                keys=AnnotationDbi::keys(org.Hs.eg.db), 
                                columns=c("ENTREZID","SYMBOL"), 
                                keytype="ENTREZID"))
  colnames(eg) <- c("gene_id","symbol")
  
  m <- match(flat$alias,eg$symbol)
  flat$entrezid <- eg$gene_id[m]
  flat <- flat[!is.na(flat$entrezid),]
  
  # keep unique cpg by gene name annotation
  id<-paste(flat$cpg,flat$entrezid,sep=".")
  d <- duplicated(id)
  flat.u <- flat[!d,]
    flat.u
    
      }
  }

####################################################################
## Code for getting mapped enterez ids and probe to gene ratios ##
####################################################################

## code obtinaed and modified from https://github.com/Oshlack/missMethyl/blob/master/R/getMappedEntrezIDs.R
## - Lines [Line Numbers]
## - Modifications: [Describe the changes, e.g., "Modified probe annotation to support EPICv2", "Adjusted normalization function for EPICv2 probes", etc.]
   

getMappedEntrezIDs <- function(sig.cpg, all.cpg=NULL, 
                               array.type=c("EPIC2"), anno=NULL)
    
    ## map CpG list to entrez ids
    ## Andrew Palmer
    ## 10 March 2024
    ## Updated 20 March 2024
    ## Modified version of Belinda Phipson's .getMappedEntrezIDs code

    ## outputs a list of the following
    ## sig.eg = mapped entrezids for sig cpgs
    ## universe = mapped enterezID for all probes on array
    ## freq = number of probes associated with each gene
    ## equiv = take into account multi gene bias
    ## de = 1s and 0s indicating which genes are DMPs
{
    # check input
    sig.cpg <- as.character(sig.cpg)
    sig.cpg <- sig.cpg[!is.na(sig.cpg)]
      
    # Get annotaton in appropriate format
    flat.u <- flattenEPICv2()
  
# remove CpGs from annotation that are not in all.cpg
    m_all <- match(flat.u$cpg, all.cpg)
    flat.u = flat.u[!is.na(m_all),]
    
    # map CpG sites to entrez gene id's
    sig.cpg <- unique(sig.cpg)
    
    m1 <- match(flat.u$cpg,sig.cpg)
    
  eg.sig <- flat.u$entrezid[!is.na(m1)]
       
    eg.sig <- unique(eg.sig)
    if(length(eg.sig)==0) {
        stop("There are no genes annotated to the significant CpGs")
    }
    
    m2 <- match(flat.u$cpg,all.cpg)
    eg.all <- flat.u$entrezid[!is.na(m2)]
    
    freq_genes <- table(eg.all)
    eg.universe <- names(freq_genes)
    
    test.de <- as.integer(eg.universe %in% eg.sig)
    
    sorted.eg.sig <- eg.universe[test.de==1]
    
    multimap <- data.frame(table(flat.u$cpg))
    multimap$Var1 <- as.character(multimap$Var1)
    m3 <- match(flat.u$cpg, multimap$Var1)
    flat.u$multimap <- multimap$Freq[m3]
    
    flat.u$inv.multimap <- 1/flat.u$multimap
    
    equivN <- tapply(flat.u$inv.multimap,flat.u$entrezid,sum)
    mm <- match(eg.universe,names(equivN))
    equivN <- equivN[mm]
    
   sig.flat <- flat.u[!is.na(m1),]
   
    
    fract <- data.frame(weight=pmin(tapply(1/sig.flat$multimap,
                                           sig.flat$entrezid,sum),
                                    1))
    
    m4 <- match(sorted.eg.sig,rownames(fract))
    fract.counts <- fract$weight[m4]
    
    out <- list(sig.eg = sorted.eg.sig, universe = eg.universe, 
                freq = freq_genes, equiv =  equivN, de = test.de, 
                fract.counts = data.frame(sigid=sorted.eg.sig,frac=fract.counts))
    out
}


####################################################################
## estimate weighted probablities ##
####################################################################

## code obtained from https://github.com/Oshlack/missMethyl/blob/master/R/gometh.R
## - Lines [Line Numbers]
## - Modifications: [Describe the changes, e.g., "Modified probe annotation to support EPICv2", "Adjusted normalization function for EPICv2 probes", etc.]
   

estimatePWF <- function(D,bias)
    ## An alternative to goseq function nullp, which is transformation invariant
    ## map CpG list to entrez ids
    ## Andrew Palmer
    ## 10 March 2024
    ## Updated 20 March 2024
    ## Modified version of Belinda Phipson's .estimatePWF code

    ## returns a vector of weighted probablities for each gene
    ## using probablity weighting function in limma
     
{
  prior.prob <- bias
  o <- order(bias)
  prior.prob[o] <- limma::tricubeMovingAverage(D[o],span=0.5)
  prior.prob
}

####################################################################
## Hypergeometric testing with prior probabalities included of
## GO terms
####################################################################
## code obtained and modfied from https://github.com/Oshlack/missMethyl/blob/master/R/gsameth.R
## - Lines [Line Numbers]
## - Modifications: [Describe the changes, e.g., "Modified probe annotation to support EPICv2", "Adjusted normalization function for EPICv2 probes", etc.]
   

hypergeometricTest <- function(pwf = pwf,
                               collection = collection,
                               sorted.eg.sig =  sorted.eg.sig,
                               eg.universe = eg.universe,
                               freq_genes = freq_genes,
                               test.de = test.de,
                               frac = frac,
                               equiv = equiv)
    
    ## hypergeomteric test including prior probablities based on number of probes per gene
    ## Andrew Palmer
    ## 10 March 2024
    ## Updated 20 MArch 2024
    ## Modified version of Belinda Phipson's gometh and gsameth code
    
    ## requires first the running of getMappedEntrezIDs function 
    ## returns a data frame with the following
    ## GOID
    ## N = number of genes in GO set
    ## DE = number of sig,cpgs in GO set
    ## P.DE = P.value of hypergeometric test
    ## FDR = benjamani hochberg corrected p value


{

    results <- matrix(NA,ncol=4,nrow=length(collection))
    
  colnames(results) <- c("N","DE","P.DE","FDR")
  rownames(results) <- names(collection)
results[,"N"] <- unlist(lapply(collection,length))

## use fractional counting to account for cpgs that map to multiple genes
 results[,"DE"] <- unlist(lapply(collection, function(x) 
     sum((sorted.eg.sig %in% x) * frac$frac)))

Nuniverse <- length(eg.universe)
m <- length(sorted.eg.sig)

## Hypergeometric test with prior probabilities

for(i in 1:length(collection)){
      InSet <- eg.universe %in% collection[[i]]
      pw.red <- sum(pwf[InSet])/results[i,"N"]
      pw.white <- sum(pwf[!InSet])/(Nuniverse-results[i,"N"])
      odds <- pw.red/pw.white
      results[i,"P.DE"] <- BiasedUrn::pWNCHypergeo(results[i,"DE"],
                                                   results[i,"N"],
                                                   Nuniverse-results[i,"N"],
                                                   m,odds,lower.tail=FALSE) + 
          BiasedUrn::dWNCHypergeo(results[i,"DE"],
                                  results[i,"N"],
                                  Nuniverse-results[i,"N"],
                                  m,odds)
      
if(results[i,"P.DE"]==0){
        message("Pvalue of exactly zero detected. Performing hypergeometric 
                test for gene set ", rownames(results)[i])
        results[i,"P.DE"] <- stats::phyper(q=results[i,"DE"]-0.5,m=m,
                                           n=Nuniverse-m,k=results[i,"N"],
                                           lower.tail=FALSE)
}
}      
  results[,"FDR"] <- stats::p.adjust(results[,"P.DE"],method="BH")
  results[,"DE"] <- floor(results[,"DE"])
  data.frame(results)

    }

####################################################################
## obtain GO terms ##
####################################################################

## code obtained from https://github.com/Oshlack/missMethyl/blob/master/R/gometh.R
## - Lines [Line Numbers]
## - Modifications: [Describe the changes, e.g., "Modified probe annotation to support EPICv2", "Adjusted normalization function for EPICv2 probes", etc.]
   

getGO <- function()

    ## obtain GO terms as a list
    ## Andrew Palmer
    ## 10 March 2024
    ## Updated 20 MArch 2024
    ## Modified version of Belinda Phipson's .getMappedEntrezIDs code

    ## outputs the following
    ## list with GO ids, terms and ontologies

{
  if(!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("org.Hs.eg.db package required but not installed.")
  egGO2ALLEGS <- utils::getFromNamespace("org.Hs.egGO2ALLEGS", "org.Hs.eg.db")
  GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[,c("gene_id", "go_id", "Ontology")]
  d <- !duplicated(GeneID.PathID[, c("gene_id", "go_id")])
  GeneID.PathID <- GeneID.PathID[d, ]
  GOID.TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, 
                                                      keys=unique(GeneID.PathID$go_id), 
                                                      columns=c("GOID","ONTOLOGY","TERM"), 
                                                      keytype="GOID"))
  go <- tapply(GeneID.PathID$gene_id, GeneID.PathID$go_id, list)
    
  list(idList=go, idTable=GOID.TERM)
}


getKEGG <- function(){
  GeneID.PathID <- limma::getGeneKEGGLinks(species.KEGG = "hsa", convert = TRUE)
  GeneID.PathID$PathwayID <- gsub("path:", "", GeneID.PathID$PathwayID)
  isna <- rowSums(is.na(GeneID.PathID[, 1:2])) > 0.5
  GeneID.PathID <- GeneID.PathID[!isna, ]
  ID.ID <- paste(GeneID.PathID[, 1], GeneID.PathID[, 2], sep = ".")
  d <- !duplicated(ID.ID)
  GeneID.PathID <- GeneID.PathID[d, ]
  PathID.PathName <- limma::getKEGGPathwayNames(species.KEGG = "hsa", 
                                         remove.qualifier = TRUE)
  #PathID.PathName$PathwayID <- paste0("path:", PathID.PathName$PathwayID)
  GeneID.PathID <- merge(GeneID.PathID, PathID.PathName, by="PathwayID")
  kegg <- tapply(GeneID.PathID$GeneID, GeneID.PathID$PathwayID, list)
  
  list(idList = kegg, idTable = PathID.PathName)
}  
