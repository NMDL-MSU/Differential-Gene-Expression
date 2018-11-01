---
title: Gene Set Enrichment KER Glycogen Diets Over Time
author: Deborah Velez-Irizarry
date: Tue Oct 30 15:35:58 EDT 2018
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---
### Description:  
Performe gene set enrichment analysis for KER Glycogen project. Look at enrichment of  
genes differentially expressed for each diet over time with depletion timepoint as  
reference.   
  
***  
**Code:**  
Parent Directory:  

>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq  
  
Directory/File:  
 
&nbsp;&nbsp;&nbsp;&nbsp;/Enrichment/Enrichment_KER_Glycogen_Depl/Enrichment_KER_Glycogen_Depl.R  
 
**Input files:**  
Directory/File:  
  
>&nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_KER_Glycogen/results_diet_over_time_depl_ref.Rdata  
  
**Output files:**  
  
Directory:  
  
>&nbsp;&nbsp;&nbsp;&nbsp;/Enrichment/Enrichment_KER_Glycogen_Depl/  

Files:

>&nbsp;&nbsp;&nbsp;&nbsp;/enrichment_function.Rdata
>&nbsp;&nbsp;&nbsp;&nbsp;/Enrichment_KER_Glycogen_Depl.Rdata  
>&nbsp;&nbsp;&nbsp;&nbsp;/Fat.Repletion.24h/*.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;/Sweet.Repletion.24h/*.txt
>&nbsp;&nbsp;&nbsp;&nbsp;/Fat.Repletion.72h/*.txt
>&nbsp;&nbsp;&nbsp;&nbsp;/Sweet.Repletion.72h/*.txt
>&nbsp;&nbsp;&nbsp;&nbsp;/Merged/*.txt
>&nbsp;&nbsp;&nbsp;&nbsp;/upregulated/*.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;/downregulated/*.txt  
 
***  
### Code  
Load required libraries


```r
library(edgeR)
```

```
## Loading required package: limma
```

```r
library(limma)
library(DOSE)
```

```
## 
```

```
## DOSE v3.0.10  For help: https://guangchuangyu.github.io/DOSE
## 
## If you use DOSE in published research, please cite:
## Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015, 31(4):608-609
```

```r
library(GO.db)
```

```
## Loading required package: AnnotationDbi
```

```
## Loading required package: stats4
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following object is masked from 'package:limma':
## 
##     plotMA
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, cbind, colnames,
##     do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, lengths, Map, mapply,
##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
##     Position, rank, rbind, Reduce, rownames, sapply, setdiff,
##     sort, table, tapply, union, unique, unsplit, which, which.max,
##     which.min
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Loading required package: IRanges
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following objects are masked from 'package:base':
## 
##     colMeans, colSums, expand.grid, rowMeans, rowSums
```

```r
library(GSEABase)
```

```
## Loading required package: annotate
```

```
## Loading required package: XML
```

```
## Loading required package: graph
```

```
## 
## Attaching package: 'graph'
```

```
## The following object is masked from 'package:XML':
## 
##     addNode
```

```r
library(clusterProfiler)
```

```
## clusterProfiler v3.2.14  For help: https://guangchuangyu.github.io/clusterProfiler
## 
## If you use clusterProfiler in published research, please cite:
## Guangchuang Yu., Li-Gen Wang, Yanyan Han, Qing-Yu He. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.
```

```r
library(org.Hs.eg.db)
```

```
## 
```

```r
library(GenomicFeatures)
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: GenomicRanges
```

Clear Environment


```r
rm(list=ls())
```

Load DE results


```r
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq/DE/DE_KER_Glycogen"
load(paste(dir, "results_diet_over_time_depl_ref.Rdata", sep="/"))
```

### Animal Groups


```r
group <- paste(anim$Diet, anim$TimePoint, sep=":")
table(group)
```

```
## group
##         Fat:Depl          Fat:Pre       Fat:Rep24h       Fat:Rep72h 
##                5                5                5                5 
##   SweetFeed:Depl    SweetFeed:Pre SweetFeed:Rep24h SweetFeed:Rep72h 
##                5                5                5                5
```

### Differential Expression Analysis Results
Reduce results list to contain only differentially expressed genes


```r
Rrst <- lapply(Rst.Diet.Depl, function(x) x[x$adj.P.Val < 0.05,])
lapply(Rrst, nrow)
```

```
## $Fat.Repletion.24h
## [1] 551
## 
## $Sweet.Repletion.24h
## [1] 32
## 
## $Fat.Repletion.72h
## [1] 4191
## 
## $Sweet.Repletion.72h
## [1] 2135
## 
## $Fat.over.time
## [1] 7819
## 
## $Sweet.over.time
## [1] 4692
```

Write DE results to file


```r
system("mkdir DE.Results")
lapply(names(Rrst), function(x) 
    write.table(Rrst[[x]], file=paste(getwd(), "/DE.Results/", gsub("[.]", "_", x), "_depl_ref.txt", sep=""),
    col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t"))
```

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
## 
## [[3]]
## NULL
## 
## [[4]]
## NULL
## 
## [[5]]
## NULL
## 
## [[6]]
## NULL
```

Merge list of genes


```r
sigG <- do.call(rbind, lapply(Rrst, function(x) data.frame(Xloc=rownames(x), x[,c(1:5,7)])))
rownames(sigG) <- NULL
nrow(sigG)
```

```
## [1] 19420
```

Seperate genes without a name


```r
na.sigG <- sigG[is.na(sigG$genes),]
nrow(na.sigG)
```

```
## [1] 616
```

```r
sigG <- sigG[!is.na(sigG$genes),]
nrow(sigG)
```

```
## [1] 18804
```

### Background for Enrichment Analysis
Background Gene List


```r
annot <- dge$genes
```

Seperate genes without a name


```r
na.Bkg <- annot[is.na(annot$genes),]
nrow(na.Bkg)
```

```
## [1] 767
```

```r
Bkg <- annot[!is.na(annot$genes),]
nrow(Bkg)
```

```
## [1] 13366
```

Obtain human EntrezIDs for gene enrichment analysis


```r
# DE gene list
sigG.entrez <- bitr(as.character(unique(sigG$genes)), fromType="SYMBOL", 
    toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```
## Warning in bitr(as.character(unique(sigG$genes)), fromType = "SYMBOL",
## toType = c("ENTREZID"), : 13.48% of input gene IDs are fail to map...
```

```r
nrow(sigG.entrez)
```

```
## [1] 7738
```

```r
# Background gene list
bkg.entrez <- bitr(as.character(unique(Bkg$genes)), fromType="SYMBOL", 
    toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```
## Warning in bitr(as.character(unique(Bkg$genes)), fromType = "SYMBOL",
## toType = c("ENTREZID"), : 15.75% of input gene IDs are fail to map...
```

```r
nrow(bkg.entrez)
```

```
## [1] 10990
```

### Functions  
**Function to performe enrichment analysis from a gene list**  
> **INPUT**: use output from bitr   
>&nbsp;&nbsp;&nbsp;&nbsp;`lst` = data frame with gene symbols and gene IDs (colnames must be "SYMBOL" and "ENTREZID")  
>&nbsp;&nbsp;&nbsp;&nbsp;`bkg` = background gene list for enrichment analysis with gene symbol and gene IDs   
>&nbsp;&nbsp;&nbsp;&nbsp;;&nbsp(colnames must be "SYMBOL" and "ENTREZID")  


```r
enrich <- function(lst, bkg){
    lstE <- lst$ENTREZID
    bkgE <- bkg$ENTREZID
    # GO: biological processes (BP)
    go.bp <- enrichGO(gene = lstE, universe = bkgE, 
        OrgDb=org.Hs.eg.db, ont = "BP", pAdjustMethod = "fdr", qvalueCutoff = 0.05,
        readable = TRUE)
    # GO: cellular component (CC)
    go.cc <- enrichGO(gene = lstE, universe = bkgE, 
        OrgDb=org.Hs.eg.db, ont = "CC", pAdjustMethod = "fdr", qvalueCutoff = 0.05,
        readable = TRUE)
    # GO: molecular function (MF)
    go.mf <- enrichGO(gene = lstE, universe = bkgE, 
        OrgDb=org.Hs.eg.db, ont = "MF", pAdjustMethod = "fdr", qvalueCutoff = 0.05,
        readable = TRUE)
    # Kegg Enrichment
    kegg <- enrichKEGG(gene = lstE, universe = bkgE, 
        organism="hsa", pAdjustMethod = "fdr", qvalueCutoff = 0.05)
    # Merge and return results
    rst <- list(Genes=lst, GO.BP=go.bp, GO.CC=go.cc, GO.MF=go.mf, Kegg=kegg)
    return(rst)
}
```

Save function to file


```r
save(enrich, file=paste(getwd(), "enrichment_function.Rdata", sep="/"))
```

### Global Enrichment Analysis  
Global gene enrichment analysis: Merge all DE genes


```r
enrich.rst.fit3 <- enrich(lst=sigG.entrez, bkg=bkg.entrez)
unlist(lapply(enrich.rst.fit3, nrow))
```

```
## Genes GO.BP GO.CC GO.MF  Kegg 
##  7738   161    38     5     1
```

Summary of significant GO terms


```r
lapply(enrich.rst.fit3, function(x) head(x[,2]))
```

```
## $Genes
## [1] "8013"  "1303"  "57528" "7057"  "5166"  "84879"
## 
## $GO.BP
## [1] "response to cytokine"                       
## [2] "innate immune response"                     
## [3] "regulation of immune response"              
## [4] "response to biotic stimulus"                
## [5] "antigen receptor-mediated signaling pathway"
## [6] "cellular response to cytokine stimulus"     
## 
## $GO.CC
## [1] "melanosome"                   "pigment granule"             
## [3] "mitochondrial inner membrane" "cell surface"                
## [5] "organelle inner membrane"     "endoplasmic reticulum lumen" 
## 
## $GO.MF
## [1] "transmembrane transporter activity"                   
## [2] "substrate-specific transmembrane transporter activity"
## [3] "ion transmembrane transporter activity"               
## [4] "active transmembrane transporter activity"            
## [5] "calcium ion binding"                                  
## 
## $Kegg
## [1] "Proteasome"
```

Save enrichment analysis results to file


```r
# Enrichemnt analysis results
data_merged <- lapply(enrich.rst.fit3, function(x) data.frame(x))

# Add gene names to Kegg results (currently has the geneIDs)
kegg.gene <- strsplit(data_merged$Kegg$geneID, "/")[[1]]
data_merged$Kegg$geneID <- paste(data_merged$Gene$SYMBOL[data_merged$Gene$ENTREZID %in% kegg.gene], collapse="/")

# Save results to file
system("mkdir Merged")
lapply(names(data_merged), function(x) 
    write.table(data_merged[[x]], file=paste(getwd(), "/Merged/", x, "_merged.txt", sep=""),
    col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t"))
```

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
## 
## [[3]]
## NULL
## 
## [[4]]
## NULL
## 
## [[5]]
## NULL
```

### Enrichment Analysis Per Contrast  
Obtain human EntrezIDs per contrast:


```r
enrichAll.rst.list <- lapply(Rrst, function(x) 
    do.call(rbind, lapply(as.character(x$genes), 
    function(y) sigG.entrez[sigG.entrez$SYMBOL %in% y,])))
names(enrichAll.rst.list) <- names(Rrst)
unlist(lapply(enrichAll.rst.list, nrow))
```

```
##   Fat.Repletion.24h Sweet.Repletion.24h   Fat.Repletion.72h 
##                 448                  27                3580 
## Sweet.Repletion.72h       Fat.over.time     Sweet.over.time 
##                1738                6618                3917
```

Gene enrichment analysis per contrast:


```r
enrichCont.rst.fit3 <- lapply(enrichAll.rst.list, 
    function(x) enrich(lst=x, bkg=bkg.entrez))
lapply(enrichCont.rst.fit3, function(x)
    unlist(lapply(x, function(y) nrow(y))))
```

```
## $Fat.Repletion.24h
## Genes GO.BP GO.CC GO.MF  Kegg 
##   448     0     0     0     0 
## 
## $Sweet.Repletion.24h
## Genes GO.BP GO.CC GO.MF  Kegg 
##    27     6     0     0     0 
## 
## $Fat.Repletion.72h
## Genes GO.BP GO.CC GO.MF  Kegg 
##  3580   177     6     0     0 
## 
## $Sweet.Repletion.72h
## Genes GO.BP GO.CC GO.MF  Kegg 
##  1738   102    40     5     3 
## 
## $Fat.over.time
## Genes GO.BP GO.CC GO.MF  Kegg 
##  6618    51    19     4     3 
## 
## $Sweet.over.time
## Genes GO.BP GO.CC GO.MF  Kegg 
##  3917   147    28    14    10
```

Gene Ontology: Biological Processes (first 10)


```r
lapply(enrichCont.rst.fit3, function(x) 
    head(x[["GO.BP"]][,2]))
```

```
## $Fat.Repletion.24h
## character(0)
## 
## $Sweet.Repletion.24h
## [1] "glycolytic process"                               
## [2] "ATP generation from ADP"                          
## [3] "osteoblast development"                           
## [4] "negative regulation of glycolytic process"        
## [5] "negative regulation of cofactor metabolic process"
## [6] "negative regulation of coenzyme metabolic process"
## 
## $Fat.Repletion.72h
## [1] "hemopoiesis"                                
## [2] "cell activation"                            
## [3] "leukocyte migration"                        
## [4] "hematopoietic or lymphoid organ development"
## [5] "leukocyte differentiation"                  
## [6] "response to external biotic stimulus"       
## 
## $Sweet.Repletion.72h
## [1] "regulation of cellular amino acid metabolic process"                                            
## [2] "regulation of cellular ketone metabolic process"                                                
## [3] "antigen processing and presentation of exogenous peptide antigen via MHC class I, TAP-dependent"
## [4] "antigen processing and presentation of peptide antigen via MHC class I"                         
## [5] "antigen processing and presentation of exogenous peptide antigen via MHC class I"               
## [6] "cellular ketone metabolic process"                                                              
## 
## $Fat.over.time
## [1] "extracellular matrix organization"      
## [2] "extracellular structure organization"   
## [3] "innate immune response"                 
## [4] "collagen metabolic process"             
## [5] "response to wounding"                   
## [6] "water-soluble vitamin metabolic process"
## 
## $Sweet.over.time
## [1] "generation of precursor metabolites and energy"     
## [2] "cellular respiration"                               
## [3] "energy derivation by oxidation of organic compounds"
## [4] "glycosyl compound metabolic process"                
## [5] "ribonucleoside triphosphate metabolic process"      
## [6] "nucleoside metabolic process"
```

Gene Ontology: Cellular Processes (first 10)


```r
lapply(enrichCont.rst.fit3, function(x) 
    head(x[["GO.CC"]][,2]))
```

```
## $Fat.Repletion.24h
## character(0)
## 
## $Sweet.Repletion.24h
## character(0)
## 
## $Fat.Repletion.72h
## [1] "vacuolar lumen"                     
## [2] "lysosomal lumen"                    
## [3] "lytic vacuole"                      
## [4] "lysosome"                           
## [5] "side of membrane"                   
## [6] "cytoplasmic side of plasma membrane"
## 
## $Sweet.Repletion.72h
## [1] "organelle inner membrane"     "mitochondrial inner membrane"
## [3] "mitochondrial outer membrane" "mitochondrial membrane part" 
## [5] "peptidase complex"            "proteasome complex"          
## 
## $Fat.over.time
## [1] "endoplasmic reticulum lumen"                  
## [2] "sarcoplasmic reticulum"                       
## [3] "collagen trimer"                              
## [4] "cell surface"                                 
## [5] "sarcoplasmic reticulum membrane"              
## [6] "intrinsic component of mitochondrial membrane"
## 
## $Sweet.over.time
## [1] "mitochondrial inner membrane"                
## [2] "organelle inner membrane"                    
## [3] "mitochondrial matrix"                        
## [4] "mitochondrial membrane part"                 
## [5] "mitochondrial protein complex"               
## [6] "inner mitochondrial membrane protein complex"
```

Gene Ontology: Molecular Processes (first 10)


```r
lapply(enrichCont.rst.fit3, function(x) 
    head(x[["GO.MP"]][,2]))
```

```
## $Fat.Repletion.24h
## NULL
## 
## $Sweet.Repletion.24h
## NULL
## 
## $Fat.Repletion.72h
## NULL
## 
## $Sweet.Repletion.72h
## NULL
## 
## $Fat.over.time
## NULL
## 
## $Sweet.over.time
## NULL
```

Kegg Pathways:


```r
lapply(enrichCont.rst.fit3, function(x) 
    head(x[["Kegg"]][,2]))
```

```
## $Fat.Repletion.24h
## character(0)
## 
## $Sweet.Repletion.24h
## character(0)
## 
## $Fat.Repletion.72h
## character(0)
## 
## $Sweet.Repletion.72h
## [1] "Proteasome"                   "Ribosome"                    
## [3] "Epstein-Barr virus infection"
## 
## $Fat.over.time
## [1] "ECM-receptor interaction"                   
## [2] "Protein processing in endoplasmic reticulum"
## [3] "Hypertrophic cardiomyopathy (HCM)"          
## 
## $Sweet.over.time
## [1] "Oxidative phosphorylation"                
## [2] "Huntington disease"                       
## [3] "Parkinson disease"                        
## [4] "Non-alcoholic fatty liver disease (NAFLD)"
## [5] "Alzheimer disease"                        
## [6] "Thermogenesis"
```

Save results


```r
data_contrast <- lapply(names(enrichCont.rst.fit3), function(x) 
    lapply(names(enrichCont.rst.fit3[[x]]), function(y) 
        data.frame(enrichCont.rst.fit3[[x]][[y]])))
names(data_contrast) <- names(enrichCont.rst.fit3)
for(i in names(data_contrast)){
    # Add names to enrichment analysis per contrast
    names(data_contrast[[i]]) <- names(data_merged)
    kegg.gene <- data_contrast[[i]]$Kegg$geneID
    # Add gene names to kegg geneIDs
    if(length(kegg.gene) > 0){
        for(j in 1:length(kegg.gene)){
            idx <- strsplit(kegg.gene[j], "/")[[1]]
            data_contrast[[i]]$Kegg$geneID[j] <- paste(data_contrast[[i]]$Gene$SYMBOL
                [data_contrast[[i]]$Gene$ENTREZID %in% idx], collapse="/")
        }
    }
}

# Save results to file
for(i in names(data_contrast)){
    system(paste("mkdir", i, sep=" "))
    for (j in names(data_contrast[[i]])){
        write.table(data_contrast[[i]][[j]], file=paste(getwd(), 
            "/", i, "/", gsub("[.]", "_", i), "_", j, ".txt", sep=""),
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
    }
}
```

### Enrichment of Up-regulated Genes  
Gene set enrichment for upregulated genes


```r
up <- lapply(Rrst, function(x) x[x$logFC > 0,])
up <- lapply(up, function(x) 
    do.call(rbind, lapply(as.character(x$genes), 
    function(y) sigG.entrez[sigG.entrez$SYMBOL %in% y,])))
up <- up[!unlist(lapply(up, is.null))]
unlist(lapply(up, nrow))
```

```
##   Fat.Repletion.24h Sweet.Repletion.24h   Fat.Repletion.72h 
##                 156                   3                2044 
## Sweet.Repletion.72h 
##                 827
```

Gene enrichment analysis per contrast for upregulated genes:


```r
enrichCont.up <- lapply(up, 
    function(x) enrich(lst=x, bkg=bkg.entrez))
lapply(enrichCont.up, function(x)
    unlist(lapply(x, function(y) nrow(y))))
```

```
## $Fat.Repletion.24h
## Genes GO.BP GO.CC GO.MF  Kegg 
##   156     0     0     1     0 
## 
## $Sweet.Repletion.24h
## Genes GO.BP GO.CC GO.MF  Kegg 
##     3     0     1     0    22 
## 
## $Fat.Repletion.72h
## Genes GO.BP GO.CC GO.MF  Kegg 
##  2044   777    67    50    30 
## 
## $Sweet.Repletion.72h
## Genes GO.BP GO.CC GO.MF  Kegg 
##   827   151    25     7     4
```

Save results


```r
data_upreg <- lapply(names(enrichCont.up), function(x) 
    lapply(names(enrichCont.up[[x]]), function(y) 
        data.frame(enrichCont.up[[x]][[y]])))
names(data_upreg) <- names(enrichCont.up)
for(i in names(data_upreg)){
    # Add names to enrichment analysis per contrast
    names(data_upreg[[i]]) <- names(data_merged)
    kegg.gene <- data_upreg[[i]]$Kegg$geneID
    # Add gene names to kegg geneIDs
    if(length(kegg.gene) > 0){
        for(j in 1:length(kegg.gene)){
            idx <- strsplit(kegg.gene[j], "/")[[1]]
            data_upreg[[i]]$Kegg$geneID[j] <- paste(data_upreg[[i]]$Gene$SYMBOL
                [data_upreg[[i]]$Gene$ENTREZID %in% idx], collapse="/")
        }
    }
}

# Save results to file
system("mkdir upregulated")
for(i in names(data_upreg)){
    system(paste("mkdir ", getwd(), "/upregulated/", i, sep=""))
    for (j in names(data_upreg[[i]])){
        write.table(data_upreg[[i]][[j]], file=paste(getwd(), 
            "/upregulated/", i, "/", gsub("[.]", "_", i), "_", j, ".txt", sep=""),
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
    }
}
```

### Enrichment of Down-regulated Genes  
Gene set enrichment for downregulated genes


```r
down <- lapply(Rrst, function(x) x[x$logFC < 0,])
down <- lapply(down, function(x) 
    do.call(rbind, lapply(as.character(x$genes), 
    function(y) sigG.entrez[sigG.entrez$SYMBOL %in% y,])))
down <- down[!unlist(lapply(down, is.null))]
unlist(lapply(down, nrow))
```

```
##   Fat.Repletion.24h Sweet.Repletion.24h   Fat.Repletion.72h 
##                 292                  24                1536 
## Sweet.Repletion.72h 
##                 911
```

Gene enrichment analysis per contrast for downregulated genes:


```r
enrichCont.down <- lapply(down, 
    function(x) enrich(lst=x, bkg=bkg.entrez))
lapply(enrichCont.down, function(x)
    unlist(lapply(x, function(y) nrow(y))))
```

```
## $Fat.Repletion.24h
## Genes GO.BP GO.CC GO.MF  Kegg 
##   292     0     0     0     0 
## 
## $Sweet.Repletion.24h
## Genes GO.BP GO.CC GO.MF  Kegg 
##    24    22     0     0     0 
## 
## $Fat.Repletion.72h
## Genes GO.BP GO.CC GO.MF  Kegg 
##  1536    52    35    17     2 
## 
## $Sweet.Repletion.72h
## Genes GO.BP GO.CC GO.MF  Kegg 
##   911   107    33     1     1
```

Save results


```r
data_downreg <- lapply(names(enrichCont.down), function(x) 
    lapply(names(enrichCont.down[[x]]), function(y) 
        data.frame(enrichCont.down[[x]][[y]])))
names(data_downreg) <- names(enrichCont.down)
for(i in names(data_downreg)){
    # Add names to enrichment analysis per contrast
    names(data_downreg[[i]]) <- names(data_merged)
    kegg.gene <- data_downreg[[i]]$Kegg$geneID
    # Add gene names to kegg geneIDs
    if(length(kegg.gene) > 0){
        for(j in 1:length(kegg.gene)){
            idx <- strsplit(kegg.gene[j], "/")[[1]]
            data_downreg[[i]]$Kegg$geneID[j] <- paste(data_downreg[[i]]$Gene$SYMBOL
                [data_downreg[[i]]$Gene$ENTREZID %in% idx], collapse="/")
        }
    }
}

# Save results to file
system("mkdir downregulated")
for(i in names(data_downreg)){
    system(paste("mkdir ", getwd(), "/downregulated/", i, sep=""))
    for (j in names(data_downreg[[i]])){
        write.table(data_downreg[[i]][[j]], file=paste(getwd(), 
            "/downregulated/", i, "/", gsub("[.]", "_", i), "_", j, ".txt", sep=""),
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
    }
}
```

### Save gene set enrichment results  


```r
save(Rrst, sigG.entrez, bkg.entrez, enrich.rst.fit3, data_merged, 
    enrichCont.rst.fit3, data_contrast, enrichCont.up, data_upreg,
    enrichCont.down, data_downreg, file=paste(getwd(), 
    "Enrichment_KER_Glycogen_Depl.Rdata", sep="/"))
```

### Run R Script


```r
htmlRunR
Enrichment_KER_Glycogen_Depl.R nodes=1,cpus-per-task=1,time=03:00:00,mem=50G \
+Gene Set Enrichment KER Glycogen Diets Over Time
```

