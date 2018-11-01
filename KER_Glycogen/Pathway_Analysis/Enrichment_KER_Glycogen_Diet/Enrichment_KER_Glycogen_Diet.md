---
title: Gene Set Enrichment KER Glycogen Diets Over Time
author: Deborah Velez-Irizarry
date: Tue Oct 30 16:24:33 EDT 2018
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---
### Description:  
Performe gene set enrichment analysis for KER Glycogen project. Look at enrichment of  
genes differentially expressed between diets over time with the Fat diet as reference.  
  
***  
**Code:**  
Parent Directory:  

>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq  
  
Directory/File:  
 
&nbsp;&nbsp;&nbsp;&nbsp;/Enrichment/Enrichment_KER_Glycogen_Diet/Enrichment_KER_Glycogen_Diet.R  
 
**Input files:**  
Directory/File:  
  
>&nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_KER_Glycogen/results_diff_between_diet_over_time.Rdata  
>&nbsp;&nbsp;&nbsp;&nbsp;/Enrichment/Enrichment_KER_Glycogen_Depl/enrichment_function.Rdata
  
**Output files:**  
  
Directory:  
  
>&nbsp;&nbsp;&nbsp;&nbsp;/Enrichment/Enrichment_KER_Glycogen_Diet/  

Files:

>&nbsp;&nbsp;&nbsp;&nbsp;/Enrichment_KER_Glycogen_Diet.Rdata  
>&nbsp;&nbsp;&nbsp;&nbsp;/Pre.Depletion/*.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;/Repletion.24h/*.txt
>&nbsp;&nbsp;&nbsp;&nbsp;/Repletion.72h/*.txt
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
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq"
load(paste(dir, "DE/DE_KER_Glycogen", 
    "results_diff_between_diet_over_time.Rdata", sep="/"))
```

Load enrichment function


```r
load(paste(dir, "Enrichment/Enrichment_KER_Glycogen_Depl", 
    "enrichment_function.Rdata", sep="/"))
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
Rrst <- lapply(Rst.Timepoint, function(x) x[x$adj.P.Val < 0.05,])
lapply(Rrst, nrow)
```

```
## $Pre.Depletion
## [1] 1
## 
## $Depletion
## [1] 0
## 
## $Repletion.24h
## [1] 3
## 
## $Repletion.72h
## [1] 742
```

Write DE results to file


```r
system("mkdir DE.Results")
lapply(names(Rrst), function(x) 
    write.table(Rrst[[x]], file=paste(getwd(), "/DE.Results/", gsub("[.]", "_", x), "_fat_ref.txt", sep=""),
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
```

Merge list of genes


```r
sigG <- do.call(rbind, lapply(Rrst, function(x) data.frame(Xloc=rownames(x), x[,c(1:5,7)])))
rownames(sigG) <- NULL
nrow(sigG)
```

```
## [1] 746
```

Seperate genes without a name


```r
na.sigG <- sigG[is.na(sigG$genes),]
nrow(na.sigG)
```

```
## [1] 13
```

```r
sigG <- sigG[!is.na(sigG$genes),]
nrow(sigG)
```

```
## [1] 733
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
## toType = c("ENTREZID"), : 10.25% of input gene IDs are fail to map...
```

```r
nrow(sigG.entrez)
```

```
## [1] 658
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

### Global Enrichment Analysis  
Global gene enrichment analysis: Merge all DE genes


```r
enrich.rst.fit1 <- enrich(lst=sigG.entrez, bkg=bkg.entrez)
unlist(lapply(enrich.rst.fit1, nrow))
```

```
## Genes GO.BP GO.CC GO.MF  Kegg 
##   658   209     9     0     2
```

Summary of significant GO terms


```r
lapply(enrich.rst.fit1, function(x) head(x[,2]))
```

```
## $Genes
## [1] "51129"  "23329"  "54625"  "57561"  "6696"   "219972"
## 
## $GO.BP
## [1] "antigen processing and presentation of exogenous peptide antigen via MHC class I, TAP-dependent"
## [2] "antigen processing and presentation of exogenous peptide antigen via MHC class I"               
## [3] "regulation of cellular amino acid metabolic process"                                            
## [4] "innate immune response"                                                                         
## [5] "antigen processing and presentation of peptide antigen via MHC class I"                         
## [6] "regulation of cellular amine metabolic process"                                                 
## 
## $GO.CC
## [1] "proteasome complex"                             
## [2] "endopeptidase complex"                          
## [3] "peptidase complex"                              
## [4] "proteasome regulatory particle"                 
## [5] "proteasome accessory complex"                   
## [6] "proteasome regulatory particle, base subcomplex"
## 
## $GO.MF
## character(0)
## 
## $Kegg
## [1] "Proteasome"                   "Epstein-Barr virus infection"
```

Save enrichment analysis results to file


```r
# Enrichemnt analysis results
data_merged <- lapply(enrich.rst.fit1, function(x) data.frame(x))

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
enrichAll.rst.list <- enrichAll.rst.list[!unlist(lapply(enrichAll.rst.list, is.null))]
enrichAll.rst.list <- enrichAll.rst.list[-1]
unlist(lapply(enrichAll.rst.list, nrow))
```

```
## Repletion.24h Repletion.72h 
##             3           654
```

Gene enrichment analysis per contrast:


```r
enrichCont.rst.fit1 <- lapply(enrichAll.rst.list, 
    function(x) enrich(lst=x, bkg=bkg.entrez))
```

```
## No gene can be mapped....
```

```
## --> return NULL...
```

```r
lapply(enrichCont.rst.fit1, function(x)
    unlist(lapply(x, function(y) nrow(y))))
```

```
## $Repletion.24h
## Genes GO.BP GO.CC GO.MF 
##     3     0     0     2 
## 
## $Repletion.72h
## Genes GO.BP GO.CC GO.MF  Kegg 
##   654   210     9     0     2
```

Gene Ontology: Biological Processes (first 10)


```r
lapply(enrichCont.rst.fit1, function(x) 
    head(x[["GO.BP"]][,2]))
```

```
## $Repletion.24h
## character(0)
## 
## $Repletion.72h
## [1] "antigen processing and presentation of exogenous peptide antigen via MHC class I, TAP-dependent"
## [2] "antigen processing and presentation of exogenous peptide antigen via MHC class I"               
## [3] "regulation of cellular amino acid metabolic process"                                            
## [4] "innate immune response"                                                                         
## [5] "antigen processing and presentation of peptide antigen via MHC class I"                         
## [6] "regulation of cellular amine metabolic process"
```

Gene Ontology: Cellular Processes (first 10)


```r
lapply(enrichCont.rst.fit1, function(x) 
    head(x[["GO.CC"]][,2]))
```

```
## $Repletion.24h
## character(0)
## 
## $Repletion.72h
## [1] "proteasome complex"                             
## [2] "endopeptidase complex"                          
## [3] "peptidase complex"                              
## [4] "proteasome regulatory particle"                 
## [5] "proteasome accessory complex"                   
## [6] "proteasome regulatory particle, base subcomplex"
```

Gene Ontology: Molecular Processes (first 10)


```r
lapply(enrichCont.rst.fit1, function(x) 
    head(x[["GO.MP"]][,2]))
```

```
## $Repletion.24h
## NULL
## 
## $Repletion.72h
## NULL
```

Kegg Pathways:


```r
lapply(enrichCont.rst.fit1, function(x) 
    head(x[["Kegg"]][,2]))
```

```
## $Repletion.24h
## NULL
## 
## $Repletion.72h
## [1] "Proteasome"                   "Epstein-Barr virus infection"
```

Save results


```r
data_contrast <- lapply(names(enrichCont.rst.fit1), function(x) 
    lapply(names(enrichCont.rst.fit1[[x]]), function(y) 
        data.frame(enrichCont.rst.fit1[[x]][[y]])))
names(data_contrast) <- names(enrichCont.rst.fit1)
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
## Repletion.24h Repletion.72h 
##             2           270
```

Gene enrichment analysis per contrast for upregulated genes:


```r
enrichCont.up <- lapply(up, 
    function(x) enrich(lst=x, bkg=bkg.entrez))
```

```
## No gene can be mapped....
```

```
## --> return NULL...
```

```r
lapply(enrichCont.up, function(x)
    unlist(lapply(x, function(y) nrow(y))))
```

```
## $Repletion.24h
## Genes GO.BP GO.CC GO.MF 
##     2    18     0     5 
## 
## $Repletion.72h
## Genes GO.BP GO.CC GO.MF  Kegg 
##   270   128     8     9     2
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
## Pre.Depletion Repletion.24h Repletion.72h 
##             1             1           384
```

Gene enrichment analysis per contrast for downregulated genes:


```r
enrichCont.down <- lapply(down, 
    function(x) enrich(lst=x, bkg=bkg.entrez))
```

```
## No gene can be mapped....
```

```
## --> return NULL...
```

```r
lapply(enrichCont.down, function(x)
    unlist(lapply(x, function(y) nrow(y))))
```

```
## $Pre.Depletion
## Genes GO.BP GO.CC GO.MF  Kegg 
##     1    25     3     1     2 
## 
## $Repletion.24h
## Genes GO.BP GO.CC GO.MF 
##     1    27     4     8 
## 
## $Repletion.72h
## Genes GO.BP GO.CC GO.MF  Kegg 
##   384   202     2     2     5
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
save(Rrst, sigG.entrez, bkg.entrez, enrich.rst.fit1, data_merged, 
    enrichCont.rst.fit1, data_contrast, enrichCont.up, data_upreg,
    enrichCont.down, data_downreg, file=paste(getwd(), 
    "Enrichment_KER_Glycogen_Diet.Rdata", sep="/"))
```

### Run R Script


```r
htmlRunR
Enrichment_KER_Glycogen_Diet.R nodes=1,cpus-per-task=1,time=03:00:00,mem=50G \
+Gene Set Enrichment KER Glycogen Diets Over Time
```

