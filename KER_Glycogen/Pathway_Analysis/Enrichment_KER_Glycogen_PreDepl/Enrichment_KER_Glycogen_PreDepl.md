---
title: Gene Set Enrichment KER Glycogen Diets Over Time
author: Deborah Velez-Irizarry
date: Tue Oct 30 15:36:10 EDT 2018
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---
### Description:  
Performe gene set enrichment analysis for KER Glycogen project. Look at enrichment of  
genes differentially expressed for each diet over time with pre depletion timepoint as  
reference.   
  
***  
**Code:**  
Parent Directory:  

>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq  
  
Directory/File:  
 
&nbsp;&nbsp;&nbsp;&nbsp;/Enrichment/Enrichment_KER_Glycogen_PreDepl/Enrichment_KER_Glycogen_PreDepl.R  
 
**Input files:**  
Directory/File:  
  
>&nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_KER_Glycogen/results_diet_over_time_pre_ref.Rdata  
>&nbsp;&nbsp;&nbsp;&nbsp;/Enrichment/Enrichment_KER_Glycogen_Depl/enrichment_function.Rdata
  
**Output files:**  
  
Directory:  
  
>&nbsp;&nbsp;&nbsp;&nbsp;/Enrichment/Enrichment_KER_Glycogen_PreDepl/  

Files:

>&nbsp;&nbsp;&nbsp;&nbsp;/Enrichment_KER_Glycogen_PreDepl.Rdata  
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
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq"
load(paste(dir, "DE/DE_KER_Glycogen", 
    "results_diet_over_time_pre_ref.Rdata", sep="/"))
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
Rrst <- lapply(Rst.Diet.Pre, function(x) x[x$adj.P.Val < 0.05,])
lapply(Rrst, nrow)
```

```
## $Fat.Depletion
## [1] 1197
## 
## $Sweet.Depletion
## [1] 918
## 
## $Fat.Repletion.24h
## [1] 72
## 
## $Sweet.Repletion.24h
## [1] 1612
## 
## $Fat.Repletion.72h
## [1] 6135
## 
## $Sweet.Repletion.72h
## [1] 3861
## 
## $Fat.over.time
## [1] 8232
## 
## $Sweet.over.time
## [1] 4016
```

Write DE results to file


```r
system("mkdir DE.Results")
lapply(names(Rrst), function(x) 
    write.table(Rrst[[x]], file=paste(getwd(), "/DE.Results/", gsub("[.]", "_", x), "_pre_ref.txt", sep=""),
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
## 
## [[7]]
## NULL
## 
## [[8]]
## NULL
```

Merge list of genes


```r
sigG <- do.call(rbind, lapply(Rrst, function(x) data.frame(Xloc=rownames(x), x[,c(1:5,7)])))
rownames(sigG) <- NULL
nrow(sigG)
```

```
## [1] 26043
```

Seperate genes without a name


```r
na.sigG <- sigG[is.na(sigG$genes),]
nrow(na.sigG)
```

```
## [1] 841
```

```r
sigG <- sigG[!is.na(sigG$genes),]
nrow(sigG)
```

```
## [1] 25202
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
## toType = c("ENTREZID"), : 13.67% of input gene IDs are fail to map...
```

```r
nrow(sigG.entrez)
```

```
## [1] 8084
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
enrich.rst.fit2 <- enrich(lst=sigG.entrez, bkg=bkg.entrez)
unlist(lapply(enrich.rst.fit2, nrow))
```

```
## Genes GO.BP GO.CC GO.MF  Kegg 
##  8084   126    32     0     0
```

Summary of significant GO terms


```r
lapply(enrich.rst.fit2, function(x) head(x[,2]))
```

```
## $Genes
## [1] "8013"  "5166"  "84879" "51129" "1734"  "1264" 
## 
## $GO.BP
## [1] "ATP metabolic process"                               
## [2] "ribonucleoside triphosphate metabolic process"       
## [3] "ribonucleoside monophosphate metabolic process"      
## [4] "oxidoreduction coenzyme metabolic process"           
## [5] "purine ribonucleoside triphosphate metabolic process"
## [6] "purine nucleoside monophosphate metabolic process"   
## 
## $GO.CC
## [1] "melanosome"                   "pigment granule"             
## [3] "oxidoreductase complex"       "mitochondrial membrane part" 
## [5] "mitochondrial inner membrane" "endoplasmic reticulum lumen" 
## 
## $GO.MF
## character(0)
## 
## $Kegg
## character(0)
```

Save enrichment analysis results to file


```r
# Enrichemnt analysis results
data_merged <- lapply(enrich.rst.fit2, function(x) data.frame(x))

# Add gene names to Kegg results (currently has the geneIDs)
kegg.gene <- strsplit(data_merged$Kegg$geneID, "/")[[1]]
```

```
## Error in strsplit(data_merged$Kegg$geneID, "/")[[1]]: subscript out of bounds
```

```r
data_merged$Kegg$geneID <- paste(data_merged$Gene$SYMBOL[data_merged$Gene$ENTREZID %in% kegg.gene], collapse="/")
```

```
## Error in data_merged$Gene$ENTREZID %in% kegg.gene: object 'kegg.gene' not found
```

```r
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
##       Fat.Depletion     Sweet.Depletion   Fat.Repletion.24h 
##                 980                 766                  63 
## Sweet.Repletion.24h   Fat.Repletion.72h Sweet.Repletion.72h 
##                1394                5269                3242 
##       Fat.over.time     Sweet.over.time 
##                6931                3335
```

Gene enrichment analysis per contrast:


```r
enrichCont.rst.fit2 <- lapply(enrichAll.rst.list, 
    function(x) enrich(lst=x, bkg=bkg.entrez))
lapply(enrichCont.rst.fit2, function(x)
    unlist(lapply(x, function(y) nrow(y))))
```

```
## $Fat.Depletion
## Genes GO.BP GO.CC GO.MF  Kegg 
##   980     0     0     0     0 
## 
## $Sweet.Depletion
## Genes GO.BP GO.CC GO.MF  Kegg 
##   766     3     0     0     0 
## 
## $Fat.Repletion.24h
## Genes GO.BP GO.CC GO.MF  Kegg 
##    63   137    12     2     5 
## 
## $Sweet.Repletion.24h
## Genes GO.BP GO.CC GO.MF  Kegg 
##  1394   118     8    18     1 
## 
## $Fat.Repletion.72h
## Genes GO.BP GO.CC GO.MF  Kegg 
##  5269   143    57    24    15 
## 
## $Sweet.Repletion.72h
## Genes GO.BP GO.CC GO.MF  Kegg 
##  3242   140    38    30    16 
## 
## $Fat.over.time
## Genes GO.BP GO.CC GO.MF  Kegg 
##  6931    48    14     2     1 
## 
## $Sweet.over.time
## Genes GO.BP GO.CC GO.MF  Kegg 
##  3335   177    25     4    10
```

Gene Ontology: Biological Processes (first 10)


```r
lapply(enrichCont.rst.fit2, function(x) 
    head(x[["GO.BP"]][,2]))
```

```
## $Fat.Depletion
## character(0)
## 
## $Sweet.Depletion
## [1] "regulation of fatty acid oxidation"         
## [2] "positive regulation of fatty acid oxidation"
## [3] "cofactor metabolic process"                 
## 
## $Fat.Repletion.24h
## [1] "regulation of skeletal muscle adaptation"
## [2] "muscle filament sliding"                 
## [3] "actin-myosin filament sliding"           
## [4] "regulation of muscle contraction"        
## [5] "skeletal muscle adaptation"              
## [6] "regulation of muscle system process"     
## 
## $Sweet.Repletion.24h
## [1] "response to cytokine"                    
## [2] "skeletal system development"             
## [3] "protein activation cascade"              
## [4] "innate immune response"                  
## [5] "complement activation, classical pathway"
## [6] "extracellular matrix organization"       
## 
## $Fat.Repletion.72h
## [1] "extracellular matrix organization"                       
## [2] "extracellular structure organization"                    
## [3] "response to wounding"                                    
## [4] "wound healing"                                           
## [5] "multicellular organismal macromolecule metabolic process"
## [6] "collagen metabolic process"                              
## 
## $Sweet.Repletion.72h
## [1] "cellular respiration"                               
## [2] "electron transport chain"                           
## [3] "generation of precursor metabolites and energy"     
## [4] "respiratory electron transport chain"               
## [5] "energy derivation by oxidation of organic compounds"
## [6] "ATP synthesis coupled electron transport"           
## 
## $Fat.over.time
## [1] "extracellular matrix organization"      
## [2] "extracellular structure organization"   
## [3] "innate immune response"                 
## [4] "collagen metabolic process"             
## [5] "response to cytokine"                   
## [6] "water-soluble vitamin metabolic process"
## 
## $Sweet.over.time
## [1] "generation of precursor metabolites and energy"
## [2] "ribonucleoside triphosphate metabolic process" 
## [3] "nucleoside triphosphate metabolic process"     
## [4] "glycosyl compound metabolic process"           
## [5] "ribonucleoside metabolic process"              
## [6] "ATP metabolic process"
```

Gene Ontology: Cellular Processes (first 10)


```r
lapply(enrichCont.rst.fit2, function(x) 
    head(x[["GO.CC"]][,2]))
```

```
## $Fat.Depletion
## character(0)
## 
## $Sweet.Depletion
## character(0)
## 
## $Fat.Repletion.24h
## [1] "proteasome complex"            "endopeptidase complex"        
## [3] "sarcomere"                     "striated muscle thin filament"
## [5] "peptidase complex"             "contractile fiber part"       
## 
## $Sweet.Repletion.24h
## [1] "proteinaceous extracellular matrix"
## [2] "extracellular matrix"              
## [3] "cell surface"                      
## [4] "collagen trimer"                   
## [5] "endoplasmic reticulum lumen"       
## [6] "secretory granule lumen"           
## 
## $Fat.Repletion.72h
## [1] "endoplasmic reticulum lumen"              
## [2] "cell surface"                             
## [3] "plasma membrane protein complex"          
## [4] "mitochondrial inner membrane"             
## [5] "protein complex involved in cell adhesion"
## [6] "integrin complex"                         
## 
## $Sweet.Repletion.72h
## [1] "mitochondrial inner membrane"  "organelle inner membrane"     
## [3] "mitochondrial membrane part"   "mitochondrial protein complex"
## [5] "oxidoreductase complex"        "respiratory chain"            
## 
## $Fat.over.time
## [1] "endoplasmic reticulum lumen"      "cell surface"                    
## [3] "sarcoplasmic reticulum"           "collagen trimer"                 
## [5] "external side of plasma membrane" "membrane region"                 
## 
## $Sweet.over.time
## [1] "mitochondrial inner membrane"  "organelle inner membrane"     
## [3] "mitochondrial membrane part"   "mitochondrial protein complex"
## [5] "mitochondrial matrix"          "mitochondrial outer membrane"
```

Gene Ontology: Molecular Processes (first 10)


```r
lapply(enrichCont.rst.fit2, function(x) 
    head(x[["GO.MP"]][,2]))
```

```
## $Fat.Depletion
## NULL
## 
## $Sweet.Depletion
## NULL
## 
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
lapply(enrichCont.rst.fit2, function(x) 
    head(x[["Kegg"]][,2]))
```

```
## $Fat.Depletion
## character(0)
## 
## $Sweet.Depletion
## character(0)
## 
## $Fat.Repletion.24h
## [1] "Proteasome"                            
## [2] "Cardiac muscle contraction"            
## [3] "Hypertrophic cardiomyopathy (HCM)"     
## [4] "Dilated cardiomyopathy (DCM)"          
## [5] "Adrenergic signaling in cardiomyocytes"
## 
## $Sweet.Repletion.24h
## [1] "Complement and coagulation cascades"
## 
## $Fat.Repletion.72h
## [1] "ECM-receptor interaction"                 
## [2] "Primary immunodeficiency"                 
## [3] "Carbon metabolism"                        
## [4] "Natural killer cell mediated cytotoxicity"
## [5] "Focal adhesion"                           
## [6] "Glycolysis / Gluconeogenesis"             
## 
## $Sweet.Repletion.72h
## [1] "Oxidative phosphorylation"                
## [2] "Parkinson disease"                        
## [3] "Huntington disease"                       
## [4] "Non-alcoholic fatty liver disease (NAFLD)"
## [5] "Alzheimer disease"                        
## [6] "Thermogenesis"                            
## 
## $Fat.over.time
## [1] "ECM-receptor interaction"
## 
## $Sweet.over.time
## [1] "Huntington disease"                       
## [2] "Non-alcoholic fatty liver disease (NAFLD)"
## [3] "Complement and coagulation cascades"      
## [4] "Oxidative phosphorylation"                
## [5] "Proteasome"                               
## [6] "Epstein-Barr virus infection"
```

Save results


```r
data_contrast <- lapply(names(enrichCont.rst.fit2), function(x) 
    lapply(names(enrichCont.rst.fit2[[x]]), function(y) 
        data.frame(enrichCont.rst.fit2[[x]][[y]])))
names(data_contrast) <- names(enrichCont.rst.fit2)
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
##       Fat.Depletion     Sweet.Depletion   Fat.Repletion.24h 
##                 574                 362                  36 
## Sweet.Repletion.24h   Fat.Repletion.72h Sweet.Repletion.72h 
##                 811                2962                1737
```

Gene enrichment analysis per contrast for upregulated genes:


```r
enrichCont.up <- lapply(up, 
    function(x) enrich(lst=x, bkg=bkg.entrez))
lapply(enrichCont.up, function(x)
    unlist(lapply(x, function(y) nrow(y))))
```

```
## $Fat.Depletion
## Genes GO.BP GO.CC GO.MF  Kegg 
##   574     0     1     0     1 
## 
## $Sweet.Depletion
## Genes GO.BP GO.CC GO.MF  Kegg 
##   362    12     2     0     4 
## 
## $Fat.Repletion.24h
## Genes GO.BP GO.CC GO.MF  Kegg 
##    36   112     9     4     1 
## 
## $Sweet.Repletion.24h
## Genes GO.BP GO.CC GO.MF  Kegg 
##   811   245    35    31     8 
## 
## $Fat.Repletion.72h
## Genes GO.BP GO.CC GO.MF  Kegg 
##  2962   881   100    85    69 
## 
## $Sweet.Repletion.72h
## Genes GO.BP GO.CC GO.MF  Kegg 
##  1737   333    68    58    36
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
##       Fat.Depletion     Sweet.Depletion   Fat.Repletion.24h 
##                 406                 404                  27 
## Sweet.Repletion.24h   Fat.Repletion.72h Sweet.Repletion.72h 
##                 583                2307                1505
```

Gene enrichment analysis per contrast for downregulated genes:


```r
enrichCont.down <- lapply(down, 
    function(x) enrich(lst=x, bkg=bkg.entrez))
lapply(enrichCont.down, function(x)
    unlist(lapply(x, function(y) nrow(y))))
```

```
## $Fat.Depletion
## Genes GO.BP GO.CC GO.MF  Kegg 
##   406     0     0     0     0 
## 
## $Sweet.Depletion
## Genes GO.BP GO.CC GO.MF  Kegg 
##   404     0     0     0    16 
## 
## $Fat.Repletion.24h
## Genes GO.BP GO.CC GO.MF  Kegg 
##    27    60    16     1     5 
## 
## $Sweet.Repletion.24h
## Genes GO.BP GO.CC GO.MF  Kegg 
##   583    96    24    12    26 
## 
## $Fat.Repletion.72h
## Genes GO.BP GO.CC GO.MF  Kegg 
##  2307   119    64    20    18 
## 
## $Sweet.Repletion.72h
## Genes GO.BP GO.CC GO.MF  Kegg 
##  1505   151    39    30    23
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
save(Rrst, sigG.entrez, bkg.entrez, enrich.rst.fit2, data_merged, 
    enrichCont.rst.fit2, data_contrast, enrichCont.up, data_upreg,
    enrichCont.down, data_downreg, file=paste(getwd(), 
    "Enrichment_KER_Glycogen_PreDepl.Rdata", sep="/"))
```

### Run R Script


```r
htmlRunR
Enrichment_KER_Glycogen_PreDepl.R nodes=1,cpus-per-task=1,time=03:00:00,mem=50G \
+Gene Set Enrichment KER Glycogen Diets Over Time
```

