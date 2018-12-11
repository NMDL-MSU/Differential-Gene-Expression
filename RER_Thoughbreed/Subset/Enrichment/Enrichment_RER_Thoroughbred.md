---
title: Gene Set Enrichment RER Thoroughbred Project
author: Deborah Velez-Irizarry
date: Tue Dec 11 11:39:17 EST 2018
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---


```r
### Description:  
```

Performe gene set enrichment analysis for RER Thoroughbred project using the subset of animals with proteomics  
data. Look at enrichment of genes differentially expressed in horses presenting recurrent exertional  
rhabdomyolysis (RER) compared to controls.  
  
***  
**Code:**  
Parent Directory:  

>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq  
  
Directory/File:  
 
&nbsp;&nbsp;&nbsp;&nbsp;/Enrichment/Enrichment_RER_Thoroughbred/Enrichment_RER_Thoroughbred.R  
 
**Input files:**  
Directory/File:  
  
>&nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_RER_Thoroughbred/DE_RER_Thoroughbred_Protein_Subset/Protein_Subset.Rdata  
>&nbsp;&nbsp;&nbsp;&nbsp;/Enrichment/Enrichment_KER_Glycogen_Depl/enrichment_function.Rdata
  
**Output files:**  
  
Directory:  
  
>&nbsp;&nbsp;&nbsp;&nbsp;/Enrichment/Enrichment_RER_Thoroughbred/  

Files:

>&nbsp;&nbsp;&nbsp;&nbsp;/Merged/*.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;/UpRegulated/*.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;/DownRegulated/*.txt  
 
***  
### R Environment    
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
## DOSE v3.8.0  For help: https://guangchuangyu.github.io/DOSE
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
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind,
##     colMeans, colnames, colSums, dirname, do.call, duplicated,
##     eval, evalq, Filter, Find, get, grep, grepl, intersect,
##     is.unsorted, lapply, lengths, Map, mapply, match, mget, order,
##     paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind,
##     Reduce, rowMeans, rownames, rowSums, sapply, setdiff, sort,
##     table, tapply, union, unique, unsplit, which, which.max,
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
## The following object is masked from 'package:base':
## 
##     expand.grid
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
## clusterProfiler v3.10.0  For help: https://guangchuangyu.github.io/software/clusterProfiler
## 
## If you use clusterProfiler in published research, please cite:
## Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.
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

Session Information


```r
sessionInfo()
```

```
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS/LAPACK: /opt/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib/libopenblas_sandybridgep-r0.3.1.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] GenomicFeatures_1.34.1 GenomicRanges_1.34.0   GenomeInfoDb_1.18.1   
##  [4] org.Hs.eg.db_3.7.0     clusterProfiler_3.10.0 GSEABase_1.44.0       
##  [7] graph_1.60.0           annotate_1.60.0        XML_3.98-1.16         
## [10] GO.db_3.7.0            AnnotationDbi_1.44.0   IRanges_2.16.0        
## [13] S4Vectors_0.20.1       Biobase_2.42.0         BiocGenerics_0.28.0   
## [16] DOSE_3.8.0             edgeR_3.24.1           limma_3.38.3          
## [19] knitr_1.20            
## 
## loaded via a namespace (and not attached):
##  [1] fgsea_1.8.0                 colorspace_1.3-2           
##  [3] ggridges_0.5.1              qvalue_2.14.0              
##  [5] XVector_0.22.0              farver_1.1.0               
##  [7] urltools_1.7.0              ggrepel_0.8.0              
##  [9] bit64_0.9-7                 xml2_1.2.0                 
## [11] splines_3.5.1               GOSemSim_2.8.0             
## [13] jsonlite_1.5                Rsamtools_1.34.0           
## [15] ggforce_0.1.3               compiler_3.5.1             
## [17] httr_1.3.1                  rvcheck_0.1.3              
## [19] assertthat_0.2.0            Matrix_1.2-14              
## [21] lazyeval_0.2.1              tweenr_1.0.0               
## [23] prettyunits_1.0.2           tools_3.5.1                
## [25] bindrcpp_0.2.2              igraph_1.2.1               
## [27] gtable_0.2.0                glue_1.3.0                 
## [29] GenomeInfoDbData_1.2.0      reshape2_1.4.3             
## [31] DO.db_2.9                   dplyr_0.7.8                
## [33] fastmatch_1.1-0             Rcpp_1.0.0                 
## [35] enrichplot_1.2.0            Biostrings_2.50.1          
## [37] rtracklayer_1.42.1          ggraph_1.0.2               
## [39] stringr_1.3.1               europepmc_0.3              
## [41] MASS_7.3-51.1               zlibbioc_1.28.0            
## [43] scales_1.0.0                hms_0.4.2                  
## [45] SummarizedExperiment_1.12.0 RColorBrewer_1.1-2         
## [47] memoise_1.1.0               gridExtra_2.3              
## [49] ggplot2_3.1.0.9000          UpSetR_1.3.3               
## [51] biomaRt_2.38.0              triebeard_0.3.0            
## [53] stringi_1.2.3               RSQLite_2.1.1              
## [55] BiocParallel_1.16.2         rlang_0.3.0.1              
## [57] pkgconfig_2.0.2             bitops_1.0-6               
## [59] matrixStats_0.54.0          evaluate_0.12              
## [61] lattice_0.20-38             purrr_0.2.5                
## [63] bindr_0.1.1                 GenomicAlignments_1.18.0   
## [65] cowplot_0.9.2               bit_1.1-14                 
## [67] tidyselect_0.2.5            plyr_1.8.4                 
## [69] magrittr_1.5                R6_2.3.0                   
## [71] DelayedArray_0.8.0          DBI_1.0.0                  
## [73] pillar_1.2.3                units_0.6-2                
## [75] RCurl_1.95-4.11             tibble_1.4.2               
## [77] crayon_1.3.4                viridis_0.5.1              
## [79] progress_1.2.0              locfit_1.5-9.1             
## [81] grid_3.5.1                  data.table_1.11.8          
## [83] blob_1.1.1                  digest_0.6.18              
## [85] xtable_1.8-2                tidyr_0.8.1                
## [87] gridGraphics_0.3-0          munsell_0.5.0              
## [89] viridisLite_0.3.0           ggplotify_0.0.3
```

### Load Data for Pathway Analysis
Load DE results


```r
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq"
load(paste(dir, 
    "/DE/DE_RER_Thoroughbred/DE_RER_Thoroughbred_Protein_Subset/Protein_Subset.Rdata", 
    sep=""))
```

Load enrichment function


```r
load(paste(dir, "Enrichment/Enrichment_KER_Glycogen_Depl", 
    "enrichment_function.Rdata", sep="/"))
```

### Animal Groups


```r
group <- anim$Dx
table(group)
```

```
## group
## Control     RER 
##       5       5
```

### Differential Expression Analysis Results
Reduce results list to contain only differentially expressed genes


```r
Rrst <- rst[rst$Q.Val < 0.05, c(1:5, 7, 9, 11:17)]
nrow(Rrst)
```

```
## [1] 812
```

Write DE results to file


```r
write.table(Rrst, file=,"DE_RER_Thoroughbred_Protein_Subset.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
```

Significant gene names


```r
sigG <- as.character(Rrst$genes)
names(sigG) <- rownames(Rrst)
```

Seperate genes without a name


```r
# Number of unknown gene transcripts
na.sigG <- sigG[is.na(sigG)]
length(na.sigG)
```

```
## [1] 32
```

```r
# Number of known genes
sigG <- sigG[!is.na(sigG)]
length(sigG)
```

```
## [1] 780
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
## [1] 789
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
sigG.entrez <- bitr(unique(sigG), fromType="SYMBOL", 
    toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
```

```
## 'select()' returned 1:1 mapping between keys and columns
```

```
## Warning in bitr(unique(sigG), fromType = "SYMBOL", toType =
## c("ENTREZID"), : 14.62% of input gene IDs are fail to map...
```

```r
nrow(sigG.entrez)
```

```
## [1] 666
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
## toType = c("ENTREZID"), : 14.77% of input gene IDs are fail to map...
```

```r
nrow(bkg.entrez)
```

```
## [1] 11095
```

### Enrichment Analysis  
Global gene enrichment analysis: Merge all DE genes


```r
enrich.rst <- enrich(lst=sigG.entrez, bkg=bkg.entrez)
unlist(lapply(enrich.rst, nrow))
```

```
## Genes GO.BP GO.CC GO.MF  Kegg 
##   666    88    46    15     8
```

Summary of significant GO terms


```r
lapply(enrich.rst, function(x) head(x[,2]))
```

```
## $Genes
## [1] "581"    "4647"   "158135" "10439"  "22980"  "55785" 
## 
## $GO.BP
## [1] "SRP-dependent cotranslational protein targeting to membrane"        
## [2] "cotranslational protein targeting to membrane"                      
## [3] "protein targeting to ER"                                            
## [4] "establishment of protein localization to endoplasmic reticulum"     
## [5] "nuclear-transcribed mRNA catabolic process, nonsense-mediated decay"
## [6] "protein localization to endoplasmic reticulum"                      
## 
## $GO.CC
## [1] "cytosolic ribosome"                "ribosomal subunit"                
## [3] "ribosome"                          "cytosolic part"                   
## [5] "cytosolic large ribosomal subunit" "large ribosomal subunit"          
## 
## $GO.MF
## [1] "structural constituent of ribosome"      
## [2] "structural molecule activity"            
## [3] "rRNA binding"                            
## [4] "NADH dehydrogenase (ubiquinone) activity"
## [5] "NADH dehydrogenase (quinone) activity"   
## [6] "NADH dehydrogenase activity"             
## 
## $Kegg
## [1] "Ribosome"                  "Oxidative phosphorylation"
## [3] "Huntington disease"        "Parkinson disease"        
## [5] "Alzheimer disease"         "Thermogenesis"
```

Save enrichment analysis results to file


```r
# Enrichemnt analysis results
data_merged <- lapply(enrich.rst, function(x) data.frame(x))

# Add gene names to Kegg results (currently has the geneIDs)
kegg.gene <- sapply(data_merged$Kegg$geneID, function(x) strsplit(x, "/")[[1]])
names(kegg.gene) <- NULL
data_merged$Kegg$geneID <- sapply(kegg.gene, function(x) 
    paste(data_merged$Gene$SYMBOL[data_merged$Gene$ENTREZID %in% x], collapse="/"))

# Save results to file
system("mkdir Merged")
z <- lapply(names(data_merged), function(x) 
    write.table(data_merged[[x]], file=paste(getwd(), "/Merged/", x, "_merged.txt", sep=""),
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t"))
```

### Enrichment of Up-regulated Genes  
Gene set enrichment for upregulated genes


```r
up <- Rrst[Rrst$logFC > 0,]
up <- lapply(as.character(up$genes), 
    function(x) sigG.entrez[sigG.entrez$SYMBOL %in% x,])
up <- do.call(rbind, up[unlist(lapply(up, nrow)) > 0])
nrow(up)
```

```
## [1] 400
```

Gene enrichment analysis for upregulated genes:


```r
enrich.upreg <- enrich(lst=up, bkg=bkg.entrez)
unlist(lapply(enrich.upreg, nrow))
```

```
## Genes GO.BP GO.CC GO.MF  Kegg 
##   400   104    48    22     8
```

Summary of significant GO terms


```r
lapply(enrich.upreg, function(x) head(x[,2]))
```

```
## $Genes
## [1] "158135" "10439"  "27129"  "64077"  "8839"   "3491"  
## 
## $GO.BP
## [1] "SRP-dependent cotranslational protein targeting to membrane"        
## [2] "cotranslational protein targeting to membrane"                      
## [3] "protein targeting to ER"                                            
## [4] "establishment of protein localization to endoplasmic reticulum"     
## [5] "nuclear-transcribed mRNA catabolic process, nonsense-mediated decay"
## [6] "protein localization to endoplasmic reticulum"                      
## 
## $GO.CC
## [1] "ribosomal subunit"                 "cytosolic ribosome"               
## [3] "ribosome"                          "cytosolic part"                   
## [5] "cytosolic large ribosomal subunit" "large ribosomal subunit"          
## 
## $GO.MF
## [1] "structural constituent of ribosome"      
## [2] "structural molecule activity"            
## [3] "rRNA binding"                            
## [4] "NADH dehydrogenase (ubiquinone) activity"
## [5] "NADH dehydrogenase (quinone) activity"   
## [6] "NADH dehydrogenase activity"             
## 
## $Kegg
## [1] "Ribosome"                  "Oxidative phosphorylation"
## [3] "Huntington disease"        "Parkinson disease"        
## [5] "Alzheimer disease"         "Thermogenesis"
```

Save results


```r
data_upreg <- lapply(enrich.upreg, function(x)  
        data.frame(x))
names(data_upreg) <- names(enrich.upreg)

# Add gene names to Kegg results (currently has the geneIDs)
kegg.gene <- sapply(data_upreg$Kegg$geneID, function(x) strsplit(x, "/")[[1]])
names(kegg.gene) <- NULL
data_upreg$Kegg$geneID <- sapply(kegg.gene, function(x) 
    paste(data_upreg$Gene$SYMBOL[data_upreg$Gene$ENTREZID %in% x], collapse="/"))


# Save results to file
system("mkdir UpRegulated")
z <- lapply(names(data_upreg), function(x) 
    write.table(data_upreg[[x]], file=paste(getwd(), "/UpRegulated/", x, "_upregulated.txt", sep=""),
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t"))
```

### Enrichment of Down-regulated Genes  
Gene set enrichment for down-regulated genes


```r
down <- Rrst[Rrst$logFC < 0,]
down <- lapply(as.character(down$genes), 
    function(x) sigG.entrez[sigG.entrez$SYMBOL %in% x,])
down <- do.call(rbind, down[unlist(lapply(down, nrow)) > 0])
nrow(down)
```

```
## [1] 266
```

Gene enrichment analysis for down-regulated genes:


```r
enrich.downreg <- enrich(lst=down, bkg=bkg.entrez)
unlist(lapply(enrich.downreg, nrow))
```

```
## Genes GO.BP GO.CC GO.MF  Kegg 
##   266     0     0     1     1
```

Summary of significant GO terms


```r
lapply(enrich.downreg, function(x) head(x[,2]))
```

```
## $Genes
## [1] "581"   "4647"  "22980" "55785" "10964" "57674"
## 
## $GO.BP
## character(0)
## 
## $GO.CC
## character(0)
## 
## $GO.MF
## [1] "histone methyltransferase activity (H3-K4 specific)"
## 
## $Kegg
## [1] "Lysine degradation"
```

Save results


```r
data_downreg <- lapply(enrich.downreg, function(x)  
        data.frame(x))
names(data_downreg) <- names(enrich.downreg)
data_downreg <- data_downreg[unlist(lapply(data_downreg, nrow)) > 0]

# Add gene names to Kegg results (currently has the geneIDs)
kegg.gene <- strsplit(data_downreg$Kegg$geneID, "/")[[1]]
data_downreg$Kegg$geneID <- paste(data_downreg$Gene$SYMBOL[data_downreg$Gene$ENTREZID %in% kegg.gene], collapse="/")

# Save results to file
system("mkdir DownRegulated")
z <- lapply(names(data_downreg), function(x) 
    write.table(data_downreg[[x]], file=paste(getwd(), "/DownRegulated/", x, "_downregulated.txt", sep=""),
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t"))
```

### Run R Script


```r
htmlRunR
Enrichment_RER_Thoroughbred.R nodes=1,cpus-per-task=1,time=03:00:00,mem=50G \
+Gene Set Enrichment RER Thoroughbred Project
```

