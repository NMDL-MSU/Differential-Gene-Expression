---
title: Differential Gene Expression on Subset
author: Deborah Velez-Irizarry
date: Wed Dec 5 16:38:29 EST 2018
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---
### Description:  
Differential gene expression analysis for subset of animals used in proteomic study.  
  
***  
**Code:**  
Parent Directory:  

> &nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq  
  
Directory/File:  
 
> &nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_KER_Glycogen/DE_RER_Thoroughbred_Protein_Subset/DE_RER_Thoroughbred_Protein_Subset.R  
 
**Input files:**  
Directory/File:  
  
> &nbsp;&nbsp;&nbsp;&nbsp;/HTSeq/htseq_counts_RER_Thoroughbred.txt  
> &nbsp;&nbsp;&nbsp;&nbsp;/RER_Thoroughbred_Animal_Information.txt  
> &nbsp;&nbsp;&nbsp;Cufflinks/MergedGTF/Annotation/annotation.txt  
  
**Output files:**  
  
Directory:  
  
> &nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_KER_Glycogen/DE_RER_Thoroughbred_Protein_Subset  
  
Files:  
  
> &nbsp;&nbsp;&nbsp;&nbsp;Protein_Subset.txt  
> &nbsp;&nbsp;&nbsp;&nbsp;Protein_Subset.Rdata  
 
Render R Script  
  
> &nbsp;&nbsp;&nbsp;&nbsp;DE_RER_Thoroughbred_Protein_Subset.qsub  
 
***  
### R Environment  
Clear Environment


```r
rm(list=ls())
```

Required Packages 


```r
library (limma)
library (edgeR)
library(qvalue)
```

**Session Information**


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
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] qvalue_2.14.0 edgeR_3.24.0  limma_3.38.2  knitr_1.20   
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.0         bindr_0.1.1        magrittr_1.5      
##  [4] splines_3.5.1      tidyselect_0.2.5   munsell_0.5.0     
##  [7] colorspace_1.3-2   lattice_0.20-38    R6_2.3.0          
## [10] rlang_0.3.0.1      stringr_1.3.1      plyr_1.8.4        
## [13] dplyr_0.7.8        tools_3.5.1        grid_3.5.1        
## [16] gtable_0.2.0       lazyeval_0.2.1     assertthat_0.2.0  
## [19] tibble_1.4.2       bindrcpp_0.2.2     reshape2_1.4.3    
## [22] purrr_0.2.5        ggplot2_3.1.0.9000 glue_1.3.0        
## [25] evaluate_0.12      stringi_1.2.3      compiler_3.5.1    
## [28] pillar_1.2.3       scales_1.0.0       locfit_1.5-9.1    
## [31] pkgconfig_2.0.2
```

### Load required R objects
Gene Counts


```r
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq"
counts <- read.table(paste(dir, "HTSeq", "htseq_counts_RER_Thoroughbred.txt", sep="/"))
colnames(counts) <- unlist(lapply(strsplit(colnames(counts), "X"), function(x) x[2]))
dim(counts)
```

```
## [1] 14155    23
```

Annotation


```r
annot <- read.table(paste(dir, "Cufflinks/MergedGTF/Annotation/annotation.txt", sep="/"),
    header=TRUE, row.names=8)
dim(annot)
```

```
## [1] 37870     8
```

Animal Information


```r
anim <- read.table(paste(dir, "RER_Thoroughbred_Animal_Information.txt", sep="/"),
    header=TRUE, row.names=2, sep="\t")[,-1]
anim$Sex <- rep("F", nrow(anim))
anim$Dx <- as.factor(anim$Dx)
anim$Age <- as.factor(anim$Age)
```

 ### Summary Function


```r
summSD <- function(x, dig=4) round(c(summary(x),
     Std.Dev.=sd(x)), dig)[c("Min.", "1st Qu.", "Median", "Mean", 
    "Std.Dev.", "3rd Qu.", "Max.")]
```

> **Data Check**: Animal IDs match between count matrix and animal matrix


```r
sum(!rownames(anim) == colnames(counts))
```

```
## [1] 0
```

### Prepare Data for DE Analysis:  
Retain gene expression on subset of animals used in proteomic analysis


```r
# Animals to retain
idx <- c("12613", "12910", "12915", "12916", "12918", "12401", "12402", "12403", "12620", "12913")

# Subset animal data frame
anim <- anim[idx,]
anim
```

```
##       Age Sex Last_Work AST  CK Biopsy_Changes      Dx Owner_Vet
## 12613   4   F       3.0 263 236              0 Control    Fenger
## 12910   4   F       1.5 297 172              0 Control    Fenger
## 12915   3   F       4.0 267 342              0 Control    Fenger
## 12916   2   F       5.0 430 243              0 Control    Fenger
## 12918   2   F       5.0 537 208              0 Control    Fenger
## 12401   7   F       2.0 921 663              2     RER Slaughter
## 12402   3   F       3.0 448 521              2     RER     Tores
## 12403   3   F       3.0 280 181              3     RER Slaughter
## 12620   4   F       3.0 219 131              2     RER    Fenger
## 12913   3   F       0.0 919 251              0     RER  Marshall
##                  Tx               Name
## 12613          <NA>       SkyHighSugar
## 12910          <NA>        Echological
## 12915          <NA>         VitaLevaEu
## 12916          <NA>            Carlexa
## 12918          <NA>         RisingFire
## 12401 no_dantrolene   DeterminedYankee
## 12402 no_dantrolene ChoppyChoppyChoppy
## 12403 no_dantrolene            Basheba
## 12620 no_dantrolene        PhyllisFlag
## 12913 no_dantrolene       MizzenColony
```

```r
# Subset count data frame
counts <- counts[, idx]
dim(counts)
```

```
## [1] 14155    10
```

Model Design


```r
design <- model.matrix(~Dx, data=anim)
```

Calculate Log Counts per Million


```r
# Create DGE object using edgeR
dge <- DGEList(counts=counts, group=anim$Dx,genes=annot[rownames(counts),], )

# Apply TMM normalization
dge <- calcNormFactors(dge)

# > Apply voom transformation
vwts <- voomWithQualityWeights(dge, design=design, plot=TRUE)
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.svg)

###  Differential Expression Analysis: Limma


```r
fit <- lmFit(vwts, design)
fit <- treat(fit, lfc=log2(1.1))
rst <- topTable(fit, coef="DxRER", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
```

Proportion of true null hypothesis (pi0)


```r
summary(qvalue(rst$P.Value))
```

```
## 
## Call:
## qvalue(p = rst$P.Value)
## 
## pi0:	0.6520466	
## 
## Cumulative number of significant calls:
## 
##           <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
## p-value       58    410  1603   2482  3397 4608 14155
## q-value        0      0     0    534  1373 2663 14155
## local FDR      0      0    16    282   721 1441  9586
```

Assume the proportion of true null hypothesis is 100%


```r
qval <- qvalue(rst$P.Value, lambda=0)
qv <- qval$qvalues
names(qv) <- rownames(rst)
summary(qval)
```

```
## 
## Call:
## qvalue(p = rst$P.Value, lambda = 0)
## 
## pi0:	1	
## 
## Cumulative number of significant calls:
## 
##           <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
## p-value       58    410  1603   2482  3397 4608 14155
## q-value        0      0     0    129   812 1813 14155
## local FDR      0      0     0    122   434  948  6060
```

Fitted values from model


```r
FV <- fitted(fit)
```

Calculate mean gene expression counts (fitted values) for RER animals


```r
mean.RER <- unlist(lapply(rownames(rst), function(x) 
    mean(FV[x, rownames(anim)[anim$Dx == "RER"]])))
names(mean.RER) <- rownames(rst)
```

Calculate mean gene expression counts (fitted values) for control animals


```r
mean.Cont <- unlist(lapply(rownames(rst), function(x) 
    mean(FV[x, rownames(anim)[!anim$Dx == "RER"]])))
names(mean.Cont ) <- rownames(rst)
```

Merge average fitted values for RER and Controls, results of treat and gene information:


```r
rst <- data.frame(Mean.RER=mean.RER, Mean.Cont=mean.Cont, rst, Q.Value=qv)
dim(rst)
```

```
## [1] 14155    18
```

```r
sum(rst$Q.Value < 0.05)
```

```
## [1] 812
```

### Save DE results to file


```r
write.table(rst, file=paste(getwd(), "Protein_Subset.txt", sep="/"), 
    col.names=TRUE, row.names=TRUE, sep="\t")
```

Save results to R data file


```r
save(dge, vwts, anim, fit, FV, rst, file=paste(getwd(), "Protein_Subset.Rdata", sep="/"))
```

### Run R Script


```r
htmlRunR
DE_RER_Thoroughbred_Protein_Subset.R nodes=1,cpus-per-task=1,time=03:00:00,mem=10G \
+Differential Gene Expression on Subset
```

