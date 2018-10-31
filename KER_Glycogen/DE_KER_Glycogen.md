---
title: Differential Gene Expression Analysis
author: Deborah Velez-Irizarry
date: Thu Oct 25 09:08:17 EDT 2018
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---
### Description:  
Differential gene expression analysis for RER glycogen study.  
  
***  
  
**Code:**  
Parent Directory:  

> &nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred  
  
Directory/File:  
 
> &nbsp;&nbsp;&nbsp;&nbsp;/RNA_Seq/DE/DE_KER_Glycogen/DE_KER_Glycogen.R  
 
**Input files:**  
Directory/File:  
  
> &nbsp;&nbsp;&nbsp;&nbsp;/RNA_Seq/HTSeq/htseq_counts_KER_Glycogen.txt  
> &nbsp;&nbsp;&nbsp;&nbsp;Proteomics/Glycogen_Kennedy_20180622/Glycogen_Project_Information.txt  
> &nbsp;&nbsp;&nbsp;/RNA_SeqCufflinks/MergedGTF/Annotation/annotation.txt  
  
**Output files:**  
  
Directory/File:  
  
> &nbsp;&nbsp;&nbsp;&nbsp;/RNA_Seq/DE/DE_KER_Glycogen/results_diff_between_diet_over_time.Rdata  
> &nbsp;&nbsp;&nbsp;&nbsp;/RNA_Seq/DE/DE_KER_Glycogen/results_diet_over_time_pre_ref.Rdata  
> &nbsp;&nbsp;&nbsp;&nbsp;/RNA_Seq/DE/DE_KER_Glycogen/results_diet_over_time_depl_ref.Rdata  
  
Render R Script  

> &nbsp;&nbsp;&nbsp;&nbsp;/RNA_Seq/DE/DE_KER_Glycogen/DE_KER_Glycogen.qsub  
 
***  
### Code  
Clear Environment


```r
rm(list=ls())
```

### Code  
Required Packages 


```r
#library(DESeq2)
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
## BLAS/LAPACK: /opt/software/OpenBLAS/0.2.20-GCC-6.4.0-2.28/lib/libopenblas_haswellp-r0.2.20.so
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
## [1] qvalue_2.12.0 edgeR_3.22.5  limma_3.36.5  knitr_1.20   
## 
## loaded via a namespace (and not attached):
##  [1] locfit_1.5-9.1   Rcpp_0.12.17     lattice_0.20-35  plyr_1.8.4      
##  [5] grid_3.5.1       gtable_0.2.0     magrittr_1.5     evaluate_0.10.1 
##  [9] scales_0.5.0     pillar_1.2.3     ggplot2_2.2.1    rlang_0.2.1     
## [13] stringi_1.2.3    reshape2_1.4.3   lazyeval_0.2.1   splines_3.5.1   
## [17] tools_3.5.1      stringr_1.3.1    munsell_0.5.0    compiler_3.5.1  
## [21] colorspace_1.3-2 tibble_1.4.2
```

### Load required R objects
> Gene Counts


```r
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred"
counts <- read.table(paste(dir, "RNA_Seq/HTSeq/htseq_counts_KER_Glycogen.txt", sep="/"))
dim(counts)
```

```
## [1] 14133    40
```

> Annotation


```r
annot <- read.table(paste(dir, "RNA_Seq/Cufflinks/MergedGTF/Annotation/annotation.txt", sep="/"),
    header=TRUE, row.names=8)
dim(annot)
```

```
## [1] 37870     8
```

> Animal Information


```r
anim <- read.table(paste(dir, 
    "Proteomics/Glycogen_Kennedy_20180622/Glycogen_Project_Information.txt", sep="/"),
    header=TRUE, sep="\t")
```

> Retain information on sequenced animals 


```r
anim <- anim[!is.na(anim$MSMS_Plate),]
anim$TimePoint <- factor(as.character(anim$TimePoint), exclude="Rep48h", 
    levels=c("Pre", "Depl", "Rep24h", "Rep72h"))
anim$Diet <- factor(as.character(anim$Diet), exclude="DietGoldenMax")
rownames(anim) <- paste("G", anim$MSMS_ID, sep="")
head(anim)
```

```
##    Animal Period      Diet DietStarch Horse TimePoint DateTrial MSMS_Plate
## G1  A7934      1       Fat        Low    Pi       Pre 3/27/2012          1
## G2  A7934      1       Fat        Low    Pi      Depl 3/30/2012          1
## G3  A7934      1       Fat        Low    Pi    Rep24h 3/31/2012          1
## G4  A7934      1       Fat        Low    Pi    Rep72h  4/2/2012          1
## G5  A7932      1 SweetFeed       High  King       Pre 3/27/2012          1
## G6  A7932      1 SweetFeed       High  King      Depl 3/30/2012          1
##    MSMS_ID GlycoMN GlycoKA
## G1       1 128.855     130
## G2       2  96.115     138
## G3       3 104.645     176
## G4       4  58.921     112
## G5       5  95.094     108
## G6       6  89.002     142
```

> **Data Check**: Animal IDs match between count matrix and animal matrix  


```r
# Should be zero
sum(!rownames(anim) == colnames(counts))
```

```
## [1] 0
```

### Prepare Data for DE Analysis: 
> Create DGE object using edgeR


```r
dge <- DGEList(counts=counts,genes=annot[rownames(counts),], )
```

> Apply TMM normalization


```r
dge <- calcNormFactors(dge)


### Look for differentially expressed genes between diet at different timepoints
```

> Model Design 


```r
design <- model.matrix(~TimePoint + TimePoint:Horse + TimePoint:Diet, data=anim)
colnames(design)
```

```
##  [1] "(Intercept)"                   "TimePointDepl"                
##  [3] "TimePointRep24h"               "TimePointRep72h"              
##  [5] "TimePointPre:HorseNash"        "TimePointDepl:HorseNash"      
##  [7] "TimePointRep24h:HorseNash"     "TimePointRep72h:HorseNash"    
##  [9] "TimePointPre:HorsePeppe"       "TimePointDepl:HorsePeppe"     
## [11] "TimePointRep24h:HorsePeppe"    "TimePointRep72h:HorsePeppe"   
## [13] "TimePointPre:HorsePi"          "TimePointDepl:HorsePi"        
## [15] "TimePointRep24h:HorsePi"       "TimePointRep72h:HorsePi"      
## [17] "TimePointPre:HorseRalph"       "TimePointDepl:HorseRalph"     
## [19] "TimePointRep24h:HorseRalph"    "TimePointRep72h:HorseRalph"   
## [21] "TimePointPre:DietSweetFeed"    "TimePointDepl:DietSweetFeed"  
## [23] "TimePointRep24h:DietSweetFeed" "TimePointRep72h:DietSweetFeed"
```

> Apply voom transformation


```r
vwts <- voomWithQualityWeights(dge, design=design, plot=TRUE)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.svg)

> Differential Expression Analysis: Limma


```r
fit <- lmFit(vwts, design)
fit <- eBayes(fit)
```

Effect of Diet in PreDepletion


```r
Dpre <- topTable(fit, coef="TimePointPre:DietSweetFeed", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(Dpre$adj.P.Val < 0.05)
```

```
## [1] 1
```

Effect of Diet in Depletion


```r
Ddep <- topTable(fit, coef="TimePointDepl:DietSweetFeed", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(Ddep$adj.P.Val < 0.05)
```

```
## [1] 0
```

Effect of Diet in 24h Repletion


```r
D24hR <- topTable(fit, coef="TimePointRep24h:DietSweetFeed", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(D24hR$adj.P.Val < 0.05)
```

```
## [1] 3
```

Effect of Diet in 72h Repletion


```r
D72hR <- topTable(fit, coef="TimePointRep72h:DietSweetFeed", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(D72hR$adj.P.Val < 0.05)
```

```
## [1] 742
```

Effect of SweetFeed diet over time


```r
DoverT <- topTable(fit, coef=21:24, 
    number=nrow(counts), sort.by="F", confint=TRUE)
sum(DoverT$adj.P.Val < 0.05)
```

```
## [1] 2
```

Merge result to list


```r
# Differentially expressed genes between diet over time
Rst.Timepoint <- list(Pre.Depletion=Dpre, Depletion=Ddep, Repletion.24h=D24hR, Repletion.72h=D72hR)
```

Save results to R data file


```r
save(anim, dge, design, vwts, fit, Rst.Timepoint,
    file=paste(getwd(), "results_diff_between_diet_over_time.Rdata", sep="/"))



### Look for changes in gene expression over time with Pre-Depletion timepoint as reference
```

> Model Design 


```r
design2 <- model.matrix(~Diet + Diet:Horse + Diet:TimePoint, data=anim)
colnames(design2)
```

```
##  [1] "(Intercept)"                   "DietSweetFeed"                
##  [3] "DietFat:HorseNash"             "DietSweetFeed:HorseNash"      
##  [5] "DietFat:HorsePeppe"            "DietSweetFeed:HorsePeppe"     
##  [7] "DietFat:HorsePi"               "DietSweetFeed:HorsePi"        
##  [9] "DietFat:HorseRalph"            "DietSweetFeed:HorseRalph"     
## [11] "DietFat:TimePointDepl"         "DietSweetFeed:TimePointDepl"  
## [13] "DietFat:TimePointRep24h"       "DietSweetFeed:TimePointRep24h"
## [15] "DietFat:TimePointRep72h"       "DietSweetFeed:TimePointRep72h"
```

> Apply voom transformation


```r
vwts2 <- voomWithQualityWeights(dge, design=design2, plot=TRUE)
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22-1.svg)

> Differential Expression Analysis: Limma


```r
fit2 <- lmFit(vwts2, design2)
fit2 <- eBayes(fit2)
```

> Genes differentially expressed between Pre-depletion and Depletion
Depletion - PreDepletion: Animals on Fat diet


```r
Fdepl <- topTable(fit2, coef="DietFat:TimePointDepl", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(Fdepl$adj.P.Val < 0.05)
```

```
## [1] 1197
```

Depletion - PreDepletion: Animals on SweetFeed diet


```r
Sdepl <- topTable(fit2, coef="DietSweetFeed:TimePointDepl", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(Fdepl$adj.P.Val < 0.05)
```

```
## [1] 1197
```

> Genes differentially expressed between Pre-depletion and Repletion at 24h
Repletion 24h - PreDepletion: Animals on Fat diet


```r
F24hRepl <- topTable(fit2, coef="DietFat:TimePointRep24h", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(F24hRepl$adj.P.Val < 0.05)
```

```
## [1] 72
```

Repletion 24h - PreDepletion: Animals on SweetFeed diet


```r
S24hRepl <- topTable(fit2, coef="DietSweetFeed:TimePointRep24h", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(S24hRepl$adj.P.Val < 0.05)
```

```
## [1] 1612
```

> Genes differentially expressed between Pre-depletion and Repletion at 72h
Repletion 72h - PreDepletion: Animals on Fat diet


```r
F72hRepl <- topTable(fit2, coef="DietFat:TimePointRep72h", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(F72hRepl$adj.P.Val < 0.05)
```

```
## [1] 6135
```

Repletion 72h - PreDepletion: Animals on SweetFeed diet


```r
S72hRepl <- topTable(fit2, coef="DietSweetFeed:TimePointRep72h", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(S72hRepl$adj.P.Val < 0.05)
```

```
## [1] 3861
```

> Genes differentially expressed between Pre-depletion over time
Change in gene expression over time compared to reference timepoint 
(Pre-Depletion) for Fat diet


```r
FoverT <- topTable(fit2, coef=c(11, 13, 15),
    number=nrow(counts), sort.by="F", confint=TRUE)
sum(FoverT$adj.P.Val < 0.05)
```

```
## [1] 8232
```

Change in gene expression over time compared to reference timepoint 
(Pre-Depletion) for SweetFeed diet


```r
SoverT <- topTable(fit2, coef=c(12, 14, 16),
    number=nrow(counts), sort.by="F", confint=TRUE)
sum(SoverT$adj.P.Val < 0.05)
```

```
## [1] 4016
```

> Save results
Merge result to list


```r
# Differentially expressed genes between diet over time
Rst.Diet.Pre <- list(Fat.Depletion=Fdepl, Sweet.Depletion=Sdepl, 
    Fat.Repletion.24h=F24hRepl, Sweet.Repletion.24h=S24hRepl,
    Fat.Repletion.72h=F72hRepl, Sweet.Repletion.72h=S72hRepl,
    Fat.over.time=FoverT, Sweet.over.time=SoverT)
```

Save results to R data file


```r
save(anim, dge, design2, vwts2, fit2, Rst.Diet.Pre,
    file=paste(getwd(), "results_diet_over_time_pre_ref.Rdata", sep="/"))



### Look for changes in gene expression over time with Depletion timepoint as reference
```

> Change reference to depletion timepoint


```r
anim$TimePoint <- factor(as.character(anim$TimePoint), 
    levels=c("Depl", "Pre", "Rep24h", "Rep72h"))
```

> Model Design 


```r
design3 <- model.matrix(~Diet + Diet:Horse + Diet:TimePoint, data=anim)
colnames(design3)
```

```
##  [1] "(Intercept)"                   "DietSweetFeed"                
##  [3] "DietFat:HorseNash"             "DietSweetFeed:HorseNash"      
##  [5] "DietFat:HorsePeppe"            "DietSweetFeed:HorsePeppe"     
##  [7] "DietFat:HorsePi"               "DietSweetFeed:HorsePi"        
##  [9] "DietFat:HorseRalph"            "DietSweetFeed:HorseRalph"     
## [11] "DietFat:TimePointPre"          "DietSweetFeed:TimePointPre"   
## [13] "DietFat:TimePointRep24h"       "DietSweetFeed:TimePointRep24h"
## [15] "DietFat:TimePointRep72h"       "DietSweetFeed:TimePointRep72h"
```

> Apply voom transformation


```r
vwts3 <- voomWithQualityWeights(dge, design=design3, plot=TRUE)
```

![plot of chunk unnamed-chunk-36](figure/unnamed-chunk-36-1.svg)

> Differential Expression Analysis: Limma


```r
fit3 <- lmFit(vwts3, design3)
fit3 <- eBayes(fit3)
```

> Genes differentially expressed between Depletion and Repletion at 24h
Repletion 24h - Depletion: Animals on Fat diet


```r
F24hRepl <- topTable(fit3, coef="DietFat:TimePointRep24h", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(F24hRepl$adj.P.Val < 0.05)
```

```
## [1] 551
```

Repletion 24h - Depletion: Animals on SweetFeed diet


```r
S24hRepl <- topTable(fit3, coef="DietSweetFeed:TimePointRep24h", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(S24hRepl$adj.P.Val < 0.05)
```

```
## [1] 32
```

> Genes differentially expressed between Depletion and Repletion at 72h
Repletion 72h - Depletion: Animals on Fat diet


```r
F72hRepl <- topTable(fit3, coef="DietFat:TimePointRep72h", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(F72hRepl$adj.P.Val < 0.05)
```

```
## [1] 4191
```

Repletion 72h - Depletion: Animals on SweetFeed diet


```r
S72hRepl <- topTable(fit3, coef="DietSweetFeed:TimePointRep72h", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(S72hRepl$adj.P.Val < 0.05)
```

```
## [1] 2135
```

> Genes differentially expressed between Pre-depletion over time
Change in gene expression over time compared to reference timepoint 
(Depletion) for Fat diet


```r
FoverT <- topTable(fit3, coef=c(11, 13, 15),
    number=nrow(counts), sort.by="F", confint=TRUE)
sum(FoverT$adj.P.Val < 0.05)
```

```
## [1] 7819
```

Change in gene expression over time compared to reference timepoint 
(Depletion) for SweetFeed diet


```r
SoverT <- topTable(fit3, coef=c(12, 14, 16),
    number=nrow(counts), sort.by="F", confint=TRUE)
sum(SoverT$adj.P.Val < 0.05)
```

```
## [1] 4692
```

> Save results
Merge result to list


```r
# Differentially expressed genes between diet over time
Rst.Diet.Depl <- list(Fat.Repletion.24h=F24hRepl, Sweet.Repletion.24h=S24hRepl,
    Fat.Repletion.72h=F72hRepl, Sweet.Repletion.72h=S72hRepl,
    Fat.over.time=FoverT, Sweet.over.time=SoverT)
```

Save results to R data file


```r
save(anim, dge, design3, vwts3, fit3, Rst.Diet.Depl,
    file=paste(getwd(), "results_diet_over_time_depl_ref.Rdata", sep="/"))
```

### Run R Script


```r
htmlRunR
DE_KER_Glycogen.R nodes=1,cpus-per-task=1,time=03:00:00,mem=10G \
+Differential Gene Expression Analysis
```

