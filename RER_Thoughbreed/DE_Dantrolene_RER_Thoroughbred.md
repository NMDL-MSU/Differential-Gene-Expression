---
title: Differential Gene Expression Analysis
author: Deborah Velez-Irizarry
date: Tue Oct 23 12:33:34 EDT 2018
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---
### Description:  
Differential gene expression analysis accounting for Dantrolene treatment.  
  
***  
  
**Code:**  
Parent Directory:  

> &nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq  
  
Directory/File:  
 
> &nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_RER_Thoroughbred/DE_Dantrolene_RER_Thoroughbred.R  
 
**Input files:**  
Directory/File:  
  
> &nbsp;&nbsp;&nbsp;&nbsp;/HTSeq/htseq_counts_RER_Thoroughbred.txt  
> &nbsp;&nbsp;&nbsp;&nbsp;/RER_Thoroughbred_Animal_Information.txt  
> &nbsp;&nbsp;&nbsp;Cufflinks/MergedGTF/Annotation/annotation.txt  
  
**Output files:**  
  
Directory/File:  
  
> &nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_Dantrolene_RER_Thoroughbred/DE_Dantrolene_RER_Thoroughbred.txt  
> &nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_Dantrolene_RER_Thoroughbred/DE_Dantrolene_RER_Thoroughbred.Rdata  
  
Render R Script  

> &nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_RER_Thoroughbred/DE_Dantrolene_RER_Thoroughbred.qsub  
 
***  
### Code  
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

**Clear environment**  


```r
rm(list=ls())
```

#### Load required R objects
> Gene Counts


```r
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq"
counts <- read.table(paste(dir, "HTSeq", "htseq_counts_RER_Thoroughbred.txt", sep="/"))
colnames(counts) <- unlist(lapply(strsplit(colnames(counts), "X"), function(x) x[2]))
dim(counts)
```

```
## [1] 14155    23
```

> Annotation


```r
annot <- read.table(paste(dir, "Cufflinks/MergedGTF/Annotation/annotation.txt", sep="/"),
    header=TRUE, row.names=8)
annot <- annot[rownames(counts),]
dim(annot)
```

```
## [1] 14155     8
```

> Animal Information


```r
anim <- read.table(paste(dir, "RER_Thoroughbred_Animal_Information.txt", sep="/"),
    header=TRUE, row.names=2, sep="\t")
anim$Sex <- rep("F", nrow(anim))
anim$Dx <- as.factor(anim$Dx)
anim$Age <- as.factor(anim$Age)
anim$Tx[is.na(anim$Tx)] <- "no_dantrolene"
anim$Tx <- as.factor(anim$Tx)
anim
```

```
##        X Age Sex Last_Work  AST   CK Biopsy_Changes      Dx Owner_Vet
## 12401  1   7   F       2.0  921  663              2     RER Slaughter
## 12402  2   3   F       3.0  448  521              2     RER     Tores
## 12403  3   3   F       3.0  280  181              3     RER Slaughter
## 12610  4   3   F       1.5  381  262              0     RER Slaughter
## 12611  5   3   F       2.0  654 1086              1     RER Slaughter
## 12612  6   2   F       2.3  503  543              0     RER Slaughter
## 12613  7   4   F       3.0  263  236              0 Control    Fenger
## 12614  8   3   F       3.0  254  200              1     RER    Fenger
## 12616  9   3   F       1.0  840 2445              1     RER Slaughter
## 12617 10   3   F       2.5  264  208              0 Control     Pagan
## 12618 11   3   F       2.0  715  735              1     RER     Pagan
## 12619 12   6   F       0.0 5357 2554              0     RER    Fenger
## 12620 13   4   F       3.0  219  131              2     RER    Fenger
## 12621 14   5   F       4.0  520  185              1     RER    Fenger
## 12622 15   2   F       3.0  541  429              0 Control    Fenger
## 12910 16   4   F       1.5  297  172              0 Control    Fenger
## 12913 17   3   F       0.0  919  251              0     RER  Marshall
## 12915 18   3   F       4.0  267  342              0 Control    Fenger
## 12916 19   2   F       5.0  430  243              0 Control    Fenger
## 12917 20   3   F       5.0  658  703              0     RER    Fenger
## 12918 21   2   F       5.0  537  208              0 Control    Fenger
## 12921 22   4   F       5.5  571  790              0     RER    Fenger
## 12924 23   3   F       2.5  384  331              0 Control  Bob Hunt
##                  Tx               Name
## 12401 no_dantrolene   DeterminedYankee
## 12402 no_dantrolene ChoppyChoppyChoppy
## 12403 no_dantrolene            Basheba
## 12610 no_dantrolene     PromiseOfPeace
## 12611    dantrolene            LaDonia
## 12612 no_dantrolene             LaDama
## 12613 no_dantrolene       SkyHighSugar
## 12614    dantrolene       HangOnSloopy
## 12616    dantrolene     TickleMyChrome
## 12617 no_dantrolene          RunAmyRun
## 12618 no_dantrolene        Bodelicious
## 12619 no_dantrolene      BessiesBullet
## 12620 no_dantrolene        PhyllisFlag
## 12621    dantrolene     AmberlyVillage
## 12622 no_dantrolene            BlueAsh
## 12910 no_dantrolene        Echological
## 12913 no_dantrolene       MizzenColony
## 12915 no_dantrolene         VitaLevaEu
## 12916 no_dantrolene            Carlexa
## 12917 no_dantrolene        BlipSaysBye
## 12918 no_dantrolene         RisingFire
## 12921    dantrolene         MinnieBlip
## 12924 no_dantrolene         Toleration
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
> Model Design


```r
design <- model.matrix(~Tx + Dx, data=anim)
head(design)
```

```
##       (Intercept) Txno_dantrolene DxRER
## 12401           1               1     1
## 12402           1               1     1
## 12403           1               1     1
## 12610           1               1     1
## 12611           1               0     1
## 12612           1               1     1
```

> Calculate Log Counts per Million


```r
# Create DGE object using edgeR
dge <- DGEList(counts=counts, group= anim$Dx,genes=annot[rownames(counts),], )

# Apply TMM normalization
dge <- calcNormFactors(dge)
```

Estimate dispersion


```r
dge <- estimateDisp(dge, design) 
```

Perform DE analysis: likelihood ratio test


```r
fit <- glmFit(dge,design)
lrt <- glmLRT(fit,coef=3)
rst <- lrt$table
toprst <- topTags(lrt)
toprst
```

```
## Coefficient:  DxRER 
##               chr     start       end  width strand            ID
## XLOC_006402 chr11  23552431  23621253  68823      + RLOC_00006752
## XLOC_033050  chr7  81718056  81763846  45791      + RLOC_00033991
## XLOC_010827 chr14  20691379  20706550  15172      - RLOC_00010367
## XLOC_030626  chr6  24085103  24096031  10929      + RLOC_00030759
## XLOC_020070 chr22  34783006  34796100  13095      + RLOC_00020320
## XLOC_016122  chr2  24116281  24123780   7500      + RLOC_00016367
## XLOC_010706 chr14   1879969   1896858  16890      - RLOC_00010127
## XLOC_031957  chr6  74809327  74826691  17365      - RLOC_00031811
## XLOC_007296 chr11  33067522  33199265 131744      - RLOC_00006942
## XLOC_025773  chr3 100918028 100959832  41805      + RLOC_00026384
##                    genes         locus      logFC   logCPM       LR
## XLOC_006402       SRCIN1 RLOC_00006752 -0.8833807 2.202644 19.90934
## XLOC_033050        AMPD3 RLOC_00033991  2.4053247 4.047958 19.58591
## XLOC_010827        LSM11 RLOC_00010367 -1.1607821 3.022872 17.98442
## XLOC_030626         ERFE RLOC_00030759  1.5636767 0.489217 16.94833
## XLOC_020070        WISP2 RLOC_00020320  0.9738592 2.250732 16.01892
## XLOC_016122        FNDC5 RLOC_00016367  1.1760330 4.615894 15.61453
## XLOC_010706        MRNIP RLOC_00010127  3.4170347 3.013458 15.12241
## XLOC_031957      ANKRD52 RLOC_00031811 -0.4971832 5.413894 14.81282
## XLOC_007296        TEX14 RLOC_00006942  2.6728058 1.971052 14.45452
## XLOC_025773 LOC102150003 RLOC_00026384 -0.9660633 3.959065 14.34099
##                   PValue        FDR
## XLOC_006402 8.120248e-06 0.06806841
## XLOC_033050 9.617578e-06 0.06806841
## XLOC_010827 2.227201e-05 0.10508676
## XLOC_030626 3.841110e-05 0.13592728
## XLOC_020070 6.271250e-05 0.17753910
## XLOC_016122 7.765539e-05 0.18320202
## XLOC_010706 1.007602e-04 0.20375156
## XLOC_031957 1.187259e-04 0.21007068
## XLOC_007296 1.435846e-04 0.21587442
## XLOC_025773 1.525075e-04 0.21587442
```

Calculate qvalues (FDR)


```r
qval <- qvalue(rst$PValue)
names(qval$qvalues) <- rownames(lrt$table)
summary(qval)
```

```
## 
## Call:
## qvalue(p = rst$PValue)
## 
## pi0:	0.9047015	
## 
## Cumulative number of significant calls:
## 
##           <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
## p-value        6     41   219    520   972 1776 14155
## q-value        0      0     0      0     0    3 14155
## local FDR      0      0     0      0     0    2 14127
```

Calculate mean gene expression counts (fitted values) for RER animals


```r
mean.RER <- do.call(rbind, lapply(rownames(rst), function(x) 
    summSD(lrt$fitted.values[x, rownames(anim)[anim$Dx == "RER"]])[c("Mean", "Std.Dev.")]))
rownames(mean.RER) <- rownames(rst)
colnames(mean.RER) <- c("Mean.RER", "Std.Dev.RER")
dim(mean.RER)
```

```
## [1] 14155     2
```

Calculate mean gene expression counts (fitted values) for control animals


```r
mean.Cont <- do.call(rbind, lapply(rownames(rst), function(x) 
    summSD(lrt$fitted.values[x, rownames(anim)[!anim$Dx == "RER"]])[c("Mean", "Std.Dev.")]))
rownames(mean.Cont ) <- rownames(rst)
colnames(mean.Cont) <- c("Mean.Cont", "Std.Dev.Cont")
dim(mean.Cont)
```

```
## [1] 14155     2
```

Merge average fitted values for RER and Controls, results of LRT and gene information:


```r
rst <- data.frame(mean.RER, mean.Cont, rst, FDR=qval$qvalue[rownames(rst)], lrt$genes[,c(1:5,7)])
dim(rst)
```

```
## [1] 14155    15
```

```r
rst[rst$FDR < 0.12,]
```

```
##             Mean.RER Std.Dev.RER Mean.Cont Std.Dev.Cont      logFC
## XLOC_006402  91.2481     20.7868  177.5494      20.7935 -0.8833807
## XLOC_010827 142.8525     31.2731  356.4857      41.7494 -1.1607821
## XLOC_033050 571.3017    198.4258  141.0059      16.5137  2.4053247
##               logCPM       LR       PValue        FDR   chr    start
## XLOC_006402 2.202644 19.90934 8.120248e-06 0.06158159 chr11 23552431
## XLOC_010827 3.022872 17.98442 2.227201e-05 0.09507215 chr14 20691379
## XLOC_033050 4.047958 19.58591 9.617578e-06 0.06158159  chr7 81718056
##                  end width strand  genes
## XLOC_006402 23621253 68823      + SRCIN1
## XLOC_010827 20706550 15172      -  LSM11
## XLOC_033050 81763846 45791      +  AMPD3
```

### Save DE results to file


```r
write.table(rst, file=paste(getwd(), "DE_Dantrolene_RER_Thoroughbred.txt", sep="/"), 
    col.names=TRUE, row.names=TRUE, sep="\t")
```

Save results to R data file


```r
save(fit, lrt, rst, file=paste(getwd(), "DE_Dantrolene_RER_Thoroughbred.Rdata", sep="/"))
```

### Run R Script


```r
htmlRunR
DE_Dantrolene_RER_Thoroughbred.R nodes=1,cpus-per-task=1,time=03:00:00,mem=10G \
+Differential Gene Expression Analysis
```

