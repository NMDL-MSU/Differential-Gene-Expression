---
title: Differential Gene Expression Analysis
author: Deborah Velez-Irizarry
date: Tue Oct 23 12:31:26 EDT 2018
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---
### Description:  
Differential Gene Expression Analysis  
  
***  
**Code:**  
Parent Directory:  

> &nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq  
  
Directory/File:  
 
> &nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_RER_Thoroughbred/DE_RER_Thoroughbred.R  
 
**Input files:**  
Directory/File:  
  
> &nbsp;&nbsp;&nbsp;&nbsp;/HTSeq/htseq_counts_RER_Thoroughbred.txt  
> &nbsp;&nbsp;&nbsp;&nbsp;/RER_Thoroughbred_Animal_Information.txt  
> &nbsp;&nbsp;&nbsp;Cufflinks/MergedGTF/Annotation/annotation.txt  
  
**Output files:**  
  
Directory/File:  
  
> &nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_RER_Thoroughbred/RER_Thoroughbred_DE_results.txt  
> &nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_RER_Thoroughbred/RER_Thoroughbred_DE_results.Rdata  
 
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
dim(annot)
```

```
## [1] 37870     8
```

> Animal Information


```r
anim <- read.table(paste(dir, "RER_Thoroughbred_Animal_Information.txt", sep="/"),
    header=TRUE, row.names=2, sep="\t")
anim$Sex <- rep("F", nrow(anim))
anim$Dx <- as.factor(anim$Dx)
anim$Age <- as.factor(anim$Age)
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
## 12610          <NA>     PromiseOfPeace
## 12611    dantrolene            LaDonia
## 12612 no_dantrolene             LaDama
## 12613          <NA>       SkyHighSugar
## 12614    dantrolene       HangOnSloopy
## 12616    dantrolene     TickleMyChrome
## 12617          <NA>          RunAmyRun
## 12618 no_dantrolene        Bodelicious
## 12619 no_dantrolene      BessiesBullet
## 12620 no_dantrolene        PhyllisFlag
## 12621    dantrolene     AmberlyVillage
## 12622          <NA>            BlueAsh
## 12910          <NA>        Echological
## 12913 no_dantrolene       MizzenColony
## 12915          <NA>         VitaLevaEu
## 12916          <NA>            Carlexa
## 12917 no_dantrolene        BlipSaysBye
## 12918          <NA>         RisingFire
## 12921    dantrolene         MinnieBlip
## 12924          <NA>         Toleration
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
design <- model.matrix(~Dx, data=anim)
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
lrt <- glmLRT(fit,coef=2)
rst <- lrt$table
toprst <- topTags(lrt)[[1]]
toprst
```

```
##               chr    start      end width strand            ID
## XLOC_010827 chr14 20691379 20706550 15172      - RLOC_00010367
## XLOC_020070 chr22 34783006 34796100 13095      + RLOC_00020320
## XLOC_006402 chr11 23552431 23621253 68823      + RLOC_00006752
## XLOC_030626  chr6 24085103 24096031 10929      + RLOC_00030759
## XLOC_019721 chr21 44474416 44474783   368      - RLOC_00019643
## XLOC_033050  chr7 81718056 81763846 45791      + RLOC_00033991
## XLOC_016193  chr2 29609705 29616396  6692      + RLOC_00016528
## XLOC_033575  chr7 35980980 36075052 94073      - RLOC_00032881
## XLOC_024504 chr28 35196997 35215514 18518      - RLOC_00024326
## XLOC_016297  chr2 37210431 37216904  6474      + RLOC_00016727
##                    genes         locus      logFC      logCPM       LR
## XLOC_010827        LSM11 RLOC_00010367 -1.1509367  3.02285256 22.71266
## XLOC_020070        WISP2 RLOC_00020320  1.1180041  2.25075894 22.58962
## XLOC_006402       SRCIN1 RLOC_00006752 -0.7920323  2.20261840 19.44269
## XLOC_030626         ERFE RLOC_00030759  1.5042219  0.48925072 18.00873
## XLOC_019721         <NA> RLOC_00019643 -2.4007706  0.61354545 17.57619
## XLOC_033050        AMPD3 RLOC_00033991  2.1859877  4.04794433 17.46333
## XLOC_016193 LOC111772405 RLOC_00016528  1.0679414 -0.07572105 15.96048
## XLOC_033575         CDON RLOC_00032881 -0.9782888  2.63162310 15.90188
## XLOC_024504        PVALB RLOC_00024326 -2.8595884  3.12555151 15.89577
## XLOC_016297        HSPB7 RLOC_00016727  1.3787859  6.91235500 15.51638
##                   PValue        FDR
## XLOC_010827 1.881248e-06 0.01419497
## XLOC_020070 2.005648e-06 0.01419497
## XLOC_006402 1.036639e-05 0.04891208
## XLOC_030626 2.198944e-05 0.06910076
## XLOC_019721 2.760226e-05 0.06910076
## XLOC_033050 2.929032e-05 0.06910076
## XLOC_016193 6.467883e-05 0.10526299
## XLOC_033575 6.671227e-05 0.10526299
## XLOC_024504 6.692807e-05 0.10526299
## XLOC_016297 8.179333e-05 0.11577845
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
## pi0:	0.9757908	
## 
## Cumulative number of significant calls:
## 
##           <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
## p-value       11     36   193    407   739 1444 14155
## q-value        0      0     0      2     3    6 14155
## local FDR      0      0     0      0     2    3  4720
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
rst[rst$FDR < 0.1,]
```

```
##             Mean.RER Std.Dev.RER Mean.Cont Std.Dev.Cont      logFC
## XLOC_006402  91.2362     20.0613  177.5486      20.7934 -0.7920323
## XLOC_010827 142.8544     31.4113  356.4877      41.7497 -1.1509367
## XLOC_019721  14.3697      3.1597   85.8176      10.0504 -2.4007706
## XLOC_020070 146.3643     32.1831   75.6747       8.8626  1.1180041
## XLOC_030626  43.6888      9.6065   17.2128       2.0159  1.5042219
## XLOC_033050 571.6372    125.6936  141.0048      16.5136  2.1859877
##                logCPM       LR       PValue        FDR   chr    start
## XLOC_006402 2.2026184 19.44269 1.036639e-05 0.04772796 chr11 23552431
## XLOC_010827 3.0228526 22.71266 1.881248e-06 0.01385132 chr14 20691379
## XLOC_019721 0.6135455 17.57619 2.760226e-05 0.06742788 chr21 44474416
## XLOC_020070 2.2507589 22.58962 2.005648e-06 0.01385132 chr22 34783006
## XLOC_030626 0.4892507 18.00873 2.198944e-05 0.06742788  chr6 24085103
## XLOC_033050 4.0479443 17.46333 2.929032e-05 0.06742788  chr7 81718056
##                  end width strand  genes
## XLOC_006402 23621253 68823      + SRCIN1
## XLOC_010827 20706550 15172      -  LSM11
## XLOC_019721 44474783   368      -   <NA>
## XLOC_020070 34796100 13095      +  WISP2
## XLOC_030626 24096031 10929      +   ERFE
## XLOC_033050 81763846 45791      +  AMPD3
```

### Save DE results to file


```r
write.table(rst, file=paste(getwd(), "RER_Thoroughbred_DE_results.txt", sep="/"), 
    col.names=TRUE, row.names=TRUE, sep="\t")
```

Save results to R data file


```r
save(fit, lrt, rst, file=paste(getwd(), "DE_RER_Thoroughbred.Rdata", sep="/"))
```

### Run R Script


```r
htmlRunR
DE_RER_Thoroughbred.R nodes=1,cpus-per-task=1,time=03:00:00,mem=10G \
+Differential Gene Expression Analysis
```

