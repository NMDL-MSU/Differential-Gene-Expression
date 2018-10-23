#' ### Description:  
#' Differential gene expression analysis accounting for Dantrolene treatment.  
#'   
#' ***  
#'   
#' **Code:**  
#' Parent Directory:  
#' 
#' > &nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq  
#'   
#' Directory/File:  
#'  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_RER_Thoroughbred/DE_Dantrolene_RER_Thoroughbred.R  
#'  
#' **Input files:**  
#' Directory/File:  
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;/HTSeq/htseq_counts_RER_Thoroughbred.txt  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/RER_Thoroughbred_Animal_Information.txt  
#' > &nbsp;&nbsp;&nbsp;Cufflinks/MergedGTF/Annotation/annotation.txt  
#'   
#' **Output files:**  
#'   
#' Directory/File:  
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_Dantrolene_RER_Thoroughbred/DE_Dantrolene_RER_Thoroughbred.txt  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_Dantrolene_RER_Thoroughbred/DE_Dantrolene_RER_Thoroughbred.Rdata  
#'   
#' Render R Script  
#' 
#' > &nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_RER_Thoroughbred/DE_Dantrolene_RER_Thoroughbred.qsub  
#'  
#' ***  

#' ### Code  
#' Clear Environment
rm(list=ls())
 
#' Required Packages 
library (limma)
library (edgeR)
library(qvalue)

#' **Session Information**
sessionInfo()

#' **Clear environment**  
rm(list=ls())

#' #### Load required R objects

#' > Gene Counts
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq"
counts <- read.table(paste(dir, "HTSeq", "htseq_counts_RER_Thoroughbred.txt", sep="/"))
colnames(counts) <- unlist(lapply(strsplit(colnames(counts), "X"), function(x) x[2]))
dim(counts)

#' > Annotation
annot <- read.table(paste(dir, "Cufflinks/MergedGTF/Annotation/annotation.txt", sep="/"),
    header=TRUE, row.names=8)
annot <- annot[rownames(counts),]
dim(annot)

#' > Animal Information
anim <- read.table(paste(dir, "RER_Thoroughbred_Animal_Information.txt", sep="/"),
    header=TRUE, row.names=2, sep="\t")
anim$Sex <- rep("F", nrow(anim))
anim$Dx <- as.factor(anim$Dx)
anim$Age <- as.factor(anim$Age)
anim$Tx[is.na(anim$Tx)] <- "no_dantrolene"
anim$Tx <- as.factor(anim$Tx)
anim

#'  ### Summary Function
summSD <- function(x, dig=4) round(c(summary(x),
     Std.Dev.=sd(x)), dig)[c("Min.", "1st Qu.", "Median", "Mean", 
    "Std.Dev.", "3rd Qu.", "Max.")]

#' > **Data Check**: Animal IDs match between count matrix and animal matrix
sum(!rownames(anim) == colnames(counts))

#' ### Prepare Data for DE Analysis: 

#' > Model Design
design <- model.matrix(~Tx + Dx, data=anim)
head(design)

#' > Calculate Log Counts per Million

# Create DGE object using edgeR
dge <- DGEList(counts=counts, group= anim$Dx,genes=annot[rownames(counts),], )

# Apply TMM normalization
dge <- calcNormFactors(dge)

#' Estimate dispersion
dge <- estimateDisp(dge, design) 


#' Perform DE analysis: likelihood ratio test
fit <- glmFit(dge,design)
lrt <- glmLRT(fit,coef=3)
rst <- lrt$table
toprst <- topTags(lrt)
toprst

#' Calculate qvalues (FDR)
qval <- qvalue(rst$PValue)
names(qval$qvalues) <- rownames(lrt$table)
summary(qval)

#' Calculate mean gene expression counts (fitted values) for RER animals
mean.RER <- do.call(rbind, lapply(rownames(rst), function(x) 
    summSD(lrt$fitted.values[x, rownames(anim)[anim$Dx == "RER"]])[c("Mean", "Std.Dev.")]))
rownames(mean.RER) <- rownames(rst)
colnames(mean.RER) <- c("Mean.RER", "Std.Dev.RER")
dim(mean.RER)

#' Calculate mean gene expression counts (fitted values) for control animals
mean.Cont <- do.call(rbind, lapply(rownames(rst), function(x) 
    summSD(lrt$fitted.values[x, rownames(anim)[!anim$Dx == "RER"]])[c("Mean", "Std.Dev.")]))
rownames(mean.Cont ) <- rownames(rst)
colnames(mean.Cont) <- c("Mean.Cont", "Std.Dev.Cont")
dim(mean.Cont)

#' Merge average fitted values for RER and Controls, results of LRT and gene information:
rst <- data.frame(mean.RER, mean.Cont, rst, FDR=qval$qvalue[rownames(rst)], lrt$genes[,c(1:5,7)])
dim(rst)
rst[rst$FDR < 0.12,]

#' ### Save DE results to file
write.table(rst, file=paste(getwd(), "DE_Dantrolene_RER_Thoroughbred.txt", sep="/"), 
    col.names=TRUE, row.names=TRUE, sep="\t")

#' Save results to R data file
save(fit, lrt, rst, file=paste(getwd(), "DE_Dantrolene_RER_Thoroughbred.Rdata", sep="/"))

#' ### Run R Script
#+ eval = FALSE
htmlRunR
DE_Dantrolene_RER_Thoroughbred.R nodes=1,cpus-per-task=1,time=03:00:00,mem=10G \
+Differential Gene Expression Analysis

