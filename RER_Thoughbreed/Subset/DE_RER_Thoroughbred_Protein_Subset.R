#' ### Description:  
#' Differential gene expression analysis for subset of animals used in proteomic study.  
#'   
#' ***  
#' **Code:**  
#' Parent Directory:  
#' 
#' > &nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq  
#'   
#' Directory/File:  
#'  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_KER_Glycogen/DE_RER_Thoroughbred_Protein_Subset/DE_RER_Thoroughbred_Protein_Subset.R  
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
#' Directory:  
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_KER_Glycogen/DE_RER_Thoroughbred_Protein_Subset  
#'   
#' Files:  
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;Protein_Subset.txt  
#' > &nbsp;&nbsp;&nbsp;&nbsp;Protein_Subset.Rdata  
#'  
#' Render R Script  
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;DE_RER_Thoroughbred_Protein_Subset.qsub  
#'  
#' ***  

#' ### R Environment  
#' Clear Environment
rm(list=ls())
 
#' Required Packages 
library (limma)
library (edgeR)
library(qvalue)

#' **Session Information**
sessionInfo()


#' ### Load required R objects

#' Gene Counts
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq"
counts <- read.table(paste(dir, "HTSeq", "htseq_counts_RER_Thoroughbred.txt", sep="/"))
colnames(counts) <- unlist(lapply(strsplit(colnames(counts), "X"), function(x) x[2]))
dim(counts)

#' Annotation
annot <- read.table(paste(dir, "Cufflinks/MergedGTF/Annotation/annotation.txt", sep="/"),
    header=TRUE, row.names=8)
dim(annot)

#' Animal Information
anim <- read.table(paste(dir, "RER_Thoroughbred_Animal_Information.txt", sep="/"),
    header=TRUE, row.names=2, sep="\t")[,-1]
anim$Sex <- rep("F", nrow(anim))
anim$Dx <- as.factor(anim$Dx)
anim$Age <- as.factor(anim$Age)

#'  ### Summary Function
summSD <- function(x, dig=4) round(c(summary(x),
     Std.Dev.=sd(x)), dig)[c("Min.", "1st Qu.", "Median", "Mean", 
    "Std.Dev.", "3rd Qu.", "Max.")]

#' > **Data Check**: Animal IDs match between count matrix and animal matrix
sum(!rownames(anim) == colnames(counts))


#' ### Prepare Data for DE Analysis:  

#' Retain gene expression on subset of animals used in proteomic analysis
# Animals to retain
idx <- c("12613", "12910", "12915", "12916", "12918", "12401", "12402", "12403", "12620", "12913")

# Subset animal data frame
anim <- anim[idx,]
anim

# Subset count data frame
counts <- counts[, idx]
dim(counts)

#' Model Design
design <- model.matrix(~Dx, data=anim)

#' Calculate Log Counts per Million

# Create DGE object using edgeR
dge <- DGEList(counts=counts, group=anim$Dx,genes=annot[rownames(counts),], )

# Apply TMM normalization
dge <- calcNormFactors(dge)

# > Apply voom transformation
vwts <- voomWithQualityWeights(dge, design=design, plot=TRUE)


#' ###  Differential Expression Analysis: Limma
fit <- lmFit(vwts, design)
fit <- treat(fit, lfc=log2(1.1))
rst <- topTable(fit, coef="DxRER", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)

#' Proportion of true null hypothesis (pi0)
summary(qvalue(rst$P.Value))

#' Assume the proportion of true null hypothesis is 100%
qval <- qvalue(rst$P.Value, lambda=0)
qv <- qval$qvalues
names(qv) <- rownames(rst)
summary(qval)

#' Fitted values from model
FV <- fitted(fit)

#' Calculate mean gene expression counts (fitted values) for RER animals
mean.RER <- unlist(lapply(rownames(rst), function(x) 
    mean(FV[x, rownames(anim)[anim$Dx == "RER"]])))
names(mean.RER) <- rownames(rst)

#' Calculate mean gene expression counts (fitted values) for control animals
mean.Cont <- unlist(lapply(rownames(rst), function(x) 
    mean(FV[x, rownames(anim)[!anim$Dx == "RER"]])))
names(mean.Cont ) <- rownames(rst)


#' Merge average fitted values for RER and Controls, results of treat and gene information:
rst <- data.frame(Mean.RER=mean.RER, Mean.Cont=mean.Cont, rst, Q.Value=qv)
dim(rst)
sum(rst$Q.Value < 0.05)

#' ### Save DE results to file
write.table(rst, file=paste(getwd(), "Protein_Subset.txt", sep="/"), 
    col.names=TRUE, row.names=TRUE, sep="\t")

#' Save results to R data file
save(dge, vwts, anim, fit, FV, rst, file=paste(getwd(), "Protein_Subset.Rdata", sep="/"))


#' ### Run R Script
#+ eval = FALSE
htmlRunR
DE_RER_Thoroughbred_Protein_Subset.R nodes=1,cpus-per-task=1,time=03:00:00,mem=10G \
+Differential Gene Expression on Subset

