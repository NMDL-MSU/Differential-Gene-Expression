#' ### Description:  
#' Differential gene expression analysis for RER glycogen study.  
#'   
#' ***  
#'   
#' **Code:**  
#' Parent Directory:  
#' 
#' > &nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred  
#'   
#' Directory/File:  
#'  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/RNA_Seq/DE/DE_KER_Glycogen/DE_KER_Glycogen.R  
#'  
#' **Input files:**  
#' Directory/File:  
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;/RNA_Seq/HTSeq/htseq_counts_KER_Glycogen.txt  
#' > &nbsp;&nbsp;&nbsp;&nbsp;Proteomics/Glycogen_Kennedy_20180622/Glycogen_Project_Information.txt  
#' > &nbsp;&nbsp;&nbsp;/RNA_SeqCufflinks/MergedGTF/Annotation/annotation.txt  
#'   
#' **Output files:**  
#'   
#' Directory/File:  
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;/RNA_Seq/DE/DE_KER_Glycogen/results_diff_between_diet_over_time.Rdata  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/RNA_Seq/DE/DE_KER_Glycogen/results_diet_over_time_pre_ref.Rdata  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/RNA_Seq/DE/DE_KER_Glycogen/results_diet_over_time_depl_ref.Rdata  
#'   
#' Render R Script  
#' 
#' > &nbsp;&nbsp;&nbsp;&nbsp;/RNA_Seq/DE/DE_KER_Glycogen/DE_KER_Glycogen.qsub  
#'  
#' ***  

#' ### Code  
#' Clear Environment
rm(list=ls())


#' ### Code  
#' Required Packages 
#library(DESeq2)
library (limma)
library (edgeR)
library(qvalue)

#' **Session Information**
sessionInfo()


#' ### Load required R objects
#' > Gene Counts
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred"
counts <- read.table(paste(dir, "RNA_Seq/HTSeq/htseq_counts_KER_Glycogen.txt", sep="/"))
dim(counts)

#' > Annotation
annot <- read.table(paste(dir, "RNA_Seq/Cufflinks/MergedGTF/Annotation/annotation.txt", sep="/"),
    header=TRUE, row.names=8)
dim(annot)

#' > Animal Information
anim <- read.table(paste(dir, 
    "Proteomics/Glycogen_Kennedy_20180622/Glycogen_Project_Information.txt", sep="/"),
    header=TRUE, sep="\t")

#' > Retain information on sequenced animals 
anim <- anim[!is.na(anim$MSMS_Plate),]
anim$TimePoint <- factor(as.character(anim$TimePoint), exclude="Rep48h", 
    levels=c("Pre", "Depl", "Rep24h", "Rep72h"))
anim$Diet <- factor(as.character(anim$Diet), exclude="DietGoldenMax")
rownames(anim) <- paste("G", anim$MSMS_ID, sep="")
head(anim)

#' > **Data Check**: Animal IDs match between count matrix and animal matrix  
# Should be zero
sum(!rownames(anim) == colnames(counts))

#' ### Prepare Data for DE Analysis: 

#' > Create DGE object using edgeR
dge <- DGEList(counts=counts,genes=annot[rownames(counts),], )

#' > Apply TMM normalization
dge <- calcNormFactors(dge)


### Look for differentially expressed genes between diet at different timepoints

#' > Model Design 
design <- model.matrix(~TimePoint + TimePoint:Horse + TimePoint:Diet, data=anim)
colnames(design)

#' > Apply voom transformation
vwts <- voomWithQualityWeights(dge, design=design, plot=TRUE)

#' > Differential Expression Analysis: Limma
fit <- lmFit(vwts, design)
fit <- eBayes(fit)

#' Effect of Diet in PreDepletion
Dpre <- topTable(fit, coef="TimePointPre:DietSweetFeed", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(Dpre$adj.P.Val < 0.05)

#' Effect of Diet in Depletion
Ddep <- topTable(fit, coef="TimePointDepl:DietSweetFeed", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(Ddep$adj.P.Val < 0.05)

#' Effect of Diet in 24h Repletion
D24hR <- topTable(fit, coef="TimePointRep24h:DietSweetFeed", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(D24hR$adj.P.Val < 0.05)

#' Effect of Diet in 72h Repletion
D72hR <- topTable(fit, coef="TimePointRep72h:DietSweetFeed", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(D72hR$adj.P.Val < 0.05)

#' Effect of SweetFeed diet over time
DoverT <- topTable(fit, coef=21:24, 
    number=nrow(counts), sort.by="F", confint=TRUE)
sum(DoverT$adj.P.Val < 0.05)

#' Merge result to list
# Differentially expressed genes between diet over time
Rst.Timepoint <- list(Pre.Depletion=Dpre, Depletion=Ddep, Repletion.24h=D24hR, Repletion.72h=D72hR)

#' Save results to R data file
save(anim, dge, design, vwts, fit, Rst.Timepoint,
    file=paste(getwd(), "results_diff_between_diet_over_time.Rdata", sep="/"))



### Look for changes in gene expression over time with Pre-Depletion timepoint as reference

#' > Model Design 
design2 <- model.matrix(~Diet + Diet:Horse + Diet:TimePoint, data=anim)
colnames(design2)

#' > Apply voom transformation
vwts2 <- voomWithQualityWeights(dge, design=design2, plot=TRUE)

#' > Differential Expression Analysis: Limma
fit2 <- lmFit(vwts2, design2)
fit2 <- eBayes(fit2)

#' > Genes differentially expressed between Pre-depletion and Depletion

#' Depletion - PreDepletion: Animals on Fat diet
Fdepl <- topTable(fit2, coef="DietFat:TimePointDepl", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(Fdepl$adj.P.Val < 0.05)

#' Depletion - PreDepletion: Animals on SweetFeed diet
Sdepl <- topTable(fit2, coef="DietSweetFeed:TimePointDepl", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(Fdepl$adj.P.Val < 0.05)


#' > Genes differentially expressed between Pre-depletion and Repletion at 24h

#' Repletion 24h - PreDepletion: Animals on Fat diet
F24hRepl <- topTable(fit2, coef="DietFat:TimePointRep24h", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(F24hRepl$adj.P.Val < 0.05)

#' Repletion 24h - PreDepletion: Animals on SweetFeed diet
S24hRepl <- topTable(fit2, coef="DietSweetFeed:TimePointRep24h", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(S24hRepl$adj.P.Val < 0.05)


#' > Genes differentially expressed between Pre-depletion and Repletion at 72h

#' Repletion 72h - PreDepletion: Animals on Fat diet
F72hRepl <- topTable(fit2, coef="DietFat:TimePointRep72h", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(F72hRepl$adj.P.Val < 0.05)

#' Repletion 72h - PreDepletion: Animals on SweetFeed diet
S72hRepl <- topTable(fit2, coef="DietSweetFeed:TimePointRep72h", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(S72hRepl$adj.P.Val < 0.05)


#' > Genes differentially expressed between Pre-depletion over time

#' Change in gene expression over time compared to reference timepoint 
#' (Pre-Depletion) for Fat diet
FoverT <- topTable(fit2, coef=c(11, 13, 15),
    number=nrow(counts), sort.by="F", confint=TRUE)
sum(FoverT$adj.P.Val < 0.05)

#' Change in gene expression over time compared to reference timepoint 
#' (Pre-Depletion) for SweetFeed diet
SoverT <- topTable(fit2, coef=c(12, 14, 16),
    number=nrow(counts), sort.by="F", confint=TRUE)
sum(SoverT$adj.P.Val < 0.05)


#' > Save results

#' Merge result to list
# Differentially expressed genes between diet over time
Rst.Diet.Pre <- list(Fat.Depletion=Fdepl, Sweet.Depletion=Sdepl, 
    Fat.Repletion.24h=F24hRepl, Sweet.Repletion.24h=S24hRepl,
    Fat.Repletion.72h=F72hRepl, Sweet.Repletion.72h=S72hRepl,
    Fat.over.time=FoverT, Sweet.over.time=SoverT)

#' Save results to R data file
save(anim, dge, design2, vwts2, fit2, Rst.Diet.Pre,
    file=paste(getwd(), "results_diet_over_time_pre_ref.Rdata", sep="/"))



### Look for changes in gene expression over time with Depletion timepoint as reference

#' > Change reference to depletion timepoint
anim$TimePoint <- factor(as.character(anim$TimePoint), 
    levels=c("Depl", "Pre", "Rep24h", "Rep72h"))

#' > Model Design 
design3 <- model.matrix(~Diet + Diet:Horse + Diet:TimePoint, data=anim)
colnames(design3)

#' > Apply voom transformation
vwts3 <- voomWithQualityWeights(dge, design=design3, plot=TRUE)

#' > Differential Expression Analysis: Limma
fit3 <- lmFit(vwts3, design3)
fit3 <- eBayes(fit3)

#' > Genes differentially expressed between Depletion and Repletion at 24h

#' Repletion 24h - Depletion: Animals on Fat diet
F24hRepl <- topTable(fit3, coef="DietFat:TimePointRep24h", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(F24hRepl$adj.P.Val < 0.05)

#' Repletion 24h - Depletion: Animals on SweetFeed diet
S24hRepl <- topTable(fit3, coef="DietSweetFeed:TimePointRep24h", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(S24hRepl$adj.P.Val < 0.05)


#' > Genes differentially expressed between Depletion and Repletion at 72h

#' Repletion 72h - Depletion: Animals on Fat diet
F72hRepl <- topTable(fit3, coef="DietFat:TimePointRep72h", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(F72hRepl$adj.P.Val < 0.05)

#' Repletion 72h - Depletion: Animals on SweetFeed diet
S72hRepl <- topTable(fit3, coef="DietSweetFeed:TimePointRep72h", 
    number=nrow(counts), sort.by="logFC", confint=TRUE)
sum(S72hRepl$adj.P.Val < 0.05)


#' > Genes differentially expressed between Pre-depletion over time

#' Change in gene expression over time compared to reference timepoint 
#' (Depletion) for Fat diet
FoverT <- topTable(fit3, coef=c(11, 13, 15),
    number=nrow(counts), sort.by="F", confint=TRUE)
sum(FoverT$adj.P.Val < 0.05)

#' Change in gene expression over time compared to reference timepoint 
#' (Depletion) for SweetFeed diet
SoverT <- topTable(fit3, coef=c(12, 14, 16),
    number=nrow(counts), sort.by="F", confint=TRUE)
sum(SoverT$adj.P.Val < 0.05)


#' > Save results

#' Merge result to list
# Differentially expressed genes between diet over time
Rst.Diet.Depl <- list(Fat.Repletion.24h=F24hRepl, Sweet.Repletion.24h=S24hRepl,
    Fat.Repletion.72h=F72hRepl, Sweet.Repletion.72h=S72hRepl,
    Fat.over.time=FoverT, Sweet.over.time=SoverT)

#' Save results to R data file
save(anim, dge, design3, vwts3, fit3, Rst.Diet.Depl,
    file=paste(getwd(), "results_diet_over_time_depl_ref.Rdata", sep="/"))


#' ### Run R Script
#+ eval = FALSE
htmlRunR
DE_KER_Glycogen.R nodes=1,cpus-per-task=1,time=03:00:00,mem=10G \
+Differential Gene Expression Analysis

