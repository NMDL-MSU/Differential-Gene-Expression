### Description:  
#' Performe gene set enrichment analysis for RER Thoroughbred project using the subset of animals with proteomics  
#' data. Look at enrichment of genes differentially expressed in horses presenting recurrent exertional  
#' rhabdomyolysis (RER) compared to controls.  
#'   
#' ***  
#' **Code:**  
#' Parent Directory:  
#' 
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq  
#'   
#' Directory/File:  
#'  
#' &nbsp;&nbsp;&nbsp;&nbsp;/Enrichment/Enrichment_RER_Thoroughbred/Enrichment_RER_Thoroughbred.R  
#'  
#' **Input files:**  
#' Directory/File:  
#'   
#' >&nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_RER_Thoroughbred/DE_RER_Thoroughbred_Protein_Subset/Protein_Subset.Rdata  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/Enrichment/Enrichment_KER_Glycogen_Depl/enrichment_function.Rdata
#'   
#' **Output files:**  
#'   
#' Directory:  
#'   
#' >&nbsp;&nbsp;&nbsp;&nbsp;/Enrichment/Enrichment_RER_Thoroughbred/  
#' 
#' Files:
#' 
#' >&nbsp;&nbsp;&nbsp;&nbsp;/Merged/*.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/UpRegulated/*.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/DownRegulated/*.txt  
#'  
#' ***  

#' ### R Environment    

#' Load required libraries
library(edgeR)
library(limma)
library(DOSE)
library(GO.db)
library(GSEABase)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GenomicFeatures)

#' Clear Environment
rm(list=ls())

#' Session Information
sessionInfo()


#' ### Load Data for Pathway Analysis

#' Load DE results
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq"
load(paste(dir, 
    "/DE/DE_RER_Thoroughbred/DE_RER_Thoroughbred_Protein_Subset/Protein_Subset.Rdata", 
    sep=""))

#' Load enrichment function
load(paste(dir, "Enrichment/Enrichment_KER_Glycogen_Depl", 
    "enrichment_function.Rdata", sep="/"))

#' ### Animal Groups
group <- anim$Dx
table(group)


#' ### Differential Expression Analysis Results

#' Reduce results list to contain only differentially expressed genes
Rrst <- rst[rst$Q.Val < 0.05, c(1:5, 7, 9, 11:17)]
nrow(Rrst)

#' Write DE results to file
write.table(Rrst, file=,"DE_RER_Thoroughbred_Protein_Subset.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")

#' Significant gene names
sigG <- as.character(Rrst$genes)
names(sigG) <- rownames(Rrst)

#' Seperate genes without a name
# Number of unknown gene transcripts
na.sigG <- sigG[is.na(sigG)]
length(na.sigG)

# Number of known genes
sigG <- sigG[!is.na(sigG)]
length(sigG)


#' ### Background for Enrichment Analysis

#' Background Gene List
annot <- dge$genes

#' Seperate genes without a name
na.Bkg <- annot[is.na(annot$genes),]
nrow(na.Bkg)

Bkg <- annot[!is.na(annot$genes),]
nrow(Bkg)

#' Obtain human EntrezIDs for gene enrichment analysis
# DE gene list
sigG.entrez <- bitr(unique(sigG), fromType="SYMBOL", 
    toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
nrow(sigG.entrez)

# Background gene list
bkg.entrez <- bitr(as.character(unique(Bkg$genes)), fromType="SYMBOL", 
    toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
nrow(bkg.entrez)


#' ### Enrichment Analysis  

#' Global gene enrichment analysis: Merge all DE genes
enrich.rst <- enrich(lst=sigG.entrez, bkg=bkg.entrez)
unlist(lapply(enrich.rst, nrow))

#' Summary of significant GO terms
lapply(enrich.rst, function(x) head(x[,2]))

#' Save enrichment analysis results to file
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


#' ### Enrichment of Up-regulated Genes  

#' Gene set enrichment for upregulated genes
up <- Rrst[Rrst$logFC > 0,]
up <- lapply(as.character(up$genes), 
    function(x) sigG.entrez[sigG.entrez$SYMBOL %in% x,])
up <- do.call(rbind, up[unlist(lapply(up, nrow)) > 0])
nrow(up)


#' Gene enrichment analysis for upregulated genes:
enrich.upreg <- enrich(lst=up, bkg=bkg.entrez)
unlist(lapply(enrich.upreg, nrow))

#' Summary of significant GO terms
lapply(enrich.upreg, function(x) head(x[,2]))

#' Save results
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



#' ### Enrichment of Down-regulated Genes  

#' Gene set enrichment for down-regulated genes
down <- Rrst[Rrst$logFC < 0,]
down <- lapply(as.character(down$genes), 
    function(x) sigG.entrez[sigG.entrez$SYMBOL %in% x,])
down <- do.call(rbind, down[unlist(lapply(down, nrow)) > 0])
nrow(down)


#' Gene enrichment analysis for down-regulated genes:
enrich.downreg <- enrich(lst=down, bkg=bkg.entrez)
unlist(lapply(enrich.downreg, nrow))

#' Summary of significant GO terms
lapply(enrich.downreg, function(x) head(x[,2]))

#' Save results
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



#' ### Run R Script
#+ eval = FALSE
htmlRunR
Enrichment_RER_Thoroughbred.R nodes=1,cpus-per-task=1,time=03:00:00,mem=50G \
+Gene Set Enrichment RER Thoroughbred Project

