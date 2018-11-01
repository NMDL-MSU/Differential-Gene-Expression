#' ### Description:  
#' Performe gene set enrichment analysis for KER Glycogen project. Look at enrichment of  
#' genes differentially expressed for each diet over time with depletion timepoint as  
#' reference.   
#'   
#' ***  
#' **Code:**  
#' Parent Directory:  
#' 
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq  
#'   
#' Directory/File:  
#'  
#' &nbsp;&nbsp;&nbsp;&nbsp;/Enrichment/Enrichment_KER_Glycogen_Depl/Enrichment_KER_Glycogen_Depl.R  
#'  
#' **Input files:**  
#' Directory/File:  
#'   
#' >&nbsp;&nbsp;&nbsp;&nbsp;/DE/DE_KER_Glycogen/results_diet_over_time_depl_ref.Rdata  
#'   
#' **Output files:**  
#'   
#' Directory:  
#'   
#' >&nbsp;&nbsp;&nbsp;&nbsp;/Enrichment/Enrichment_KER_Glycogen_Depl/  
#' 
#' Files:
#' 
#' >&nbsp;&nbsp;&nbsp;&nbsp;/enrichment_function.Rdata
#' >&nbsp;&nbsp;&nbsp;&nbsp;/Enrichment_KER_Glycogen_Depl.Rdata  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/Fat.Repletion.24h/*.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/Sweet.Repletion.24h/*.txt
#' >&nbsp;&nbsp;&nbsp;&nbsp;/Fat.Repletion.72h/*.txt
#' >&nbsp;&nbsp;&nbsp;&nbsp;/Sweet.Repletion.72h/*.txt
#' >&nbsp;&nbsp;&nbsp;&nbsp;/Merged/*.txt
#' >&nbsp;&nbsp;&nbsp;&nbsp;/upregulated/*.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/downregulated/*.txt  
#'  
#' ***  

#' ### Code  

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

#' Load DE results
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq/DE/DE_KER_Glycogen"
load(paste(dir, "results_diet_over_time_depl_ref.Rdata", sep="/"))


#' ### Animal Groups
group <- paste(anim$Diet, anim$TimePoint, sep=":")
table(group)


#' ### Differential Expression Analysis Results

#' Reduce results list to contain only differentially expressed genes
Rrst <- lapply(Rst.Diet.Depl, function(x) x[x$adj.P.Val < 0.05,])
lapply(Rrst, nrow)

#' Write DE results to file
system("mkdir DE.Results")
lapply(names(Rrst), function(x) 
    write.table(Rrst[[x]], file=paste(getwd(), "/DE.Results/", gsub("[.]", "_", x), "_depl_ref.txt", sep=""),
    col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t"))

#' Merge list of genes
sigG <- do.call(rbind, lapply(Rrst, function(x) data.frame(Xloc=rownames(x), x[,c(1:5,7)])))
rownames(sigG) <- NULL
nrow(sigG)

#' Seperate genes without a name
na.sigG <- sigG[is.na(sigG$genes),]
nrow(na.sigG)

sigG <- sigG[!is.na(sigG$genes),]
nrow(sigG)


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
sigG.entrez <- bitr(as.character(unique(sigG$genes)), fromType="SYMBOL", 
    toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
nrow(sigG.entrez)

# Background gene list
bkg.entrez <- bitr(as.character(unique(Bkg$genes)), fromType="SYMBOL", 
    toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
nrow(bkg.entrez)


#' ### Functions  
#' **Function to performe enrichment analysis from a gene list**  

#' > **INPUT**: use output from bitr   
#' >&nbsp;&nbsp;&nbsp;&nbsp;`lst` = data frame with gene symbols and gene IDs (colnames must be "SYMBOL" and "ENTREZID")  
#' >&nbsp;&nbsp;&nbsp;&nbsp;`bkg` = background gene list for enrichment analysis with gene symbol and gene IDs   
#' >&nbsp;&nbsp;&nbsp;&nbsp;;&nbsp(colnames must be "SYMBOL" and "ENTREZID")  
enrich <- function(lst, bkg){
    lstE <- lst$ENTREZID
    bkgE <- bkg$ENTREZID
    # GO: biological processes (BP)
    go.bp <- enrichGO(gene = lstE, universe = bkgE, 
        OrgDb=org.Hs.eg.db, ont = "BP", pAdjustMethod = "fdr", qvalueCutoff = 0.05,
        readable = TRUE)
    # GO: cellular component (CC)
    go.cc <- enrichGO(gene = lstE, universe = bkgE, 
        OrgDb=org.Hs.eg.db, ont = "CC", pAdjustMethod = "fdr", qvalueCutoff = 0.05,
        readable = TRUE)
    # GO: molecular function (MF)
    go.mf <- enrichGO(gene = lstE, universe = bkgE, 
        OrgDb=org.Hs.eg.db, ont = "MF", pAdjustMethod = "fdr", qvalueCutoff = 0.05,
        readable = TRUE)
    # Kegg Enrichment
    kegg <- enrichKEGG(gene = lstE, universe = bkgE, 
        organism="hsa", pAdjustMethod = "fdr", qvalueCutoff = 0.05)
    # Merge and return results
    rst <- list(Genes=lst, GO.BP=go.bp, GO.CC=go.cc, GO.MF=go.mf, Kegg=kegg)
    return(rst)
}

#' Save function to file
save(enrich, file=paste(getwd(), "enrichment_function.Rdata", sep="/"))

#' ### Global Enrichment Analysis  

#' Global gene enrichment analysis: Merge all DE genes
enrich.rst.fit3 <- enrich(lst=sigG.entrez, bkg=bkg.entrez)
unlist(lapply(enrich.rst.fit3, nrow))

#' Summary of significant GO terms
lapply(enrich.rst.fit3, function(x) head(x[,2]))

#' Save enrichment analysis results to file
# Enrichemnt analysis results
data_merged <- lapply(enrich.rst.fit3, function(x) data.frame(x))

# Add gene names to Kegg results (currently has the geneIDs)
kegg.gene <- strsplit(data_merged$Kegg$geneID, "/")[[1]]
data_merged$Kegg$geneID <- paste(data_merged$Gene$SYMBOL[data_merged$Gene$ENTREZID %in% kegg.gene], collapse="/")

# Save results to file
system("mkdir Merged")
lapply(names(data_merged), function(x) 
    write.table(data_merged[[x]], file=paste(getwd(), "/Merged/", x, "_merged.txt", sep=""),
    col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t"))


#' ### Enrichment Analysis Per Contrast  

#' Obtain human EntrezIDs per contrast:
enrichAll.rst.list <- lapply(Rrst, function(x) 
    do.call(rbind, lapply(as.character(x$genes), 
    function(y) sigG.entrez[sigG.entrez$SYMBOL %in% y,])))
names(enrichAll.rst.list) <- names(Rrst)
unlist(lapply(enrichAll.rst.list, nrow))

#' Gene enrichment analysis per contrast:
enrichCont.rst.fit3 <- lapply(enrichAll.rst.list, 
    function(x) enrich(lst=x, bkg=bkg.entrez))
lapply(enrichCont.rst.fit3, function(x)
    unlist(lapply(x, function(y) nrow(y))))

#' Gene Ontology: Biological Processes (first 10)
lapply(enrichCont.rst.fit3, function(x) 
    head(x[["GO.BP"]][,2]))

#' Gene Ontology: Cellular Processes (first 10)
lapply(enrichCont.rst.fit3, function(x) 
    head(x[["GO.CC"]][,2]))

#' Gene Ontology: Molecular Processes (first 10)
lapply(enrichCont.rst.fit3, function(x) 
    head(x[["GO.MP"]][,2]))

#' Kegg Pathways:
lapply(enrichCont.rst.fit3, function(x) 
    head(x[["Kegg"]][,2]))

#' Save results
data_contrast <- lapply(names(enrichCont.rst.fit3), function(x) 
    lapply(names(enrichCont.rst.fit3[[x]]), function(y) 
        data.frame(enrichCont.rst.fit3[[x]][[y]])))
names(data_contrast) <- names(enrichCont.rst.fit3)
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


#' ### Enrichment of Up-regulated Genes  

#' Gene set enrichment for upregulated genes
up <- lapply(Rrst, function(x) x[x$logFC > 0,])
up <- lapply(up, function(x) 
    do.call(rbind, lapply(as.character(x$genes), 
    function(y) sigG.entrez[sigG.entrez$SYMBOL %in% y,])))
up <- up[!unlist(lapply(up, is.null))]
unlist(lapply(up, nrow))


#' Gene enrichment analysis per contrast for upregulated genes:
enrichCont.up <- lapply(up, 
    function(x) enrich(lst=x, bkg=bkg.entrez))
lapply(enrichCont.up, function(x)
    unlist(lapply(x, function(y) nrow(y))))

#' Save results
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




#' ### Enrichment of Down-regulated Genes  

#' Gene set enrichment for downregulated genes
down <- lapply(Rrst, function(x) x[x$logFC < 0,])
down <- lapply(down, function(x) 
    do.call(rbind, lapply(as.character(x$genes), 
    function(y) sigG.entrez[sigG.entrez$SYMBOL %in% y,])))
down <- down[!unlist(lapply(down, is.null))]
unlist(lapply(down, nrow))


#' Gene enrichment analysis per contrast for downregulated genes:
enrichCont.down <- lapply(down, 
    function(x) enrich(lst=x, bkg=bkg.entrez))
lapply(enrichCont.down, function(x)
    unlist(lapply(x, function(y) nrow(y))))

#' Save results
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


#' ### Save gene set enrichment results  
save(Rrst, sigG.entrez, bkg.entrez, enrich.rst.fit3, data_merged, 
    enrichCont.rst.fit3, data_contrast, enrichCont.up, data_upreg,
    enrichCont.down, data_downreg, file=paste(getwd(), 
    "Enrichment_KER_Glycogen_Depl.Rdata", sep="/"))



#' ### Run R Script
#+ eval = FALSE
htmlRunR
Enrichment_KER_Glycogen_Depl.R nodes=1,cpus-per-task=1,time=03:00:00,mem=50G \
+Gene Set Enrichment KER Glycogen Diets Over Time
