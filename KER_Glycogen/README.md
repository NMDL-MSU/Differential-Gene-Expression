# Differential-Gene-Expression
## Project KER Glycogen

**Critical pathways limiting glycogen repletion and performance in the horse defined by
proteomic and transcriptomic analyses**

*Objective:*  
To determine the time course of differential expression specific genes and biological pathways during
glycogen depletion and repletion in fit thoroughbred horses fed a high starch (Sweet Feed) compared to an
isocaloric low starch, high fat (Fat) diet. 

### Study design
Five Thoroughbred geldings were given two diets (low and high starch) with appropriate washout periods between diet. 
Muscle biopsies were taken from the gluteal muscle at four-time points: pre-depletion of glycogen, depletion, repletion at 24h and repletion at 72h. 


### Differential Gene Expression Analysis
Determined differentially expressed genes were identified for the following comparisons:  
* Between diets per time point  
* Between diets over time  
* Between time points per diet  
* Between time points per diet over time  

#### &nbsp;&nbsp;&nbsp;&nbsp;DE Analysis:
The scripts used to run the DE analysis can be viewed in detail [here](https://github.com/NMDL-MSU/Differential-Gene-Expression/blob/master/KER_Glycogen/DE_KER_Glycogen.md)

#### &nbsp;&nbsp;&nbsp;&nbsp;DE Results:
The DE results for all comparisons can be downloaded [here](https://github.com/NMDL-MSU/Differential-Gene-Expression/blob/master/KER_Glycogen/Results_DE_Analysis_KER_Glycogen.xlsx?raw=true)


### Gene Set Enrichment Analysis of Differentially Expressed Genes
Pathway analysis for differentially expressed genes was performed with ['clusterprofiler'](https://guangchuangyu.github.io/software/clusterProfiler/). The scripts used to evaluate enriched pathways for each comparison evaluated can be accessed through the following links:
* [Enriched pathways for genes differentially expressed between diet over time](https://htmlpreview.github.io/?https://github.com/NMDL-MSU/Differential-Gene-Expression/blob/master/KER_Glycogen/Pathway_Analysis/Enrichment_KER_Glycogen_Diet/Enrichment_KER_Glycogen_Diet.html)
* Enriched pathways for genes differentially expressed in each diet over time:
    * [Pre-depletion time point as reference](https://htmlpreview.github.io/?https://github.com/NMDL-MSU/Differential-Gene-Expression/blob/master/KER_Glycogen/Pathway_Analysis/Enrichment_KER_Glycogen_PreDepl/Enrichment_KER_Glycogen_PreDepl.html)
    * [Depletion time point as reference](https://htmlpreview.github.io/?https://github.com/NMDL-MSU/Differential-Gene-Expression/blob/master/KER_Glycogen/Pathway_Analysis/Enrichment_KER_Glycogen_Depl/Enrichment_KER_Glycogen_Depl.html)
    

#### &nbsp;&nbsp;&nbsp;&nbsp;Enriched pathway results:
The enriched pathways can be downloaded using the following links:
* [Pathways for genes differentially expressed between diet over time](https://github.com/NMDL-MSU/Differential-Gene-Expression/blob/master/KER_Glycogen/Pathway_Analysis/Enrichment_KER_Glycogen_Diet/KER_Glycogen_GO-Kegg_Results_Diet.xlsx?raw=true)
* Pathways for genes differentially expressed in each diet over time:
  * [Pre-depletion time point as reference](https://github.com/NMDL-MSU/Differential-Gene-Expression/blob/master/KER_Glycogen/Pathway_Analysis/Enrichment_KER_Glycogen_PreDepl/KER_Glycogen_GO-Kegg_Results_Timepoint_Pre-Depletion_Ref.xlsx?raw=true)
  * [Depletion time point as reference](https://github.com/NMDL-MSU/Differential-Gene-Expression/blob/master/KER_Glycogen/Pathway_Analysis/Enrichment_KER_Glycogen_Depl/KER_Glycogen_GO-Kegg_Results_Timepoint_Depletion_Ref.xlsx?raw=true)
