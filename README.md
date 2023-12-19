## Description 
This repository contains analysis rundown of [GeneSelectR](https://github.com/dzhakparov/GeneSelectR) applied to [TCGA-BRCA](https://portal.gdc.cancer.gov/projects/TCGA-BRCA) RNA sequencing data. 

## Repository Structure Description 
Repository contains following folders: 
- *bin* - contains scripts that were used to get the TCGA-BRCA dataset 
- *docker* - contains Dockerfile to reproduce the analysis 
- *docs* - contains R Markdown of the analysis and its rendered versions in different formats 
- *raw-data* - contains raw data files for the analysis 
- *results* - contains outputs and result files produced during the analysis

## Aim of the Benchmarking and Outcomes 
The goal of the analysis is to showcase the GeneSelectR package applied to a common problem and compare how it performs in comparison to traditional differential gene expression analysis (DGE). According to the results, boruta feature selection method, a default method of GeneSelectR, showed a better performance in classification task with higher mean accuracy compared to other methods. Additionally, it yielded a list of 106 transcripts that were biologically relevant to cancer metabolic and immune processes. On the other hand DGE-derived list results in more than a thousand of possible transcripts. Please refer to the R Mardown file for more details.  

## Contact details 
If you have any questions regarding the analysis or the GeneSelectR package please email to damir.zhakparov@uzh.ch. If you have an issue with GeneSelectR package please make an issue in the official [repository](https://github.com/dzhakparov/GeneSelectR). 
