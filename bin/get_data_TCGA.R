library(TCGAbiolinks)
work_dir <- getwd()
data_dir <- file.path(work_dir,'raw-data')
path_to_folder <- file.path(data_dir)
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")

GDCdownload(query, directory = path_to_folder)
data <- GDCprepare(query, directory = path_to_folder)

# get the raw data (?)
raw_counts <- SummarizedExperiment::assay(data, "unstranded")
sample_metadata <- as.data.frame(SummarizedExperiment::colData(data))
