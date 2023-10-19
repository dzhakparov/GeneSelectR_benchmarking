GeneSelectR::set_reticulate_python()
library(dplyr)
library(GeneSelectR)
library(DESeq2)
library(ggplot2)

################################################################################
# 1. Prepare the data
################################################################################
work_dir <- getwd()
data_dir <- file.path(work_dir,'raw-data')
output_dir <- file.path(work_dir, 'results')
clinical <- read.csv(file.path(data_dir, 'clinical.csv'))

# create a new column with extarcted substring
sample_metadata$short_id <- gsub(".*-(A\\w{3})$", "\\1", sample_metadata$patient)
sample_metadata <- sample_metadata[sample_metadata$short_id %in% clinical$X, ]
sample_metadata <- sample_metadata %>% dplyr::filter(!is.na(paper_BRCA_Subtype_PAM50))

# subset the counts matrix
raw_counts <- raw_counts[,colnames(raw_counts) %in% sample_metadata$barcode]
View(sample_metadata[duplicated(sample_metadata$short_id),c('barcode','short_id', 'paper_BRCA_Subtype_PAM50')])

################################################################################
# 2. DESeq2 Analysis
# Basal  Her2  LumA  LumB
# 80    38   188    74
################################################################################
# Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = sample_metadata,
                              design = ~ paper_BRCA_Subtype_PAM50)

# filter the low counts
table(sample_metadata$paper_BRCA_Subtype_PAM50)
smallestGroupSize <- 38 #smallest group is 38 samples
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds)

# Get DE results
res <- results(dds)
filtered_res <- subset(res, padj < 0.001 & abs(log2FoldChange) > log2(5))
filtered_df <- as.data.frame(filtered_res)
deg_list <- rownames(filtered_df)
# Extract VST (Variance-Stabilized Transformed) Data
vsd <- vst(dds, blind = FALSE)  # blind=FALSE uses sample information in transformation
vsd_matrix <- assay(vsd)
vsd_matrix <- t(vsd_matrix)

# do a PCA of VST counts
pca_res <- prcomp(vsd_matrix)
# Calculate explained variance
explained_var <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)
# Plot using ggplot2
ggplot() +
  geom_point(aes(x = pca_res$x[, 1], y = pca_res$x[, 2], color = sample_metadata$paper_BRCA_Subtype_PAM50)) +
  labs(
    x = paste0("PC1: ", explained_var[1], "% variance"),
    y = paste0("PC2: ", explained_var[2], "% variance"),
    title = "PCA of DESeq2 Normalized Counts"
  ) +
  #scale_color_manual(values = c("red", "blue",'darkgreen','purple')) +  # Replace with your color scheme
  theme_minimal()

# do PCA on deseq2-normalized log transformed counts
norm_counts <- counts(dds, normalized=TRUE)
log_norm_counts <- log2(norm_counts + 1)
pca_res <- prcomp(t(log_norm_counts), center = TRUE, scale. = FALSE)
plot(pca_res$x[, 1], pca_res$x[, 2],
     xlab = "PC1",
     ylab = "PC2",
     main = "PCA of DESeq2 Normalized Counts",
     col = sample_metadata$paper_BRCA_Subtype_PAM50)  # Replace with your color vector

write.csv(vsd_matrix, file = file.path(data_dir, 'normalized_vst_matrix.csv'))
saveRDS(sample_metadata, file = file.path(data_dir, 'sample_metadata.rds'))
################################################################################
# 3. GeneSelectR Analysis
################################################################################
GeneSelectR::set_reticulate_python()

exp <- read.csv(file.path(data_dir, 'normalized_vst_matrix.csv'))
rownames(exp) <- exp$X
exp$X <- NULL
sample_metadata <- readRDS(file.path(data_dir, 'sample_metadata.rds'))
sample_metadata <- sample_metadata[match(rownames(exp),rownames(sample_metadata)), ]
sample_metadata$paper_BRCA_Subtype_PAM50 <- as.factor(sample_metadata$paper_BRCA_Subtype_PAM50)
sample_metadata$num_label <- as.integer(sample_metadata$paper_BRCA_Subtype_PAM50)
X <-  exp
y <- sample_metadata %>% select(num_label)

selection_results <- GeneSelectR(X,
                                 y,
                                 njobs = -1,
                                 n_splits = 5,
                                 max_features = 250,
                                 perform_test_split = TRUE,
                                 scoring = 'accuracy',
                                 calculate_permutation_importance = TRUE,
                                 search_type = 'random',
                                 n_iter = 150L)
saveRDS(selection_results, file = file.path(output_dir,'selection_results.rds'))

# plot performance metrics
plot_metrics(selection_results)

# plot feature importance
feature_imps <- plot_feature_importance(selection_results, top_n_features = 20)

# show if there is overlap
overlap <- calculate_overlap_coefficients(selection_results)
plot_overlap_heatmaps(overlap)

#inspect overlap with DEGs
overlap_degs <- calculate_overlap_coefficients(selection_results, custom_lists = list('DEGS' = deg_list))
plot_overlap_heatmaps(overlap_degs)

# annotate the lists
# fecth the backfround list
background <- as.character(colnames(vsd_matrix))

# remove version number from the ensembl ids
background <- gsub("\\..*", "", background)

# annotate the lists
ah <- AnnotationHub::AnnotationHub()
human_ens <- AnnotationHub::query(ah, c("Homo sapiens", "EnsDb"))
human_ens <- human_ens[['AH98047']]
annotations_ahb <- ensembldb::genes(human_ens, return.type = "data.frame") %>%
  dplyr::select(gene_id,gene_name,entrezid,gene_biotype)

custom_list <- list('background' = background,
                    'DEGs' = deg_list)

annotations_df <- annotate_gene_lists(pipeline_results = selection_results,
                                      annotations_ahb = annotations_ahb,
                                      format = 'ensembl_id',
                                      custom_lists = custom_list)

# perform GO Analysis
annotated_GO_inbuilt <- GO_enrichment_analysis(annotations_df,
                                       list_type = 'inbuilt', #run GO enrichment on permutation based selected features
                                       keyType = 'ENSEMBL', # run analysis with ENSEMBLIDs
                                       background = background,
                                       ont = 'BP') # run BP ontology

annotated_GO_permutation <- GO_enrichment_analysis(annotations_df,
                                               list_type = 'permutation', #run GO enrichment on permutation based selected features
                                               keyType = 'ENSEMBL', # run analysis with ENSEMBLIDs
                                               background = background,
                                               ont = 'BP') # run BP ontology

# perform SS analysis
# simplifyEnrichment produces a heatmap
hmap_inbuilt <- run_simplify_enrichment(annotated_GO_inbuilt,
                                method = 'louvain',
                                measure = 'Rel',
                                ont = 'BP',
                                padj_column = 'pvalue',
                                padj_cutoff = 0.01)

hmap_permutation <- run_simplify_enrichment(annotated_GO_permutation,
                                        method = 'louvain',
                                        measure = 'Jiang',
                                        ont = 'BP',
                                        padj_column = 'pvalue',
                                        padj_cutoff = 0.01)


# the best so far:
# hmap_inbuilt <- run_simplify_enrichment(annotated_GO_inbuilt,
#                                         method = 'louvain',
#                                         measure = 'Jiang',
#                                         ont = 'BP',
#                                         padj_column = 'pvalue',
#                                         padj_cutoff = 0.01)
#
# hmap_permutation <- run_simplify_enrichment(annotated_GO_permutation,
#                                             method = 'louvain',
#                                             measure = 'Jiang',
#                                             ont = 'BP',
#                                             padj_column = 'pvalue',
#                                             padj_cutoff = 0.01)
