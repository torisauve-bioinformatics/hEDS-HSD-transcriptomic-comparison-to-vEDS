library(tidyverse)
library(dplyr)
library(tidyr)
library(tibble)
library(limma)
library(edgeR)
library(pheatmap)

# Read all count files
files <- list.files(pattern="\\.txt$", full.names=TRUE)
print(files)

# Read each file
read_counts <- function(f) {
  df <- read.table(f, header=TRUE, skip=1, sep="\t", row.names=1)
  # Last column is the counts
  counts <- df[, ncol(df), drop=FALSE]
  # Clean up sample name from filename
  colnames(counts) <- basename(f) |> 
    gsub("counts_trimmed_|_Aligned.txt", "", x=_)
  return(counts)
}

# Apply read_counts fxn to each file in files list
count_list <- lapply(files, read_counts)

# Merge counts from all of the files
count_matrix <- do.call(cbind, count_list)

# Read in fetchngs sample sheet
sample_sheet <- read.csv("samplesheet.csv") 
ss_needed <- sample_sheet %>% dplyr::select(fastq_2, run_accession, scientific_name) # Keep needed columns from metadata
ss_needed <- ss_needed %>%
  unite("experiment", run_accession, fastq_2, sep = "_", remove = TRUE) # Create name that matches names in count_matrix 

# Create a cleaned meta data
meta_cleaned <- ss_needed %>%
  extract(
    col = scientific_name, # Break scientific name into useful information
    into = c("condition", "patient_id", "rep_type"),
    # Regex Breakdown:
    # (.*)        -> Capture everything up to the last space (Condition)
    # \\s([CP]\\d+) -> Capture a space followed by C or P and digits (Patient ID)
    # _?(.*)?     -> Capture an optional underscore and everything after (Replicate)
    regex = "(.*)\\s([CPH]\\d+)_?(.*)?",
    remove = FALSE
  ) %>%
  mutate(rep_type = ifelse(rep_type == "" | is.na(rep_type), "base", rep_type)) #if rep_type empty, call base


meta_cleaned <- meta_cleaned %>%
  mutate(rep_type = ifelse(rep_type == "" | is.na(rep_type), "base", rep_type)) %>%
  # Add clean condition label based on what's in the condition column
  mutate(condition_clean = case_when(
    grepl("hEDS", condition, ignore.case = TRUE)    ~ "hEDS",
    grepl("HSD",  condition, ignore.case = TRUE)    ~ "HSD",
    grepl("Control", condition, ignore.case = TRUE) ~ "control",
    TRUE ~ NA_character_  # Catches anything unexpected
  ))


# Add labels to differentiate between replicate types 
meta_cleaned <- meta_cleaned %>%
  mutate(collapse_id = case_when(
    rep_type == "technical replicate" ~ paste0(patient_id, "_A"),
    rep_type == "base"                ~ paste0(patient_id, "_A"),
    rep_type == "biological replicate" ~ paste0(patient_id, "_B"),
    TRUE                               ~ patient_id
  ))


meta_cleaned <- meta_cleaned %>%
  mutate(unique_column_name = paste(patient_id, rep_type, sep = "_")) # Create Unique ID for each

# Make it so column names in count_matrix match simplified unique column names 
meta_cleaned <- meta_cleaned %>% mutate(counts_file_name = paste0(experiment, "_Aligned.sortedByCoord.out.txt"))
mapping <- setNames(meta_cleaned$counts_file_name, meta_cleaned$unique_column_name) 
genecounts <- count_matrix %>% rename(any_of(mapping))


write.csv(genecounts, "genecounts.csv", row.names = TRUE) # Convert genecounts to csv


# Align columns to metadata
genecounts_numeric <- genecounts[, meta_cleaned$unique_column_name]
rownames(meta_cleaned) <- meta_cleaned$unique_column_name
genecounts_numeric <- genecounts_numeric[, rownames(meta_cleaned)]

# Log-transform raw counts (add 1 to avoid log(0))
log_counts <- log2(as.matrix(genecounts_numeric) + 1)


gene_vars <- apply(log_counts, 1, var)
log_counts_filtered <- log_counts[gene_vars > 0, ]


# Build DGEList
dge <- DGEList(counts = round(as.matrix(genecounts_numeric)))
# Sum technical replicates (A) biological replicates (B) are left separate
collapse_factor <- factor(meta_cleaned$collapse_id)
dge_collapsed <- sumTechReps(dge, ID = collapse_factor)

# Rebuild metadata after collapse
meta_collapsed <- meta_cleaned %>%
  distinct(collapse_id, .keep_all = TRUE) %>%
  arrange(match(collapse_id, colnames(dge_collapsed$counts)))

rownames(meta_collapsed) <- meta_collapsed$collapse_id

# Set condition as factor with control as reference
meta_collapsed$condition <- factor(meta_collapsed$condition_clean, 
                                   levels = c("control", "hEDS", "HSD"))
# Filter Low Expression Genes
keep <- filterByExpr(dge_collapsed, 
                     group = meta_collapsed$condition_clean,
                     min.count = 10)

dge_filtered <- dge_collapsed[keep, ]

# Normalize 
dge_filtered <- calcNormFactors(dge_filtered, method = "TMM")
dge_filtered$samples


# Define Design matrix
design <- model.matrix(~ 0 + condition, data = meta_collapsed)

colnames(design) <- levels(meta_collapsed$condition) # Rename columns for clarity 



# VOOM and Duplicate correlation (accounts for biological replicates)
v1 <- voom(dge_filtered, design, plot = TRUE)

# Take biological replicates into account
corfit1 <- duplicateCorrelation(v1, design, 
                                block = meta_collapsed$patient_id)

corfit1$consensus
# Second pass with voom correlation estimates
v2 <- voom(dge_filtered, design, plot = TRUE,
           block = meta_collapsed$patient_id,
           correlation = corfit1$consensus)


v2_counts <- v2$E
gene_vars <- apply(v2_counts, 1, var)

top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:500]

# PCA Analysis
pca_res <- prcomp(t(v2_counts[top_genes, ]), scale. = TRUE)
pca_df <- as.data.frame(pca_res$x)
pca_df$SampleID <- rownames(pca_df)

# Calculate variance explained for axis labels
var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
pca_df <- pca_df %>% 
  left_join(meta_collapsed, by = c("SampleID" = "collapse_id"))

pca_df <- pca_df %>%
  mutate(srx_numeric = as.numeric(gsub("SRX([0-9]+)_.*", "\\1", experiment)))


# Plot using ggplot2
ggplot(pca_df, aes(x = PC1, y = PC2, color=condition_clean, label = SampleID)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5) +
  labs(title = "PCA of Top 500 Most Variable Genes",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_minimal()+
  theme(
    plot.title   = element_text(size = 16),    # title
    axis.title   = element_text(size = 13),    # x and y axis labels
    axis.text    = element_text(size = 11),    # axis tick labels
    legend.title = element_text(size = 12),    # legend title
    legend.text  = element_text(size = 10)     # legend item labels
  )

ggplot(pca_df, aes(x = PC1, y = PC2, color=srx_numeric, label = SampleID)) +
  geom_point(size = 3) +
  geom_text(vjust = 2) +
  labs(title = "PCA of Top 500 Most Variable Genes",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_minimal() +
  theme(
    plot.title   = element_text(size = 16),    # title
    axis.title   = element_text(size = 13),    # x and y axis labels
    axis.text    = element_text(size = 11),    # axis tick labels
    legend.title = element_text(size = 12),    # legend title
    legend.text  = element_text(size = 10)     # legend item labels
  )

ggplot(pca_df, aes(x = PC3, y = PC4, color=condition_clean, label = SampleID)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5) +
  labs(title = "PCA of Top 500 Most Variable Genes",
       x = paste0("PC3 (", var_explained[3], "%)"),
       y = paste0("PC4 (", var_explained[4], "%)")) +
  theme_minimal() +
  theme(
    plot.title   = element_text(size = 16),    # title
    axis.title   = element_text(size = 13),    # x and y axis labels
    axis.text    = element_text(size = 11),    # axis tick labels
    legend.title = element_text(size = 12),    # legend title
    legend.text  = element_text(size = 10)     # legend item labels
  )


ggplot(pca_df, aes(x = PC1, y = PC3, color=condition_clean, label = SampleID)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5) +
  labs(title = "PCA of Top 500 Most Variable Genes",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC3 (", var_explained[3], "%)")) +
  theme_minimal() +
  theme(
    plot.title   = element_text(size = 16),    # title
    axis.title   = element_text(size = 13),    # x and y axis labels
    axis.text    = element_text(size = 11),    # axis tick labels
    legend.title = element_text(size = 12),    # legend title
    legend.text  = element_text(size = 10)     # legend item labels
  )



set.seed(123) # For reproducibility
k_val <- 3

# Run k-means on the first 5 PCs (which capture most variance)
km_res <- kmeans(pca_res$x[, 1:5], centers = k_val, nstart = 25)

# Add the cluster assignments back to your plotting dataframe
pca_df$cluster <- as.factor(km_res$cluster)


wss <- sapply(1:10, function(k){
  kmeans(pca_res$x[, 1:5], centers = k, nstart = 25)$tot.withinss
})

ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster, shape = condition_clean)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(title = paste0("K-means Clustering (k=", k_val, ")"),
       subtitle = "Colors = K-means Clusters | Shapes = Actual Diagnosis",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_minimal()



genecounts_numeric <- genecounts_numeric[, meta_cleaned$unique_column_name]
# After sumTechReps, update the metadata order:
meta_collapsed <- meta_collapsed[match(colnames(dge_collapsed), meta_collapsed$collapse_id), ]
# Second pass duplicateCorrection with updated voom object
corfit2 <- duplicateCorrelation(v2, design,
                                block = meta_collapsed$patient_id)

# Fit Model
fit <- lmFit(v2, design,
             block = meta_collapsed$patient_id,
             correlation = corfit2$consensus)

contrasts <- makeContrasts(
  HED_vs_control = hEDS - control,
  HSD_vs_control = HSD - control,
  HED_vs_HSD     = hEDS - HSD,
  levels = design
)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)


get_DEGs <- function(fit, coef, fdr = 0.05, lfc = 1) {
  # Extract results
  res <- topTable(fit, coef = coef, n = Inf, adjust.method = "BH")
  
  # Capture gene_id from rownames explicitly
  res$gene_id <- rownames(res)
  
  
  # Restore the row names if they were lost during merge
  rownames(res) <- res$gene_id
  
  # Filter and Sort
  res_filtered <- res %>%
    filter(adj.P.Val < fdr, abs(logFC) > lfc) %>%
    arrange(adj.P.Val)
  
  return(res_filtered)
}

# Apply function for each comparison
DEGs_HED_vs_ctrl <- get_DEGs(fit2, "HED_vs_control")
nrow(DEGs_HED_vs_ctrl)
DEGs_HSD_vs_ctrl <- get_DEGs(fit2, "HSD_vs_control")
nrow(DEGs_HSD_vs_ctrl)
DEGs_HED_vs_HSD  <- get_DEGs(fit2, "HED_vs_HSD")
nrow(DEGs_HED_vs_HSD)

# Convert to CSV file
write.csv(DEGs_HED_vs_ctrl, "DEGs_HEDS_vs_control.csv", row.names = FALSE)
write.csv(DEGs_HSD_vs_ctrl, "DEGs_HSD_vs_control.csv", row.names = FALSE)
write.csv(DEGs_HED_vs_HSD,  "DEGs_HEDS_vs_HSD.csv",    row.names = FALSE)

# Function to make volcano plots for each comparison
make_volcano <- function(fit, coef, title, p_thresh=0.05, lfc_thresh=1){
  res <- topTable(fit, coef = coef, n = Inf)
  res$gene_id <- rownames(res)
  
  # Add column for coloring
  res <- res %>%
    mutate(sig = case_when(
      adj.P.Val < p_thresh & logFC > lfc_thresh ~ "Up",
      adj.P.Val < p_thresh & logFC < -lfc_thresh ~ "Down",
      TRUE ~ "Not Sig"
    ))
  # Plot
  ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Down" = "royalblue", "Not Sig" = "grey", "Up" = "firebrick")) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed") +
    labs(title = title,
         x = "Log2 Fold Change",
         y = "-log10 Adjusted P-value") +
    theme_minimal()
  
}

# Use function to create volcano plots
make_volcano(fit2, "HED_vs_control", "Volcano Plot: hEDS vs Control")
make_volcano(fit2, "HSD_vs_control", "Volcano Plot: HSD vs Control")
make_volcano(fit2, "HED_vs_HSD",     "Volcano Plot: hEDS vs HSD")

# Load packages for enrichment 
library(clusterProfiler)
library(org.Hs.eg.db)

run_GO_with_conversion <- function(deg_results) {
  
  # Convert Symbols to ENSEMBL or ENTREZ IDs
  # This function maps your symbols to the standard database
  ids <- bitr(deg_results$gene_id, 
              fromType = "SYMBOL", 
              toType   = c("ENSEMBL", "ENTREZID"), 
              OrgDb    = org.Hs.eg.db)
  
  # Run the Enrichment using the new IDs
  ego <- enrichGO(gene          = ids$ENSEMBL,
                  universe      = NULL, 
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENSEMBL",
                  ont           = "ALL", 
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2)
  
  # Add the readable symbols back to the output for easier reading
  ego <- setReadable(ego, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
  return(ego)
}

# GO results for each comparison 
heds_GO_results = run_GO_with_conversion(DEGs_HED_vs_ctrl)
nrow(heds_GO_results)
hsd_GO_results = run_GO_with_conversion(DEGs_HSD_vs_ctrl)
heds_hsd_GO_results = run_GO_with_conversion(DEGs_HED_vs_HSD)

write.csv(heds_GO_results, "GO_terms_hEDS_vs_control.csv", row.names = FALSE)
write.csv(hsd_GO_results,  "GO_terms_HSD_vs_control.csv",  row.names = FALSE)
write.csv(heds_hsd_GO_results,  "GO_terms_hEDS_vs_HSD.csv",     row.names = FALSE)

# Create bar plots for each
plot_go_from_combined <- function(ego_all, n = 10, title = "Top GO Terms by Ontology") {
  
  df_all <- as.data.frame(ego_all)
  
  # Split by ontology and take top n each
  df_combined <- df_all %>%
    group_by(ONTOLOGY) %>%
    slice_min(order_by = p.adjust, n = n) %>%
    ungroup() %>%
    mutate(
      log10p = -log10(p.adjust),
      Description = reorder(Description, log10p)  # order bars by significance
    )
  
  ggplot(df_combined, aes(x = log10p, y = Description, fill = ONTOLOGY)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Set1", name = "GO Category") +
    labs(
      title = title,
      x = "-log10(p.adjust)",
      y = ""
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title  = element_text(size = 11, face = "bold"),
      legend.position = "right"
    )
}

plot_go_from_combined(heds_GO_results,     title = "Top GO Terms — hEDS vs Control")


# Create barplot based on GO category 
GO_Type_Distribution <- function(result){
  res_df <- as.data.frame(result) # Convert result to df
  
  go_type_counts <- res_df %>%
    group_by(ONTOLOGY) %>%
    summarise(count = n()) %>%
    mutate(percentage = count / sum(count))
  
  ggplot(go_type_counts, aes(x = "", y = count, fill = ONTOLOGY)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    theme_void() + # Remove background, grid, and labels
    labs(title = "Distribution of GO Terms: hEDS",
         fill = "GO Type",
         x = NULL,
         y = NULL) +
    geom_text(aes(label = paste0(round(percentage * 100), "%")),
              position = position_stack(vjust = 0.5)) +
    scale_fill_brewer(palette = "Set1") +
    theme(
      plot.title   = element_text(size = 16, hjust = 0.5),
      legend.title = element_text(size = 14),
      legend.text  = element_text(size = 12))
}

GO_Type_Distribution(heds_GO_results)
GO_Type_Distribution
GO_Type_Distribution

