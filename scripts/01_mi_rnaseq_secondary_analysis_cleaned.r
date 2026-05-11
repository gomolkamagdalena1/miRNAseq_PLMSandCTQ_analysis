## ---
## Purpose: Secondary analysis of human small RNA sequencing data
## Project: Small RNA sequencing of serum and sperm samples associated with adverse childhood experiences
## Original author: Adria-Jaume Roura
## Adapted for publication repository by: Magdalena Gomolka
## Notes:
##   - Raw sequencing data and processed count files are deposited in ArrayExpress/BioStudies.
##   - Sensitive participant-level metadata are not included in this repository.
##   - Before running, place input files in the local input directory and adjust the paths below.
## ---

## =============================
## 0. User settings
## =============================

# Run this script from the project root directory.
# Expected local structure:
#   input/
#   results/
#   figures/

input_dir <- "input"
results_dir <- "results"
figures_dir <- "figures"

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Input files. Use local, non-identifiable filenames.
counts_file <- file.path(input_dir, "primary_UMI_counts.csv")
registry_file <- file.path(input_dir, "illumina_sample_registry.tsv")
metadata_file <- file.path(input_dir, "sample_metadata.tsv")
tarbase_file <- file.path(input_dir, "tarbase_data.tsv")

# Samples excluded before downstream analysis.
# Reasons included low sequencing quality or library preparation issues.
excluded_sample_files <- c(
  "AJM_F54_0004_P28_V11",
  "AJM_F54_0007_P28_V11",
  "ALI_F53_0126_P28_V11",
  "ALI_F54_0058_P28_V11",
  "ALI_F54_0045_P28_V11"
)

# Minimum total UMI count across all samples for retaining a small RNA feature.
minimum_feature_sum <- 20

# Set to TRUE only if the target/enrichment input file is available locally.
RUN_ENRICHMENT <- FALSE


## =============================
## 1. Libraries
## =============================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(corrplot)
  library(pheatmap)
  library(RColorBrewer)
  library(EnhancedVolcano)
  library(BiocParallel)
})

# Optional packages used only for exploratory target/enrichment analysis.
if (RUN_ENRICHMENT) {
  suppressPackageStartupMessages({
    library(biomaRt)
    library(clusterProfiler)
    library(ReactomePA)
    library(org.Hs.eg.db)
    library(enrichplot)
  })
}

# More portable parallelization. MulticoreParam is not available on Windows.
if (.Platform$OS.type == "windows") {
  register(SerialParam())
} else {
  register(MulticoreParam(workers = 8))
}


## =============================
## 2. Helper functions
## =============================

clean_feature_names <- function(x) {
  gsub(pattern = "/gb/.*", replacement = "", x = x)
}

save_pdf <- function(filename, plot_expr, width = 10, height = 8) {
  pdf(file = filename, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  force(plot_expr)
}

make_coldata <- function(metadata) {
  coldata <- data.frame(
    row.names = metadata$file,
    tissue = metadata$sperm_serum,
    sample_name = metadata$name_of_sample,
    sample_id = metadata$sample,
    gender = metadata$gender,
    group = metadata$group,
    stringsAsFactors = FALSE
  )
  coldata$group <- factor(coldata$group)
  coldata$gender <- factor(coldata$gender)
  coldata$tissue <- factor(coldata$tissue)
  coldata
}

summarize_deseq_results <- function(res_df) {
  message(
    "Differentially UP-regulated features with q < 0.05: ",
    sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE)
  )
  message(
    "Differentially UP-regulated features with q < 0.10: ",
    sum(res_df$padj < 0.10 & res_df$log2FoldChange > 0, na.rm = TRUE)
  )
  message(
    "Differentially DOWN-regulated features with q < 0.05: ",
    sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)
  )
  message(
    "Differentially DOWN-regulated features with q < 0.10: ",
    sum(res_df$padj < 0.10 & res_df$log2FoldChange < 0, na.rm = TRUE)
  )
}

run_deseq2 <- function(count_matrix, metadata, group_1, group_2, tissue_type) {
  message("Running DESeq2: ", tissue_type, ", ", group_1, " vs ", group_2)

  coldata <- make_coldata(metadata)
  coldata <- coldata[coldata$tissue == tissue_type, , drop = FALSE]

  # Keep only the two groups used in the contrast.
  coldata <- coldata[coldata$group %in% c(group_1, group_2), , drop = FALSE]
  coldata$group <- droplevels(coldata$group)
  coldata$gender <- droplevels(coldata$gender)

  counts_subset <- count_matrix[, rownames(coldata), drop = FALSE]
  counts_subset <- round(as.matrix(counts_subset))
  storage.mode(counts_subset) <- "integer"

  message("Number of samples: ", nrow(coldata))
  message("Number of features: ", nrow(counts_subset))

  if (tissue_type == "Serum") {
    dds <- DESeqDataSetFromMatrix(
      countData = counts_subset,
      colData = coldata,
      design = ~ gender + group
    )
  } else if (tissue_type == "Sperm") {
    dds <- DESeqDataSetFromMatrix(
      countData = counts_subset,
      colData = coldata,
      design = ~ group
    )
  } else {
    stop("Unknown tissue_type. Expected 'Serum' or 'Sperm'.")
  }

  dds <- DESeq(dds)
  vst_obj <- varianceStabilizingTransformation(dds, blind = FALSE)

  res <- results(
    dds,
    contrast = c("group", group_1, group_2),
    independentFiltering = TRUE
  )
  res_df <- as.data.frame(res[order(res$padj), ])
  res_df$feature <- rownames(res_df)
  res_df <- res_df[, c("feature", setdiff(colnames(res_df), "feature"))]

  summarize_deseq_results(res_df)

  out_file <- file.path(results_dir, paste0("DESeq2_", group_1, "_vs_", group_2, ".tsv"))
  write.table(res_df, out_file, sep = "\t", quote = FALSE, row.names = FALSE)

  list(
    dds = dds,
    results = res_df,
    vst = vst_obj,
    group_1 = group_1,
    group_2 = group_2,
    tissue_type = tissue_type
  )
}

plot_pca_pdf <- function(vst_obj, intgroup, output_file, color_values = NULL, width = 10, height = 8) {
  p <- plotPCA(vst_obj, intgroup = intgroup) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 14, face = "plain", colour = "black"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      legend.title = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    )

  if (!is.null(color_values)) {
    p <- p + scale_color_manual(values = color_values)
  }

  pdf(output_file, width = width, height = height)
  print(p)
  dev.off()
}

plot_volcano_pdf <- function(comp_object, output_file, width = 12, height = 10) {
  res <- comp_object$results
  rownames(res) <- clean_feature_names(res$feature)

  p <- EnhancedVolcano(
    res,
    lab = rownames(res),
    x = "log2FoldChange",
    y = "padj",
    pCutoff = 0.05,
    FCcutoff = 1,
    xlab = bquote(~Log[2]~ "fold change"),
    ylab = "Statistical significance (-log10 q-value)",
    labSize = 4,
    labCol = "black",
    labFace = "bold",
    boxedLabels = FALSE,
    colAlpha = 4 / 5,
    legendPosition = "right",
    legendLabSize = 10,
    legendIconSize = 3.0,
    drawConnectors = TRUE,
    widthConnectors = 0.1,
    colConnectors = "red2",
    title = "",
    subtitle = "",
    axisLabSize = 13,
    caption = NULL,
    xlim = c(-7, 7),
    ylim = c(0, 10),
    legendVisible = FALSE,
    col = c("grey30", "grey30", "royalblue", "red2"),
    border = "full",
    borderWidth = 1.5,
    borderColour = "black"
  ) +
    theme(
      axis.text.x = element_text(size = 18, color = "black"),
      axis.title = element_text(size = 18, color = "black", face = "plain"),
      axis.text.y = element_text(vjust = 0.5, size = 18, color = "black")
    )

  pdf(output_file, width = width, height = height)
  print(p)
  dev.off()
}

plot_deg_heatmap_pdf <- function(comp_object, output_file, annotation_colors = NULL, width = 12, height = 10) {
  res <- comp_object$results
  deg <- res[!is.na(res$padj) & res$padj < 0.10, , drop = FALSE]

  if (nrow(deg) == 0) {
    warning("No features with padj < 0.10 for ", output_file)
    return(invisible(NULL))
  }

  mat <- assay(comp_object$vst)[deg$feature, , drop = FALSE]
  anno <- as.data.frame(colData(comp_object$dds)[, "group", drop = FALSE])
  colnames(anno) <- "group"
  anno <- anno[order(anno$group), , drop = FALSE]

  mat <- mat[, rownames(anno), drop = FALSE]

  pdf(output_file, width = width, height = height)
  pheatmap(
    mat,
    annotation_col = anno,
    show_colnames = FALSE,
    clustering_method = "ward.D2",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    show_rownames = TRUE,
    scale = "row",
    fontsize_col = 7,
    angle_col = 90,
    fontsize_row = 8.5,
    cluster_cols = FALSE,
    annotation_colors = annotation_colors,
    border_color = NA
  )
  dev.off()
}

plot_normalized_counts_pdf <- function(comp_object, output_file) {
  res <- comp_object$results
  deg <- res[!is.na(res$padj) & res$padj < 0.10, , drop = FALSE]

  if (nrow(deg) == 0) {
    warning("No features with padj < 0.10 for ", output_file)
    return(invisible(NULL))
  }

  pdf(output_file)
  for (feature_id in deg$feature) {
    d <- plotCounts(
      comp_object$dds,
      gene = feature_id,
      intgroup = "group",
      returnData = TRUE
    )

    p <- ggplot(d, aes(x = group, y = count, fill = group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_point(position = position_jitter(width = 0.1, height = 0), shape = 20, size = 4) +
      stat_summary(fun = mean, geom = "point", shape = 4, size = 5) +
      labs(title = paste("small RNA-seq:", clean_feature_names(feature_id)), x = "", y = "Normalized counts") +
      theme_bw() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold")
      )

    print(p)
  }
  dev.off()
}


## =============================
## 3. Load input data
## =============================

# Primary QIAGEN GeneGlobe output containing UMI counts.
miRNA_reads <- read.csv(
  counts_file,
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Discard raw read count columns and retain UMI columns.
miRNA_reads <- miRNA_reads %>%
  dplyr::select(-contains("READs")) %>%
  mutate(across(everything(), ~ as.numeric(as.character(.x))))

# Remove sample suffix produced by sequencing/core facility naming, if present.
colnames(miRNA_reads) <- gsub(
  pattern = "_S[0-9]+\\.UMIs$",
  replacement = "",
  x = colnames(miRNA_reads)
)

# Remove excluded samples if present in the count matrix.
miRNA_reads <- miRNA_reads[, !colnames(miRNA_reads) %in% excluded_sample_files, drop = FALSE]
miRNA_UMIs <- miRNA_reads

# Illumina registry mapping sequencing file IDs to sample IDs.
illumina_registry <- read.delim(
  registry_file,
  header = FALSE,
  stringsAsFactors = FALSE
)
colnames(illumina_registry) <- c("file", "sample")

# Sample metadata. Expected columns:
# sample, sperm_serum, type_of_sample, gender, name_of_sample
metadata <- read.delim(
  metadata_file,
  header = FALSE,
  stringsAsFactors = FALSE
)
colnames(metadata) <- c("sample", "sperm_serum", "type_of_sample", "gender", "name_of_sample")

metadata <- merge(metadata, illumina_registry, by = "sample")
metadata$type_of_sample <- gsub(" ", "_", metadata$type_of_sample)
metadata$name_of_sample <- gsub(" ", "-", metadata$name_of_sample)
metadata$group <- paste(metadata$sperm_serum, metadata$type_of_sample, sep = "_")

metadata <- metadata[!metadata$file %in% excluded_sample_files, , drop = FALSE]

metadata$group <- gsub("Sperm_Control", "CTQ0", metadata$group)
metadata$group <- gsub("Sperm_1_Stress_event", "CTQ1", metadata$group)
metadata$group <- gsub("Sperm_2_Stress_events", "CTQ2", metadata$group)
metadata$group <- gsub("Serum_Control", "CTRL", metadata$group)
metadata$group <- gsub("Serum_PLMS", "PLMS", metadata$group)

message("Group summary:")
print(table(metadata$group))

# Ensure count matrix and metadata are aligned.
missing_from_counts <- setdiff(metadata$file, colnames(miRNA_UMIs))
if (length(missing_from_counts) > 0) {
  stop("The following metadata samples are missing from the count matrix: ", paste(missing_from_counts, collapse = ", "))
}

miRNA_UMIs <- miRNA_UMIs[, metadata$file, drop = FALSE]


## =============================
## 4. Quality assessment
## =============================

# Correlation matrix across all retained samples.
cor_matrix <- cor(as.matrix(miRNA_UMIs), use = "complete.obs")
rownames(cor_matrix) <- metadata$group
colnames(cor_matrix) <- metadata$group

cohort_colors <- c(
  CTRL = "chartreuse4",
  CTQ0 = "aquamarine4",
  PLMS = "brown1",
  CTQ1 = "firebrick1",
  CTQ2 = "firebrick4"
)

gender_colors <- c(
  Female = "blueviolet",
  Male = "orange"
)

label_colors_group <- rownames(cor_matrix)
label_colors_group <- gsub("CTRL", cohort_colors["CTRL"], label_colors_group)
label_colors_group <- gsub("CTQ0", cohort_colors["CTQ0"], label_colors_group)
label_colors_group <- gsub("PLMS", cohort_colors["PLMS"], label_colors_group)
label_colors_group <- gsub("CTQ1", cohort_colors["CTQ1"], label_colors_group)
label_colors_group <- gsub("CTQ2", cohort_colors["CTQ2"], label_colors_group)

pdf(file = file.path(figures_dir, "cohort_correlation_groups.pdf"), width = 15, height = 15)
corrplot(
  cor_matrix,
  is.corr = TRUE,
  type = "upper",
  method = "square",
  tl.cex = 0.6,
  tl.col = label_colors_group,
  tl.srt = 45
)
dev.off()


## =============================
## 5. Pre-filtering
## =============================

keep <- rowSums(miRNA_UMIs, na.rm = TRUE) >= minimum_feature_sum
miRNA_UMIs <- miRNA_UMIs[keep, , drop = FALSE]
miRNA_UMIs <- miRNA_UMIs[complete.cases(miRNA_UMIs), , drop = FALSE]

message("Number of features retained after pre-filtering: ", nrow(miRNA_UMIs))


## =============================
## 6. Differential expression analyses
## =============================

comp_serum_plms_vs_ctrl <- run_deseq2(
  count_matrix = miRNA_UMIs,
  metadata = metadata,
  group_1 = "PLMS",
  group_2 = "CTRL",
  tissue_type = "Serum"
)

comp_sperm_ctq1_vs_ctq0 <- run_deseq2(
  count_matrix = miRNA_UMIs,
  metadata = metadata,
  group_1 = "CTQ1",
  group_2 = "CTQ0",
  tissue_type = "Sperm"
)

comp_sperm_ctq2_vs_ctq0 <- run_deseq2(
  count_matrix = miRNA_UMIs,
  metadata = metadata,
  group_1 = "CTQ2",
  group_2 = "CTQ0",
  tissue_type = "Sperm"
)

comp_sperm_ctq2_vs_ctq1 <- run_deseq2(
  count_matrix = miRNA_UMIs,
  metadata = metadata,
  group_1 = "CTQ2",
  group_2 = "CTQ1",
  tissue_type = "Sperm"
)


## =============================
## 7. PCA plots
## =============================

plot_pca_pdf(
  comp_serum_plms_vs_ctrl$vst,
  intgroup = "group",
  output_file = file.path(figures_dir, "PCA_serum_group.pdf"),
  color_values = cohort_colors
)

plot_pca_pdf(
  comp_serum_plms_vs_ctrl$vst,
  intgroup = "gender",
  output_file = file.path(figures_dir, "PCA_serum_gender.pdf"),
  color_values = gender_colors
)

plot_pca_pdf(
  comp_sperm_ctq1_vs_ctq0$vst,
  intgroup = "group",
  output_file = file.path(figures_dir, "PCA_sperm_CTQ1_vs_CTQ0_group.pdf"),
  color_values = cohort_colors
)

plot_pca_pdf(
  comp_sperm_ctq2_vs_ctq0$vst,
  intgroup = "group",
  output_file = file.path(figures_dir, "PCA_sperm_CTQ2_vs_CTQ0_group.pdf"),
  color_values = cohort_colors
)


## =============================
## 8. Heatmaps and volcano plots
## =============================

plot_deg_heatmap_pdf(
  comp_serum_plms_vs_ctrl,
  output_file = file.path(figures_dir, "heatmap_DE_miRNA_q01_PLMS_vs_CTRL.pdf"),
  annotation_colors = list(group = cohort_colors)
)

plot_deg_heatmap_pdf(
  comp_sperm_ctq1_vs_ctq0,
  output_file = file.path(figures_dir, "heatmap_DE_miRNA_q01_CTQ1_vs_CTQ0.pdf"),
  annotation_colors = list(group = cohort_colors),
  width = 8,
  height = 8
)

plot_deg_heatmap_pdf(
  comp_sperm_ctq2_vs_ctq0,
  output_file = file.path(figures_dir, "heatmap_DE_miRNA_q01_CTQ2_vs_CTQ0.pdf"),
  annotation_colors = list(group = cohort_colors)
)

plot_volcano_pdf(
  comp_serum_plms_vs_ctrl,
  output_file = file.path(figures_dir, "volcano_DE_miRNA_PLMS_vs_CTRL.pdf")
)

plot_volcano_pdf(
  comp_sperm_ctq1_vs_ctq0,
  output_file = file.path(figures_dir, "volcano_DE_miRNA_CTQ1_vs_CTQ0.pdf")
)

plot_volcano_pdf(
  comp_sperm_ctq2_vs_ctq0,
  output_file = file.path(figures_dir, "volcano_DE_miRNA_CTQ2_vs_CTQ0.pdf")
)


## =============================
## 9. Normalized count plots
## =============================

plot_normalized_counts_pdf(
  comp_serum_plms_vs_ctrl,
  output_file = file.path(figures_dir, "plotCounts_PLMS_vs_CTRL_q01.pdf")
)

plot_normalized_counts_pdf(
  comp_sperm_ctq1_vs_ctq0,
  output_file = file.path(figures_dir, "plotCounts_CTQ1_vs_CTQ0_q01.pdf")
)

plot_normalized_counts_pdf(
  comp_sperm_ctq2_vs_ctq0,
  output_file = file.path(figures_dir, "plotCounts_CTQ2_vs_CTQ0_q01.pdf")
)


## =============================
## 10. Exploratory overlap analysis
## =============================

deg_serum <- comp_serum_plms_vs_ctrl$results %>%
  filter(!is.na(padj), padj < 0.10)

deg_sperm <- comp_sperm_ctq2_vs_ctq0$results %>%
  filter(!is.na(padj), padj < 0.10)

deg_serum_sperm <- merge(deg_serum, deg_sperm, by = "feature", suffixes = c("_serum", "_sperm"))

write.table(
  deg_serum_sperm,
  file.path(results_dir, "DEG_serum_sperm_intersection.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


## =============================
## 11. Optional exploratory target/enrichment analysis
## =============================

if (RUN_ENRICHMENT) {
  mart <- useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "hsapiens_gene_ensembl",
    host = "www.ensembl.org",
    ensemblRedirect = FALSE
  )

  tarbase_data <- read.table(
    tarbase_file,
    header = TRUE,
    sep = "\t",
    na.strings = "NA",
    dec = ".",
    strip.white = TRUE,
    stringsAsFactors = FALSE
  )

  tarbase_data <- subset(
    tarbase_data,
    species == "Homo sapiens" & direct_indirect == "DIRECT" & up_down == "DOWN"
  )

  run_target_enrichment <- function(comp_object, comparison_label) {
    res <- comp_object$results
    res$feature_clean <- clean_feature_names(res$feature)

    up_features <- res %>%
      filter(!is.na(padj), padj < 0.05, log2FoldChange > 0) %>%
      pull(feature_clean)

    down_features <- res %>%
      filter(!is.na(padj), padj < 0.05, log2FoldChange < 0) %>%
      pull(feature_clean)

    feature_sets <- list(up = up_features, down = down_features)

    for (direction in names(feature_sets)) {
      selected_features <- feature_sets[[direction]]
      if (length(selected_features) == 0) next

      targets <- tarbase_data[tarbase_data$mirna %in% selected_features, ]
      targets <- targets[!duplicated(targets$geneName), ]

      if (nrow(targets) == 0) next

      gene_info <- getBM(
        values = targets$geneName,
        filters = "hgnc_symbol",
        mart = mart,
        attributes = c("hgnc_symbol", "entrezgene_id")
      )

      entrez_ids <- unique(na.omit(gene_info$entrezgene_id))
      if (length(entrez_ids) == 0) next

      go <- enrichGO(
        gene = entrez_ids,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP"
      )

      reactome <- enrichPathway(
        gene = entrez_ids,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        organism = "human"
      )

      pdf(file.path(figures_dir, paste0("GO_", direction, "_miRNA_targets_", comparison_label, ".pdf")), width = 10, height = 10)
      print(dotplot(go, showCategory = 20, title = paste("GO:", direction, "miRNA targets", comparison_label)))
      dev.off()

      pdf(file.path(figures_dir, paste0("Reactome_", direction, "_miRNA_targets_", comparison_label, ".pdf")), width = 10, height = 10)
      print(dotplot(reactome, showCategory = 20, title = paste("Reactome:", direction, "miRNA targets", comparison_label)))
      dev.off()
    }
  }

  run_target_enrichment(comp_serum_plms_vs_ctrl, "PLMS_vs_CTRL")
  run_target_enrichment(comp_sperm_ctq2_vs_ctq0, "CTQ2_vs_CTQ0")
}


## =============================
## 12. Session information
## =============================

writeLines(capture.output(sessionInfo()), file.path(results_dir, "sessionInfo.txt"))
