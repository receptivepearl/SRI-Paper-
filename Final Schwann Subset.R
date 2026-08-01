# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)

target_genes <- c("ERBB3", "PMP22")
expr_threshold <- 0 

# Check that both genes exist
missing_genes <- setdiff(target_genes, rownames(Merged_Joined))

if (length(missing_genes) > 0) {
  stop(paste("Missing genes in dataset:", paste(missing_genes, collapse = ", ")))
}

print("Both ERBB3 and PMP22 found in dataset.")

expr_data <- FetchData(
  Merged_Joined,
  vars = target_genes,
  layer = "data"
)


# DGE and GSEA Analysis and Visuals

# ============================================================
# SCHWANN CELL SEX-BIASED DGE + GSEA ANALYSIS
# Male vs Female PDAC Schwann Cells
# ERBB3+ / PMP22+ True Schwann subset
# ============================================================


# ============================================================
# 0. LOAD / INSTALL REQUIRED PACKAGES
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c(
  "Seurat", "dplyr", "ggplot2", "ggrepel",
  "msigdbr", "tibble", "tidyr", "patchwork"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("fgsea", quietly = TRUE)) {
  BiocManager::install("fgsea", ask = FALSE, update = TRUE)
}

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(msigdbr)
library(fgsea)
library(tibble)
library(tidyr)
library(patchwork)


# ============================================================
# 1. LOAD TRUE SCHWANN SUBSET
# ============================================================

true_schwann_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/true_schwann_ERBB3_PMP22_subset.rds"
)

DefaultAssay(true_schwann_subset) <- "RNA"

# Optional for Seurat v5 if layers are split
# true_schwann_subset <- JoinLayers(true_schwann_subset)

print(true_schwann_subset)


# ============================================================
# 2. ADD SEX METADATA TO SCHWANN SUBSET
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",   "Sample_02" = "Male",   "Sample_03" = "Male",
  "Sample_04" = "Female", "Sample_05" = "Female", "Sample_06" = "Female",
  "Sample_07" = "Male",   "Sample_08" = "Female", "Sample_09" = "Male",
  "Sample_10" = "Female", "Sample_11" = "Female", "Sample_12" = "Male",
  "Sample_13" = "Female", "Sample_14" = "Male",   "Sample_15" = "Female",
  "Sample_16" = "Male",   "Sample_17" = "Female"
)

true_schwann_subset$Sex <- unname(
  sex_mapping[as.character(true_schwann_subset$orig.ident)]
)

# Keep only cells with clear sex labels
true_schwann_subset <- subset(
  true_schwann_subset,
  subset = Sex %in% c("Male", "Female")
)

true_schwann_subset$Sex <- factor(
  true_schwann_subset$Sex,
  levels = c("Female", "Male")
)

print("Schwann cells by sex:")
print(table(true_schwann_subset$Sex))

print("Schwann cells by patient and sex:")
print(table(true_schwann_subset$orig.ident, true_schwann_subset$Sex))


# ============================================================
# 3. SAFETY CHECK: MAKE SURE BOTH SEXES HAVE CELLS
# ============================================================

sex_counts <- table(true_schwann_subset$Sex)

if (!all(c("Male", "Female") %in% names(sex_counts))) {
  stop("Both Male and Female Schwann cells are required for DGE.")
}

if (any(sex_counts[c("Male", "Female")] < 3)) {
  stop("One sex has fewer than 3 Schwann cells. DGE is not reliable or may fail.")
}


# ============================================================
# 4. FILTER GENES BY AVERAGE EXPRESSION > 0.5 IN BOTH SEXES
# ============================================================

print("Filtering genes by average expression > 0.5 in both male and female Schwann cells...")

avg_expression_by_sex <- AverageExpression(
  true_schwann_subset,
  assays = "RNA",
  group.by = "Sex",
  slot = "data"
)$RNA

print(colnames(avg_expression_by_sex))

genes_expr_gt_0.5 <- rownames(avg_expression_by_sex)[
  avg_expression_by_sex[, "Male"] > 0.5 &
    avg_expression_by_sex[, "Female"] > 0.5
]

cat(
  "Number of genes with average expression > 0.5 in BOTH Male and Female Schwann cells:",
  length(genes_expr_gt_0.5),
  "\n"
)

head(genes_expr_gt_0.5, 20)

if (length(genes_expr_gt_0.5) < 50) {
  warning("Very few genes passed the expression > 0.5 filter. GSEA may be underpowered.")
}


# ============================================================
# 5. RUN DIFFERENTIAL GENE EXPRESSION: MALE VS FEMALE SCHWANN
# ============================================================

print("Running DGE: Male vs Female True Schwann cells...")

schwann_sex_markers <- FindMarkers(
  object = true_schwann_subset,
  ident.1 = "Male",
  ident.2 = "Female",
  group.by = "Sex",
  features = genes_expr_gt_0.5,
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.05,
  only.pos = FALSE
)

schwann_sex_markers$Gene <- rownames(schwann_sex_markers)

cat("Number of genes returned by FindMarkers:",
    nrow(schwann_sex_markers), "\n")


# ============================================================
# 6. LABEL SIGNIFICANT GENES FOR VOLCANO PLOT
# ============================================================

schwann_sex_markers$Significance <- "Not Significant"

p_val_cutoff <- 0.05
log2fc_cutoff <- 0.5

# Positive avg_log2FC = higher in Male Schwann cells
# Negative avg_log2FC = higher in Female Schwann cells

schwann_sex_markers$Significance[
  schwann_sex_markers$p_val_adj < p_val_cutoff &
    schwann_sex_markers$avg_log2FC > log2fc_cutoff
] <- "Upregulated in Males"

schwann_sex_markers$Significance[
  schwann_sex_markers$p_val_adj < p_val_cutoff &
    schwann_sex_markers$avg_log2FC < -log2fc_cutoff
] <- "Upregulated in Females"

schwann_sex_markers$Significance <- factor(
  schwann_sex_markers$Significance,
  levels = c("Upregulated in Females", "Not Significant", "Upregulated in Males")
)

print(table(schwann_sex_markers$Significance))


# ============================================================
# 7. HANDLE ZERO ADJUSTED P-VALUES SAFELY
# ============================================================

smallest_nonzero_padj <- min(
  schwann_sex_markers$p_val_adj[schwann_sex_markers$p_val_adj > 0],
  na.rm = TRUE
)

schwann_sex_markers <- schwann_sex_markers %>%
  mutate(
    p_val_adj_safe = ifelse(p_val_adj == 0, smallest_nonzero_padj, p_val_adj),
    neg_log10_padj = -log10(p_val_adj_safe)
  )


# ============================================================
# 8. IDENTIFY TOP GENES TO LABEL
# ============================================================

top_schwann_genes <- schwann_sex_markers %>%
  filter(Significance != "Not Significant") %>%
  group_by(Significance) %>%
  slice_min(order_by = p_val_adj_safe, n = 15) %>%
  ungroup()


# ============================================================
# 9. CREATE OUTPUT FOLDER
# ============================================================

schwann_outdir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis/Schwann_DGE_GSEA_Figures"
dir.create(schwann_outdir, showWarnings = FALSE, recursive = TRUE)

write.csv(
  schwann_sex_markers,
  file.path(schwann_outdir, "Schwann_Male_vs_Female_DGE_ExpressionGT0.5.csv"),
  row.names = FALSE
)


# ============================================================
# 10. VOLCANO PLOT: MALE VS FEMALE SCHWANN CELLS
# ============================================================

schwann_volcano <- ggplot(
  schwann_sex_markers,
  aes(x = avg_log2FC, y = neg_log10_padj, color = Significance)
) +
  geom_point(alpha = 0.8, size = 1.5) +
  
  scale_color_manual(
    values = c(
      "Upregulated in Females" = "#8B0000",
      "Not Significant" = "grey80",
      "Upregulated in Males" = "#0072B2"
    ),
    labels = c(
      "Upregulated in Females" = "Schwann cells from females",
      "Not Significant" = "Not Significant",
      "Upregulated in Males" = "Schwann cells from males"
    )
  ) +
  
  geom_vline(
    xintercept = c(-log2fc_cutoff, log2fc_cutoff),
    linetype = "dashed",
    color = "black",
    alpha = 0.5
  ) +
  
  geom_hline(
    yintercept = -log10(p_val_cutoff),
    linetype = "dashed",
    color = "black",
    alpha = 0.5
  ) +
  
  geom_text_repel(
    data = top_schwann_genes,
    aes(label = Gene),
    size = 4,
    fontface = "italic",
    color = "black",
    box.padding = 0.5,
    point.padding = 0.4,
    max.overlaps = 30,
    min.segment.length = 0
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.03, 0.18))
  ) +
  
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  coord_cartesian(clip = "off") +
  
  ggtitle("Schwann Cell Differential Expression: Male vs. Female") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted P-Value)"
  ) +
  
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16,
      margin = margin(b = 15)
    ),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, face = "bold"),
    plot.margin = margin(t = 25, r = 25, b = 20, l = 20)
  )

print(schwann_volcano)

ggsave(
  file.path(schwann_outdir, "Figure_1_Schwann_Male_vs_Female_Volcano.pdf"),
  plot = schwann_volcano,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(schwann_outdir, "Figure_1_Schwann_Male_vs_Female_Volcano.png"),
  plot = schwann_volcano,
  width = 9,
  height = 7,
  dpi = 300
)


# ============================================================
# 11. PREPARE RANKED GENE LIST FOR GSEA
# ============================================================

print("Preparing ranked gene list for Schwann GSEA...")

schwann_markers_ranked <- schwann_sex_markers %>%
  filter(!is.na(avg_log2FC)) %>%
  filter(!is.na(p_val_adj_safe)) %>%
  mutate(
    rank_score = sign(avg_log2FC) * -log10(p_val_adj_safe)
  ) %>%
  filter(!is.na(rank_score), is.finite(rank_score)) %>%
  arrange(desc(rank_score))

schwann_gene_ranks <- schwann_markers_ranked$rank_score
names(schwann_gene_ranks) <- schwann_markers_ranked$Gene

schwann_gene_ranks <- schwann_gene_ranks[!duplicated(names(schwann_gene_ranks))]
schwann_gene_ranks <- sort(schwann_gene_ranks, decreasing = TRUE)

cat("Number of Schwann genes going into GSEA:",
    length(schwann_gene_ranks), "\n")

head(schwann_gene_ranks, 10)
tail(schwann_gene_ranks, 10)


# ============================================================
# 12. FETCH HALLMARK PATHWAYS
# ============================================================

print("Fetching Hallmark biological pathways...")

hallmark_sets <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

pathways_list <- split(
  x = hallmark_sets$gene_symbol,
  f = hallmark_sets$gs_name
)


# ============================================================
# 13. RUN fgsea
# ============================================================

print("Running fgsea for Schwann cells...")

set.seed(42)

schwann_fgsea_results <- fgsea(
  pathways = pathways_list,
  stats = schwann_gene_ranks,
  minSize = 10,
  maxSize = 500,
  eps = 0
)

schwann_fgsea_results <- schwann_fgsea_results %>%
  arrange(padj) %>%
  mutate(
    pathway_clean = gsub("HALLMARK_", "", pathway),
    pathway_clean = gsub("_", " ", pathway_clean),
    Direction = ifelse(
      NES > 0,
      "Enriched in Male Schwann Cells",
      "Enriched in Female Schwann Cells"
    ),
    neg_log10_FDR = -log10(padj)
  )

print(head(schwann_fgsea_results, 20))


# ============================================================
# 14. SAVE fgsea RESULTS SAFELY
# ============================================================
# fgsea has a list-column called leadingEdge, so collapse it before write.csv.

schwann_fgsea_results_save <- schwann_fgsea_results %>%
  as.data.frame() %>%
  mutate(
    leadingEdge = sapply(leadingEdge, paste, collapse = ";")
  )

write.csv(
  schwann_fgsea_results_save,
  file.path(schwann_outdir, "Schwann_Male_vs_Female_fgsea_ExpressionGT0.5_Results.csv"),
  row.names = FALSE
)


# ============================================================
# 15. GSEA FIGURE 1: NES DOTPLOT
# ============================================================

schwann_plot_gsea <- schwann_fgsea_results %>%
  filter(!is.na(NES), !is.na(padj)) %>%
  filter(padj < 0.25) %>%
  arrange(NES)

if (nrow(schwann_plot_gsea) < 5) {
  schwann_plot_gsea <- schwann_fgsea_results %>%
    filter(!is.na(NES), !is.na(padj)) %>%
    arrange(padj) %>%
    head(25) %>%
    arrange(NES)
}

schwann_plot_gsea$pathway_clean <- factor(
  schwann_plot_gsea$pathway_clean,
  levels = schwann_plot_gsea$pathway_clean
)

schwann_gsea_dotplot <- ggplot(
  schwann_plot_gsea,
  aes(x = NES, y = pathway_clean)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black",
    alpha = 0.5
  ) +
  geom_point(
    aes(size = size, color = neg_log10_FDR),
    alpha = 0.9
  ) +
  scale_color_gradient(
    low = "lightgrey",
    high = "darkred",
    name = "-log10(FDR)"
  ) +
  labs(
    title = "Hallmark GSEA: Male vs Female Schwann Cells",
    subtitle = "Positive NES = male-enriched; negative NES = female-enriched",
    x = "Normalized Enrichment Score",
    y = "Hallmark Pathway",
    size = "Gene Set Size"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 8)),
    plot.subtitle = element_text(hjust = 0.5, size = 11, margin = margin(b = 12)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

print(schwann_gsea_dotplot)

ggsave(
  file.path(schwann_outdir, "Figure_2_Schwann_GSEA_NES_Dotplot.pdf"),
  plot = schwann_gsea_dotplot,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(schwann_outdir, "Figure_2_Schwann_GSEA_NES_Dotplot.png"),
  plot = schwann_gsea_dotplot,
  width = 9,
  height = 7,
  dpi = 300
)


# ============================================================
# 16. GSEA FIGURE 2: TOP MALE VS FEMALE PATHWAY BARPLOT
# ============================================================

top_male_schwann <- schwann_fgsea_results %>%
  filter(NES > 0) %>%
  arrange(padj) %>%
  head(10)

top_female_schwann <- schwann_fgsea_results %>%
  filter(NES < 0) %>%
  arrange(padj) %>%
  head(10)

top_schwann_combined <- bind_rows(
  top_female_schwann,
  top_male_schwann
) %>%
  arrange(NES)

top_schwann_combined$pathway_clean <- factor(
  top_schwann_combined$pathway_clean,
  levels = top_schwann_combined$pathway_clean
)

schwann_gsea_barplot <- ggplot(
  top_schwann_combined,
  aes(x = pathway_clean, y = NES, fill = Direction)
) +
  geom_col(width = 0.75, color = "black", linewidth = 0.25) +
  coord_flip() +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = c(
    "Enriched in Female Schwann Cells" = "#8B0000",
    "Enriched in Male Schwann Cells" = "#0072B2"
  )) +
  labs(
    title = "Top Sex-Biased Pathways in PDAC Schwann Cells",
    x = "Hallmark Pathway",
    y = "Normalized Enrichment Score",
    fill = ""
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 11, face = "bold")
  )

print(schwann_gsea_barplot)

ggsave(
  file.path(schwann_outdir, "Figure_3_Schwann_GSEA_Barplot.pdf"),
  plot = schwann_gsea_barplot,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(schwann_outdir, "Figure_3_Schwann_GSEA_Barplot.png"),
  plot = schwann_gsea_barplot,
  width = 9,
  height = 7,
  dpi = 300
)


# ============================================================
# 17. GSEA FIGURE 3: TOP ENRICHMENT CURVES
# ============================================================

top_male_schwann_pathway <- schwann_fgsea_results %>%
  filter(NES > 0) %>%
  arrange(padj) %>%
  slice(1) %>%
  pull(pathway)

top_female_schwann_pathway <- schwann_fgsea_results %>%
  filter(NES < 0) %>%
  arrange(padj) %>%
  slice(1) %>%
  pull(pathway)

if (length(top_male_schwann_pathway) > 0) {
  
  male_schwann_curve <- plotEnrichment(
    pathways_list[[top_male_schwann_pathway]],
    schwann_gene_ranks
  ) +
    labs(
      title = paste(
        "Male-Enriched:",
        gsub("_", " ", gsub("HALLMARK_", "", top_male_schwann_pathway))
      ),
      x = "Ranked Genes",
      y = "Enrichment Score"
    ) +
    theme_classic(base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      axis.title = element_text(face = "bold")
    )
  
  print(male_schwann_curve)
  
  ggsave(
    file.path(schwann_outdir, "Figure_4A_Top_Male_Schwann_Enrichment_Curve.pdf"),
    plot = male_schwann_curve,
    width = 8,
    height = 5,
    dpi = 300
  )
}

if (length(top_female_schwann_pathway) > 0) {
  
  female_schwann_curve <- plotEnrichment(
    pathways_list[[top_female_schwann_pathway]],
    schwann_gene_ranks
  ) +
    labs(
      title = paste(
        "Female-Enriched:",
        gsub("_", " ", gsub("HALLMARK_", "", top_female_schwann_pathway))
      ),
      x = "Ranked Genes",
      y = "Enrichment Score"
    ) +
    theme_classic(base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      axis.title = element_text(face = "bold")
    )
  
  print(female_schwann_curve)
  
  ggsave(
    file.path(schwann_outdir, "Figure_4B_Top_Female_Schwann_Enrichment_Curve.pdf"),
    plot = female_schwann_curve,
    width = 8,
    height = 5,
    dpi = 300
  )
}


# ============================================================
# 18. GSEA FIGURE 4: LEADING-EDGE GENE HEATMAP
# ============================================================

leading_edge_schwann_genes <- schwann_fgsea_results %>%
  filter(pathway %in% c(top_male_schwann_pathway, top_female_schwann_pathway)) %>%
  pull(leadingEdge) %>%
  unlist() %>%
  unique()

leading_edge_schwann_genes <- intersect(
  leading_edge_schwann_genes,
  rownames(true_schwann_subset)
)

leading_edge_schwann_genes <- head(leading_edge_schwann_genes, 50)

if (length(leading_edge_schwann_genes) >= 2) {
  
  schwann_avg_exp <- AverageExpression(
    true_schwann_subset,
    features = leading_edge_schwann_genes,
    group.by = "Sex",
    assays = "RNA",
    slot = "data"
  )$RNA
  
  schwann_avg_exp_scaled <- t(scale(t(schwann_avg_exp)))
  schwann_avg_exp_scaled[is.na(schwann_avg_exp_scaled)] <- 0
  
  schwann_heatmap_df <- as.data.frame(schwann_avg_exp_scaled) %>%
    rownames_to_column("Gene") %>%
    pivot_longer(
      cols = -Gene,
      names_to = "Sex",
      values_to = "ScaledExpression"
    )
  
  schwann_heatmap_df$Sex <- factor(
    schwann_heatmap_df$Sex,
    levels = c("Female", "Male")
  )
  
  schwann_leading_edge_heatmap <- ggplot(
    schwann_heatmap_df,
    aes(x = Sex, y = Gene, fill = ScaledExpression)
  ) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low = "#8B0000",
      mid = "white",
      high = "#0072B2",
      midpoint = 0,
      name = "Scaled\nExpression"
    ) +
    labs(
      title = "Leading-Edge Genes from Top Sex-Biased Schwann Pathways",
      x = "Patient Sex",
      y = "Leading-Edge Gene"
    ) +
    theme_classic(base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
      axis.title = element_text(face = "bold", size = 13),
      axis.text.x = element_text(size = 12, face = "bold", color = "black"),
      axis.text.y = element_text(size = 8, face = "italic", color = "black"),
      legend.title = element_text(face = "bold")
    )
  
  print(schwann_leading_edge_heatmap)
  
  ggsave(
    file.path(schwann_outdir, "Figure_5_Schwann_Leading_Edge_Heatmap.pdf"),
    plot = schwann_leading_edge_heatmap,
    width = 6,
    height = 10,
    dpi = 300
  )
  
  ggsave(
    file.path(schwann_outdir, "Figure_5_Schwann_Leading_Edge_Heatmap.png"),
    plot = schwann_leading_edge_heatmap,
    width = 6,
    height = 10,
    dpi = 300
  )
}


# ============================================================
# 19. PRINT FINAL SUMMARY
# ============================================================

cat("\n======================================================\n")
cat("SCHWANN DGE + GSEA COMPLETE\n")
cat("Only genes with average expression > 0.5 in BOTH sexes were included.\n")
cat("Positive avg_log2FC = higher in Male Schwann cells.\n")
cat("Negative avg_log2FC = higher in Female Schwann cells.\n")
cat("Positive NES = pathway enriched in Male Schwann cells.\n")
cat("Negative NES = pathway enriched in Female Schwann cells.\n")
cat("Results saved to:\n")
cat(schwann_outdir, "\n")
cat("======================================================\n")

cat("\nTop Male-Upregulated Schwann Genes:\n")
print(
  schwann_sex_markers %>%
    filter(Significance == "Upregulated in Males") %>%
    arrange(desc(avg_log2FC)) %>%
    dplyr::select(Gene, avg_log2FC, p_val_adj, pct.1, pct.2) %>%
    head(20)
)

cat("\nTop Female-Upregulated Schwann Genes:\n")
print(
  schwann_sex_markers %>%
    filter(Significance == "Upregulated in Females") %>%
    arrange(avg_log2FC) %>%
    dplyr::select(Gene, avg_log2FC, p_val_adj, pct.1, pct.2) %>%
    head(20)
)

cat("\nTop Male-Enriched Schwann Pathways:\n")
print(
  schwann_fgsea_results %>%
    filter(NES > 0) %>%
    arrange(padj) %>%
    dplyr::select(pathway_clean, NES, pval, padj, size) %>%
    head(10)
)

cat("\nTop Female-Enriched Schwann Pathways:\n")
print(
  schwann_fgsea_results %>%
    filter(NES < 0) %>%
    arrange(padj) %>%
    dplyr::select(pathway_clean, NES, pval, padj, size) %>%
    head(10)
)

# Identify cells that co-express BOTH ERBB3 and PMP22
true_schwann_cells <- rownames(expr_data)[
  expr_data$ERBB3 > expr_threshold &
    expr_data$PMP22 > expr_threshold
]

print(paste("Number of ERBB3+ / PMP22+ true Schwann cells:", length(true_schwann_cells)))

true_schwann_label <- rep("Other", ncol(Merged_Joined))
names(true_schwann_label) <- colnames(Merged_Joined)

true_schwann_label[true_schwann_cells] <- "True Schwann"

Merged_Joined <- AddMetaData(
  Merged_Joined,
  metadata = true_schwann_label,
  col.name = "true_schwann_subset"
)

# Make it a factor so plotting order is consistent
Merged_Joined$true_schwann_subset <- factor(
  Merged_Joined$true_schwann_subset,
  levels = c("Other", "True Schwann")
)

Idents(Merged_Joined) <- "true_schwann_subset"

true_schwann_subset <- subset(
  Merged_Joined,
  subset = true_schwann_subset == "True Schwann"
)

print(true_schwann_subset)

saveRDS(
  true_schwann_subset,
  file = "/Users/sudharshiniram/Desktop/SRI Data Analysis/true_schwann_ERBB3_PMP22_subset.rds"
)

summary_data <- FetchData(
  Merged_Joined,
  vars = c("orig.ident", "true_schwann_subset")
)

schwann_summary <- summary_data %>%
  group_by(orig.ident) %>%
  summarise(
    Total_Cells = n(),
    True_Schwann_Cells = sum(true_schwann_subset == "True Schwann"),
    Percentage_True_Schwann = (True_Schwann_Cells / Total_Cells) * 100,
    .groups = "drop"
  )

patient_mapping <- c(
  "Sample_01" = "P01", "Sample_02" = "P02", "Sample_03" = "P03",
  "Sample_04" = "P04", "Sample_05" = "P05", "Sample_06" = "P06",
  "Sample_07" = "P07", "Sample_08" = "P08", "Sample_09" = "P09",
  "Sample_10" = "P10", "Sample_11" = "P11", "Sample_12" = "P12",
  "Sample_13" = "P13", "Sample_14" = "P14", "Sample_15" = "P15",
  "Sample_16" = "P16", "Sample_17" = "P17"
)

schwann_summary$PatientID <- patient_mapping[as.character(schwann_summary$orig.ident)]
schwann_summary$PatientID <- factor(
  schwann_summary$PatientID,
  levels = sprintf("P%02d", 1:17)
)

print("FINAL TRUE SCHWANN SUBSET: ERBB3+ / PMP22+ Co-Expression")
print("======================================================")

print(
  as.data.frame(
    schwann_summary[, c(
      "PatientID",
      "Total_Cells",
      "True_Schwann_Cells",
      "Percentage_True_Schwann"
    )]
  )
)

# Save summary table
write.csv(
  schwann_summary,
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/True_Schwann_ERBB3_PMP22_Per_Sample.csv",
  row.names = FALSE
)

# ============================================================
# 8. BAR CHART OF TRUE SCHWANN PERCENTAGE BY SAMPLE
# ============================================================
schwann_graph <- ggplot(
  schwann_summary,
  aes(x = PatientID, y = Percentage_True_Schwann)
) +
  geom_bar(stat = "identity", fill = "navy", color = "black", linewidth = 0.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  ggtitle("Proportion of ERBB3+ / PMP22+ True Schwann Cells") +
  labs(
    x = "Patient Sample",
    y = "Percentage of Total Cells (%)"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16,
      margin = margin(b = 15)
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 11,
      face = "bold",
      color = "black"
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 13,
      margin = margin(t = 10)
    ),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(
      face = "bold",
      size = 13,
      margin = margin(r = 10)
    ),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
  )

print(schwann_graph)



