library(dplyr)
library(ggplot2)
library(ggVennDiagram)

# ============================================================
# 1. MAKE SURE sex_markers HAS A GENE COLUMN
# ============================================================
sex_markers$Gene <- rownames(sex_markers)

# ============================================================
# 2. SET LOG FOLD-CHANGE THRESHOLD
# ============================================================
logfc_cutoff <- 0.5

# ============================================================
# 3. DEFINE MALE-UP AND FEMALE-UP GENES WITHOUT SIGNIFICANCE COLUMN
# ============================================================
# Since FindMarkers used ident.1 = "Male" and ident.2 = "Female":
# avg_log2FC > 0 means higher in males
# avg_log2FC < 0 means higher in females

male_up_genes <- sex_markers %>%
  filter(avg_log2FC > logfc_cutoff) %>%
  pull(Gene) %>%
  unique()

female_up_genes <- sex_markers %>%
  filter(avg_log2FC < -logfc_cutoff) %>%
  pull(Gene) %>%
  unique()

# ============================================================
# 4. CALCULATE VENN COUNTS WITHOUT DOUBLE COUNTING
# ============================================================
male_only_genes <- setdiff(male_up_genes, female_up_genes)
female_only_genes <- setdiff(female_up_genes, male_up_genes)
both_genes <- intersect(male_up_genes, female_up_genes)

total_unique_upregulated_genes <- length(union(male_up_genes, female_up_genes))

cat("======================================================\n")
cat("Male-only upregulated genes:", length(male_only_genes), "\n")
cat("Both male and female upregulated genes:", length(both_genes), "\n")
cat("Female-only upregulated genes:", length(female_only_genes), "\n")
cat("Total unique upregulated genes:", total_unique_upregulated_genes, "\n")
cat("======================================================\n")

cat("\nMale-only genes:\n")
print(male_only_genes)

cat("\nBoth male and female genes:\n")
print(both_genes)

cat("\nFemale-only genes:\n")
print(female_only_genes)

# ============================================================
# 5. MAKE VENN DIAGRAM
# ============================================================
venn_list <- list(
  "Male Upregulated" = male_up_genes,
  "Female Upregulated" = female_up_genes
)

venn_plot <- ggVennDiagram(
  venn_list,
  label = "count",
  label_alpha = 0
) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ggtitle("Male vs Female Upregulated Genes | log2FC > 0.5") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    legend.position = "none"
  )

print(venn_plot)

# ============================================================
# 6. SAVE VENN DIAGRAM
# ============================================================
ggsave(
  filename = "/Users/sudharshiniram/Desktop/SRI Data Analysis/Male_Female_Upregulated_Genes_Venn_logFC0.5.png",
  plot = venn_plot,
  width = 7,
  height = 6,
  dpi = 300
)

# ============================================================
# 7. SAVE GENE LISTS
# ============================================================
venn_gene_summary <- data.frame(
  Category = c(
    rep("Male Only Upregulated", length(male_only_genes)),
    rep("Both Male and Female Upregulated", length(both_genes)),
    rep("Female Only Upregulated", length(female_only_genes))
  ),
  Gene = c(
    male_only_genes,
    both_genes,
    female_only_genes
  )
)

write.csv(
  venn_gene_summary,
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Male_Female_Upregulated_Genes_Venn_GeneLists_logFC0.5.csv",
  row.names = FALSE
)


library(Seurat)
library(dplyr)
library(ggplot2)
library(ggVennDiagram)

# ============================================================
# 1. CHECK SEX COLUMN
# ============================================================
print(colnames(fibroblast_subset@meta.data))
print(table(fibroblast_subset$Sex))

# ============================================================
# 2. SPLIT FIBROBLASTS BY SEX
# ============================================================
male_fibroblasts <- subset(fibroblast_subset, subset = Sex == "Male")
female_fibroblasts <- subset(fibroblast_subset, subset = Sex == "Female")

cat("Male fibroblast cells:", ncol(male_fibroblasts), "\n")
cat("Female fibroblast cells:", ncol(female_fibroblasts), "\n")

# ============================================================
# 3. CALCULATE AVERAGE EXPRESSION PER GENE IN EACH SEX
# ============================================================
default_assay <- DefaultAssay(fibroblast_subset)

male_avg <- AverageExpression(
  male_fibroblasts,
  assays = default_assay,
  layer = "data"
)[[default_assay]]

female_avg <- AverageExpression(
  female_fibroblasts,
  assays = default_assay,
  layer = "data"
)[[default_assay]]

male_avg_df <- data.frame(
  Gene = rownames(male_avg),
  Male_Avg_Expression = male_avg[, 1]
)

female_avg_df <- data.frame(
  Gene = rownames(female_avg),
  Female_Avg_Expression = female_avg[, 1]
)

avg_expr_df <- inner_join(male_avg_df, female_avg_df, by = "Gene")

# ============================================================
# 4. DEFINE GENES EXPRESSED ABOVE 0.5
# ============================================================
expr_cutoff <- 0.5

male_up_genes <- avg_expr_df %>%
  filter(Male_Avg_Expression > expr_cutoff) %>%
  pull(Gene) %>%
  unique()

female_up_genes <- avg_expr_df %>%
  filter(Female_Avg_Expression > expr_cutoff) %>%
  pull(Gene) %>%
  unique()

# ============================================================
# 5. CALCULATE VENN COUNTS WITHOUT DOUBLE COUNTING
# ============================================================
male_only_genes <- setdiff(male_up_genes, female_up_genes)
female_only_genes <- setdiff(female_up_genes, male_up_genes)
both_genes <- intersect(male_up_genes, female_up_genes)

total_unique_upregulated_genes <- length(union(male_up_genes, female_up_genes))

cat("======================================================\n")
cat("Male-only genes expressed above 0.5:", length(male_only_genes), "\n")
cat("Both male and female genes expressed above 0.5:", length(both_genes), "\n")
cat("Female-only genes expressed above 0.5:", length(female_only_genes), "\n")
cat("Total unique genes expressed above 0.5:", total_unique_upregulated_genes, "\n")
cat("======================================================\n")

cat("\nMale-only genes:\n")
print(male_only_genes)

cat("\nBoth male and female genes:\n")
print(both_genes)

cat("\nFemale-only genes:\n")
print(female_only_genes)

# ============================================================
# 6. MAKE VENN DIAGRAM
# ============================================================
venn_list <- list(
  "Male Fibroblasts" = male_up_genes,
  "Female Fibroblasts" = female_up_genes
)

venn_plot <- ggVennDiagram(
  venn_list,
  label = "count",
  label_alpha = 0
) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ggtitle("Genes Expressed Above 0.5 in Male and Female Fibroblasts") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15
    ),
    legend.position = "none"
  )

print(venn_plot)

# ============================================================
# 7. SAVE VENN DIAGRAM
# ============================================================
ggsave(
  filename = "/Users/sudharshiniram/Desktop/SRI Data Analysis/Male_Female_Fibroblast_Expression_Venn.png",
  plot = venn_plot,
  width = 7,
  height = 6,
  dpi = 300
)

# ============================================================
# 8. SAVE GENE LISTS
# ============================================================
venn_gene_summary <- data.frame(
  Category = c(
    rep("Male Only Above 0.5", length(male_only_genes)),
    rep("Both Male and Female Above 0.5", length(both_genes)),
    rep("Female Only Above 0.5", length(female_only_genes))
  ),
  Gene = c(
    male_only_genes,
    both_genes,
    female_only_genes
  )
)

write.csv(
  venn_gene_summary,
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Male_Female_Fibroblast_Expression_Venn_GeneLists.csv",
  row.names = FALSE
)

# ============================================================
# 9. SAVE FULL AVERAGE EXPRESSION TABLE TOO
# ============================================================
write.csv(
  avg_expr_df,
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Male_Female_Fibroblast_AverageExpression.csv",
  row.names = FALSE
)

# ============================================================
# ADD SEX TO Merged_Joined USING orig.ident
# ============================================================
sample_sex_map <- fibroblast_subset@meta.data %>%
  select(orig.ident, Sex) %>%
  distinct()

print(sample_sex_map)

Merged_Joined@meta.data <- Merged_Joined@meta.data %>%
  left_join(sample_sex_map, by = "orig.ident")

print(table(Merged_Joined$Sex, useNA = "ifany"))