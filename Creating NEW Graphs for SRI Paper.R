

### 1B: Dotplot of Overall Umap and all Genes
library(Seurat)
library(ggplot2)
library(dplyr)

DefaultAssay(Merged_Joined) <- "RNA"

Merged_Joined$master_celltype_revised <- case_when(
  as.character(Merged_Joined$seurat_clusters) == "0" ~ "Cancer 1",
  as.character(Merged_Joined$seurat_clusters) == "9" ~ "Cancer 2",
  as.character(Merged_Joined$seurat_clusters) == "2" ~ "Ductal",
  as.character(Merged_Joined$seurat_clusters) == "1" ~ "Fibroblasts",
  as.character(Merged_Joined$seurat_clusters) == "10" ~ "Smooth muscle",
  as.character(Merged_Joined$seurat_clusters) == "3" ~ "Lymphoid",
  as.character(Merged_Joined$seurat_clusters) == "8" ~ "Lymphoid",
  as.character(Merged_Joined$seurat_clusters) == "4" ~ "Myeloid",
  as.character(Merged_Joined$seurat_clusters) == "7" ~ "Endothelial",
  as.character(Merged_Joined$seurat_clusters) == "5" ~ "Acinar",
  as.character(Merged_Joined$seurat_clusters) == "6" ~ "Endocrine",
  as.character(Merged_Joined$seurat_clusters) %in% c("11", "12", "13", "14") ~ "Other",
  TRUE ~ "Other"
)

celltype_order <- c(
  "Cancer 1",
  "Cancer 2",
  "Ductal",
  "Fibroblasts",
  "Smooth muscle",
  "Lymphoid",
  "Myeloid",
  "Endothelial",
  "Acinar",
  "Endocrine",
  "Other"
)

Merged_Joined$master_celltype_revised <- factor(
  Merged_Joined$master_celltype_revised,
  levels = celltype_order
)

gene_order <- c(
  "EPCAM", "KRT19", "CEACAM6", "MUC1",
  "CFTR", "SLC4A4", "SCTR",
  "COL1A1", "DCN", "PDGFRA", "FAP",
  "RGS5", "PDGFRB", "NOTCH3",
  "MS4A1", "BLK", "CD79A", "CD19",
  "CD3D", "CD8A", "NKG7", "IL7R",
  "CD68", "CD14", "C1QA", "IRF8", "CLEC9A",
  "PECAM1", "VWF", "PLVAP",
  "PRSS1", "PNLIP", "AMY2A",
  "CHGA", "INS", "GCG",
  "SOX10", "S100B", "PLP1"
)

gene_order <- unique(gene_order)

gene_order_present <- gene_order[
  gene_order %in% rownames(Merged_Joined)
]

missing_genes <- setdiff(gene_order, gene_order_present)

print("Missing genes:")
print(missing_genes)

dotplot_base <- DotPlot(
  Merged_Joined,
  features = gene_order_present,
  group.by = "master_celltype_revised"
)

dotplot_data <- dotplot_base$data

dotplot_data$id <- factor(
  dotplot_data$id,
  levels = celltype_order
)

dotplot_data$features.plot <- factor(
  dotplot_data$features.plot,
  levels = gene_order_present
)

dotplot_data$avg.exp.scaled_plot <- ifelse(
  dotplot_data$id == "Other",
  NA,
  dotplot_data$avg.exp.scaled
)

master_dotplot_revised_ordered <- ggplot(
  dotplot_data,
  aes(
    x = id,
    y = features.plot,
    size = pct.exp,
    color = avg.exp.scaled_plot
  )
) +
  geom_point() +
  scale_color_gradient(
    low = "#E6E6E6",
    high = "navy",
    na.value = "grey70",
    name = "Average Expression"
  ) +
  scale_size(
    range = c(0, 6),
    name = "Percent Expressed"
  ) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 11,
      color = "black"
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 13,
      margin = margin(t = 10)
    ),
    axis.text.y = element_text(
      face = "italic",
      size = 10,
      color = "black"
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 13,
      margin = margin(r = 10)
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16,
      margin = margin(b = 15)
    ),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("Master Cell Type Marker Expression") +
  labs(
    x = "Identified Cell Type",
    y = "Marker Gene"
  )

print(master_dotplot_revised_ordered)




#### 1D: Cell Type Percentage for all Samples

library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)

plot_data <- Merged_Joined@meta.data

plot_data$CellType <- case_when(
  as.character(plot_data$seurat_clusters) == "0" ~ "Cancer 1",
  as.character(plot_data$seurat_clusters) == "1" ~ "Fibroblast",
  as.character(plot_data$seurat_clusters) == "2" ~ "Ductal",
  as.character(plot_data$seurat_clusters) == "3" ~ "B Cells",
  as.character(plot_data$seurat_clusters) == "4" ~ "Myeloid",
  as.character(plot_data$seurat_clusters) == "5" ~ "Acinar",
  as.character(plot_data$seurat_clusters) == "6" ~ "Endocrine",
  as.character(plot_data$seurat_clusters) == "7" ~ "Endothelial",
  as.character(plot_data$seurat_clusters) == "8" ~ "Lymphoid",
  as.character(plot_data$seurat_clusters) == "9" ~ "Cancer 2",
  as.character(plot_data$seurat_clusters) == "10" ~ "Smooth muscle",
  as.character(plot_data$seurat_clusters) %in% c("11", "12", "13", "14") ~ "Other",
  TRUE ~ "Other"
)

celltype_order <- c(
  "Cancer 1",
  "Fibroblast",
  "Ductal",
  "B Cells",
  "Myeloid",
  "Acinar",
  "Endocrine",
  "Endothelial",
  "Lymphoid",
  "Cancer 2",
  "Smooth muscle",
  "Other"
)

plot_data$CellType <- factor(
  plot_data$CellType,
  levels = celltype_order
)

patient_mapping <- c(
  "Sample_01" = "P01",
  "Sample_02" = "P02",
  "Sample_03" = "P03",
  "Sample_04" = "P04",
  "Sample_05" = "P05",
  "Sample_06" = "P06",
  "Sample_07" = "P07",
  "Sample_08" = "P08",
  "Sample_09" = "P09",
  "Sample_10" = "P10",
  "Sample_11" = "P11",
  "Sample_12" = "P12",
  "Sample_13" = "P13",
  "Sample_14" = "P14",
  "Sample_15" = "P15",
  "Sample_16" = "P16",
  "Sample_17" = "P17"
)

plot_data$PatientID <- patient_mapping[as.character(plot_data$orig.ident)]

plot_data$PatientID <- factor(
  plot_data$PatientID,
  levels = sprintf("P%02d", 1:17)
)

main_celltypes <- setdiff(celltype_order, "Other")

celltype_colors <- hue_pal()(length(main_celltypes))
names(celltype_colors) <- main_celltypes

celltype_colors <- c(
  celltype_colors,
  "Other" = "grey70"
)

proportion_plot_revised <- ggplot(
  plot_data,
  aes(x = PatientID, fill = CellType)
) +
  geom_bar(
    position = "fill",
    color = "white",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = celltype_colors,
    breaks = celltype_order,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.01),
    expand = c(0, 0)
  ) +
  ggtitle("Cell Type Proportions per Patient Sample") +
  labs(
    x = "Patient Sample",
    y = "Proportion of Cells",
    fill = ""
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 10)),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 10)),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm")
  ) +
  guides(
    fill = guide_legend(nrow = 3, byrow = TRUE)
  )

print(proportion_plot_revised)



### 2A: Fibroblast Subset UMAP
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)

fibroblast_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS"
)

fibroblast_subset$CAF_cluster_number <- as.character(fibroblast_subset$seurat_clusters)

fibroblast_subset$CAF_cluster_number <- factor(
  fibroblast_subset$CAF_cluster_number,
  levels = c("0", "1", "2", "3", "4")
)

caf_legend_labels <- c(
  "0" = "0 - Fibroblasts",
  "1" = "1 - myCAF",
  "2" = "2 - iCAF",
  "3" = "3 - Pancreatic stellate",
  "4" = "4 - apCAF"
)

caf_colors <- hue_pal()(5)
names(caf_colors) <- c("0", "1", "2", "3", "4")

fibroblast_umap <- DimPlot(
  fibroblast_subset,
  reduction = "sub_umap",
  group.by = "CAF_cluster_number",
  label = TRUE,
  repel = TRUE,
  label.size = 6,
  label.color = "black",
  label.box = TRUE,
  pt.size = 0.4
) +
  scale_color_manual(
    values = caf_colors,
    breaks = names(caf_legend_labels),
    labels = caf_legend_labels
  ) +
  labs(
    x = "umap_1",
    y = "umap_2"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )

print(fibroblast_umap)

### 2B: Fibroblast Dotplot

library(Seurat)
library(ggplot2)
library(dplyr)

fibroblast_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS"
)

DefaultAssay(fibroblast_subset) <- "RNA"

fibroblast_subset$CAF_type <- case_when(
  as.character(fibroblast_subset$seurat_clusters) == "0" ~ "Fibroblasts",
  as.character(fibroblast_subset$seurat_clusters) == "1" ~ "myCAF",
  as.character(fibroblast_subset$seurat_clusters) == "2" ~ "iCAF",
  as.character(fibroblast_subset$seurat_clusters) == "3" ~ "Pancreatic stellate",
  as.character(fibroblast_subset$seurat_clusters) == "4" ~ "apCAF",
  TRUE ~ "Other"
)

fibroblast_subset$CAF_type <- factor(
  fibroblast_subset$CAF_type,
  levels = c(
    "myCAF",
    "Pancreatic stellate",
    "Fibroblasts",
    "iCAF",
    "apCAF"
  )
)

caf_markers_clean <- c(
  "MMP11", "POSTN", "COL1A1", "COL11A1",
  "IGFBP5", "CCL2", "APOE", "C2orf40",
  "ACTA2",
  "COL4A1", "COL4A2",
  "COL1A2", "PDGFRA", "C3"
)

caf_markers_present <- caf_markers_clean[
  caf_markers_clean %in% rownames(fibroblast_subset)
]

missing_genes <- setdiff(caf_markers_clean, caf_markers_present)

print("Genes included in DotPlot:")
print(caf_markers_present)

print("Missing genes:")
print(missing_genes)

caf_dotplot_ordered <- DotPlot(
  fibroblast_subset,
  features = caf_markers_present,
  group.by = "CAF_type"
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "navy") +
  ggtitle("CAF Subtype Marker Expression") +
  labs(
    x = "Marker genes",
    y = "CAF subtype"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 11, color = "black"),
    axis.text.y = element_text(face = "bold", size = 11, color = "black"),
    axis.title = element_text(face = "bold", size = 13),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  )

print(caf_dotplot_ordered)




### 2C: Fibroblast UMAP Split by Sex

library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)

fibroblast_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS"
)

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

fibroblast_subset$Sex <- unname(
  sex_mapping[as.character(fibroblast_subset$orig.ident)]
)

fibroblast_subset$CAF_cluster <- as.character(fibroblast_subset$seurat_clusters)

fibroblast_subset$CAF_cluster <- factor(
  fibroblast_subset$CAF_cluster,
  levels = c("0", "1", "2", "3", "4")
)

caf_legend_labels <- c(
  "0" = "0 - Fibroblasts",
  "1" = "1 - myCAF",
  "2" = "2 - iCAF",
  "3" = "3 - Pancreatic stellate",
  "4" = "4 - apCAF"
)

caf_colors <- hue_pal()(5)
names(caf_colors) <- c("0", "1", "2", "3", "4")

female_fibros <- subset(
  fibroblast_subset,
  subset = Sex == "Female"
)

male_fibros <- subset(
  fibroblast_subset,
  subset = Sex == "Male"
)

female_umap <- DimPlot(
  female_fibros,
  reduction = "sub_umap",
  group.by = "CAF_cluster",
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 0.5
) +
  scale_color_manual(
    values = caf_colors,
    breaks = names(caf_legend_labels),
    labels = caf_legend_labels
  ) +
  ggtitle("Female PDAC Fibroblasts") +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = "none"
  )

male_umap <- DimPlot(
  male_fibros,
  reduction = "sub_umap",
  group.by = "CAF_cluster",
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 0.5
) +
  scale_color_manual(
    values = caf_colors,
    breaks = names(caf_legend_labels),
    labels = caf_legend_labels
  ) +
  ggtitle("Male PDAC Fibroblasts") +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = "right"
  )

combined_fibro_plot <- female_umap + male_umap + plot_layout(guides = "collect")

print(combined_fibro_plot)


### 2D: Part 2 Fibroblast Cell Type Percentage with male versus female Samples

library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)

fibroblast_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS"
)

fibroblast_subset$CAF_type <- case_when(
  as.character(fibroblast_subset$seurat_clusters) == "0" ~ "Fibroblasts",
  as.character(fibroblast_subset$seurat_clusters) == "1" ~ "myCAF",
  as.character(fibroblast_subset$seurat_clusters) == "2" ~ "iCAF",
  as.character(fibroblast_subset$seurat_clusters) == "3" ~ "Pancreatic Stellate",
  as.character(fibroblast_subset$seurat_clusters) == "4" ~ "apCAF",
  TRUE ~ "Other"
)

caf_order <- c(
  "Fibroblasts",
  "myCAF",
  "iCAF",
  "Pancreatic Stellate",
  "apCAF"
)

fibroblast_subset$CAF_type <- factor(
  fibroblast_subset$CAF_type,
  levels = caf_order
)

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

patient_mapping <- c(
  "Sample_01" = "P01",
  "Sample_02" = "P02",
  "Sample_03" = "P03",
  "Sample_04" = "P04",
  "Sample_05" = "P05",
  "Sample_06" = "P06",
  "Sample_07" = "P07",
  "Sample_08" = "P08",
  "Sample_09" = "P09",
  "Sample_10" = "P10",
  "Sample_11" = "P11",
  "Sample_12" = "P12",
  "Sample_13" = "P13",
  "Sample_14" = "P14",
  "Sample_15" = "P15",
  "Sample_16" = "P16",
  "Sample_17" = "P17"
)

plot_data <- fibroblast_subset@meta.data

plot_data$CellType <- fibroblast_subset$CAF_type
plot_data$Sex <- unname(sex_mapping[as.character(plot_data$orig.ident)])
plot_data$PatientID <- unname(patient_mapping[as.character(plot_data$orig.ident)])

plot_data$CellType <- factor(
  plot_data$CellType,
  levels = caf_order
)

plot_data$PatientID <- factor(
  plot_data$PatientID,
  levels = sprintf("P%02d", 1:17)
)

plot_data <- plot_data %>%
  filter(!is.na(Sex), !is.na(PatientID), !is.na(CellType))

caf_colors <- hue_pal()(length(caf_order))
names(caf_colors) <- caf_order

plot_data_female <- plot_data %>%
  filter(Sex == "Female")

plot_data_male <- plot_data %>%
  filter(Sex == "Male")

female_caf_plot <- ggplot(
  plot_data_female,
  aes(x = PatientID, fill = CellType)
) +
  geom_bar(
    position = "fill",
    color = "white",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = caf_colors,
    breaks = caf_order,
    drop = FALSE
  ) +
  scale_x_discrete(drop = TRUE) +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.01),
    expand = c(0, 0)
  ) +
  ggtitle("Female Samples") +
  labs(
    x = "Patient Sample",
    y = "Proportion",
    fill = "Subpopulation"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  )

male_caf_plot <- ggplot(
  plot_data_male,
  aes(x = PatientID, fill = CellType)
) +
  geom_bar(
    position = "fill",
    color = "white",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = caf_colors,
    breaks = caf_order,
    drop = FALSE
  ) +
  scale_x_discrete(drop = TRUE) +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.01),
    expand = c(0, 0)
  ) +
  ggtitle("Male Samples") +
  labs(
    x = "Patient Sample",
    y = "Proportion",
    fill = "Subpopulation"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  )

caf_sex_comparison_plot <- female_caf_plot + male_caf_plot +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "CAF Subpopulation Proportions by Sex",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 18,
        margin = margin(b = 15)
      )
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )
print(caf_sex_comparison_plot)



### 3B: 1 male vs female GSEA Fibroblast
# NO AVERAGE EXPRESSION > 0.5 FILTER

library(Seurat)
library(dplyr)
library(ggplot2)
library(msigdbr)
library(fgsea)
library(tibble)
library(tidyr)
library(patchwork)

fibroblast_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS"
)

DefaultAssay(fibroblast_subset) <- "RNA"

# Optional but helpful for Seurat v5 if layers are split
fibroblast_subset <- tryCatch(
  JoinLayers(fibroblast_subset),
  error = function(e) fibroblast_subset
)


sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

fibroblast_subset$Sex <- unname(
  sex_mapping[as.character(fibroblast_subset$orig.ident)]
)

print(table(fibroblast_subset$orig.ident, fibroblast_subset$Sex, useNA = "ifany"))
print(table(fibroblast_subset$Sex, useNA = "ifany"))

fibroblast_subset <- subset(
  fibroblast_subset,
  subset = Sex %in% c("Male", "Female")
)

print(table(fibroblast_subset$Sex))

sex_markers <- FindMarkers(
  object = fibroblast_subset,
  ident.1 = "Male",
  ident.2 = "Female",
  group.by = "Sex",
  logfc.threshold = 0,
  min.pct = 0,
  only.pos = FALSE
)

sex_markers$Gene <- rownames(sex_markers)

cat(
  "Number of genes returned by FindMarkers without average expression > 0.5 filter:",
  nrow(sex_markers),
  "\n"
)


print("Fetching Hallmark biological pathways...")

hallmark_sets <- tryCatch(
  msigdbr(
    species = "Homo sapiens",
    category = "H"
  ),
  error = function(e) {
    msigdbr(
      species = "Homo sapiens",
      collection = "H"
    )
  }
)

pathways_list <- split(
  x = hallmark_sets$gene_symbol,
  f = hallmark_sets$gs_name
)

sex_markers_clean <- sex_markers %>%
  filter(!is.na(p_val_adj)) %>%
  filter(!is.na(avg_log2FC))

nonzero_padj <- sex_markers_clean$p_val_adj[
  sex_markers_clean$p_val_adj > 0
]

smallest_nonzero_padj <- ifelse(
  length(nonzero_padj) > 0,
  min(nonzero_padj, na.rm = TRUE),
  .Machine$double.xmin
)

sex_markers_clean <- sex_markers_clean %>%
  mutate(
    p_val_adj_safe = ifelse(
      p_val_adj == 0,
      smallest_nonzero_padj,
      p_val_adj
    ),
    rank_score = sign(avg_log2FC) * -log10(p_val_adj_safe)
  ) %>%
  filter(!is.na(rank_score), is.finite(rank_score)) %>%
  arrange(desc(abs(rank_score))) %>%
  group_by(Gene) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(desc(rank_score))

gene_ranks <- sex_markers_clean$rank_score
names(gene_ranks) <- sex_markers_clean$Gene

gene_ranks <- sort(gene_ranks, decreasing = TRUE)

cat(
  "Number of genes going into fgsea without average expression > 0.5 filter:",
  length(gene_ranks),
  "\n"
)

head(gene_ranks, 10)
tail(gene_ranks, 10)

print("Running fgsea...")

set.seed(42)

fgsea_results <- fgsea(
  pathways = pathways_list,
  stats = gene_ranks,
  minSize = 10,
  maxSize = 500,
  eps = 0
)

nonzero_padj_gsea <- fgsea_results$padj[
  fgsea_results$padj > 0
]

smallest_nonzero_gsea_padj <- ifelse(
  length(nonzero_padj_gsea) > 0,
  min(nonzero_padj_gsea, na.rm = TRUE),
  .Machine$double.xmin
)

fgsea_results <- fgsea_results %>%
  arrange(padj) %>%
  mutate(
    pathway_clean = gsub("HALLMARK_", "", pathway),
    pathway_clean = gsub("_", " ", pathway_clean),
    Direction = ifelse(
      NES > 0,
      "Enriched in Male Fibroblasts",
      "Enriched in Female Fibroblasts"
    ),
    padj_safe = ifelse(
      padj == 0,
      smallest_nonzero_gsea_padj,
      padj
    ),
    neg_log10_FDR = -log10(padj_safe)
  )

print(head(fgsea_results, 20))

gsea_outdir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis/GSEA_Figures_NoExpressionFilter"

dir.create(
  gsea_outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

fgsea_results_save <- fgsea_results %>%
  as.data.frame() %>%
  mutate(
    leadingEdge = sapply(leadingEdge, paste, collapse = ";")
  )

write.csv(
  fgsea_results_save,
  file.path(
    gsea_outdir,
    "Fibroblast_Male_vs_Female_fgsea_NoExpressionFilter_Results.csv"
  ),
  row.names = FALSE
)

plot_gsea <- fgsea_results %>%
  filter(!is.na(NES), !is.na(padj)) %>%
  filter(padj < 0.25) %>%
  arrange(NES)

if (nrow(plot_gsea) < 5) {
  plot_gsea <- fgsea_results %>%
    filter(!is.na(NES), !is.na(padj)) %>%
    arrange(padj) %>%
    head(25) %>%
    arrange(NES)
}

plot_gsea$pathway_clean <- factor(
  plot_gsea$pathway_clean,
  levels = plot_gsea$pathway_clean
)

gsea_nes_dotplot <- ggplot(
  plot_gsea,
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
    title = "Hallmark GSEA: Male vs Female PDAC Fibroblasts",
    subtitle = "Positive NES = male-enriched; negative NES = female-enriched",
    x = "Normalized Enrichment Score",
    y = "Hallmark Pathway",
    size = "Gene Set Size"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16,
      margin = margin(b = 8)
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 11,
      margin = margin(b = 12)
    ),
    axis.title = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

print(gsea_nes_dotplot)

ggsave(
  file.path(gsea_outdir, "Figure_1_fgsea_NES_Dotplot_NoExpressionFilter.pdf"),
  plot = gsea_nes_dotplot,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(gsea_outdir, "Figure_1_fgsea_NES_Dotplot_NoExpressionFilter.png"),
  plot = gsea_nes_dotplot,
  width = 9,
  height = 7,
  dpi = 300
)

sig_gsea <- fgsea_results %>%
  filter(!is.na(NES), !is.na(padj))

sig_gsea_for_bar <- sig_gsea %>%
  filter(padj < 0.25)

if (nrow(sig_gsea_for_bar) < 10) {
  sig_gsea_for_bar <- sig_gsea
}

top_male_pathways <- sig_gsea_for_bar %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 10)

top_female_pathways <- sig_gsea_for_bar %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  slice_head(n = 10)

top_sex_biased <- bind_rows(
  top_female_pathways,
  top_male_pathways
) %>%
  arrange(NES)

top_sex_biased$pathway_clean <- factor(
  top_sex_biased$pathway_clean,
  levels = top_sex_biased$pathway_clean
)

sex_biased_barplot <- ggplot(
  top_sex_biased,
  aes(x = NES, y = pathway_clean, fill = Direction)
) +
  geom_col(
    color = "black",
    linewidth = 0.25
  ) +
  geom_vline(
    xintercept = 0,
    color = "black",
    linewidth = 0.4
  ) +
  scale_fill_manual(
    values = c(
      "Enriched in Female Fibroblasts" = "darkred",
      "Enriched in Male Fibroblasts" = "#0072B2"
    )
  ) +
  labs(
    title = "Top Sex-Biased Pathways in PDAC Fibroblasts",
    x = "Normalized Enrichment Score",
    y = "Hallmark Pathway",
    fill = ""
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
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 11, face = "bold")
  )

print(sex_biased_barplot)

ggsave(
  file.path(gsea_outdir, "Figure_2_Top_Sex_Biased_Pathways_NoExpressionFilter.pdf"),
  plot = sex_biased_barplot,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(gsea_outdir, "Figure_2_Top_Sex_Biased_Pathways_NoExpressionFilter.png"),
  plot = sex_biased_barplot,
  width = 10,
  height = 7,
  dpi = 300
)


#Recreating 3B: Part 1 with the spermatogensis
library(dplyr)
library(ggplot2)
library(scales)

pathways_to_plot <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_HYPOXIA",
  "HALLMARK_SPERMATOGENESIS",
  "HALLMARK_PANCREAS_BETA_CELLS",
  "HALLMARK_APICAL_SURFACE",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING"
)

plot_df <- gsea_df %>%
  filter(ID %in% pathways_to_plot) %>%
  mutate(
    pathway = gsub("HALLMARK_", "", ID),
    pathway = gsub("_", " ", pathway),
    pathway = factor(
      pathway,
      levels = rev(gsub("_", " ", gsub("HALLMARK_", "", pathways_to_plot)))
    )
  )

ggplot(plot_df, aes(x = NES, y = pathway)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.8) +
  geom_point(aes(size = setSize, color = p.adjust)) +
  scale_color_gradient(
    low = "red",
    high = "blue",
    name = "Adjusted P-value",
    breaks = c(0.05, 0.10, 0.15, 0.20),
    limits = c(0.025, 0.225),
    oob = scales::squish
  ) +
  scale_size_continuous(
    name = "Gene Set Size",
    range = c(2, 9),
    breaks = c(50, 100, 150)
  ) +
  labs(
    title = "Hallmark GSEA: Male vs Female PDAC Fibroblasts",
    x = "Normalized Enrichment Score",
    y = "Hallmark Pathway"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12)
  )

plot_df[order(plot_df$NES, decreasing = TRUE), c("ID", "NES", "p.adjust", "setSize")]



### 3F: Venn Diagram

# ============================================================
# VENN-STYLE DIAGRAM:
# FEMALE-UPREGULATED, MALE-UPREGULATED, AND SHARED EXPRESSED GENES
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)

if (!requireNamespace("ggforce", quietly = TRUE)) {
  install.packages("ggforce")
}
library(ggforce)

fibroblast_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS"
)

DefaultAssay(fibroblast_subset) <- "RNA"

fibroblast_subset <- tryCatch(
  JoinLayers(fibroblast_subset),
  error = function(e) fibroblast_subset
)

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

fibroblast_subset$Sex <- unname(
  sex_mapping[as.character(fibroblast_subset$orig.ident)]
)

fibroblast_subset <- subset(
  fibroblast_subset,
  subset = Sex %in% c("Male", "Female")
)
print(table(fibroblast_subset$Sex))
sex_markers <- FindMarkers(
  object = fibroblast_subset,
  ident.1 = "Male",
  ident.2 = "Female",
  group.by = "Sex",
  logfc.threshold = 0,
  min.pct = 0,
  only.pos = FALSE
)

sex_markers$Gene <- rownames(sex_markers)

p_val_cutoff <- 0.05
log2fc_cutoff <- 0.5

female_up_genes <- sex_markers %>%
  filter(
    p_val_adj < p_val_cutoff,
    avg_log2FC < -log2fc_cutoff
  ) %>%
  pull(Gene) %>%
  unique()

male_up_genes <- sex_markers %>%
  filter(
    p_val_adj < p_val_cutoff,
    avg_log2FC > log2fc_cutoff
  ) %>%
  pull(Gene) %>%
  unique()

female_only_up_genes <- setdiff(female_up_genes, male_up_genes)
male_only_up_genes <- setdiff(male_up_genes, female_up_genes)

# ============================================================
# 2. DEFINE GENES EXPRESSED IN BOTH MALE AND FEMALE FIBROBLASTS
# ============================================================

expr_mat <- tryCatch(
  GetAssayData(fibroblast_subset, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(fibroblast_subset, assay = "RNA", slot = "data")
)

male_cells <- colnames(fibroblast_subset)[fibroblast_subset$Sex == "Male"]
female_cells <- colnames(fibroblast_subset)[fibroblast_subset$Sex == "Female"]

pct_male_expressed <- Matrix::rowMeans(expr_mat[, male_cells, drop = FALSE] > 0)
pct_female_expressed <- Matrix::rowMeans(expr_mat[, female_cells, drop = FALSE] > 0)

min_pct_shared_expression <- 0.05

shared_expressed_genes <- names(which(
  pct_male_expressed >= min_pct_shared_expression &
    pct_female_expressed >= min_pct_shared_expression
))

# ============================================================
# 3. PRINT COUNTS
# ============================================================

cat("Female-only upregulated genes:", length(female_only_up_genes), "\n")
cat("Male-only upregulated genes:", length(male_only_up_genes), "\n")
cat("Genes expressed in both sexes:", length(shared_expressed_genes), "\n")

# ============================================================
# 4. SAVE GENE LISTS
# ============================================================

venn_gene_lists <- bind_rows(
  data.frame(
    Category = "Female-only upregulated",
    Gene = female_only_up_genes
  ),
  data.frame(
    Category = "Expressed in both sexes",
    Gene = shared_expressed_genes
  ),
  data.frame(
    Category = "Male-only upregulated",
    Gene = male_only_up_genes
  )
)

write.csv(
  venn_gene_lists,
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Female_Male_Upregulated_SharedExpressed_GeneLists.csv",
  row.names = FALSE
)

# ============================================================
# 5. CREATE VENN-STYLE FIGURE
# ============================================================

circle_data <- data.frame(
  group = c("Female", "Male"),
  x0 = c(-0.75, 0.75),
  y0 = c(0, 0),
  r = c(1.25, 1.25)
)

venn_plot <- ggplot() +
  geom_circle(
    data = circle_data,
    aes(x0 = x0, y0 = y0, r = r, fill = group, color = group),
    alpha = 0.45,
    linewidth = 1.3
  ) +
  annotate(
    "text",
    x = -1.25,
    y = 0,
    label = paste0(
      length(female_only_up_genes),
      "\nFemale-only\nupregulated"
    ),
    size = 5,
    fontface = "bold",
    color = "#8B0000"
  ) +
  annotate(
    "text",
    x = 0,
    y = 0,
    label = paste0(
      length(shared_expressed_genes),
      "\nExpressed\nin both"
    ),
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  annotate(
    "text",
    x = 1.25,
    y = 0,
    label = paste0(
      length(male_only_up_genes),
      "\nMale-only\nupregulated"
    ),
    size = 5,
    fontface = "bold",
    color = "#004C99"
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#D55E00",
      "Male" = "#0072B2"
    )
  ) +
  scale_color_manual(
    values = c(
      "Female" = "#8B0000",
      "Male" = "#004C99"
    )
  ) +
  coord_fixed() +
  xlim(-2.3, 2.3) +
  ylim(-1.5, 1.5) +
  ggtitle("Sex-Biased and Shared Expressed Genes in PDAC Fibroblasts") +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16,
      margin = margin(b = 15)
    ),
    legend.position = "none"
  )

print(venn_plot)

### 3F: Venn Diagram for General Cells

# ============================================================
# VENN-STYLE DIAGRAM:
# ALL GENERAL CELL TYPES
# FEMALE-UPREGULATED, MALE-UPREGULATED, AND SHARED EXPRESSED GENES
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)

if (!requireNamespace("ggforce", quietly = TRUE)) {
  install.packages("ggforce")
}
library(ggforce)

Merged_Joined <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Joined.RDS"
)

DefaultAssay(Merged_Joined) <- "RNA"

Merged_Joined <- tryCatch(
  JoinLayers(Merged_Joined),
  error = function(e) Merged_Joined
)

# ============================================================
# 1. ADD SEX METADATA
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

Merged_Joined$Sex <- unname(
  sex_mapping[as.character(Merged_Joined$orig.ident)]
)

Merged_Joined <- subset(
  Merged_Joined,
  subset = Sex %in% c("Male", "Female")
)

print(table(Merged_Joined$Sex))

# ============================================================
# 2. RUN MALE VS FEMALE DGE ACROSS ALL GENERAL CELL TYPES
# Positive avg_log2FC = higher in males
# Negative avg_log2FC = higher in females
# ============================================================

sex_markers_all_celltypes <- FindMarkers(
  object = Merged_Joined,
  ident.1 = "Male",
  ident.2 = "Female",
  group.by = "Sex",
  logfc.threshold = 0,
  min.pct = 0,
  only.pos = FALSE
)

sex_markers_all_celltypes$Gene <- rownames(sex_markers_all_celltypes)

p_val_cutoff <- 0.05
log2fc_cutoff <- 0.5

female_up_genes_all <- sex_markers_all_celltypes %>%
  filter(
    p_val_adj < p_val_cutoff,
    avg_log2FC < -log2fc_cutoff
  ) %>%
  pull(Gene) %>%
  unique()

male_up_genes_all <- sex_markers_all_celltypes %>%
  filter(
    p_val_adj < p_val_cutoff,
    avg_log2FC > log2fc_cutoff
  ) %>%
  pull(Gene) %>%
  unique()

female_only_up_genes_all <- setdiff(
  female_up_genes_all,
  male_up_genes_all
)

male_only_up_genes_all <- setdiff(
  male_up_genes_all,
  female_up_genes_all
)

# ============================================================
# 3. DEFINE GENES EXPRESSED IN BOTH SEXES ACROSS ALL GENERAL CELL TYPES
# ============================================================

expr_mat <- tryCatch(
  GetAssayData(Merged_Joined, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(Merged_Joined, assay = "RNA", slot = "data")
)

male_cells <- colnames(Merged_Joined)[Merged_Joined$Sex == "Male"]
female_cells <- colnames(Merged_Joined)[Merged_Joined$Sex == "Female"]

pct_male_expressed <- Matrix::rowMeans(
  expr_mat[, male_cells, drop = FALSE] > 0
)

pct_female_expressed <- Matrix::rowMeans(
  expr_mat[, female_cells, drop = FALSE] > 0
)

min_pct_shared_expression <- 0.05

shared_expressed_genes_all <- names(which(
  pct_male_expressed >= min_pct_shared_expression &
    pct_female_expressed >= min_pct_shared_expression
))

# ============================================================
# 4. PRINT COUNTS
# ============================================================

cat("Female-only upregulated genes across all cell types:",
    length(female_only_up_genes_all), "\n")

cat("Male-only upregulated genes across all cell types:",
    length(male_only_up_genes_all), "\n")

cat("Genes expressed in both sexes across all cell types:",
    length(shared_expressed_genes_all), "\n")

# ============================================================
# 5. SAVE GENE LISTS
# ============================================================

venn_gene_lists_all_celltypes <- bind_rows(
  data.frame(
    Category = "Female-only upregulated",
    Gene = female_only_up_genes_all
  ),
  data.frame(
    Category = "Expressed in both sexes",
    Gene = shared_expressed_genes_all
  ),
  data.frame(
    Category = "Male-only upregulated",
    Gene = male_only_up_genes_all
  )
)

write.csv(
  venn_gene_lists_all_celltypes,
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/AllCellTypes_Female_Male_Upregulated_SharedExpressed_GeneLists.csv",
  row.names = FALSE
)

write.csv(
  sex_markers_all_celltypes,
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/AllCellTypes_Male_vs_Female_DGE.csv",
  row.names = FALSE
)

# ============================================================
# 6. CREATE VENN-STYLE FIGURE
# ============================================================

circle_data <- data.frame(
  group = c("Female", "Male"),
  x0 = c(-0.75, 0.75),
  y0 = c(0, 0),
  r = c(1.25, 1.25)
)

venn_plot_all_celltypes <- ggplot() +
  geom_circle(
    data = circle_data,
    aes(
      x0 = x0,
      y0 = y0,
      r = r,
      fill = group,
      color = group
    ),
    alpha = 0.45,
    linewidth = 1.3
  ) +
  annotate(
    "text",
    x = -1.25,
    y = 0,
    label = paste0(
      length(female_only_up_genes_all),
      "\nFemale-only\nupregulated"
    ),
    size = 5,
    fontface = "bold",
    color = "#8B0000"
  ) +
  annotate(
    "text",
    x = 0,
    y = 0,
    label = paste0(
      length(shared_expressed_genes_all),
      "\nExpressed\nin both"
    ),
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  annotate(
    "text",
    x = 1.25,
    y = 0,
    label = paste0(
      length(male_only_up_genes_all),
      "\nMale-only\nupregulated"
    ),
    size = 5,
    fontface = "bold",
    color = "#004C99"
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#D55E00",
      "Male" = "#0072B2"
    )
  ) +
  scale_color_manual(
    values = c(
      "Female" = "#8B0000",
      "Male" = "#004C99"
    )
  ) +
  coord_fixed() +
  xlim(-2.3, 2.3) +
  ylim(-1.5, 1.5) +
  ggtitle("Sex-Biased and Shared Expressed Genes Across All Cell Types") +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16,
      margin = margin(b = 15)
    ),
    legend.position = "none"
  )

print(venn_plot_all_celltypes)

ggsave(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/AllCellTypes_Female_Male_Upregulated_SharedExpressed_Venn.pdf",
  plot = venn_plot_all_celltypes,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/AllCellTypes_Female_Male_Upregulated_SharedExpressed_Venn.png",
  plot = venn_plot_all_celltypes,
  width = 8,
  height = 5,
  dpi = 300
)

### Part 3 and 4

library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)

if (!requireNamespace("ggforce", quietly = TRUE)) {
  install.packages("ggforce")
}
library(ggforce)

fibroblast_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS"
)

DefaultAssay(fibroblast_subset) <- "RNA"

fibroblast_subset <- tryCatch(
  JoinLayers(fibroblast_subset),
  error = function(e) fibroblast_subset
)

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

fibroblast_subset$Sex <- unname(
  sex_mapping[as.character(fibroblast_subset$orig.ident)]
)

fibroblast_subset$CAF_type <- case_when(
  as.character(fibroblast_subset$seurat_clusters) == "0" ~ "Fibroblasts",
  as.character(fibroblast_subset$seurat_clusters) == "1" ~ "myCAF",
  as.character(fibroblast_subset$seurat_clusters) == "2" ~ "iCAF",
  as.character(fibroblast_subset$seurat_clusters) == "3" ~ "Pancreatic stellate",
  as.character(fibroblast_subset$seurat_clusters) == "4" ~ "apCAF",
  TRUE ~ "Other"
)

fibroblast_subset <- subset(
  fibroblast_subset,
  subset = Sex %in% c("Male", "Female") & CAF_type %in% c("iCAF", "myCAF")
)

make_caf_venn <- function(seurat_obj, caf_name, output_prefix) {
  caf_subset <- subset(
    seurat_obj,
    subset = CAF_type == caf_name
  )
  
  sex_markers <- FindMarkers(
    object = caf_subset,
    ident.1 = "Male",
    ident.2 = "Female",
    group.by = "Sex",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = FALSE
  )
  
  sex_markers$Gene <- rownames(sex_markers)
  
  p_val_cutoff <- 0.05
  log2fc_cutoff <- 0.5
  
  female_up_genes <- sex_markers %>%
    filter(p_val_adj < p_val_cutoff, avg_log2FC < -log2fc_cutoff) %>%
    pull(Gene) %>%
    unique()
  
  male_up_genes <- sex_markers %>%
    filter(p_val_adj < p_val_cutoff, avg_log2FC > log2fc_cutoff) %>%
    pull(Gene) %>%
    unique()
  
  female_only_up_genes <- setdiff(female_up_genes, male_up_genes)
  male_only_up_genes <- setdiff(male_up_genes, female_up_genes)
  
  expr_mat <- tryCatch(
    GetAssayData(caf_subset, assay = "RNA", layer = "data"),
    error = function(e) GetAssayData(caf_subset, assay = "RNA", slot = "data")
  )
  
  male_cells <- colnames(caf_subset)[caf_subset$Sex == "Male"]
  female_cells <- colnames(caf_subset)[caf_subset$Sex == "Female"]
  
  pct_male_expressed <- Matrix::rowMeans(expr_mat[, male_cells, drop = FALSE] > 0)
  pct_female_expressed <- Matrix::rowMeans(expr_mat[, female_cells, drop = FALSE] > 0)
  
  min_pct_shared_expression <- 0.05
  
  shared_expressed_genes <- names(which(
    pct_male_expressed >= min_pct_shared_expression &
      pct_female_expressed >= min_pct_shared_expression
  ))
  
  venn_gene_lists <- bind_rows(
    data.frame(Category = "Female-only upregulated", Gene = female_only_up_genes),
    data.frame(Category = "Expressed in both sexes", Gene = shared_expressed_genes),
    data.frame(Category = "Male-only upregulated", Gene = male_only_up_genes)
  )
  
  write.csv(
    venn_gene_lists,
    paste0("/Users/sudharshiniram/Desktop/SRI Data Analysis/", output_prefix, "_Venn_GeneLists.csv"),
    row.names = FALSE
  )
  
  write.csv(
    sex_markers,
    paste0("/Users/sudharshiniram/Desktop/SRI Data Analysis/", output_prefix, "_Male_vs_Female_DGE.csv"),
    row.names = FALSE
  )
  
  circle_data <- data.frame(
    group = c("Female", "Male"),
    x0 = c(-0.75, 0.75),
    y0 = c(0, 0),
    r = c(1.25, 1.25)
  )
  
  venn_plot <- ggplot() +
    geom_circle(
      data = circle_data,
      aes(x0 = x0, y0 = y0, r = r, fill = group, color = group),
      alpha = 0.45,
      linewidth = 1.3
    ) +
    annotate(
      "text",
      x = -1.25,
      y = 0,
      label = paste0(length(female_only_up_genes), "\nFemale-only\nupregulated"),
      size = 5,
      fontface = "bold",
      color = "#8B0000"
    ) +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = paste0(length(shared_expressed_genes), "\nExpressed\nin both"),
      size = 5,
      fontface = "bold",
      color = "black"
    ) +
    annotate(
      "text",
      x = 1.25,
      y = 0,
      label = paste0(length(male_only_up_genes), "\nMale-only\nupregulated"),
      size = 5,
      fontface = "bold",
      color = "#004C99"
    ) +
    scale_fill_manual(values = c("Female" = "#D55E00", "Male" = "#0072B2")) +
    scale_color_manual(values = c("Female" = "#8B0000", "Male" = "#004C99")) +
    coord_fixed() +
    xlim(-2.3, 2.3) +
    ylim(-1.5, 1.5) +
    ggtitle(paste0(caf_name, ": Sex-Biased and Shared Expressed Genes")) +
    theme_void(base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
      legend.position = "none"
    )
  
  print(venn_plot)
  
  ggsave(
    paste0("/Users/sudharshiniram/Desktop/SRI Data Analysis/", output_prefix, "_Venn.pdf"),
    plot = venn_plot,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  ggsave(
    paste0("/Users/sudharshiniram/Desktop/SRI Data Analysis/", output_prefix, "_Venn.png"),
    plot = venn_plot,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  return(venn_plot)
}

icaf_venn <- make_caf_venn(
  fibroblast_subset,
  caf_name = "iCAF",
  output_prefix = "iCAF_Female_Male_Upregulated_SharedExpressed"
)

mycaf_venn <- make_caf_venn(
  fibroblast_subset,
  caf_name = "myCAF",
  output_prefix = "myCAF_Female_Male_Upregulated_SharedExpressed"
)



### 3F: Part 5 Schwann Cells

library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)

if (!requireNamespace("ggforce", quietly = TRUE)) {
  install.packages("ggforce")
}
library(ggforce)

schwann_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/true_schwann_ERBB3_PMP22_subset.RDS"
)

DefaultAssay(schwann_subset) <- "RNA"

schwann_subset <- tryCatch(
  JoinLayers(schwann_subset),
  error = function(e) schwann_subset
)

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

schwann_subset$Sex <- unname(
  sex_mapping[as.character(schwann_subset$orig.ident)]
)

schwann_subset <- subset(
  schwann_subset,
  subset = Sex %in% c("Male", "Female")
)

print(table(schwann_subset$Sex))

sex_markers <- FindMarkers(
  object = schwann_subset,
  ident.1 = "Male",
  ident.2 = "Female",
  group.by = "Sex",
  logfc.threshold = 0,
  min.pct = 0,
  only.pos = FALSE
)

sex_markers$Gene <- rownames(sex_markers)

p_val_cutoff <- 0.05
log2fc_cutoff <- 0.5

female_up_genes <- sex_markers %>%
  filter(p_val_adj < p_val_cutoff, avg_log2FC < -log2fc_cutoff) %>%
  pull(Gene) %>%
  unique()

male_up_genes <- sex_markers %>%
  filter(p_val_adj < p_val_cutoff, avg_log2FC > log2fc_cutoff) %>%
  pull(Gene) %>%
  unique()

female_only_up_genes <- setdiff(female_up_genes, male_up_genes)
male_only_up_genes <- setdiff(male_up_genes, female_up_genes)

expr_mat <- tryCatch(
  GetAssayData(schwann_subset, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(schwann_subset, assay = "RNA", slot = "data")
)

male_cells <- colnames(schwann_subset)[schwann_subset$Sex == "Male"]
female_cells <- colnames(schwann_subset)[schwann_subset$Sex == "Female"]

pct_male_expressed <- Matrix::rowMeans(expr_mat[, male_cells, drop = FALSE] > 0)
pct_female_expressed <- Matrix::rowMeans(expr_mat[, female_cells, drop = FALSE] > 0)

min_pct_shared_expression <- 0.05

shared_expressed_genes <- names(which(
  pct_male_expressed >= min_pct_shared_expression &
    pct_female_expressed >= min_pct_shared_expression
))

cat("Female-only upregulated genes:", length(female_only_up_genes), "\n")
cat("Male-only upregulated genes:", length(male_only_up_genes), "\n")
cat("Genes expressed in both sexes:", length(shared_expressed_genes), "\n")

venn_gene_lists <- bind_rows(
  data.frame(Category = "Female-only upregulated", Gene = female_only_up_genes),
  data.frame(Category = "Expressed in both sexes", Gene = shared_expressed_genes),
  data.frame(Category = "Male-only upregulated", Gene = male_only_up_genes)
)

write.csv(
  venn_gene_lists,
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/TrueSchwann_ERBB3_PMP22_Female_Male_Upregulated_SharedExpressed_Venn_GeneLists.csv",
  row.names = FALSE
)

write.csv(
  sex_markers,
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/TrueSchwann_ERBB3_PMP22_Male_vs_Female_DGE.csv",
  row.names = FALSE
)

circle_data <- data.frame(
  group = c("Female", "Male"),
  x0 = c(-0.75, 0.75),
  y0 = c(0, 0),
  r = c(1.25, 1.25)
)

schwann_venn <- ggplot() +
  geom_circle(
    data = circle_data,
    aes(x0 = x0, y0 = y0, r = r, fill = group, color = group),
    alpha = 0.45,
    linewidth = 1.3
  ) +
  annotate(
    "text",
    x = -1.25,
    y = 0,
    label = paste0(length(female_only_up_genes), "\nFemale-only\nupregulated"),
    size = 5,
    fontface = "bold",
    color = "#8B0000"
  ) +
  annotate(
    "text",
    x = 0,
    y = 0,
    label = paste0(length(shared_expressed_genes), "\nExpressed\nin both"),
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  annotate(
    "text",
    x = 1.25,
    y = 0,
    label = paste0(length(male_only_up_genes), "\nMale-only\nupregulated"),
    size = 5,
    fontface = "bold",
    color = "#004C99"
  ) +
  scale_fill_manual(values = c("Female" = "#D55E00", "Male" = "#0072B2")) +
  scale_color_manual(values = c("Female" = "#8B0000", "Male" = "#004C99")) +
  coord_fixed() +
  xlim(-2.3, 2.3) +
  ylim(-1.5, 1.5) +
  ggtitle("Schwann Cells: Sex-Biased and Shared Expressed Genes") +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    legend.position = "none"
  )

print(schwann_venn)

ggsave(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/TrueSchwann_ERBB3_PMP22_Female_Male_Upregulated_SharedExpressed_Venn.pdf",
  plot = schwann_venn,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/TrueSchwann_ERBB3_PMP22_Female_Male_Upregulated_SharedExpressed_Venn.png",
  plot = schwann_venn,
  width = 8,
  height = 5,
  dpi = 300
)


### 4A: Part 2 CellChat Interactions between Cell Types

library(Seurat)
library(CellChat)
library(dplyr)
library(circlize)
library(patchwork)

Merged_Joined <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Joined.RDS")
myeloid_subset <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/myeloid_subset.RDS")
BCell_Subset <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/BCell_Subset.RDS")
TCell_Subset <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/TCell_Subset.RDS")

DefaultAssay(Merged_Joined) <- "RNA"

Merged_Joined <- tryCatch(
  JoinLayers(Merged_Joined),
  error = function(e) Merged_Joined
)

Merged_Joined$cellchat_celltype <- "Other"

fibroblast_cells <- colnames(Merged_Joined)[as.character(Merged_Joined$seurat_clusters) == "1"]

if ("true_schwann_subset" %in% colnames(Merged_Joined@meta.data)) {
  schwann_cells <- colnames(Merged_Joined)[
    as.character(Merged_Joined$true_schwann_subset) %in% c("True Schwann", "TRUE", "True", "true")
  ]
} else {
  schwann_cells <- colnames(Merged_Joined)[as.character(Merged_Joined$seurat_clusters) == "11"]
}

b_cells <- intersect(colnames(BCell_Subset), colnames(Merged_Joined))
t_cells <- intersect(colnames(TCell_Subset), colnames(Merged_Joined))

if ("myeloid_subcluster" %in% colnames(myeloid_subset@meta.data)) {
  myeloid_clusters <- as.character(myeloid_subset$myeloid_subcluster)
} else if ("seurat_clusters" %in% colnames(myeloid_subset@meta.data)) {
  myeloid_clusters <- as.character(myeloid_subset$seurat_clusters)
} else {
  myeloid_clusters <- as.character(Idents(myeloid_subset))
}

macrophage_cells <- colnames(myeloid_subset)[myeloid_clusters %in% c("0", "1")]
dendritic_cells <- colnames(myeloid_subset)[myeloid_clusters == "3"]

macrophage_cells <- intersect(macrophage_cells, colnames(Merged_Joined))
dendritic_cells <- intersect(dendritic_cells, colnames(Merged_Joined))

Merged_Joined$cellchat_celltype[colnames(Merged_Joined) %in% schwann_cells] <- "Schwann cells"
Merged_Joined$cellchat_celltype[colnames(Merged_Joined) %in% fibroblast_cells] <- "Fibroblasts"
Merged_Joined$cellchat_celltype[colnames(Merged_Joined) %in% b_cells] <- "B cells"
Merged_Joined$cellchat_celltype[colnames(Merged_Joined) %in% t_cells] <- "T cells"
Merged_Joined$cellchat_celltype[colnames(Merged_Joined) %in% macrophage_cells] <- "Macrophages"
Merged_Joined$cellchat_celltype[colnames(Merged_Joined) %in% dendritic_cells] <- "Dendritic cells"

cellchat_groups <- c(
  "Schwann cells",
  "Fibroblasts",
  "B cells",
  "T cells",
  "Macrophages",
  "Dendritic cells"
)

cellchat_seurat <- subset(
  Merged_Joined,
  subset = cellchat_celltype %in% cellchat_groups
)

cellchat_seurat$cellchat_celltype <- factor(
  cellchat_seurat$cellchat_celltype,
  levels = cellchat_groups
)

print(table(cellchat_seurat$cellchat_celltype))

if (any(table(cellchat_seurat$cellchat_celltype) < 10)) {
  warning("One or more groups has fewer than 10 cells. CellChat may filter that group out.")
}

cellchat_obj <- createCellChat(
  object = cellchat_seurat,
  group.by = "cellchat_celltype",
  assay = "RNA"
)

CellChatDB <- CellChatDB.human
cellchat_obj@DB <- CellChatDB

cellchat_obj <- subsetData(cellchat_obj)
cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)

cellchat_obj <- computeCommunProb(
  cellchat_obj,
  type = "triMean",
  raw.use = TRUE
)

cellchat_obj <- filterCommunication(
  cellchat_obj,
  min.cells = 10
)

cellchat_obj <- computeCommunProbPathway(cellchat_obj)
cellchat_obj <- aggregateNet(cellchat_obj)

group_size <- as.numeric(table(cellchat_obj@idents))

cellchat_colors <- c(
  "Schwann cells" = "#E69F00",
  "Fibroblasts" = "#009E73",
  "B cells" = "#CC79A7",
  "T cells" = "#0072B2",
  "Macrophages" = "#D55E00",
  "Dendritic cells" = "#56B4E9"
)

cellchat_colors <- cellchat_colors[cellchat_groups]

pdf(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/CellChat_Circle_Number_of_Interactions.pdf",
  width = 8,
  height = 8
)

netVisual_circle(
  cellchat_obj@net$count,
  vertex.weight = group_size,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = cellchat_colors,
  title.name = "Number of CellChat Interactions"
)

dev.off()

pdf(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/CellChat_Circle_Interaction_Strength.pdf",
  width = 8,
  height = 8
)

netVisual_circle(
  cellchat_obj@net$weight,
  vertex.weight = group_size,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = cellchat_colors,
  title.name = "CellChat Interaction Strength"
)

dev.off()

pdf(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/CellChat_Chord_All_Ligand_Receptor_Interactions.pdf",
  width = 10,
  height = 10
)

netVisual_chord_gene(
  cellchat_obj,
  sources.use = cellchat_groups,
  targets.use = cellchat_groups,
  slot.name = "net",
  color.use = cellchat_colors,
  lab.cex = 0.7,
  legend.pos.y = 30
)

dev.off()

netVisual_circle(
  cellchat_obj@net$count,
  vertex.weight = group_size,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = cellchat_colors,
  title.name = "Number of CellChat Interactions"
)

netVisual_circle(
  cellchat_obj@net$weight,
  vertex.weight = group_size,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = cellchat_colors,
  title.name = "CellChat Interaction Strength"
)

netVisual_chord_gene(
  cellchat_obj,
  sources.use = cellchat_groups,
  targets.use = cellchat_groups,
  slot.name = "net",
  color.use = cellchat_colors,
  lab.cex = 0.7,
  legend.pos.y = 30
)

saveRDS(
  cellchat_obj,
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/CellChat_Schwann_Fibroblast_B_T_Macrophage_DC.rds"
)


### 4A: Part 3 CellChat Male vs Female

library(Seurat)
library(CellChat)
library(dplyr)
library(circlize)

outdir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis/CellChat_Publication_Figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

Merged_Joined <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Joined.RDS")
myeloid_subset <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/myeloid_subset.RDS")
BCell_Subset <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/BCell_Subset.RDS")
TCell_Subset <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/TCell_Subset.RDS")

DefaultAssay(Merged_Joined) <- "RNA"

Merged_Joined <- tryCatch(
  JoinLayers(Merged_Joined),
  error = function(e) Merged_Joined
)

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

Merged_Joined$Sex <- unname(sex_mapping[as.character(Merged_Joined$orig.ident)])
Merged_Joined$cellchat_celltype <- "Other"

fibroblast_cells <- colnames(Merged_Joined)[as.character(Merged_Joined$seurat_clusters) == "1"]

if ("true_schwann_subset" %in% colnames(Merged_Joined@meta.data)) {
  schwann_cells <- colnames(Merged_Joined)[
    as.character(Merged_Joined$true_schwann_subset) %in% c("True Schwann", "TRUE", "True", "true")
  ]
} else {
  schwann_cells <- colnames(Merged_Joined)[as.character(Merged_Joined$seurat_clusters) == "11"]
}

b_cells <- intersect(colnames(BCell_Subset), colnames(Merged_Joined))
t_cells <- intersect(colnames(TCell_Subset), colnames(Merged_Joined))

if ("myeloid_subcluster" %in% colnames(myeloid_subset@meta.data)) {
  myeloid_clusters <- as.character(myeloid_subset$myeloid_subcluster)
} else if ("seurat_clusters" %in% colnames(myeloid_subset@meta.data)) {
  myeloid_clusters <- as.character(myeloid_subset$seurat_clusters)
} else {
  myeloid_clusters <- as.character(Idents(myeloid_subset))
}

macrophage_cells <- intersect(
  colnames(myeloid_subset)[myeloid_clusters %in% c("0", "1")],
  colnames(Merged_Joined)
)

dendritic_cells <- intersect(
  colnames(myeloid_subset)[myeloid_clusters == "3"],
  colnames(Merged_Joined)
)

Merged_Joined$cellchat_celltype[colnames(Merged_Joined) %in% schwann_cells] <- "Schwann cells"
Merged_Joined$cellchat_celltype[colnames(Merged_Joined) %in% fibroblast_cells] <- "Fibroblasts"
Merged_Joined$cellchat_celltype[colnames(Merged_Joined) %in% b_cells] <- "B cells"
Merged_Joined$cellchat_celltype[colnames(Merged_Joined) %in% t_cells] <- "T cells"
Merged_Joined$cellchat_celltype[colnames(Merged_Joined) %in% macrophage_cells] <- "Macrophages"
Merged_Joined$cellchat_celltype[colnames(Merged_Joined) %in% dendritic_cells] <- "Dendritic cells"

cellchat_groups <- c(
  "Schwann cells",
  "Fibroblasts",
  "B cells",
  "T cells",
  "Macrophages",
  "Dendritic cells"
)

cellchat_colors <- c(
  "Schwann cells" = "#E69F00",
  "Fibroblasts" = "#009E73",
  "B cells" = "#CC79A7",
  "T cells" = "#0072B2",
  "Macrophages" = "#D55E00",
  "Dendritic cells" = "#56B4E9"
)

prepare_cellchat_seurat <- function(seurat_obj, sex_filter = NULL) {
  if (is.null(sex_filter)) {
    obj <- subset(
      seurat_obj,
      subset = cellchat_celltype %in% cellchat_groups & Sex %in% c("Male", "Female")
    )
  } else {
    obj <- subset(
      seurat_obj,
      subset = cellchat_celltype %in% cellchat_groups & Sex == sex_filter
    )
  }
  
  obj$cellchat_celltype <- factor(
    obj$cellchat_celltype,
    levels = cellchat_groups
  )
  
  obj$cellchat_celltype <- droplevels(obj$cellchat_celltype)
  return(obj)
}

run_cellchat <- function(seurat_obj, min_cells_filter = 10) {
  print(table(seurat_obj$cellchat_celltype))
  
  if (any(table(seurat_obj$cellchat_celltype) < min_cells_filter)) {
    warning("One or more groups has fewer than the min_cells_filter. Consider lowering min_cells_filter.")
  }
  
  cellchat_obj <- createCellChat(
    object = seurat_obj,
    group.by = "cellchat_celltype",
    assay = "RNA"
  )
  
  cellchat_obj@DB <- CellChatDB.human
  cellchat_obj <- subsetData(cellchat_obj)
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  
  cellchat_obj <- computeCommunProb(
    cellchat_obj,
    type = "triMean",
    raw.use = TRUE
  )
  
  cellchat_obj <- filterCommunication(
    cellchat_obj,
    min.cells = min_cells_filter
  )
  
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)
  cellchat_obj <- aggregateNet(cellchat_obj)
  return(cellchat_obj)
}

save_circle_plot <- function(cellchat_obj, net_type, title_text, file_prefix) {
  net_matrix <- if (net_type == "count") {
    cellchat_obj@net$count
  } else {
    cellchat_obj@net$weight
  }
  
  plot_groups <- rownames(net_matrix)
  group_size <- as.numeric(table(cellchat_obj@idents)[plot_groups])
  color_use <- cellchat_colors[plot_groups]
  
  pdf(
    file.path(outdir, paste0(file_prefix, ".pdf")),
    width = 8,
    height = 8,
    useDingbats = FALSE
  )
  
  par(mar = c(1, 1, 4, 1), xpd = TRUE)
  
  netVisual_circle(
    net_matrix,
    vertex.weight = group_size,
    vertex.label.cex = 1.15,
    vertex.label.color = "black",
    vertex.size.max = 18,
    weight.scale = TRUE,
    edge.width.max = 8,
    alpha.edge = 0.65,
    label.edge = FALSE,
    color.use = color_use,
    title.name = title_text
  )
  
  dev.off()
  
  png(
    file.path(outdir, paste0(file_prefix, ".png")),
    width = 8,
    height = 8,
    units = "in",
    res = 600
  )
  
  par(mar = c(1, 1, 4, 1), xpd = TRUE)
  
  netVisual_circle(
    net_matrix,
    vertex.weight = group_size,
    vertex.label.cex = 1.15,
    vertex.label.color = "black",
    vertex.size.max = 18,
    weight.scale = TRUE,
    edge.width.max = 8,
    alpha.edge = 0.65,
    label.edge = FALSE,
    color.use = color_use,
    title.name = title_text
  )
  
  dev.off()
}

save_chord_plot <- function(cellchat_obj, title_text, file_prefix) {
  plot_groups <- levels(cellchat_obj@idents)
  color_use <- cellchat_colors[plot_groups]
  
  pdf(
    file.path(outdir, paste0(file_prefix, ".pdf")),
    width = 12,
    height = 10,
    useDingbats = FALSE
  )
  
  par(mar = c(1, 1, 4, 8), xpd = TRUE)
  circos.clear()
  
  netVisual_chord_gene(
    cellchat_obj,
    sources.use = plot_groups,
    targets.use = plot_groups,
    slot.name = "net",
    color.use = color_use,
    lab.cex = 0.45,
    legend.pos.x = 8,
    legend.pos.y = 4
  )
  
  title(main = title_text, cex.main = 1.4, font.main = 2)
  circos.clear()
  dev.off()
  
  png(
    file.path(outdir, paste0(file_prefix, ".png")),
    width = 12,
    height = 10,
    units = "in",
    res = 600
  )
  
  par(mar = c(1, 1, 4, 8), xpd = TRUE)
  circos.clear()
  
  netVisual_chord_gene(
    cellchat_obj,
    sources.use = plot_groups,
    targets.use = plot_groups,
    slot.name = "net",
    color.use = color_use,
    lab.cex = 0.45,
    legend.pos.x = 8,
    legend.pos.y = 4
  )
  
  title(main = title_text, cex.main = 1.4, font.main = 2)
  circos.clear()
  dev.off()
}

cellchat_seurat_all <- prepare_cellchat_seurat(Merged_Joined)
cellchat_all <- run_cellchat(cellchat_seurat_all, min_cells_filter = 10)

save_circle_plot(
  cellchat_all,
  net_type = "count",
  title_text = "CellChat Network: Number of Ligand-Receptor Interactions",
  file_prefix = "Figure_CellChat_AllSamples_Number_of_Interactions"
)

save_circle_plot(
  cellchat_all,
  net_type = "weight",
  title_text = "CellChat Network: Overall Interaction Strength",
  file_prefix = "Figure_CellChat_AllSamples_Interaction_Strength"
)

save_chord_plot(
  cellchat_all,
  title_text = "CellChat Ligand-Receptor Interactions Across Cell Types",
  file_prefix = "Figure_CellChat_AllSamples_Ligand_Receptor_Chord"
)

write.csv(
  subsetCommunication(cellchat_all),
  file.path(outdir, "CellChat_AllSamples_All_Interactions.csv"),
  row.names = FALSE
)

saveRDS(
  cellchat_all,
  file.path(outdir, "CellChat_AllSamples_Object.rds")
)

cellchat_seurat_female <- prepare_cellchat_seurat(Merged_Joined, sex_filter = "Female")
cellchat_seurat_male <- prepare_cellchat_seurat(Merged_Joined, sex_filter = "Male")

cellchat_female <- run_cellchat(cellchat_seurat_female, min_cells_filter = 10)
cellchat_male <- run_cellchat(cellchat_seurat_male, min_cells_filter = 10)

saveRDS(
  cellchat_female,
  file.path(outdir, "CellChat_Female_Object.rds")
)

saveRDS(
  cellchat_male,
  file.path(outdir, "CellChat_Male_Object.rds")
)

write.csv(
  subsetCommunication(cellchat_female),
  file.path(outdir, "CellChat_Female_All_Interactions.csv"),
  row.names = FALSE
)

write.csv(
  subsetCommunication(cellchat_male),
  file.path(outdir, "CellChat_Male_All_Interactions.csv"),
  row.names = FALSE
)

female_net <- cellchat_female@net$weight
male_net <- cellchat_male@net$weight

female_groups <- rownames(female_net)
male_groups <- rownames(male_net)

female_group_size <- as.numeric(table(cellchat_female@idents)[female_groups])
male_group_size <- as.numeric(table(cellchat_male@idents)[male_groups])

edge_max <- max(female_net, male_net, na.rm = TRUE)
vertex_max <- max(female_group_size, male_group_size, na.rm = TRUE)

pdf(
  file.path(outdir, "Figure_CellChat_Female_vs_Male_Interaction_Strength.pdf"),
  width = 14,
  height = 7,
  useDingbats = FALSE
)

par(mfrow = c(1, 2), mar = c(1, 1, 4, 1), xpd = TRUE)

netVisual_circle(
  female_net,
  vertex.weight = female_group_size,
  vertex.weight.max = vertex_max,
  vertex.label.cex = 0.95,
  vertex.label.color = "black",
  vertex.size.max = 18,
  weight.scale = TRUE,
  edge.weight.max = edge_max,
  edge.width.max = 8,
  alpha.edge = 0.65,
  label.edge = FALSE,
  color.use = cellchat_colors[female_groups],
  title.name = "Female Samples: Interaction Strength"
)

netVisual_circle(
  male_net,
  vertex.weight = male_group_size,
  vertex.weight.max = vertex_max,
  vertex.label.cex = 0.95,
  vertex.label.color = "black",
  vertex.size.max = 18,
  weight.scale = TRUE,
  edge.weight.max = edge_max,
  edge.width.max = 8,
  alpha.edge = 0.65,
  label.edge = FALSE,
  color.use = cellchat_colors[male_groups],
  title.name = "Male Samples: Interaction Strength"
)

dev.off()

png(
  file.path(outdir, "Figure_CellChat_Female_vs_Male_Interaction_Strength.png"),
  width = 14,
  height = 7,
  units = "in",
  res = 600
)

par(mfrow = c(1, 2), mar = c(1, 1, 4, 1), xpd = TRUE)

netVisual_circle(
  female_net,
  vertex.weight = female_group_size,
  vertex.weight.max = vertex_max,
  vertex.label.cex = 0.95,
  vertex.label.color = "black",
  vertex.size.max = 18,
  weight.scale = TRUE,
  edge.weight.max = edge_max,
  edge.width.max = 8,
  alpha.edge = 0.65,
  label.edge = FALSE,
  color.use = cellchat_colors[female_groups],
  title.name = "Female Samples: Interaction Strength"
)

netVisual_circle(
  male_net,
  vertex.weight = male_group_size,
  vertex.weight.max = vertex_max,
  vertex.label.cex = 0.95,
  vertex.label.color = "black",
  vertex.size.max = 18,
  weight.scale = TRUE,
  edge.weight.max = edge_max,
  edge.width.max = 8,
  alpha.edge = 0.65,
  label.edge = FALSE,
  color.use = cellchat_colors[male_groups],
  title.name = "Male Samples: Interaction Strength"
)

dev.off()



### Bar Graph of Cell Chat Interactions


library(CellChat)
library(dplyr)
library(ggplot2)

outdir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis/CellChat_Publication_Figures"

cellchat_all <- readRDS(file.path(outdir, "CellChat_AllSamples_Object.rds"))
cellchat_male <- readRDS(file.path(outdir, "CellChat_Male_Object.rds"))
cellchat_female <- readRDS(file.path(outdir, "CellChat_Female_Object.rds"))

cellchat_colors <- c(
  "Schwann cells" = "#E69F00",
  "Fibroblasts" = "#009E73",
  "B cells" = "#CC79A7",
  "T cells" = "#0072B2",
  "Macrophages" = "#D55E00",
  "Dendritic cells" = "#56B4E9"
)

make_cellchat_barplot <- function(cellchat_obj, dataset_name) {
  net_matrix <- cellchat_obj@net$weight
  
  interaction_df <- as.data.frame(as.table(net_matrix))
  colnames(interaction_df) <- c("Source", "Target", "Interaction_Strength")
  
  interaction_df <- interaction_df %>%
    filter(Interaction_Strength > 0) %>%
    mutate(
      Source = as.character(Source),
      Target = as.character(Target),
      Interaction = paste(Source, "→", Target)
    ) %>%
    arrange(Interaction_Strength)
  
  interaction_df$Interaction <- factor(
    interaction_df$Interaction,
    levels = interaction_df$Interaction
  )
  
  ggplot(
    interaction_df,
    aes(x = Interaction, y = Interaction_Strength, fill = Source)
  ) +
    geom_col(width = 0.75, color = "black", linewidth = 0.25) +
    coord_flip() +
    scale_fill_manual(values = cellchat_colors) +
    labs(
      title = paste0(dataset_name, ": CellChat Interaction Strength"),
      x = "Predicted signaling direction",
      y = "CellChat interaction strength",
      fill = "Source cell type"
    ) +
    theme_classic(base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
      axis.title = element_text(face = "bold", size = 13),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.text.x = element_text(size = 11, color = "black"),
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      legend.position = "right"
    )
}

barplot_all <- make_cellchat_barplot(
  cellchat_all,
  dataset_name = "All Samples"
)

barplot_male <- make_cellchat_barplot(
  cellchat_male,
  dataset_name = "Male Samples"
)

barplot_female <- make_cellchat_barplot(
  cellchat_female,
  dataset_name = "Female Samples"
)

print(barplot_all)
print(barplot_male)
print(barplot_female)


### 1A: Revising title and color to grey

library(Seurat)
library(ggplot2)
library(dplyr)

Merged_Joined <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Joined.RDS"
)

Merged_Joined$CellType_Plot <- case_when(
  as.character(Merged_Joined$seurat_clusters) == "0" ~ "Cancer 1",
  as.character(Merged_Joined$seurat_clusters) == "1" ~ "Fibroblast",
  as.character(Merged_Joined$seurat_clusters) == "2" ~ "Ductal",
  as.character(Merged_Joined$seurat_clusters) == "3" ~ "B Cells",
  as.character(Merged_Joined$seurat_clusters) == "4" ~ "Myeloid",
  as.character(Merged_Joined$seurat_clusters) == "5" ~ "Acinar",
  as.character(Merged_Joined$seurat_clusters) == "6" ~ "Endocrine",
  as.character(Merged_Joined$seurat_clusters) == "7" ~ "Endothelial",
  as.character(Merged_Joined$seurat_clusters) == "8" ~ "T/NK Cells",
  as.character(Merged_Joined$seurat_clusters) == "9" ~ "Cancer 2",
  as.character(Merged_Joined$seurat_clusters) == "10" ~ "Smooth Muscle",
  as.character(Merged_Joined$seurat_clusters) == "12" ~ "Adipocytes",
  as.character(Merged_Joined$seurat_clusters) %in% c("11", "13", "14") ~ "Unlabeled",
  TRUE ~ "Unlabeled"
)

celltype_order <- c(
  "Cancer 1",
  "Fibroblast",
  "Ductal",
  "B Cells",
  "Myeloid",
  "Acinar",
  "Endocrine",
  "Endothelial",
  "T/NK Cells",
  "Cancer 2",
  "Smooth Muscle",
  "Adipocytes",
  "Unlabeled"
)

Merged_Joined$CellType_Plot <- factor(
  Merged_Joined$CellType_Plot,
  levels = celltype_order
)

celltype_colors <- c(
  "Cancer 1" = "#F8766D",
  "Fibroblast" = "#E68613",
  "Ductal" = "#B79F00",
  "B Cells" = "#7CAE00",
  "Myeloid" = "#00A600",
  "Acinar" = "#00BF7D",
  "Endocrine" = "#00BFC4",
  "Endothelial" = "#00A9CF",
  "T/NK Cells" = "#00A6FF",
  "Cancer 2" = "#8494FF",
  "Smooth Muscle" = "#C77CFF",
  "Adipocytes" = "#FF61C3",
  "Unlabeled" = "grey75"
)

umap_coords <- as.data.frame(Embeddings(Merged_Joined, reduction = "umap"))
colnames(umap_coords)[1:2] <- c("umap_1", "umap_2")
umap_coords$CellType_Plot <- Merged_Joined$CellType_Plot

label_data <- umap_coords %>%
  filter(CellType_Plot != "Unlabeled") %>%
  group_by(CellType_Plot) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2),
    .groups = "drop"
  )

umap_plot <- DimPlot(
  Merged_Joined,
  reduction = "umap",
  group.by = "CellType_Plot",
  label = FALSE,
  pt.size = 0.5
) +
  geom_text(
    data = label_data,
    aes(x = umap_1, y = umap_2, label = CellType_Plot),
    color = "black",
    size = 5
  ) +
  scale_color_manual(
    values = celltype_colors,
    breaks = setdiff(celltype_order, "Unlabeled")
  ) +
  ggtitle("UMAP") +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = "right"
  )

print(umap_plot)


### 1C: UMAP by Sex

female_umap <- DimPlot(
  female_obj,
  reduction = "umap",
  group.by = "Cluster_Plot",
  label = FALSE,
  pt.size = 0.5
) +
  geom_text(
    data = female_label_df,
    aes(x = umap_1, y = umap_2, label = Cluster_Plot),
    color = "black",
    size = 5
  ) +
  scale_color_manual(
    values = cluster_colors,
    breaks = legend_breaks,
    labels = legend_labels,
    drop = FALSE
  ) +
  ggtitle("Female") +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )

male_umap <- DimPlot(
  male_obj,
  reduction = "umap",
  group.by = "Cluster_Plot",
  label = FALSE,
  pt.size = 0.5
) +
  geom_text(
    data = male_label_df,
    aes(x = umap_1, y = umap_2, label = Cluster_Plot),
    color = "black",
    size = 5
  ) +
  scale_color_manual(
    values = cluster_colors,
    breaks = legend_breaks,
    labels = legend_labels,
    drop = FALSE
  ) +
  ggtitle("Male") +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )

combined_sex_umap <- female_umap + male_umap +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )

print(combined_sex_umap)


### 1D: Cell Type Percentage for all Samples

library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)

Merged_Joined <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Joined.RDS"
)

plot_data <- Merged_Joined@meta.data

plot_data$CellType <- case_when(
  as.character(plot_data$seurat_clusters) == "0" ~ "Cancer 1",
  as.character(plot_data$seurat_clusters) == "1" ~ "Fibroblast",
  as.character(plot_data$seurat_clusters) == "2" ~ "Ductal",
  as.character(plot_data$seurat_clusters) == "3" ~ "B Cells",
  as.character(plot_data$seurat_clusters) == "4" ~ "Myeloid",
  as.character(plot_data$seurat_clusters) == "5" ~ "Acinar",
  as.character(plot_data$seurat_clusters) == "6" ~ "Endocrine",
  as.character(plot_data$seurat_clusters) == "7" ~ "Endothelial",
  as.character(plot_data$seurat_clusters) == "8" ~ "T / NK Cells",
  as.character(plot_data$seurat_clusters) == "9" ~ "Cancer 2",
  as.character(plot_data$seurat_clusters) == "10" ~ "Smooth Muscle",
  as.character(plot_data$seurat_clusters) == "12" ~ "Adipocytes",
  as.character(plot_data$seurat_clusters) %in% c("11", "13", "14") ~ "Other",
  TRUE ~ "Other"
)

celltype_order <- c(
  "Cancer 1",
  "Fibroblast",
  "Ductal",
  "B Cells",
  "Myeloid",
  "Acinar",
  "Endocrine",
  "Endothelial",
  "T / NK Cells",
  "Cancer 2",
  "Smooth Muscle",
  "Adipocytes",
  "Other"
)

plot_data$CellType <- factor(
  plot_data$CellType,
  levels = celltype_order
)

patient_mapping <- c(
  "Sample_01" = "P01",
  "Sample_02" = "P02",
  "Sample_03" = "P03",
  "Sample_04" = "P04",
  "Sample_05" = "P05",
  "Sample_06" = "P06",
  "Sample_07" = "P07",
  "Sample_08" = "P08",
  "Sample_09" = "P09",
  "Sample_10" = "P10",
  "Sample_11" = "P11",
  "Sample_12" = "P12",
  "Sample_13" = "P13",
  "Sample_14" = "P14",
  "Sample_15" = "P15",
  "Sample_16" = "P16",
  "Sample_17" = "P17"
)

plot_data$PatientID <- patient_mapping[
  as.character(plot_data$orig.ident)
]

plot_data$PatientID <- factor(
  plot_data$PatientID,
  levels = sprintf("P%02d", 1:17)
)

celltype_colors <- c(
  "Cancer 1" = "#F8766D",
  "Fibroblast" = "#E68613",
  "Ductal" = "#B79F00",
  "B Cells" = "#7CAE00",
  "Myeloid" = "#00A600",
  "Acinar" = "#00BF7D",
  "Endocrine" = "#00BFC4",
  "Endothelial" = "#00A9CF",
  "T / NK Cells" = "#00A6FF",
  "Cancer 2" = "#8494FF",
  "Smooth Muscle" = "#C77CFF",
  "Adipocytes" = "#FF61C3",
  "Other" = "grey75"
)

proportion_plot_revised <- ggplot(
  plot_data,
  aes(x = PatientID, fill = CellType)
) +
  geom_bar(
    position = "fill",
    color = "white",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = celltype_colors,
    breaks = celltype_order,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.01),
    expand = c(0, 0)
  ) +
  ggtitle("Cell Type Proportions per Patient Sample") +
  labs(
    x = "Patient Sample",
    y = "Proportion of Cells",
    fill = ""
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 10)),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 10)),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm")
  ) +
  guides(
    fill = guide_legend(nrow = 3, byrow = TRUE)
  )

print(proportion_plot_revised)


### 2B

library(Seurat)
library(ggplot2)
library(dplyr)

fibroblast_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS"
)

DefaultAssay(fibroblast_subset) <- "RNA"

fibroblast_subset$CAF_type <- case_when(
  as.character(fibroblast_subset$seurat_clusters) == "0" ~ "Fibroblasts",
  as.character(fibroblast_subset$seurat_clusters) == "1" ~ "myCAF",
  as.character(fibroblast_subset$seurat_clusters) == "2" ~ "iCAF",
  as.character(fibroblast_subset$seurat_clusters) == "3" ~ "Pancreatic stellate",
  as.character(fibroblast_subset$seurat_clusters) == "4" ~ "apCAF",
  TRUE ~ "Other"
)

fibroblast_subset$CAF_type <- factor(
  fibroblast_subset$CAF_type,
  levels = c(
    "myCAF",
    "iCAF",
    "Pancreatic stellate",
    "Fibroblasts",
    "apCAF"
  )
)

caf_markers_clean <- c(
  "MMP11", "POSTN", "COL1A1", "COL11A1",
  "IGFBP5", "CCL2", "APOE", "C2orf40",
  "ACTA2",
  "COL4A1", "COL4A2",
  "COL1A2", "PDGFRA", "C3"
)

caf_markers_present <- caf_markers_clean[
  caf_markers_clean %in% rownames(fibroblast_subset)
]

missing_genes <- setdiff(caf_markers_clean, caf_markers_present)

print("Genes included in DotPlot:")
print(caf_markers_present)

print("Missing genes:")
print(missing_genes)

caf_dotplot_ordered <- DotPlot(
  fibroblast_subset,
  features = caf_markers_present,
  group.by = "CAF_type"
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "navy") +
  labs(
    x = "Marker genes",
    y = "CAF subtype"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 11, color = "black"),
    axis.text.y = element_text(face = "bold", size = 11, color = "black"),
    axis.title = element_text(face = "bold", size = 13),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  )

print(caf_dotplot_ordered)


### 2C: Fibroblast UMAP Split by Sex 


library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)

fibroblast_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS"
)

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

fibroblast_subset$Sex <- unname(
  sex_mapping[as.character(fibroblast_subset$orig.ident)]
)

fibroblast_subset$CAF_cluster <- as.character(fibroblast_subset$seurat_clusters)

fibroblast_subset$CAF_cluster <- factor(
  fibroblast_subset$CAF_cluster,
  levels = c("0", "1", "2", "3", "4")
)

caf_legend_labels <- c(
  "0" = "0 - Fibroblasts",
  "1" = "1 - myCAF",
  "2" = "2 - iCAF",
  "3" = "3 - Pancreatic stellate",
  "4" = "4 - apCAF"
)

caf_colors <- hue_pal()(5)
names(caf_colors) <- c("0", "1", "2", "3", "4")

female_fibros <- subset(
  fibroblast_subset,
  subset = Sex == "Female"
)

male_fibros <- subset(
  fibroblast_subset,
  subset = Sex == "Male"
)

female_umap <- DimPlot(
  female_fibros,
  reduction = "sub_umap",
  group.by = "CAF_cluster",
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 0.5
) +
  scale_color_manual(
    values = caf_colors,
    breaks = names(caf_legend_labels),
    labels = caf_legend_labels,
    drop = FALSE
  ) +
  labs(
    title = "Female",
    x = "umap_1",
    y = "umap_2"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )

male_umap <- DimPlot(
  male_fibros,
  reduction = "sub_umap",
  group.by = "CAF_cluster",
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 0.5
) +
  scale_color_manual(
    values = caf_colors,
    breaks = names(caf_legend_labels),
    labels = caf_legend_labels,
    drop = FALSE
  ) +
  labs(
    title = "Male",
    x = "umap_1",
    y = "umap_2"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )

combined_fibro_plot <- female_umap + male_umap +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.title = element_blank()
  )

print(combined_fibro_plot)


### 2D: Part 2 Fibroblast Cell Type Percentage with male versus female Samples

library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)

fibroblast_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS"
)

fibroblast_subset$CAF_type <- case_when(
  as.character(fibroblast_subset$seurat_clusters) == "0" ~ "Fibroblasts",
  as.character(fibroblast_subset$seurat_clusters) == "1" ~ "myCAF",
  as.character(fibroblast_subset$seurat_clusters) == "2" ~ "iCAF",
  as.character(fibroblast_subset$seurat_clusters) == "3" ~ "Pancreatic Stellate",
  as.character(fibroblast_subset$seurat_clusters) == "4" ~ "apCAF",
  TRUE ~ "Other"
)

caf_order <- c(
  "Fibroblasts",
  "myCAF",
  "iCAF",
  "Pancreatic Stellate",
  "apCAF"
)

fibroblast_subset$CAF_type <- factor(
  fibroblast_subset$CAF_type,
  levels = caf_order
)

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

patient_mapping <- c(
  "Sample_01" = "P01",
  "Sample_02" = "P02",
  "Sample_03" = "P03",
  "Sample_04" = "P04",
  "Sample_05" = "P05",
  "Sample_06" = "P06",
  "Sample_07" = "P07",
  "Sample_08" = "P08",
  "Sample_09" = "P09",
  "Sample_10" = "P10",
  "Sample_11" = "P11",
  "Sample_12" = "P12",
  "Sample_13" = "P13",
  "Sample_14" = "P14",
  "Sample_15" = "P15",
  "Sample_16" = "P16",
  "Sample_17" = "P17"
)

plot_data <- fibroblast_subset@meta.data

plot_data$CellType <- fibroblast_subset$CAF_type
plot_data$Sex <- unname(sex_mapping[as.character(plot_data$orig.ident)])
plot_data$PatientID <- unname(patient_mapping[as.character(plot_data$orig.ident)])

plot_data$CellType <- factor(
  plot_data$CellType,
  levels = caf_order
)

plot_data$PatientID <- factor(
  plot_data$PatientID,
  levels = sprintf("P%02d", 1:17)
)

plot_data <- plot_data %>%
  filter(!is.na(Sex), !is.na(PatientID), !is.na(CellType))

caf_colors <- hue_pal()(length(caf_order))
names(caf_colors) <- caf_order

plot_data_female <- plot_data %>%
  filter(Sex == "Female")

plot_data_male <- plot_data %>%
  filter(Sex == "Male")

female_caf_plot <- ggplot(
  plot_data_female,
  aes(x = PatientID, fill = CellType)
) +
  geom_bar(
    position = "fill",
    color = "white",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = caf_colors,
    breaks = caf_order,
    drop = FALSE
  ) +
  scale_x_discrete(drop = TRUE) +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.01),
    expand = c(0, 0)
  ) +
  labs(
    title = "Female",
    x = "Patient Sample",
    y = "Proportion",
    fill = ""
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )

male_caf_plot <- ggplot(
  plot_data_male,
  aes(x = PatientID, fill = CellType)
) +
  geom_bar(
    position = "fill",
    color = "white",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = caf_colors,
    breaks = caf_order,
    drop = FALSE
  ) +
  scale_x_discrete(drop = TRUE) +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.01),
    expand = c(0, 0)
  ) +
  labs(
    title = "Male",
    x = "Patient Sample",
    y = "Proportion",
    fill = ""
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )

caf_sex_comparison_plot <- female_caf_plot + male_caf_plot +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

print(caf_sex_comparison_plot)


### 3A: Volcano Plot of DGE of fibroblast populations Male versus Female.png

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)

fibroblast_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS"
)

DefaultAssay(fibroblast_subset) <- "RNA"

fibroblast_subset <- tryCatch(
  JoinLayers(fibroblast_subset),
  error = function(e) fibroblast_subset
)

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

sex_vector <- unname(
  sex_mapping[as.character(fibroblast_subset$orig.ident)]
)

names(sex_vector) <- colnames(fibroblast_subset)

fibroblast_subset <- AddMetaData(
  object = fibroblast_subset,
  metadata = sex_vector,
  col.name = "Sex"
)

print(colnames(fibroblast_subset@meta.data))
print(table(fibroblast_subset@meta.data$Sex, useNA = "ifany"))

fibroblast_subset <- subset(
  fibroblast_subset,
  subset = Sex %in% c("Male", "Female")
)

sex_markers <- FindMarkers(
  object = fibroblast_subset,
  ident.1 = "Male",
  ident.2 = "Female",
  group.by = "Sex",
  logfc.threshold = 0.25,
  min.pct = 0.10
)

sex_markers$Gene <- rownames(sex_markers)

sex_markers$Significance <- "Not Significant"

p_val_cutoff <- 0.05
log2fc_cutoff <- 0.5

sex_markers$Significance[
  sex_markers$p_val_adj < p_val_cutoff &
    sex_markers$avg_log2FC > log2fc_cutoff
] <- "Male"

sex_markers$Significance[
  sex_markers$p_val_adj < p_val_cutoff &
    sex_markers$avg_log2FC < -log2fc_cutoff
] <- "Female"

sex_markers$Significance <- factor(
  sex_markers$Significance,
  levels = c("Female", "Not Significant", "Male")
)

smallest_nonzero_padj <- min(
  sex_markers$p_val_adj[sex_markers$p_val_adj > 0],
  na.rm = TRUE
)

sex_markers <- sex_markers %>%
  mutate(
    p_val_adj_safe = ifelse(p_val_adj == 0, smallest_nonzero_padj, p_val_adj),
    neg_log10_padj = -log10(p_val_adj_safe)
  )

top_genes <- sex_markers %>%
  filter(Significance != "Not Significant") %>%
  group_by(Significance) %>%
  slice_min(order_by = p_val_adj_safe, n = 15) %>%
  ungroup()

print("======================================================")
print("DIFFERENTIAL EXPRESSION: MALE VS FEMALE FIBROBLASTS")
print("======================================================")
print(table(sex_markers$Significance))

write.csv(
  sex_markers,
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Male_vs_Female_Fibroblasts_DGE.csv",
  row.names = FALSE
)

volcano_plot <- ggplot(
  sex_markers,
  aes(x = avg_log2FC, y = neg_log10_padj, color = Significance)
) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(
    values = c(
      "Female" = "#D62728",
      "Not Significant" = "grey80",
      "Male" = "#0072B2"
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
    data = top_genes,
    aes(label = Gene),
    size = 4,
    fontface = "italic",
    color = "black",
    box.padding = 0.5,
    max.overlaps = 20
  ) +
  labs(
    title = "Differential Gene Expression",
    x = expression(Log[2]~"Fold Change"),
    y = expression(-Log[10]~"(Adjusted P-Value)"),
    color = ""
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )

print(volcano_plot)





### 3B Part 2: 

library(Seurat)
library(dplyr)
library(ggplot2)
library(msigdbr)
library(fgsea)
library(tibble)
library(tidyr)
library(patchwork)

fibroblast_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS"
)

DefaultAssay(fibroblast_subset) <- "RNA"

fibroblast_subset <- tryCatch(
  JoinLayers(fibroblast_subset),
  error = function(e) fibroblast_subset
)

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

sex_vector <- unname(
  sex_mapping[as.character(fibroblast_subset$orig.ident)]
)

names(sex_vector) <- colnames(fibroblast_subset)

fibroblast_subset <- AddMetaData(
  object = fibroblast_subset,
  metadata = sex_vector,
  col.name = "Sex"
)

fibroblast_subset <- subset(
  fibroblast_subset,
  subset = Sex %in% c("Male", "Female")
)

print(table(fibroblast_subset$Sex))

sex_markers <- FindMarkers(
  object = fibroblast_subset,
  ident.1 = "Male",
  ident.2 = "Female",
  group.by = "Sex",
  logfc.threshold = 0,
  min.pct = 0,
  only.pos = FALSE
)

sex_markers$Gene <- rownames(sex_markers)

sex_markers_clean <- sex_markers %>%
  filter(!is.na(p_val_adj)) %>%
  filter(!is.na(avg_log2FC))

smallest_nonzero_padj <- min(
  sex_markers_clean$p_val_adj[sex_markers_clean$p_val_adj > 0],
  na.rm = TRUE
)

sex_markers_clean <- sex_markers_clean %>%
  mutate(
    p_val_adj_safe = ifelse(
      p_val_adj == 0,
      smallest_nonzero_padj,
      p_val_adj
    ),
    rank_score = sign(avg_log2FC) * -log10(p_val_adj_safe)
  ) %>%
  filter(!is.na(rank_score), is.finite(rank_score)) %>%
  arrange(desc(abs(rank_score))) %>%
  group_by(Gene) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(desc(rank_score))

gene_ranks <- sex_markers_clean$rank_score
names(gene_ranks) <- sex_markers_clean$Gene
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

cat("Number of genes used for fgsea:", length(gene_ranks), "\n")

hallmark_sets <- tryCatch(
  msigdbr(
    species = "Homo sapiens",
    category = "H"
  ),
  error = function(e) {
    msigdbr(
      species = "Homo sapiens",
      collection = "H"
    )
  }
)

pathways_list <- split(
  x = hallmark_sets$gene_symbol,
  f = hallmark_sets$gs_name
)

set.seed(42)

fgsea_results <- fgsea(
  pathways = pathways_list,
  stats = gene_ranks,
  minSize = 10,
  maxSize = 500,
  eps = 0
)

nonzero_padj_gsea <- fgsea_results$padj[
  fgsea_results$padj > 0
]

smallest_nonzero_gsea_padj <- ifelse(
  length(nonzero_padj_gsea) > 0,
  min(nonzero_padj_gsea, na.rm = TRUE),
  .Machine$double.xmin
)

fgsea_results <- fgsea_results %>%
  arrange(padj) %>%
  mutate(
    pathway_clean = gsub("HALLMARK_", "", pathway),
    pathway_clean = gsub("_", " ", pathway_clean),
    Direction = ifelse(
      NES > 0,
      "Male",
      "Female"
    ),
    padj_safe = ifelse(
      padj == 0,
      smallest_nonzero_gsea_padj,
      padj
    ),
    neg_log10_FDR = -log10(padj_safe)
  )

print(head(fgsea_results, 20))

gsea_outdir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis/GSEA_Figures_NoExpressionFilter"
dir.create(gsea_outdir, showWarnings = FALSE, recursive = TRUE)

fgsea_results_save <- fgsea_results %>%
  as.data.frame() %>%
  mutate(
    leadingEdge = sapply(leadingEdge, paste, collapse = ";")
  )

write.csv(
  fgsea_results_save,
  file.path(
    gsea_outdir,
    "Fibroblast_Male_vs_Female_fgsea_NoExpressionFilter_Results.csv"
  ),
  row.names = FALSE
)

plot_gsea <- fgsea_results %>%
  filter(!is.na(NES), !is.na(padj)) %>%
  filter(padj < 0.25) %>%
  arrange(NES)

if (nrow(plot_gsea) < 5) {
  plot_gsea <- fgsea_results %>%
    filter(!is.na(NES), !is.na(padj)) %>%
    arrange(padj) %>%
    head(25) %>%
    arrange(NES)
}

plot_gsea$pathway_clean <- factor(
  plot_gsea$pathway_clean,
  levels = plot_gsea$pathway_clean
)

gsea_nes_dotplot <- ggplot(
  plot_gsea,
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
    title = NULL,
    x = "Normalized Enrichment Score",
    y = "Hallmark Pathway",
    size = "Gene Set Size"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    axis.title = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

print(gsea_nes_dotplot)

sig_gsea <- fgsea_results %>%
  filter(!is.na(NES), !is.na(padj))

sig_gsea_for_bar <- sig_gsea %>%
  filter(padj < 0.25)

if (nrow(sig_gsea_for_bar) < 10) {
  sig_gsea_for_bar <- sig_gsea
}

top_male <- sig_gsea_for_bar %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 10)

top_female <- sig_gsea_for_bar %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  slice_head(n = 10)

top_combined <- bind_rows(
  top_female,
  top_male
) %>%
  arrange(NES)

top_combined$pathway_clean <- factor(
  top_combined$pathway_clean,
  levels = top_combined$pathway_clean
)

gsea_barplot <- ggplot(
  top_combined,
  aes(x = pathway_clean, y = NES, fill = Direction)
) +
  geom_col(width = 0.75, color = "black", linewidth = 0.25) +
  coord_flip() +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  scale_fill_manual(
    values = c(
      "Female" = "#8B0000",
      "Male" = "#0072B2"
    )
  ) +
  labs(
    title = NULL,
    x = "Hallmark Pathway",
    y = "Normalized Enrichment Score",
    fill = ""
  ) +
  theme_classic(base_family = "sans") +
  theme(
    axis.title = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, face = "bold")
  )

print(gsea_barplot)


### 3C: DGE schwann cells Male verus Female


library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)

true_schwann_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/true_schwann_ERBB3_PMP22_subset.rds"
)

DefaultAssay(true_schwann_subset) <- "RNA"

true_schwann_subset <- tryCatch(
  JoinLayers(true_schwann_subset),
  error = function(e) true_schwann_subset
)

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

true_schwann_subset$Sex <- unname(
  sex_mapping[as.character(true_schwann_subset$orig.ident)]
)

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

sex_counts <- table(true_schwann_subset$Sex)

if (!all(c("Male", "Female") %in% names(sex_counts))) {
  stop("Both Male and Female Schwann cells are required for DGE.")
}

if (any(sex_counts[c("Male", "Female")] < 3)) {
  stop("One sex has fewer than 3 Schwann cells. DGE is not reliable or may fail.")
}

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
  warning("Very few genes passed the expression > 0.5 filter. DGE may be underpowered.")
}

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

cat("Number of genes returned by FindMarkers:", nrow(schwann_sex_markers), "\n")

schwann_sex_markers$Significance <- "Not Significant"

p_val_cutoff <- 0.05
log2fc_cutoff <- 0.5

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

smallest_nonzero_padj <- min(
  schwann_sex_markers$p_val_adj[schwann_sex_markers$p_val_adj > 0],
  na.rm = TRUE
)

schwann_sex_markers <- schwann_sex_markers %>%
  mutate(
    p_val_adj_safe = ifelse(p_val_adj == 0, smallest_nonzero_padj, p_val_adj),
    neg_log10_padj = -log10(p_val_adj_safe)
  )

top_schwann_genes <- schwann_sex_markers %>%
  filter(Significance != "Not Significant") %>%
  group_by(Significance) %>%
  slice_min(order_by = p_val_adj_safe, n = 15) %>%
  ungroup()

write.csv(
  schwann_sex_markers,
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Schwann_Male_vs_Female_DGE_ExpressionGT0.5.csv",
  row.names = FALSE
)

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
      "Upregulated in Females" = "Female",
      "Not Significant" = "Not Significant",
      "Upregulated in Males" = "Male"
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
  labs(
    title = "Differential Gene Expression",
    x = expression(Log[2]~"Fold Change"),
    y = expression(-Log[10]~"(Adjusted P-Value)"),
    color = ""
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



### 3D: Part 2 GSEA schwann cells Male Versus Female

library(Seurat)
library(dplyr)
library(ggplot2)
library(msigdbr)
library(fgsea)
library(tibble)
library(tidyr)

true_schwann_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/true_schwann_ERBB3_PMP22_subset.rds"
)

DefaultAssay(true_schwann_subset) <- "RNA"

true_schwann_subset <- tryCatch(
  JoinLayers(true_schwann_subset),
  error = function(e) true_schwann_subset
)

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

true_schwann_subset$Sex <- unname(
  sex_mapping[as.character(true_schwann_subset$orig.ident)]
)

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

sex_counts <- table(true_schwann_subset$Sex)

if (!all(c("Male", "Female") %in% names(sex_counts))) {
  stop("Both Male and Female Schwann cells are required for DGE.")
}

if (any(sex_counts[c("Male", "Female")] < 3)) {
  stop("One sex has fewer than 3 Schwann cells. DGE is not reliable or may fail.")
}

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

if (length(genes_expr_gt_0.5) < 50) {
  warning("Very few genes passed the expression > 0.5 filter. GSEA may be underpowered.")
}

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

cat("Number of genes returned by FindMarkers:", nrow(schwann_sex_markers), "\n")

smallest_nonzero_padj <- min(
  schwann_sex_markers$p_val_adj[schwann_sex_markers$p_val_adj > 0],
  na.rm = TRUE
)

schwann_sex_markers <- schwann_sex_markers %>%
  mutate(
    p_val_adj_safe = ifelse(p_val_adj == 0, smallest_nonzero_padj, p_val_adj),
    rank_score = sign(avg_log2FC) * -log10(p_val_adj_safe)
  ) %>%
  filter(!is.na(rank_score), is.finite(rank_score)) %>%
  arrange(desc(rank_score))

schwann_gene_ranks <- schwann_sex_markers$rank_score
names(schwann_gene_ranks) <- schwann_sex_markers$Gene

schwann_gene_ranks <- schwann_gene_ranks[!duplicated(names(schwann_gene_ranks))]
schwann_gene_ranks <- sort(schwann_gene_ranks, decreasing = TRUE)

cat("Number of Schwann genes going into GSEA:", length(schwann_gene_ranks), "\n")

hallmark_sets <- tryCatch(
  msigdbr(
    species = "Homo sapiens",
    category = "H"
  ),
  error = function(e) {
    msigdbr(
      species = "Homo sapiens",
      collection = "H"
    )
  }
)

pathways_list <- split(
  x = hallmark_sets$gene_symbol,
  f = hallmark_sets$gs_name
)

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
      "Male",
      "Female"
    )
  )

print(head(schwann_fgsea_results, 20))

schwann_outdir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis/Schwann_DGE_GSEA_Figures"
dir.create(schwann_outdir, showWarnings = FALSE, recursive = TRUE)

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
  scale_fill_manual(
    values = c(
      "Female" = "#8B0000",
      "Male" = "#0072B2"
    )
  ) +
  labs(
    title = NULL,
    x = "Hallmark Pathway",
    y = "Normalized Enrichment Score",
    fill = ""
  ) +
  theme_classic(base_family = "sans") +
  theme(
    axis.title = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, face = "bold")
  )

print(schwann_gsea_barplot)

ggsave(
  file.path(schwann_outdir, "Figure_3_Schwann_GSEA_Barplot_NoTitle.pdf"),
  plot = schwann_gsea_barplot,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(schwann_outdir, "Figure_3_Schwann_GSEA_Barplot_NoTitle.png"),
  plot = schwann_gsea_barplot,
  width = 9,
  height = 7,
  dpi = 300
)


### DGE of apCAF

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)

fibroblast_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS"
)

DefaultAssay(fibroblast_subset) <- "RNA"

fibroblast_subset <- tryCatch(
  JoinLayers(fibroblast_subset),
  error = function(e) fibroblast_subset
)

apcaf_subset <- subset(
  fibroblast_subset,
  subset = as.character(seurat_clusters) == "4"
)

saveRDS(
  apcaf_subset,
  file = "/Users/sudharshiniram/Desktop/SRI Data Analysis/apCAF_Subset.RDS"
)

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

apcaf_subset$Sex <- unname(
  sex_mapping[as.character(apcaf_subset$orig.ident)]
)

apcaf_subset <- subset(
  apcaf_subset,
  subset = Sex %in% c("Male", "Female")
)

apcaf_subset$Sex <- factor(
  apcaf_subset$Sex,
  levels = c("Female", "Male")
)

print("apCAF cells by sex:")
print(table(apcaf_subset$Sex))

print("apCAF cells by patient and sex:")
print(table(apcaf_subset$orig.ident, apcaf_subset$Sex))

sex_counts <- table(apcaf_subset$Sex)

if (!all(c("Male", "Female") %in% names(sex_counts))) {
  stop("Both Male and Female apCAF cells are required for DGE.")
}

if (any(sex_counts[c("Male", "Female")] < 3)) {
  stop("One sex has fewer than 3 apCAF cells. DGE is not reliable or may fail.")
}

apcaf_markers <- FindMarkers(
  object = apcaf_subset,
  ident.1 = "Male",
  ident.2 = "Female",
  group.by = "Sex",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.05,
  only.pos = FALSE
)

apcaf_markers$Gene <- rownames(apcaf_markers)

p_val_cutoff <- 0.05
log2fc_cutoff <- 0.5

apcaf_markers$Significance <- "Not Significant"

apcaf_markers$Significance[
  apcaf_markers$p_val_adj < p_val_cutoff &
    apcaf_markers$avg_log2FC > log2fc_cutoff
] <- "Male"

apcaf_markers$Significance[
  apcaf_markers$p_val_adj < p_val_cutoff &
    apcaf_markers$avg_log2FC < -log2fc_cutoff
] <- "Female"

apcaf_markers$Significance <- factor(
  apcaf_markers$Significance,
  levels = c("Female", "Not Significant", "Male")
)

smallest_nonzero_padj <- min(
  apcaf_markers$p_val_adj[apcaf_markers$p_val_adj > 0],
  na.rm = TRUE
)

apcaf_markers <- apcaf_markers %>%
  mutate(
    p_val_adj_safe = ifelse(p_val_adj == 0, smallest_nonzero_padj, p_val_adj),
    neg_log10_padj = -log10(p_val_adj_safe)
  )

top_40_genes <- apcaf_markers %>%
  filter(Significance != "Not Significant") %>%
  arrange(p_val_adj_safe, desc(abs(avg_log2FC))) %>%
  slice_head(n = 40)

if (nrow(top_40_genes) < 40) {
  top_40_genes <- apcaf_markers %>%
    arrange(p_val_adj_safe, desc(abs(avg_log2FC))) %>%
    slice_head(n = 40)
}

print("Top 40 most differentially expressed apCAF genes:")
print(top_40_genes[, c("Gene", "avg_log2FC", "p_val_adj", "pct.1", "pct.2")])

write.csv(
  apcaf_markers,
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/apCAF_Male_vs_Female_DGE.csv",
  row.names = FALSE
)

apcaf_volcano <- ggplot(
  apcaf_markers,
  aes(x = avg_log2FC, y = neg_log10_padj, color = Significance)
) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(
    values = c(
      "Female" = "#8B0000",
      "Not Significant" = "grey80",
      "Male" = "#0072B2"
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
    data = top_40_genes,
    aes(label = Gene),
    size = 3.6,
    fontface = "italic",
    color = "black",
    box.padding = 0.5,
    point.padding = 0.4,
    max.overlaps = 40,
    min.segment.length = 0
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.03, 0.18))
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = NULL,
    x = expression(Log[2]~"Fold Change"),
    y = expression(-Log[10]~"(Adjusted P-Value)"),
    color = ""
  ) +
  theme_classic(base_family = "sans") +
  theme(
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, face = "bold"),
    plot.margin = margin(t = 25, r = 25, b = 20, l = 20)
  )

print(apcaf_volcano)

ggsave(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/apCAF_Male_vs_Female_Top40_DGE_Volcano.pdf",
  plot = apcaf_volcano,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/apCAF_Male_vs_Female_Top40_DGE_Volcano.png",
  plot = apcaf_volcano,
  width = 9,
  height = 7,
  dpi = 300
)


### apCAF only genes DGE


library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)

fibroblast_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS"
)

DefaultAssay(fibroblast_subset) <- "RNA"

fibroblast_subset <- tryCatch(
  JoinLayers(fibroblast_subset),
  error = function(e) fibroblast_subset
)

apcaf_subset <- subset(
  fibroblast_subset,
  subset = as.character(seurat_clusters) == "4"
)

saveRDS(
  apcaf_subset,
  file = "/Users/sudharshiniram/Desktop/SRI Data Analysis/apCAF_Subset.RDS"
)

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

apcaf_subset$Sex <- unname(
  sex_mapping[as.character(apcaf_subset$orig.ident)]
)

apcaf_subset <- subset(
  apcaf_subset,
  subset = Sex %in% c("Male", "Female")
)

apcaf_subset$Sex <- factor(
  apcaf_subset$Sex,
  levels = c("Female", "Male")
)

print("apCAF cells by sex:")
print(table(apcaf_subset$Sex))

print("apCAF cells by patient and sex:")
print(table(apcaf_subset$orig.ident, apcaf_subset$Sex))

sex_counts <- table(apcaf_subset$Sex)

if (!all(c("Male", "Female") %in% names(sex_counts))) {
  stop("Both Male and Female apCAF cells are required for DGE.")
}

if (any(sex_counts[c("Male", "Female")] < 3)) {
  stop("One sex has fewer than 3 apCAF cells. DGE is not reliable or may fail.")
}

apcaf_gene_set <- unique(c(
  "CD74",
  "CIITA",
  "HLA-DRA",
  "HLA-DRB1",
  "HLA-DRB3",
  "HLA-DRB4",
  "HLA-DRB5",
  "HLA-DPA1",
  "HLA-DPB1",
  "HLA-DQA1",
  "HLA-DQA2",
  "HLA-DQB1",
  "HLA-DQB2",
  "HLA-DMA",
  "HLA-DMB",
  "HLA-DOA",
  "HLA-DOB",
  "IFI30",
  "CTSB",
  "CTSS",
  "CTSL",
  "LAMP1",
  "LAMP2",
  "B2M",
  "HLA-A",
  "HLA-B",
  "HLA-C",
  "TAP1",
  "TAP2",
  "TAPBP",
  "CALR",
  "CANX",
  "PDIA3",
  "PSMB8",
  "PSMB9",
  "PSMB10",
  "NLRC5",
  "IRF1",
  "STAT1",
  "GBP1",
  "GBP2",
  "GBP4",
  "GBP5",
  "IFIT1",
  "IFIT2",
  "IFIT3",
  "ISG15",
  "MX1",
  "OAS1",
  "OAS2",
  "OAS3",
  "CXCL9",
  "CXCL10",
  "CXCL11"
))

apcaf_genes_present <- apcaf_gene_set[apcaf_gene_set %in% rownames(apcaf_subset)]
missing_apcaf_genes <- setdiff(apcaf_gene_set, apcaf_genes_present)

print("apCAF genes included in DGE:")
print(apcaf_genes_present)

print("apCAF genes missing from dataset:")
print(missing_apcaf_genes)

if (length(apcaf_genes_present) < 2) {
  stop("Fewer than 2 apCAF genes were found in the dataset.")
}

apcaf_markers <- FindMarkers(
  object = apcaf_subset,
  ident.1 = "Male",
  ident.2 = "Female",
  group.by = "Sex",
  features = apcaf_genes_present,
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0,
  only.pos = FALSE
)

apcaf_markers$Gene <- rownames(apcaf_markers)

p_val_cutoff <- 0.05
log2fc_cutoff <- 0.5

apcaf_markers$Significance <- "Not Significant"

apcaf_markers$Significance[
  apcaf_markers$p_val_adj < p_val_cutoff &
    apcaf_markers$avg_log2FC > log2fc_cutoff
] <- "Male"

apcaf_markers$Significance[
  apcaf_markers$p_val_adj < p_val_cutoff &
    apcaf_markers$avg_log2FC < -log2fc_cutoff
] <- "Female"

apcaf_markers$Significance <- factor(
  apcaf_markers$Significance,
  levels = c("Female", "Not Significant", "Male")
)

smallest_nonzero_padj <- min(
  apcaf_markers$p_val_adj[apcaf_markers$p_val_adj > 0],
  na.rm = TRUE
)

apcaf_markers <- apcaf_markers %>%
  mutate(
    p_val_adj_safe = ifelse(p_val_adj == 0, smallest_nonzero_padj, p_val_adj),
    neg_log10_padj = -log10(p_val_adj_safe)
  )

top_40_apcaf_genes <- apcaf_markers %>%
  arrange(p_val_adj_safe, desc(abs(avg_log2FC))) %>%
  slice_head(n = min(40, nrow(apcaf_markers)))

print("Top apCAF genes shown on volcano plot:")
print(top_40_apcaf_genes[, c("Gene", "avg_log2FC", "p_val_adj", "pct.1", "pct.2")])

write.csv(
  apcaf_markers,
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/apCAF_OnlyGenes_Male_vs_Female_DGE.csv",
  row.names = FALSE
)

apcaf_volcano <- ggplot(
  apcaf_markers,
  aes(x = avg_log2FC, y = neg_log10_padj, color = Significance)
) +
  geom_point(alpha = 0.85, size = 2) +
  scale_color_manual(
    values = c(
      "Female" = "#8B0000",
      "Not Significant" = "grey80",
      "Male" = "#0072B2"
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
    data = top_40_apcaf_genes,
    aes(label = Gene),
    size = 3.7,
    fontface = "italic",
    color = "black",
    box.padding = 0.5,
    point.padding = 0.4,
    max.overlaps = 40,
    min.segment.length = 0
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.03, 0.18))
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = NULL,
    x = expression(Log[2]~"Fold Change"),
    y = expression(-Log[10]~"(Adjusted P-Value)"),
    color = ""
  ) +
  theme_classic(base_family = "sans") +
  theme(
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, face = "bold"),
    plot.margin = margin(t = 25, r = 25, b = 20, l = 20)
  )

print(apcaf_volcano)

ggsave(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/apCAF_OnlyGenes_Male_vs_Female_Top40_Volcano.pdf",
  plot = apcaf_volcano,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/apCAF_OnlyGenes_Male_vs_Female_Top40_Volcano.png",
  plot = apcaf_volcano,
  width = 9,
  height = 7,
  dpi = 300
)



## doesn't plot nonapCAF genes

apcaf_markers_for_plot <- apcaf_markers

top_40_genes <- apcaf_markers %>%
  filter(!(Gene %in% non_apcaf_genes_to_exclude)) %>%
  arrange(p_val_adj_safe, desc(abs(avg_log2FC))) %>%
  slice_head(n = 40)

if (nrow(top_40_genes) < 40) {
  warning(
    paste(
      "Only",
      nrow(top_40_genes),
      "genes available after removing excluded genes. Fewer than 40 will be labeled."
    )
  )
}

print("Top 40 plotted genes after removing excluded genes:")
print(top_40_genes[, c("Gene", "avg_log2FC", "p_val_adj", "pct.1", "pct.2")])

print("Excluded genes that would have appeared but were removed from the top 40 labels:")
print(
  apcaf_markers %>%
    filter(Gene %in% non_apcaf_genes_to_exclude) %>%
    arrange(p_val_adj_safe, desc(abs(avg_log2FC))) %>%
    dplyr::select(Gene, avg_log2FC, p_val_adj, pct.1, pct.2)
)

write.csv(
  apcaf_markers,
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/apCAF_Male_vs_Female_DGE_AllGenes.csv",
  row.names = FALSE
)

write.csv(
  top_40_genes,
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/apCAF_Male_vs_Female_Top40_FilteredGenes.csv",
  row.names = FALSE
)

apcaf_volcano <- ggplot(
  apcaf_markers_for_plot,
  aes(x = avg_log2FC, y = neg_log10_padj, color = Significance)
) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(
    values = c(
      "Female" = "#8B0000",
      "Not Significant" = "grey80",
      "Male" = "#0072B2"
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
    data = top_40_genes,
    aes(label = Gene),
    size = 3.6,
    fontface = "italic",
    color = "black",
    box.padding = 0.5,
    point.padding = 0.4,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.03, 0.18))
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = NULL,
    x = expression(Log[2]~"Fold Change"),
    y = expression(-Log[10]~"(Adjusted P-Value)"),
    color = ""
  ) +
  theme_classic(base_family = "sans") +
  theme(
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, face = "bold"),
    plot.margin = margin(t = 25, r = 25, b = 20, l = 20)
  )

print(apcaf_volcano)

ggsave(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/apCAF_Male_vs_Female_Top40_DGE_Volcano_Filtered.pdf",
  plot = apcaf_volcano,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/apCAF_Male_vs_Female_Top40_DGE_Volcano_Filtered.png",
  plot = apcaf_volcano,
  width = 9,
  height = 7,
  dpi = 300
)


### Dotplot apCAFs

library(Seurat)
library(ggplot2)

fibroblast_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS"
)

DefaultAssay(fibroblast_subset) <- "RNA"

fibroblast_subset$CAF_type <- case_when(
  as.character(fibroblast_subset$seurat_clusters) == "0" ~ "Fibroblasts",
  as.character(fibroblast_subset$seurat_clusters) == "1" ~ "myCAF",
  as.character(fibroblast_subset$seurat_clusters) == "2" ~ "iCAF",
  as.character(fibroblast_subset$seurat_clusters) == "3" ~ "Pancreatic stellate",
  as.character(fibroblast_subset$seurat_clusters) == "4" ~ "apCAF",
  TRUE ~ "Other"
)

fibroblast_subset$CAF_type <- factor(
  fibroblast_subset$CAF_type,
  levels = rev(c(
    "Fibroblasts",
    "myCAF",
    "iCAF",
    "Pancreatic stellate",
    "apCAF"
  ))
)

apcaf_dotplot_genes <- c(
  "CD74",
  "HLA-DRA",
  "CIITA",
  "CTSS",
  "LAMP2",
  "TAP1",
  "TAP2"
)


apcaf_dotplot_genes_present <- apcaf_dotplot_genes[
  apcaf_dotplot_genes %in% rownames(fibroblast_subset)
]

print("Genes included:")
print(apcaf_dotplot_genes_present)

print("Missing genes:")
print(setdiff(apcaf_dotplot_genes, apcaf_dotplot_genes_present))

apcaf_dotplot <- DotPlot(
  fibroblast_subset,
  features = apcaf_dotplot_genes_present,
  group.by = "CAF_type"
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "navy") +
  labs(
    x = "apCAF marker genes",
    y = "CAF subtype"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 11, color = "black"),
    axis.text.y = element_text(face = "bold", size = 11, color = "black"),
    axis.title = element_text(face = "bold", size = 13),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  )

print(apcaf_dotplot)



### 2B: Fibroblast Dotplot (with new apCAF markers)

library(Seurat)
library(ggplot2)
library(dplyr)

fibroblast_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS"
)

DefaultAssay(fibroblast_subset) <- "RNA"

fibroblast_subset$CAF_type <- case_when(
  as.character(fibroblast_subset$seurat_clusters) == "0" ~ "Fibroblasts",
  as.character(fibroblast_subset$seurat_clusters) == "1" ~ "myCAFs",
  as.character(fibroblast_subset$seurat_clusters) == "2" ~ "iCAFs",
  as.character(fibroblast_subset$seurat_clusters) == "3" ~ "PSCs",
  as.character(fibroblast_subset$seurat_clusters) == "4" ~ "apCAFs",
  TRUE ~ "Other"
)

fibroblast_subset$CAF_type <- factor(
  fibroblast_subset$CAF_type,
  levels = c(
    "Fibroblasts",
    "myCAFs",
    "iCAFs",
    "apCAFs",
    "PSCs"
  )
)

caf_markers_clean <- c(
  "IGFBP5",
  "MMP11",
  "COL1A2",
  "COL1A1",
  "POSTN",
  "COL11A1",
  "ACTA2",
  "COL4A1",
  "COL4A2",
  "PDGFRA",
  "C3",
  "CD74",
  "HLA-DRA",
  "CTSS",
  "LAMP2",
  "VIT",
  "MEOX2",
  "APOD"
)

caf_markers_present <- caf_markers_clean[
  caf_markers_clean %in% rownames(fibroblast_subset)
]

missing_genes <- setdiff(caf_markers_clean, caf_markers_present)

print("Genes included in DotPlot:")
print(caf_markers_present)

print("Missing genes:")
print(missing_genes)

caf_dotplot_ordered <- DotPlot(
  fibroblast_subset,
  features = caf_markers_present,
  group.by = "CAF_type",
  cols = c("lightgrey", "navy")
) +
  coord_flip() +
  scale_size_continuous(
    name = "Percent Expressed",
    breaks = c(10, 15, 25, 50, 75),
    range = c(0, 7)
  ) +
  labs(
    x = "CAF subtype",
    y = "Marker genes",
    color = "Average Expression"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 11,
      color = "black"
    ),
    axis.text.y = element_text(
      face = "italic",
      size = 11,
      color = "black"
    ),
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    legend.text = element_text(size = 9)
  )

print(caf_dotplot_ordered)





### 3G: CAF sustypes male versus female RATIOS

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

fibroblast_subset <- readRDS(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS"
)

DefaultAssay(fibroblast_subset) <- "RNA"

# Add sex metadata
sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

fibroblast_subset$Sex <- unname(
  sex_mapping[as.character(fibroblast_subset$orig.ident)]
)

fibroblast_subset <- subset(
  fibroblast_subset,
  subset = Sex %in% c("Male", "Female")
)

fibroblast_subset$Sex <- factor(
  fibroblast_subset$Sex,
  levels = c("Female", "Male")
)

# Add CAF subtype labels
fibroblast_subset$CAF_type <- case_when(
  as.character(fibroblast_subset$seurat_clusters) == "0" ~ "Fibroblasts",
  as.character(fibroblast_subset$seurat_clusters) == "1" ~ "myCAF",
  as.character(fibroblast_subset$seurat_clusters) == "2" ~ "iCAF",
  as.character(fibroblast_subset$seurat_clusters) == "3" ~ "PSC",
  as.character(fibroblast_subset$seurat_clusters) == "4" ~ "apCAF",
  TRUE ~ "Other"
)

fibroblast_subset$CAF_type <- factor(
  fibroblast_subset$CAF_type,
  levels = c("Fibroblasts", "myCAF", "iCAF", "PSC", "apCAF")
)

# Make metadata dataframe
meta_df <- as.data.frame(fibroblast_subset@meta.data)

# Count iCAFs and myCAFs by sex
icaf_mycaf_counts <- meta_df %>%
  dplyr::filter(CAF_type %in% c("iCAF", "myCAF")) %>%
  dplyr::count(Sex, CAF_type, name = "Cell_Count") %>%
  tidyr::complete(
    Sex,
    CAF_type,
    fill = list(Cell_Count = 0)
  ) %>%
  dplyr::arrange(Sex, CAF_type)

print("iCAF and myCAF counts by sex:")
print(icaf_mycaf_counts)

# Convert to wide format
icaf_mycaf_wide <- icaf_mycaf_counts %>%
  tidyr::pivot_wider(
    names_from = CAF_type,
    values_from = Cell_Count,
    values_fill = 0
  )

# Calculate iCAF:myCAF ratio
icaf_mycaf_ratio_table <- icaf_mycaf_wide %>%
  dplyr::mutate(
    iCAF_to_myCAF_ratio = ifelse(myCAF == 0, NA, iCAF / myCAF),
    count_ratio = paste0(iCAF, ":", myCAF),
    plot_label = paste0(
      "iCAF:myCAF = ", count_ratio,
      "\nRatio = ", round(iCAF_to_myCAF_ratio, 2)
    )
  ) %>%
  dplyr::select(
    Sex,
    iCAF,
    myCAF,
    count_ratio,
    iCAF_to_myCAF_ratio,
    plot_label
  )

print("iCAF to myCAF ratio table:")
print(icaf_mycaf_ratio_table)

# Bar graph: iCAF:myCAF ratio in female versus male
icaf_mycaf_ratio_barplot <- ggplot(
  icaf_mycaf_ratio_table,
  aes(x = Sex, y = iCAF_to_myCAF_ratio, fill = Sex)
) +
  geom_col(
    width = 0.6,
    color = "black",
    linewidth = 0.4
  ) +
  geom_text(
    aes(label = plot_label),
    vjust = -0.4,
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#8B0000",
      "Male" = "#0072B2"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    title = "iCAF:myCAF Ratio by Sex",
    x = "Sex",
    y = "iCAF:myCAF ratio"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    axis.text = element_text(
      face = "bold",
      size = 11,
      color = "black"
    ),
    legend.position = "none"
  )

print(icaf_mycaf_ratio_barplot)



### 4D: CellChat on CAF populations and schwann cell interactions  (iCAF, myCAF, apCAF) 
# ============================================================
# REVISED CELLCHAT CIRCLE PLOTS:
# Properly labeled All Samples and Female vs Male comparison
# ============================================================

# Output folder
cellchat_outdir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis/CellChat_CAF_Subtypes_Schwann_Circle"

dir.create(
  cellchat_outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

# Groups to include in circle plots
target_groups <- c(
  "myCAFs",
  "iCAFs",
  "apCAFs",
  "Schwann cells"
)

# ============================================================
# 1. HELPER FUNCTION: EXTRACT COUNT OR WEIGHT MATRIX
# ============================================================

get_cellchat_matrix <- function(
    cellchat_obj,
    measure = "count",
    remove_self_loops = FALSE
) {
  
  if (measure == "count") {
    mat <- cellchat_obj@net$count
  } else if (measure == "weight") {
    mat <- cellchat_obj@net$weight
  } else {
    stop("measure must be either 'count' or 'weight'")
  }
  
  keep_groups <- target_groups[target_groups %in% rownames(mat)]
  
  if (length(keep_groups) < 2) {
    stop("Could not find at least two target groups in the CellChat matrix.")
  }
  
  mat_use <- mat[keep_groups, keep_groups, drop = FALSE]
  
  if (remove_self_loops) {
    diag(mat_use) <- 0
  }
  
  return(mat_use)
}

# ============================================================
# 2. HELPER FUNCTION: GET NODE SIZES
# ============================================================

get_group_size <- function(cellchat_obj, mat_use) {
  
  group_size <- as.numeric(table(cellchat_obj@idents)[rownames(mat_use)])
  names(group_size) <- rownames(mat_use)
  
  group_size[is.na(group_size)] <- 1
  
  return(group_size)
}

# ============================================================
# 3. HELPER FUNCTION: PLOT ONE CIRCLE
# ============================================================

plot_one_cellchat_circle <- function(
    cellchat_obj,
    panel_title,
    measure = "count",
    remove_self_loops = FALSE,
    shared_edge_max = NULL
) {
  
  mat_use <- get_cellchat_matrix(
    cellchat_obj = cellchat_obj,
    measure = measure,
    remove_self_loops = remove_self_loops
  )
  
  group_size <- get_group_size(
    cellchat_obj = cellchat_obj,
    mat_use = mat_use
  )
  
  if (is.null(shared_edge_max)) {
    edge_max <- max(mat_use, na.rm = TRUE)
  } else {
    edge_max <- shared_edge_max
  }
  
  if (!is.finite(edge_max) || edge_max == 0) {
    edge_max <- 1
  }
  
  netVisual_circle(
    mat_use,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = TRUE,
    edge.weight.max = edge_max,
    edge.width.max = 8,
    vertex.label.cex = 1.3,
    title.name = panel_title
  )
}

# ============================================================
# 4. SAVE SINGLE ALL-SAMPLES CIRCLE PLOT
# ============================================================

save_all_samples_circle_plot <- function(
    cellchat_obj,
    measure = "count",
    remove_self_loops = FALSE,
    filename_prefix
) {
  
  if (measure == "count") {
    general_title <- "CellChat Network: Number of Ligand-Receptor Interactions"
  } else {
    general_title <- "CellChat Network: Interaction Strength"
  }
  
  if (remove_self_loops) {
    general_title <- paste0(general_title, "\nCross-Talk Only")
  }
  
  # PDF
  pdf(
    file.path(cellchat_outdir, paste0(filename_prefix, ".pdf")),
    width = 8,
    height = 8
  )
  
  par(
    oma = c(0, 0, 3.5, 0),
    mar = c(1, 1, 3, 1)
  )
  
  plot_one_cellchat_circle(
    cellchat_obj = cellchat_obj,
    panel_title = "All Samples",
    measure = measure,
    remove_self_loops = remove_self_loops
  )
  
  mtext(
    general_title,
    outer = TRUE,
    side = 3,
    line = 1,
    cex = 1.4,
    font = 2
  )
  
  dev.off()
  
  # PNG
  png(
    file.path(cellchat_outdir, paste0(filename_prefix, ".png")),
    width = 2400,
    height = 2400,
    res = 300
  )
  
  par(
    oma = c(0, 0, 3.5, 0),
    mar = c(1, 1, 3, 1)
  )
  
  plot_one_cellchat_circle(
    cellchat_obj = cellchat_obj,
    panel_title = "All Samples",
    measure = measure,
    remove_self_loops = remove_self_loops
  )
  
  mtext(
    general_title,
    outer = TRUE,
    side = 3,
    line = 1,
    cex = 1.4,
    font = 2
  )
  
  dev.off()
}

# ============================================================
# 5. SAVE FEMALE VS MALE CIRCLE COMPARISON PLOT
# ============================================================

save_female_male_circle_comparison <- function(
    cellchat_female,
    cellchat_male,
    measure = "count",
    remove_self_loops = FALSE,
    filename_prefix
) {
  
  if (measure == "count") {
    general_title <- "CellChat Network: Number of Ligand-Receptor Interactions"
  } else {
    general_title <- "CellChat Network: Interaction Strength"
  }
  
  if (remove_self_loops) {
    general_title <- paste0(general_title, "\nCross-Talk Only")
  }
  
  # Use same edge scale for Female and Male so line widths are comparable
  female_mat <- get_cellchat_matrix(
    cellchat_obj = cellchat_female,
    measure = measure,
    remove_self_loops = remove_self_loops
  )
  
  male_mat <- get_cellchat_matrix(
    cellchat_obj = cellchat_male,
    measure = measure,
    remove_self_loops = remove_self_loops
  )
  
  shared_edge_max <- max(
    c(female_mat, male_mat),
    na.rm = TRUE
  )
  
  if (!is.finite(shared_edge_max) || shared_edge_max == 0) {
    shared_edge_max <- 1
  }
  
  # PDF
  pdf(
    file.path(cellchat_outdir, paste0(filename_prefix, ".pdf")),
    width = 14,
    height = 8
  )
  
  par(
    mfrow = c(1, 2),
    oma = c(0, 0, 4.5, 0),
    mar = c(1, 1, 4, 1)
  )
  
  plot_one_cellchat_circle(
    cellchat_obj = cellchat_female,
    panel_title = "Female",
    measure = measure,
    remove_self_loops = remove_self_loops,
    shared_edge_max = shared_edge_max
  )
  
  plot_one_cellchat_circle(
    cellchat_obj = cellchat_male,
    panel_title = "Male",
    measure = measure,
    remove_self_loops = remove_self_loops,
    shared_edge_max = shared_edge_max
  )
  
  mtext(
    general_title,
    outer = TRUE,
    side = 3,
    line = 1.2,
    cex = 1.6,
    font = 2
  )
  
  dev.off()
  
  # PNG
  png(
    file.path(cellchat_outdir, paste0(filename_prefix, ".png")),
    width = 4200,
    height = 2400,
    res = 300
  )
  
  par(
    mfrow = c(1, 2),
    oma = c(0, 0, 4.5, 0),
    mar = c(1, 1, 4, 1)
  )
  
  plot_one_cellchat_circle(
    cellchat_obj = cellchat_female,
    panel_title = "Female",
    measure = measure,
    remove_self_loops = remove_self_loops,
    shared_edge_max = shared_edge_max
  )
  
  plot_one_cellchat_circle(
    cellchat_obj = cellchat_male,
    panel_title = "Male",
    measure = measure,
    remove_self_loops = remove_self_loops,
    shared_edge_max = shared_edge_max
  )
  
  mtext(
    general_title,
    outer = TRUE,
    side = 3,
    line = 1.2,
    cex = 1.6,
    font = 2
  )
  
  dev.off()
}

# ============================================================
# 6. SAVE ALL SAMPLES CIRCLE PLOTS
# ============================================================

save_all_samples_circle_plot(
  cellchat_obj = cellchat_all,
  measure = "count",
  remove_self_loops = FALSE,
  filename_prefix = "AllSamples_CAF_Subtypes_Schwann_Circle_Count"
)

save_all_samples_circle_plot(
  cellchat_obj = cellchat_all,
  measure = "weight",
  remove_self_loops = FALSE,
  filename_prefix = "AllSamples_CAF_Subtypes_Schwann_Circle_Strength"
)

save_all_samples_circle_plot(
  cellchat_obj = cellchat_all,
  measure = "count",
  remove_self_loops = TRUE,
  filename_prefix = "AllSamples_CAF_Subtypes_Schwann_Circle_Count_CrossTalkOnly"
)

save_all_samples_circle_plot(
  cellchat_obj = cellchat_all,
  measure = "weight",
  remove_self_loops = TRUE,
  filename_prefix = "AllSamples_CAF_Subtypes_Schwann_Circle_Strength_CrossTalkOnly"
)

# ============================================================
# 7. SAVE FEMALE VS MALE COMPARISON CIRCLE PLOTS
# ============================================================

save_female_male_circle_comparison(
  cellchat_female = cellchat_female,
  cellchat_male = cellchat_male,
  measure = "count",
  remove_self_loops = FALSE,
  filename_prefix = "Female_vs_Male_CAF_Subtypes_Schwann_Circle_Count"
)

save_female_male_circle_comparison(
  cellchat_female = cellchat_female,
  cellchat_male = cellchat_male,
  measure = "weight",
  remove_self_loops = FALSE,
  filename_prefix = "Female_vs_Male_CAF_Subtypes_Schwann_Circle_Strength"
)

save_female_male_circle_comparison(
  cellchat_female = cellchat_female,
  cellchat_male = cellchat_male,
  measure = "count",
  remove_self_loops = TRUE,
  filename_prefix = "Female_vs_Male_CAF_Subtypes_Schwann_Circle_Count_CrossTalkOnly"
)

save_female_male_circle_comparison(
  cellchat_female = cellchat_female,
  cellchat_male = cellchat_male,
  measure = "weight",
  remove_self_loops = TRUE,
  filename_prefix = "Female_vs_Male_CAF_Subtypes_Schwann_Circle_Strength_CrossTalkOnly"
)

# ============================================================
# 8. FINAL MESSAGE
# ============================================================

cat("\n======================================================\n")
cat("REVISED CELLCHAT CIRCLE PLOTS SAVED\n")
cat("Folder:\n")
cat(cellchat_outdir, "\n")
cat("Saved plot types:\n")
cat("- All Samples: count and strength\n")
cat("- All Samples cross-talk only: count and strength\n")
cat("- Female vs Male comparison: count and strength\n")
cat("- Female vs Male cross-talk only: count and strength\n")
cat("======================================================\n")


### 4D: CellChat on CAF subtypes and schwann interactions without edge labels

# ============================================================
# FULL STANDALONE SCRIPT:
# CELLCHAT CIRCLE PLOTS WITHOUT EDGE NUMBERS
# Titles = "Female" or "Male"
# Prevents overlapping/cut-off labels
# myCAFs, iCAFs, apCAFs, and Schwann cells
# ============================================================

# ============================================================
# 0. LOAD REQUIRED LIBRARIES
# ============================================================

library(CellChat)
library(Seurat)
library(dplyr)

options(stringsAsFactors = FALSE)

# Close any open plotting devices
graphics.off()

# ============================================================
# 1. SET INPUT AND OUTPUT FOLDER
# ============================================================

cellchat_outdir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis/CellChat_CAF_Subtypes_Schwann_Circle"

dir.create(
  cellchat_outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# 2. READ SAVED CELLCHAT OBJECTS
# ============================================================

cellchat_female <- readRDS(
  file.path(cellchat_outdir, "cellchat_female_myCAF_iCAF_apCAF_schwann.rds")
)

cellchat_male <- readRDS(
  file.path(cellchat_outdir, "cellchat_male_myCAF_iCAF_apCAF_schwann.rds")
)

cellchat_female <- updateCellChat(cellchat_female)
cellchat_male <- updateCellChat(cellchat_male)

# ============================================================
# 3. GROUPS TO INCLUDE IN CIRCLE PLOTS
# ============================================================

target_groups <- c(
  "myCAFs",
  "iCAFs",
  "apCAFs",
  "Schwann cells"
)

print("Groups available in Female CellChat object:")
print(levels(cellchat_female@idents))

print("Groups available in Male CellChat object:")
print(levels(cellchat_male@idents))

# ============================================================
# 4. FIXED LAYOUT TO PREVENT LABEL OVERLAP
# ============================================================
# This manually places the nodes farther apart and away from edges.

circle_layout <- matrix(
  c(
    0.70,  0.00,   # myCAFs: right
    0.00,  0.75,   # iCAFs: top
    -0.70,  0.00,   # apCAFs: left
    0.00, -0.75    # Schwann cells: bottom
  ),
  ncol = 2,
  byrow = TRUE
)

rownames(circle_layout) <- target_groups

# ============================================================
# 5. HELPER FUNCTION: EXTRACT COUNT OR WEIGHT MATRIX
# ============================================================

get_cellchat_matrix <- function(
    cellchat_obj,
    measure = "count",
    remove_self_loops = FALSE
) {
  
  if (measure == "count") {
    mat <- cellchat_obj@net$count
  } else if (measure == "weight") {
    mat <- cellchat_obj@net$weight
  } else {
    stop("measure must be either 'count' or 'weight'")
  }
  
  keep_groups <- target_groups[target_groups %in% rownames(mat)]
  
  if (length(keep_groups) < 2) {
    stop("Could not find at least two target groups in the CellChat matrix.")
  }
  
  mat_use <- mat[keep_groups, keep_groups, drop = FALSE]
  
  if (remove_self_loops) {
    diag(mat_use) <- 0
  }
  
  return(mat_use)
}

# ============================================================
# 6. HELPER FUNCTION: GET NODE SIZES
# ============================================================

get_group_size <- function(cellchat_obj, mat_use) {
  
  group_size <- as.numeric(table(cellchat_obj@idents)[rownames(mat_use)])
  names(group_size) <- rownames(mat_use)
  
  group_size[is.na(group_size)] <- 1
  
  return(group_size)
}

# ============================================================
# 7. HELPER FUNCTION: GET SHARED EDGE SCALE
# ============================================================

get_shared_edge_max <- function(
    cellchat_female,
    cellchat_male,
    measure = "count",
    remove_self_loops = FALSE
) {
  
  female_mat <- get_cellchat_matrix(
    cellchat_obj = cellchat_female,
    measure = measure,
    remove_self_loops = remove_self_loops
  )
  
  male_mat <- get_cellchat_matrix(
    cellchat_obj = cellchat_male,
    measure = measure,
    remove_self_loops = remove_self_loops
  )
  
  shared_edge_max <- max(
    c(female_mat, male_mat),
    na.rm = TRUE
  )
  
  if (!is.finite(shared_edge_max) || shared_edge_max == 0) {
    shared_edge_max <- 1
  }
  
  return(shared_edge_max)
}

# ============================================================
# 8. HELPER FUNCTION: PLOT ONE CIRCLE WITHOUT EDGE NUMBERS
# ============================================================

plot_one_cellchat_circle <- function(
    cellchat_obj,
    title_name,
    measure = "count",
    remove_self_loops = FALSE,
    shared_edge_max = NULL
) {
  
  mat_use <- get_cellchat_matrix(
    cellchat_obj = cellchat_obj,
    measure = measure,
    remove_self_loops = remove_self_loops
  )
  
  group_size <- get_group_size(
    cellchat_obj = cellchat_obj,
    mat_use = mat_use
  )
  
  if (is.null(shared_edge_max)) {
    edge_max <- max(mat_use, na.rm = TRUE)
  } else {
    edge_max <- shared_edge_max
  }
  
  if (!is.finite(edge_max) || edge_max == 0) {
    edge_max <- 1
  }
  
  layout_use <- circle_layout[rownames(mat_use), , drop = FALSE]
  
  netVisual_circle(
    mat_use,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.label.cex = 0,
    edge.weight.max = edge_max,
    edge.width.max = 7,
    vertex.label.cex = 1.1,
    title.name = title_name,
    layout = layout_use,
    margin = 0.35
  )
}

# ============================================================
# 9. HELPER FUNCTION: SAVE ONE INDIVIDUAL CIRCLE PLOT
# ============================================================

save_individual_circle_plot <- function(
    cellchat_obj,
    title_name,
    measure = "count",
    remove_self_loops = FALSE,
    filename_prefix,
    shared_edge_max = NULL
) {
  
  # PDF version
  pdf(
    file.path(cellchat_outdir, paste0(filename_prefix, ".pdf")),
    width = 9.5,
    height = 9
  )
  
  par(
    oma = c(0, 0, 0, 0),
    mar = c(1, 1, 4, 1),
    xpd = TRUE
  )
  
  plot_one_cellchat_circle(
    cellchat_obj = cellchat_obj,
    title_name = title_name,
    measure = measure,
    remove_self_loops = remove_self_loops,
    shared_edge_max = shared_edge_max
  )
  
  dev.off()
  
  # PNG version
  png(
    file.path(cellchat_outdir, paste0(filename_prefix, ".png")),
    width = 3000,
    height = 2850,
    res = 300
  )
  
  par(
    oma = c(0, 0, 0, 0),
    mar = c(1, 1, 4, 1),
    xpd = TRUE
  )
  
  plot_one_cellchat_circle(
    cellchat_obj = cellchat_obj,
    title_name = title_name,
    measure = measure,
    remove_self_loops = remove_self_loops,
    shared_edge_max = shared_edge_max
  )
  
  dev.off()
}

# ============================================================
# 10. SHARED EDGE SCALES FOR COMPARABLE FEMALE AND MALE PLOTS
# ============================================================

shared_count_max <- get_shared_edge_max(
  cellchat_female = cellchat_female,
  cellchat_male = cellchat_male,
  measure = "count",
  remove_self_loops = FALSE
)

shared_weight_max <- get_shared_edge_max(
  cellchat_female = cellchat_female,
  cellchat_male = cellchat_male,
  measure = "weight",
  remove_self_loops = FALSE
)

shared_count_crosstalk_max <- get_shared_edge_max(
  cellchat_female = cellchat_female,
  cellchat_male = cellchat_male,
  measure = "count",
  remove_self_loops = TRUE
)

shared_weight_crosstalk_max <- get_shared_edge_max(
  cellchat_female = cellchat_female,
  cellchat_male = cellchat_male,
  measure = "weight",
  remove_self_loops = TRUE
)

# ============================================================
# 11. SAVE FEMALE AND MALE COUNT CIRCLE PLOTS
# ============================================================

save_individual_circle_plot(
  cellchat_obj = cellchat_female,
  title_name = "Female",
  measure = "count",
  remove_self_loops = FALSE,
  filename_prefix = "Female_CAF_Subtypes_Schwann_Circle_Count_Final",
  shared_edge_max = shared_count_max
)

save_individual_circle_plot(
  cellchat_obj = cellchat_male,
  title_name = "Male",
  measure = "count",
  remove_self_loops = FALSE,
  filename_prefix = "Male_CAF_Subtypes_Schwann_Circle_Count_Final",
  shared_edge_max = shared_count_max
)

# ============================================================
# 12. SAVE FEMALE AND MALE STRENGTH CIRCLE PLOTS
# ============================================================

save_individual_circle_plot(
  cellchat_obj = cellchat_female,
  title_name = "Female",
  measure = "weight",
  remove_self_loops = FALSE,
  filename_prefix = "Female_CAF_Subtypes_Schwann_Circle_Strength_Final",
  shared_edge_max = shared_weight_max
)

save_individual_circle_plot(
  cellchat_obj = cellchat_male,
  title_name = "Male",
  measure = "weight",
  remove_self_loops = FALSE,
  filename_prefix = "Male_CAF_Subtypes_Schwann_Circle_Strength_Final",
  shared_edge_max = shared_weight_max
)

# ============================================================
# 13. SAVE FEMALE AND MALE CROSS-TALK ONLY COUNT PLOTS
# ============================================================
# These remove self-loops like iCAF -> iCAF, myCAF -> myCAF, etc.

save_individual_circle_plot(
  cellchat_obj = cellchat_female,
  title_name = "Female",
  measure = "count",
  remove_self_loops = TRUE,
  filename_prefix = "Female_CAF_Subtypes_Schwann_Circle_Count_CrossTalkOnly_Final",
  shared_edge_max = shared_count_crosstalk_max
)

save_individual_circle_plot(
  cellchat_obj = cellchat_male,
  title_name = "Male",
  measure = "count",
  remove_self_loops = TRUE,
  filename_prefix = "Male_CAF_Subtypes_Schwann_Circle_Count_CrossTalkOnly_Final",
  shared_edge_max = shared_count_crosstalk_max
)

# ============================================================
# 14. SAVE FEMALE AND MALE CROSS-TALK ONLY STRENGTH PLOTS
# ============================================================

save_individual_circle_plot(
  cellchat_obj = cellchat_female,
  title_name = "Female",
  measure = "weight",
  remove_self_loops = TRUE,
  filename_prefix = "Female_CAF_Subtypes_Schwann_Circle_Strength_CrossTalkOnly_Final",
  shared_edge_max = shared_weight_crosstalk_max
)

save_individual_circle_plot(
  cellchat_obj = cellchat_male,
  title_name = "Male",
  measure = "weight",
  remove_self_loops = TRUE,
  filename_prefix = "Male_CAF_Subtypes_Schwann_Circle_Strength_CrossTalkOnly_Final",
  shared_edge_max = shared_weight_crosstalk_max
)

# ============================================================
# 15. FINAL MESSAGE
# ============================================================

cat("\n======================================================\n")
cat("FINAL CELLCHAT CIRCLE PLOTS SAVED\n")
cat("Folder:\n")
cat(cellchat_outdir, "\n")
cat("Updates applied:\n")
cat("- Titles added to every graph: Female or Male\n")
cat("- Edge numbers removed\n")
cat("- Manual node layout added to prevent label overlap\n")
cat("- Larger plotting canvas added to prevent cut-off labels\n")
cat("======================================================\n")


### Running T Test Statistics on General UMAP

# ============================================================
# GENERAL UMAP CLUSTER COMPOSITION ANALYSIS
# Cell count + proportion differences by sex
# T-test per cluster/cell type
# No individual dots/points on plots
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

options(stringsAsFactors = FALSE)

# ============================================================
# 1. LOAD OBJECT
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

Merged_Joined <- readRDS(
  file.path(base_dir, "Merged_Joined.RDS")
)

DefaultAssay(Merged_Joined) <- "RNA"

Merged_Joined <- tryCatch(
  JoinLayers(Merged_Joined),
  error = function(e) Merged_Joined
)

# ============================================================
# 2. ADD SEX METADATA
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

Merged_Joined$Sex <- unname(
  sex_mapping[as.character(Merged_Joined$orig.ident)]
)

Merged_Joined$Sex <- factor(
  Merged_Joined$Sex,
  levels = c("Female", "Male")
)

# ============================================================
# 3. CREATE GENERAL UMAP CLUSTER/CELL TYPE LABELS
# ============================================================

Merged_Joined$general_celltype <- dplyr::case_when(
  as.character(Merged_Joined$seurat_clusters) == "0"  ~ "Cancer 1",
  as.character(Merged_Joined$seurat_clusters) == "1"  ~ "Fibroblast",
  as.character(Merged_Joined$seurat_clusters) == "2"  ~ "Ductal",
  as.character(Merged_Joined$seurat_clusters) == "3"  ~ "B Cells",
  as.character(Merged_Joined$seurat_clusters) == "4"  ~ "Myeloid",
  as.character(Merged_Joined$seurat_clusters) == "5"  ~ "Acinar",
  as.character(Merged_Joined$seurat_clusters) == "6"  ~ "Endocrine",
  as.character(Merged_Joined$seurat_clusters) == "7"  ~ "Endothelial",
  as.character(Merged_Joined$seurat_clusters) == "8"  ~ "T / NK Cells",
  as.character(Merged_Joined$seurat_clusters) == "9"  ~ "Cancer 2",
  as.character(Merged_Joined$seurat_clusters) == "10" ~ "Smooth Muscle",
  as.character(Merged_Joined$seurat_clusters) == "11" ~ "Schwann Cells",
  as.character(Merged_Joined$seurat_clusters) == "12" ~ "Adipocytes",
  TRUE ~ "Other"
)

celltype_order <- c(
  "Cancer 1",
  "Cancer 2",
  "Ductal",
  "Acinar",
  "Endocrine",
  "Fibroblast",
  "Schwann Cells",
  "Myeloid",
  "B Cells",
  "T / NK Cells",
  "Endothelial",
  "Smooth Muscle",
  "Adipocytes",
  "Other"
)

Merged_Joined$general_celltype <- factor(
  Merged_Joined$general_celltype,
  levels = celltype_order
)

# ============================================================
# 4. CREATE OUTPUT FOLDER
# ============================================================

composition_outdir <- file.path(
  base_dir,
  "General_UMAP_Cluster_Composition_Ttests"
)

dir.create(
  composition_outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# 5. COUNT CELLS PER SAMPLE PER CLUSTER
# ============================================================
# Each patient/sample is treated as the statistical unit.

cell_counts_per_sample <- Merged_Joined@meta.data %>%
  filter(
    !is.na(Sex),
    !is.na(orig.ident),
    !is.na(general_celltype)
  ) %>%
  dplyr::count(
    orig.ident,
    Sex,
    general_celltype,
    name = "Cell_Count"
  ) %>%
  complete(
    nesting(orig.ident, Sex),
    general_celltype = factor(celltype_order, levels = celltype_order),
    fill = list(Cell_Count = 0)
  ) %>%
  mutate(
    general_celltype = factor(
      general_celltype,
      levels = celltype_order
    )
  )

sample_totals <- cell_counts_per_sample %>%
  group_by(orig.ident, Sex) %>%
  summarise(
    Total_Cells_Per_Sample = sum(Cell_Count),
    .groups = "drop"
  )

composition_per_sample <- cell_counts_per_sample %>%
  left_join(
    sample_totals,
    by = c("orig.ident", "Sex")
  ) %>%
  mutate(
    Percent_of_Sample = 100 * Cell_Count / Total_Cells_Per_Sample
  )

write.csv(
  composition_per_sample,
  file.path(composition_outdir, "General_UMAP_Cell_Counts_and_Proportions_Per_Sample.csv"),
  row.names = FALSE
)

# ============================================================
# 6. T-TEST FUNCTION
# ============================================================

run_t_test_by_cluster <- function(data, value_col) {
  
  results <- data %>%
    group_by(general_celltype) %>%
    group_modify(~{
      
      female_values <- .x[[value_col]][.x$Sex == "Female"]
      male_values <- .x[[value_col]][.x$Sex == "Male"]
      
      p_val <- tryCatch(
        {
          if (
            length(na.omit(female_values)) >= 2 &&
            length(na.omit(male_values)) >= 2
          ) {
            t.test(female_values, male_values)$p.value
          } else {
            NA_real_
          }
        },
        error = function(e) NA_real_
      )
      
      tibble(
        Female_Mean = mean(female_values, na.rm = TRUE),
        Male_Mean = mean(male_values, na.rm = TRUE),
        Female_SD = sd(female_values, na.rm = TRUE),
        Male_SD = sd(male_values, na.rm = TRUE),
        Female_N = length(na.omit(female_values)),
        Male_N = length(na.omit(male_values)),
        p_value = p_val
      )
    }) %>%
    ungroup() %>%
    mutate(
      p_adj_BH = p.adjust(p_value, method = "BH"),
      Significance = case_when(
        is.na(p_adj_BH) ~ "NA",
        p_adj_BH < 0.001 ~ "***",
        p_adj_BH < 0.01 ~ "**",
        p_adj_BH < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  
  return(results)
}

count_ttest_results <- run_t_test_by_cluster(
  data = composition_per_sample,
  value_col = "Cell_Count"
)

proportion_ttest_results <- run_t_test_by_cluster(
  data = composition_per_sample,
  value_col = "Percent_of_Sample"
)

write.csv(
  count_ttest_results,
  file.path(composition_outdir, "Ttest_Cell_Counts_Per_Cluster_By_Sex.csv"),
  row.names = FALSE
)

write.csv(
  proportion_ttest_results,
  file.path(composition_outdir, "Ttest_Proportions_Per_Cluster_By_Sex.csv"),
  row.names = FALSE
)

print("T-test results for raw cell counts:")
print(count_ttest_results)

print("T-test results for proportions:")
print(proportion_ttest_results)

# ============================================================
# 7. SUMMARY DATA FOR BAR PLOTS
# ============================================================

plot_summary <- composition_per_sample %>%
  group_by(general_celltype, Sex) %>%
  summarise(
    Mean_Cell_Count = mean(Cell_Count, na.rm = TRUE),
    SEM_Cell_Count = sd(Cell_Count, na.rm = TRUE) / sqrt(n()),
    Mean_Percent = mean(Percent_of_Sample, na.rm = TRUE),
    SEM_Percent = sd(Percent_of_Sample, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

write.csv(
  plot_summary,
  file.path(composition_outdir, "General_UMAP_Composition_Barplot_Summary.csv"),
  row.names = FALSE
)

# ============================================================
# 8. CREATE P-VALUE LABELS
# ============================================================

count_p_labels <- plot_summary %>%
  group_by(general_celltype) %>%
  summarise(
    y_position = max(Mean_Cell_Count + SEM_Cell_Count, na.rm = TRUE) * 1.2,
    .groups = "drop"
  ) %>%
  left_join(
    count_ttest_results %>%
      dplyr::select(general_celltype, p_value, p_adj_BH, Significance),
    by = "general_celltype"
  ) %>%
  mutate(
    y_position = ifelse(
      !is.finite(y_position) | y_position == 0,
      1,
      y_position
    ),
    label = ifelse(
      is.na(p_value),
      "p = NA",
      paste0("p = ", signif(p_value, 2))
    )
  )

proportion_p_labels <- plot_summary %>%
  group_by(general_celltype) %>%
  summarise(
    y_position = max(Mean_Percent + SEM_Percent, na.rm = TRUE) * 1.2,
    .groups = "drop"
  ) %>%
  left_join(
    proportion_ttest_results %>%
      dplyr::select(general_celltype, p_value, p_adj_BH, Significance),
    by = "general_celltype"
  ) %>%
  mutate(
    y_position = ifelse(
      !is.finite(y_position) | y_position == 0,
      0.05,
      y_position
    ),
    label = ifelse(
      is.na(p_value),
      "p = NA",
      paste0("p = ", signif(p_value, 2))
    )
  )

# ============================================================
# 9. RAW CELL COUNT BAR GRAPH
# ============================================================

count_plot <- ggplot(
  plot_summary,
  aes(
    x = Sex,
    y = Mean_Cell_Count,
    fill = Sex
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_Cell_Count - SEM_Cell_Count,
      ymax = Mean_Cell_Count + SEM_Cell_Count
    ),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_text(
    data = count_p_labels,
    aes(
      x = 1.5,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    size = 3.2,
    fontface = "bold"
  ) +
  facet_wrap(
    ~general_celltype,
    scales = "free_y",
    ncol = 4
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    title = "General UMAP Cluster Cell Counts by Sex",
    x = "Sex",
    y = "Number of cells per sample",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    legend.position = "none"
  )

print(count_plot)

ggsave(
  file.path(composition_outdir, "General_UMAP_Cell_Counts_By_Sex_Barplot_NoDots.png"),
  plot = count_plot,
  width = 13,
  height = 9,
  dpi = 300
)

ggsave(
  file.path(composition_outdir, "General_UMAP_Cell_Counts_By_Sex_Barplot_NoDots.pdf"),
  plot = count_plot,
  width = 13,
  height = 9,
  dpi = 300
)

# ============================================================
# 10. PROPORTION BAR GRAPH
# ============================================================

proportion_plot <- ggplot(
  plot_summary,
  aes(
    x = Sex,
    y = Mean_Percent,
    fill = Sex
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_Percent - SEM_Percent,
      ymax = Mean_Percent + SEM_Percent
    ),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_text(
    data = proportion_p_labels,
    aes(
      x = 1.5,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    size = 3.2,
    fontface = "bold"
  ) +
  facet_wrap(
    ~general_celltype,
    scales = "free_y",
    ncol = 4
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    title = "General UMAP Cluster Proportions by Sex",
    x = "Sex",
    y = "Percent of total cells per sample",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    legend.position = "none"
  )

print(proportion_plot)

ggsave(
  file.path(composition_outdir, "General_UMAP_Proportions_By_Sex_Barplot_NoDots.png"),
  plot = proportion_plot,
  width = 13,
  height = 9,
  dpi = 300
)

ggsave(
  file.path(composition_outdir, "General_UMAP_Proportions_By_Sex_Barplot_NoDots.pdf"),
  plot = proportion_plot,
  width = 13,
  height = 9,
  dpi = 300
)

# ============================================================
# 11. SIGNIFICANT CLUSTERS ONLY - PROPORTIONS
# ============================================================

significant_proportion_clusters <- proportion_ttest_results %>%
  filter(!is.na(p_value), p_value < 0.05) %>%
  pull(general_celltype)

if (length(significant_proportion_clusters) > 0) {
  
  proportion_plot_sig <- ggplot(
    plot_summary %>%
      filter(general_celltype %in% significant_proportion_clusters),
    aes(
      x = Sex,
      y = Mean_Percent,
      fill = Sex
    )
  ) +
    geom_col(
      color = "black",
      linewidth = 0.3,
      width = 0.65
    ) +
    geom_errorbar(
      aes(
        ymin = Mean_Percent - SEM_Percent,
        ymax = Mean_Percent + SEM_Percent
      ),
      width = 0.2,
      linewidth = 0.4
    ) +
    geom_text(
      data = proportion_p_labels %>%
        filter(general_celltype %in% significant_proportion_clusters),
      aes(
        x = 1.5,
        y = y_position,
        label = label
      ),
      inherit.aes = FALSE,
      size = 3.5,
      fontface = "bold"
    ) +
    facet_wrap(
      ~general_celltype,
      scales = "free_y",
      ncol = 3
    ) +
    scale_fill_manual(
      values = c(
        "Female" = "#F8766D",
        "Male" = "#00BFC4"
      )
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.25))
    ) +
    labs(
      title = "Significant General UMAP Proportion Differences by Sex",
      x = "Sex",
      y = "Percent of total cells per sample",
      fill = "Sex"
    ) +
    theme_classic(base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 10, color = "black"),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "none"
    )
  
  print(proportion_plot_sig)
  
  ggsave(
    file.path(composition_outdir, "Significant_General_UMAP_Proportion_Differences_By_Sex_NoDots.png"),
    plot = proportion_plot_sig,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  ggsave(
    file.path(composition_outdir, "Significant_General_UMAP_Proportion_Differences_By_Sex_NoDots.pdf"),
    plot = proportion_plot_sig,
    width = 10,
    height = 6,
    dpi = 300
  )
}

# ============================================================
# 12. SIGNIFICANT CLUSTERS ONLY - CELL COUNTS
# ============================================================

significant_count_clusters <- count_ttest_results %>%
  filter(!is.na(p_value), p_value < 0.05) %>%
  pull(general_celltype)

if (length(significant_count_clusters) > 0) {
  
  count_plot_sig <- ggplot(
    plot_summary %>%
      filter(general_celltype %in% significant_count_clusters),
    aes(
      x = Sex,
      y = Mean_Cell_Count,
      fill = Sex
    )
  ) +
    geom_col(
      color = "black",
      linewidth = 0.3,
      width = 0.65
    ) +
    geom_errorbar(
      aes(
        ymin = Mean_Cell_Count - SEM_Cell_Count,
        ymax = Mean_Cell_Count + SEM_Cell_Count
      ),
      width = 0.2,
      linewidth = 0.4
    ) +
    geom_text(
      data = count_p_labels %>%
        filter(general_celltype %in% significant_count_clusters),
      aes(
        x = 1.5,
        y = y_position,
        label = label
      ),
      inherit.aes = FALSE,
      size = 3.5,
      fontface = "bold"
    ) +
    facet_wrap(
      ~general_celltype,
      scales = "free_y",
      ncol = 3
    ) +
    scale_fill_manual(
      values = c(
        "Female" = "#F8766D",
        "Male" = "#00BFC4"
      )
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.25))
    ) +
    labs(
      title = "Significant General UMAP Cell Count Differences by Sex",
      x = "Sex",
      y = "Number of cells per sample",
      fill = "Sex"
    ) +
    theme_classic(base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 10, color = "black"),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "none"
    )
  
  print(count_plot_sig)
  
  ggsave(
    file.path(composition_outdir, "Significant_General_UMAP_Cell_Count_Differences_By_Sex_NoDots.png"),
    plot = count_plot_sig,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  ggsave(
    file.path(composition_outdir, "Significant_General_UMAP_Cell_Count_Differences_By_Sex_NoDots.pdf"),
    plot = count_plot_sig,
    width = 10,
    height = 6,
    dpi = 300
  )
}

# ============================================================
# 13. FINAL MESSAGE
# ============================================================

cat("\n======================================================\n")
cat("GENERAL UMAP CLUSTER COMPOSITION ANALYSIS COMPLETE\n")
cat("Results saved to:\n")
cat(composition_outdir, "\n\n")
cat("Main outputs:\n")
cat("- General_UMAP_Cell_Counts_By_Sex_Barplot_NoDots.png/pdf\n")
cat("- General_UMAP_Proportions_By_Sex_Barplot_NoDots.png/pdf\n")
cat("- Ttest_Cell_Counts_Per_Cluster_By_Sex.csv\n")
cat("- Ttest_Proportions_Per_Cluster_By_Sex.csv\n")
cat("- General_UMAP_Cell_Counts_and_Proportions_Per_Sample.csv\n")
cat("- General_UMAP_Composition_Barplot_Summary.csv\n")
cat("======================================================\n")



### T Test on Schwann-Fibroblast

# ============================================================
# SCHWANN-FIBROBLAST CELLCHAT T-TEST ANALYSIS
# Interaction count + interaction strength by sex
# Statistical unit = patient/sample
# No individual dots/points on plots
# ============================================================

library(CellChat)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

options(stringsAsFactors = FALSE)

# ============================================================
# 1. SET PATHS
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

merged_path <- file.path(base_dir, "Merged_Joined.RDS")

cellchat_outdir <- file.path(
  base_dir,
  "CellChat_Schwann_Fibroblast_PerSample_Ttests"
)

dir.create(
  cellchat_outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

per_sample_rds_dir <- file.path(
  cellchat_outdir,
  "Per_Sample_CellChat_RDS"
)

dir.create(
  per_sample_rds_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

summary_csv_path <- file.path(
  cellchat_outdir,
  "Schwann_Fibroblast_CellChat_PerSample_Count_Strength.csv"
)

# ============================================================
# 2. LOAD SEURAT OBJECT
# ============================================================

Merged_Joined <- readRDS(merged_path)

DefaultAssay(Merged_Joined) <- "RNA"

Merged_Joined <- tryCatch(
  JoinLayers(Merged_Joined),
  error = function(e) Merged_Joined
)

# ============================================================
# 3. ADD SEX METADATA
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

Merged_Joined$Sex <- unname(
  sex_mapping[as.character(Merged_Joined$orig.ident)]
)

Merged_Joined$Sex <- factor(
  Merged_Joined$Sex,
  levels = c("Female", "Male")
)

# ============================================================
# 4. CREATE CELLCHAT LABELS
# ============================================================

Merged_Joined$cellchat_celltype <- NA_character_

# Fibroblasts = broad cluster 1
Merged_Joined$cellchat_celltype[
  as.character(Merged_Joined$seurat_clusters) == "1"
] <- "Fibroblasts"

# True Schwann cells
if ("true_schwann_subset" %in% colnames(Merged_Joined@meta.data)) {
  
  Merged_Joined$cellchat_celltype[
    as.character(Merged_Joined$true_schwann_subset) %in% c(
      "True Schwann",
      "TRUE",
      "True",
      "true"
    )
  ] <- "Schwann cells"
  
} else {
  
  warning("true_schwann_subset not found. Using seurat cluster 11 as Schwann fallback.")
  
  Merged_Joined$cellchat_celltype[
    as.character(Merged_Joined$seurat_clusters) == "11"
  ] <- "Schwann cells"
}

target_groups <- c(
  "Fibroblasts",
  "Schwann cells"
)

target_cells <- subset(
  Merged_Joined,
  subset = cellchat_celltype %in% target_groups &
    Sex %in% c("Female", "Male")
)

target_cells$cellchat_celltype <- factor(
  target_cells$cellchat_celltype,
  levels = target_groups
)

target_cells$Sex <- factor(
  target_cells$Sex,
  levels = c("Female", "Male")
)

print("Cell counts by sample, sex, and CellChat group:")
print(table(target_cells$orig.ident, target_cells$Sex, target_cells$cellchat_celltype))

# ============================================================
# 5. HELPER FUNCTION TO GET RNA DATA
# ============================================================

get_rna_data <- function(seurat_obj) {
  
  data_input <- tryCatch(
    GetAssayData(
      seurat_obj,
      assay = "RNA",
      layer = "data"
    ),
    error = function(e) {
      GetAssayData(
        seurat_obj,
        assay = "RNA",
        slot = "data"
      )
    }
  )
  
  return(data_input)
}

# ============================================================
# 6. RUN CELLCHAT FOR ONE SAMPLE
# ============================================================

run_cellchat_one_sample <- function(seurat_obj, sample_id, min_cells_filter = 10) {
  
  print(paste("Running CellChat for", sample_id, "..."))
  
  DefaultAssay(seurat_obj) <- "RNA"
  
  meta <- seurat_obj@meta.data
  
  meta$cellchat_celltype <- droplevels(
    factor(
      meta$cellchat_celltype,
      levels = target_groups
    )
  )
  
  meta$samples <- as.character(meta$orig.ident)
  
  print(table(meta$cellchat_celltype))
  
  data.input <- get_rna_data(seurat_obj)
  
  cellchat <- createCellChat(
    object = data.input,
    meta = meta,
    group.by = "cellchat_celltype"
  )
  
  cellchat@DB <- CellChatDB.human
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(
    cellchat,
    type = "triMean",
    raw.use = TRUE
  )
  
  cellchat <- filterCommunication(
    cellchat,
    min.cells = min_cells_filter
  )
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  cellchat <- updateCellChat(cellchat)
  
  return(cellchat)
}

# ============================================================
# 7. RUN OR READ PER-SAMPLE CELLCHAT RESULTS
# ============================================================
# Important:
# A t-test needs patient/sample-level values.
# This does NOT t-test one pooled female CellChat object vs one pooled male object.
# It runs CellChat separately per sample, then compares sample-level
# interaction count and strength between Female and Male samples.

min_cells_required <- 10

if (file.exists(summary_csv_path)) {
  
  print("Reading saved per-sample CellChat summary:")
  print(summary_csv_path)
  
  per_sample_summary <- read.csv(summary_csv_path)
  
} else {
  
  sample_ids <- sort(unique(as.character(target_cells$orig.ident)))
  
  per_sample_summary_list <- list()
  
  for (sample_id in sample_ids) {
    
    sample_obj <- subset(
      target_cells,
      subset = orig.ident == sample_id
    )
    
    sample_sex <- unique(as.character(sample_obj$Sex))
    sample_sex <- sample_sex[!is.na(sample_sex)][1]
    
    group_counts <- table(sample_obj$cellchat_celltype)
    
    fibroblast_n <- ifelse(
      "Fibroblasts" %in% names(group_counts),
      as.numeric(group_counts["Fibroblasts"]),
      0
    )
    
    schwann_n <- ifelse(
      "Schwann cells" %in% names(group_counts),
      as.numeric(group_counts["Schwann cells"]),
      0
    )
    
    print("======================================================")
    print(paste("Sample:", sample_id))
    print(paste("Sex:", sample_sex))
    print(paste("Fibroblast cells:", fibroblast_n))
    print(paste("Schwann cells:", schwann_n))
    
    if (
      fibroblast_n < min_cells_required ||
      schwann_n < min_cells_required
    ) {
      
      warning(
        paste(
          sample_id,
          "does not have enough Fibroblast or Schwann cells for CellChat. Marking as NA."
        )
      )
      
      sample_summary <- tibble(
        orig.ident = sample_id,
        Sex = sample_sex,
        Direction = c(
          "Fibroblasts → Schwann cells",
          "Schwann cells → Fibroblasts",
          "Both directions"
        ),
        Fibroblast_Cell_Count = fibroblast_n,
        Schwann_Cell_Count = schwann_n,
        Interaction_Count = NA_real_,
        Total_Interaction_Strength = NA_real_,
        Mean_Interaction_Strength = NA_real_
      )
      
      per_sample_summary_list[[sample_id]] <- sample_summary
      next
    }
    
    rds_path <- file.path(
      per_sample_rds_dir,
      paste0(sample_id, "_fibroblast_schwann_cellchat.rds")
    )
    
    if (file.exists(rds_path)) {
      
      print(paste("Reading saved CellChat object:", rds_path))
      cellchat_sample <- readRDS(rds_path)
      cellchat_sample <- updateCellChat(cellchat_sample)
      
    } else {
      
      cellchat_sample <- run_cellchat_one_sample(
        seurat_obj = sample_obj,
        sample_id = sample_id,
        min_cells_filter = min_cells_required
      )
      
      saveRDS(
        cellchat_sample,
        rds_path
      )
    }
    
    df_sample <- subsetCommunication(cellchat_sample)
    
    fib_schwann_sample <- df_sample %>%
      filter(
        (source == "Fibroblasts" & target == "Schwann cells") |
          (source == "Schwann cells" & target == "Fibroblasts")
      ) %>%
      mutate(
        Direction = paste(source, "→", target)
      )
    
    if (nrow(fib_schwann_sample) == 0) {
      
      directional_summary <- tibble(
        Direction = c(
          "Fibroblasts → Schwann cells",
          "Schwann cells → Fibroblasts"
        ),
        Interaction_Count = 0,
        Total_Interaction_Strength = 0,
        Mean_Interaction_Strength = 0
      )
      
    } else {
      
      directional_summary <- fib_schwann_sample %>%
        group_by(Direction) %>%
        summarise(
          Interaction_Count = n(),
          Total_Interaction_Strength = sum(prob, na.rm = TRUE),
          Mean_Interaction_Strength = mean(prob, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        complete(
          Direction = c(
            "Fibroblasts → Schwann cells",
            "Schwann cells → Fibroblasts"
          ),
          fill = list(
            Interaction_Count = 0,
            Total_Interaction_Strength = 0,
            Mean_Interaction_Strength = 0
          )
        )
    }
    
    both_directions_summary <- directional_summary %>%
      summarise(
        Direction = "Both directions",
        Interaction_Count = sum(Interaction_Count, na.rm = TRUE),
        Total_Interaction_Strength = sum(Total_Interaction_Strength, na.rm = TRUE),
        Mean_Interaction_Strength = ifelse(
          sum(Interaction_Count, na.rm = TRUE) > 0,
          Total_Interaction_Strength / Interaction_Count,
          0
        )
      )
    
    sample_summary <- bind_rows(
      directional_summary,
      both_directions_summary
    ) %>%
      mutate(
        orig.ident = sample_id,
        Sex = sample_sex,
        Fibroblast_Cell_Count = fibroblast_n,
        Schwann_Cell_Count = schwann_n
      ) %>%
      dplyr::select(
        orig.ident,
        Sex,
        Direction,
        Fibroblast_Cell_Count,
        Schwann_Cell_Count,
        Interaction_Count,
        Total_Interaction_Strength,
        Mean_Interaction_Strength
      )
    
    per_sample_summary_list[[sample_id]] <- sample_summary
  }
  
  per_sample_summary <- bind_rows(per_sample_summary_list)
  
  write.csv(
    per_sample_summary,
    summary_csv_path,
    row.names = FALSE
  )
}

print("Per-sample Schwann-Fibroblast CellChat summary:")
print(per_sample_summary)

# ============================================================
# 8. T-TEST FUNCTION
# ============================================================

run_t_test_by_direction <- function(data, value_col) {
  
  results <- data %>%
    group_by(Direction) %>%
    group_modify(~{
      
      female_values <- .x[[value_col]][.x$Sex == "Female"]
      male_values <- .x[[value_col]][.x$Sex == "Male"]
      
      p_val <- tryCatch(
        {
          if (
            length(na.omit(female_values)) >= 2 &&
            length(na.omit(male_values)) >= 2
          ) {
            t.test(female_values, male_values)$p.value
          } else {
            NA_real_
          }
        },
        error = function(e) NA_real_
      )
      
      tibble(
        Female_Mean = mean(female_values, na.rm = TRUE),
        Male_Mean = mean(male_values, na.rm = TRUE),
        Female_SD = sd(female_values, na.rm = TRUE),
        Male_SD = sd(male_values, na.rm = TRUE),
        Female_N = length(na.omit(female_values)),
        Male_N = length(na.omit(male_values)),
        p_value = p_val
      )
    }) %>%
    ungroup() %>%
    mutate(
      p_adj_BH = p.adjust(p_value, method = "BH"),
      Significance = case_when(
        is.na(p_value) ~ "NA",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      Significance_BH = case_when(
        is.na(p_adj_BH) ~ "NA",
        p_adj_BH < 0.001 ~ "***",
        p_adj_BH < 0.01 ~ "**",
        p_adj_BH < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  
  return(results)
}

count_ttest_results <- run_t_test_by_direction(
  data = per_sample_summary,
  value_col = "Interaction_Count"
)

strength_ttest_results <- run_t_test_by_direction(
  data = per_sample_summary,
  value_col = "Total_Interaction_Strength"
)

write.csv(
  count_ttest_results,
  file.path(cellchat_outdir, "Ttest_Schwann_Fibroblast_Interaction_Count_By_Sex.csv"),
  row.names = FALSE
)

write.csv(
  strength_ttest_results,
  file.path(cellchat_outdir, "Ttest_Schwann_Fibroblast_Interaction_Strength_By_Sex.csv"),
  row.names = FALSE
)

print("T-test results for CellChat interaction count:")
print(count_ttest_results)

print("T-test results for CellChat interaction strength:")
print(strength_ttest_results)

# ============================================================
# 9. SUMMARY DATA FOR BAR PLOTS
# ============================================================

plot_summary <- per_sample_summary %>%
  group_by(Direction, Sex) %>%
  summarise(
    Mean_Interaction_Count = mean(Interaction_Count, na.rm = TRUE),
    SEM_Interaction_Count = sd(Interaction_Count, na.rm = TRUE) / sqrt(sum(!is.na(Interaction_Count))),
    Mean_Interaction_Strength = mean(Total_Interaction_Strength, na.rm = TRUE),
    SEM_Interaction_Strength = sd(Total_Interaction_Strength, na.rm = TRUE) / sqrt(sum(!is.na(Total_Interaction_Strength))),
    N = n(),
    .groups = "drop"
  )

write.csv(
  plot_summary,
  file.path(cellchat_outdir, "Schwann_Fibroblast_CellChat_Barplot_Summary.csv"),
  row.names = FALSE
)

# ============================================================
# 10. P-VALUE LABELS
# ============================================================

count_p_labels <- plot_summary %>%
  group_by(Direction) %>%
  summarise(
    y_position = max(
      Mean_Interaction_Count + SEM_Interaction_Count,
      na.rm = TRUE
    ) * 1.2,
    .groups = "drop"
  ) %>%
  left_join(
    count_ttest_results %>%
      dplyr::select(Direction, p_value, p_adj_BH, Significance),
    by = "Direction"
  ) %>%
  mutate(
    y_position = ifelse(
      !is.finite(y_position) | y_position == 0,
      1,
      y_position
    ),
    label = ifelse(
      is.na(p_value),
      "p = NA",
      paste0("p = ", signif(p_value, 2))
    )
  )

strength_p_labels <- plot_summary %>%
  group_by(Direction) %>%
  summarise(
    y_position = max(
      Mean_Interaction_Strength + SEM_Interaction_Strength,
      na.rm = TRUE
    ) * 1.2,
    .groups = "drop"
  ) %>%
  left_join(
    strength_ttest_results %>%
      dplyr::select(Direction, p_value, p_adj_BH, Significance),
    by = "Direction"
  ) %>%
  mutate(
    y_position = ifelse(
      !is.finite(y_position) | y_position == 0,
      0.01,
      y_position
    ),
    label = ifelse(
      is.na(p_value),
      "p = NA",
      paste0("p = ", signif(p_value, 2))
    )
  )

# ============================================================
# 11. BAR PLOT: INTERACTION COUNT
# ============================================================

count_plot <- ggplot(
  plot_summary,
  aes(
    x = Sex,
    y = Mean_Interaction_Count,
    fill = Sex
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_Interaction_Count - SEM_Interaction_Count,
      ymax = Mean_Interaction_Count + SEM_Interaction_Count
    ),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_text(
    data = count_p_labels,
    aes(
      x = 1.5,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  ) +
  facet_wrap(
    ~Direction,
    scales = "free_y",
    ncol = 3
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    title = "Schwann-Fibroblast CellChat Interaction Count by Sex",
    x = "Sex",
    y = "Mean number of inferred interactions per sample",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "none"
  )

print(count_plot)

ggsave(
  file.path(cellchat_outdir, "Schwann_Fibroblast_CellChat_Interaction_Count_Ttest_Barplot.png"),
  plot = count_plot,
  width = 12,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(cellchat_outdir, "Schwann_Fibroblast_CellChat_Interaction_Count_Ttest_Barplot.pdf"),
  plot = count_plot,
  width = 12,
  height = 5.5,
  dpi = 300
)

# ============================================================
# 12. BAR PLOT: INTERACTION STRENGTH
# ============================================================

strength_plot <- ggplot(
  plot_summary,
  aes(
    x = Sex,
    y = Mean_Interaction_Strength,
    fill = Sex
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_Interaction_Strength - SEM_Interaction_Strength,
      ymax = Mean_Interaction_Strength + SEM_Interaction_Strength
    ),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_text(
    data = strength_p_labels,
    aes(
      x = 1.5,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  ) +
  facet_wrap(
    ~Direction,
    scales = "free_y",
    ncol = 3
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    title = "Schwann-Fibroblast CellChat Interaction Strength by Sex",
    x = "Sex",
    y = "Mean total interaction strength per sample",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "none"
  )

print(strength_plot)

ggsave(
  file.path(cellchat_outdir, "Schwann_Fibroblast_CellChat_Interaction_Strength_Ttest_Barplot.png"),
  plot = strength_plot,
  width = 12,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(cellchat_outdir, "Schwann_Fibroblast_CellChat_Interaction_Strength_Ttest_Barplot.pdf"),
  plot = strength_plot,
  width = 12,
  height = 5.5,
  dpi = 300
)

# ============================================================
# 13. FINAL MESSAGE
# ============================================================

cat("\n======================================================\n")
cat("SCHWANN-FIBROBLAST CELLCHAT T-TEST ANALYSIS COMPLETE\n")
cat("Results saved to:\n")
cat(cellchat_outdir, "\n\n")
cat("Main outputs:\n")
cat("- Schwann_Fibroblast_CellChat_PerSample_Count_Strength.csv\n")
cat("- Ttest_Schwann_Fibroblast_Interaction_Count_By_Sex.csv\n")
cat("- Ttest_Schwann_Fibroblast_Interaction_Strength_By_Sex.csv\n")
cat("- Schwann_Fibroblast_CellChat_Interaction_Count_Ttest_Barplot.png/pdf\n")
cat("- Schwann_Fibroblast_CellChat_Interaction_Strength_Ttest_Barplot.png/pdf\n")
cat("======================================================\n")


### T Test for iCAF: myCAF ratio

# ============================================================
# iCAF:myCAF RATIO T-TEST ANALYSIS
# Female vs Male samples
# Statistical unit = patient/sample
# No individual dots/points on plot
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

options(stringsAsFactors = FALSE)

# ============================================================
# 1. LOAD FIBROBLAST SUBSET
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

fibroblast_subset <- readRDS(
  file.path(base_dir, "Fibroblast_Subset.RDS")
)

DefaultAssay(fibroblast_subset) <- "RNA"

fibroblast_subset <- tryCatch(
  JoinLayers(fibroblast_subset),
  error = function(e) fibroblast_subset
)

# ============================================================
# 2. ADD SEX METADATA
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

fibroblast_subset$Sex <- unname(
  sex_mapping[as.character(fibroblast_subset$orig.ident)]
)

fibroblast_subset$Sex <- factor(
  fibroblast_subset$Sex,
  levels = c("Female", "Male")
)

# ============================================================
# 3. LABEL CAF SUBTYPES
# ============================================================
# Based on your fibroblast subset clustering:
# cluster 1 = myCAF
# cluster 2 = iCAF

fibroblast_subset$CAF_type <- dplyr::case_when(
  as.character(fibroblast_subset$seurat_clusters) == "1" ~ "myCAF",
  as.character(fibroblast_subset$seurat_clusters) == "2" ~ "iCAF",
  TRUE ~ "Other"
)

fibroblast_subset$CAF_type <- factor(
  fibroblast_subset$CAF_type,
  levels = c("myCAF", "iCAF", "Other")
)

# ============================================================
# 4. CREATE OUTPUT FOLDER
# ============================================================

ratio_outdir <- file.path(
  base_dir,
  "iCAF_myCAF_Ratio_Ttest"
)

dir.create(
  ratio_outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# 5. COUNT iCAFs AND myCAFs PER SAMPLE
# ============================================================

icaf_mycaf_counts <- fibroblast_subset@meta.data %>%
  filter(
    !is.na(Sex),
    !is.na(orig.ident),
    CAF_type %in% c("iCAF", "myCAF")
  ) %>%
  dplyr::count(
    orig.ident,
    Sex,
    CAF_type,
    name = "Cell_Count"
  ) %>%
  complete(
    nesting(orig.ident, Sex),
    CAF_type = c("iCAF", "myCAF"),
    fill = list(Cell_Count = 0)
  )

icaf_mycaf_wide <- icaf_mycaf_counts %>%
  pivot_wider(
    names_from = CAF_type,
    values_from = Cell_Count,
    values_fill = 0
  ) %>%
  mutate(
    Total_iCAF_myCAF = iCAF + myCAF,
    
    # Raw iCAF:myCAF ratio
    iCAF_myCAF_Ratio = ifelse(
      myCAF > 0,
      iCAF / myCAF,
      NA_real_
    ),
    
    # Percent iCAF within iCAF + myCAF compartment
    iCAF_Percent_of_iCAF_myCAF = ifelse(
      Total_iCAF_myCAF > 0,
      100 * iCAF / Total_iCAF_myCAF,
      NA_real_
    ),
    
    # Log2 ratio with small pseudocount to avoid division by zero
    Log2_iCAF_myCAF_Ratio = log2(
      (iCAF + 0.5) / (myCAF + 0.5)
    )
  )

write.csv(
  icaf_mycaf_wide,
  file.path(ratio_outdir, "iCAF_myCAF_Counts_and_Ratios_Per_Sample.csv"),
  row.names = FALSE
)

print("iCAF and myCAF counts/ratios per sample:")
print(icaf_mycaf_wide)

# ============================================================
# 6. T-TEST FUNCTION
# ============================================================

run_ratio_ttest <- function(data, value_col) {
  
  female_values <- data[[value_col]][data$Sex == "Female"]
  male_values <- data[[value_col]][data$Sex == "Male"]
  
  p_val <- tryCatch(
    {
      if (
        length(na.omit(female_values)) >= 2 &&
        length(na.omit(male_values)) >= 2
      ) {
        t.test(female_values, male_values)$p.value
      } else {
        NA_real_
      }
    },
    error = function(e) NA_real_
  )
  
  result <- tibble(
    Metric = value_col,
    Female_Mean = mean(female_values, na.rm = TRUE),
    Male_Mean = mean(male_values, na.rm = TRUE),
    Female_SD = sd(female_values, na.rm = TRUE),
    Male_SD = sd(male_values, na.rm = TRUE),
    Female_N = length(na.omit(female_values)),
    Male_N = length(na.omit(male_values)),
    p_value = p_val,
    Significance = case_when(
      is.na(p_val) ~ "NA",
      p_val < 0.001 ~ "***",
      p_val < 0.01 ~ "**",
      p_val < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
  
  return(result)
}

ratio_ttest_results <- bind_rows(
  run_ratio_ttest(icaf_mycaf_wide, "iCAF_myCAF_Ratio"),
  run_ratio_ttest(icaf_mycaf_wide, "iCAF_Percent_of_iCAF_myCAF"),
  run_ratio_ttest(icaf_mycaf_wide, "Log2_iCAF_myCAF_Ratio")
) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    Significance_BH = case_when(
      is.na(p_adj_BH) ~ "NA",
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01 ~ "**",
      p_adj_BH < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

write.csv(
  ratio_ttest_results,
  file.path(ratio_outdir, "Ttest_iCAF_myCAF_Ratio_Female_vs_Male.csv"),
  row.names = FALSE
)

print("T-test results:")
print(ratio_ttest_results)

# ============================================================
# 7. BARPLOT SUMMARY FOR RAW iCAF:myCAF RATIO
# ============================================================

ratio_plot_summary <- icaf_mycaf_wide %>%
  group_by(Sex) %>%
  summarise(
    Mean_Ratio = mean(iCAF_myCAF_Ratio, na.rm = TRUE),
    SEM_Ratio = sd(iCAF_myCAF_Ratio, na.rm = TRUE) /
      sqrt(sum(!is.na(iCAF_myCAF_Ratio))),
    N = sum(!is.na(iCAF_myCAF_Ratio)),
    .groups = "drop"
  )

ratio_p_value <- ratio_ttest_results %>%
  filter(Metric == "iCAF_myCAF_Ratio") %>%
  pull(p_value)

ratio_p_label <- ifelse(
  is.na(ratio_p_value),
  "p = NA",
  paste0("p = ", signif(ratio_p_value, 2))
)

ratio_y_position <- max(
  ratio_plot_summary$Mean_Ratio + ratio_plot_summary$SEM_Ratio,
  na.rm = TRUE
) * 1.25

if (!is.finite(ratio_y_position) || ratio_y_position == 0) {
  ratio_y_position <- 1
}

# ============================================================
# 8. RAW iCAF:myCAF RATIO BARPLOT
# ============================================================

ratio_plot <- ggplot(
  ratio_plot_summary,
  aes(
    x = Sex,
    y = Mean_Ratio,
    fill = Sex
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_Ratio - SEM_Ratio,
      ymax = Mean_Ratio + SEM_Ratio
    ),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      x = 1.5,
      y = ratio_y_position,
      label = ratio_p_label
    ),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    title = "iCAF:myCAF Ratio by Sex",
    x = "Sex",
    y = "Mean iCAF:myCAF ratio per sample",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "none"
  )

print(ratio_plot)

ggsave(
  file.path(ratio_outdir, "iCAF_myCAF_Ratio_Female_vs_Male_Barplot.png"),
  plot = ratio_plot,
  width = 6,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(ratio_outdir, "iCAF_myCAF_Ratio_Female_vs_Male_Barplot.pdf"),
  plot = ratio_plot,
  width = 6,
  height = 5,
  dpi = 300
)

# ============================================================
# 9. BARPLOT SUMMARY FOR LOG2 iCAF:myCAF RATIO
# ============================================================

log2_plot_summary <- icaf_mycaf_wide %>%
  group_by(Sex) %>%
  summarise(
    Mean_Log2_Ratio = mean(Log2_iCAF_myCAF_Ratio, na.rm = TRUE),
    SEM_Log2_Ratio = sd(Log2_iCAF_myCAF_Ratio, na.rm = TRUE) /
      sqrt(sum(!is.na(Log2_iCAF_myCAF_Ratio))),
    N = sum(!is.na(Log2_iCAF_myCAF_Ratio)),
    .groups = "drop"
  )

log2_p_value <- ratio_ttest_results %>%
  filter(Metric == "Log2_iCAF_myCAF_Ratio") %>%
  pull(p_value)

log2_p_label <- ifelse(
  is.na(log2_p_value),
  "p = NA",
  paste0("p = ", signif(log2_p_value, 2))
)

log2_y_position <- max(
  log2_plot_summary$Mean_Log2_Ratio + log2_plot_summary$SEM_Log2_Ratio,
  na.rm = TRUE
) * 1.25

if (!is.finite(log2_y_position)) {
  log2_y_position <- 1
}

# ============================================================
# 10. LOG2 iCAF:myCAF RATIO BARPLOT
# ============================================================

log2_ratio_plot <- ggplot(
  log2_plot_summary,
  aes(
    x = Sex,
    y = Mean_Log2_Ratio,
    fill = Sex
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_Log2_Ratio - SEM_Log2_Ratio,
      ymax = Mean_Log2_Ratio + SEM_Log2_Ratio
    ),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      x = 1.5,
      y = log2_y_position,
      label = log2_p_label
    ),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.15, 0.25))
  ) +
  labs(
    title = "Log2 iCAF:myCAF Ratio by Sex",
    x = "Sex",
    y = expression("Mean log"[2] * "((iCAF + 0.5)/(myCAF + 0.5)) per sample"),
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "none"
  )

print(log2_ratio_plot)

ggsave(
  file.path(ratio_outdir, "Log2_iCAF_myCAF_Ratio_Female_vs_Male_Barplot.png"),
  plot = log2_ratio_plot,
  width = 6.5,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(ratio_outdir, "Log2_iCAF_myCAF_Ratio_Female_vs_Male_Barplot.pdf"),
  plot = log2_ratio_plot,
  width = 6.5,
  height = 5,
  dpi = 300
)

# ============================================================
# 11. FINAL MESSAGE
# ============================================================

cat("\n======================================================\n")
cat("iCAF:myCAF RATIO T-TEST ANALYSIS COMPLETE\n")
cat("Results saved to:\n")
cat(ratio_outdir, "\n\n")
cat("Main outputs:\n")
cat("- iCAF_myCAF_Counts_and_Ratios_Per_Sample.csv\n")
cat("- Ttest_iCAF_myCAF_Ratio_Female_vs_Male.csv\n")
cat("- iCAF_myCAF_Ratio_Female_vs_Male_Barplot.png/pdf\n")
cat("- Log2_iCAF_myCAF_Ratio_Female_vs_Male_Barplot.png/pdf\n")
cat("======================================================\n")



### T Test on CAF Subsetting

# ============================================================
# FIBROBLAST SUBTYPE COMPOSITION T-TEST ANALYSIS
# iCAF, myCAF, apCAF, PSCs, Fibroblasts
# Female vs Male samples
# Statistical unit = patient/sample
# Bar graphs for cell counts and proportions
# No individual dots/points
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

options(stringsAsFactors = FALSE)

# ============================================================
# 1. LOAD FIBROBLAST SUBSET
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

fibroblast_subset <- readRDS(
  file.path(base_dir, "Fibroblast_Subset.RDS")
)

DefaultAssay(fibroblast_subset) <- "RNA"

fibroblast_subset <- tryCatch(
  JoinLayers(fibroblast_subset),
  error = function(e) fibroblast_subset
)

# ============================================================
# 2. ADD SEX METADATA
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

fibroblast_subset$Sex <- unname(
  sex_mapping[as.character(fibroblast_subset$orig.ident)]
)

fibroblast_subset$Sex <- factor(
  fibroblast_subset$Sex,
  levels = c("Female", "Male")
)

# ============================================================
# 3. LABEL FIBROBLAST SUBTYPES
# ============================================================
# Based on your fibroblast subset clustering:
# 0 = Fibroblasts
# 1 = myCAF
# 2 = iCAF
# 3 = Pancreatic stellate / PSC
# 4 = apCAF

fibroblast_subset$CAF_type <- dplyr::case_when(
  as.character(fibroblast_subset$seurat_clusters) == "0" ~ "Fibroblasts",
  as.character(fibroblast_subset$seurat_clusters) == "1" ~ "myCAFs",
  as.character(fibroblast_subset$seurat_clusters) == "2" ~ "iCAFs",
  as.character(fibroblast_subset$seurat_clusters) == "3" ~ "PSCs",
  as.character(fibroblast_subset$seurat_clusters) == "4" ~ "apCAFs",
  TRUE ~ "Other"
)

caf_order <- c(
  "Fibroblasts",
  "myCAFs",
  "iCAFs",
  "PSCs",
  "apCAFs"
)

fibroblast_subset$CAF_type <- factor(
  fibroblast_subset$CAF_type,
  levels = c(caf_order, "Other")
)

# ============================================================
# 4. CREATE OUTPUT FOLDER
# ============================================================

caf_outdir <- file.path(
  base_dir,
  "Fibroblast_Subtype_Composition_Ttests"
)

dir.create(
  caf_outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# 5. COUNT CAF SUBTYPES PER SAMPLE
# ============================================================

caf_counts_per_sample <- fibroblast_subset@meta.data %>%
  filter(
    !is.na(Sex),
    !is.na(orig.ident),
    CAF_type %in% caf_order
  ) %>%
  dplyr::count(
    orig.ident,
    Sex,
    CAF_type,
    name = "Cell_Count"
  ) %>%
  complete(
    nesting(orig.ident, Sex),
    CAF_type = factor(caf_order, levels = caf_order),
    fill = list(Cell_Count = 0)
  ) %>%
  mutate(
    CAF_type = factor(CAF_type, levels = caf_order)
  )

sample_totals <- caf_counts_per_sample %>%
  group_by(orig.ident, Sex) %>%
  summarise(
    Total_Fibroblast_Subset_Cells = sum(Cell_Count),
    .groups = "drop"
  )

caf_composition_per_sample <- caf_counts_per_sample %>%
  left_join(
    sample_totals,
    by = c("orig.ident", "Sex")
  ) %>%
  mutate(
    Percent_of_Fibroblast_Subset = ifelse(
      Total_Fibroblast_Subset_Cells > 0,
      100 * Cell_Count / Total_Fibroblast_Subset_Cells,
      NA_real_
    )
  )

write.csv(
  caf_composition_per_sample,
  file.path(caf_outdir, "CAF_Subtype_Counts_and_Proportions_Per_Sample.csv"),
  row.names = FALSE
)

print("CAF subtype counts and proportions per sample:")
print(caf_composition_per_sample)

# ============================================================
# 6. T-TEST FUNCTION
# ============================================================

run_t_test_by_caf_type <- function(data, value_col) {
  
  results <- data %>%
    group_by(CAF_type) %>%
    group_modify(~{
      
      female_values <- .x[[value_col]][.x$Sex == "Female"]
      male_values <- .x[[value_col]][.x$Sex == "Male"]
      
      p_val <- tryCatch(
        {
          if (
            length(na.omit(female_values)) >= 2 &&
            length(na.omit(male_values)) >= 2
          ) {
            t.test(female_values, male_values)$p.value
          } else {
            NA_real_
          }
        },
        error = function(e) NA_real_
      )
      
      tibble(
        Female_Mean = mean(female_values, na.rm = TRUE),
        Male_Mean = mean(male_values, na.rm = TRUE),
        Female_SD = sd(female_values, na.rm = TRUE),
        Male_SD = sd(male_values, na.rm = TRUE),
        Female_N = length(na.omit(female_values)),
        Male_N = length(na.omit(male_values)),
        p_value = p_val
      )
    }) %>%
    ungroup() %>%
    mutate(
      p_adj_BH = p.adjust(p_value, method = "BH"),
      Significance = case_when(
        is.na(p_value) ~ "NA",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      Significance_BH = case_when(
        is.na(p_adj_BH) ~ "NA",
        p_adj_BH < 0.001 ~ "***",
        p_adj_BH < 0.01 ~ "**",
        p_adj_BH < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  
  return(results)
}

count_ttest_results <- run_t_test_by_caf_type(
  data = caf_composition_per_sample,
  value_col = "Cell_Count"
)

proportion_ttest_results <- run_t_test_by_caf_type(
  data = caf_composition_per_sample,
  value_col = "Percent_of_Fibroblast_Subset"
)

write.csv(
  count_ttest_results,
  file.path(caf_outdir, "Ttest_CAF_Subtype_Cell_Counts_By_Sex.csv"),
  row.names = FALSE
)

write.csv(
  proportion_ttest_results,
  file.path(caf_outdir, "Ttest_CAF_Subtype_Proportions_By_Sex.csv"),
  row.names = FALSE
)

print("T-test results for CAF subtype cell counts:")
print(count_ttest_results)

print("T-test results for CAF subtype proportions:")
print(proportion_ttest_results)

# ============================================================
# 7. SUMMARY DATA FOR BAR PLOTS
# ============================================================

plot_summary <- caf_composition_per_sample %>%
  group_by(CAF_type, Sex) %>%
  summarise(
    Mean_Cell_Count = mean(Cell_Count, na.rm = TRUE),
    SEM_Cell_Count = sd(Cell_Count, na.rm = TRUE) /
      sqrt(sum(!is.na(Cell_Count))),
    Mean_Percent = mean(Percent_of_Fibroblast_Subset, na.rm = TRUE),
    SEM_Percent = sd(Percent_of_Fibroblast_Subset, na.rm = TRUE) /
      sqrt(sum(!is.na(Percent_of_Fibroblast_Subset))),
    N = n(),
    .groups = "drop"
  )

write.csv(
  plot_summary,
  file.path(caf_outdir, "CAF_Subtype_Barplot_Summary.csv"),
  row.names = FALSE
)

# ============================================================
# 8. CREATE P-VALUE LABELS
# ============================================================

count_p_labels <- plot_summary %>%
  group_by(CAF_type) %>%
  summarise(
    y_position = max(Mean_Cell_Count + SEM_Cell_Count, na.rm = TRUE) * 1.25,
    .groups = "drop"
  ) %>%
  left_join(
    count_ttest_results %>%
      dplyr::select(CAF_type, p_value, p_adj_BH, Significance),
    by = "CAF_type"
  ) %>%
  mutate(
    y_position = ifelse(
      !is.finite(y_position) | y_position == 0,
      1,
      y_position
    ),
    label = ifelse(
      is.na(p_value),
      "p = NA",
      paste0("p = ", signif(p_value, 2))
    )
  )

proportion_p_labels <- plot_summary %>%
  group_by(CAF_type) %>%
  summarise(
    y_position = max(Mean_Percent + SEM_Percent, na.rm = TRUE) * 1.25,
    .groups = "drop"
  ) %>%
  left_join(
    proportion_ttest_results %>%
      dplyr::select(CAF_type, p_value, p_adj_BH, Significance),
    by = "CAF_type"
  ) %>%
  mutate(
    y_position = ifelse(
      !is.finite(y_position) | y_position == 0,
      0.05,
      y_position
    ),
    label = ifelse(
      is.na(p_value),
      "p = NA",
      paste0("p = ", signif(p_value, 2))
    )
  )

# ============================================================
# 9. BARPLOT: CAF SUBTYPE CELL COUNTS
# ============================================================

count_plot <- ggplot(
  plot_summary,
  aes(
    x = Sex,
    y = Mean_Cell_Count,
    fill = Sex
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_Cell_Count - SEM_Cell_Count,
      ymax = Mean_Cell_Count + SEM_Cell_Count
    ),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_text(
    data = count_p_labels,
    aes(
      x = 1.5,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    size = 3.8,
    fontface = "bold"
  ) +
  facet_wrap(
    ~CAF_type,
    scales = "free_y",
    ncol = 3
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.28))
  ) +
  labs(
    title = "Fibroblast Subtype Cell Counts by Sex",
    x = "Sex",
    y = "Mean number of cells per sample",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    axis.text = element_text(
      size = 11,
      color = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    legend.position = "none"
  )

print(count_plot)

ggsave(
  file.path(caf_outdir, "CAF_Subtype_Cell_Counts_By_Sex_Ttest_Barplot.png"),
  plot = count_plot,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(caf_outdir, "CAF_Subtype_Cell_Counts_By_Sex_Ttest_Barplot.pdf"),
  plot = count_plot,
  width = 10,
  height = 7,
  dpi = 300
)

# ============================================================
# 10. BARPLOT: CAF SUBTYPE PROPORTIONS
# ============================================================

proportion_plot <- ggplot(
  plot_summary,
  aes(
    x = Sex,
    y = Mean_Percent,
    fill = Sex
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_Percent - SEM_Percent,
      ymax = Mean_Percent + SEM_Percent
    ),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_text(
    data = proportion_p_labels,
    aes(
      x = 1.5,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    size = 3.8,
    fontface = "bold"
  ) +
  facet_wrap(
    ~CAF_type,
    scales = "free_y",
    ncol = 3
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.28))
  ) +
  labs(
    title = "Fibroblast Subtype Proportions by Sex",
    x = "Sex",
    y = "Mean percent of fibroblast subset per sample",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    axis.text = element_text(
      size = 11,
      color = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    legend.position = "none"
  )

print(proportion_plot)

ggsave(
  file.path(caf_outdir, "CAF_Subtype_Proportions_By_Sex_Ttest_Barplot.png"),
  plot = proportion_plot,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(caf_outdir, "CAF_Subtype_Proportions_By_Sex_Ttest_Barplot.pdf"),
  plot = proportion_plot,
  width = 10,
  height = 7,
  dpi = 300
)

# ============================================================
# 11. OPTIONAL: SIGNIFICANT SUBTYPES ONLY - PROPORTIONS
# ============================================================

significant_proportion_subtypes <- proportion_ttest_results %>%
  filter(!is.na(p_value), p_value < 0.05) %>%
  pull(CAF_type)

if (length(significant_proportion_subtypes) > 0) {
  
  proportion_plot_sig <- ggplot(
    plot_summary %>%
      filter(CAF_type %in% significant_proportion_subtypes),
    aes(
      x = Sex,
      y = Mean_Percent,
      fill = Sex
    )
  ) +
    geom_col(
      color = "black",
      linewidth = 0.3,
      width = 0.65
    ) +
    geom_errorbar(
      aes(
        ymin = Mean_Percent - SEM_Percent,
        ymax = Mean_Percent + SEM_Percent
      ),
      width = 0.2,
      linewidth = 0.4
    ) +
    geom_text(
      data = proportion_p_labels %>%
        filter(CAF_type %in% significant_proportion_subtypes),
      aes(
        x = 1.5,
        y = y_position,
        label = label
      ),
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold"
    ) +
    facet_wrap(
      ~CAF_type,
      scales = "free_y",
      ncol = 3
    ) +
    scale_fill_manual(
      values = c(
        "Female" = "#F8766D",
        "Male" = "#00BFC4"
      )
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.28))
    ) +
    labs(
      title = "Significant Fibroblast Subtype Proportion Differences by Sex",
      x = "Sex",
      y = "Mean percent of fibroblast subset per sample",
      fill = "Sex"
    ) +
    theme_classic(base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 11, color = "black"),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "none"
    )
  
  print(proportion_plot_sig)
  
  ggsave(
    file.path(caf_outdir, "Significant_CAF_Subtype_Proportion_Differences_By_Sex.png"),
    plot = proportion_plot_sig,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  ggsave(
    file.path(caf_outdir, "Significant_CAF_Subtype_Proportion_Differences_By_Sex.pdf"),
    plot = proportion_plot_sig,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# ============================================================
# 12. OPTIONAL: SIGNIFICANT SUBTYPES ONLY - COUNTS
# ============================================================

significant_count_subtypes <- count_ttest_results %>%
  filter(!is.na(p_value), p_value < 0.05) %>%
  pull(CAF_type)

if (length(significant_count_subtypes) > 0) {
  
  count_plot_sig <- ggplot(
    plot_summary %>%
      filter(CAF_type %in% significant_count_subtypes),
    aes(
      x = Sex,
      y = Mean_Cell_Count,
      fill = Sex
    )
  ) +
    geom_col(
      color = "black",
      linewidth = 0.3,
      width = 0.65
    ) +
    geom_errorbar(
      aes(
        ymin = Mean_Cell_Count - SEM_Cell_Count,
        ymax = Mean_Cell_Count + SEM_Cell_Count
      ),
      width = 0.2,
      linewidth = 0.4
    ) +
    geom_text(
      data = count_p_labels %>%
        filter(CAF_type %in% significant_count_subtypes),
      aes(
        x = 1.5,
        y = y_position,
        label = label
      ),
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold"
    ) +
    facet_wrap(
      ~CAF_type,
      scales = "free_y",
      ncol = 3
    ) +
    scale_fill_manual(
      values = c(
        "Female" = "#F8766D",
        "Male" = "#00BFC4"
      )
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.28))
    ) +
    labs(
      title = "Significant Fibroblast Subtype Cell Count Differences by Sex",
      x = "Sex",
      y = "Mean number of cells per sample",
      fill = "Sex"
    ) +
    theme_classic(base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 11, color = "black"),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "none"
    )
  
  print(count_plot_sig)
  
  ggsave(
    file.path(caf_outdir, "Significant_CAF_Subtype_Cell_Count_Differences_By_Sex.png"),
    plot = count_plot_sig,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  ggsave(
    file.path(caf_outdir, "Significant_CAF_Subtype_Cell_Count_Differences_By_Sex.pdf"),
    plot = count_plot_sig,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# ============================================================
# 13. FINAL MESSAGE
# ============================================================

cat("\n======================================================\n")
cat("FIBROBLAST SUBTYPE COMPOSITION T-TEST ANALYSIS COMPLETE\n")
cat("Results saved to:\n")
cat(caf_outdir, "\n\n")
cat("Main outputs:\n")
cat("- CAF_Subtype_Counts_and_Proportions_Per_Sample.csv\n")
cat("- Ttest_CAF_Subtype_Cell_Counts_By_Sex.csv\n")
cat("- Ttest_CAF_Subtype_Proportions_By_Sex.csv\n")
cat("- CAF_Subtype_Cell_Counts_By_Sex_Ttest_Barplot.png/pdf\n")
cat("- CAF_Subtype_Proportions_By_Sex_Ttest_Barplot.png/pdf\n")
cat("- CAF_Subtype_Barplot_Summary.csv\n")
cat("======================================================\n")


### T Test of CAF Subtypes


# ============================================================
# FIBROBLAST SUBTYPE COMPOSITION T-TEST ANALYSIS
# General Fibroblasts, myCAFs, iCAFs, PSCs, apCAFs
# Female vs Male samples
# Statistical unit = patient/sample
# Bar graphs for counts and proportions
# No individual dots/points
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

options(stringsAsFactors = FALSE)

# ============================================================
# 1. LOAD FIBROBLAST SUBSET
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

fibroblast_path <- file.path(
  base_dir,
  "Fibroblast_Subset.RDS"
)

fibroblast_subset <- readRDS(fibroblast_path)

DefaultAssay(fibroblast_subset) <- "RNA"

fibroblast_subset <- tryCatch(
  JoinLayers(fibroblast_subset),
  error = function(e) fibroblast_subset
)

print("Number of fibroblast subset cells:")
print(ncol(fibroblast_subset))

# ============================================================
# 2. ADD SEX METADATA
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

fibroblast_subset$Sex <- unname(
  sex_mapping[as.character(fibroblast_subset$orig.ident)]
)

fibroblast_subset$Sex <- factor(
  fibroblast_subset$Sex,
  levels = c("Female", "Male")
)

# ============================================================
# 3. LABEL FIBROBLAST SUBTYPES
# ============================================================
# Your fibroblast subset clustering:
# 0 = General Fibroblasts
# 1 = myCAFs
# 2 = iCAFs
# 3 = Pancreatic stellate cells / PSCs
# 4 = apCAFs

fibroblast_subset$Fibroblast_Subtype <- dplyr::case_when(
  as.character(fibroblast_subset$seurat_clusters) == "0" ~ "General Fibroblasts",
  as.character(fibroblast_subset$seurat_clusters) == "1" ~ "myCAFs",
  as.character(fibroblast_subset$seurat_clusters) == "2" ~ "iCAFs",
  as.character(fibroblast_subset$seurat_clusters) == "3" ~ "PSCs",
  as.character(fibroblast_subset$seurat_clusters) == "4" ~ "apCAFs",
  TRUE ~ "Other"
)

fibroblast_order <- c(
  "General Fibroblasts",
  "myCAFs",
  "iCAFs",
  "PSCs",
  "apCAFs"
)

fibroblast_subset$Fibroblast_Subtype <- factor(
  fibroblast_subset$Fibroblast_Subtype,
  levels = c(fibroblast_order, "Other")
)

print("Fibroblast subtype counts:")
print(table(fibroblast_subset$Fibroblast_Subtype))

print("Fibroblast subtype counts by sex:")
print(table(fibroblast_subset$Fibroblast_Subtype, fibroblast_subset$Sex))

# ============================================================
# 4. CREATE OUTPUT FOLDER
# ============================================================

fibroblast_outdir <- file.path(
  base_dir,
  "Fibroblast_Subtype_Composition_Ttests"
)

dir.create(
  fibroblast_outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# 5. COUNT FIBROBLAST SUBTYPES PER SAMPLE
# ============================================================

fibroblast_counts_per_sample <- fibroblast_subset@meta.data %>%
  filter(
    !is.na(Sex),
    !is.na(orig.ident),
    Fibroblast_Subtype %in% fibroblast_order
  ) %>%
  dplyr::count(
    orig.ident,
    Sex,
    Fibroblast_Subtype,
    name = "Cell_Count"
  ) %>%
  complete(
    nesting(orig.ident, Sex),
    Fibroblast_Subtype = factor(
      fibroblast_order,
      levels = fibroblast_order
    ),
    fill = list(Cell_Count = 0)
  ) %>%
  mutate(
    Fibroblast_Subtype = factor(
      Fibroblast_Subtype,
      levels = fibroblast_order
    )
  )

sample_totals <- fibroblast_counts_per_sample %>%
  group_by(orig.ident, Sex) %>%
  summarise(
    Total_Fibroblast_Cells = sum(Cell_Count),
    .groups = "drop"
  )

fibroblast_composition_per_sample <- fibroblast_counts_per_sample %>%
  left_join(
    sample_totals,
    by = c("orig.ident", "Sex")
  ) %>%
  mutate(
    Percent_of_Fibroblast_Subset = ifelse(
      Total_Fibroblast_Cells > 0,
      100 * Cell_Count / Total_Fibroblast_Cells,
      NA_real_
    )
  )

write.csv(
  fibroblast_composition_per_sample,
  file.path(fibroblast_outdir, "Fibroblast_Subtype_Counts_and_Proportions_Per_Sample.csv"),
  row.names = FALSE
)

print("Fibroblast subtype counts and proportions per sample:")
print(fibroblast_composition_per_sample)

# ============================================================
# 6. T-TEST FUNCTION
# ============================================================

run_t_test_by_fibroblast_subtype <- function(data, value_col) {
  
  results <- data %>%
    group_by(Fibroblast_Subtype) %>%
    group_modify(~{
      
      female_values <- .x[[value_col]][.x$Sex == "Female"]
      male_values <- .x[[value_col]][.x$Sex == "Male"]
      
      p_val <- tryCatch(
        {
          if (
            length(na.omit(female_values)) >= 2 &&
            length(na.omit(male_values)) >= 2
          ) {
            t.test(female_values, male_values)$p.value
          } else {
            NA_real_
          }
        },
        error = function(e) NA_real_
      )
      
      tibble(
        Female_Mean = mean(female_values, na.rm = TRUE),
        Male_Mean = mean(male_values, na.rm = TRUE),
        Female_SD = sd(female_values, na.rm = TRUE),
        Male_SD = sd(male_values, na.rm = TRUE),
        Female_N = length(na.omit(female_values)),
        Male_N = length(na.omit(male_values)),
        p_value = p_val
      )
    }) %>%
    ungroup() %>%
    mutate(
      p_adj_BH = p.adjust(p_value, method = "BH"),
      Significance = case_when(
        is.na(p_value) ~ "NA",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      Significance_BH = case_when(
        is.na(p_adj_BH) ~ "NA",
        p_adj_BH < 0.001 ~ "***",
        p_adj_BH < 0.01 ~ "**",
        p_adj_BH < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  
  return(results)
}

count_ttest_results <- run_t_test_by_fibroblast_subtype(
  data = fibroblast_composition_per_sample,
  value_col = "Cell_Count"
)

proportion_ttest_results <- run_t_test_by_fibroblast_subtype(
  data = fibroblast_composition_per_sample,
  value_col = "Percent_of_Fibroblast_Subset"
)

write.csv(
  count_ttest_results,
  file.path(fibroblast_outdir, "Ttest_Fibroblast_Subtype_Cell_Counts_By_Sex.csv"),
  row.names = FALSE
)

write.csv(
  proportion_ttest_results,
  file.path(fibroblast_outdir, "Ttest_Fibroblast_Subtype_Proportions_By_Sex.csv"),
  row.names = FALSE
)

print("T-test results for fibroblast subtype cell counts:")
print(count_ttest_results)

print("T-test results for fibroblast subtype proportions:")
print(proportion_ttest_results)

# ============================================================
# 7. SUMMARY DATA FOR BAR PLOTS
# ============================================================

plot_summary <- fibroblast_composition_per_sample %>%
  group_by(Fibroblast_Subtype, Sex) %>%
  summarise(
    Mean_Cell_Count = mean(Cell_Count, na.rm = TRUE),
    SEM_Cell_Count = sd(Cell_Count, na.rm = TRUE) /
      sqrt(sum(!is.na(Cell_Count))),
    Mean_Percent = mean(Percent_of_Fibroblast_Subset, na.rm = TRUE),
    SEM_Percent = sd(Percent_of_Fibroblast_Subset, na.rm = TRUE) /
      sqrt(sum(!is.na(Percent_of_Fibroblast_Subset))),
    N = n(),
    .groups = "drop"
  )

write.csv(
  plot_summary,
  file.path(fibroblast_outdir, "Fibroblast_Subtype_Barplot_Summary.csv"),
  row.names = FALSE
)

# ============================================================
# 8. CREATE P-VALUE LABELS
# ============================================================

count_p_labels <- plot_summary %>%
  group_by(Fibroblast_Subtype) %>%
  summarise(
    y_position = max(Mean_Cell_Count + SEM_Cell_Count, na.rm = TRUE) * 1.25,
    .groups = "drop"
  ) %>%
  left_join(
    count_ttest_results %>%
      dplyr::select(Fibroblast_Subtype, p_value, p_adj_BH, Significance),
    by = "Fibroblast_Subtype"
  ) %>%
  mutate(
    y_position = ifelse(
      !is.finite(y_position) | y_position == 0,
      1,
      y_position
    ),
    label = ifelse(
      is.na(p_value),
      "p = NA",
      paste0("p = ", signif(p_value, 2))
    )
  )

proportion_p_labels <- plot_summary %>%
  group_by(Fibroblast_Subtype) %>%
  summarise(
    y_position = max(Mean_Percent + SEM_Percent, na.rm = TRUE) * 1.25,
    .groups = "drop"
  ) %>%
  left_join(
    proportion_ttest_results %>%
      dplyr::select(Fibroblast_Subtype, p_value, p_adj_BH, Significance),
    by = "Fibroblast_Subtype"
  ) %>%
  mutate(
    y_position = ifelse(
      !is.finite(y_position) | y_position == 0,
      0.05,
      y_position
    ),
    label = ifelse(
      is.na(p_value),
      "p = NA",
      paste0("p = ", signif(p_value, 2))
    )
  )

# ============================================================
# 9. BARPLOT: FIBROBLAST SUBTYPE CELL COUNTS
# ============================================================

count_plot <- ggplot(
  plot_summary,
  aes(
    x = Sex,
    y = Mean_Cell_Count,
    fill = Sex
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_Cell_Count - SEM_Cell_Count,
      ymax = Mean_Cell_Count + SEM_Cell_Count
    ),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_text(
    data = count_p_labels,
    aes(
      x = 1.5,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    size = 3.8,
    fontface = "bold"
  ) +
  facet_wrap(
    ~Fibroblast_Subtype,
    scales = "free_y",
    ncol = 3
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.28))
  ) +
  labs(
    title = "Fibroblast Subtype Cell Counts by Sex",
    x = "Sex",
    y = "Mean number of cells per sample",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    axis.text = element_text(
      size = 11,
      color = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    legend.position = "none"
  )

print(count_plot)

ggsave(
  file.path(fibroblast_outdir, "Fibroblast_Subtype_Cell_Counts_By_Sex_Ttest_Barplot.png"),
  plot = count_plot,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(fibroblast_outdir, "Fibroblast_Subtype_Cell_Counts_By_Sex_Ttest_Barplot.pdf"),
  plot = count_plot,
  width = 10,
  height = 7,
  dpi = 300
)

# ============================================================
# 10. BARPLOT: FIBROBLAST SUBTYPE PROPORTIONS
# ============================================================

proportion_plot <- ggplot(
  plot_summary,
  aes(
    x = Sex,
    y = Mean_Percent,
    fill = Sex
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_Percent - SEM_Percent,
      ymax = Mean_Percent + SEM_Percent
    ),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_text(
    data = proportion_p_labels,
    aes(
      x = 1.5,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    size = 3.8,
    fontface = "bold"
  ) +
  facet_wrap(
    ~Fibroblast_Subtype,
    scales = "free_y",
    ncol = 3
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.28))
  ) +
  labs(
    title = "Fibroblast Subtype Proportions by Sex",
    x = "Sex",
    y = "Mean percent of fibroblast subset per sample",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    axis.text = element_text(
      size = 11,
      color = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    legend.position = "none"
  )

print(proportion_plot)

ggsave(
  file.path(fibroblast_outdir, "Fibroblast_Subtype_Proportions_By_Sex_Ttest_Barplot.png"),
  plot = proportion_plot,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(fibroblast_outdir, "Fibroblast_Subtype_Proportions_By_Sex_Ttest_Barplot.pdf"),
  plot = proportion_plot,
  width = 10,
  height = 7,
  dpi = 300
)

# ============================================================
# 11. OPTIONAL: SIGNIFICANT SUBTYPES ONLY - PROPORTIONS
# ============================================================

significant_proportion_subtypes <- proportion_ttest_results %>%
  filter(!is.na(p_value), p_value < 0.05) %>%
  pull(Fibroblast_Subtype)

if (length(significant_proportion_subtypes) > 0) {
  
  proportion_plot_sig <- ggplot(
    plot_summary %>%
      filter(Fibroblast_Subtype %in% significant_proportion_subtypes),
    aes(
      x = Sex,
      y = Mean_Percent,
      fill = Sex
    )
  ) +
    geom_col(
      color = "black",
      linewidth = 0.3,
      width = 0.65
    ) +
    geom_errorbar(
      aes(
        ymin = Mean_Percent - SEM_Percent,
        ymax = Mean_Percent + SEM_Percent
      ),
      width = 0.2,
      linewidth = 0.4
    ) +
    geom_text(
      data = proportion_p_labels %>%
        filter(Fibroblast_Subtype %in% significant_proportion_subtypes),
      aes(
        x = 1.5,
        y = y_position,
        label = label
      ),
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold"
    ) +
    facet_wrap(
      ~Fibroblast_Subtype,
      scales = "free_y",
      ncol = 3
    ) +
    scale_fill_manual(
      values = c(
        "Female" = "#F8766D",
        "Male" = "#00BFC4"
      )
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.28))
    ) +
    labs(
      title = "Significant Fibroblast Subtype Proportion Differences by Sex",
      x = "Sex",
      y = "Mean percent of fibroblast subset per sample",
      fill = "Sex"
    ) +
    theme_classic(base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 11, color = "black"),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "none"
    )
  
  print(proportion_plot_sig)
  
  ggsave(
    file.path(fibroblast_outdir, "Significant_Fibroblast_Subtype_Proportion_Differences_By_Sex.png"),
    plot = proportion_plot_sig,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  ggsave(
    file.path(fibroblast_outdir, "Significant_Fibroblast_Subtype_Proportion_Differences_By_Sex.pdf"),
    plot = proportion_plot_sig,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# ============================================================
# 12. OPTIONAL: SIGNIFICANT SUBTYPES ONLY - COUNTS
# ============================================================

significant_count_subtypes <- count_ttest_results %>%
  filter(!is.na(p_value), p_value < 0.05) %>%
  pull(Fibroblast_Subtype)

if (length(significant_count_subtypes) > 0) {
  
  count_plot_sig <- ggplot(
    plot_summary %>%
      filter(Fibroblast_Subtype %in% significant_count_subtypes),
    aes(
      x = Sex,
      y = Mean_Cell_Count,
      fill = Sex
    )
  ) +
    geom_col(
      color = "black",
      linewidth = 0.3,
      width = 0.65
    ) +
    geom_errorbar(
      aes(
        ymin = Mean_Cell_Count - SEM_Cell_Count,
        ymax = Mean_Cell_Count + SEM_Cell_Count
      ),
      width = 0.2,
      linewidth = 0.4
    ) +
    geom_text(
      data = count_p_labels %>%
        filter(Fibroblast_Subtype %in% significant_count_subtypes),
      aes(
        x = 1.5,
        y = y_position,
        label = label
      ),
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold"
    ) +
    facet_wrap(
      ~Fibroblast_Subtype,
      scales = "free_y",
      ncol = 3
    ) +
    scale_fill_manual(
      values = c(
        "Female" = "#F8766D",
        "Male" = "#00BFC4"
      )
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.28))
    ) +
    labs(
      title = "Significant Fibroblast Subtype Cell Count Differences by Sex",
      x = "Sex",
      y = "Mean number of cells per sample",
      fill = "Sex"
    ) +
    theme_classic(base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 11, color = "black"),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "none"
    )
  
  print(count_plot_sig)
  
  ggsave(
    file.path(fibroblast_outdir, "Significant_Fibroblast_Subtype_Cell_Count_Differences_By_Sex.png"),
    plot = count_plot_sig,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  ggsave(
    file.path(fibroblast_outdir, "Significant_Fibroblast_Subtype_Cell_Count_Differences_By_Sex.pdf"),
    plot = count_plot_sig,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# ============================================================
# 13. FINAL MESSAGE
# ============================================================

cat("\n======================================================\n")
cat("FIBROBLAST SUBTYPE COMPOSITION T-TEST ANALYSIS COMPLETE\n")
cat("Analysis used:\n")
cat(fibroblast_path, "\n\n")
cat("Results saved to:\n")
cat(fibroblast_outdir, "\n\n")
cat("Main outputs:\n")
cat("- Fibroblast_Subtype_Counts_and_Proportions_Per_Sample.csv\n")
cat("- Ttest_Fibroblast_Subtype_Cell_Counts_By_Sex.csv\n")
cat("- Ttest_Fibroblast_Subtype_Proportions_By_Sex.csv\n")
cat("- Fibroblast_Subtype_Cell_Counts_By_Sex_Ttest_Barplot.png/pdf\n")
cat("- Fibroblast_Subtype_Proportions_By_Sex_Ttest_Barplot.png/pdf\n")
cat("- Fibroblast_Subtype_Barplot_Summary.csv\n")
cat("======================================================\n")


### Checking Collagen Type 1 and 3 Levels Across Male versus Female Samples

# ============================================================
# COLLAGEN TYPE I AND TYPE III EXPRESSION ANALYSIS
# Female vs Male samples
# Statistical unit = patient/sample
# Genes:
# Type I collagen = COL1A1 + COL1A2
# Type III collagen = COL3A1
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)

options(stringsAsFactors = FALSE)

# ============================================================
# 1. SET PATHS
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

# Default: fibroblast subset
object_path <- file.path(base_dir, "Fibroblast_Subset.RDS")

# If you want to run across the entire main UMAP instead, use this:
# object_path <- file.path(base_dir, "Merged_Joined.RDS")

collagen_outdir <- file.path(
  base_dir,
  "Collagen_Type1_Type3_Expression_Ttests"
)

dir.create(
  collagen_outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# 2. LOAD SEURAT OBJECT
# ============================================================

seurat_obj <- readRDS(object_path)

DefaultAssay(seurat_obj) <- "RNA"

seurat_obj <- tryCatch(
  JoinLayers(seurat_obj),
  error = function(e) seurat_obj
)

print("Loaded object:")
print(object_path)

print("Number of cells analyzed:")
print(ncol(seurat_obj))

# ============================================================
# 3. ADD SEX METADATA
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

seurat_obj$Sex <- unname(
  sex_mapping[as.character(seurat_obj$orig.ident)]
)

seurat_obj$Sex <- factor(
  seurat_obj$Sex,
  levels = c("Female", "Male")
)

print("Cells by sex:")
print(table(seurat_obj$Sex))

print("Cells by sample:")
print(table(seurat_obj$orig.ident, seurat_obj$Sex))

# ============================================================
# 4. DEFINE COLLAGEN GENE SETS
# ============================================================

collagen_sets <- list(
  "Type I Collagen" = c("COL1A1", "COL1A2"),
  "Type III Collagen" = c("COL3A1")
)

# ============================================================
# 5. GET NORMALIZED RNA EXPRESSION DATA
# ============================================================

get_rna_data <- function(seurat_obj) {
  
  data_input <- tryCatch(
    GetAssayData(
      seurat_obj,
      assay = "RNA",
      layer = "data"
    ),
    error = function(e) {
      GetAssayData(
        seurat_obj,
        assay = "RNA",
        slot = "data"
      )
    }
  )
  
  return(data_input)
}

expr_mat <- get_rna_data(seurat_obj)

all_collagen_genes <- unique(unlist(collagen_sets))
genes_present <- intersect(all_collagen_genes, rownames(expr_mat))
genes_missing <- setdiff(all_collagen_genes, rownames(expr_mat))

print("Collagen genes present:")
print(genes_present)

print("Collagen genes missing:")
print(genes_missing)

if (length(genes_present) == 0) {
  stop("None of the collagen genes were found in the RNA expression matrix.")
}

# ============================================================
# 6. CREATE CELL-LEVEL COLLAGEN EXPRESSION TABLE
# ============================================================

meta <- seurat_obj@meta.data
meta$Cell_ID <- rownames(meta)

meta <- meta[colnames(expr_mat), ]

cell_expression_list <- list()

for (set_name in names(collagen_sets)) {
  
  genes_use <- intersect(
    collagen_sets[[set_name]],
    rownames(expr_mat)
  )
  
  if (length(genes_use) == 0) {
    warning(paste("No genes found for:", set_name))
    next
  }
  
  set_expr <- expr_mat[genes_use, , drop = FALSE]
  
  total_set_expression_per_cell <- Matrix::colSums(set_expr)
  
  cell_expression_list[[set_name]] <- tibble(
    Cell_ID = colnames(expr_mat),
    orig.ident = meta$orig.ident,
    Sex = meta$Sex,
    Collagen_Set = set_name,
    Genes_Used = paste(genes_use, collapse = ", "),
    Total_Set_Expression_Per_Cell = as.numeric(total_set_expression_per_cell),
    Expressing = Total_Set_Expression_Per_Cell > 0
  )
}

cell_expression_df <- bind_rows(cell_expression_list)

write.csv(
  cell_expression_df,
  file.path(collagen_outdir, "Cell_Level_Collagen_Type1_Type3_Expression.csv"),
  row.names = FALSE
)

# ============================================================
# 7. SUMMARIZE COLLAGEN EXPRESSION PER SAMPLE
# ============================================================
# Total_Expression = summed expression across all analyzed cells in each sample
# Mean_Expression_Per_Cell = average expression per analyzed cell
# Percent_Expressing = percent of analyzed cells expressing the gene set

collagen_sample_summary <- cell_expression_df %>%
  filter(!is.na(Sex)) %>%
  group_by(orig.ident, Sex, Collagen_Set, Genes_Used) %>%
  summarise(
    Cells_Analyzed = n(),
    Total_Expression = sum(Total_Set_Expression_Per_Cell, na.rm = TRUE),
    Mean_Expression_Per_Cell = mean(Total_Set_Expression_Per_Cell, na.rm = TRUE),
    Median_Expression_Per_Cell = median(Total_Set_Expression_Per_Cell, na.rm = TRUE),
    Percent_Expressing = 100 * mean(Expressing, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  collagen_sample_summary,
  file.path(collagen_outdir, "Sample_Level_Collagen_Type1_Type3_Expression.csv"),
  row.names = FALSE
)

print("Sample-level collagen expression summary:")
print(collagen_sample_summary)

# ============================================================
# 8. T-TEST FUNCTION
# ============================================================

run_collagen_ttest <- function(data, value_col) {
  
  results <- data %>%
    group_by(Collagen_Set, Genes_Used) %>%
    group_modify(~{
      
      female_values <- .x[[value_col]][.x$Sex == "Female"]
      male_values <- .x[[value_col]][.x$Sex == "Male"]
      
      p_val <- tryCatch(
        {
          if (
            length(na.omit(female_values)) >= 2 &&
            length(na.omit(male_values)) >= 2
          ) {
            t.test(female_values, male_values)$p.value
          } else {
            NA_real_
          }
        },
        error = function(e) NA_real_
      )
      
      tibble(
        Metric = value_col,
        Female_Mean = mean(female_values, na.rm = TRUE),
        Male_Mean = mean(male_values, na.rm = TRUE),
        Female_SD = sd(female_values, na.rm = TRUE),
        Male_SD = sd(male_values, na.rm = TRUE),
        Female_N = length(na.omit(female_values)),
        Male_N = length(na.omit(male_values)),
        p_value = p_val
      )
    }) %>%
    ungroup()
  
  return(results)
}

collagen_ttest_results <- bind_rows(
  run_collagen_ttest(
    collagen_sample_summary,
    "Total_Expression"
  ),
  run_collagen_ttest(
    collagen_sample_summary,
    "Mean_Expression_Per_Cell"
  ),
  run_collagen_ttest(
    collagen_sample_summary,
    "Percent_Expressing"
  )
) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    Significance = case_when(
      is.na(p_value) ~ "NA",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    Significance_BH = case_when(
      is.na(p_adj_BH) ~ "NA",
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01 ~ "**",
      p_adj_BH < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

write.csv(
  collagen_ttest_results,
  file.path(collagen_outdir, "Ttest_Collagen_Type1_Type3_Expression_By_Sex.csv"),
  row.names = FALSE
)

print("Collagen expression t-test results:")
print(collagen_ttest_results)

# ============================================================
# 9. SUMMARY DATA FOR TOTAL EXPRESSION BARPLOT
# ============================================================

total_plot_summary <- collagen_sample_summary %>%
  group_by(Collagen_Set, Sex) %>%
  summarise(
    Mean_Total_Expression = mean(Total_Expression, na.rm = TRUE),
    SEM_Total_Expression = sd(Total_Expression, na.rm = TRUE) /
      sqrt(sum(!is.na(Total_Expression))),
    N = sum(!is.na(Total_Expression)),
    .groups = "drop"
  )

total_p_labels <- total_plot_summary %>%
  group_by(Collagen_Set) %>%
  summarise(
    y_position = max(
      Mean_Total_Expression + SEM_Total_Expression,
      na.rm = TRUE
    ) * 1.25,
    .groups = "drop"
  ) %>%
  left_join(
    collagen_ttest_results %>%
      filter(Metric == "Total_Expression") %>%
      dplyr::select(Collagen_Set, p_value, p_adj_BH, Significance),
    by = "Collagen_Set"
  ) %>%
  mutate(
    y_position = ifelse(
      !is.finite(y_position) | y_position == 0,
      1,
      y_position
    ),
    label = ifelse(
      is.na(p_value),
      "p = NA",
      paste0("p = ", signif(p_value, 2))
    )
  )

# ============================================================
# 10. BARPLOT: TOTAL COLLAGEN EXPRESSION
# ============================================================

total_expression_plot <- ggplot(
  total_plot_summary,
  aes(
    x = Collagen_Set,
    y = Mean_Total_Expression,
    fill = Sex
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_Total_Expression - SEM_Total_Expression,
      ymax = Mean_Total_Expression + SEM_Total_Expression
    ),
    position = position_dodge(width = 0.75),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_text(
    data = total_p_labels,
    aes(
      x = Collagen_Set,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    size = 4.5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.28))
  ) +
  labs(
    title = "Collagen Type I and Type III Total Expression by Sex",
    x = "Collagen gene set",
    y = "Mean total normalized expression per sample",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text = element_text(
      size = 11
    )
  )

print(total_expression_plot)

ggsave(
  file.path(collagen_outdir, "Collagen_Type1_Type3_Total_Expression_By_Sex_Barplot.png"),
  plot = total_expression_plot,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(collagen_outdir, "Collagen_Type1_Type3_Total_Expression_By_Sex_Barplot.pdf"),
  plot = total_expression_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# 11. SUMMARY DATA FOR MEAN PER-CELL EXPRESSION BARPLOT
# ============================================================
# This is useful because total expression can be affected by
# how many fibroblast cells were captured per sample.

mean_plot_summary <- collagen_sample_summary %>%
  group_by(Collagen_Set, Sex) %>%
  summarise(
    Mean_Per_Cell_Expression = mean(Mean_Expression_Per_Cell, na.rm = TRUE),
    SEM_Per_Cell_Expression = sd(Mean_Expression_Per_Cell, na.rm = TRUE) /
      sqrt(sum(!is.na(Mean_Expression_Per_Cell))),
    N = sum(!is.na(Mean_Expression_Per_Cell)),
    .groups = "drop"
  )

mean_p_labels <- mean_plot_summary %>%
  group_by(Collagen_Set) %>%
  summarise(
    y_position = max(
      Mean_Per_Cell_Expression + SEM_Per_Cell_Expression,
      na.rm = TRUE
    ) * 1.25,
    .groups = "drop"
  ) %>%
  left_join(
    collagen_ttest_results %>%
      filter(Metric == "Mean_Expression_Per_Cell") %>%
      dplyr::select(Collagen_Set, p_value, p_adj_BH, Significance),
    by = "Collagen_Set"
  ) %>%
  mutate(
    y_position = ifelse(
      !is.finite(y_position) | y_position == 0,
      0.05,
      y_position
    ),
    label = ifelse(
      is.na(p_value),
      "p = NA",
      paste0("p = ", signif(p_value, 2))
    )
  )

# ============================================================
# 12. BARPLOT: MEAN PER-CELL COLLAGEN EXPRESSION
# ============================================================

mean_expression_plot <- ggplot(
  mean_plot_summary,
  aes(
    x = Collagen_Set,
    y = Mean_Per_Cell_Expression,
    fill = Sex
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_Per_Cell_Expression - SEM_Per_Cell_Expression,
      ymax = Mean_Per_Cell_Expression + SEM_Per_Cell_Expression
    ),
    position = position_dodge(width = 0.75),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_text(
    data = mean_p_labels,
    aes(
      x = Collagen_Set,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    size = 4.5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.28))
  ) +
  labs(
    title = "Collagen Type I and Type III Mean Expression by Sex",
    x = "Collagen gene set",
    y = "Mean normalized expression per cell",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text = element_text(
      size = 11
    )
  )

print(mean_expression_plot)

ggsave(
  file.path(collagen_outdir, "Collagen_Type1_Type3_Mean_Per_Cell_Expression_By_Sex_Barplot.png"),
  plot = mean_expression_plot,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(collagen_outdir, "Collagen_Type1_Type3_Mean_Per_Cell_Expression_By_Sex_Barplot.pdf"),
  plot = mean_expression_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# 13. FINAL MESSAGE
# ============================================================

cat("\n======================================================\n")
cat("COLLAGEN TYPE I / TYPE III EXPRESSION ANALYSIS COMPLETE\n")
cat("Analysis used:\n")
cat(object_path, "\n\n")
cat("Results saved to:\n")
cat(collagen_outdir, "\n\n")
cat("Main outputs:\n")
cat("- Sample_Level_Collagen_Type1_Type3_Expression.csv\n")
cat("- Ttest_Collagen_Type1_Type3_Expression_By_Sex.csv\n")
cat("- Collagen_Type1_Type3_Total_Expression_By_Sex_Barplot.png/pdf\n")
cat("- Collagen_Type1_Type3_Mean_Per_Cell_Expression_By_Sex_Barplot.png/pdf\n")
cat("======================================================\n")



### Overall Collagen Expression Male versus Female

# ============================================================
# OVERALL COLLAGEN GENE EXPRESSION ANALYSIS
# Female vs Male samples
# Statistical unit = patient/sample
# Default object = Fibroblast_Subset.RDS
# Overall collagen expression = summed normalized RNA expression
# across all collagen genes present in the object
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)

options(stringsAsFactors = FALSE)

# ============================================================
# 1. SET PATHS
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

# Use fibroblast subset by default because collagen is mainly CAF/fibroblast-associated
object_path <- file.path(base_dir, "Fibroblast_Subset.RDS")

# If you want to run this across the whole UMAP instead, use:
# object_path <- file.path(base_dir, "Merged_Joined.RDS")

collagen_outdir <- file.path(
  base_dir,
  "Overall_Collagen_Expression_Ttests"
)

dir.create(
  collagen_outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# 2. LOAD SEURAT OBJECT
# ============================================================

seurat_obj <- readRDS(object_path)

DefaultAssay(seurat_obj) <- "RNA"

seurat_obj <- tryCatch(
  JoinLayers(seurat_obj),
  error = function(e) seurat_obj
)

print("Loaded object:")
print(object_path)

print("Number of cells analyzed:")
print(ncol(seurat_obj))

# ============================================================
# 3. ADD SEX METADATA
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

seurat_obj$Sex <- unname(
  sex_mapping[as.character(seurat_obj$orig.ident)]
)

seurat_obj$Sex <- factor(
  seurat_obj$Sex,
  levels = c("Female", "Male")
)

print("Cells by sex:")
print(table(seurat_obj$Sex))

print("Cells by sample:")
print(table(seurat_obj$orig.ident, seurat_obj$Sex))

# ============================================================
# 4. DEFINE OVERALL COLLAGEN GENE SET
# ============================================================
# This includes the main human collagen alpha-chain genes.
# The code automatically keeps only genes present in your RNA matrix.

all_collagen_genes <- c(
  "COL1A1", "COL1A2",
  "COL2A1",
  "COL3A1",
  "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6",
  "COL5A1", "COL5A2", "COL5A3",
  "COL6A1", "COL6A2", "COL6A3", "COL6A4", "COL6A5", "COL6A6",
  "COL7A1",
  "COL8A1", "COL8A2",
  "COL9A1", "COL9A2", "COL9A3",
  "COL10A1",
  "COL11A1", "COL11A2",
  "COL12A1",
  "COL13A1",
  "COL14A1",
  "COL15A1",
  "COL16A1",
  "COL17A1",
  "COL18A1",
  "COL19A1",
  "COL20A1",
  "COL21A1",
  "COL22A1",
  "COL23A1",
  "COL24A1",
  "COL25A1",
  "COL26A1",
  "COL27A1",
  "COL28A1"
)

# ============================================================
# 5. GET NORMALIZED RNA EXPRESSION DATA
# ============================================================

get_rna_data <- function(seurat_obj) {
  
  data_input <- tryCatch(
    GetAssayData(
      seurat_obj,
      assay = "RNA",
      layer = "data"
    ),
    error = function(e) {
      GetAssayData(
        seurat_obj,
        assay = "RNA",
        slot = "data"
      )
    }
  )
  
  return(data_input)
}

expr_mat <- get_rna_data(seurat_obj)

collagen_genes_present <- intersect(
  all_collagen_genes,
  rownames(expr_mat)
)

collagen_genes_missing <- setdiff(
  all_collagen_genes,
  rownames(expr_mat)
)

print("Collagen genes present in object:")
print(collagen_genes_present)

print("Collagen genes missing from object:")
print(collagen_genes_missing)

if (length(collagen_genes_present) == 0) {
  stop("No collagen genes were found in the RNA expression matrix.")
}

write.csv(
  data.frame(
    Collagen_Genes_Present = collagen_genes_present
  ),
  file.path(collagen_outdir, "Collagen_Genes_Used.csv"),
  row.names = FALSE
)

# ============================================================
# 6. CALCULATE CELL-LEVEL OVERALL COLLAGEN EXPRESSION
# ============================================================

meta <- seurat_obj@meta.data
meta$Cell_ID <- rownames(meta)

meta <- meta[colnames(expr_mat), ]

overall_collagen_expression_per_cell <- Matrix::colSums(
  expr_mat[collagen_genes_present, , drop = FALSE]
)

cell_expression_df <- tibble(
  Cell_ID = colnames(expr_mat),
  orig.ident = meta$orig.ident,
  Sex = meta$Sex,
  Overall_Collagen_Expression = as.numeric(overall_collagen_expression_per_cell),
  Expressing_Collagen = Overall_Collagen_Expression > 0
)

write.csv(
  cell_expression_df,
  file.path(collagen_outdir, "Cell_Level_Overall_Collagen_Expression.csv"),
  row.names = FALSE
)

# ============================================================
# 7. SUMMARIZE OVERALL COLLAGEN EXPRESSION PER SAMPLE
# ============================================================
# Total_Expression = summed collagen expression across all analyzed cells in each sample
# Mean_Expression_Per_Cell = average collagen expression per analyzed cell
# Percent_Collagen_Expressing = percent of analyzed cells with collagen expression > 0

collagen_sample_summary <- cell_expression_df %>%
  filter(!is.na(Sex)) %>%
  group_by(orig.ident, Sex) %>%
  summarise(
    Cells_Analyzed = n(),
    Total_Collagen_Expression = sum(Overall_Collagen_Expression, na.rm = TRUE),
    Mean_Collagen_Expression_Per_Cell = mean(Overall_Collagen_Expression, na.rm = TRUE),
    Median_Collagen_Expression_Per_Cell = median(Overall_Collagen_Expression, na.rm = TRUE),
    Percent_Collagen_Expressing = 100 * mean(Expressing_Collagen, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  collagen_sample_summary,
  file.path(collagen_outdir, "Sample_Level_Overall_Collagen_Expression.csv"),
  row.names = FALSE
)

print("Sample-level overall collagen expression:")
print(collagen_sample_summary)

# ============================================================
# 8. RUN STATISTICAL T-TESTS
# ============================================================

run_collagen_ttest <- function(data, value_col) {
  
  female_values <- data[[value_col]][data$Sex == "Female"]
  male_values <- data[[value_col]][data$Sex == "Male"]
  
  p_val <- tryCatch(
    {
      if (
        length(na.omit(female_values)) >= 2 &&
        length(na.omit(male_values)) >= 2
      ) {
        t.test(female_values, male_values)$p.value
      } else {
        NA_real_
      }
    },
    error = function(e) NA_real_
  )
  
  tibble(
    Metric = value_col,
    Female_Mean = mean(female_values, na.rm = TRUE),
    Male_Mean = mean(male_values, na.rm = TRUE),
    Female_SD = sd(female_values, na.rm = TRUE),
    Male_SD = sd(male_values, na.rm = TRUE),
    Female_N = length(na.omit(female_values)),
    Male_N = length(na.omit(male_values)),
    p_value = p_val
  )
}

collagen_ttest_results <- bind_rows(
  run_collagen_ttest(
    collagen_sample_summary,
    "Total_Collagen_Expression"
  ),
  run_collagen_ttest(
    collagen_sample_summary,
    "Mean_Collagen_Expression_Per_Cell"
  ),
  run_collagen_ttest(
    collagen_sample_summary,
    "Percent_Collagen_Expressing"
  )
) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    Significance = case_when(
      is.na(p_value) ~ "NA",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    Significance_BH = case_when(
      is.na(p_adj_BH) ~ "NA",
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01 ~ "**",
      p_adj_BH < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

write.csv(
  collagen_ttest_results,
  file.path(collagen_outdir, "Ttest_Overall_Collagen_Expression_By_Sex.csv"),
  row.names = FALSE
)

print("Overall collagen expression t-test results:")
print(collagen_ttest_results)

# ============================================================
# 9. SUMMARY DATA FOR TOTAL COLLAGEN BARPLOT
# ============================================================

total_plot_summary <- collagen_sample_summary %>%
  group_by(Sex) %>%
  summarise(
    Mean_Total_Collagen_Expression = mean(Total_Collagen_Expression, na.rm = TRUE),
    SEM_Total_Collagen_Expression = sd(Total_Collagen_Expression, na.rm = TRUE) /
      sqrt(sum(!is.na(Total_Collagen_Expression))),
    N = sum(!is.na(Total_Collagen_Expression)),
    .groups = "drop"
  )

total_p_value <- collagen_ttest_results %>%
  filter(Metric == "Total_Collagen_Expression") %>%
  pull(p_value)

total_p_label <- ifelse(
  is.na(total_p_value),
  "p = NA",
  paste0("p = ", signif(total_p_value, 2))
)

total_y_position <- max(
  total_plot_summary$Mean_Total_Collagen_Expression +
    total_plot_summary$SEM_Total_Collagen_Expression,
  na.rm = TRUE
) * 1.25

if (!is.finite(total_y_position) || total_y_position == 0) {
  total_y_position <- 1
}

# ============================================================
# 10. BARPLOT: TOTAL OVERALL COLLAGEN EXPRESSION
# ============================================================

total_collagen_plot <- ggplot(
  total_plot_summary,
  aes(
    x = Sex,
    y = Mean_Total_Collagen_Expression,
    fill = Sex
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_Total_Collagen_Expression - SEM_Total_Collagen_Expression,
      ymax = Mean_Total_Collagen_Expression + SEM_Total_Collagen_Expression
    ),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      x = 1.5,
      y = total_y_position,
      label = total_p_label
    ),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.28))
  ) +
  labs(
    title = "Overall Collagen Gene Expression by Sex",
    x = "Sex",
    y = "Mean total normalized collagen expression per sample",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "none"
  )

print(total_collagen_plot)

ggsave(
  file.path(collagen_outdir, "Overall_Collagen_Total_Expression_By_Sex_Barplot.png"),
  plot = total_collagen_plot,
  width = 6,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(collagen_outdir, "Overall_Collagen_Total_Expression_By_Sex_Barplot.pdf"),
  plot = total_collagen_plot,
  width = 6,
  height = 5,
  dpi = 300
)

# ============================================================
# 11. SUMMARY DATA FOR MEAN PER-CELL COLLAGEN BARPLOT
# ============================================================
# This second plot controls for differences in number of cells captured per sample.

mean_plot_summary <- collagen_sample_summary %>%
  group_by(Sex) %>%
  summarise(
    Mean_Collagen_Per_Cell = mean(Mean_Collagen_Expression_Per_Cell, na.rm = TRUE),
    SEM_Collagen_Per_Cell = sd(Mean_Collagen_Expression_Per_Cell, na.rm = TRUE) /
      sqrt(sum(!is.na(Mean_Collagen_Expression_Per_Cell))),
    N = sum(!is.na(Mean_Collagen_Expression_Per_Cell)),
    .groups = "drop"
  )

mean_p_value <- collagen_ttest_results %>%
  filter(Metric == "Mean_Collagen_Expression_Per_Cell") %>%
  pull(p_value)

mean_p_label <- ifelse(
  is.na(mean_p_value),
  "p = NA",
  paste0("p = ", signif(mean_p_value, 2))
)

mean_y_position <- max(
  mean_plot_summary$Mean_Collagen_Per_Cell +
    mean_plot_summary$SEM_Collagen_Per_Cell,
  na.rm = TRUE
) * 1.25

if (!is.finite(mean_y_position) || mean_y_position == 0) {
  mean_y_position <- 0.05
}

# ============================================================
# 12. BARPLOT: MEAN PER-CELL OVERALL COLLAGEN EXPRESSION
# ============================================================

mean_collagen_plot <- ggplot(
  mean_plot_summary,
  aes(
    x = Sex,
    y = Mean_Collagen_Per_Cell,
    fill = Sex
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_Collagen_Per_Cell - SEM_Collagen_Per_Cell,
      ymax = Mean_Collagen_Per_Cell + SEM_Collagen_Per_Cell
    ),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      x = 1.5,
      y = mean_y_position,
      label = mean_p_label
    ),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.28))
  ) +
  labs(
    title = "Mean Overall Collagen Expression by Sex",
    x = "Sex",
    y = "Mean normalized collagen expression per cell",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "none"
  )

print(mean_collagen_plot)

ggsave(
  file.path(collagen_outdir, "Overall_Collagen_Mean_Per_Cell_Expression_By_Sex_Barplot.png"),
  plot = mean_collagen_plot,
  width = 6,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(collagen_outdir, "Overall_Collagen_Mean_Per_Cell_Expression_By_Sex_Barplot.pdf"),
  plot = mean_collagen_plot,
  width = 6,
  height = 5,
  dpi = 300
)

# ============================================================
# 13. FINAL MESSAGE
# ============================================================

cat("\n======================================================\n")
cat("OVERALL COLLAGEN EXPRESSION ANALYSIS COMPLETE\n")
cat("Analysis used:\n")
cat(object_path, "\n\n")
cat("Results saved to:\n")
cat(collagen_outdir, "\n\n")
cat("Collagen genes used:\n")
cat(paste(collagen_genes_present, collapse = ", "), "\n\n")
cat("Main outputs:\n")
cat("- Collagen_Genes_Used.csv\n")
cat("- Cell_Level_Overall_Collagen_Expression.csv\n")
cat("- Sample_Level_Overall_Collagen_Expression.csv\n")
cat("- Ttest_Overall_Collagen_Expression_By_Sex.csv\n")
cat("- Overall_Collagen_Total_Expression_By_Sex_Barplot.png/pdf\n")
cat("- Overall_Collagen_Mean_Per_Cell_Expression_By_Sex_Barplot.png/pdf\n")
cat("======================================================\n")



### Collagen Type 3 Male verus Female

# ============================================================
# COLLAGEN TYPE III / COL3A1 EXPRESSION ANALYSIS
# Female vs Male samples
# Statistical unit = patient/sample
# Default object = Fibroblast_Subset.RDS
# Gene tested = COL3A1
# Includes t-test p-value on barplots
# No individual dots/points
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)

options(stringsAsFactors = FALSE)

# ============================================================
# 1. SET PATHS
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

# Use fibroblast subset by default because COL3A1 is mainly fibroblast/CAF-associated
object_path <- file.path(base_dir, "Fibroblast_Subset.RDS")

# If you want to run this across the whole UMAP instead, use:
# object_path <- file.path(base_dir, "Merged_Joined.RDS")

col3a1_outdir <- file.path(
  base_dir,
  "COL3A1_Collagen_Type3_Expression_Ttests"
)

dir.create(
  col3a1_outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# 2. LOAD SEURAT OBJECT
# ============================================================

seurat_obj <- readRDS(object_path)

DefaultAssay(seurat_obj) <- "RNA"

seurat_obj <- tryCatch(
  JoinLayers(seurat_obj),
  error = function(e) seurat_obj
)

print("Loaded object:")
print(object_path)

print("Number of cells analyzed:")
print(ncol(seurat_obj))

# ============================================================
# 3. ADD SEX METADATA
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

seurat_obj$Sex <- unname(
  sex_mapping[as.character(seurat_obj$orig.ident)]
)

seurat_obj$Sex <- factor(
  seurat_obj$Sex,
  levels = c("Female", "Male")
)

print("Cells by sex:")
print(table(seurat_obj$Sex))

print("Cells by sample:")
print(table(seurat_obj$orig.ident, seurat_obj$Sex))

# ============================================================
# 4. GET NORMALIZED RNA EXPRESSION DATA
# ============================================================

get_rna_data <- function(seurat_obj) {
  
  data_input <- tryCatch(
    GetAssayData(
      seurat_obj,
      assay = "RNA",
      layer = "data"
    ),
    error = function(e) {
      GetAssayData(
        seurat_obj,
        assay = "RNA",
        slot = "data"
      )
    }
  )
  
  return(data_input)
}

expr_mat <- get_rna_data(seurat_obj)

# ============================================================
# 5. CHECK COL3A1 PRESENCE
# ============================================================

gene_use <- "COL3A1"

if (!(gene_use %in% rownames(expr_mat))) {
  stop("COL3A1 was not found in the RNA expression matrix.")
}

print("Gene used:")
print(gene_use)

# ============================================================
# 6. CREATE CELL-LEVEL COL3A1 EXPRESSION TABLE
# ============================================================

meta <- seurat_obj@meta.data
meta$Cell_ID <- rownames(meta)

meta <- meta[colnames(expr_mat), ]

col3a1_expression_per_cell <- as.numeric(
  expr_mat[gene_use, ]
)

cell_expression_df <- tibble(
  Cell_ID = colnames(expr_mat),
  orig.ident = meta$orig.ident,
  Sex = meta$Sex,
  COL3A1_Expression = col3a1_expression_per_cell,
  Expressing_COL3A1 = COL3A1_Expression > 0
)

write.csv(
  cell_expression_df,
  file.path(col3a1_outdir, "Cell_Level_COL3A1_Expression.csv"),
  row.names = FALSE
)

# ============================================================
# 7. SUMMARIZE COL3A1 EXPRESSION PER SAMPLE
# ============================================================
# Total_COL3A1_Expression = summed COL3A1 expression across cells in each sample
# Mean_COL3A1_Expression_Per_Cell = average COL3A1 expression per analyzed cell
# Percent_COL3A1_Expressing = percent of cells with COL3A1 expression > 0

col3a1_sample_summary <- cell_expression_df %>%
  filter(!is.na(Sex)) %>%
  group_by(orig.ident, Sex) %>%
  summarise(
    Cells_Analyzed = n(),
    Total_COL3A1_Expression = sum(COL3A1_Expression, na.rm = TRUE),
    Mean_COL3A1_Expression_Per_Cell = mean(COL3A1_Expression, na.rm = TRUE),
    Median_COL3A1_Expression_Per_Cell = median(COL3A1_Expression, na.rm = TRUE),
    Percent_COL3A1_Expressing = 100 * mean(Expressing_COL3A1, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  col3a1_sample_summary,
  file.path(col3a1_outdir, "Sample_Level_COL3A1_Expression.csv"),
  row.names = FALSE
)

print("Sample-level COL3A1 expression:")
print(col3a1_sample_summary)

# ============================================================
# 8. RUN STATISTICAL T-TESTS
# ============================================================

run_col3a1_ttest <- function(data, value_col) {
  
  female_values <- data[[value_col]][data$Sex == "Female"]
  male_values <- data[[value_col]][data$Sex == "Male"]
  
  p_val <- tryCatch(
    {
      if (
        length(na.omit(female_values)) >= 2 &&
        length(na.omit(male_values)) >= 2
      ) {
        t.test(female_values, male_values)$p.value
      } else {
        NA_real_
      }
    },
    error = function(e) NA_real_
  )
  
  tibble(
    Metric = value_col,
    Female_Mean = mean(female_values, na.rm = TRUE),
    Male_Mean = mean(male_values, na.rm = TRUE),
    Female_SD = sd(female_values, na.rm = TRUE),
    Male_SD = sd(male_values, na.rm = TRUE),
    Female_N = length(na.omit(female_values)),
    Male_N = length(na.omit(male_values)),
    p_value = p_val
  )
}

col3a1_ttest_results <- bind_rows(
  run_col3a1_ttest(
    col3a1_sample_summary,
    "Total_COL3A1_Expression"
  ),
  run_col3a1_ttest(
    col3a1_sample_summary,
    "Mean_COL3A1_Expression_Per_Cell"
  ),
  run_col3a1_ttest(
    col3a1_sample_summary,
    "Percent_COL3A1_Expressing"
  )
) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    Significance = case_when(
      is.na(p_value) ~ "NA",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    Significance_BH = case_when(
      is.na(p_adj_BH) ~ "NA",
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01 ~ "**",
      p_adj_BH < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

write.csv(
  col3a1_ttest_results,
  file.path(col3a1_outdir, "Ttest_COL3A1_Expression_By_Sex.csv"),
  row.names = FALSE
)

print("COL3A1 expression t-test results:")
print(col3a1_ttest_results)

# ============================================================
# 9. SUMMARY DATA FOR TOTAL COL3A1 BARPLOT
# ============================================================

total_plot_summary <- col3a1_sample_summary %>%
  group_by(Sex) %>%
  summarise(
    Mean_Total_COL3A1_Expression = mean(Total_COL3A1_Expression, na.rm = TRUE),
    SEM_Total_COL3A1_Expression = sd(Total_COL3A1_Expression, na.rm = TRUE) /
      sqrt(sum(!is.na(Total_COL3A1_Expression))),
    N = sum(!is.na(Total_COL3A1_Expression)),
    .groups = "drop"
  )

total_p_value <- col3a1_ttest_results %>%
  filter(Metric == "Total_COL3A1_Expression") %>%
  pull(p_value)

total_p_label <- ifelse(
  is.na(total_p_value),
  "p = NA",
  paste0("p = ", signif(total_p_value, 2))
)

total_y_position <- max(
  total_plot_summary$Mean_Total_COL3A1_Expression +
    total_plot_summary$SEM_Total_COL3A1_Expression,
  na.rm = TRUE
) * 1.25

if (!is.finite(total_y_position) || total_y_position == 0) {
  total_y_position <- 1
}

# ============================================================
# 10. BARPLOT: TOTAL COL3A1 EXPRESSION
# ============================================================

total_col3a1_plot <- ggplot(
  total_plot_summary,
  aes(
    x = Sex,
    y = Mean_Total_COL3A1_Expression,
    fill = Sex
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_Total_COL3A1_Expression - SEM_Total_COL3A1_Expression,
      ymax = Mean_Total_COL3A1_Expression + SEM_Total_COL3A1_Expression
    ),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      x = 1.5,
      y = total_y_position,
      label = total_p_label
    ),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.28))
  ) +
  labs(
    title = "Collagen Type III Expression by Sex",
    x = "Sex",
    y = "Mean total normalized COL3A1 expression per sample",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "none"
  )

print(total_col3a1_plot)

ggsave(
  file.path(col3a1_outdir, "COL3A1_Total_Expression_By_Sex_Barplot.png"),
  plot = total_col3a1_plot,
  width = 6,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(col3a1_outdir, "COL3A1_Total_Expression_By_Sex_Barplot.pdf"),
  plot = total_col3a1_plot,
  width = 6,
  height = 5,
  dpi = 300
)

# ============================================================
# 11. SUMMARY DATA FOR MEAN PER-CELL COL3A1 BARPLOT
# ============================================================
# This second plot controls for differences in the number of cells captured per sample.

mean_plot_summary <- col3a1_sample_summary %>%
  group_by(Sex) %>%
  summarise(
    Mean_COL3A1_Per_Cell = mean(Mean_COL3A1_Expression_Per_Cell, na.rm = TRUE),
    SEM_COL3A1_Per_Cell = sd(Mean_COL3A1_Expression_Per_Cell, na.rm = TRUE) /
      sqrt(sum(!is.na(Mean_COL3A1_Expression_Per_Cell))),
    N = sum(!is.na(Mean_COL3A1_Expression_Per_Cell)),
    .groups = "drop"
  )

mean_p_value <- col3a1_ttest_results %>%
  filter(Metric == "Mean_COL3A1_Expression_Per_Cell") %>%
  pull(p_value)

mean_p_label <- ifelse(
  is.na(mean_p_value),
  "p = NA",
  paste0("p = ", signif(mean_p_value, 2))
)

mean_y_position <- max(
  mean_plot_summary$Mean_COL3A1_Per_Cell +
    mean_plot_summary$SEM_COL3A1_Per_Cell,
  na.rm = TRUE
) * 1.25

if (!is.finite(mean_y_position) || mean_y_position == 0) {
  mean_y_position <- 0.05
}

# ============================================================
# 12. BARPLOT: MEAN PER-CELL COL3A1 EXPRESSION
# ============================================================

mean_col3a1_plot <- ggplot(
  mean_plot_summary,
  aes(
    x = Sex,
    y = Mean_COL3A1_Per_Cell,
    fill = Sex
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_COL3A1_Per_Cell - SEM_COL3A1_Per_Cell,
      ymax = Mean_COL3A1_Per_Cell + SEM_COL3A1_Per_Cell
    ),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      x = 1.5,
      y = mean_y_position,
      label = mean_p_label
    ),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.28))
  ) +
  labs(
    title = "Mean Collagen Type III Expression by Sex",
    x = "Sex",
    y = "Mean normalized COL3A1 expression per cell",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "none"
  )

print(mean_col3a1_plot)

ggsave(
  file.path(col3a1_outdir, "COL3A1_Mean_Per_Cell_Expression_By_Sex_Barplot.png"),
  plot = mean_col3a1_plot,
  width = 6,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(col3a1_outdir, "COL3A1_Mean_Per_Cell_Expression_By_Sex_Barplot.pdf"),
  plot = mean_col3a1_plot,
  width = 6,
  height = 5,
  dpi = 300
)

# ============================================================
# 13. FINAL MESSAGE
# ============================================================

cat("\n======================================================\n")
cat("COLLAGEN TYPE III / COL3A1 EXPRESSION ANALYSIS COMPLETE\n")
cat("Analysis used:\n")
cat(object_path, "\n\n")
cat("Gene used: COL3A1\n\n")
cat("Results saved to:\n")
cat(col3a1_outdir, "\n\n")
cat("Main outputs:\n")
cat("- Cell_Level_COL3A1_Expression.csv\n")
cat("- Sample_Level_COL3A1_Expression.csv\n")
cat("- Ttest_COL3A1_Expression_By_Sex.csv\n")
cat("- COL3A1_Total_Expression_By_Sex_Barplot.png/pdf\n")
cat("- COL3A1_Mean_Per_Cell_Expression_By_Sex_Barplot.png/pdf\n")
cat("======================================================\n")


### Type 3 Collagen Per Total Collagen P Value

# ============================================================
# TYPE III COLLAGEN AS A FRACTION OF TOTAL COLLAGEN EXPRESSION
# Female vs Male samples
# Metric = COL3A1 expression / total collagen gene expression
# Statistical unit = patient/sample
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)

options(stringsAsFactors = FALSE)

# ============================================================
# 1. SET PATHS
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

# Use fibroblast subset by default
object_path <- file.path(base_dir, "Fibroblast_Subset.RDS")

# If you want whole UMAP instead, use:
# object_path <- file.path(base_dir, "Merged_Joined.RDS")

outdir <- file.path(
  base_dir,
  "COL3A1_Per_Total_Collagen_Ttest"
)

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# 2. LOAD OBJECT
# ============================================================

seurat_obj <- readRDS(object_path)

DefaultAssay(seurat_obj) <- "RNA"

seurat_obj <- tryCatch(
  JoinLayers(seurat_obj),
  error = function(e) seurat_obj
)

print("Loaded object:")
print(object_path)

print("Number of cells analyzed:")
print(ncol(seurat_obj))

# ============================================================
# 3. ADD SEX METADATA
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

seurat_obj$Sex <- unname(
  sex_mapping[as.character(seurat_obj$orig.ident)]
)

seurat_obj$Sex <- factor(
  seurat_obj$Sex,
  levels = c("Female", "Male")
)

# ============================================================
# 4. DEFINE TOTAL COLLAGEN GENE SET
# ============================================================

total_collagen_genes <- c(
  "COL1A1", "COL1A2",
  "COL2A1",
  "COL3A1",
  "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6",
  "COL5A1", "COL5A2", "COL5A3",
  "COL6A1", "COL6A2", "COL6A3", "COL6A4", "COL6A5", "COL6A6",
  "COL7A1",
  "COL8A1", "COL8A2",
  "COL9A1", "COL9A2", "COL9A3",
  "COL10A1",
  "COL11A1", "COL11A2",
  "COL12A1",
  "COL13A1",
  "COL14A1",
  "COL15A1",
  "COL16A1",
  "COL17A1",
  "COL18A1",
  "COL19A1",
  "COL20A1",
  "COL21A1",
  "COL22A1",
  "COL23A1",
  "COL24A1",
  "COL25A1",
  "COL26A1",
  "COL27A1",
  "COL28A1"
)

type3_gene <- "COL3A1"

# ============================================================
# 5. GET NORMALIZED RNA DATA
# ============================================================

get_rna_data <- function(seurat_obj) {
  
  data_input <- tryCatch(
    GetAssayData(
      seurat_obj,
      assay = "RNA",
      layer = "data"
    ),
    error = function(e) {
      GetAssayData(
        seurat_obj,
        assay = "RNA",
        slot = "data"
      )
    }
  )
  
  return(data_input)
}

expr_mat <- get_rna_data(seurat_obj)

collagen_genes_present <- intersect(
  total_collagen_genes,
  rownames(expr_mat)
)

if (!(type3_gene %in% rownames(expr_mat))) {
  stop("COL3A1 was not found in the RNA expression matrix.")
}

if (length(collagen_genes_present) == 0) {
  stop("No collagen genes were found in the RNA expression matrix.")
}

print("Total collagen genes used:")
print(collagen_genes_present)

write.csv(
  data.frame(Collagen_Genes_Used = collagen_genes_present),
  file.path(outdir, "Total_Collagen_Genes_Used.csv"),
  row.names = FALSE
)

# ============================================================
# 6. CALCULATE CELL-LEVEL EXPRESSION
# ============================================================

meta <- seurat_obj@meta.data
meta$Cell_ID <- rownames(meta)

meta <- meta[colnames(expr_mat), ]

col3a1_expr_per_cell <- as.numeric(
  expr_mat[type3_gene, ]
)

total_collagen_expr_per_cell <- Matrix::colSums(
  expr_mat[collagen_genes_present, , drop = FALSE]
)

cell_level_df <- tibble(
  Cell_ID = colnames(expr_mat),
  orig.ident = meta$orig.ident,
  Sex = meta$Sex,
  COL3A1_Expression = col3a1_expr_per_cell,
  Total_Collagen_Expression = as.numeric(total_collagen_expr_per_cell),
  COL3A1_Fraction_of_Total_Collagen = ifelse(
    Total_Collagen_Expression > 0,
    COL3A1_Expression / Total_Collagen_Expression,
    NA_real_
  ),
  COL3A1_Percent_of_Total_Collagen = 100 * COL3A1_Fraction_of_Total_Collagen
)

write.csv(
  cell_level_df,
  file.path(outdir, "Cell_Level_COL3A1_Per_Total_Collagen.csv"),
  row.names = FALSE
)

# ============================================================
# 7. SUMMARIZE PER SAMPLE
# ============================================================
# This is the key step:
# For each sample:
# COL3A1_Percent_of_Total_Collagen =
# 100 * total COL3A1 expression / total collagen expression

sample_level_df <- cell_level_df %>%
  filter(!is.na(Sex)) %>%
  group_by(orig.ident, Sex) %>%
  summarise(
    Cells_Analyzed = n(),
    Total_COL3A1_Expression = sum(COL3A1_Expression, na.rm = TRUE),
    Total_Collagen_Expression = sum(Total_Collagen_Expression, na.rm = TRUE),
    COL3A1_Fraction_of_Total_Collagen = ifelse(
      Total_Collagen_Expression > 0,
      Total_COL3A1_Expression / Total_Collagen_Expression,
      NA_real_
    ),
    COL3A1_Percent_of_Total_Collagen = 100 * COL3A1_Fraction_of_Total_Collagen,
    .groups = "drop"
  )

write.csv(
  sample_level_df,
  file.path(outdir, "Sample_Level_COL3A1_Per_Total_Collagen.csv"),
  row.names = FALSE
)

print("Sample-level COL3A1 per total collagen:")
print(sample_level_df)

# ============================================================
# 8. RUN T-TEST
# ============================================================

female_values <- sample_level_df$COL3A1_Percent_of_Total_Collagen[
  sample_level_df$Sex == "Female"
]

male_values <- sample_level_df$COL3A1_Percent_of_Total_Collagen[
  sample_level_df$Sex == "Male"
]

col3a1_ratio_ttest <- t.test(
  female_values,
  male_values
)

ttest_results <- tibble(
  Metric = "COL3A1_Percent_of_Total_Collagen",
  Female_Mean = mean(female_values, na.rm = TRUE),
  Male_Mean = mean(male_values, na.rm = TRUE),
  Female_SD = sd(female_values, na.rm = TRUE),
  Male_SD = sd(male_values, na.rm = TRUE),
  Female_N = length(na.omit(female_values)),
  Male_N = length(na.omit(male_values)),
  p_value = col3a1_ratio_ttest$p.value,
  Significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "ns"
  )
)

write.csv(
  ttest_results,
  file.path(outdir, "Ttest_COL3A1_Per_Total_Collagen_By_Sex.csv"),
  row.names = FALSE
)

print("T-test result:")
print(ttest_results)

print("Raw t-test output:")
print(col3a1_ratio_ttest)

# ============================================================
# 9. BARPLOT WITH P-VALUE
# ============================================================

plot_summary <- sample_level_df %>%
  group_by(Sex) %>%
  summarise(
    Mean_COL3A1_Percent = mean(COL3A1_Percent_of_Total_Collagen, na.rm = TRUE),
    SEM_COL3A1_Percent = sd(COL3A1_Percent_of_Total_Collagen, na.rm = TRUE) /
      sqrt(sum(!is.na(COL3A1_Percent_of_Total_Collagen))),
    N = sum(!is.na(COL3A1_Percent_of_Total_Collagen)),
    .groups = "drop"
  )

p_label <- paste0(
  "p = ",
  signif(ttest_results$p_value, 2)
)

y_position <- max(
  plot_summary$Mean_COL3A1_Percent + plot_summary$SEM_COL3A1_Percent,
  na.rm = TRUE
) * 1.25

if (!is.finite(y_position) || y_position == 0) {
  y_position <- 0.05
}

col3a1_ratio_plot <- ggplot(
  plot_summary,
  aes(
    x = Sex,
    y = Mean_COL3A1_Percent,
    fill = Sex
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_COL3A1_Percent - SEM_COL3A1_Percent,
      ymax = Mean_COL3A1_Percent + SEM_COL3A1_Percent
    ),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      x = 1.5,
      y = y_position,
      label = p_label
    ),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.28))
  ) +
  labs(
    title = "Type III Collagen Fraction of Total Collagen by Sex",
    x = "Sex",
    y = "COL3A1 expression as percent of total collagen expression",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 17
    ),
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "none"
  )

print(col3a1_ratio_plot)

ggsave(
  file.path(outdir, "COL3A1_Per_Total_Collagen_By_Sex_Barplot.png"),
  plot = col3a1_ratio_plot,
  width = 6.5,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(outdir, "COL3A1_Per_Total_Collagen_By_Sex_Barplot.pdf"),
  plot = col3a1_ratio_plot,
  width = 6.5,
  height = 5,
  dpi = 300
)

# ============================================================
# 10. FINAL MESSAGE
# ============================================================

cat("\n======================================================\n")
cat("COL3A1 PER TOTAL COLLAGEN ANALYSIS COMPLETE\n")
cat("Analysis used:\n")
cat(object_path, "\n\n")
cat("Metric tested:\n")
cat("100 * Total COL3A1 expression / Total collagen expression per sample\n\n")
cat("P-value:\n")
cat(ttest_results$p_value, "\n\n")
cat("Results saved to:\n")
cat(outdir, "\n")




### Collagen Table

# ============================================================
# RAW COLLAGEN VALUES PER SAMPLE + P-VALUE TABLE
# Saves CSV files to Downloads
#
# Output 1:
# Raw collagen values for each of the 17 samples
#
# Output 2:
# Female vs Male t-test p-values for each collagen metric
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(Matrix)

options(stringsAsFactors = FALSE)

# ============================================================
# 1. LOAD OBJECT
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

# Default: fibroblast subset
object_path <- file.path(base_dir, "Fibroblast_Subset.RDS")

# To run on the full UMAP instead, use:
# object_path <- file.path(base_dir, "Merged_Joined.RDS")

seurat_obj <- readRDS(object_path)

DefaultAssay(seurat_obj) <- "RNA"

seurat_obj <- tryCatch(
  JoinLayers(seurat_obj),
  error = function(e) seurat_obj
)

cat("\nLoaded object:\n")
cat(object_path, "\n")

cat("\nNumber of cells analyzed:\n")
print(ncol(seurat_obj))

# ============================================================
# 2. SET DOWNLOADS OUTPUT PATH
# ============================================================

downloads_dir <- path.expand("~/Downloads")

raw_sample_csv <- file.path(
  downloads_dir,
  "Raw_Collagen_Values_Per_Sample.csv"
)

pvalue_csv <- file.path(
  downloads_dir,
  "Collagen_Pvalue_Summary_Table.csv"
)

# ============================================================
# 3. ADD SEX METADATA
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

seurat_obj$Sex <- unname(
  sex_mapping[as.character(seurat_obj$orig.ident)]
)

seurat_obj$Sex <- factor(
  seurat_obj$Sex,
  levels = c("Female", "Male")
)

cat("\nCells by sex:\n")
print(table(seurat_obj$Sex))

cat("\nCells by sample and sex:\n")
print(table(seurat_obj$orig.ident, seurat_obj$Sex))

# ============================================================
# 4. DEFINE COLLAGEN GENE SETS
# ============================================================

type1_genes <- c("COL1A1", "COL1A2")
type3_genes <- c("COL3A1")

total_collagen_genes <- c(
  "COL1A1", "COL1A2",
  "COL2A1",
  "COL3A1",
  "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6",
  "COL5A1", "COL5A2", "COL5A3",
  "COL6A1", "COL6A2", "COL6A3", "COL6A4", "COL6A5", "COL6A6",
  "COL7A1",
  "COL8A1", "COL8A2",
  "COL9A1", "COL9A2", "COL9A3",
  "COL10A1",
  "COL11A1", "COL11A2",
  "COL12A1",
  "COL13A1",
  "COL14A1",
  "COL15A1",
  "COL16A1",
  "COL17A1",
  "COL18A1",
  "COL19A1",
  "COL20A1",
  "COL21A1",
  "COL22A1",
  "COL23A1",
  "COL24A1",
  "COL25A1",
  "COL26A1",
  "COL27A1",
  "COL28A1"
)

# ============================================================
# 5. GET NORMALIZED RNA DATA
# ============================================================

get_rna_data <- function(seurat_obj) {
  
  data_input <- tryCatch(
    GetAssayData(
      seurat_obj,
      assay = "RNA",
      layer = "data"
    ),
    error = function(e) {
      GetAssayData(
        seurat_obj,
        assay = "RNA",
        slot = "data"
      )
    }
  )
  
  return(data_input)
}

expr_mat <- get_rna_data(seurat_obj)

type1_genes_present <- intersect(type1_genes, rownames(expr_mat))
type3_genes_present <- intersect(type3_genes, rownames(expr_mat))
total_collagen_genes_present <- intersect(total_collagen_genes, rownames(expr_mat))

cat("\nType I collagen genes used:\n")
print(type1_genes_present)

cat("\nType III collagen genes used:\n")
print(type3_genes_present)

cat("\nTotal collagen genes used:\n")
print(total_collagen_genes_present)

if (length(type1_genes_present) == 0) {
  stop("No Type I collagen genes were found.")
}

if (length(type3_genes_present) == 0) {
  stop("COL3A1 was not found.")
}

if (length(total_collagen_genes_present) == 0) {
  stop("No total collagen genes were found.")
}

# ============================================================
# 6. CALCULATE CELL-LEVEL EXPRESSION
# ============================================================

meta <- seurat_obj@meta.data
meta$Cell_ID <- rownames(meta)

meta <- meta[colnames(expr_mat), ]

type1_expr_per_cell <- Matrix::colSums(
  expr_mat[type1_genes_present, , drop = FALSE]
)

type3_expr_per_cell <- Matrix::colSums(
  expr_mat[type3_genes_present, , drop = FALSE]
)

total_collagen_expr_per_cell <- Matrix::colSums(
  expr_mat[total_collagen_genes_present, , drop = FALSE]
)

cell_level_df <- tibble(
  Cell_ID = colnames(expr_mat),
  orig.ident = meta$orig.ident,
  Sex = meta$Sex,
  Type_I_Collagen_Expression = as.numeric(type1_expr_per_cell),
  Type_III_Collagen_Expression = as.numeric(type3_expr_per_cell),
  Total_Collagen_Expression = as.numeric(total_collagen_expr_per_cell)
)

# ============================================================
# 7. RAW VALUES PER SAMPLE
# ============================================================

raw_sample_collagen_table <- cell_level_df %>%
  filter(!is.na(Sex)) %>%
  group_by(orig.ident, Sex) %>%
  summarise(
    Cells_Analyzed = n(),
    
    Type_I_Collagen_Expression = sum(
      Type_I_Collagen_Expression,
      na.rm = TRUE
    ),
    
    Type_III_Collagen_Expression = sum(
      Type_III_Collagen_Expression,
      na.rm = TRUE
    ),
    
    Total_Collagen_Expression = sum(
      Total_Collagen_Expression,
      na.rm = TRUE
    ),
    
    Type_III_Percent_of_Total_Collagen = ifelse(
      Total_Collagen_Expression > 0,
      100 * Type_III_Collagen_Expression / Total_Collagen_Expression,
      NA_real_
    ),
    
    Type_I_Percent_of_Total_Collagen = ifelse(
      Total_Collagen_Expression > 0,
      100 * Type_I_Collagen_Expression / Total_Collagen_Expression,
      NA_real_
    ),
    
    .groups = "drop"
  ) %>%
  arrange(orig.ident) %>%
  mutate(
    Type_I_Collagen_Expression = round(Type_I_Collagen_Expression, 4),
    Type_III_Collagen_Expression = round(Type_III_Collagen_Expression, 4),
    Total_Collagen_Expression = round(Total_Collagen_Expression, 4),
    Type_III_Percent_of_Total_Collagen = round(Type_III_Percent_of_Total_Collagen, 4),
    Type_I_Percent_of_Total_Collagen = round(Type_I_Percent_of_Total_Collagen, 4)
  )

cat("\n======================================================\n")
cat("RAW COLLAGEN VALUES PER SAMPLE\n")
cat("======================================================\n")
print(raw_sample_collagen_table, n = Inf)

write.csv(
  raw_sample_collagen_table,
  raw_sample_csv,
  row.names = FALSE
)

# ============================================================
# 8. T-TEST FUNCTION
# ============================================================

run_metric_ttest <- function(data, metric_col, metric_label, genes_used) {
  
  female_values <- data[[metric_col]][data$Sex == "Female"]
  male_values <- data[[metric_col]][data$Sex == "Male"]
  
  p_val <- tryCatch(
    {
      if (
        length(na.omit(female_values)) >= 2 &&
        length(na.omit(male_values)) >= 2
      ) {
        t.test(female_values, male_values)$p.value
      } else {
        NA_real_
      }
    },
    error = function(e) NA_real_
  )
  
  tibble(
    Metric = metric_label,
    Genes_Used = genes_used,
    Female_Mean = mean(female_values, na.rm = TRUE),
    Male_Mean = mean(male_values, na.rm = TRUE),
    Female_SD = sd(female_values, na.rm = TRUE),
    Male_SD = sd(male_values, na.rm = TRUE),
    Female_N = length(na.omit(female_values)),
    Male_N = length(na.omit(male_values)),
    p_value = p_val,
    Significance = case_when(
      is.na(p_val) ~ "NA",
      p_val < 0.001 ~ "***",
      p_val < 0.01 ~ "**",
      p_val < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
}

# ============================================================
# 9. P-VALUE SUMMARY TABLE
# ============================================================

collagen_pvalue_summary_table <- bind_rows(
  
  run_metric_ttest(
    data = raw_sample_collagen_table,
    metric_col = "Type_I_Collagen_Expression",
    metric_label = "Type I collagen expression",
    genes_used = paste(type1_genes_present, collapse = ", ")
  ),
  
  run_metric_ttest(
    data = raw_sample_collagen_table,
    metric_col = "Type_III_Collagen_Expression",
    metric_label = "Type III collagen expression",
    genes_used = paste(type3_genes_present, collapse = ", ")
  ),
  
  run_metric_ttest(
    data = raw_sample_collagen_table,
    metric_col = "Total_Collagen_Expression",
    metric_label = "Total collagen expression",
    genes_used = paste(total_collagen_genes_present, collapse = ", ")
  ),
  
  run_metric_ttest(
    data = raw_sample_collagen_table,
    metric_col = "Type_III_Percent_of_Total_Collagen",
    metric_label = "Type III collagen / total collagen (%)",
    genes_used = paste0(
      paste(type3_genes_present, collapse = ", "),
      " / total collagen"
    )
  ),
  
  run_metric_ttest(
    data = raw_sample_collagen_table,
    metric_col = "Type_I_Percent_of_Total_Collagen",
    metric_label = "Type I collagen / total collagen (%)",
    genes_used = paste0(
      paste(type1_genes_present, collapse = ", "),
      " / total collagen"
    )
  )
  
) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    Female_Mean = round(Female_Mean, 4),
    Male_Mean = round(Male_Mean, 4),
    Female_SD = round(Female_SD, 4),
    Male_SD = round(Male_SD, 4),
    p_value = signif(p_value, 4),
    p_adj_BH = signif(p_adj_BH, 4)
  ) %>%
  dplyr::select(
    Metric,
    Genes_Used,
    Female_Mean,
    Male_Mean,
    Female_SD,
    Male_SD,
    Female_N,
    Male_N,
    p_value,
    p_adj_BH,
    Significance
  )

cat("\n======================================================\n")
cat("P-VALUE SUMMARY TABLE\n")
cat("======================================================\n")
print(collagen_pvalue_summary_table, n = Inf)

write.csv(
  collagen_pvalue_summary_table,
  pvalue_csv,
  row.names = FALSE
)

# ============================================================
# 10. FINAL MESSAGE
# ============================================================

cat("\n======================================================\n")
cat("RAW SAMPLE COLLAGEN TABLE COMPLETE\n")
cat("Analysis used:\n")
cat(object_path, "\n\n")
cat("CSV files saved to Downloads:\n")
cat(raw_sample_csv, "\n")
cat(pvalue_csv, "\n")
cat("======================================================\n")



### PSC count

# ============================================================
# PSC PROPORTION IN FIBROBLAST SUBSET BY SEX
# Statistical unit = patient/sample
#
# PSCs = fibroblast_subset cluster 3
# Proportion = PSC cells / total fibroblast subset cells per sample
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

options(stringsAsFactors = FALSE)

# ============================================================
# 1. LOAD FIBROBLAST SUBSET
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

fibroblast_subset <- readRDS(
  file.path(base_dir, "Fibroblast_Subset.RDS")
)

DefaultAssay(fibroblast_subset) <- "RNA"

fibroblast_subset <- tryCatch(
  JoinLayers(fibroblast_subset),
  error = function(e) fibroblast_subset
)

# ============================================================
# 2. ADD SEX METADATA
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

fibroblast_subset$Sex <- unname(
  sex_mapping[as.character(fibroblast_subset$orig.ident)]
)

fibroblast_subset$Sex <- factor(
  fibroblast_subset$Sex,
  levels = c("Female", "Male")
)

# ============================================================
# 3. LABEL PSCs
# ============================================================
# Fibroblast subtype mapping:
# 0 = General Fibroblasts
# 1 = myCAFs
# 2 = iCAFs
# 3 = PSCs
# 4 = apCAFs

fibroblast_subset$CAF_type <- case_when(
  as.character(fibroblast_subset$seurat_clusters) == "0" ~ "General Fibroblasts",
  as.character(fibroblast_subset$seurat_clusters) == "1" ~ "myCAFs",
  as.character(fibroblast_subset$seurat_clusters) == "2" ~ "iCAFs",
  as.character(fibroblast_subset$seurat_clusters) == "3" ~ "PSCs",
  as.character(fibroblast_subset$seurat_clusters) == "4" ~ "apCAFs",
  TRUE ~ "Other"
)

fibroblast_subset$CAF_type <- factor(
  fibroblast_subset$CAF_type,
  levels = c(
    "General Fibroblasts",
    "myCAFs",
    "iCAFs",
    "PSCs",
    "apCAFs",
    "Other"
  )
)

cat("\nCell counts by CAF subtype and sex:\n")
print(table(fibroblast_subset$CAF_type, fibroblast_subset$Sex))

# ============================================================
# 4. RAW PSC COUNTS AND PROPORTIONS PER SAMPLE
# ============================================================

fibroblast_meta <- fibroblast_subset@meta.data

sample_totals <- fibroblast_meta %>%
  filter(!is.na(Sex)) %>%
  group_by(orig.ident, Sex) %>%
  summarise(
    Total_Fibroblast_Subset_Cells = n(),
    .groups = "drop"
  )

psc_counts <- fibroblast_meta %>%
  filter(!is.na(Sex), CAF_type == "PSCs") %>%
  group_by(orig.ident, Sex) %>%
  summarise(
    PSC_Cell_Count = n(),
    .groups = "drop"
  )

psc_sample_table <- sample_totals %>%
  left_join(
    psc_counts,
    by = c("orig.ident", "Sex")
  ) %>%
  mutate(
    PSC_Cell_Count = ifelse(
      is.na(PSC_Cell_Count),
      0,
      PSC_Cell_Count
    ),
    PSC_Proportion = PSC_Cell_Count / Total_Fibroblast_Subset_Cells,
    PSC_Percent = 100 * PSC_Proportion
  ) %>%
  arrange(orig.ident)

cat("\n======================================================\n")
cat("RAW PSC VALUES PER SAMPLE\n")
cat("======================================================\n")
print(psc_sample_table, n = Inf)

# ============================================================
# 5. MALE VS FEMALE SUMMARY: EXACT MEAN AND SD
# ============================================================

psc_sex_summary <- psc_sample_table %>%
  group_by(Sex) %>%
  summarise(
    Number_of_Samples = n(),
    
    Mean_PSC_Cell_Count = mean(PSC_Cell_Count, na.rm = TRUE),
    SD_PSC_Cell_Count = sd(PSC_Cell_Count, na.rm = TRUE),
    
    Mean_Total_Fibroblast_Cells = mean(Total_Fibroblast_Subset_Cells, na.rm = TRUE),
    SD_Total_Fibroblast_Cells = sd(Total_Fibroblast_Subset_Cells, na.rm = TRUE),
    
    Mean_PSC_Proportion = mean(PSC_Proportion, na.rm = TRUE),
    SD_PSC_Proportion = sd(PSC_Proportion, na.rm = TRUE),
    
    Mean_PSC_Percent = mean(PSC_Percent, na.rm = TRUE),
    SD_PSC_Percent = sd(PSC_Percent, na.rm = TRUE),
    
    .groups = "drop"
  )

cat("\n======================================================\n")
cat("MALE VS FEMALE PSC SUMMARY\n")
cat("======================================================\n")
print(psc_sex_summary, n = Inf)

# ============================================================
# 6. T-TESTS
# ============================================================

female_psc_percent <- psc_sample_table$PSC_Percent[
  psc_sample_table$Sex == "Female"
]

male_psc_percent <- psc_sample_table$PSC_Percent[
  psc_sample_table$Sex == "Male"
]

female_psc_count <- psc_sample_table$PSC_Cell_Count[
  psc_sample_table$Sex == "Female"
]

male_psc_count <- psc_sample_table$PSC_Cell_Count[
  psc_sample_table$Sex == "Male"
]

ttest_psc_percent <- t.test(
  female_psc_percent,
  male_psc_percent
)

ttest_psc_count <- t.test(
  female_psc_count,
  male_psc_count
)

psc_ttest_table <- tibble(
  Metric = c(
    "PSC percent of fibroblast subset",
    "PSC cell count"
  ),
  Female_Mean = c(
    mean(female_psc_percent, na.rm = TRUE),
    mean(female_psc_count, na.rm = TRUE)
  ),
  Male_Mean = c(
    mean(male_psc_percent, na.rm = TRUE),
    mean(male_psc_count, na.rm = TRUE)
  ),
  Female_SD = c(
    sd(female_psc_percent, na.rm = TRUE),
    sd(female_psc_count, na.rm = TRUE)
  ),
  Male_SD = c(
    sd(male_psc_percent, na.rm = TRUE),
    sd(male_psc_count, na.rm = TRUE)
  ),
  Female_N = c(
    length(na.omit(female_psc_percent)),
    length(na.omit(female_psc_count))
  ),
  Male_N = c(
    length(na.omit(male_psc_percent)),
    length(na.omit(male_psc_count))
  ),
  p_value = c(
    ttest_psc_percent$p.value,
    ttest_psc_count$p.value
  )
) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

cat("\n======================================================\n")
cat("PSC FEMALE VS MALE T-TEST RESULTS\n")
cat("======================================================\n")
print(psc_ttest_table, n = Inf)

cat("\nRaw t-test for PSC percent:\n")
print(ttest_psc_percent)

cat("\nRaw t-test for PSC cell count:\n")
print(ttest_psc_count)

# ============================================================
# 7. CLEAN BARPLOT: PSC PROPORTION BY SEX
# ============================================================

p_label <- paste0(
  "p = ",
  signif(ttest_psc_percent$p.value, 2)
)

plot_summary <- psc_sample_table %>%
  group_by(Sex) %>%
  summarise(
    Mean_PSC_Percent = mean(PSC_Percent, na.rm = TRUE),
    SEM_PSC_Percent = sd(PSC_Percent, na.rm = TRUE) /
      sqrt(sum(!is.na(PSC_Percent))),
    .groups = "drop"
  )

y_position <- max(
  plot_summary$Mean_PSC_Percent + plot_summary$SEM_PSC_Percent,
  na.rm = TRUE
) * 1.25

if (!is.finite(y_position) || y_position == 0) {
  y_position <- 1
}

psc_percent_plot <- ggplot(
  plot_summary,
  aes(
    x = Sex,
    y = Mean_PSC_Percent,
    fill = Sex
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_errorbar(
    aes(
      ymin = Mean_PSC_Percent - SEM_PSC_Percent,
      ymax = Mean_PSC_Percent + SEM_PSC_Percent
    ),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      x = 1.5,
      y = y_position,
      label = p_label
    ),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.28))
  ) +
  labs(
    title = "PSC Proportion Within Fibroblast Subset by Sex",
    x = "Sex",
    y = "Mean PSC percent of fibroblast subset per sample"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 17
    ),
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "none"
  )

print(psc_percent_plot)

# ============================================================
# 8. OPTIONAL: SAVE TABLES TO DOWNLOADS
# ============================================================

write.csv(
  psc_sample_table,
  file.path(path.expand("~/Downloads"), "PSC_Raw_Values_Per_Sample.csv"),
  row.names = FALSE
)

write.csv(
  psc_sex_summary,
  file.path(path.expand("~/Downloads"), "PSC_Sex_Summary_Mean_SD.csv"),
  row.names = FALSE
)

write.csv(
  psc_ttest_table,
  file.path(path.expand("~/Downloads"), "PSC_Ttest_By_Sex.csv"),
  row.names = FALSE
)

cat("\n======================================================\n")
cat("PSC PROPORTION ANALYSIS COMPLETE\n")
cat("CSV files saved to Downloads:\n")
cat("- PSC_Raw_Values_Per_Sample.csv\n")
cat("- PSC_Sex_Summary_Mean_SD.csv\n")
cat("- PSC_Ttest_By_Sex.csv\n")
cat("======================================================\n")


### Finding Log2 Fold Change Values

# ============================================================
# LOG2 FOLD CHANGE VALUES:
# Fibroblasts: COL3A1, COL5A2, COL12A1
# Schwann cells: INPP4B
#
# Comparison:
# Male vs Female
#
# Interpretation:
# avg_log2FC > 0 = higher in Male
# avg_log2FC < 0 = higher in Female
# ============================================================

library(Seurat)
library(dplyr)
library(tibble)

options(stringsAsFactors = FALSE)

# ============================================================
# 1. SET PATHS
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

fibroblast_path <- file.path(
  base_dir,
  "Fibroblast_Subset.RDS"
)

schwann_path <- file.path(
  base_dir,
  "true_schwann_ERBB3_PMP22_subset.RDS"
)

# If your Schwann file is saved with lowercase .rds
if (!file.exists(schwann_path)) {
  schwann_path <- file.path(
    base_dir,
    "true_schwann_ERBB3_PMP22_subset.rds"
  )
}

if (!file.exists(fibroblast_path)) {
  stop("Fibroblast_Subset.RDS was not found.")
}

if (!file.exists(schwann_path)) {
  stop("true_schwann_ERBB3_PMP22_subset.RDS was not found.")
}

# ============================================================
# 2. SEX MAPPING
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

# ============================================================
# 3. PREPARE SEURAT OBJECT
# ============================================================

prepare_object <- function(seurat_obj) {
  
  DefaultAssay(seurat_obj) <- "RNA"
  
  seurat_obj <- tryCatch(
    JoinLayers(seurat_obj),
    error = function(e) seurat_obj
  )
  
  seurat_obj$Sex <- unname(
    sex_mapping[as.character(seurat_obj$orig.ident)]
  )
  
  seurat_obj <- subset(
    seurat_obj,
    subset = Sex %in% c("Male", "Female")
  )
  
  seurat_obj$Sex <- factor(
    seurat_obj$Sex,
    levels = c("Female", "Male")
  )
  
  return(seurat_obj)
}

# ============================================================
# 4. FUNCTION TO RUN GENE-SPECIFIC LOG2FC
# ============================================================

run_gene_log2fc <- function(seurat_obj, genes, cell_type_label) {
  
  genes_present <- intersect(
    genes,
    rownames(seurat_obj)
  )
  
  genes_missing <- setdiff(
    genes,
    rownames(seurat_obj)
  )
  
  cat("\n======================================================\n")
  cat(cell_type_label, "\n")
  cat("Genes present:\n")
  print(genes_present)
  cat("Genes missing:\n")
  print(genes_missing)
  cat("======================================================\n")
  
  if (length(genes_present) == 0) {
    stop(
      paste(
        "None of the requested genes were found in",
        cell_type_label
      )
    )
  }
  
  markers <- FindMarkers(
    object = seurat_obj,
    ident.1 = "Male",
    ident.2 = "Female",
    group.by = "Sex",
    features = genes_present,
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = FALSE
  )
  
  if (nrow(markers) == 0) {
    warning(
      paste(
        "FindMarkers returned no rows for",
        cell_type_label
      )
    )
    
    return(
      tibble(
        Cell_Type = cell_type_label,
        Gene = genes_present,
        Comparison = "Male vs Female",
        avg_log2FC = NA_real_,
        Direction = NA_character_,
        p_val = NA_real_,
        p_val_adj = NA_real_,
        pct.1 = NA_real_,
        pct.2 = NA_real_,
        pct.1_label = "Percent expressed in Male",
        pct.2_label = "Percent expressed in Female"
      )
    )
  }
  
  markers$Gene <- rownames(markers)
  
  result <- markers %>%
    as_tibble() %>%
    dplyr::select(
      Gene,
      avg_log2FC,
      p_val,
      p_val_adj,
      pct.1,
      pct.2
    ) %>%
    mutate(
      Cell_Type = cell_type_label,
      Comparison = "Male vs Female",
      Direction = case_when(
        avg_log2FC > 0 ~ "Higher in Male",
        avg_log2FC < 0 ~ "Higher in Female",
        TRUE ~ "No difference"
      ),
      pct.1_label = "Percent expressed in Male",
      pct.2_label = "Percent expressed in Female"
    ) %>%
    dplyr::select(
      Cell_Type,
      Gene,
      Comparison,
      avg_log2FC,
      Direction,
      p_val,
      p_val_adj,
      pct.1,
      pct.2,
      pct.1_label,
      pct.2_label
    )
  
  return(result)
}

# ============================================================
# 5. FIBROBLAST LOG2FC VALUES
# ============================================================

fibroblast_subset <- readRDS(fibroblast_path)

fibroblast_subset <- prepare_object(
  fibroblast_subset
)

cat("\n======================================================\n")
cat("FIBROBLAST CELLS BY SEX\n")
cat("======================================================\n")
print(table(fibroblast_subset$Sex))

cat("\n======================================================\n")
cat("FIBROBLAST CELLS BY SAMPLE AND SEX\n")
cat("======================================================\n")
print(table(fibroblast_subset$orig.ident, fibroblast_subset$Sex))

fibroblast_genes <- c(
  "COL3A1",
  "COL5A2",
  "COL12A1"
)

fibroblast_log2fc <- run_gene_log2fc(
  seurat_obj = fibroblast_subset,
  genes = fibroblast_genes,
  cell_type_label = "Fibroblasts"
)

# ============================================================
# 6. SCHWANN LOG2FC VALUE
# ============================================================

true_schwann_subset <- readRDS(schwann_path)

true_schwann_subset <- prepare_object(
  true_schwann_subset
)

cat("\n======================================================\n")
cat("SCHWANN CELLS BY SEX\n")
cat("======================================================\n")
print(table(true_schwann_subset$Sex))

cat("\n======================================================\n")
cat("SCHWANN CELLS BY SAMPLE AND SEX\n")
cat("======================================================\n")
print(table(true_schwann_subset$orig.ident, true_schwann_subset$Sex))

schwann_genes <- c(
  "INPP4B"
)

schwann_log2fc <- run_gene_log2fc(
  seurat_obj = true_schwann_subset,
  genes = schwann_genes,
  cell_type_label = "Schwann cells"
)

# ============================================================
# 7. COMBINE RESULTS
# ============================================================

target_gene_log2fc_results <- dplyr::bind_rows(
  fibroblast_log2fc,
  schwann_log2fc
) %>%
  arrange(Cell_Type, Gene) %>%
  as_tibble()

cat("\n======================================================\n")
cat("REQUESTED GENE LOG2 FOLD CHANGE RESULTS\n")
cat("======================================================\n")

# Safe print version
print(
  as.data.frame(target_gene_log2fc_results),
  row.names = FALSE
)

# ============================================================
# 8. PRINT ONLY KEY VALUES
# ============================================================

key_results <- target_gene_log2fc_results %>%
  dplyr::select(
    Cell_Type,
    Gene,
    avg_log2FC,
    Direction,
    p_val,
    p_val_adj,
    pct.1,
    pct.2
  )

cat("\n======================================================\n")
cat("KEY LOG2FC VALUES ONLY\n")
cat("======================================================\n")

print(
  as.data.frame(key_results),
  row.names = FALSE
)

# ============================================================
# 9. SAVE TO DOWNLOADS
# ============================================================

output_file <- file.path(
  path.expand("~/Downloads"),
  "Requested_Gene_Log2FC_Fibroblast_Schwann_Male_vs_Female.csv"
)

write.csv(
  as.data.frame(target_gene_log2fc_results),
  output_file,
  row.names = FALSE
)

key_output_file <- file.path(
  path.expand("~/Downloads"),
  "Requested_Gene_Log2FC_Key_Values.csv"
)

write.csv(
  as.data.frame(key_results),
  key_output_file,
  row.names = FALSE
)

cat("\n======================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("CSV files saved to Downloads:\n")
cat(output_file, "\n")
cat(key_output_file, "\n")
cat("======================================================\n")

cat("\nINTERPRETATION:\n")
cat("avg_log2FC > 0 means higher in Male.\n")
cat("avg_log2FC < 0 means higher in Female.\n")
cat("pct.1 = percent expressed in Male.\n")
cat("pct.2 = percent expressed in Female.\n")
cat("======================================================\n")



### Patient Level Summary Table

# ============================================================
# PATIENT-LEVEL SUMMARY CSV - FIXED VERSION
#
# Includes only:
# - Patient ID
# - Sex
# - Total nuclei
# - Total fibroblast nuclei
# - Total Schwann nuclei
# - Percentage Schwann of total nuclei
# - Percentage iCAF of fibroblast nuclei
# - Percentage myCAF of fibroblast nuclei
# - Percentage apCAF of fibroblast nuclei
# - PSC-like proportion of fibroblast nuclei
#
# Output:
# ~/Downloads/Patient_Level_Simplified_Summary.csv
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)

options(stringsAsFactors = FALSE)

# ============================================================
# 1. PATHS
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

merged_path <- file.path(
  base_dir,
  "Merged_Joined.RDS"
)

fibroblast_path <- file.path(
  base_dir,
  "Fibroblast_Subset.RDS"
)

schwann_path <- file.path(
  base_dir,
  "true_schwann_ERBB3_PMP22_subset.RDS"
)

if (!file.exists(schwann_path)) {
  schwann_path <- file.path(
    base_dir,
    "true_schwann_ERBB3_PMP22_subset.rds"
  )
}

output_file <- file.path(
  path.expand("~/Downloads"),
  "Patient_Level_Simplified_Summary.csv"
)

# ============================================================
# 2. CHECK FILES EXIST
# ============================================================

if (!file.exists(merged_path)) {
  stop("Merged_Joined.RDS was not found.")
}

if (!file.exists(fibroblast_path)) {
  stop("Fibroblast_Subset.RDS was not found.")
}

if (!file.exists(schwann_path)) {
  stop("true_schwann_ERBB3_PMP22_subset.RDS was not found.")
}

# ============================================================
# 3. SEX AND PATIENT MAPPING
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

patient_mapping <- c(
  "Sample_01" = "P01",
  "Sample_02" = "P02",
  "Sample_03" = "P03",
  "Sample_04" = "P04",
  "Sample_05" = "P05",
  "Sample_06" = "P06",
  "Sample_07" = "P07",
  "Sample_08" = "P08",
  "Sample_09" = "P09",
  "Sample_10" = "P10",
  "Sample_11" = "P11",
  "Sample_12" = "P12",
  "Sample_13" = "P13",
  "Sample_14" = "P14",
  "Sample_15" = "P15",
  "Sample_16" = "P16",
  "Sample_17" = "P17"
)

sample_base <- tibble::tibble(
  orig.ident = names(sex_mapping),
  Patient_ID = unname(patient_mapping[names(sex_mapping)]),
  Sex = unname(sex_mapping[names(sex_mapping)])
) %>%
  dplyr::mutate(
    Patient_ID = factor(Patient_ID, levels = sprintf("P%02d", 1:17)),
    Sex = factor(Sex, levels = c("Female", "Male"))
  ) %>%
  dplyr::arrange(Patient_ID)

# ============================================================
# 4. LOAD OBJECTS
# ============================================================

Merged_Joined <- readRDS(merged_path)
fibroblast_subset <- readRDS(fibroblast_path)
true_schwann_subset <- readRDS(schwann_path)

DefaultAssay(Merged_Joined) <- "RNA"
DefaultAssay(fibroblast_subset) <- "RNA"
DefaultAssay(true_schwann_subset) <- "RNA"

Merged_Joined <- tryCatch(
  JoinLayers(Merged_Joined),
  error = function(e) Merged_Joined
)

fibroblast_subset <- tryCatch(
  JoinLayers(fibroblast_subset),
  error = function(e) fibroblast_subset
)

true_schwann_subset <- tryCatch(
  JoinLayers(true_schwann_subset),
  error = function(e) true_schwann_subset
)

# ============================================================
# 5. CREATE METADATA TABLES
# ============================================================

merged_meta <- Merged_Joined@meta.data %>%
  tibble::rownames_to_column("Cell_ID") %>%
  dplyr::mutate(
    orig.ident = as.character(orig.ident)
  ) %>%
  dplyr::filter(
    orig.ident %in% names(sex_mapping)
  )

fibroblast_meta <- fibroblast_subset@meta.data %>%
  tibble::rownames_to_column("Cell_ID") %>%
  dplyr::mutate(
    orig.ident = as.character(orig.ident)
  ) %>%
  dplyr::filter(
    orig.ident %in% names(sex_mapping)
  )

schwann_meta <- true_schwann_subset@meta.data %>%
  tibble::rownames_to_column("Cell_ID") %>%
  dplyr::mutate(
    orig.ident = as.character(orig.ident)
  ) %>%
  dplyr::filter(
    orig.ident %in% names(sex_mapping)
  )

# ============================================================
# 6. TOTAL NUCLEI PER PATIENT
# ============================================================

total_nuclei_table <- merged_meta %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(
    Total_Nuclei = dplyr::n(),
    .groups = "drop"
  )

# ============================================================
# 7. TOTAL FIBROBLAST NUCLEI PER PATIENT
# ============================================================

fibroblast_nuclei_table <- fibroblast_meta %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(
    Total_Fibroblast_Nuclei = dplyr::n(),
    .groups = "drop"
  )

# ============================================================
# 8. TOTAL SCHWANN NUCLEI PER PATIENT
# ============================================================

schwann_nuclei_table <- schwann_meta %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(
    Total_Schwann_Nuclei = dplyr::n(),
    .groups = "drop"
  )

# ============================================================
# 9. FIBROBLAST SUBTYPE PERCENTAGES
# ============================================================
# Fibroblast subset cluster mapping:
# 1 = myCAF
# 2 = iCAF
# 3 = PSC-like
# 4 = apCAF

fibroblast_meta <- fibroblast_meta %>%
  dplyr::mutate(
    Fibroblast_Subtype = dplyr::case_when(
      as.character(seurat_clusters) == "1" ~ "myCAF",
      as.character(seurat_clusters) == "2" ~ "iCAF",
      as.character(seurat_clusters) == "3" ~ "PSC_like",
      as.character(seurat_clusters) == "4" ~ "apCAF",
      TRUE ~ "Other"
    )
  )

subtype_order <- c(
  "myCAF",
  "iCAF",
  "apCAF",
  "PSC_like"
)

fibroblast_subtype_counts <- fibroblast_meta %>%
  dplyr::filter(
    Fibroblast_Subtype %in% subtype_order
  ) %>%
  dplyr::group_by(
    orig.ident,
    Fibroblast_Subtype
  ) %>%
  dplyr::summarise(
    Subtype_Nuclei = dplyr::n(),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    orig.ident = names(sex_mapping),
    Fibroblast_Subtype = subtype_order,
    fill = list(Subtype_Nuclei = 0)
  ) %>%
  dplyr::left_join(
    fibroblast_nuclei_table,
    by = "orig.ident"
  ) %>%
  dplyr::mutate(
    Total_Fibroblast_Nuclei = ifelse(
      is.na(Total_Fibroblast_Nuclei),
      0,
      Total_Fibroblast_Nuclei
    ),
    Subtype_Percent_of_Fibroblast_Nuclei = ifelse(
      Total_Fibroblast_Nuclei > 0,
      100 * Subtype_Nuclei / Total_Fibroblast_Nuclei,
      NA_real_
    ),
    Subtype_Proportion_of_Fibroblast_Nuclei = ifelse(
      Total_Fibroblast_Nuclei > 0,
      Subtype_Nuclei / Total_Fibroblast_Nuclei,
      NA_real_
    )
  )

# ============================================================
# 10. WIDE TABLE FOR iCAF, myCAF, apCAF PERCENTAGES
# ============================================================

fibroblast_subtype_percent_wide <- fibroblast_subtype_counts %>%
  dplyr::filter(
    Fibroblast_Subtype %in% c("iCAF", "myCAF", "apCAF")
  ) %>%
  dplyr::select(
    orig.ident,
    Fibroblast_Subtype,
    Subtype_Percent_of_Fibroblast_Nuclei
  ) %>%
  dplyr::mutate(
    Column_Name = paste0(
      Fibroblast_Subtype,
      "_Percent_of_Fibroblast_Nuclei"
    )
  ) %>%
  dplyr::select(
    orig.ident,
    Column_Name,
    Subtype_Percent_of_Fibroblast_Nuclei
  ) %>%
  tidyr::pivot_wider(
    names_from = Column_Name,
    values_from = Subtype_Percent_of_Fibroblast_Nuclei
  )

# ============================================================
# 11. PSC-LIKE PROPORTION TABLE
# ============================================================

psc_proportion_table <- fibroblast_subtype_counts %>%
  dplyr::filter(
    Fibroblast_Subtype == "PSC_like"
  ) %>%
  dplyr::select(
    orig.ident,
    PSC_like_Proportion = Subtype_Proportion_of_Fibroblast_Nuclei
  )

# ============================================================
# 12. COMBINE FINAL PATIENT-LEVEL TABLE
# ============================================================

patient_level_summary <- sample_base %>%
  dplyr::left_join(
    total_nuclei_table,
    by = "orig.ident"
  ) %>%
  dplyr::left_join(
    fibroblast_nuclei_table,
    by = "orig.ident"
  ) %>%
  dplyr::left_join(
    schwann_nuclei_table,
    by = "orig.ident"
  ) %>%
  dplyr::mutate(
    Total_Nuclei = ifelse(
      is.na(Total_Nuclei),
      0,
      Total_Nuclei
    ),
    Total_Fibroblast_Nuclei = ifelse(
      is.na(Total_Fibroblast_Nuclei),
      0,
      Total_Fibroblast_Nuclei
    ),
    Total_Schwann_Nuclei = ifelse(
      is.na(Total_Schwann_Nuclei),
      0,
      Total_Schwann_Nuclei
    ),
    Percent_Schwann_of_Total_Nuclei = ifelse(
      Total_Nuclei > 0,
      100 * Total_Schwann_Nuclei / Total_Nuclei,
      NA_real_
    )
  ) %>%
  dplyr::left_join(
    fibroblast_subtype_percent_wide,
    by = "orig.ident"
  ) %>%
  dplyr::left_join(
    psc_proportion_table,
    by = "orig.ident"
  ) %>%
  dplyr::mutate(
    iCAF_Percent_of_Fibroblast_Nuclei = ifelse(
      is.na(iCAF_Percent_of_Fibroblast_Nuclei),
      0,
      iCAF_Percent_of_Fibroblast_Nuclei
    ),
    myCAF_Percent_of_Fibroblast_Nuclei = ifelse(
      is.na(myCAF_Percent_of_Fibroblast_Nuclei),
      0,
      myCAF_Percent_of_Fibroblast_Nuclei
    ),
    apCAF_Percent_of_Fibroblast_Nuclei = ifelse(
      is.na(apCAF_Percent_of_Fibroblast_Nuclei),
      0,
      apCAF_Percent_of_Fibroblast_Nuclei
    ),
    PSC_like_Proportion = ifelse(
      is.na(PSC_like_Proportion),
      0,
      PSC_like_Proportion
    )
  ) %>%
  dplyr::select(
    Patient_ID,
    orig.ident,
    Sex,
    Total_Nuclei,
    Total_Fibroblast_Nuclei,
    Total_Schwann_Nuclei,
    Percent_Schwann_of_Total_Nuclei,
    iCAF_Percent_of_Fibroblast_Nuclei,
    myCAF_Percent_of_Fibroblast_Nuclei,
    apCAF_Percent_of_Fibroblast_Nuclei,
    PSC_like_Proportion
  ) %>%
  dplyr::arrange(Patient_ID)

# ============================================================
# 13. ROUND NUMERIC VALUES
# ============================================================

patient_level_summary <- patient_level_summary %>%
  dplyr::mutate(
    Percent_Schwann_of_Total_Nuclei = round(
      Percent_Schwann_of_Total_Nuclei,
      4
    ),
    iCAF_Percent_of_Fibroblast_Nuclei = round(
      iCAF_Percent_of_Fibroblast_Nuclei,
      4
    ),
    myCAF_Percent_of_Fibroblast_Nuclei = round(
      myCAF_Percent_of_Fibroblast_Nuclei,
      4
    ),
    apCAF_Percent_of_Fibroblast_Nuclei = round(
      apCAF_Percent_of_Fibroblast_Nuclei,
      4
    ),
    PSC_like_Proportion = round(
      PSC_like_Proportion,
      6
    )
  )

# ============================================================
# 14. PRINT FINAL TABLE
# ============================================================

cat("\n======================================================\n")
cat("PATIENT-LEVEL SIMPLIFIED SUMMARY TABLE\n")
cat("======================================================\n")

print(
  as.data.frame(patient_level_summary),
  row.names = FALSE
)

# ============================================================
# 15. SAVE CSV
# ============================================================

write.csv(
  as.data.frame(patient_level_summary),
  output_file,
  row.names = FALSE
)

cat("\n======================================================\n")
cat("CSV SAVED TO DOWNLOADS:\n")
cat(output_file, "\n")
cat("======================================================\n")


### Patient-Level Fibroblast DGE

# ============================================================
# PATIENT-LEVEL PSEUDOBULK DIFFERENTIAL EXPRESSION
# Male versus Female PDAC Fibroblasts
#
# Statistical unit: PATIENT, not nucleus
#
# This version:
# - does NOT install/update packages
# - does NOT save files
# - prints tables in console
# - prints volcano plot in RStudio
#
# Positive log2FC = higher in male patients
# Negative log2FC = higher in female patients
# ============================================================

# ============================================================
# 0. LOAD PACKAGES WITHOUT INSTALLING/UPDATING
# ============================================================

required_packages <- c(
  "Seurat",
  "Matrix",
  "edgeR",
  "dplyr",
  "tibble",
  "ggplot2",
  "ggrepel"
)

missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "These packages are missing: ",
      paste(missing_packages, collapse = ", "),
      "\nInstall them manually first, then rerun this script."
    )
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
})

options(stringsAsFactors = FALSE)

# ============================================================
# 1. LOAD FIBROBLAST SEURAT OBJECT
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

fibroblast_path <- file.path(
  base_dir,
  "Fibroblast_Subset.RDS"
)

if (!file.exists(fibroblast_path)) {
  stop("Fibroblast_Subset.RDS was not found.")
}

fibroblast_subset <- readRDS(fibroblast_path)

DefaultAssay(fibroblast_subset) <- "RNA"

fibroblast_subset <- tryCatch(
  JoinLayers(
    fibroblast_subset,
    assay = "RNA"
  ),
  error = function(e) {
    message(
      "JoinLayers was not required or could not be applied: ",
      e$message
    )
    fibroblast_subset
  }
)

# ============================================================
# 2. ADD PATIENT-LEVEL SEX INFORMATION
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

fibroblast_subset$Sex <- unname(
  sex_mapping[as.character(fibroblast_subset$orig.ident)]
)

fibroblast_subset <- subset(
  fibroblast_subset,
  subset = !is.na(orig.ident) &
    Sex %in% c("Female", "Male")
)

fibroblast_subset$orig.ident <- as.character(
  fibroblast_subset$orig.ident
)

fibroblast_subset$Sex <- factor(
  fibroblast_subset$Sex,
  levels = c("Female", "Male")
)

cat("\nFibroblast nuclei by sex:\n")
print(table(fibroblast_subset$Sex))

cat("\nFibroblast nuclei by patient:\n")
print(table(fibroblast_subset$orig.ident))

# ============================================================
# 3. EXTRACT RAW RNA COUNTS
# ============================================================
# Pseudobulk DE must use raw counts.

raw_counts <- tryCatch(
  LayerData(
    object = fibroblast_subset,
    assay = "RNA",
    layer = "counts"
  ),
  error = function(e) {
    GetAssayData(
      object = fibroblast_subset,
      assay = "RNA",
      slot = "counts"
    )
  }
)

raw_counts <- as(
  raw_counts,
  "dgCMatrix"
)

fibroblast_metadata <- fibroblast_subset@meta.data[
  colnames(raw_counts),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(fibroblast_metadata),
    colnames(raw_counts)
  )
)

# ============================================================
# 4. VERIFY EACH PATIENT HAS ONE SEX LABEL
# ============================================================

patient_sex_check <- fibroblast_metadata %>%
  tibble::rownames_to_column("Cell_ID") %>%
  dplyr::transmute(
    Patient = as.character(orig.ident),
    Sex = as.character(Sex)
  ) %>%
  dplyr::group_by(Patient) %>%
  dplyr::summarise(
    Number_of_Sex_Labels = dplyr::n_distinct(Sex),
    Sex_Label = paste(unique(Sex), collapse = ";"),
    .groups = "drop"
  )

cat("\nPatient sex label check:\n")
print(
  as.data.frame(patient_sex_check),
  row.names = FALSE
)

if (any(patient_sex_check$Number_of_Sex_Labels != 1)) {
  stop(
    "At least one patient has more than one sex label. Check orig.ident and Sex metadata."
  )
}

# ============================================================
# 5. CREATE PATIENT-LEVEL PSEUDOBULK COUNTS
# ============================================================
# Original matrix: genes x nuclei
# New matrix: genes x patients

patient_factor <- factor(
  as.character(fibroblast_metadata$orig.ident),
  levels = sort(
    unique(
      as.character(fibroblast_metadata$orig.ident)
    )
  )
)

aggregation_matrix <- Matrix::sparse.model.matrix(
  ~ 0 + patient_factor
)

colnames(aggregation_matrix) <- levels(patient_factor)

pseudobulk_counts <- raw_counts %*% aggregation_matrix

colnames(pseudobulk_counts) <- levels(patient_factor)

cat("\nPseudobulk count matrix dimensions:\n")
cat(
  nrow(pseudobulk_counts),
  "genes x",
  ncol(pseudobulk_counts),
  "patients\n"
)

cat("\nPseudobulk count columns:\n")
print(colnames(pseudobulk_counts))

# ============================================================
# 6. CREATE ONE METADATA ROW PER PATIENT - FIXED ALIGNMENT
# ============================================================

patient_metadata <- fibroblast_metadata %>%
  tibble::rownames_to_column("Cell_ID") %>%
  dplyr::transmute(
    Patient = as.character(orig.ident),
    Sex = as.character(Sex)
  ) %>%
  dplyr::group_by(Patient, Sex) %>%
  dplyr::summarise(
    Fibroblast_Nuclei = dplyr::n(),
    .groups = "drop"
  ) %>%
  as.data.frame()

cat("\nPatient metadata before matching:\n")
print(patient_metadata)

# Check that every pseudobulk column has matching metadata
missing_metadata <- setdiff(
  colnames(pseudobulk_counts),
  patient_metadata$Patient
)

if (length(missing_metadata) > 0) {
  stop(
    paste0(
      "These patients are in pseudobulk_counts but missing from patient_metadata: ",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

# Reorder metadata to exactly match pseudobulk count columns
patient_metadata <- patient_metadata[
  match(
    colnames(pseudobulk_counts),
    patient_metadata$Patient
  ),
  ,
  drop = FALSE
]

# Force rownames to match pseudobulk count columns
rownames(patient_metadata) <- patient_metadata$Patient

patient_metadata$Sex <- factor(
  patient_metadata$Sex,
  levels = c("Female", "Male")
)

cat("\nChecking patient metadata alignment before filtering:\n")

cat("\npseudobulk count columns:\n")
print(colnames(pseudobulk_counts))

cat("\npatient metadata rownames:\n")
print(rownames(patient_metadata))

if (!all(rownames(patient_metadata) == colnames(pseudobulk_counts))) {
  stop(
    "Patient metadata rownames do not match pseudobulk count columns before filtering."
  )
}

cat("\nPatient-level pseudobulk metadata:\n")
print(patient_metadata)

cat("\nNumber of patients by sex before filtering:\n")
print(table(patient_metadata$Sex))

# ============================================================
# 7. FILTER PATIENTS WITH TOO FEW FIBROBLAST NUCLEI
# ============================================================

min_fibroblast_nuclei <- 50L

eligible_patients <- patient_metadata$Patient[
  patient_metadata$Fibroblast_Nuclei >= min_fibroblast_nuclei
]

excluded_patients <- setdiff(
  patient_metadata$Patient,
  eligible_patients
)

if (length(excluded_patients) > 0) {
  message(
    "Excluding patients with fewer than ",
    min_fibroblast_nuclei,
    " fibroblast nuclei: ",
    paste(excluded_patients, collapse = ", ")
  )
}

# Subset pseudobulk matrix
pseudobulk_counts <- pseudobulk_counts[
  ,
  eligible_patients,
  drop = FALSE
]

# Subset metadata using Patient column
patient_metadata <- patient_metadata[
  patient_metadata$Patient %in% eligible_patients,
  ,
  drop = FALSE
]

# Reorder metadata again after filtering
patient_metadata <- patient_metadata[
  match(
    colnames(pseudobulk_counts),
    patient_metadata$Patient
  ),
  ,
  drop = FALSE
]

rownames(patient_metadata) <- patient_metadata$Patient

patient_metadata$Sex <- factor(
  patient_metadata$Sex,
  levels = c("Female", "Male")
)

cat("\nChecking patient metadata alignment after filtering:\n")

cat("\npseudobulk count columns:\n")
print(colnames(pseudobulk_counts))

cat("\npatient metadata rownames:\n")
print(rownames(patient_metadata))

if (!all(rownames(patient_metadata) == colnames(pseudobulk_counts))) {
  stop(
    "After filtering, patient metadata rownames do not match pseudobulk count columns."
  )
}

cat("\nPatients retained by sex after filtering:\n")
print(table(patient_metadata$Sex))

if (any(table(patient_metadata$Sex) < 2)) {
  stop(
    "Too few patients remain in one sex group for patient-level differential expression."
  )
}

# ============================================================
# 8. CREATE edgeR OBJECT
# ============================================================

dge <- edgeR::DGEList(
  counts = pseudobulk_counts,
  samples = patient_metadata
)

# Female is the reference group.
# SexMale = Male - Female.
design <- model.matrix(
  ~ Sex,
  data = patient_metadata
)

rownames(design) <- rownames(patient_metadata)

cat("\nDesign matrix:\n")
print(design)

if (!"SexMale" %in% colnames(design)) {
  stop("SexMale coefficient was not found in the design matrix.")
}

# ============================================================
# 9. FILTER LOWLY EXPRESSED GENES
# ============================================================

keep_genes <- edgeR::filterByExpr(
  dge,
  design = design
)

cat("\nGenes before expression filtering:", nrow(dge), "\n")
cat("Genes retained after expression filtering:", sum(keep_genes), "\n")

dge <- dge[
  keep_genes,
  ,
  keep.lib.sizes = FALSE
]

# ============================================================
# 10. NORMALIZE PATIENT-LEVEL LIBRARY SIZES
# ============================================================

dge <- edgeR::calcNormFactors(
  dge,
  method = "TMM"
)

cat("\nedgeR patient-level library information:\n")
print(dge$samples)

# ============================================================
# 11. FIT PATIENT-LEVEL DIFFERENTIAL EXPRESSION MODEL
# ============================================================

dge <- edgeR::estimateDisp(
  dge,
  design,
  robust = TRUE
)

fit <- edgeR::glmQLFit(
  dge,
  design,
  robust = TRUE
)

qlf_test <- edgeR::glmQLFTest(
  fit,
  coef = "SexMale"
)

# ============================================================
# 12. EXTRACT ALL RESULTS
# ============================================================

fibroblast_patient_level_DE <- edgeR::topTags(
  qlf_test,
  n = Inf,
  sort.by = "PValue"
)$table %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::rename(
    log2FC_Male_vs_Female = logFC,
    Average_LogCPM = logCPM,
    Raw_P_Value = PValue,
    FDR = FDR
  ) %>%
  dplyr::mutate(
    Direction = dplyr::case_when(
      FDR < 0.05 &
        log2FC_Male_vs_Female >= 0.5 ~
        "Male-upregulated",
      
      FDR < 0.05 &
        log2FC_Male_vs_Female <= -0.5 ~
        "Female-upregulated",
      
      TRUE ~ "Not significant"
    ),
    GSEA_Rank = sign(log2FC_Male_vs_Female) * sqrt(F)
  ) %>%
  dplyr::arrange(
    FDR,
    dplyr::desc(abs(log2FC_Male_vs_Female))
  )

# ============================================================
# 13. PRINT RESULT SUMMARY
# ============================================================

cat("\nPatient-level differential-expression summary:\n")
print(table(fibroblast_patient_level_DE$Direction))

cat("\nTop 25 patient-level differential genes:\n")

top_results_print <- fibroblast_patient_level_DE %>%
  dplyr::select(
    Gene,
    log2FC_Male_vs_Female,
    Average_LogCPM,
    F,
    Raw_P_Value,
    FDR,
    Direction
  ) %>%
  head(25)

print(
  as.data.frame(top_results_print),
  row.names = FALSE
)

# ============================================================
# 14. CHECK SPECIFIC GENES OF INTEREST
# ============================================================

genes_of_interest <- c(
  "COL3A1",
  "COL5A2",
  "COL12A1",
  "INPP4B"
)

genes_of_interest_results <- fibroblast_patient_level_DE %>%
  dplyr::filter(
    Gene %in% genes_of_interest
  ) %>%
  dplyr::select(
    Gene,
    log2FC_Male_vs_Female,
    Average_LogCPM,
    Raw_P_Value,
    FDR,
    Direction
  )

cat("\nGenes of interest in patient-level fibroblast DE:\n")
print(
  as.data.frame(genes_of_interest_results),
  row.names = FALSE
)

# ============================================================
# 15. PATIENT-LEVEL VOLCANO PLOT
# ============================================================

volcano_data <- fibroblast_patient_level_DE %>%
  dplyr::mutate(
    FDR_for_plot = pmax(
      FDR,
      .Machine$double.xmin
    ),
    Negative_Log10_FDR = -log10(FDR_for_plot)
  )

top_label_genes <- volcano_data %>%
  dplyr::filter(
    Direction != "Not significant"
  ) %>%
  dplyr::arrange(FDR) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::pull(Gene)

volcano_data <- volcano_data %>%
  dplyr::mutate(
    Label = ifelse(
      Gene %in% top_label_genes,
      Gene,
      NA_character_
    )
  )

patient_level_volcano <- ggplot(
  volcano_data,
  aes(
    x = log2FC_Male_vs_Female,
    y = Negative_Log10_FDR
  )
) +
  geom_point(
    aes(color = Direction),
    alpha = 0.65,
    size = 1.7
  ) +
  geom_vline(
    xintercept = c(-0.5, 0.5),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  ggrepel::geom_text_repel(
    aes(label = Label),
    max.overlaps = 20,
    size = 3.4,
    box.padding = 0.4,
    point.padding = 0.3,
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = c(
      "Female-upregulated" = "#B2182B",
      "Male-upregulated" = "#2166AC",
      "Not significant" = "grey70"
    )
  ) +
  labs(
    title = "Patient-Level Pseudobulk Differential Expression",
    subtitle = paste0(
      "PDAC fibroblasts: ",
      sum(patient_metadata$Sex == "Male"),
      " male versus ",
      sum(patient_metadata$Sex == "Female"),
      " female patients"
    ),
    x = "log2 Fold Change: Male versus Female",
    y = "-log10(FDR)",
    color = ""
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 11
    ),
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    legend.position = "bottom",
    legend.text = element_text(
      face = "bold"
    )
  )

# Print plot directly in RStudio
print(patient_level_volcano)

# ============================================================
# 16. FINAL MESSAGE
# ============================================================

cat("\n======================================================\n")
cat("PATIENT-LEVEL PSEUDOBULK FIBROBLAST DE COMPLETE\n")
cat("======================================================\n")
cat("Statistical unit used: PATIENT\n")
cat("Number of male patients:", sum(patient_metadata$Sex == "Male"), "\n")
cat("Number of female patients:", sum(patient_metadata$Sex == "Female"), "\n")
cat("No files were saved. Plot was printed in RStudio.\n")
cat("======================================================\n")


### PSC Analysis Patient Level



# ============================================================
# DONOR-LEVEL PSC-LIKE FIBROBLAST COMPOSITION ANALYSIS
#
# PSC-like cells = fibroblast subset cluster 3
#
# Includes:
# 1. Donor-level compositional model
# 2. PSC-like percent of all nuclei per patient
# 3. PSC-like percent within fibroblasts per patient
# 4. Raw patient points
# 5. Error bars = mean ± SEM
# 6. BH multiplicity correction
# 7. Leave-one-patient-out sensitivity analysis
# 8. Separate plot of PSC-like nuclei count per patient
#
# This version prints plots in RStudio and does not save files.
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

options(stringsAsFactors = FALSE)

# ============================================================
# 1. PATHS
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

merged_path <- file.path(
  base_dir,
  "Merged_Joined.RDS"
)

fibroblast_path <- file.path(
  base_dir,
  "Fibroblast_Subset.RDS"
)

if (!file.exists(merged_path)) {
  stop("Merged_Joined.RDS was not found.")
}

if (!file.exists(fibroblast_path)) {
  stop("Fibroblast_Subset.RDS was not found.")
}

# ============================================================
# 2. SEX AND PATIENT MAPPING
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

patient_mapping <- c(
  "Sample_01" = "P01",
  "Sample_02" = "P02",
  "Sample_03" = "P03",
  "Sample_04" = "P04",
  "Sample_05" = "P05",
  "Sample_06" = "P06",
  "Sample_07" = "P07",
  "Sample_08" = "P08",
  "Sample_09" = "P09",
  "Sample_10" = "P10",
  "Sample_11" = "P11",
  "Sample_12" = "P12",
  "Sample_13" = "P13",
  "Sample_14" = "P14",
  "Sample_15" = "P15",
  "Sample_16" = "P16",
  "Sample_17" = "P17"
)

sample_base <- tibble::tibble(
  orig.ident = names(sex_mapping),
  Patient_ID = unname(patient_mapping[names(sex_mapping)]),
  Sex = unname(sex_mapping[names(sex_mapping)])
) %>%
  dplyr::mutate(
    Patient_ID = factor(Patient_ID, levels = sprintf("P%02d", 1:17)),
    Sex = factor(Sex, levels = c("Female", "Male"))
  ) %>%
  dplyr::arrange(Patient_ID)

# ============================================================
# 3. LOAD OBJECTS
# ============================================================

Merged_Joined <- readRDS(merged_path)
fibroblast_subset <- readRDS(fibroblast_path)

DefaultAssay(Merged_Joined) <- "RNA"
DefaultAssay(fibroblast_subset) <- "RNA"

Merged_Joined <- tryCatch(
  JoinLayers(Merged_Joined),
  error = function(e) Merged_Joined
)

fibroblast_subset <- tryCatch(
  JoinLayers(fibroblast_subset),
  error = function(e) fibroblast_subset
)

# ============================================================
# 4. CREATE METADATA TABLES
# ============================================================

merged_meta <- Merged_Joined@meta.data %>%
  tibble::rownames_to_column("Cell_ID") %>%
  dplyr::mutate(
    orig.ident = as.character(orig.ident)
  ) %>%
  dplyr::filter(
    orig.ident %in% names(sex_mapping)
  )

fibroblast_meta <- fibroblast_subset@meta.data %>%
  tibble::rownames_to_column("Cell_ID") %>%
  dplyr::mutate(
    orig.ident = as.character(orig.ident),
    seurat_clusters = as.character(seurat_clusters)
  ) %>%
  dplyr::filter(
    orig.ident %in% names(sex_mapping)
  )

# ============================================================
# 5. DONOR-LEVEL PSC-LIKE COUNTS AND PERCENTAGES
# ============================================================

total_nuclei_table <- merged_meta %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(
    Total_Nuclei = dplyr::n(),
    .groups = "drop"
  )

fibroblast_nuclei_table <- fibroblast_meta %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(
    Total_Fibroblast_Nuclei = dplyr::n(),
    .groups = "drop"
  )

psc_nuclei_table <- fibroblast_meta %>%
  dplyr::filter(
    seurat_clusters == "3"
  ) %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(
    PSC_like_Nuclei = dplyr::n(),
    .groups = "drop"
  )

psc_patient_table <- sample_base %>%
  dplyr::left_join(
    total_nuclei_table,
    by = "orig.ident"
  ) %>%
  dplyr::left_join(
    fibroblast_nuclei_table,
    by = "orig.ident"
  ) %>%
  dplyr::left_join(
    psc_nuclei_table,
    by = "orig.ident"
  ) %>%
  dplyr::mutate(
    Total_Nuclei = ifelse(is.na(Total_Nuclei), 0, Total_Nuclei),
    Total_Fibroblast_Nuclei = ifelse(
      is.na(Total_Fibroblast_Nuclei),
      0,
      Total_Fibroblast_Nuclei
    ),
    PSC_like_Nuclei = ifelse(is.na(PSC_like_Nuclei), 0, PSC_like_Nuclei),
    
    Percent_PSC_of_All_Nuclei = ifelse(
      Total_Nuclei > 0,
      100 * PSC_like_Nuclei / Total_Nuclei,
      NA_real_
    ),
    
    Percent_PSC_within_Fibroblasts = ifelse(
      Total_Fibroblast_Nuclei > 0,
      100 * PSC_like_Nuclei / Total_Fibroblast_Nuclei,
      NA_real_
    ),
    
    PSC_Proportion_of_All_Nuclei = ifelse(
      Total_Nuclei > 0,
      PSC_like_Nuclei / Total_Nuclei,
      NA_real_
    ),
    
    PSC_Proportion_within_Fibroblasts = ifelse(
      Total_Fibroblast_Nuclei > 0,
      PSC_like_Nuclei / Total_Fibroblast_Nuclei,
      NA_real_
    )
  ) %>%
  dplyr::arrange(Patient_ID)

cat("\n======================================================\n")
cat("DONOR-LEVEL PSC-LIKE PATIENT TABLE\n")
cat("======================================================\n")

print(
  as.data.frame(psc_patient_table),
  row.names = FALSE
)

# ============================================================
# 6. SUMMARY TABLE FOR ERROR BARS
# Error bars are explicitly mean ± SEM.
# ============================================================

psc_long <- psc_patient_table %>%
  dplyr::select(
    Patient_ID,
    orig.ident,
    Sex,
    Percent_PSC_of_All_Nuclei,
    Percent_PSC_within_Fibroblasts
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      Percent_PSC_of_All_Nuclei,
      Percent_PSC_within_Fibroblasts
    ),
    names_to = "Metric",
    values_to = "PSC_Percent"
  ) %>%
  dplyr::mutate(
    Metric = dplyr::case_when(
      Metric == "Percent_PSC_of_All_Nuclei" ~
        "Percent of PSC-like Fibroblasts from all nuclei",
      Metric == "Percent_PSC_within_Fibroblasts" ~
        "Percent of PSC-like Fibroblasts Within fibroblasts",
      TRUE ~ Metric
    )
  )

psc_summary <- psc_long %>%
  dplyr::group_by(Metric, Sex) %>%
  dplyr::summarise(
    N_Patients = dplyr::n(),
    Mean = mean(PSC_Percent, na.rm = TRUE),
    SD = sd(PSC_Percent, na.rm = TRUE),
    SEM = SD / sqrt(N_Patients),
    Lower_SEM = Mean - SEM,
    Upper_SEM = Mean + SEM,
    .groups = "drop"
  )

cat("\n======================================================\n")
cat("PSC-LIKE SUMMARY BY SEX\n")
cat("Error bars = mean ± SEM\n")
cat("======================================================\n")

print(
  as.data.frame(psc_summary),
  row.names = FALSE
)

# ============================================================
# 7. DONOR-LEVEL COMPOSITIONAL MODELS
# ============================================================
# Uses empirical logit transformation:
# logit((PSC_like_Nuclei + 0.5) / (denominator + 1))
#
# This keeps the patient as the statistical unit.
# SexMale coefficient = Male - Female on logit scale.
# Effect size reported below = Female mean - Male mean in percentage points.
# ============================================================

fit_compositional_model <- function(
    data,
    metric_name,
    numerator_col,
    denominator_col,
    percent_col
) {
  
  model_data <- data %>%
    dplyr::mutate(
      Numerator = .data[[numerator_col]],
      Denominator = .data[[denominator_col]],
      Percent = .data[[percent_col]],
      Empirical_Proportion = (Numerator + 0.5) / (Denominator + 1),
      Empirical_Logit = log(
        Empirical_Proportion / (1 - Empirical_Proportion)
      )
    ) %>%
    dplyr::filter(
      !is.na(Empirical_Logit),
      is.finite(Empirical_Logit),
      !is.na(Sex)
    )
  
  fit <- lm(
    Empirical_Logit ~ Sex,
    data = model_data
  )
  
  fit_summary <- summary(fit)
  
  p_value <- fit_summary$coefficients["SexMale", "Pr(>|t|)"]
  beta_male_vs_female <- fit_summary$coefficients["SexMale", "Estimate"]
  
  female_mean_percent <- mean(
    model_data$Percent[model_data$Sex == "Female"],
    na.rm = TRUE
  )
  
  male_mean_percent <- mean(
    model_data$Percent[model_data$Sex == "Male"],
    na.rm = TRUE
  )
  
  female_minus_male_percent <- female_mean_percent - male_mean_percent
  
  out <- data.frame(
    Metric = metric_name,
    Model = "lm(empirical logit donor-level proportion ~ Sex)",
    N_Female_Patients = sum(model_data$Sex == "Female"),
    N_Male_Patients = sum(model_data$Sex == "Male"),
    Female_Mean_Percent = female_mean_percent,
    Male_Mean_Percent = male_mean_percent,
    Effect_Size_Female_Minus_Male_Percentage_Points =
      female_minus_male_percent,
    Beta_Male_vs_Female_Logit_Scale = beta_male_vs_female,
    Raw_P_Value = p_value,
    stringsAsFactors = FALSE
  )
  
  return(out)
}

model_results <- dplyr::bind_rows(
  fit_compositional_model(
    data = psc_patient_table,
    metric_name = "PSC-like % of all nuclei",
    numerator_col = "PSC_like_Nuclei",
    denominator_col = "Total_Nuclei",
    percent_col = "Percent_PSC_of_All_Nuclei"
  ),
  fit_compositional_model(
    data = psc_patient_table,
    metric_name = "PSC-like % within fibroblasts",
    numerator_col = "PSC_like_Nuclei",
    denominator_col = "Total_Fibroblast_Nuclei",
    percent_col = "Percent_PSC_within_Fibroblasts"
  )
) %>%
  dplyr::mutate(
    BH_Adjusted_P_Value = p.adjust(
      Raw_P_Value,
      method = "BH"
    )
  )

cat("\n======================================================\n")
cat("DONOR-LEVEL COMPOSITIONAL MODEL RESULTS\n")
cat("Multiplicity correction: BH across the two PSC composition tests\n")
cat("======================================================\n")

print(
  as.data.frame(model_results),
  row.names = FALSE
)

# ============================================================
# 8. LEAVE-ONE-PATIENT-OUT ANALYSIS
# ============================================================

run_leave_one_out <- function(
    data,
    metric_name,
    numerator_col,
    denominator_col,
    percent_col
) {
  
  patient_ids <- as.character(data$Patient_ID)
  
  loo_results <- lapply(
    patient_ids,
    function(left_out_patient) {
      
      loo_data <- data %>%
        dplyr::filter(
          as.character(Patient_ID) != left_out_patient
        ) %>%
        dplyr::mutate(
          Numerator = .data[[numerator_col]],
          Denominator = .data[[denominator_col]],
          Percent = .data[[percent_col]],
          Empirical_Proportion = (Numerator + 0.5) / (Denominator + 1),
          Empirical_Logit = log(
            Empirical_Proportion / (1 - Empirical_Proportion)
          )
        ) %>%
        dplyr::filter(
          !is.na(Empirical_Logit),
          is.finite(Empirical_Logit),
          !is.na(Sex)
        )
      
      if (length(unique(loo_data$Sex)) < 2) {
        return(NULL)
      }
      
      fit <- lm(
        Empirical_Logit ~ Sex,
        data = loo_data
      )
      
      fit_summary <- summary(fit)
      
      p_value <- fit_summary$coefficients["SexMale", "Pr(>|t|)"]
      beta_male_vs_female <- fit_summary$coefficients["SexMale", "Estimate"]
      
      female_mean_percent <- mean(
        loo_data$Percent[loo_data$Sex == "Female"],
        na.rm = TRUE
      )
      
      male_mean_percent <- mean(
        loo_data$Percent[loo_data$Sex == "Male"],
        na.rm = TRUE
      )
      
      effect_female_minus_male <- female_mean_percent - male_mean_percent
      
      data.frame(
        Metric = metric_name,
        Left_Out_Patient = left_out_patient,
        N_Female_Patients = sum(loo_data$Sex == "Female"),
        N_Male_Patients = sum(loo_data$Sex == "Male"),
        Female_Mean_Percent = female_mean_percent,
        Male_Mean_Percent = male_mean_percent,
        Effect_Size_Female_Minus_Male_Percentage_Points =
          effect_female_minus_male,
        Beta_Male_vs_Female_Logit_Scale = beta_male_vs_female,
        Raw_P_Value = p_value,
        Female_Mean_Remains_Higher = effect_female_minus_male > 0,
        Result_Reversed = effect_female_minus_male < 0,
        stringsAsFactors = FALSE
      )
    }
  )
  
  loo_results <- dplyr::bind_rows(loo_results)
  
  loo_summary <- loo_results %>%
    dplyr::group_by(Metric) %>%
    dplyr::summarise(
      Number_of_Leave_One_Out_Iterations = dplyr::n(),
      Effect_Size_Range_Min = min(
        Effect_Size_Female_Minus_Male_Percentage_Points,
        na.rm = TRUE
      ),
      Effect_Size_Range_Max = max(
        Effect_Size_Female_Minus_Male_Percentage_Points,
        na.rm = TRUE
      ),
      Female_Mean_Remains_Higher_N_Iterations = sum(
        Female_Mean_Remains_Higher,
        na.rm = TRUE
      ),
      P_Value_Range_Min = min(Raw_P_Value, na.rm = TRUE),
      P_Value_Range_Max = max(Raw_P_Value, na.rm = TRUE),
      Any_Patient_Reverses_Result = any(Result_Reversed, na.rm = TRUE),
      Patients_That_Reverse_Result = paste(
        Left_Out_Patient[Result_Reversed],
        collapse = ", "
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Patients_That_Reverse_Result = ifelse(
        Patients_That_Reverse_Result == "",
        "None",
        Patients_That_Reverse_Result
      )
    )
  
  return(
    list(
      detailed = loo_results,
      summary = loo_summary
    )
  )
}

loo_all_nuclei <- run_leave_one_out(
  data = psc_patient_table,
  metric_name = "PSC-like % of all nuclei",
  numerator_col = "PSC_like_Nuclei",
  denominator_col = "Total_Nuclei",
  percent_col = "Percent_PSC_of_All_Nuclei"
)

loo_within_fibroblasts <- run_leave_one_out(
  data = psc_patient_table,
  metric_name = "PSC-like % within fibroblasts",
  numerator_col = "PSC_like_Nuclei",
  denominator_col = "Total_Fibroblast_Nuclei",
  percent_col = "Percent_PSC_within_Fibroblasts"
)

loo_detailed <- dplyr::bind_rows(
  loo_all_nuclei$detailed,
  loo_within_fibroblasts$detailed
)

loo_summary <- dplyr::bind_rows(
  loo_all_nuclei$summary,
  loo_within_fibroblasts$summary
)

cat("\n======================================================\n")
cat("LEAVE-ONE-PATIENT-OUT SUMMARY\n")
cat("======================================================\n")

print(
  as.data.frame(loo_summary),
  row.names = FALSE
)

cat("\n======================================================\n")
cat("LEAVE-ONE-PATIENT-OUT DETAILED RESULTS\n")
cat("======================================================\n")

print(
  as.data.frame(loo_detailed),
  row.names = FALSE
)

# ============================================================
# 9. PLOT 1:
# PSC-LIKE COMPOSITION USING BOTH DENOMINATORS
# Raw patient points + mean ± SEM
# ============================================================

model_labels <- model_results %>%
  dplyr::mutate(
    Metric = as.character(Metric),
    Label = paste0(
      "BH p = ",
      formatC(
        BH_Adjusted_P_Value,
        format = "e",
        digits = 2
      ),
      "\nError bars = mean \u00B1 SEM"
    )
  ) %>%
  dplyr::select(
    Metric,
    Label
  )

label_positions <- psc_long %>%
  dplyr::group_by(Metric) %>%
  dplyr::summarise(
    y_position = max(PSC_Percent, na.rm = TRUE) * 1.18,
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    model_labels,
    by = "Metric"
  )

composition_plot <- ggplot(
  psc_long,
  aes(
    x = Sex,
    y = PSC_Percent,
    color = Sex
  )
) +
  geom_jitter(
    width = 0.08,
    height = 0,
    size = 3,
    alpha = 0.9
  ) +
  geom_errorbar(
    data = psc_summary,
    aes(
      x = Sex,
      ymin = Lower_SEM,
      ymax = Upper_SEM
    ),
    inherit.aes = FALSE,
    width = 0.18,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_point(
    data = psc_summary,
    aes(
      x = Sex,
      y = Mean
    ),
    inherit.aes = FALSE,
    size = 4,
    shape = 18,
    color = "black"
  ) +
  geom_text(
    data = label_positions,
    aes(
      x = 1.5,
      y = y_position,
      label = Label
    ),
    inherit.aes = FALSE,
    size = 3.4,
    color = "black"
  ) +
  facet_wrap(
    ~ Metric,
    scales = "free_y"
  ) +
  scale_color_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  labs(
    title = "Donor-Level PSC-like Fibroblast Composition",
    subtitle = "Each point is one patient; error bars show mean \u00B1 SEM",
    x = "",
    y = "PSC-like nuclei (%)",
    color = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 11
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    axis.title = element_text(
      face = "bold",
      size = 12
    ),
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    legend.position = "bottom"
  )

print(composition_plot)

# ============================================================
# 10. PLOT 2:
# PSC-LIKE NUCLEI COUNT IN EVERY PATIENT
# ============================================================

count_plot <- ggplot(
  psc_patient_table,
  aes(
    x = Patient_ID,
    y = PSC_like_Nuclei,
    fill = Sex
  )
) +
  geom_col(
    width = 0.75,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(
      label = PSC_like_Nuclei
    ),
    vjust = -0.25,
    size = 3
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  labs(
    title = "PSC-like Nuclei Count Per Patient",
    subtitle = "PSC-like cells defined as fibroblast subset cluster 3",
    x = "Patient",
    y = "Number of PSC-like nuclei",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 11
    ),
    axis.title = element_text(
      face = "bold",
      size = 12
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 10,
      color = "black"
    ),
    legend.position = "bottom"
  )

print(count_plot)

# ============================================================
# 11. FINAL CONSOLE SUMMARY
# ============================================================

cat("\n======================================================\n")
cat("PSC-LIKE DONOR-LEVEL ANALYSIS COMPLETE\n")
cat("======================================================\n")
cat("Statistical unit: patient/sample\n")
cat("PSC-like definition: fibroblast subset seurat cluster 3\n")
cat("Composition model: empirical logit donor-level proportion ~ Sex\n")
cat("Error bars in composition plot: mean ± SEM\n")
cat("Multiplicity correction: BH correction across two composition tests\n")
cat("Plots printed in RStudio. No files were saved.\n")
cat("======================================================\n")


# ============================================================
# CLEAN LEAVE-ONE-PATIENT-OUT REPORT
# ============================================================

cat("\n======================================================\n")
cat("LEAVE-ONE-PATIENT-OUT REPORT\n")
cat("======================================================\n")

for (i in seq_len(nrow(loo_summary))) {
  
  cat("\nMetric:", loo_summary$Metric[i], "\n")
  cat("Number of leave-one-patient-out iterations:",
      loo_summary$Number_of_Leave_One_Out_Iterations[i], "\n")
  
  cat("Range of effect-size estimates, female minus male percentage points:",
      round(loo_summary$Effect_Size_Range_Min[i], 4),
      "to",
      round(loo_summary$Effect_Size_Range_Max[i], 4),
      "\n")
  
  cat("Number of iterations in which female mean remains higher:",
      loo_summary$Female_Mean_Remains_Higher_N_Iterations[i],
      "of",
      loo_summary$Number_of_Leave_One_Out_Iterations[i],
      "\n")
  
  cat("Range of p-values:",
      formatC(loo_summary$P_Value_Range_Min[i], format = "e", digits = 2),
      "to",
      formatC(loo_summary$P_Value_Range_Max[i], format = "e", digits = 2),
      "\n")
  
  cat("Does leaving out one patient reverse the result?:",
      loo_summary$Any_Patient_Reverses_Result[i],
      "\n")
  
  cat("Patients that reverse the result:",
      loo_summary$Patients_That_Reverse_Result[i],
      "\n")
}

cat("\n======================================================\n")



### PSC counting

# ============================================================
# PSC-LIKE NUCLEI COUNT: FEMALE VS MALE AVERAGE
#
# Goal:
# - Calculate PSC-like nuclei count per patient
# - Compare average PSC-like nuclei count between females and males
# - Plot mean PSC-like nuclei count by sex
# - Show raw patient points
# - Add p-value
#
# PSC-like cells = fibroblast subset cluster 3
# Female = red, Male = blue
# Statistical unit = patient/sample
# ============================================================

library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)

options(stringsAsFactors = FALSE)

# ============================================================
# 1. LOAD FIBROBLAST SUBSET
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

fibroblast_path <- file.path(
  base_dir,
  "Fibroblast_Subset.RDS"
)

if (!file.exists(fibroblast_path)) {
  stop("Fibroblast_Subset.RDS was not found.")
}

fibroblast_subset <- readRDS(fibroblast_path)

DefaultAssay(fibroblast_subset) <- "RNA"

fibroblast_subset <- tryCatch(
  JoinLayers(fibroblast_subset),
  error = function(e) fibroblast_subset
)

# ============================================================
# 2. SEX AND PATIENT MAPPING
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

patient_mapping <- c(
  "Sample_01" = "P01",
  "Sample_02" = "P02",
  "Sample_03" = "P03",
  "Sample_04" = "P04",
  "Sample_05" = "P05",
  "Sample_06" = "P06",
  "Sample_07" = "P07",
  "Sample_08" = "P08",
  "Sample_09" = "P09",
  "Sample_10" = "P10",
  "Sample_11" = "P11",
  "Sample_12" = "P12",
  "Sample_13" = "P13",
  "Sample_14" = "P14",
  "Sample_15" = "P15",
  "Sample_16" = "P16",
  "Sample_17" = "P17"
)

# ============================================================
# 3. CREATE PATIENT-LEVEL PSC COUNT TABLE
# ============================================================

fibroblast_meta <- fibroblast_subset@meta.data %>%
  tibble::rownames_to_column("Cell_ID") %>%
  dplyr::mutate(
    orig.ident = as.character(orig.ident),
    seurat_clusters = as.character(seurat_clusters),
    Sex = unname(sex_mapping[orig.ident]),
    Patient_ID = unname(patient_mapping[orig.ident])
  ) %>%
  dplyr::filter(
    !is.na(Sex),
    !is.na(Patient_ID)
  )

sample_base <- tibble::tibble(
  orig.ident = names(sex_mapping),
  Patient_ID = unname(patient_mapping[names(sex_mapping)]),
  Sex = unname(sex_mapping[names(sex_mapping)])
) %>%
  dplyr::mutate(
    Patient_ID = factor(Patient_ID, levels = sprintf("P%02d", 1:17)),
    Sex = factor(Sex, levels = c("Female", "Male"))
  )

psc_count_per_patient <- fibroblast_meta %>%
  dplyr::filter(
    seurat_clusters == "3"
  ) %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(
    PSC_like_Nuclei = dplyr::n(),
    .groups = "drop"
  )

psc_patient_table <- sample_base %>%
  dplyr::left_join(
    psc_count_per_patient,
    by = "orig.ident"
  ) %>%
  dplyr::mutate(
    PSC_like_Nuclei = ifelse(
      is.na(PSC_like_Nuclei),
      0,
      PSC_like_Nuclei
    )
  ) %>%
  dplyr::arrange(Patient_ID)

cat("\n======================================================\n")
cat("PSC-LIKE NUCLEI COUNT PER PATIENT\n")
cat("======================================================\n")

print(
  as.data.frame(psc_patient_table),
  row.names = FALSE
)

# ============================================================
# 4. SEX-LEVEL SUMMARY: MEAN ± SEM
# ============================================================

psc_summary <- psc_patient_table %>%
  dplyr::group_by(Sex) %>%
  dplyr::summarise(
    N_Patients = dplyr::n(),
    Mean_PSC_like_Nuclei = mean(PSC_like_Nuclei, na.rm = TRUE),
    SD_PSC_like_Nuclei = sd(PSC_like_Nuclei, na.rm = TRUE),
    SEM_PSC_like_Nuclei = SD_PSC_like_Nuclei / sqrt(N_Patients),
    Lower_SEM = Mean_PSC_like_Nuclei - SEM_PSC_like_Nuclei,
    Upper_SEM = Mean_PSC_like_Nuclei + SEM_PSC_like_Nuclei,
    .groups = "drop"
  )

cat("\n======================================================\n")
cat("PSC-LIKE NUCLEI COUNT SUMMARY BY SEX\n")
cat("Error bars = mean ± SEM\n")
cat("======================================================\n")

print(
  as.data.frame(psc_summary),
  row.names = FALSE
)

# ============================================================
# 5. T-TEST: FEMALE VS MALE PSC-LIKE NUCLEI COUNT
# ============================================================

psc_ttest <- t.test(
  PSC_like_Nuclei ~ Sex,
  data = psc_patient_table
)

p_value <- psc_ttest$p.value
p_label <- paste0(
  "Welch t-test p = ",
  signif(p_value, 3)
)

cat("\n======================================================\n")
cat("PSC-LIKE NUCLEI COUNT T-TEST\n")
cat("======================================================\n")
print(psc_ttest)

# ============================================================
# 6. PLOT: AVERAGE PSC-LIKE NUCLEI COUNT BY SEX
# ============================================================

y_max <- max(
  psc_patient_table$PSC_like_Nuclei,
  psc_summary$Upper_SEM,
  na.rm = TRUE
)

psc_average_count_plot <- ggplot(
  psc_summary,
  aes(
    x = Sex,
    y = Mean_PSC_like_Nuclei,
    fill = Sex
  )
) +
  geom_col(
    width = 0.65,
    color = "black",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = Lower_SEM,
      ymax = Upper_SEM
    ),
    width = 0.18,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_jitter(
    data = psc_patient_table,
    aes(
      x = Sex,
      y = PSC_like_Nuclei,
      color = Sex
    ),
    inherit.aes = FALSE,
    width = 0.08,
    size = 3,
    alpha = 0.9
  ) +
  annotate(
    "text",
    x = 1.5,
    y = y_max * 1.15,
    label = p_label,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  scale_color_manual(
    values = c(
      "Female" = "#B2182B",
      "Male" = "#0072B2"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title = "Average PSC-like Nuclei Count by Sex",
    subtitle = "Each point represents one patient; error bars show mean ± SEM",
    x = "",
    y = "Average number of PSC-like nuclei",
    fill = "Sex",
    color = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 17
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 12
    ),
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(
      size = 12,
      color = "black",
      face = "bold"
    ),
    legend.position = "bottom",
    legend.title = element_text(
      face = "bold"
    ),
    legend.text = element_text(
      size = 11
    )
  )

print(psc_average_count_plot)

# ============================================================
# 7. FINAL MESSAGE
# ============================================================

cat("\n======================================================\n")
cat("PSC-LIKE AVERAGE COUNT PLOT COMPLETE\n")
cat("======================================================\n")
cat("Statistical unit: patient/sample\n")
cat("PSC-like cells defined as fibroblast subset cluster 3\n")
cat("Female bar: red\n")
cat("Male bar: blue\n")
cat("Error bars: mean ± SEM\n")
cat("P-value: Welch two-sample t-test\n")
cat("======================================================\n")


### Cell Composition FINAL FINAL
# ============================================================
# AVERAGE CELL-TYPE COMPOSITION BY SEX
# Full revised code with attached legend color scheme
#
# Each bar = average patient-level cell-type proportion
# Female and Male bars each sum to 100%
# Statistical unit for averaging = patient/sample
#
# Prints plot in RStudio.
# Does not save files unless optional ggsave lines are uncommented.
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scales)

options(stringsAsFactors = FALSE)

# ============================================================
# 1. LOAD MAIN SEURAT OBJECT
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

merged_path <- file.path(
  base_dir,
  "Merged_Joined.RDS"
)

if (!file.exists(merged_path)) {
  stop("Merged_Joined.RDS was not found.")
}

Merged_Joined <- readRDS(merged_path)

DefaultAssay(Merged_Joined) <- "RNA"

Merged_Joined <- tryCatch(
  JoinLayers(Merged_Joined),
  error = function(e) Merged_Joined
)

# ============================================================
# 2. SEX AND PATIENT MAPPING
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

patient_mapping <- c(
  "Sample_01" = "P01",
  "Sample_02" = "P02",
  "Sample_03" = "P03",
  "Sample_04" = "P04",
  "Sample_05" = "P05",
  "Sample_06" = "P06",
  "Sample_07" = "P07",
  "Sample_08" = "P08",
  "Sample_09" = "P09",
  "Sample_10" = "P10",
  "Sample_11" = "P11",
  "Sample_12" = "P12",
  "Sample_13" = "P13",
  "Sample_14" = "P14",
  "Sample_15" = "P15",
  "Sample_16" = "P16",
  "Sample_17" = "P17"
)

sample_base <- tibble::tibble(
  orig.ident = names(sex_mapping),
  Patient_ID = unname(patient_mapping[names(sex_mapping)]),
  Sex = unname(sex_mapping[names(sex_mapping)])
) %>%
  dplyr::mutate(
    Patient_ID = factor(Patient_ID, levels = sprintf("P%02d", 1:17)),
    Sex = factor(Sex, levels = c("Female", "Male"))
  )

# ============================================================
# 3. DEFINE CELL TYPES
# ============================================================
# This matches the attached legend:
# Cancer 1, Fibroblast, Ductal, B Cells, Myeloid,
# Acinar, Endocrine, Endothelial, T / NK Cells, Cancer 2,
# Smooth Muscle, Adipocytes, Other
#
# Clusters 11, 13, and 14 are grouped as Other here.
# ============================================================

celltype_order <- c(
  "Cancer 1",
  "Fibroblast",
  "Ductal",
  "B Cells",
  "Myeloid",
  "Acinar",
  "Endocrine",
  "Endothelial",
  "T / NK Cells",
  "Cancer 2",
  "Smooth Muscle",
  "Adipocytes",
  "Other"
)

plot_data <- Merged_Joined@meta.data %>%
  tibble::rownames_to_column("Cell_ID") %>%
  dplyr::mutate(
    orig.ident = as.character(orig.ident),
    seurat_clusters = as.character(seurat_clusters),
    Sex = unname(sex_mapping[orig.ident]),
    Patient_ID = unname(patient_mapping[orig.ident]),
    
    CellType = dplyr::case_when(
      seurat_clusters == "0" ~ "Cancer 1",
      seurat_clusters == "1" ~ "Fibroblast",
      seurat_clusters == "2" ~ "Ductal",
      seurat_clusters == "3" ~ "B Cells",
      seurat_clusters == "4" ~ "Myeloid",
      seurat_clusters == "5" ~ "Acinar",
      seurat_clusters == "6" ~ "Endocrine",
      seurat_clusters == "7" ~ "Endothelial",
      seurat_clusters == "8" ~ "T / NK Cells",
      seurat_clusters == "9" ~ "Cancer 2",
      seurat_clusters == "10" ~ "Smooth Muscle",
      seurat_clusters == "12" ~ "Adipocytes",
      seurat_clusters %in% c("11", "13", "14") ~ "Other",
      TRUE ~ "Other"
    ),
    
    Sex = factor(Sex, levels = c("Female", "Male")),
    Patient_ID = factor(Patient_ID, levels = sprintf("P%02d", 1:17)),
    CellType = factor(CellType, levels = celltype_order)
  ) %>%
  dplyr::filter(
    !is.na(Sex),
    !is.na(Patient_ID),
    !is.na(CellType)
  )

# ============================================================
# 4. CALCULATE PATIENT-LEVEL CELL-TYPE PROPORTIONS
# ============================================================

patient_totals <- plot_data %>%
  dplyr::group_by(
    orig.ident,
    Patient_ID,
    Sex
  ) %>%
  dplyr::summarise(
    Total_Nuclei = dplyr::n(),
    .groups = "drop"
  )

patient_celltype_counts <- plot_data %>%
  dplyr::group_by(
    orig.ident,
    Patient_ID,
    Sex,
    CellType
  ) %>%
  dplyr::summarise(
    Cell_Count = dplyr::n(),
    .groups = "drop"
  )

patient_celltype_proportions <- tidyr::expand_grid(
  sample_base,
  CellType = factor(celltype_order, levels = celltype_order)
) %>%
  dplyr::left_join(
    patient_celltype_counts,
    by = c("orig.ident", "Patient_ID", "Sex", "CellType")
  ) %>%
  dplyr::left_join(
    patient_totals,
    by = c("orig.ident", "Patient_ID", "Sex")
  ) %>%
  dplyr::mutate(
    Cell_Count = ifelse(is.na(Cell_Count), 0, Cell_Count),
    Proportion = ifelse(
      Total_Nuclei > 0,
      Cell_Count / Total_Nuclei,
      NA_real_
    ),
    Percent = 100 * Proportion
  ) %>%
  dplyr::arrange(
    Patient_ID,
    CellType
  )

cat("\n======================================================\n")
cat("PATIENT-LEVEL CELL-TYPE PROPORTIONS\n")
cat("======================================================\n")

print(
  as.data.frame(patient_celltype_proportions),
  row.names = FALSE
)

# ============================================================
# 5. AVERAGE PATIENT-LEVEL PROPORTION BY SEX
# ============================================================

average_celltype_by_sex <- patient_celltype_proportions %>%
  dplyr::group_by(
    Sex,
    CellType
  ) %>%
  dplyr::summarise(
    Mean_Percent = mean(Percent, na.rm = TRUE),
    SD_Percent = sd(Percent, na.rm = TRUE),
    N_Patients = sum(!is.na(Percent)),
    SEM_Percent = SD_Percent / sqrt(N_Patients),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Sex = factor(Sex, levels = c("Female", "Male")),
    CellType = factor(CellType, levels = celltype_order)
  )

cat("\n======================================================\n")
cat("AVERAGE CELL-TYPE PROPORTIONS BY SEX\n")
cat("======================================================\n")

print(
  as.data.frame(average_celltype_by_sex),
  row.names = FALSE
)

cat("\nCheck that each sex sums to approximately 100%:\n")

print(
  average_celltype_by_sex %>%
    dplyr::group_by(Sex) %>%
    dplyr::summarise(
      Total_Percent = sum(Mean_Percent, na.rm = TRUE),
      .groups = "drop"
    )
)

# ============================================================
# 6. COLORS - MATCHING ATTACHED LEGEND
# ============================================================

celltype_colors <- c(
  "Cancer 1" = "#F8766D",
  "Fibroblast" = "#E68613",
  "Ductal" = "#B79F00",
  "B Cells" = "#7CAE00",
  "Myeloid" = "#00A600",
  "Acinar" = "#00BF7D",
  "Endocrine" = "#00BFC4",
  "Endothelial" = "#00A9CF",
  "T / NK Cells" = "#00A6FF",
  "Cancer 2" = "#8494FF",
  "Smooth Muscle" = "#C77CFF",
  "Adipocytes" = "#FF61C3",
  "Other" = "grey75"
)

# ============================================================
# 7. STACKED BAR PLOT
# ============================================================

average_celltype_plot <- ggplot(
  average_celltype_by_sex,
  aes(
    x = Sex,
    y = Mean_Percent,
    fill = CellType
  )
) +
  geom_col(
    width = 0.6,
    color = "white",
    linewidth = 0.45
  ) +
  scale_fill_manual(
    values = celltype_colors,
    breaks = celltype_order,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +
  labs(
    title = "Average Cell-Type Composition by Sex",
    x = "Patient Sex",
    y = "Average Proportion (%)",
    fill = "Cell Type"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 14,
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      face = "bold",
      size = 13,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text = element_text(
      size = 10
    ),
    legend.position = "right",
    legend.key.size = unit(0.55, "cm"),
    panel.border = element_blank()
  ) +
  guides(
    fill = guide_legend(
      ncol = 1,
      byrow = TRUE
    )
  )

print(average_celltype_plot)

# ============================================================
# 8. OPTIONAL: SAVE FIGURE
# Uncomment if you want to save the plot.
# ============================================================

# ggsave(
#   file.path(base_dir, "Average_CellType_Composition_By_Sex.pdf"),
#   plot = average_celltype_plot,
#   width = 10,
#   height = 6,
#   dpi = 300
# )
#
# ggsave(
#   file.path(base_dir, "Average_CellType_Composition_By_Sex.png"),
#   plot = average_celltype_plot,
#   width = 10,
#   height = 6,
#   dpi = 300
# )

# ============================================================
# 9. FINAL MESSAGE
# ============================================================

cat("\n======================================================\n")
cat("AVERAGE CELL-TYPE COMPOSITION GRAPH COMPLETE\n")
cat("======================================================\n")
cat("Statistical unit for averaging: patient/sample\n")
cat("Each bar shows the average patient-level cell-type proportion by sex.\n")
cat("Color scheme matches the attached legend.\n")
cat("Plot printed in RStudio.\n")
cat("======================================================\n")

### Finding Log2fold Change

# ============================================================
# FIND LOG2 FOLD CHANGE OF INPP4B IN FIBROBLAST DGE
#
# Comparison:
# Male fibroblasts vs Female fibroblasts
#
# Interpretation:
# avg_log2FC > 0 = higher in Male fibroblasts
# avg_log2FC < 0 = higher in Female fibroblasts
#
# Object:
# Fibroblast_Subset.RDS
# ============================================================

library(Seurat)
library(dplyr)
library(tibble)

options(stringsAsFactors = FALSE)

# ============================================================
# 1. LOAD FIBROBLAST SUBSET
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

fibroblast_path <- file.path(
  base_dir,
  "Fibroblast_Subset.RDS"
)

if (!file.exists(fibroblast_path)) {
  stop("Fibroblast_Subset.RDS was not found.")
}

fibroblast_subset <- readRDS(fibroblast_path)

DefaultAssay(fibroblast_subset) <- "RNA"

fibroblast_subset <- tryCatch(
  JoinLayers(fibroblast_subset),
  error = function(e) fibroblast_subset
)

# ============================================================
# 2. ADD SEX METADATA
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",
  "Sample_02" = "Male",
  "Sample_03" = "Male",
  "Sample_04" = "Female",
  "Sample_05" = "Female",
  "Sample_06" = "Female",
  "Sample_07" = "Male",
  "Sample_08" = "Female",
  "Sample_09" = "Male",
  "Sample_10" = "Female",
  "Sample_11" = "Female",
  "Sample_12" = "Male",
  "Sample_13" = "Female",
  "Sample_14" = "Male",
  "Sample_15" = "Female",
  "Sample_16" = "Male",
  "Sample_17" = "Female"
)

fibroblast_subset$Sex <- unname(
  sex_mapping[as.character(fibroblast_subset$orig.ident)]
)

fibroblast_subset$Sex <- factor(
  fibroblast_subset$Sex,
  levels = c("Female", "Male")
)

fibroblast_subset <- subset(
  fibroblast_subset,
  subset = Sex %in% c("Female", "Male")
)

cat("\nFibroblast nuclei by sex:\n")
print(table(fibroblast_subset$Sex))

cat("\nFibroblast nuclei by patient and sex:\n")
print(table(fibroblast_subset$orig.ident, fibroblast_subset$Sex))

# ============================================================
# 3. CHECK THAT INPP4B EXISTS
# ============================================================

gene_to_check <- "INPP4B"

if (!gene_to_check %in% rownames(fibroblast_subset)) {
  stop("INPP4B was not found in the fibroblast subset.")
}

# ============================================================
# 4. RUN FIBROBLAST DGE FOR INPP4B ONLY
# ============================================================
# ident.1 = Male
# ident.2 = Female
#
# Therefore:
# positive avg_log2FC = higher in Male fibroblasts
# negative avg_log2FC = higher in Female fibroblasts
# ============================================================

inpp4b_dge <- FindMarkers(
  object = fibroblast_subset,
  ident.1 = "Male",
  ident.2 = "Female",
  group.by = "Sex",
  features = gene_to_check,
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0,
  only.pos = FALSE
)

inpp4b_dge$Gene <- rownames(inpp4b_dge)

inpp4b_result <- inpp4b_dge %>%
  dplyr::select(
    Gene,
    avg_log2FC,
    p_val,
    p_val_adj,
    pct.1,
    pct.2
  ) %>%
  dplyr::mutate(
    Comparison = "Male vs Female fibroblasts",
    Direction = dplyr::case_when(
      avg_log2FC > 0 ~ "Higher in Male fibroblasts",
      avg_log2FC < 0 ~ "Higher in Female fibroblasts",
      TRUE ~ "No log2FC difference"
    ),
    pct.1_label = "Percent expressed in Male fibroblasts",
    pct.2_label = "Percent expressed in Female fibroblasts"
  ) %>%
  dplyr::select(
    Comparison,
    Gene,
    avg_log2FC,
    Direction,
    p_val,
    p_val_adj,
    pct.1,
    pct.1_label,
    pct.2,
    pct.2_label
  )

# ============================================================
# 5. PRINT RESULT
# ============================================================

cat("\n======================================================\n")
cat("INPP4B LOG2 FOLD CHANGE IN FIBROBLASTS\n")
cat("======================================================\n")
cat("Comparison: Male fibroblasts vs Female fibroblasts\n")
cat("Positive avg_log2FC = higher in Male\n")
cat("Negative avg_log2FC = higher in Female\n")
cat("======================================================\n\n")

print(
  as.data.frame(inpp4b_result),
  row.names = FALSE
)

cat("\nINPP4B avg_log2FC value:\n")
print(inpp4b_result$avg_log2FC)

cat("\nDirection:\n")
print(inpp4b_result$Direction)

# ============================================================
# 6. OPTIONAL: SAVE RESULT TO DOWNLOADS
# ============================================================

write.csv(
  as.data.frame(inpp4b_result),
  file = "~/Downloads/INPP4B_Fibroblast_Log2FC_Male_vs_Female.csv",
  row.names = FALSE
)

cat("\nSaved CSV to:\n")
cat("~/Downloads/INPP4B_Fibroblast_Log2FC_Male_vs_Female.csv\n")
cat("======================================================\n")