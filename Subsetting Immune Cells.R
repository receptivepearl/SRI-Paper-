library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)

Merged_Joined <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Joined.RDS")
DefaultAssay(Merged_Joined) <- "RNA"
Idents(Merged_Joined) <- "seurat_clusters"
myeloid_subset <- subset(Merged_Joined, idents = "4")
myeloid_subset$parent_cluster <- myeloid_subset$seurat_clusters
myeloid_subset
table(myeloid_subset$parent_cluster)
DefaultAssay(myeloid_subset) <- "RNA"
myeloid_subset <- JoinLayers(myeloid_subset)
myeloid_subset[["RNA"]] <- split(
  myeloid_subset[["RNA"]],
  f = myeloid_subset$orig.ident
)
myeloid_subset <- NormalizeData(myeloid_subset)
myeloid_subset <- FindVariableFeatures(
  myeloid_subset,
  selection.method = "vst",
  nfeatures = 3000
)
myeloid_subset <- ScaleData(myeloid_subset)
myeloid_subset <- RunPCA(
  myeloid_subset,
  features = VariableFeatures(object = myeloid_subset)
)
myeloid_subset <- IntegrateLayers(
  object = myeloid_subset,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE
)
myeloid_subset <- JoinLayers(myeloid_subset)
myeloid_subset <- FindNeighbors(
  myeloid_subset,
  reduction = "harmony",
  dims = 1:15
)
myeloid_subset <- FindClusters(
  myeloid_subset,
  resolution = 0.1
)
myeloid_subset$myeloid_subcluster <- Idents(myeloid_subset)
myeloid_subset <- RunUMAP(
  myeloid_subset,
  reduction = "harmony",
  dims = 1:15,
  reduction.name = "myeloid_umap"
)
DimPlot(
  myeloid_subset,
  reduction = "myeloid_umap",
  group.by = "myeloid_subcluster",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Harmony-Reintegrated Myeloid Sub-Clusters") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
table(myeloid_subset$myeloid_subcluster)

saveRDS(myeloid_subset, file = '/Users/sudharshiniram/Desktop/SRI Data Analysis/Myeloid_Subset.RDS')

# Broad myeloid subtype markers
DotPlot(
  myeloid_subset,
  features = c(
    "CSF3R", "FCGR3B", "CEACAM8", "S100A8",
    "C1QA", "C1QB", "APOE", "SPP1",
    "CD1C", "CLEC9A", "XCR1", "LILRA4"
  ),
  group.by = "myeloid_subcluster"
) +
  RotatedAxis() +
  ggtitle("Myeloid Subtype Markers") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


# Macrophages / TAMs - Cluster 0 and 1
macrophage_markers <- c(
  "C1QA", "C1QB", "C1QC", "APOE", "APOC1",
  "MSR1", "MRC1", "CD163", "SPP1", "GPNMB", "TREM2"
)
macrophage_markers <- macrophage_markers[
  macrophage_markers %in% rownames(myeloid_subset)
]
DotPlot(
  myeloid_subset,
  features = macrophage_markers,
  group.by = "myeloid_subcluster"
) +
  RotatedAxis() +
  ggtitle("Macrophage / TAM Markers") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


# Dendritic cells - Cluster 3
VlnPlot(
  myeloid_subset,
  features = c("CLEC9A", "XCR1", "BATF3", "IRF8", "CD1C", "FCER1A"),
  group.by = "myeloid_subcluster",
  pt.size = 0
)
VlnPlot(
  myeloid_subset,
  features = c("CLEC10A", "LILRA4", "CLEC4C", "IL3RA", "TCF4", "CCR7"),
  group.by = "myeloid_subcluster",
  pt.size = 0
)
FeaturePlot(
  myeloid_subset,
  reduction = "myeloid_umap",
  features = c("CLEC9A", "XCR1", "BATF3", "IRF8", "CD1C", "LILRA4"),
  cols = c("lightgrey", "blue")
)
DotPlot(
  myeloid_subset,
  features = c("CLEC9A", "XCR1", "BATF3", "IRF8", "CD1C", "FCER1A",
               "CLEC10A", "LILRA4", "CLEC4C", "IL3RA", "TCF4", "CCR7"),
  group.by = "myeloid_subcluster"
) +
  RotatedAxis() +
  ggtitle("Additional Human Dendritic Markers") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


# Neutrophils / granulocytes - ?
DotPlot(
  myeloid_subset,
  features = c(
    "CSF3R",
    "FCGR3B",
    "CEACAM8",
    "CXCR2",
    "MPO",
    "ELANE",
    "FPR1",
    "MMP8",
    "MMP9",
    "FUT4"
  ),
  group.by = "myeloid_subcluster"
) +
  RotatedAxis() +
  ggtitle("Neutrophil / Granulocyte Markers") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

DotPlot(
  myeloid_subset,
  features = c(
    "FCGR3B",
    "NAMPT",
    "IFITM2",
    "TNFAIP2",
    "PI3",
    "ITGAM",   # CD11b
    "FUT4",    # CD15
    "FCGR3A",  # CD16
    "IL7R"     # CD127
  ),
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  ggtitle("Neutrophil Marker DotPlot in Myeloid Subset") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Monocytes / inflammatory myeloid
VlnPlot(
  myeloid_subset,
  features = c("FCN1", "VCAN", "LYZ", "CCR2", "S100A8", "S100A9"),
  group.by = "myeloid_subcluster",
  pt.size = 0
)
FeaturePlot(
  myeloid_subset,
  reduction = "myeloid_umap",
  features = c("FCN1", "VCAN", "LYZ", "CCR2", "S100A8", "S100A9"),
  cols = c("lightgrey", "blue")
)

Idents(myeloid_subset) <- "myeloid_subcluster"
myeloid_markers <- FindAllMarkers(
  myeloid_subset,
  only.pos = TRUE,
  min.pct = 0.15,
  logfc.threshold = 0.25
)
top20_myeloid <- myeloid_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20)
View(top20_myeloid)
output_path <- "/Users/sudharshiniram/Desktop/top20_myeloid_markers_by_cluster.csv"
write.csv(
  top20_myeloid,
  file = output_path,
  row.names = FALSE
)
file.exists(output_path)

# Myeloid / epithelial / fibroblast / acinar contamination check
DotPlot(
  myeloid_subset,
  features = c(
    "PTPRC", "LST1", "TYROBP", "FCER1G", "AIF1", "CD68",
    "EPCAM", "KRT19", "KRT8", "KRT7", "CFTR", "MUC1",
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL11A1",
    "PRSS1", "PRSS2", "PNLIP", "SPINK1"
  ),
  group.by = "myeloid_subcluster"
) +
  RotatedAxis()

Idents(myeloid_subset) <- "myeloid_subcluster"
myeloid_subset_clean <- subset(
  myeloid_subset,
  idents = c("1", "2", "3")
)
table(myeloid_subset_clean$myeloid_subcluster)
DefaultAssay(myeloid_subset_clean) <- "RNA"
myeloid_subset_clean <- JoinLayers(myeloid_subset_clean)
counts_mat <- GetAssayData(
  myeloid_subset_clean,
  assay = "RNA",
  layer = "counts"
)
myeloid_subset_clean[["RNA"]] <- CreateAssay5Object(counts = counts_mat)
DefaultAssay(myeloid_subset_clean) <- "RNA"
myeloid_subset_clean <- NormalizeData(myeloid_subset_clean)
myeloid_subset_clean <- FindVariableFeatures(
  myeloid_subset_clean,
  selection.method = "vst",
  nfeatures = 3000
)
myeloid_subset_clean <- ScaleData(myeloid_subset_clean)
myeloid_subset_clean <- RunPCA(
  myeloid_subset_clean,
  features = VariableFeatures(myeloid_subset_clean)
)
myeloid_subset_clean <- harmony::RunHarmony(
  object = myeloid_subset_clean,
  group.by.vars = "orig.ident",
  reduction.use = "pca",
  dims.use = 1:20,
  reduction.save = "harmony"
)
myeloid_subset_clean <- FindNeighbors(
  myeloid_subset_clean,
  reduction = "harmony",
  dims = 1:20
)
myeloid_subset_clean <- FindClusters(
  myeloid_subset_clean,
  resolution = 0.05
)
myeloid_subset_clean$myeloid_subcluster_clean <- Idents(myeloid_subset_clean)
myeloid_subset_clean <- RunUMAP(
  myeloid_subset_clean,
  reduction = "harmony",
  dims = 1:20,
  reduction.name = "myeloid_umap_clean"
)
DimPlot(
  myeloid_subset_clean,
  reduction = "myeloid_umap_clean",
  group.by = "myeloid_subcluster_clean",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Cleaned Harmony-Reintegrated Myeloid Subclusters") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Cleaned myeloid subtype markers
DotPlot(
  myeloid_subset_clean,
  features = c(
    "PTPRC", "LST1", "TYROBP", "FCER1G", "AIF1", "CD68",
    "C1QA", "C1QB", "APOE", "SPP1",
    "FCN1", "VCAN", "CSF3R",
    "CLEC9A", "XCR1", "BATF3", "IRF8"
  ),
  group.by = "myeloid_subcluster_clean"
) +
  RotatedAxis()

# Cleaned myeloid / epithelial / fibroblast / acinar contamination check
DotPlot(
  myeloid_subset_clean,
  features = c(
    "PTPRC", "LST1", "TYROBP", "FCER1G", "AIF1", "CD68",
    "EPCAM", "KRT19", "KRT8", "KRT7", "CFTR", "MUC1",
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL11A1",
    "PRSS1", "PRSS2", "PNLIP", "SPINK1"
  ),
  group.by = "myeloid_subcluster_clean"
) +
  RotatedAxis()

# Macrophages / TAMs
macrophage_markers <- c(
  "C1QA", "C1QB", "C1QC", "APOE", "APOC1",
  "MSR1", "MRC1", "CD163", "SPP1", "GPNMB", "TREM2"
)
macrophage_markers <- macrophage_markers[
  macrophage_markers %in% rownames(myeloid_subset)
]
DotPlot(
  myeloid_subset_clean,
  features = macrophage_markers,
  group.by = "myeloid_subcluster"
) +
  RotatedAxis() +
  ggtitle("Macrophage / TAM Markers") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#Dendritic
DotPlot(
  myeloid_subset_clean,
  features = c("CLEC9A", "XCR1", "BATF3", "IRF8", "CD1C", "FCER1A",
               "CLEC10A", "LILRA4", "CLEC4C", "IL3RA", "TCF4", "CCR7"),
  group.by = "myeloid_subcluster"
) +
  RotatedAxis() +
  ggtitle("Additional Human Dendritic Markers") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Neutrophils / granulocytes
DotPlot(
  myeloid_subset,
  features = c(
    "CSF3R",
    "FCGR3B",
    "CEACAM8",
    "CXCR2",
    "MPO",
    "ELANE",
    "FPR1",
    "MMP8",
    "MMP9",
    "FUT4"
  ),
  group.by = "myeloid_subcluster"
) +
  RotatedAxis() +
  ggtitle("Neutrophil / Granulocyte Markers") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))



### TESTING IMMUNE CELLS (B & T CELLS) ####

#B Cells - Cluster 3
FeaturePlot(Merged_Labeled, features=c("CD19", "PTPRC", "CD27", "IGHD", "TNFRSF17", "CXCR4", "SDC1", "SLAMF7"), cols = c('lightgrey', 'blue'))

# T cells - Cluster 8
FeaturePlot(Merged_Labeled, features=c("CD27", "CD8A", "SELL", "CD44", "CD3E", "CD4"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged_Labeled, features=c("CD8A", "CD3E", "CD4"), cols = c('lightgrey', 'blue'))

### TESTING IMMUNE CELLS: SEPARATE B CELL AND T CELL SUBSETS ####

DefaultAssay(Merged_Joined) <- "RNA"
Idents(Merged_Joined) <- "seurat_clusters"

# B cells
bcell_subset <- subset(
  Merged_Joined,
  idents = "3"
)

bcell_markers <- c(
  "JCHAIN", "MZB1", "CD38", "BLK", "CD22", "MS4A1"
)

bcell_markers_present <- bcell_markers[
  bcell_markers %in% rownames(bcell_subset)
]

print(bcell_markers_present)

DotPlot(
  bcell_subset,
  features = bcell_markers_present,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  ggtitle("B Cell Marker Expression - Cluster 3") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

saveRDS(bcell_subset, file = '/Users/sudharshiniram/Desktop/SRI Data Analysis/BCell_Subset.RDS')


# T cells
tcell_subset <- subset(
  Merged_Joined,
  idents = "8"
)

tcell_markers <- c(
  "CD27", "CD8A", "SELL",
  "CD44", "CD3E", "CD4", "IL7R"
)

tcell_markers_present <- tcell_markers[
  tcell_markers %in% rownames(tcell_subset)
]

print(tcell_markers_present)

DotPlot(
  tcell_subset,
  features = tcell_markers_present,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  ggtitle("T Cell Marker Expression - Cluster 8") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

saveRDS(tcell_subset, file = '/Users/sudharshiniram/Desktop/SRI Data Analysis/TCell_Subset.RDS')
