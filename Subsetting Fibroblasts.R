library(Seurat)
library(ggplot2)

rm(Merged.markers, top20_genes, fibroblast_subset, Merged, Merged_Integrated, Merged_Joined, Merged_small, catch_obj, expr_matrix, catch_cluster_ids, Merged_sccatch)
gc()

#Testing Feature Plots
Merged <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Joined.RDS")
DimPlot(Merged, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("UMAP of Merged Dataset")


# Fibroblasts / Cancer-Associated Fibroblasts - CAFs
# Broad stromal / fibroblast markers.
# ------------------------------------------------------------
FeaturePlot(Merged, features=c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "SPARC"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("FN1", "ACTA2", "PDGFRA", "FAP", "THY1", "MMP2"), cols = c('lightgrey', 'blue'))

# Subsetting Fibroblasts
fibroblast_subset <- subset(Merged, idents = "1") 

# 1. Clear old reductions
fibroblast_subset@reductions$pca <- NULL
fibroblast_subset@reductions$harmony <- NULL
fibroblast_subset@reductions$umap <- NULL

# 2. Split layers by original sample so Harmony can detect batches
fibroblast_subset[["RNA"]] <- split(fibroblast_subset[["RNA"]], f = fibroblast_subset$orig.ident)

# 3. Standard re-processing
fibroblast_subset <- FindVariableFeatures(fibroblast_subset, selection.method = "vst", nfeatures = 2000)
fibroblast_subset <- ScaleData(fibroblast_subset)
fibroblast_subset <- RunPCA(fibroblast_subset, features = VariableFeatures(object = fibroblast_subset))

# 4. Re-integrate with Harmony
fibroblast_subset <- IntegrateLayers(
  object = fibroblast_subset, 
  method = HarmonyIntegration, 
  orig.reduction = "pca",
  new.reduction = "harmony", 
  verbose = FALSE
)

# 5. Re-join layers so FeaturePlot and DotPlot function correctly
fibroblast_subset <- JoinLayers(fibroblast_subset)

# 6. Cluster and UMAP using the new Harmony reduction
fibroblast_subset <- FindNeighbors(fibroblast_subset, reduction = "harmony", dims = 1:15)
fibroblast_subset <- FindClusters(fibroblast_subset, resolution = 0.2)
fibroblast_subset <- RunUMAP(fibroblast_subset, reduction = "harmony", dims = 1:15, reduction.name = "sub_umap")

DimPlot(fibroblast_subset, reduction = "sub_umap", label = TRUE, repel = TRUE) + 
  ggtitle("Fibroblast Sub-Clusters (Cluster 1)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#cleaning out immune cells
fibroblast_subset <- subset(fibroblast_subset, idents = c("3", "4"), invert = TRUE)
fibroblast_subset <- subset(
  fibroblast_subset, 
  subset = seurat_clusters != "2" | (seurat_clusters == "2" & PTPRC < 1)
)
fibroblast_subset@reductions$pca <- NULL
fibroblast_subset@reductions$harmony <- NULL
fibroblast_subset@reductions$umap <- NULL

#checking fibroblast subset has no immune cells
VlnPlot(fibroblast_subset, 
        features = c("PTPRC", "EPCAM", "ACTA2", "PDGFRA", "DCN"), 
        pt.size = 0, # Hides the individual dots to clearly show the bulk distribution
        stack = TRUE, 
        flip = TRUE)


# myCAFs
FeaturePlot(fibroblast_subset, features=c("ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2", "POSTN"), reduction = "sub_umap", cols = c('lightgrey', 'blue'))

# iCAFs - inflammatory CAFs
FeaturePlot(fibroblast_subset, features=c("IL6", "CXCL12", "CXCL14", "LIF", "CFD", "PDGFRA"), reduction = "sub_umap", cols = c('lightgrey', 'blue'))

# apCAFs - antigen-presenting CAFs
FeaturePlot(fibroblast_subset, features=c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA"), reduction = "sub_umap", cols = c('lightgrey', 'blue'))

# ECM-remodeling / matrix CAFs (FIXED: Switched from Merged to fibroblast_subset)
FeaturePlot(fibroblast_subset, features=c("COL11A1", "COL10A1", "POSTN", "MMP11", "THBS2", "INHBA"), reduction = "sub_umap", cols = c('lightgrey', 'blue'))

# Group by your new sub-clusters to see exactly which cluster number lights up
DotPlot(fibroblast_subset, 
        features = c("IL6", "CXCL12", "CXCL14", "LIF", "CFD", "PDGFRA"), 
        group.by = "seurat_clusters") + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("iCAF Marker Expression by Sub-cluster")

DotPlot(fibroblast_subset, 
        features = c("ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2", "POSTN"), 
        group.by = "seurat_clusters") + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("myCAF Marker Expression by Sub-cluster")

DotPlot(fibroblast_subset, 
        features = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA"), 
        group.by = "seurat_clusters") + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("apCAF Marker Expression by Sub-cluster")





# General fibroblast / CAF markers
FeaturePlot(
  fibroblast_subset,
  features = c("LIF", "CD34", "CD248", "CXCL12", "IL6", "ACTA2", "FGF18", "BCAM"),
  cols = c("lightgrey", "blue")
)

# apCAFs / antigen-presenting CAFs
FeaturePlot(
  fibroblast_subset_clean,
  features = c("CD74", "PDPN", "FAP", "C7", "PDGFRB"), 
  cols = c("lightgrey", "blue")
)

FeaturePlot(
  fibroblast_subset,
  features = c("RGS5", "NOTCH3", "ADIRF", "CRIP1"),
  cols = c("lightgrey", "blue")
)

# iCAFs / inflammatory CAFs
FeaturePlot(
  fibroblast_subset_clean,
  features = c("CXCL12", "IL6", "C3", "PLA2G2A", "STAT1", "IRF1", "IRF9"),
  cols = c("lightgrey", "blue")
)

# myCAFs / myofibroblastic CAFs
FeaturePlot(
  fibroblast_subset_clean,
  features = c("TAGLN", "MMP11", "HOPX", "PDGFRA", "ACTA2", "POSTN"),
  cols = c("lightgrey", "blue")
)

# matrix CAFs / ECM-remodeling CAFs
FeaturePlot(
  fibroblast_subset,
  features = c("MMP11", "COL1A1", "COL3A1", "POSTN", "CDH11"),
  cols = c("lightgrey", "blue")
)

# senCAFs / senescence-like CAFs
FeaturePlot(
  fibroblast_subset,
  features = c("COL12A1", "TAGLN", "COL8A1", "C3", "PDGFRA", "HAS1", "PI16", "DPT", "CD74"),
  cols = c("lightgrey", "blue")
)

# Normal-like fibroblast / fibroblast-associated markers
FeaturePlot(
  fibroblast_subset,
  features = c("ENG", "PDPN", "FAP"),
  cols = c("lightgrey", "blue")
)


# 1. Find all markers for the new fibroblast sub-clusters
# (This step takes a moment to compute)
fibroblast.markers <- FindAllMarkers(fibroblast_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 2. Extract the top 10 significant genes per cluster
top10_fibro_genes <- fibroblast.markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%                  # Keep only statistically significant genes
  slice_max(n = 10, order_by = avg_log2FC)      # Grab the top 10 based on highest expression fold change

# 3. Save this clean list to a new CSV
write.csv(top10_fibro_genes, 
          file = "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Top10_Markers.csv", 
          row.names = FALSE)


print(top10_fibro_genes %>% slice_max(n = 3, order_by = avg_log2FC))