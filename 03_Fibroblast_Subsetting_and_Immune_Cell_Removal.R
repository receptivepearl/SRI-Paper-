library(Seurat)
library(ggplot2)
library(dplyr)

rm(fibroblast_subset_clean, Merged.markers, top20_genes, fibroblast_subset, Merged, Merged_Integrated, Merged_Joined, Merged_small, catch_obj, expr_matrix, catch_cluster_ids, Merged_sccatch)
gc()
# ------------------------------------------------------------
# INITIALIZATION
# ------------------------------------------------------------
Merged <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Joined.RDS")

# Subsetting Fibroblasts (Right now, everything is just "Cluster 1")
fibroblast_subset <- subset(Merged, idents = "1") 

# ------------------------------------------------------------
# FIRST PASS: Generate Temporary Sub-Clusters
# ------------------------------------------------------------
fibroblast_subset[["RNA"]] <- split(fibroblast_subset[["RNA"]], f = fibroblast_subset$orig.ident)
fibroblast_subset <- FindVariableFeatures(fibroblast_subset, selection.method = "vst", nfeatures = 2000)
fibroblast_subset <- ScaleData(fibroblast_subset)
fibroblast_subset <- RunPCA(fibroblast_subset, features = VariableFeatures(object = fibroblast_subset))

fibroblast_subset <- IntegrateLayers(
  object = fibroblast_subset, 
  method = HarmonyIntegration, 
  orig.reduction = "pca",
  new.reduction = "harmony", 
  verbose = FALSE
)

fibroblast_subset <- JoinLayers(fibroblast_subset)

fibroblast_subset <- FindNeighbors(fibroblast_subset, reduction = "harmony", dims = 1:15)
fibroblast_subset <- FindClusters(fibroblast_subset, resolution = 0.2)

# THE FIX: You must run RunUMAP before you can plot a UMAP!
fibroblast_subset <- RunUMAP(fibroblast_subset, reduction = "harmony", dims = 1:15, reduction.name = "sub_umap")

DimPlot(fibroblast_subset, reduction = "sub_umap", label = TRUE, repel = TRUE) + 
  ggtitle("Fibroblast Sub-Clusters") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

VlnPlot(fibroblast_subset, 
        features = c("PTPRC", "EPCAM", "ACTA2", "PDGFRA", "DCN"), 
        pt.size = 0, 
        stack = TRUE, 
        flip = TRUE)

saveRDS(fibroblast_subset, file = '/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS')
fibroblast_subset <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS")

# ------------------------------------------------------------
# SURGICAL CLEANUP: Global PTPRC Removal
# ------------------------------------------------------------
# THE FIX: Added 'subset = PTPRC < 1' back in to guarantee no immune cells survive
fibroblast_subset_clean <- subset(
  fibroblast_subset, 
  invert = TRUE,
  subset = PTPRC < 1
)

# Clear the old, dirty reductions
fibroblast_subset_clean@reductions$pca <- NULL
fibroblast_subset_clean@reductions$harmony <- NULL
fibroblast_subset_clean@reductions$umap <- NULL
# ------------------------------------------------------------
# SECOND PASS: Build the Final, Clean Map
# ------------------------------------------------------------
fibroblast_subset_clean[["RNA"]] <- split(fibroblast_subset_clean[["RNA"]], f = fibroblast_subset_clean$orig.ident)

fibroblast_subset_clean <- FindVariableFeatures(fibroblast_subset_clean, selection.method = "vst", nfeatures = 2000)
fibroblast_subset_clean <- ScaleData(fibroblast_subset_clean)
fibroblast_subset_clean <- RunPCA(fibroblast_subset_clean, features = VariableFeatures(object = fibroblast_subset_clean))

fibroblast_subset_clean <- IntegrateLayers(
  object = fibroblast_subset_clean, 
  method = HarmonyIntegration, 
  orig.reduction = "pca",
  new.reduction = "harmony", 
  verbose = FALSE
)

fibroblast_subset_clean <- JoinLayers(fibroblast_subset_clean)

# Final Clustering and UMAP
fibroblast_subset_clean <- FindNeighbors(fibroblast_subset_clean, reduction = "harmony", dims = 1:20)
fibroblast_subset_clean <- FindClusters(fibroblast_subset_clean, resolution = 0.4)
fibroblast_subset_clean <- RunUMAP(fibroblast_subset_clean, reduction = "harmony", dims = 1:20, reduction.name = "sub_umap")

saveRDS(fibroblast_subset_clean, file = '/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset_Clean.RDS')
fibroblast_subset_clean <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset_Clean.RDS")


# ------------------------------------------------------------
# PLOTTING & VERIFICATION
# ------------------------------------------------------------
DimPlot(fibroblast_subset_clean, reduction = "sub_umap", label = TRUE, repel = TRUE) + 
  ggtitle("CLEAN Fibroblast Sub-Clusters") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Verification: Ensure immune cells are gone
VlnPlot(fibroblast_subset_clean, 
        features = c("PTPRC", "EPCAM", "ACTA2", "PDGFRA", "DCN"), 
        pt.size = 0, 
        stack = TRUE, 
        flip = TRUE)

# myCAFs
FeaturePlot(fibroblast_subset_clean, features=c("ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2", "POSTN"), reduction = "sub_umap", cols = c('lightgrey', 'blue'))

# iCAFs
FeaturePlot(fibroblast_subset_clean, features=c("IL6", "CXCL12", "CXCL14", "LIF", "CFD", "PDGFRA"), reduction = "sub_umap", cols = c('lightgrey', 'blue'))

# apCAFs
FeaturePlot(fibroblast_subset_clean, features=c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA"), reduction = "sub_umap", cols = c('lightgrey', 'blue'))

#iCAFs
DotPlot(fibroblast_subset_clean, features = c("IL6", "CXCL12", "CXCL14", "LIF", "CFD", "PDGFRA"), group.by = "seurat_clusters") + 
  coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("iCAF Marker Expression")

#myCAFs
DotPlot(fibroblast_subset_clean, features = c("ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2", "POSTN"), group.by = "seurat_clusters") + 
  coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("myCAF Marker Expression")

#apCAFs
DotPlot(fibroblast_subset_clean, features = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA"), group.by = "seurat_clusters") + 
  coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("myCAF Marker Expression")


# ------------------------------------------------------------
# FINDING SUB-CLUSTER MARKERS
# ------------------------------------------------------------
fibroblast.markers <- FindAllMarkers(fibroblast_subset_clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10_fibro_genes <- fibroblast.markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%                  
  slice_max(n = 10, order_by = avg_log2FC)      

write.csv(top10_fibro_genes, file = "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Top10_Markers.csv", row.names = FALSE)

print(top10_fibro_genes %>% slice_max(n = 3, order_by = avg_log2FC))
