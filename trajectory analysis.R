library(Seurat)
library(ggplot2)
library(harmony) 

Merged <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Joined.RDS")

fibroblast_subset <- subset(Merged, idents = "1") 

fibroblast_subset@reductions$pca <- NULL
fibroblast_subset@reductions$harmony <- NULL
fibroblast_subset@reductions$umap <- NULL

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
fibroblast_subset <- RunUMAP(fibroblast_subset, reduction = "harmony", dims = 1:15, reduction.name = "sub_umap")

DimPlot(fibroblast_subset, reduction = "sub_umap", label = TRUE, repel = TRUE) + 
  ggtitle("Fibroblast Sub-Clusters (Integrated)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

saveRDS(fibroblast_subset, file = '/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS')


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

rm(list = ls())
gc()
fibroblast_subset <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS")


library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(ggplot2)
library(viridis) # Required for the plasma color scale


# ------------------------------------------------------------
# 2. CONVERT TO BIOCONDUCTOR FORMAT
# ------------------------------------------------------------
# Slingshot requires a SingleCellExperiment (SCE) object instead of Seurat
# This function safely ports over your counts, PCA, UMAP, and cluster labels
sce <- as.SingleCellExperiment(fibroblast_subset)


# ------------------------------------------------------------
# 3. RUN SLINGSHOT TRAJECTORY INFERENCE
# ------------------------------------------------------------
sce.sling <- slingshot(sce, clusterLabels = 'ident', reducedDim = 'PCA')

# Extract the calculated pseudotime values for each branching path
pseudo.paths <- slingPseudotime(sce.sling)

# Calculate the average pseudotime across all paths for each cell
# and add it directly into your Seurat object's metadata for easy plotting
fibroblast_subset$pseudotime <- rowMeans(pseudo.paths, na.rm = TRUE)


# ------------------------------------------------------------
# 4. VISUALIZE THE TRAJECTORY ON YOUR UMAP
# ------------------------------------------------------------
# Base Plot: Color the cells by their pseudotime progression (dark = early, bright = late)
p1 <- FeaturePlot(fibroblast_subset, features = "pseudotime", reduction = "sub_umap") +
  scale_color_viridis_c(option = "plasma") + 
  ggtitle("Fibroblast Trajectory (Pseudotime)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Project the PCA-based trajectory curves onto your 2D "SUB_UMAP" visual space
embedded <- embedCurves(sce.sling, "SUB_UMAP") 
embedded_curves <- slingCurves(embedded)

# Loop through each branching path and draw it as a thick black line
for (path in embedded_curves) {
  path_data <- as.data.frame(path$s[path$ord, ])
    colnames(path_data) <- c("x", "y")
    p1 <- p1 + geom_path(data = path_data, aes(x = x, y = y), linewidth = 1.2, color = "black", inherit.aes = FALSE)
}

print(p1)