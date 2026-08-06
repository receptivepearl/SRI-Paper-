library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(Azimuth)

Sys.setenv('R_MAX_VSIZE'=32000000000)

Sample_01 <- CreateSeuratObject(counts=Read10X('/Users/sudharshiniram/Desktop/SRI Data/Sample_01'), project="Sample_01", min.cells=3, min.features=200)
Sample_02 <- CreateSeuratObject(counts=Read10X('/Users/sudharshiniram/Desktop/SRI Data/Sample_02'), project="Sample_02", min.cells=3, min.features=200)
Sample_03 <- CreateSeuratObject(counts=Read10X('/Users/sudharshiniram/Desktop/SRI Data/Sample_03'), project="Sample_03", min.cells=3, min.features=200)
Sample_04 <- CreateSeuratObject(counts=Read10X('/Users/sudharshiniram/Desktop/SRI Data/Sample_04'), project="Sample_04", min.cells=3, min.features=200)
Sample_05 <- CreateSeuratObject(counts=Read10X('/Users/sudharshiniram/Desktop/SRI Data/Sample_05'), project="Sample_05", min.cells=3, min.features=200)
Sample_06 <- CreateSeuratObject(counts=Read10X('/Users/sudharshiniram/Desktop/SRI Data/Sample_06'), project="Sample_06", min.cells=3, min.features=200)
Sample_07 <- CreateSeuratObject(counts=Read10X('/Users/sudharshiniram/Desktop/SRI Data/Sample_07'), project="Sample_07", min.cells=3, min.features=200)
Sample_08 <- CreateSeuratObject(counts=Read10X('/Users/sudharshiniram/Desktop/SRI Data/Sample_08'), project="Sample_08", min.cells=3, min.features=200)
Sample_09 <- CreateSeuratObject(counts=Read10X('/Users/sudharshiniram/Desktop/SRI Data/Sample_09'), project="Sample_09", min.cells=3, min.features=200)

gc() 

Sample_10 <- CreateSeuratObject(counts=Read10X('/Users/sudharshiniram/Desktop/SRI Data/Sample_10'), project="Sample_10", min.cells=3, min.features=200)
Sample_11 <- CreateSeuratObject(counts=Read10X('/Users/sudharshiniram/Desktop/SRI Data/Sample_11'), project="Sample_11", min.cells=3, min.features=200)
Sample_12 <- CreateSeuratObject(counts=Read10X('/Users/sudharshiniram/Desktop/SRI Data/Sample_12'), project="Sample_12", min.cells=3, min.features=200)
Sample_13 <- CreateSeuratObject(counts=Read10X('/Users/sudharshiniram/Desktop/SRI Data/Sample_13'), project="Sample_13", min.cells=3, min.features=200)
Sample_14 <- CreateSeuratObject(counts=Read10X('/Users/sudharshiniram/Desktop/SRI Data/Sample_14'), project="Sample_14", min.cells=3, min.features=200)
Sample_15 <- CreateSeuratObject(counts=Read10X('/Users/sudharshiniram/Desktop/SRI Data/Sample_15'), project="Sample_15", min.cells=3, min.features=200)
Sample_16 <- CreateSeuratObject(counts=Read10X('/Users/sudharshiniram/Desktop/SRI Data/Sample_16'), project="Sample_16", min.cells=3, min.features=200)
Sample_17 <- CreateSeuratObject(counts=Read10X('/Users/sudharshiniram/Desktop/SRI Data/Sample_17'), project="Sample_17", min.cells=3, min.features=200)

gc()

Merged <- merge(x = Sample_01, 
                y = c(Sample_02, Sample_03, Sample_04, Sample_05, Sample_06, 
                      Sample_07, Sample_08, Sample_09, Sample_10, Sample_11, 
                      Sample_12, Sample_13, Sample_14, Sample_15, Sample_16, Sample_17), 
                add.cell.ids = c("S01", "S02", "S03", "S04", "S05", "S06", 
                                 "S07", "S08", "S09", "S10", "S11", "S12", 
                                 "S13", "S14", "S15", "S16", "S17"))

rm(Sample_01, Sample_02, Sample_03, Sample_04, Sample_05, Sample_06, Sample_07, Sample_08, Sample_09, Sample_10, Sample_11, Sample_12, Sample_13, Sample_14, Sample_15, Sample_16, Sample_17)
gc()
saveRDS(Merged, file = '/Users/sudharshiniram/Desktop/SRI Data Analysis/Preprocessed Merge.rds')


# QC and Filtering
Merged <- PercentageFeatureSet(Merged, pattern = "^MT-", col.name = "percent.mt")
Merged <- subset(Merged, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5 & nCount_RNA > 500 & nCount_RNA < 10000)


# Normalization & Dimensionality Reduction
Merged <- NormalizeData(Merged, normalization.method = "LogNormalize", scale.factor = 10000)
Merged <- FindVariableFeatures(Merged, selection.method = "vst", nfeatures = 2000)

gc()
Merged <- ScaleData(Merged)
Merged <- RunPCA(Merged, features = VariableFeatures(object = Merged))

#Integration
install.packages("harmony")

gc()
Merged <- ScaleData(Merged)
Merged <- RunPCA(Merged, features = VariableFeatures(object = Merged))
gc()
Merged <- IntegrateLayers(
  object = Merged, 
  method = HarmonyIntegration, 
  orig.reduction = "pca", 
  new.reduction = "harmony", 
  verbose = TRUE
)
Merged <- FindNeighbors(Merged, reduction = "harmony", dims = 1:20)
Merged <- FindClusters(Merged, resolution = 0.5)
Merged <- RunUMAP(Merged, reduction = "harmony", dims = 1:20)
DimPlot(Merged, reduction = "umap", label = TRUE, repel = TRUE) + 
  ggtitle("Fully Integrated PDAC Data (Harmony)")
Merged <- JoinLayers(Merged)
saveRDS(Merged, file = '/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Integrated.rds')
#FeaturePlot(Merged, features=c("PTPRC"), cols = c('lightgray', 'blue'))

'''
#scCatch
install.packages("scCATCH")
library(scCATCH)
Merged_small <- subset(Merged, downsample = 300)
expr_matrix_small <- LayerData(Merged_small, assay = "RNA", layer = "data")
catch_obj <- createscCATCH(data = expr_matrix_small, cluster = as.character(Idents(Merged_small)))
catch_obj <- findmarkergene(
  object = catch_obj, 
  species = "Human", 
  marker = cellmatch, 
  tissue = c("Pancreas", "Pancreatic acinar tissue", "Pancreatic islet", "Blood vessel", "Fetal pancreas", "Lymphoid tissue"), 
  cancer = "Normal"
)
catch_obj <- findcelltype(catch_obj)
catch_cluster_ids <- setNames(catch_obj@celltype$cell_type, catch_obj@celltype$cluster)
Merged <- RenameIdents(Merged, catch_cluster_ids)
Merged$scCATCH_annotation <- Idents(Merged)
rm(Merged_small, expr_matrix_small, catch_obj)
gc()
scCATCH_umap <- DimPlot(Merged, reduction = "umap", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("PDAC Tumor Microenvironment by scCATCH") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(scCATCH_umap)
ggsave(filename = "/Users/sudharshiniram/Desktop/SRI Data Analysis/scCATCH_UMAP.png", plot = scCATCH_umap, width = 14, height = 10, dpi = 300)
saveRDS(Merged, file = '/Users/sudharshiniram/Desktop/SRI Data Analysis/Data_Analysis_Final_scCATCH_Labeled.rds')
'''



#Azimuth Package
library(Azimuth)
library(Seurat)
library(ggplot2)
library(patchwork)
remotes::install_version("TFMPvalue", version = "0.0.9", upgrade = "never")
BiocManager::install("TFBSTools", ask = FALSE)
remotes::install_github("satijalab/azimuth", ref = "master", upgrade = "never")
install.packages("rappdirs")
Merged <- RunAzimuth(query = Merged, reference = "pancreasref")
head(Merged@meta.data)
plot_types <- DimPlot(Merged, group.by = "predicted.annotation.l1", reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("Azimuth Predicted Cell Types") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
plot_scores <- FeaturePlot(Merged, features = "mapping.score", reduction = "umap") +
  ggtitle("Mapping Confidence (Dark = Healthy, Light = Malignant)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
final_combined_plot <- plot_types + plot_scores
ggsave(filename = "/Users/sudharshiniram/Desktop/SRI Data Analysis/Azimuth_Tumor_Map.png", 
       plot = final_combined_plot, 
       width = 16,  
       height = 8,  
       dpi = 300) 
azimuth_types <- DimPlot(Merged, group.by = "predicted.annotation.l1", reduction = "umap", label = TRUE, repel = TRUE) + 
  ggtitle("Azimuth Predicted Cell Types") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
azimuth_scores <- FeaturePlot(Merged, features = "mapping.score", reduction = "umap") + 
  ggtitle("Mapping Confidence (Dark = Healthy, Light = Malignant)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
final_azimuth_plot <- azimuth_types + azimuth_scores
print(final_azimuth_plot)

# Calculate Marker Genes
'''
gc()

Merged.markers <- FindAllMarkers(Merged, only.pos = TRUE, max.cells.per.ident = 500)

#saveRDS(Merged.markers, file = "/Users/sudharshiniram/Desktop/SRI Data Analysis/Marker_Genes.rds")
write.csv(Merged.markers, file = "/Users/sudharshiniram/Desktop/SRI Data Analysis/Marker_Genes.csv")


# Label Clusters
my_cell_types <- c(
  "Epithelial / PDAC",        # 0
  "myCAF",                    # 1
  "B-Cell",                   # 2
  "Macrophage",               # 3
  "Epithelial",               # 4
  "Unknown",                  # 5
  "myCAF",                    # 6
  "Unknown",                  # 7
  "Unknown",                  # 8
  "myCAF",                    # 9
  "Macrophage / Dendritic",   # 10
  "Endocrine (Islet)",        # 11
  "Unknown",                  # 12
  "Endothelial",              # 13
  "Acinar",                   # 14
  "Unknown",                  # 15
  "Epithelial (Metaplastic)", # 16
  "Unknown",                  # 17
  "Unknown",                  # 18
  "Acinar",                   # 19
  "Acinar",                   # 20
  "Epithelial / PDAC",        # 21
  "iCAF",                     # 22
  "Endocrine (Islet)",        # 23
  "CAF",                      # 24
  "Acinar",                   # 25
  "Macrophage (M2/TAM)",      # 26
  "Epithelial / Ductal",      # 27
  "Unknown",                  # 28
  "T-Cell / NK Cell",         # 29
  "CD8+ Cytotoxic T-Cell",    # 30
  "Endothelial",              # 31
  "NK Cell",                  # 32
  "Monocyte / Myeloid",       # 33
  "Unknown",                  # 34
  "Endothelial",              # 35
  "iCAF",                     # 36
  "Neural / Schwann Cell",    # 37
  "Lymphatic Endothelial",    # 38
  "Endothelial",              # 39
  "Unknown",                  # 40
  "Epithelial",               # 41
  "Endothelial",              # 42
  "Endocrine (Islet)"         # 43
)

names(my_cell_types) <- levels(Merged)
Merged <- RenameIdents(Merged, my_cell_types)


# Graphing
final_umap <- DimPlot(Merged, reduction = "umap", label = TRUE, label.size = 3.5, repel = TRUE) + 
  ggtitle("PDAC Tumor Microenvironment by Cell Type") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(filename = "/Users/sudharshiniram/Desktop/SRI Data Analysis/Final_Labeled_UMAP.png", 
       plot = final_umap, 
       width = 14,      
       height = 10,     
       dpi = 300)  

saveRDS(Merged, file = '/Users/sudharshiniram/Desktop/SRI Data Analysis/Data_Analysis_Final_Labeled.rds')

# Neurons
# Note: 'Tuj1' is an antibody target alias for TUBB3, 'vChat'/'Vacht' are aliases for SLC18A3, and 'Brn3a' is POU4F1.
FeaturePlot(Merged, features=c("PRPH", "CALCA", "CHAT", "TUBB3", "TH", "NEFL", "MAP2", "MAPT"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("CHAT", "TRPV1", "SLC18A3", "POU4F1"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("TUBB3", "UCHL1", "RBFOX3", "SYP"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("TUBB3"), cols = c('lightgrey', 'blue'))

# Macrophages
# Note: Humans do not have Cd209f, the closest is CD209. 
FeaturePlot(Merged, features=c("CD209", "CCR2", "MRC1", "CD163"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("CX3CR1", "CD80", "CD86", "NOS2", "CD38", "ARG1", "MRC1", "CCR2", "MERTK"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("CD68", "SIGLEC1", "ADGRE1", "CD163"), cols = c('lightgrey', 'blue'))

# Neutrophils
# Note: Humans lack Ly6g/Ly6c1. Replaced with standard human neutrophil markers CEACAM8 (CD66b) and FCGR3B (CD16b).
FeaturePlot(Merged, features=c("CEACAM8", "FCGR3B", "CXCR4", "CD33", "FUT4", "CD14"), cols = c('lightgrey', 'blue'))

# Monocytes
# Note: Replaced Ly6g/Ly6c2 with human monocyte markers CD14 and FCGR3A (CD16a).
FeaturePlot(Merged, features=c("CD14", "FCGR3A", "CXCR4", "CD33"), cols = c('lightgrey', 'blue'))

# Schwann cells
FeaturePlot(Merged, features=c("MBP", "S100B", "MPZ", "NCAM1", "JUN", "NRG1", "SCN7A", "ERBB3", "SOX2"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("S100B", "GFAP", "NGFR", "MBP", "POU3F1", "EGR2"), cols = c('lightgrey', 'blue'))

# Neuroendocrine cells
FeaturePlot(Merged, features=c("CALCB", "SCG5", "SCG2", "CHGB", "PCSK1N", "CHGA", "SCG3"), cols = c('lightgrey', 'blue'))

# Acinar cells
# Note: If this is pancreatic tissue, AMY2A is usually the prevalent amylase over AMY1A.
FeaturePlot(Merged, features=c("GP2", "GATA4", "AMY2A", "AMY1A", "PRSS1"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("LIPA", "LPL", "LIPE", "PNLIP", "PNPLA2"), cols = c('lightgrey', 'blue'))

# Duct cells
FeaturePlot(Merged, features=c("KRT19", "SOX9", "PDX1", "P2RY1", "BMPR1A"), cols = c('lightgrey', 'blue'))

# VlnPlot(Merged, features = c("NINJ1"))

# apCAFs
# Note: Hladrb1 formatted to the proper human hyphenated symbol HLA-DRB1.
FeaturePlot(Merged, features=c("CD74", "ACTA2", "HLA-DRB1"), cols = c('lightgrey', 'blue'))

# Mesothelial cells
FeaturePlot(Merged, features=c("CALB2", "KRT5", "KRT6A", "CD44", "THBD"), cols = c('lightgrey', 'blue'))

# B cells
FeaturePlot(Merged, features=c("CD19", "PTPRC", "CD27", "IGHD", "TNFRSF17", "CXCR4", "SDC1", "SLAMF7"), cols = c('lightgrey', 'blue'))

# T cells
FeaturePlot(Merged, features=c("CD27", "CD8A", "SELL", "CD44", "CD3E", "CD4"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("CD8A", "CD3E", "CD4"), cols = c('lightgrey', 'blue'))

# Endocrine cells
# Note: Humans only have one insulin gene (INS), not Ins1/Ins2.
FeaturePlot(Merged, features=c("TFRC", "GCG", "SST", "INS"), cols = c('lightgrey', 'blue'))

# Eosinophils
# Note: The functional human equivalent of mouse Siglecf is SIGLEC8.
FeaturePlot(Merged, features=c("CCR3", "ITGA2", "SIGLEC8", "CEACAM8"), cols = c('lightgrey', 'blue'))

# Basophils
# Note: Crth2 is an alias; the official HGNC symbol is PTGDR2. I included both just in case your reference genome uses the alias.
FeaturePlot(Merged, features=c("IL3RA", "ENPP3", "PTGDR2", "CRTH2", "CD63"), cols = c('lightgrey', 'blue'))

'''


