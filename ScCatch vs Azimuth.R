install.packages("SeuratObject")
install.packages(c("SeuratObject", "Seurat"))
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

Merged <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Preprocessed Merge.rds")

# QC and Filtering
#Merged <- PercentageFeatureSet(Merged, pattern = "^MT-", col.name = "percent.mt")
#Merged <- subset(Merged, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 2 & nCount_RNA > 500 & nCount_RNA < 10000)


# QC and Filtering
# 1. Join the layers so Seurat can see all the counts at once
Merged <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Preprocessed Merge.rds")
Merged <- JoinLayers(Merged)
Merged <- PercentageFeatureSet(Merged, pattern = "^MT-", col.name = "percent.mt")
head(Merged@meta.data)
Merged <- subset(Merged, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 2 & nCount_RNA > 500 & nCount_RNA < 10000)

Merged[["RNA"]] <- split(Merged[["RNA"]], f = Merged$orig.ident)

# Normalization & Dimensionality Reduction
Merged <- NormalizeData(Merged, normalization.method = "LogNormalize", scale.factor = 10000)
Merged <- FindVariableFeatures(Merged, selection.method = "vst", nfeatures = 2000)

gc()
Merged <- ScaleData(Merged)

Merged <- RunPCA(Merged, features = VariableFeatures(object = Merged))

saveRDS(Merged, file = '/Users/sudharshiniram/Desktop/SRI Data Analysis/PostProcessing Merged.RDS')
Merged <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/PostProcessing Merged.RDS")

#Data integration

Merged_Integrated <- IntegrateLayers(
  object = Merged, 
  method = HarmonyIntegration, 
  orig.reduction = "pca",
  new.reduction = "harmony", 
  verbose = FALSE
)

#Find Neighbors/clusters and UMAP
Merged_Integrated <- FindNeighbors(Merged_Integrated, reduction = "harmony", dims = 1:25)
Merged_Integrated <- FindClusters(Merged_Integrated, resolution = 0.2)
Merged_Integrated <- RunUMAP(Merged_Integrated, reduction = "harmony", dims = 1:25)
DimPlot(Merged_Integrated, reduction = "umap")

saveRDS(Merged_Integrated, file ="/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Integrated.rds")
Merged_Integrated <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Integrated.rds")
Merged_Joined <- JoinLayers(Merged_Integrated)
saveRDS(Merged_Joined, file ="/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Joined.RDS")
rm(Merged,Merged_Integrated, Merged_Joined, Merged_small, catch_obj, expr_matrix, catch_cluster_ids)
gc()

# check to see everything is integrated and joined correctly
'''
Merged_Integrated <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Integrated.rds")
Merged_Joined <-readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Joined.RDS")
head(Merged_Integrated@meta.data)
Reductions(Merged_Integrated)
Layers(Merged_Joined)
#Merged_Joined[["RNA"]]$data[1:5, 1:5]
'''

#scCatch
Sys.setenv('R_MAX_VSIZE'=32000000000)
Merged_Integrated <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Integrated.rds")
Merged_Joined <- JoinLayers(Merged_Integrated)
rm(Merged_Integrated)
gc()
install.packages("scCATCH") 
library(scCATCH)
Merged_small <- subset(Merged_Joined, downsample = 1500)
expr_matrix <- LayerData(Merged_small, assay = "RNA", layer = "data")
catch_obj <- createscCATCH(data = expr_matrix, cluster = as.character(Idents(Merged_small)))
catch_obj <- findmarkergene(
  object = catch_obj, 
  species = "Human", 
  marker = cellmatch, 
  tissue = c("Pancreas", "Pancreatic acinar tissue", "Pancreatic islet", "Blood", "Peripheral blood", "Lymph node", "Lymphoid tissue"), 
  cancer = "Normal"
)

catch_obj <- findcelltype(catch_obj)
catch_cluster_ids <- setNames(catch_obj@celltype$cell_type, catch_obj@celltype$cluster)
Merged_sccatch <- RenameIdents(Merged_Joined, catch_cluster_ids)
Merged_Joined$scCATCH_annotation <- Idents(Merged_Joined)
rm(Merged_small, expr_matrix, catch_obj)
gc()

scCATCH_umap <- DimPlot(Merged_sccatch, reduction = "umap", label = TRUE, label.size = 3.5, repel = TRUE, raster = FALSE) + 
  ggtitle("PDAC Tumor Microenvironment by scCATCH") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(scCATCH_umap)
ggsave(filename = "/Users/sudharshiniram/Desktop/SRI Data Analysis/scCATCH_UMAP.png", plot = scCATCH_umap, width = 14, height = 10, dpi = 300)
saveRDS(Merged_Joined, file = '/Users/sudharshiniram/Desktop/SRI Data Analysis/Data_Analysis_Final_scCATCH_Labeled.rds')
graphics.off()

scCATCH_umap <- DimPlot(Merged_sccatch, reduction = "umap", label = TRUE, label.size = 3.5, repel = TRUE) + 
  ggtitle("PDAC Tumor Microenvironment by scCATCH") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(scCATCH_umap)
ggsave(filename = "/Users/sudharshiniram/Desktop/SRI Data Analysis/scCATCH_UMAP.png", 
       plot = scCATCH_umap, 
       width = 14, 
       height = 10, 
       dpi = 300)

rm(Merged_Integrated, Merged_Joined, Merged_small, catch_obj, expr_matrix, catch_cluster_ids, Merged_sccatch)
gc()



# Azimuth 
library(Azimuth)
library(Seurat)
library(ggplot2)
library(patchwork)
Merged_Joined <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Joined.RDS")
Merged_Joined <- RunAzimuth(query = Merged_Joined, reference = "pancreasref")
head(Merged_Joined@meta.data)
plot_types <- DimPlot(Merged_Joined, group.by = "predicted.annotation.l1", reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("Azimuth Predicted Cell Subtypes (Level 1)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
plot_scores <- FeaturePlot(Merged_Joined, features = "mapping.score", reduction = "umap") +
  ggtitle("Mapping Confidence (Dark = Healthy, Light = Malignant)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
final_combined_plot <- plot_types + plot_scores
print(final_combined_plot)
ggsave(filename = "/Users/sudharshiniram/Desktop/SRI Data Analysis/Azimuth_Tumor_Map.png", 
       plot = final_combined_plot, 
       width = 16,  
       height = 8,  
       dpi = 300) 

Merged_Joined <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Joined.RDS")
# 1. Recalculate the markers
Merged.markers <- FindAllMarkers(Merged_Joined, only.pos = TRUE, max.cells.per.ident = 500)
# 2. Pull the top 20 genes per cluster
library(dplyr)
top20_genes <- Merged.markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%                  
  slice_max(n = 20, order_by = avg_log2FC)      

# 3. Save the filtered list
write.csv(top20_genes, 
          file = "/Users/sudharshiniram/Desktop/SRI Data Analysis/Top20_Marker_Genes.csv", 
          row.names = FALSE)
'''
raw_counts <- table(Merged_Joined$seurat_clusters, Merged_Joined$predicted.annotation.l1)
percent_table <- prop.table(raw_counts, margin = 1) * 100
percent_table_rounded <- round(percent_table, 2)
print(percent_table_rounded)


write.csv(as.data.frame.matrix(percent_table_rounded), 
          file = "/Users/sudharshiniram/Desktop/SRI Data Analysis/Cluster_Accuracy_Percentages.csv", 
          row.names = TRUE)




# Calculate Marker Genes

gc()

Merged.markers <- FindAllMarkers(Merged_Joined, only.pos = TRUE, max.cells.per.ident = 500)

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

names(my_cell_types) <- levels(Merged_Joined)
Merged_Joined <- RenameIdents(Merged_Joined, my_cell_types)


# Graphing
final_umap <- DimPlot(Merged_Joined, reduction = "umap", label = TRUE, label.size = 3.5, repel = TRUE) + 
  ggtitle("PDAC Tumor Microenvironment by Cell Type") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(filename = "/Users/sudharshiniram/Desktop/SRI Data Analysis/Final_Labeled_UMAP.png", 
       plot = final_umap, 
       width = 14,     
       height = 10,     
       dpi = 300)  

saveRDS(Merged_Joined, file = '/Users/sudharshiniram/Desktop/SRI Data Analysis/Data_Analysis_Final_Labeled.rds')



#Testing Feature Plots
Merged_Labeled <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Data_Analysis_Final_Labeled.rds")

# Neurons
FeaturePlot(Merged_Labeled, features=c("PRPH", "CALCA", "CHAT", "TUBB3", "TH", "NEFL"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged_Labeled, features=c("CHAT", "TRPV1", "SLC18A3", "POU4F1"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged_Labeled, features=c("TUBB3", "UCHL1", "RBFOX3", "SYP"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged_Labeled, features=c("MAP2", "MAPT"), cols = c('lightgrey', 'blue'))

# Macrophages
# Note: Humans do not have Cd209f, the closest is CD209. 
FeaturePlot(Merged_Labeled, features=c("CD209", "CCR2", "MRC1", "CD163"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged_Labeled, features=c("CX3CR1", "CD80", "CD86", "NOS2", "CD38", "ARG1", "MRC1", "CCR2", "MERTK"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged_Labeled, features=c("CD68", "SIGLEC1", "ADGRE1", "CD163"), cols = c('lightgrey', 'blue'))

# Neutrophils
# Note: Humans lack Ly6g/Ly6c1. Replaced with standard human neutrophil markers CEACAM8 (CD66b) and FCGR3B (CD16b).
FeaturePlot(Merged_Labeled, features=c("CEACAM8", "FCGR3B", "CXCR4", "CD33", "FUT4", "CD14"), cols = c('lightgrey', 'blue'))

# Monocytes
# Note: Replaced Ly6g/Ly6c2 with human monocyte markers CD14 and FCGR3A (CD16a).
FeaturePlot(Merged_Labeled, features=c("CD14", "FCGR3A", "CXCR4", "CD33"), cols = c('lightgrey', 'blue'))

# Schwann cells
FeaturePlot(Merged_Labeled, features=c("MBP", "S100B", "MPZ", "NCAM1", "JUN", "NRG1", "SCN7A", "ERBB3", "SOX2"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged_Labeled, features=c("S100B", "GFAP", "NGFR", "MBP", "POU3F1", "EGR2"), cols = c('lightgrey', 'blue'))

# Neuroendocrine cells
FeaturePlot(Merged_Labeled, features=c("CALCB", "SCG5", "SCG2", "CHGB", "PCSK1N", "CHGA", "SCG3"), cols = c('lightgrey', 'blue'))

# Acinar cells
# Note: If this is pancreatic tissue, AMY2A is usually the prevalent amylase over AMY1A.
FeaturePlot(Merged_Labeled, features=c("GP2", "GATA4", "AMY2A", "AMY1A", "PRSS1"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged_Labeled, features=c("LIPA", "LPL", "LIPE", "PNLIP", "PNPLA2"), cols = c('lightgrey', 'blue'))

# Duct cells
FeaturePlot(Merged_Labeled, features=c("KRT19", "SOX9", "PDX1", "P2RY1", "BMPR1A"), cols = c('lightgrey', 'blue'))

# VlnPlot(Merged_Labeled, features = c("NINJ1"))

# apCAFs
# Note: Hladrb1 formatted to the proper human hyphenated symbol HLA-DRB1.
FeaturePlot(Merged_Labeled, features=c("CD74", "ACTA2", "HLA-DRB1"), cols = c('lightgrey', 'blue'))

# Mesothelial cells
FeaturePlot(Merged_Labeled, features=c("CALB2", "KRT5", "KRT6A", "CD44", "THBD"), cols = c('lightgrey', 'blue'))

# B cells
FeaturePlot(Merged_Labeled, features=c("CD19", "PTPRC", "CD27", "IGHD", "TNFRSF17", "CXCR4", "SDC1", "SLAMF7"), cols = c('lightgrey', 'blue'))

# T cells
FeaturePlot(Merged_Labeled, features=c("CD27", "CD8A", "SELL", "CD44", "CD3E", "CD4"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged_Labeled, features=c("CD8A", "CD3E", "CD4"), cols = c('lightgrey', 'blue'))

# Endocrine cells
# Note: Humans only have one insulin gene (INS), not Ins1/Ins2.
FeaturePlot(Merged_Labeled, features=c("TFRC", "GCG", "SST", "INS"), cols = c('lightgrey', 'blue'))

# Eosinophils
# Note: The functional human equivalent of mouse Siglecf is SIGLEC8.
FeaturePlot(Merged_Labeled, features=c("CCR3", "ITGA2", "SIGLEC8", "CEACAM8"), cols = c('lightgrey', 'blue'))

# Basophils
# Note: Crth2 is an alias; the official HGNC symbol is PTGDR2. I included both just in case your reference genome uses the alias.
FeaturePlot(Merged_Labeled, features=c("IL3RA", "ENPP3", "PTGDR2", "CRTH2", "CD63"), cols = c('lightgrey', 'blue'))