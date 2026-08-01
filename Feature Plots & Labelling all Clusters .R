library(Seurat)
library(ggplot2)

rm(fibroblast_subset, Merged, Merged_Integrated, Merged_Joined, Merged_small, catch_obj, expr_matrix, catch_cluster_ids, Merged_sccatch)
gc()

#Testing Feature Plots
Merged <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Joined.RDS")
DimPlot(Merged, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("UMAP of Merged Dataset")

# ============================================================
# BROAD CELL TYPE FEATURE PLOTS - HUMAN PDAC snRNA-seq
# Each FeaturePlot has max 6 genes.
# ============================================================

FeaturePlot(Merged, features=c("PMP22"), cols = c('lightgrey', 'blue'))

# ------------------------------------------------------------
# Epithelial / Ductal / Malignant Tumor Cells
# Broad epithelial markers: use these first to find epithelial clusters.
# Malignant vs normal ductal usually requires CNV/inferCNV or tumor-enriched markers.
# ------------------------------------------------------------
FeaturePlot(Merged, features=c("EPCAM", "KRT19", "KRT8", "KRT7", "MUC1", "SOX9"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("CEACAM6", "TSPAN8", "FXYD3", "MSLN", "SPP1", "CFTR"), cols = c('lightgrey', 'blue'))


# ------------------------------------------------------------
# Acinar Cells
# Pancreatic enzyme-producing cells.
# ------------------------------------------------------------
FeaturePlot(Merged, features=c("PRSS1", "PRSS2", "CTRB2", "PNLIP", "CEL", "GP2"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("AMY2A", "REG1A", "PTF1A", "GATA4"), cols = c('lightgrey', 'blue'))


# ------------------------------------------------------------
# Endocrine / Islet Cells
# Broad endocrine/islet markers plus alpha, beta, delta cell hormones.
# ------------------------------------------------------------
FeaturePlot(Merged, features=c("CHGA", "CHGB", "INSM1", "ISL1", "SCG2", "PCSK1"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("INS", "GCG", "SST", "PPY"), cols = c('lightgrey', 'blue'))


# ------------------------------------------------------------
# Fibroblasts / Cancer-Associated Fibroblasts - CAFs
# Broad stromal / fibroblast markers.
# ------------------------------------------------------------
FeaturePlot(Merged, features=c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "SPARC"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("FN1", "ACTA2", "PDGFRA", "FAP", "THY1", "MMP2"), cols = c('lightgrey', 'blue'))

#Subsetting Fibroblasts
fibroblast_subset <- subset(Merged, idents = "1") 
fibroblast_subset@reductions$pca <- NULL
fibroblast_subset@reductions$harmony <- NULL
fibroblast_subset@reductions$umap <- NULL
fibroblast_subset <- FindVariableFeatures(fibroblast_subset, selection.method = "vst", nfeatures = 2000)
fibroblast_subset <- ScaleData(fibroblast_subset)
fibroblast_subset <- RunPCA(fibroblast_subset, features = VariableFeatures(object = fibroblast_subset))
fibroblast_subset <- FindNeighbors(fibroblast_subset, reduction = "pca", dims = 1:15)
fibroblast_subset <- FindClusters(fibroblast_subset, resolution = 0.2)
fibroblast_subset <- RunUMAP(fibroblast_subset, reduction = "pca", dims = 1:15, reduction.name = "sub_umap")
DimPlot(fibroblast_subset, reduction = "sub_umap", label = TRUE, repel = TRUE) + 
  ggtitle("Fibroblast Sub-Clusters (Cluster 1)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#myCAFs
FeaturePlot(fibroblast_subset, features=c("ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2", "POSTN"), cols = c('lightgrey', 'blue'))
# iCAFs - inflammatory CAFs
FeaturePlot(fibroblast_subset, features=c("IL6", "CXCL12", "CXCL14", "LIF", "CFD", "PDGFRA"), cols = c('lightgrey', 'blue'))
# apCAFs - antigen-presenting CAFs
FeaturePlot(fibroblast_subset, features=c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA"), cols = c('lightgrey', 'blue'))
# # ECM-remodeling / matrix CAFs
FeaturePlot(Merged, features=c("COL11A1", "COL10A1", "POSTN", "MMP11", "THBS2", "INHBA"), cols = c('lightgrey', 'blue'))

# Group by your new sub-clusters to see exactly which cluster number lights up
DotPlot(fibroblast_subset, 
        features = c("IL6", "CXCL12", "CXCL14", "LIF", "CFD", "PDGFRA"), 
        group.by = "seurat_clusters") + 
  coord_flip() + # Flips the axes so it's easier to read
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("iCAF Marker Expression by Sub-cluster")

DotPlot(fibroblast_subset, 
        features = c("ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2", "POSTN"), 
        group.by = "seurat_clusters") + 
  coord_flip() + # Flips the axes so it's easier to read
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("myCAF Marker Expression by Sub-cluster")

DotPlot(fibroblast_subset, 
        features = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA"), 
        group.by = "seurat_clusters") + 
  coord_flip() + # Flips the axes so it's easier to read
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("apCAF Marker Expression by Sub-cluster")

# ------------------------------------------------------------
# Endothelial Cells
# Vascular endothelial markers.
# ------------------------------------------------------------
FeaturePlot(Merged, features=c("PECAM1", "VWF", "PLVAP", "KDR", "FLT1", "ESAM"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("SPARCL1", "EMCN", "CD34"), cols = c('lightgrey', 'blue'))


# ------------------------------------------------------------
# Pericytes / Smooth Muscle / Mural Cells
# Often near endothelial cells but separate from true fibroblasts.
# ------------------------------------------------------------
FeaturePlot(Merged, features=c("RGS5", "PDGFRB", "CSPG4", "MCAM", "NOTCH3", "ACTA2"), cols = c('lightgrey', 'blue'))


# ------------------------------------------------------------
# Myeloid Lineage / Macrophages / Monocytes
# Broad myeloid markers, then TAM/monocyte-like markers.
# ------------------------------------------------------------
FeaturePlot(Merged, features=c("LST1", "TYROBP", "FCER1G", "CD68", "C1QA", "C1QB"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("C1QC", "FCN1", "CD14", "FCGR3A", "MRC1", "CD163"), cols = c('lightgrey', 'blue'))


# ------------------------------------------------------------
# Granulocytes
# Neutrophils, eosinophils, and basophils may be rare in snRNA-seq.
# ------------------------------------------------------------
FeaturePlot(Merged, features=c("CEACAM8", "FCGR3B", "CXCR2", "S100A8", "S100A9", "MPO"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("CCR3", "SIGLEC8", "CLC", "IL3RA", "ENPP3", "PTGDR2"), cols = c('lightgrey', 'blue'))


# ------------------------------------------------------------
# T Cells
# Broad T cell markers plus CD4/CD8/cytotoxic markers.
# ------------------------------------------------------------
FeaturePlot(Merged, features=c("CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "IL7R"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("CD4", "CD8A", "CD8B", "NKG7", "PRF1", "GZMB"), cols = c('lightgrey', 'blue'))


# ------------------------------------------------------------
# NK Cells
# NK cells can overlap with cytotoxic CD8 T cells, so compare with CD3D/CD3E.
# ------------------------------------------------------------
FeaturePlot(Merged, features=c("NKG7", "GNLY", "KLRD1", "KLRF1", "FCGR3A", "NCAM1"), cols = c('lightgrey', 'blue'))


# ------------------------------------------------------------
# B Cells / Plasma Cells
# B cell markers followed by plasma cell markers.
# ------------------------------------------------------------
FeaturePlot(Merged, features=c("MS4A1", "CD79A", "CD79B", "CD19", "CD22", "BLK"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("JCHAIN", "MZB1", "SDC1", "CD38", "SLAMF7", "TNFRSF17"), cols = c('lightgrey', 'blue'))


# ------------------------------------------------------------
# Mast Cells
# Rare but important immune population in PDAC.
# ------------------------------------------------------------
FeaturePlot(Merged, features=c("TPSAB1", "TPSB2", "CPA3", "MS4A2", "KIT", "HDC"), cols = c('lightgrey', 'blue'))


# ------------------------------------------------------------
# Schwann Cells / Neural Lineage
# Important for PDAC neural invasion and nerve-associated tumor biology.
# ------------------------------------------------------------
FeaturePlot(Merged, features=c("SOX10", "S100B", "PLP1", "CDH19", "GFAP", "NGFR"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("PMP22", "MPZ", "MBP", "ERBB3", "POU3F1", "EGR2"), cols = c('lightgrey', 'blue'))

FeaturePlot(Merged, features=c("JUN", "NCAM1", "SCN7A", "SOX2", "MBP"), cols = c('lightgrey', 'blue'))


# ------------------------------------------------------------
# Neurons / Neural Cells
# True neuronal nuclei may be rare in PDAC snRNA-seq, but these help detect neural contamination or nerve-associated clusters.
# ------------------------------------------------------------
FeaturePlot(Merged, features=c("PRPH", "CALCA", "CHAT", "TUBB3", "TH", "NEFL"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("TRPV1", "SLC18A3", "POU4F1", "UCHL1", "RBFOX3", "SYP"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("MAP2", "MAPT"), cols = c('lightgrey', 'blue'))


# ------------------------------------------------------------
# Mesothelial Cells
# Can appear in pancreatic/peritoneal-associated samples.
# ------------------------------------------------------------
FeaturePlot(Merged, features=c("CALB2", "KRT5", "KRT6A", "MSLN", "WT1", "THBD"), cols = c('lightgrey', 'blue'))

#Healthy Ductal
FeaturePlot(Merged, features=c("CFTR", "AMBP"), cols = c('lightgrey', 'blue'))
# Healthy / normal ductal-enriched cells
# Best markers to separate normal ductal cells from malignant ductal-like tumor cells
FeaturePlot(Merged, features=c("CFTR", "SLC4A4", "SCTR", "BICC1", "PKHD1", "CA2"), cols = c('lightgrey', 'blue'))
FeaturePlot(Merged, features=c("ANXA4", "SLC26A6"), cols = c('lightgrey', 'blue'))


# ============================================================
# OPTIONAL SUBTYPE FEATURE PLOTS
# Commented out so they will NOT run automatically.
# Uncomment after broad cell types are assigned.
# ============================================================


# ------------------------------------------------------------
# Epithelial Subtypes
# ------------------------------------------------------------

# # Normal ductal-like epithelial cells
# FeaturePlot(Merged, features=c("CFTR", "KRT19", "KRT7", "KRT8", "SOX9", "SPP1"), cols = c('lightgrey', 'blue'))

# # Malignant / PDAC-enriched ductal epithelial cells
# FeaturePlot(Merged, features=c("EPCAM", "KRT19", "MUC1", "CEACAM6", "TSPAN8", "FXYD3"), cols = c('lightgrey', 'blue'))
# FeaturePlot(Merged, features=c("MSLN", "KRT17", "S100A6", "S100A10"), cols = c('lightgrey', 'blue'))


# ------------------------------------------------------------
# CAF Subtypes
# ------------------------------------------------------------

# # myCAFs - myofibroblastic CAFs
# FeaturePlot(Merged, features=c("ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2", "POSTN"), cols = c('lightgrey', 'blue'))

# # iCAFs - inflammatory CAFs
# FeaturePlot(Merged, features=c("IL6", "CXCL12", "CXCL14", "LIF", "CFD", "PDGFRA"), cols = c('lightgrey', 'blue'))

# # apCAFs - antigen-presenting CAFs
# FeaturePlot(Merged, features=c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA"), cols = c('lightgrey', 'blue'))

# # ECM-remodeling / matrix CAFs
# FeaturePlot(Merged, features=c("COL11A1", "COL10A1", "POSTN", "MMP11", "THBS2", "INHBA"), cols = c('lightgrey', 'blue'))


# ------------------------------------------------------------
# Myeloid Subtypes
# ------------------------------------------------------------

# # Monocytes
#FeaturePlot(Merged, features=c("FCN1", "VCAN", "S100A8", "S100A9", "CD14", "CCR2"), cols = c('lightgrey', 'blue'))

# # Tissue-resident / complement-high TAMs
# FeaturePlot(Merged, features=c("C1QA", "C1QB", "C1QC", "APOE", "MRC1", "CD163"), cols = c('lightgrey', 'blue'))

# # Inflammatory macrophages
#FeaturePlot(Merged, features=c("IL1B", "CXCL8", "NFKBIA", "TNF", "S100A8", "S100A9"), cols = c('lightgrey', 'blue'))

# # Dendritic cells
# FeaturePlot(Merged, features=c("FCER1A", "CLEC10A", "CD1C", "LILRA4", "CLEC9A", "XCR1"), cols = c('lightgrey', 'blue'))


# ------------------------------------------------------------
# T Cell Subtypes
# ------------------------------------------------------------

# # CD4 memory/helper T cells
# FeaturePlot(Merged, features=c("CD3D", "CD3E", "CD4", "IL7R", "SELL", "CCR7"), cols = c('lightgrey', 'blue'))

# # CD8 cytotoxic T cells
# FeaturePlot(Merged, features=c("CD3D", "CD3E", "CD8A", "CD8B", "NKG7", "GZMB"), cols = c('lightgrey', 'blue'))

# # Exhausted / dysfunctional T cells
# FeaturePlot(Merged, features=c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CTLA4"), cols = c('lightgrey', 'blue'))

# # Regulatory T cells
# FeaturePlot(Merged, features=c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT", "TNFRSF18"), cols = c('lightgrey', 'blue'))


# ------------------------------------------------------------
# B Cell / Plasma Subtypes
# ------------------------------------------------------------

# # Naive / mature B cells
# FeaturePlot(Merged, features=c("MS4A1", "CD79A", "CD79B", "CD19", "IGHD", "TCL1A"), cols = c('lightgrey', 'blue'))

# # Plasma cells
# FeaturePlot(Merged, features=c("JCHAIN", "MZB1", "SDC1", "CD38", "XBP1", "PRDM1"), cols = c('lightgrey', 'blue'))


# Check DCs
#FeaturePlot(Merged, features = c("CD1C", "CLEC9A", "FLT3", "LAMP3"), cols = c('lightgrey', 'blue'))

# Check for Doublets (Compare these side-by-side)
#FeaturePlot(Merged, features = c("CD3E", "LYZ"), cols = c('lightgrey', 'blue'))

# ------------------------------------------------------------
# Endothelial Subtypes
# ------------------------------------------------------------

# # General vascular endothelial cells
# FeaturePlot(Merged, features=c("PECAM1", "VWF", "KDR", "FLT1", "ESAM", "CD34"), cols = c('lightgrey', 'blue'))

# # Lymphatic endothelial cells
# FeaturePlot(Merged, features=c("PROX1", "PDPN", "LYVE1", "FLT4", "CCL21", "MMRN1"), cols = c('lightgrey', 'blue'))


# ------------------------------------------------------------
# Schwann / Neural Subtypes
# ------------------------------------------------------------

# # Myelinating Schwann cells
# FeaturePlot(Merged, features=c("MPZ", "MBP", "PMP22", "PLP1", "EGR2", "POU3F1"), cols = c('lightgrey', 'blue'))

# # Repair-like / activated Schwann cells
# FeaturePlot(Merged, features=c("NGFR", "GFAP", "JUN", "SOX2", "ERBB3", "S100B"), cols = c('lightgrey', 'blue'))