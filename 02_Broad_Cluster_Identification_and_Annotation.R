library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

Sys.setenv("R_MAX_VSIZE" = 32000000000)

# ============================================================
# 0. Check that Merged is already loaded
# ============================================================

if (!exists("Merged")) {
  stop("Merged object is not loaded in the R environment. Open your saved RDS file first using readRDS().")
}

# Optional: save your current identities before reclustering
Merged$old_cluster_or_label <- Idents(Merged)

# ============================================================
# 1. Re-run PCA with more PCs
# ============================================================

set.seed(1234)

Merged <- RunPCA(
  object = Merged,
  features = VariableFeatures(object = Merged),
  npcs = 50
)

# Save elbow plot
elbow_plot <- ElbowPlot(Merged, ndims = 50) +
  ggtitle("Elbow Plot: First 50 PCs")

ggsave(
  filename = "/Users/sudharshiniram/Desktop/SRI Data Analysis/ElbowPlot_50PCs.png",
  plot = elbow_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# 2. Use 30 PCs for broad clustering
# ============================================================

pc_use <- 1:30

Merged <- FindNeighbors(
  object = Merged,
  dims = pc_use
)

Merged <- FindClusters(
  object = Merged,
  resolution = 0.2
)

# Save broad cluster IDs before labeling
Merged$broad_res_0.2_cluster <- Idents(Merged)

# ============================================================
# 3. Re-run UMAP with broad settings
# ============================================================

Merged <- RunUMAP(
  object = Merged,
  dims = pc_use,
  n.neighbors = 30,
  min.dist = 0.3
)

# ============================================================
# 4. Plot broad computational clusters before cell-type labels
# ============================================================

broad_cluster_umap <- DimPlot(
  object = Merged,
  reduction = "umap",
  group.by = "broad_res_0.2_cluster",
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 0.08
) +
  NoLegend() +
  ggtitle("PDAC Broad Clusters, Resolution 0.2") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

print(broad_cluster_umap)

ggsave(
  filename = "/Users/sudharshiniram/Desktop/SRI Data Analysis/Broad_Clusters_Resolution_0.2_UMAP.png",
  plot = broad_cluster_umap,
  width = 12,
  height = 9,
  dpi = 300
)

# ============================================================
# 5. Get top marker genes for the new broad clusters
# ============================================================

Merged <- JoinLayers(Merged)
gc()

'''
broad_markers <- FindAllMarkers(
  object = Merged,
  only.pos = TRUE,
  max.cells.per.ident = 500
)

write.csv(
  broad_markers,
  file = "/Users/sudharshiniram/Desktop/SRI Data Analysis/Broad_Resolution_0.2_All_Markers.csv",
  row.names = FALSE
)

top30_broad_markers <- broad_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(
  top30_broad_markers,
  file = "/Users/sudharshiniram/Desktop/SRI Data Analysis/Broad_Resolution_0.2_Top30_Markers.csv",
  row.names = FALSE
)

# ============================================================
# 6. Broad marker dot plot to help label clusters
# ============================================================

broad_marker_panel <- c(
  # epithelial / ductal / PDAC
  "EPCAM", "KRT19", "KRT8", "KRT18", "SOX9", "MUC1", "MMP7", "TSPAN8", "LCN2",
  
  # acinar
  "PRSS1", "PRSS2", "CPA1", "CTRB1", "CTRB2", "CELA3A", "PNLIP",
  
  # endocrine
  "CHGA", "CHGB", "INS", "GCG", "SST", "PPY",
  
  # CAF / fibroblast
  "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "ACTA2", "TAGLN", "IL6", "CXCL12",
  
  # endothelial / lymphatic / pericyte
  "PECAM1", "VWF", "CDH5", "PLVAP", "PROX1", "PDPN", "LYVE1", "RGS5", "PDGFRB",
  
  # immune
  "PTPRC", "LYZ", "C1QA", "C1QB", "CD68", "MRC1", "CD163",
  "CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "MS4A1", "CD79A",
  
  # plasma / mast / Schwann
  "MZB1", "JCHAIN", "TPSAB1", "CPA3", "SOX10", "S100B", "PLP1", "MPZ", "GFAP"
)

# Only keep genes that actually exist in your dataset
broad_marker_panel <- broad_marker_panel[broad_marker_panel %in% rownames(Merged)]

broad_dotplot <- DotPlot(
  object = Merged,
  features = broad_marker_panel,
  group.by = "broad_res_0.2_cluster"
) +
  RotatedAxis() +
  ggtitle("Broad Marker DotPlot, Resolution 0.2") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 10)
  )

print(broad_dotplot)

ggsave(
  filename = "/Users/sudharshiniram/Desktop/SRI Data Analysis/Broad_Resolution_0.2_DotPlot.png",
  plot = broad_dotplot,
  width = 18,
  height = 8,
  dpi = 300
)

# ============================================================
# 7. Save updated object
# ============================================================

saveRDS(
  object = Merged,
  file = "/Users/sudharshiniram/Desktop/SRI Data Analysis/Data_Analysis_Broad_Resolution_0.2.rds",
  compress = FALSE
)

'''


DefaultAssay(Merged) <- "RNA"

# ============================================================
# Labeled FeaturePlot gene sets for human PDAC snRNA-seq
# ============================================================

feature_plot_genes <- list(
  
  # 1. Broad overview
  "01_Broad_Lineage_Overview" = c(
    "EPCAM", "KRT19", "KRT8", "KRT18",
    "PRSS1", "CPA1", "CTRB1",
    "CHGA", "INS", "GCG",
    "COL1A1", "COL1A2", "DCN", "LUM",
    "PTPRC", "LYZ", "CD3E", "MS4A1",
    "PECAM1", "VWF", "CDH5",
    "SOX10", "S100B"
  ),
  
  # 2. Epithelial / ductal / PDAC
  "02_Epithelial_Ductal_PDAC" = c(
    "EPCAM", "KRT19", "KRT8", "KRT18", "KRT7",
    "SOX9", "MUC1", "MSLN", "TSPAN8", "MMP7", "LCN2"
  ),
  
  # 3. Classical / gastric-like epithelial
  "03_Classical_Gastric_Epithelial" = c(
    "GATA6", "TFF1", "TFF2", "TFF3",
    "AGR2", "CLDN18", "CTSE", "MUC5AC", "MUC6"
  ),
  
  # 4. Basal-like / squamous-like epithelial
  "04_Basal_Squamous_Epithelial" = c(
    "KRT5", "KRT6A", "KRT14", "KRT17",
    "TP63", "S100A2", "LAMC2", "ITGA6", "ITGB4"
  ),
  
  # 5. Acinar cells
  "05_Acinar_Cells" = c(
    "PRSS1", "PRSS2", "CPA1", "CPA2",
    "CTRB1", "CTRB2", "CELA3A", "CELA3B",
    "AMY2A", "AMY2B", "PNLIP", "CLPS",
    "REG1A", "REG1B", "GP2"
  ),
  
  # 6. Endocrine / islet cells
  "06_Endocrine_Islet_Cells" = c(
    "CHGA", "CHGB", "SCG2", "SCG3", "SCG5",
    "PCSK1", "PCSK2", "SLC30A8", "RFX6", "ISL1",
    "INS", "IAPP",
    "GCG", "ARX",
    "SST",
    "PPY", "PYY"
  ),
  
  # 7. CAF / fibroblast
  "07_CAF_Fibroblast" = c(
    "COL1A1", "COL1A2", "COL3A1",
    "DCN", "LUM", "COL6A1", "COL6A2", "COL6A3",
    "FAP", "PDGFRA", "PDGFRB", "THY1", "FN1"
  ),
  
  # 8. myCAF / activated CAF
  "08_myCAF_Activated_CAF" = c(
    "ACTA2", "TAGLN", "MYL9", "TPM2",
    "POSTN", "MMP11", "LRRC15",
    "COL11A1", "COL12A1", "THBS2", "INHBA"
  ),
  
  # 9. iCAF / inflammatory CAF
  "09_iCAF_Inflammatory_CAF" = c(
    "IL6", "CXCL12", "CXCL14", "CXCL1", "CXCL2",
    "CCL2", "CFD", "HAS1", "LIF", "PLA2G2A", "PDGFRA"
  ),
  
  # 10. apCAF / antigen-presenting CAF
  "10_apCAF_Markers" = c(
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
    "CD74", "CIITA", "SLPI",
    "COL1A1", "DCN", "LUM"
  ),
  
  # 11. Blood vascular endothelial
  "11_Endothelial_Vascular" = c(
    "PECAM1", "VWF", "CDH5", "KDR", "FLT1",
    "ESAM", "EMCN", "RAMP2", "RAMP3",
    "PLVAP", "ENG", "TEK", "TIE1", "EGFL7", "NOS3"
  ),
  
  # 12. Lymphatic endothelial
  "12_Lymphatic_Endothelial" = c(
    "PROX1", "PDPN", "LYVE1", "FLT4",
    "CCL21", "MMRN1", "TFF3"
  ),
  
  # 13. Pericyte / vascular smooth muscle
  "13_Pericyte_VSMC" = c(
    "RGS5", "MCAM", "CSPG4", "PDGFRB",
    "NOTCH3", "ACTA2", "MYH11", "DES"
  ),
  
  # 14. Macrophage / TAM / myeloid
  "14_Macrophage_TAM_Myeloid" = c(
    "PTPRC", "LYZ", "LST1", "TYROBP",
    "C1QA", "C1QB", "C1QC",
    "CD68", "MRC1", "CD163", "MSR1",
    "APOE", "MARCO", "SPP1", "GPNMB", "FOLR2", "CCL18"
  ),
  
  # 15. Monocytes
  "15_Monocytes" = c(
    "LYZ", "FCN1", "VCAN", "S100A8", "S100A9",
    "LST1", "CD14", "CCR2", "FCGR3A"
  ),
  
  # 16. Neutrophil / granulocyte-like
  "16_Neutrophil_Granulocyte" = c(
    "S100A8", "S100A9", "FCGR3B", "CSF3R",
    "CXCR2", "MMP9", "MNDA", "LCN2", "CEACAM8"
  ),
  
  # 17. Dendritic cells
  "17_Dendritic_Cells" = c(
    "FCER1A", "CLEC10A", "CD1C", "ITGAX",
    "HLA-DRA", "HLA-DPB1", "CD74",
    "CLEC9A", "XCR1", "LILRA4", "CLEC4C"
  ),
  
  # 18. General T cells
  "18_T_Cells_General" = c(
    "PTPRC", "CD3D", "CD3E", "CD3G", "TRAC",
    "CD2", "CD247", "IL7R", "LTB"
  ),
  
  # 19. CD4 / CD8 / cytotoxic T cells
  "19_CD4_CD8_Cytotoxic_T_Cells" = c(
    "CD4", "IL7R", "CCR7", "TCF7",
    "CD8A", "CD8B",
    "GZMA", "GZMB", "PRF1", "CCL5", "CST7", "NKG7"
  ),
  
  # 20. NK cells
  "20_NK_Cells" = c(
    "NKG7", "GNLY", "KLRD1", "KLRF1",
    "PRF1", "GZMB", "GZMH", "FCGR3A", "TYROBP"
  ),
  
  # 21. Treg / exhausted T cell states
  "21_Treg_Exhausted_T_Cells" = c(
    "FOXP3", "IL2RA", "CTLA4", "IKZF2",
    "PDCD1", "LAG3", "TIGIT", "HAVCR2", "TOX"
  ),
  
  # 22. B cells
  "22_B_Cells" = c(
    "MS4A1", "CD79A", "CD79B", "BANK1",
    "CD19", "CD37", "CD74", "HLA-DRA"
  ),
  
  # 23. Plasma cells
  "23_Plasma_Cells" = c(
    "MZB1", "JCHAIN", "XBP1",
    "IGKC", "IGHG1", "IGHG3", "SDC1", "DERL3"
  ),
  
  # 24. Mast cells
  "24_Mast_Cells" = c(
    "TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HDC"
  ),
  
  # 25. Eosinophils
  "25_Eosinophils" = c(
    "CCR3", "IL5RA", "CLC", "PRG2", "RNASE2", "SIGLEC8"
  ),
  
  # 26. Basophils
  "26_Basophils" = c(
    "IL3RA", "ENPP3", "FCER1A", "FCER1G", "MS4A2", "HDC"
  ),
  
  # 27. Schwann / neural glial cells
  "27_Schwann_Neural_Glial" = c(
    "SOX10", "S100B", "PLP1", "MPZ", "MBP",
    "GFAP", "NGFR", "ERBB3", "POU3F1",
    "NRG1", "SCN7A", "NCAM1", "CRYAB", "GPM6B", "CDH19"
  ),
  
  # 28. Neuron-like markers
  "28_Neuron_like_Markers" = c(
    "RBFOX3", "SNAP25", "TUBB3", "UCHL1",
    "SYT1", "SYP", "NEFL", "NEFM",
    "MAP2", "MAPT", "TH", "CHAT", "PRPH"
  ),
  
  # 29. Cycling / proliferating cells
  "29_Cycling_Proliferating_Cells" = c(
    "MKI67", "TOP2A", "UBE2C", "CENPF",
    "STMN1", "HMGB2", "PCNA", "TYMS"
  ),
  
  # 30. Mesothelial-like cells
  "30_Mesothelial_like_Cells" = c(
    "CALB2", "KRT5", "KRT6A", "CD44",
    "THBD", "MSLN", "WT1", "UPK3B"
  ),
  
  # 31. Mixed / doublet signal check
  "31_Mixed_Doublet_Check" = c(
    "EPCAM", "KRT19",
    "COL1A1", "COL3A1",
    "PTPRC", "LYZ",
    "PECAM1", "VWF",
    "PRSS1", "CPA1",
    "CHGA", "INS"
  )
)


featureplot_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis/Broad_Celltype_FeaturePlots"

if (!dir.exists(featureplot_dir)) {
  dir.create(featureplot_dir, recursive = TRUE)
}


make_feature_plot <- function(plot_name, genes, object = Merged, ncol = 4) {
  
  genes_present <- genes[genes %in% rownames(object)]
  genes_missing <- genes[!genes %in% rownames(object)]
  
  cat("\n==============================\n")
  cat("Plot:", plot_name, "\n")
  cat("Genes present:", paste(genes_present, collapse = ", "), "\n")
  cat("Genes missing:", paste(genes_missing, collapse = ", "), "\n")
  
  if (length(genes_present) == 0) {
    warning(paste("No genes found for:", plot_name))
    return(NULL)
  }
  
  p <- FeaturePlot(
    object = object,
    features = genes_present,
    reduction = "umap",
    cols = c("lightgrey", "blue"),
    order = TRUE,
    raster = TRUE,
    ncol = ncol
  ) +
    plot_annotation(title = gsub("_", " ", plot_name))
  
  return(p)
}


feature_plots <- list()

for (plot_name in names(feature_plot_genes)) {
  feature_plots[[plot_name]] <- make_feature_plot(
    plot_name = plot_name,
    genes = feature_plot_genes[[plot_name]],
    object = Merged,
    ncol = 4
  )
}


for (plot_name in names(feature_plots)) {
  
  p <- feature_plots[[plot_name]]
  
  if (!is.null(p)) {
    ggsave(
      filename = file.path(featureplot_dir, paste0(plot_name, ".png")),
      plot = p,
      width = 16,
      height = 10,
      dpi = 300
    )
  }
}

feature_plots[["01_Broad_Lineage_Overview"]]
feature_plots[["27_Schwann_Neural_Glial"]]
feature_plots[["31_Mixed_Doublet_Check"]]
