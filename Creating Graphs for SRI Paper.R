# Load required libraries
library(Seurat)
library(ggplot2)

# ============================================================
# INITIALIZATION & DATA LOADING
# ============================================================
# Ensure we are using a consistent object name throughout the script
Merged_Joined <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Joined.RDS")


# ============================================================
# 1. PUBLICATION-QUALITY UMAP
# ============================================================
# Generate a clean, grid-free UMAP with professional typography
umap_plot <- DimPlot(Merged_Joined, reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.5) + 
  ggtitle("UMAP Projection of PDAC snRNA-seq Samples") +
  theme_classic(base_family = "sans") + # Standardizes font to Arial/Helvetica, removes grid lines
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10)
  )

# Print the UMAP
print(umap_plot)


# ============================================================
# 2. RENAME CLUSTERS TO BIOLOGICAL NAMES
# ============================================================
# Replace numerical clusters with formal biological classifications
new_cluster_ids <- c(
  "0" = "Cancer 1",
  "1" = "Fibroblast",
  "2" = "Healthy Ductal",
  "3" = "Lymphoid B Cells",
  "4" = "Myeloid",
  "5" = "Acinar",
  "6" = "Endocrine / Islet",
  "7" = "Endothelial",
  "8" = "Lymphoid T/NK Cells",
  "9" = "Cancer 2",
  "10" = "Pericytes / Mural",
  "11" = "Schwann Cells"
)

# Apply the formal names to the object
Idents(Merged_Joined) <- "seurat_clusters" 
Merged_Joined <- RenameIdents(Merged_Joined, new_cluster_ids)


# ============================================================
# 3. DEFINE THE MASTER GENE LIST
# ============================================================
# Grouped logically by lineage to create a striking diagonal pattern
master_genes <- c(
  # Epithelial & Cancer (0, 2, 9)
  "EPCAM", "KRT19", "CEACAM6", "MUC1", "CFTR", "SLC4A4", "SCTR",
  # Stromal & Fibroblast (1, 10)
  "COL1A1", "DCN", "PDGFRA", "FAP", "RGS5", "PDGFRB", "NOTCH3",
  # Immune (3, 4, 8)
  "MS4A1", "CD79A", "CD19", "CD68", "CD14", "C1QA", "CD3D", "CD8A", "NKG7",
  # Endothelial (7)
  "PECAM1", "VWF", "PLVAP",
  # Pancreatic Functional (5, 6)
  "PRSS1", "PNLIP", "AMY2A", "CHGA", "INS", "GCG",
  # Neural / Schwann (11)
  "SOX10", "S100B", "PLP1"
)


# ============================================================
# 4. PUBLICATION-QUALITY DOTPLOT
# ============================================================
master_dotplot <- DotPlot(Merged_Joined, features = master_genes) + 
  coord_flip() + # Places genes on the Y-axis and cell types on the X-axis
  theme_minimal(base_family = "sans") +
  theme(
    # X-Axis: Cell Types
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 10)),
    
    # Y-Axis: Genes (Italicized by scientific convention)
    axis.text.y = element_text(face = "italic", size = 10, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 10)),
    
    # Title & Legend
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    
    # Clean up panel borders
    panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.minor = element_blank()
  ) +
  # Use a striking, high-contrast color palette 
  scale_color_gradient(low = "#E6E6E6", high = "navy") + 
  ggtitle("Master Cell Type Marker Expression") +
  labs(x = "Identified Cell Type", y = "Marker Gene")

# Print the final Master DotPlot
print(master_dotplot)

# ============================================================
# ALTERNATE UMAP: NAMES IN LEGEND ONLY
# ============================================================
umap_named_legend <- DimPlot(Merged_Joined, reduction = "umap", label = FALSE, pt.size = 0.5) + 
  ggtitle("UMAP Projection of PDAC snRNA-seq Samples") +
  theme_classic(base_family = "sans") + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_blank() # Keeps the legend clean by removing the default title
  )

# Print the UMAP
print(umap_named_legend)

# ============================================================
# ALTERNATE DOTPLOT: CLUSTER NUMBERS ON X-AXIS
# ============================================================
# Note the addition of: group.by = "seurat_clusters"
numbered_dotplot <- DotPlot(Merged_Joined, features = master_genes, group.by = "seurat_clusters") + 
  coord_flip() + 
  theme_minimal(base_family = "sans") +
  theme(
    # X-Axis: Cluster Numbers
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 11, color = "black"), # Angle set to 0 since numbers are short
    axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 10)),
    
    # Y-Axis: Genes
    axis.text.y = element_text(face = "italic", size = 10, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 10)),
    
    # Title & Legend
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    
    # Clean up panel borders
    panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.minor = element_blank()
  ) +
  scale_color_gradient(low = "#E6E6E6", high = "navy") + 
  ggtitle("Cluster Marker Expression") +
  labs(x = "Marker Gene", y = "Cluster Number")

# Print the numbered DotPlot
print(numbered_dotplot)

# Load required libraries
library(Seurat)
library(ggplot2)
library(scales)    # Required to extract default Seurat colors
library(patchwork) # Allows us to stitch the plots together nicely

# ============================================================
# 1. PREPARE THE LABELS AND COLORS
# ============================================================
# Set active identities back to pure numbers for the on-graph labels
Idents(Merged_Joined) <- "seurat_clusters" 

# Generate a fixed 12-color palette to guarantee identical colors across both plots
cluster_colors <- hue_pal()(12)
names(cluster_colors) <- 0:11

# Define the highly descriptive legend labels
legend_labels <- c(
  "0" = "0 - Cancer 1",
  "1" = "1 - Fibroblast",
  "2" = "2 - Healthy Ductal",
  "3" = "3 - Lymphoid B Cells",
  "4" = "4 - Myeloid",
  "5" = "5 - Acinar",
  "6" = "6 - Endocrine / Islet",
  "7" = "7 - Endothelial",
  "8" = "8 - Lymphoid T/NK Cells",
  "9" = "9 - Cancer 2",
  "10" = "10 - Pericytes / Mural",
  "11" = "11 - Schwann Cells"
)


# Load required libraries
library(Seurat)
library(ggplot2)
library(scales)    # Required to extract default ggplot colors
library(patchwork) # Allows us to stitch the plots together seamlessly



# ============================================================
# 1. INJECT CLINICAL METADATA (SEX)
# ============================================================
# Mapping each project ID to its clinical sex
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

# THE FIX: Use unname() to strip the names off the vector so Seurat accepts it
Merged_Joined$Sex <- unname(sex_mapping[as.character(Merged_Joined$orig.ident)])


# ============================================================
# 2. PREPARE IDENTITIES, LABELS, AND COLORS
# ============================================================
# Set active identities to pure numbers so only numbers appear on the graph
Idents(Merged_Joined) <- "seurat_clusters" 

# Generate a fixed 12-color palette to guarantee identical colors across both plots
cluster_colors <- hue_pal()(12)
names(cluster_colors) <- 0:11

# Define the highly descriptive legend labels for the side panel
legend_labels <- c(
  "0" = "0 - Cancer 1",
  "1" = "1 - Fibroblast",
  "2" = "2 - Healthy Ductal",
  "3" = "3 - Lymphoid B Cells",
  "4" = "4 - Myeloid",
  "5" = "5 - Acinar",
  "6" = "6 - Endocrine / Islet",
  "7" = "7 - Endothelial",
  "8" = "8 - Lymphoid T/NK Cells",
  "9" = "9 - Cancer 2",
  "10" = "10 - Pericytes / Mural",
  "11" = "11 - Schwann Cells"
)


# ============================================================
# 3. SUBSET THE DATA BY SEX
# ============================================================
female_subset <- subset(Merged_Joined, subset = Sex == "Female")
male_subset <- subset(Merged_Joined, subset = Sex == "Male")


# ============================================================
# 4. GENERATE FEMALE UMAP
# ============================================================
female_umap <- DimPlot(female_subset, reduction = "umap", label = TRUE, label.size = 5, pt.size = 0.5) + 
  scale_color_manual(values = cluster_colors, labels = legend_labels) + 
  ggtitle("UMAP Projection: Female PDAC Samples") +
  theme_classic(base_family = "sans") + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.position = "right"
  )


# ============================================================
# 5. GENERATE MALE UMAP
# ============================================================
male_umap <- DimPlot(male_subset, reduction = "umap", label = TRUE, label.size = 5, pt.size = 0.5) + 
  scale_color_manual(values = cluster_colors, labels = legend_labels) + 
  ggtitle("UMAP Projection: Male PDAC Samples") +
  theme_classic(base_family = "sans") + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.position = "right"
  )


# ============================================================
# 6. PRINT AND COMBINE THE PLOTS
# ============================================================
# Stitch them together side-by-side for a unified, publication-ready figure
combined_plot <- female_umap + male_umap + plot_layout(guides = "collect")
print(combined_plot)
# ============================================================
# SINGLE UMAP: COLORED BY SEX (FEMALE = RED, MALE = BLUE)
# ============================================================
sex_umap <- DimPlot(Merged_Joined, reduction = "umap", group.by = "Sex", pt.size = 0.5) + 
  # Assign specific colors to the exact metadata text
  scale_color_manual(values = c("Female" = "red", "Male" = "blue")) + 
  ggtitle("UMAP Projection by Patient Sex") +
  theme_classic(base_family = "sans") + 
  theme(
    # Professional typography and centering
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.text = element_text(size = 11),
    legend.title = element_blank(), # Removes the "Sex" title above the legend for a cleaner look
    legend.position = "right"
  )

# Print the plot
print(sex_umap)


# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales) 

# ============================================================
# 1. EXTRACT AND FORMAT DATA FOR PLOTTING
# ============================================================
# Extract the metadata from your Seurat object
plot_data <- Merged_Joined@meta.data

# Use the original cluster numbers so we can map them precisely to your list
plot_data$CellType <- as.factor(plot_data$seurat_clusters)

# Create a mapping key to rename the 17 samples for the publication figure
# Ensure these left-side names perfectly match your 'orig.ident' values
patient_mapping <- c(
  "Sample_01" = "P01", "Sample_02" = "P02", "Sample_03" = "P03",
  "Sample_04" = "P04", "Sample_05" = "P05", "Sample_06" = "P06",
  "Sample_07" = "P07", "Sample_08" = "P08", "Sample_09" = "P09",
  "Sample_10" = "P10", "Sample_11" = "P11", "Sample_12" = "P12",
  "Sample_13" = "P13", "Sample_14" = "P14", "Sample_15" = "P15",
  "Sample_16" = "P16", "Sample_17" = "P17"
)

# Apply the mapping to create a clean 'PatientID' column
plot_data$PatientID <- patient_mapping[as.character(plot_data$orig.ident)]

# Lock the order of the new PatientIDs so they display sequentially on the graph
plot_data$PatientID <- factor(plot_data$PatientID, levels = sprintf("P%02d", 1:17))


# ============================================================
# 2. PREPARE COLORS AND LEGEND LABELS
# ============================================================
# Generate a fixed color palette based on the number of clusters
cluster_colors <- hue_pal()(length(levels(plot_data$CellType)))
names(cluster_colors) <- levels(plot_data$CellType)

# Define your exact numbered legend labels
legend_labels <- c(
  "0" = "0 - Cancer 1",
  "1" = "1 - Fibroblast",
  "2" = "2 - Healthy Ductal",
  "3" = "3 - Lymphoid B Cells",
  "4" = "4 - Myeloid",
  "5" = "5 - Acinar",
  "6" = "6 - Endocrine / Islet",
  "7" = "7 - Endothelial",
  "8" = "8 - Lymphoid T Cells / NK Cells",
  "9" = "9 - Cancer 2",
  "10" = "10 - Pericytes / Smooth Muscle / Mural Cells",
  "11" = "11 - Schwann Cells"
)


# ============================================================
# 3. GENERATE THE STACKED PROPORTION BAR CHART
# ============================================================
# Note that we are now using x = PatientID instead of orig.ident
proportion_plot <- ggplot(plot_data, aes(x = PatientID, fill = CellType)) +
  # geom_bar with position="fill" calculates proportions (0 to 1.00)
  # color="white" adds a microscopic border between segments for visual clarity
  geom_bar(position = "fill", color = "white", linewidth = 0.2) + 
  
  # Apply the assigned colors and overwrite the legend with your custom labels
  scale_fill_manual(values = cluster_colors, labels = legend_labels) +
  
  # Ensure the y-axis strictly follows a 0.00 to 1.00 proportion scale
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), expand = c(0, 0)) +
  
  # Professional Titles and Labels
  ggtitle("Cell Type Proportions per Patient Sample") +
  labs(x = "Patient Sample", y = "Proportion of Cells", fill = "") +
  
  # Publication-standard formatting
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    
    # X-axis: Angled sample names for clarity
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 10)),
    
    # Y-axis formatting
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 10)),
    
    # Bottom legend formatting
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm")
  ) +
  # Format legend into 4 neat rows to accommodate long string lengths
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))

# Print the final plot
print(proportion_plot)


# Load required libraries
library(Seurat)
library(ggplot2)

# ============================================================
# 1. ORGANIZE YOUR GENES OF INTEREST
# ============================================================
# Combine all markers into a single vector. 
# Keeping them in this specific order ensures they group together on the x-axis.
caf_markers <- c(
  # myCAF markers
  "ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2", "POSTN",
  # iCAF markers
  "IL6", "CXCL12", "CXCL14", "LIF", "CFD", "PDGFRA",
  # apCAF markers
  "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA"
)

fibroblast_subset <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Subset.RDS")

# ============================================================
# 2. GENERATE PUBLICATION-READY DOTPLOT
# ============================================================
pub_dotplot <- DotPlot(fibroblast_subset, features = caf_markers) + 
  
  # Set a high-contrast color gradient (light grey for 0, navy for high expression)
  scale_color_gradient(low = "lightgrey", high = "navy") +
  
  # Apply a clean, white-background theme favored by journals
  theme_classic(base_family = "sans") + 
  
  # Fine-tune the typography and layout
  theme(
    # Angle gene names to 45 degrees, italicize them, and align them properly
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", color = "black", size = 11),
    
    # Make cluster labels bold and clean
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    
    # Bold axis titles
    axis.title = element_text(face = "bold", size = 13),
    
    # Center and bold the main title
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    
    # Format the legend for readability
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)
  ) +
  
  # Add descriptive titles and labels
  labs(
    title = "Expression of Fibroblast Subpopulation Markers",
    x = "CAF Marker Genes",
    y = "Cluster Identity"
  )

# ============================================================
# 3. DISPLAY AND SAVE
# ============================================================
# Print to the viewer
print(pub_dotplot)

pub_dotplot + 
  geom_vline(xintercept = c(6.5, 12.5), linetype = "dashed", color = "gray50", alpha = 0.5)
# Optional: Save as a high-resolution PDF for your manuscript
ggsave("CAF_Subset_DotPlot.pdf", plot = pub_dotplot, width = 10, height = 5, dpi = 300)

# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)

# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)

# ============================================================
# 1. INJECT CLINICAL METADATA INTO THE SUBSET
# ============================================================
# We must re-define the mapping for this specific object since it was loaded fresh from an RDS
sex_mapping <- c(
  "Sample_01" = "Male",   "Sample_02" = "Male",   "Sample_03" = "Male",
  "Sample_04" = "Female", "Sample_05" = "Female", "Sample_06" = "Female",
  "Sample_07" = "Male",   "Sample_08" = "Female", "Sample_09" = "Male",
  "Sample_10" = "Female", "Sample_11" = "Female", "Sample_12" = "Male",
  "Sample_13" = "Female", "Sample_14" = "Male",   "Sample_15" = "Female",
  "Sample_16" = "Male",   "Sample_17" = "Female"
)

# Inject the Sex column into the fibroblast object
fibroblast_subset$Sex <- unname(sex_mapping[as.character(fibroblast_subset$orig.ident)])


# ============================================================
# 2. SUBSET THE FIBROBLASTS BY SEX
# ============================================================
# Now that the 'Sex' column exists, this will work perfectly!
female_fibros <- subset(fibroblast_subset, subset = Sex == "Female")
male_fibros <- subset(fibroblast_subset, subset = Sex == "Male")


# ============================================================
# 3. PREPARE UNIFIED COLORS
# ============================================================
# Extract the active identities (Fixed: changed fibroblast_subset_clean to fibroblast_subset)
fibro_idents <- levels(Idents(fibroblast_subset))

# Generate a fixed color palette so the same cell type is the exact same color on both graphs
fibro_colors <- hue_pal()(length(fibro_idents))
names(fibro_colors) <- fibro_idents


# ============================================================
# 4. GENERATE FEMALE FIBROBLAST UMAP
# ============================================================
female_umap <- DimPlot(female_fibros, reduction = "sub_umap", label = TRUE, label.size = 5, pt.size = 0.5) + 
  scale_color_manual(values = fibro_colors) + 
  ggtitle("Female PDAC Fibroblasts") +
  theme_classic(base_family = "sans") + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.position = "none" # We will hide this legend and rely on a shared one later
  )


# ============================================================
# 5. GENERATE MALE FIBROBLAST UMAP
# ============================================================
male_umap <- DimPlot(male_fibros, reduction = "sub_umap", label = TRUE, label.size = 5, pt.size = 0.5) + 
  scale_color_manual(values = fibro_colors) + 
  ggtitle("Male PDAC Fibroblasts") +
  theme_classic(base_family = "sans") + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.position = "right" # Keep this legend so it serves both plots
  )


# ============================================================
# 6. STITCH AND PRINT THE FINAL FIGURE
# ============================================================
# Use patchwork to combine them side-by-side and collect the legend cleanly on the right
combined_fibro_plot <- female_umap + male_umap + plot_layout(guides = "collect")

# Print the final publication-ready plot
print(combined_fibro_plot)

#Differential Gene Expression Male vs Female
# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
# Install ggrepel if you don't have it: install.packages("ggrepel")
library(ggrepel) # This is the secret to clean, non-overlapping text labels on graphs!

# ============================================================
# 1. PREPARE SEURAT OBJECT FOR DGE
# ============================================================
# CRITICAL: Always perform differential expression on the RNA assay, not the integrated assay!
DefaultAssay(fibroblast_subset) <- "RNA"

# Optionally, ensure your RNA data is joined and normalized (required for Seurat v5)
fibroblast_subset <- JoinLayers(fibroblast_subset)


# ============================================================
# 2. RUN DIFFERENTIAL GENE EXPRESSION (MALE VS FEMALE)
# ============================================================
# ident.1 = Male, ident.2 = Female
# Note: Positive Log2FC values mean the gene is higher in MALES. 
#       Negative Log2FC values mean the gene is higher in FEMALES.
sex_markers <- FindMarkers(
  object = fibroblast_subset,
  ident.1 = "Male",
  ident.2 = "Female",
  group.by = "Sex",
  logfc.threshold = 0.25, # Only look at genes with at least a 0.25 log-fold change
  min.pct = 0.10          # The gene must be expressed in at least 10% of the cells
)

# Move the gene names from the row names into their own dedicated column for easier plotting
sex_markers$Gene <- rownames(sex_markers)


# ============================================================
# 3. PREPARE DATA FOR THE VOLCANO PLOT
# ============================================================
# Create a new column to label whether a gene is significant, and in which direction
sex_markers$Significance <- "Not Significant"

# Define our strict publication thresholds (Adjust these if your data is very sparse!)
p_val_cutoff <- 0.05
log2fc_cutoff <- 0.5 

# Apply the labels based on the math
sex_markers$Significance[sex_markers$p_val_adj < p_val_cutoff & sex_markers$avg_log2FC > log2fc_cutoff] <- "Upregulated in Males"
sex_markers$Significance[sex_markers$p_val_adj < p_val_cutoff & sex_markers$avg_log2FC < -log2fc_cutoff] <- "Upregulated in Females"

# Convert to factor so we can control the exact colors in ggplot
sex_markers$Significance <- factor(sex_markers$Significance, 
                                   levels = c("Upregulated in Females", "Not Significant", "Upregulated in Males"))

# Identify the Top 15 most significant genes in both directions to label on the graph
top_genes <- sex_markers %>%
  filter(Significance != "Not Significant") %>%
  group_by(Significance) %>%
  slice_min(order_by = p_val_adj, n = 15) %>% # Grabs the 15 lowest p-values per group
  ungroup()


# ============================================================
# 4. PRINT METRICS & SAVE TABLE
# ============================================================
print("======================================================")
print("DIFFERENTIAL EXPRESSION: MALE VS FEMALE FIBROBLASTS")
print("======================================================")
print(table(sex_markers$Significance))

# Save the full mathematical results to a CSV for your supplemental data!
write.csv(sex_markers, "/Users/sudharshiniram/Desktop/SRI Data Analysis/Male_vs_Female_Fibroblasts_DGE.csv", row.names = FALSE)


# ============================================================
# 5. GENERATE PUBLICATION-READY VOLCANO PLOT
# ============================================================
volcano_plot <- ggplot(sex_markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  
  # Plot all the genes as tiny dots
  geom_point(alpha = 0.8, size = 1.5) +
  
  # Apply standard Red (Female) and Blue (Male) colors, with Grey for insignificant genes
  scale_color_manual(values = c("Upregulated in Females" = "#D55E00", 
                                "Not Significant" = "grey80", 
                                "Upregulated in Males" = "#0072B2")) +
  
  # Add dashed threshold lines so the reader can see your cutoff criteria
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(p_val_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +
  
  # Add the text labels for the top genes using ggrepel so they organically avoid overlapping
  geom_text_repel(data = top_genes, aes(label = Gene), 
                  size = 4, fontface = "italic", color = "black", 
                  box.padding = 0.5, max.overlaps = 20) +
  
  # Titles and Labels
  ggtitle("Differential Expression: Male vs. Female Fibroblasts") +
  labs(x = "Log2 Fold Change", y = "-Log10(Adjusted P-Value)") +
  
  # Publication-standard clean theme
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )

# Print the beautiful plot!
print(volcano_plot)

# Load required libraries
library(Seurat)
library(ggplot2)

Merged_Joined <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Merged_Joined.RDS")

# ============================================================
# 1. RENAME CLUSTERS TO BIOLOGICAL IDENTITIES
# ============================================================
# Define the exact mapping from numerical clusters to biological names
new_cluster_ids <- c(
  "0" = "Cancer 1",
  "1" = "Fibroblast",
  "2" = "Ductal",
  "3" = "B Cells",
  "4" = "Myeloid",
  "5" = "Acinar",
  "6" = "Endocrine",
  "7" = "Endothelial",
  "8" = "T/NK Cells",
  "9" = "Cancer 2",
  "10" = "Smooth Muscle",
  "11" = " ",
  "12" = "Adipocytes",
  "13" = " ",
  "14" = " "
)

# Ensure the active identities are set to your numerical clusters first
Idents(Merged_Joined) <- "seurat_clusters" 

# Apply the formal biological names to the Seurat object
Merged_Joined <- RenameIdents(Merged_Joined, new_cluster_ids)


# ============================================================
# 2. GENERATE PUBLICATION-READY UMAP
# ============================================================
pub_umap <- DimPlot(Merged_Joined, reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.5, repel = TRUE) + 
  
  # Professional Titles
  ggtitle("UMAP Projection of Cell Types") +
  
  # Publication-standard clean theme
  theme_classic(base_family = "sans") + 
  
  # Fine-tune typography and layout
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    
    # Legend formatting
    legend.text = element_text(size = 11),
    legend.title = element_blank(), # Hides the legend title for a cleaner look
    legend.position = "right"
  )

# ============================================================
# 3. DISPLAY THE PLOT
# ============================================================
print(pub_umap)

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)

# ============================================================
# 1. RENAME CLUSTERS TO YOUR EXACT CAF LABELS
# ============================================================
# Ensure we are looking at the numerical clusters first
Idents(fibroblast_subset) <- "seurat_clusters" 

# Define your exact biological mapping
caf_cluster_ids <- c(
  "0" = "Fibroblasts",
  "1" = "myCAF",
  "2" = "iCAF",
  "3" = "Pancreatic Stellate",
  "4" = "apCAF"
)

# Apply the labels to the Seurat object
fibroblast_subset <- RenameIdents(fibroblast_subset, caf_cluster_ids)


# ============================================================
# 2. EXTRACT METADATA AND INJECT MAPPINGS
# ============================================================
# Extract metadata from the active fibroblast subset
plot_data <- fibroblast_subset@meta.data

# Capture the newly named active identities (myCAF, iCAF, etc.)
plot_data$CellType <- Idents(fibroblast_subset)

# Apply our standard P01-P17 patient naming mapping
patient_mapping <- c(
  "Sample_01" = "P01", "Sample_02" = "P02", "Sample_03" = "P03",
  "Sample_04" = "P04", "Sample_05" = "P05", "Sample_06" = "P06",
  "Sample_07" = "P07", "Sample_08" = "P08", "Sample_09" = "P09",
  "Sample_10" = "P10", "Sample_11" = "P11", "Sample_12" = "P12",
  "Sample_13" = "P13", "Sample_14" = "P14", "Sample_15" = "P15",
  "Sample_16" = "P16", "Sample_17" = "P17"
)
plot_data$PatientID <- patient_mapping[as.character(plot_data$orig.ident)]
plot_data$PatientID <- factor(plot_data$PatientID, levels = sprintf("P%02d", 1:17))

# Inject the Sex mapping directly into the plotting data frame
sex_mapping <- c(
  "Sample_01" = "Male",   "Sample_02" = "Male",   "Sample_03" = "Male",
  "Sample_04" = "Female", "Sample_05" = "Female", "Sample_06" = "Female",
  "Sample_07" = "Male",   "Sample_08" = "Female", "Sample_09" = "Male",
  "Sample_10" = "Female", "Sample_11" = "Female", "Sample_12" = "Male",
  "Sample_13" = "Female", "Sample_14" = "Male",   "Sample_15" = "Female",
  "Sample_16" = "Male",   "Sample_17" = "Female"
)
plot_data$Sex <- sex_mapping[as.character(plot_data$orig.ident)]


# ============================================================
# 3. LOCK IN THE COLOR PALETTE
# ============================================================
# Generate a fixed color palette so cell types stay identical across all 3 graphs
fibro_types <- levels(plot_data$CellType)
fibro_colors <- hue_pal()(length(fibro_types))
names(fibro_colors) <- fibro_types


# ============================================================
# 4. GENERATE GRAPH 1: ALL SAMPLES
# ============================================================
all_samples_plot <- ggplot(plot_data, aes(x = PatientID, fill = CellType)) +
  geom_bar(position = "fill", color = "white", linewidth = 0.2) + 
  scale_fill_manual(values = fibro_colors) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), expand = c(0, 0)) +
  ggtitle("CAF Subpopulation Proportions (All Samples)") +
  labs(x = "Patient Sample", y = "Proportion of Fibroblasts", fill = "Subpopulation") +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 10)),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 10)),
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    legend.title = element_text(face = "bold", size = 11)
  )

# Print the Master Graph
print(all_samples_plot)


# ============================================================
# 5. SPLIT DATA BY SEX FOR COMPARISON GRAPHS
# ============================================================
plot_data_female <- plot_data %>% filter(Sex == "Female")
plot_data_male <- plot_data %>% filter(Sex == "Male")


# ============================================================
# 6. GENERATE GRAPH 2: FEMALE SAMPLES ONLY
# ============================================================
female_plot <- ggplot(plot_data_female, aes(x = PatientID, fill = CellType)) +
  geom_bar(position = "fill", color = "white", linewidth = 0.2) + 
  scale_fill_manual(values = fibro_colors) +
  scale_x_discrete(drop = TRUE) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), expand = c(0, 0)) +
  ggtitle("Female Samples") +
  labs(x = "Patient Sample", y = "Proportion", fill = "Subpopulation") +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 10)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(face = "bold", size = 12)
  )


# ============================================================
# 7. GENERATE GRAPH 3: MALE SAMPLES ONLY
# ============================================================
male_plot <- ggplot(plot_data_male, aes(x = PatientID, fill = CellType)) +
  geom_bar(position = "fill", color = "white", linewidth = 0.2) + 
  scale_fill_manual(values = fibro_colors) +
  scale_x_discrete(drop = TRUE) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), expand = c(0, 0)) +
  ggtitle("Male Samples") +
  labs(x = "Patient Sample", y = "Proportion", fill = "Subpopulation") +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 10)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(face = "bold", size = 12)
  )


# ============================================================
# 8. STITCH AND PRINT THE SEX COMPARISON FIGURE
# ============================================================
# Use patchwork to place them side-by-side with a single, shared legend
sex_comparison_plot <- female_plot + male_plot + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", legend.title = element_text(face = "bold"))

# Print the final combined plot
print(sex_comparison_plot)

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)

# ============================================================
# 1. RENAME CLUSTERS TO YOUR EXACT CAF LABELS
# ============================================================
Idents(fibroblast_subset) <- "seurat_clusters" 

caf_cluster_ids <- c(
  "0" = "Fibroblasts",
  "1" = "myCAF",
  "2" = "iCAF",
  "3" = "Pancreatic Stellate",
  "4" = "apCAF"
)
fibroblast_subset <- RenameIdents(fibroblast_subset, caf_cluster_ids)


# ============================================================
# 2. EXTRACT METADATA AND INJECT MAPPINGS
# ============================================================
plot_data <- fibroblast_subset@meta.data
plot_data$CellType <- Idents(fibroblast_subset)

# Inject Patient and Sex mapping directly into the plotting data frame
patient_mapping <- c(
  "Sample_01" = "P01", "Sample_02" = "P02", "Sample_03" = "P03",
  "Sample_04" = "P04", "Sample_05" = "P05", "Sample_06" = "P06",
  "Sample_07" = "P07", "Sample_08" = "P08", "Sample_09" = "P09",
  "Sample_10" = "P10", "Sample_11" = "P11", "Sample_12" = "P12",
  "Sample_13" = "P13", "Sample_14" = "P14", "Sample_15" = "P15",
  "Sample_16" = "P16", "Sample_17" = "P17"
)
plot_data$PatientID <- patient_mapping[as.character(plot_data$orig.ident)]

sex_mapping <- c(
  "Sample_01" = "Male",   "Sample_02" = "Male",   "Sample_03" = "Male",
  "Sample_04" = "Female", "Sample_05" = "Female", "Sample_06" = "Female",
  "Sample_07" = "Male",   "Sample_08" = "Female", "Sample_09" = "Male",
  "Sample_10" = "Female", "Sample_11" = "Female", "Sample_12" = "Male",
  "Sample_13" = "Female", "Sample_14" = "Male",   "Sample_15" = "Female",
  "Sample_16" = "Male",   "Sample_17" = "Female"
)
plot_data$Sex <- sex_mapping[as.character(plot_data$orig.ident)]


# ============================================================
# 3. LOCK IN THE COLOR PALETTE
# ============================================================
fibro_types <- levels(plot_data$CellType)
fibro_colors <- hue_pal()(length(fibro_types))
names(fibro_colors) <- fibro_types


# ============================================================
# 4. STATISTICAL MATH: CALCULATE TRUE AVERAGES
# ============================================================
# Step A: Calculate exact proportions per patient (prevents missing zero-counts from breaking math)
patient_props <- as.data.frame(prop.table(table(plot_data$PatientID, plot_data$CellType), margin = 1))
colnames(patient_props) <- c("PatientID", "CellType", "Proportion")

# Step B: Attach the Sex data back to these patient proportions
unique_patients <- unique(plot_data[, c("PatientID", "Sex")])
patient_props <- merge(patient_props, unique_patients, by = "PatientID")

# Step C: Calculate the rigorous average proportion per Sex
avg_props <- patient_props %>%
  group_by(Sex, CellType) %>%
  summarise(Mean_Proportion = mean(Proportion), .groups = "drop") %>%
  # Lock the factor order so Female is always on the left, Male on the right
  mutate(Sex = factor(Sex, levels = c("Female", "Male")))


# ============================================================
# 5. PRINT THE EXACT AVERAGES TO THE CONSOLE
# ============================================================
print("======================================================")
print("TRUE AVERAGE FIBROBLAST PROPORTIONS BY SEX")
print("======================================================")
print(as.data.frame(avg_props))


# ============================================================
# 6. GENERATE THE SIDE-BY-SIDE GRAPH
# ============================================================
avg_sex_plot <- ggplot(avg_props, aes(x = Sex, y = Mean_Proportion, fill = CellType)) +
  
  # stat="identity" tells ggplot to use our exact mathematical averages, not its own counts
  geom_bar(stat = "identity", color = "white", linewidth = 0.5, width = 0.6) + 
  
  scale_fill_manual(values = fibro_colors) +
  
  # Format y-axis as clean percentages (0% to 100%)
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  
  ggtitle("Average Fibroblast Subpopulations by Sex") +
  labs(x = "Patient Sex", y = "Average Proportion (%)", fill = "Subpopulation") +
  
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.text.x = element_text(size = 13, face = "bold", color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 10)),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 10)),
    legend.position = "right",
    legend.text = element_text(size = 11),
    legend.title = element_text(face = "bold", size = 12)
  )

# Print the final Publication Graph
print(avg_sex_plot)


# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)

# ============================================================
# 1. RUN DIFFERENTIAL GENE EXPRESSION (MALE VS FEMALE)
#    ONLY FOR GENES WITH AVG EXPRESSION > 0.5
# ============================================================

print("Calculating differential expression for Fibroblasts...")

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Ensure you are on the RNA assay for accurate DGE math
DefaultAssay(fibroblast_subset) <- "RNA"

# Optional but recommended for Seurat v5 if layers are split
# fibroblast_subset <- JoinLayers(fibroblast_subset)

# Make sure Sex column only contains Male/Female
fibroblast_subset <- subset(fibroblast_subset, subset = Sex %in% c("Male", "Female"))

# Check cell counts
print(table(fibroblast_subset$Sex))

# ============================================================
# 0. LOAD REQUIRED LIBRARIES
# ============================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)


# ============================================================
# 1. RUN DIFFERENTIAL GENE EXPRESSION (MALE VS FEMALE)
#    ONLY GENES WITH AVERAGE EXPRESSION > 0.5 IN BOTH SEXES
# ============================================================

print("Calculating differential expression for Fibroblasts...")

# Ensure you are on the RNA assay for accurate DGE math
DefaultAssay(fibroblast_subset) <- "RNA"

# Optional but recommended for Seurat v5 if layers are split
# fibroblast_subset <- JoinLayers(fibroblast_subset)

# Keep only cells with clear Sex labels
fibroblast_subset <- subset(fibroblast_subset, subset = Sex %in% c("Male", "Female"))

# Check cell counts
print(table(fibroblast_subset$Sex))


# ============================================================
# 1A. IDENTIFY GENES WITH AVERAGE EXPRESSION > 0.5
#     IN BOTH MALE AND FEMALE FIBROBLASTS
# ============================================================

avg_expression_by_sex <- AverageExpression(
  fibroblast_subset,
  assays = "RNA",
  group.by = "Sex",
  slot = "data"
)$RNA

# Check column names to confirm Male and Female are present
print(colnames(avg_expression_by_sex))

# Keep genes expressed above 0.5 in BOTH male and female fibroblasts
genes_expr_gt_0.5 <- rownames(avg_expression_by_sex)[
  avg_expression_by_sex[, "Male"] > 0.5 &
    avg_expression_by_sex[, "Female"] > 0.5
]

cat("Number of genes with average expression > 0.5 in BOTH Male and Female fibroblasts:",
    length(genes_expr_gt_0.5), "\n")

head(genes_expr_gt_0.5, 20)


# ============================================================
# 1B. RUN DGE ONLY ON GENES PASSING EXPRESSION FILTER
# ============================================================

sex_markers <- FindMarkers(
  object = fibroblast_subset,
  ident.1 = "Male",
  ident.2 = "Female",
  group.by = "Sex",
  features = genes_expr_gt_0.5,
  logfc.threshold = 0.25,
  min.pct = 0.10,
  only.pos = FALSE
)

# Move gene names from row names into their own column
sex_markers$Gene <- rownames(sex_markers)

cat("Number of genes returned by FindMarkers after expression filter:",
    nrow(sex_markers), "\n")


# ============================================================
# 2. PREPARE THE DATA FOR PLOTTING
# ============================================================

# Create a new column to label significance and direction
sex_markers$Significance <- "Not Significant"

# Define strict publication thresholds
p_val_cutoff <- 0.05
log2fc_cutoff <- 0.5

# Positive log2FC = higher in male fibroblasts
# Negative log2FC = higher in female fibroblasts
sex_markers$Significance[
  sex_markers$p_val_adj < p_val_cutoff &
    sex_markers$avg_log2FC > log2fc_cutoff
] <- "Upregulated in Males"

sex_markers$Significance[
  sex_markers$p_val_adj < p_val_cutoff &
    sex_markers$avg_log2FC < -log2fc_cutoff
] <- "Upregulated in Females"

# Convert to factor to lock legend order
sex_markers$Significance <- factor(
  sex_markers$Significance,
  levels = c("Upregulated in Females", "Not Significant", "Upregulated in Males")
)

print(table(sex_markers$Significance))


# ============================================================
# 2A. HANDLE ZERO ADJUSTED P-VALUES SAFELY FOR PLOTTING
# ============================================================
# If p_val_adj is 0, -log10(0) becomes Inf and can distort/cut the graph.
# This replaces zeros with the smallest nonzero adjusted p-value.

smallest_nonzero_padj <- min(sex_markers$p_val_adj[sex_markers$p_val_adj > 0], na.rm = TRUE)

sex_markers <- sex_markers %>%
  mutate(
    p_val_adj_safe = ifelse(p_val_adj == 0, smallest_nonzero_padj, p_val_adj),
    neg_log10_padj = -log10(p_val_adj_safe)
  )


# ============================================================
# 2B. IDENTIFY TOP GENES TO LABEL
# ============================================================

top_genes <- sex_markers %>%
  filter(Significance != "Not Significant") %>%
  group_by(Significance) %>%
  slice_min(order_by = p_val_adj_safe, n = 15) %>%
  ungroup()


# ============================================================
# 3. GENERATE THE VOLCANO PLOT
# ============================================================

volcano_plot_fibro <- ggplot(
  sex_markers,
  aes(x = avg_log2FC, y = neg_log10_padj, color = Significance)
) +
  
  # Plot all genes
  geom_point(alpha = 0.8, size = 1.5) +
  
  # Modified legend labels
  scale_color_manual(
    values = c(
      "Upregulated in Females" = "#8B0000",
      "Not Significant" = "grey80",
      "Upregulated in Males" = "#0072B2"
    ),
    labels = c(
      "Upregulated in Females" = "Fibroblasts from females",
      "Not Significant" = "Not Significant",
      "Upregulated in Males" = "Fibroblasts from Males"
    )
  ) +
  
  # Threshold lines
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
  
  # Gene labels
  geom_text_repel(
    data = top_genes,
    aes(label = Gene),
    size = 4,
    fontface = "italic",
    color = "black",
    box.padding = 0.5,
    point.padding = 0.4,
    max.overlaps = 30,
    min.segment.length = 0
  ) +
  
  # Add extra room so top circles/labels are not clipped
  scale_y_continuous(
    expand = expansion(mult = c(0.03, 0.18))
  ) +
  
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  coord_cartesian(clip = "off") +
  
  # Titles and labels
  ggtitle("Fibroblast Differential Expression: Male vs. Female") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted P-Value)"
  ) +
  
  # Publication-style theme
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
    
    # Extra plot margin helps prevent labels from being cut off
    plot.margin = margin(t = 25, r = 25, b = 20, l = 20)
  )


# ============================================================
# 4. DISPLAY THE PLOT
# ============================================================

print(volcano_plot_fibro)


# ============================================================
# 5. OPTIONAL: SAVE THE PLOT
# ============================================================

ggsave(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Male_vs_Female_Volcano_ExpressionGT0.5.pdf",
  plot = volcano_plot_fibro,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Male_vs_Female_Volcano_ExpressionGT0.5.png",
  plot = volcano_plot_fibro,
  width = 9,
  height = 7,
  dpi = 300
)


# ============================================================
# 1. LOAD REQUIRED LIBRARIES
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(msigdbr)
library(fgsea)
library(tibble)
library(tidyr)
library(patchwork)


# ============================================================
# 2. FETCH BIOLOGICAL PATHWAYS
# ============================================================

print("Fetching Hallmark biological pathways...")

hallmark_sets <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

pathways_list <- split(
  x = hallmark_sets$gene_symbol,
  f = hallmark_sets$gs_name
)


# ============================================================
# 3. FILTER sex_markers FOR EXPRESSION > 0.5
# ============================================================

print("Filtering genes by average expression > 0.5 in both male and female fibroblasts...")

DefaultAssay(fibroblast_subset) <- "RNA"

# Optional for Seurat v5 if layers are split
# fibroblast_subset <- JoinLayers(fibroblast_subset)

fibroblast_subset <- subset(fibroblast_subset, subset = Sex %in% c("Male", "Female"))

avg_expression_by_sex <- AverageExpression(
  fibroblast_subset,
  assays = "RNA",
  group.by = "Sex",
  slot = "data"
)$RNA

print(colnames(avg_expression_by_sex))

genes_expr_gt_0.5 <- rownames(avg_expression_by_sex)[
  avg_expression_by_sex[, "Male"] > 0.5 &
    avg_expression_by_sex[, "Female"] > 0.5
]

cat("Number of genes with average expression > 0.5 in BOTH Male and Female fibroblasts:",
    length(genes_expr_gt_0.5), "\n")


# ============================================================
# 4. PREPARE RANKED GENE LIST FOR fgsea
# ============================================================

# Make sure sex_markers has a Gene column
if (!"Gene" %in% colnames(sex_markers)) {
  sex_markers$Gene <- rownames(sex_markers)
}

sex_markers_clean <- sex_markers %>%
  filter(Gene %in% genes_expr_gt_0.5) %>%
  filter(!is.na(p_val_adj)) %>%
  filter(!is.na(avg_log2FC))

# Handle p_val_adj = 0 safely
smallest_nonzero_padj <- min(
  sex_markers_clean$p_val_adj[sex_markers_clean$p_val_adj > 0],
  na.rm = TRUE
)

sex_markers_clean <- sex_markers_clean %>%
  mutate(
    p_val_adj_safe = ifelse(p_val_adj == 0, smallest_nonzero_padj, p_val_adj),
    rank_score = sign(avg_log2FC) * -log10(p_val_adj_safe)
  ) %>%
  filter(!is.na(rank_score), is.finite(rank_score)) %>%
  arrange(desc(rank_score))

gene_ranks <- sex_markers_clean$rank_score
names(gene_ranks) <- sex_markers_clean$Gene

gene_ranks <- gene_ranks[!duplicated(names(gene_ranks))]
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

cat("Number of genes going into fgsea:", length(gene_ranks), "\n")

head(gene_ranks, 10)
tail(gene_ranks, 10)


# ============================================================
# 5. RUN fgsea
# ============================================================

print("Running fgsea...")

set.seed(42)

fgsea_results <- fgsea(
  pathways = pathways_list,
  stats = gene_ranks,
  minSize = 10,
  maxSize = 500,
  eps = 0
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
    neg_log10_FDR = -log10(padj)
  )

print(head(fgsea_results, 20))


# ============================================================
# 6. SAVE RESULTS
# ============================================================
# ============================================================
# 6. SAVE RESULTS
# ============================================================

gsea_outdir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis/GSEA_Figures_ExpressionGT0.5"
dir.create(gsea_outdir, showWarnings = FALSE, recursive = TRUE)

# fgsea has a list-column called leadingEdge.
# write.csv cannot save list-columns, so we collapse leadingEdge genes into one string.
fgsea_results_save <- fgsea_results %>%
  as.data.frame() %>%
  mutate(
    leadingEdge = sapply(leadingEdge, paste, collapse = ";")
  )

write.csv(
  fgsea_results_save,
  file.path(gsea_outdir, "Fibroblast_Male_vs_Female_fgsea_ExpressionGT0.5_Results.csv"),
  row.names = FALSE
)


# ============================================================
# 7. FIGURE 1: NES DOTPLOT
# ============================================================

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

gsea_nes_dotplot <- ggplot(plot_gsea, aes(x = NES, y = pathway_clean)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_point(aes(size = size, color = neg_log10_FDR), alpha = 0.9) +
  scale_color_gradient(low = "lightgrey", high = "darkred", name = "-log10(FDR)") +
  labs(
    title = "Hallmark GSEA: Male vs Female PDAC Fibroblasts",
    subtitle = "Positive NES = male-enriched; negative NES = female-enriched",
    x = "Normalized Enrichment Score",
    y = "Hallmark Pathway",
    size = "Gene Set Size"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 8)),
    plot.subtitle = element_text(hjust = 0.5, size = 11, margin = margin(b = 12)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

print(gsea_nes_dotplot)

ggsave(
  file.path(gsea_outdir, "Figure_1_fgsea_NES_Dotplot_ExpressionGT0.5.pdf"),
  plot = gsea_nes_dotplot,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(gsea_outdir, "Figure_1_fgsea_NES_Dotplot_ExpressionGT0.5.png"),
  plot = gsea_nes_dotplot,
  width = 9,
  height = 7,
  dpi = 300
)


# ============================================================
# 8. FIGURE 2: BARPLOT
# ============================================================

top_male <- fgsea_results %>%
  filter(NES > 0) %>%
  arrange(padj) %>%
  head(10)

top_female <- fgsea_results %>%
  filter(NES < 0) %>%
  arrange(padj) %>%
  head(10)

top_combined <- bind_rows(top_female, top_male) %>%
  arrange(NES)

top_combined$pathway_clean <- factor(
  top_combined$pathway_clean,
  levels = top_combined$pathway_clean
)

gsea_barplot <- ggplot(top_combined, aes(x = pathway_clean, y = NES, fill = Direction)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.25) +
  coord_flip() +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = c(
    "Enriched in Female Fibroblasts" = "#8B0000",
    "Enriched in Male Fibroblasts" = "#0072B2"
  )) +
  labs(
    title = "Top Sex-Biased Pathways in PDAC Fibroblasts",
    x = "Hallmark Pathway",
    y = "Normalized Enrichment Score",
    fill = ""
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 11, face = "bold")
  )

print(gsea_barplot)

ggsave(
  file.path(gsea_outdir, "Figure_2_fgsea_Barplot_ExpressionGT0.5.pdf"),
  plot = gsea_barplot,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(gsea_outdir, "Figure_2_fgsea_Barplot_ExpressionGT0.5.png"),
  plot = gsea_barplot,
  width = 9,
  height = 7,
  dpi = 300
)


# ============================================================
# 9. FIGURE 3: ENRICHMENT CURVES
# ============================================================

top_male_pathway <- fgsea_results %>%
  filter(NES > 0) %>%
  arrange(padj) %>%
  slice(1) %>%
  pull(pathway)

top_female_pathway <- fgsea_results %>%
  filter(NES < 0) %>%
  arrange(padj) %>%
  slice(1) %>%
  pull(pathway)

if (length(top_male_pathway) > 0) {
  male_curve <- plotEnrichment(
    pathways_list[[top_male_pathway]],
    gene_ranks
  ) +
    labs(
      title = paste("Male-Enriched:", gsub("_", " ", gsub("HALLMARK_", "", top_male_pathway))),
      x = "Ranked Genes",
      y = "Enrichment Score"
    ) +
    theme_classic(base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      axis.title = element_text(face = "bold")
    )
  
  print(male_curve)
  
  ggsave(
    file.path(gsea_outdir, "Figure_3A_Top_Male_Enrichment_Curve_ExpressionGT0.5.pdf"),
    plot = male_curve,
    width = 8,
    height = 5,
    dpi = 300
  )
}

if (length(top_female_pathway) > 0) {
  female_curve <- plotEnrichment(
    pathways_list[[top_female_pathway]],
    gene_ranks
  ) +
    labs(
      title = paste("Female-Enriched:", gsub("_", " ", gsub("HALLMARK_", "", top_female_pathway))),
      x = "Ranked Genes",
      y = "Enrichment Score"
    ) +
    theme_classic(base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      axis.title = element_text(face = "bold")
    )
  
  print(female_curve)
  
  ggsave(
    file.path(gsea_outdir, "Figure_3B_Top_Female_Enrichment_Curve_ExpressionGT0.5.pdf"),
    plot = female_curve,
    width = 8,
    height = 5,
    dpi = 300
  )
}


# ============================================================
# 10. FINAL SUMMARY
# ============================================================

cat("\n======================================================\n")
cat("fgsea completed successfully.\n")
cat("Do not use enrichIt() unless you are using the escape package.\n")
cat("Figures saved to:\n")
cat(gsea_outdir, "\n")
cat("======================================================\n")
# Load required libraries
# Install if needed: install.packages("ggVennDiagram")
install.packages("ggVennDiagram")
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggVennDiagram)

# ============================================================
# 1. EXTRACT GENERAL FIBROBLAST GENE LISTS
# ============================================================
print("Extracting General Fibroblast signatures...")
# Using the 'sex_markers' data frame we already generated in Step 3a
gen_male_up <- sex_markers %>% filter(Significance == "Upregulated in Males") %>% pull(Gene)
gen_female_up <- sex_markers %>% filter(Significance == "Upregulated in Females") %>% pull(Gene)


# ============================================================
# 2. RUN DGE FOR iCAF SUBSET ONLY
# ============================================================
print("Running Male vs Female DGE for iCAFs...")

# Isolate just the iCAFs
icaf_subset <- subset(fibroblast_subset, idents = "iCAF")

# Run the math
icaf_markers <- FindMarkers(icaf_subset, ident.1 = "Male", ident.2 = "Female", group.by = "Sex", 
                            logfc.threshold = 0.25, min.pct = 0.05)

# Extract the significant genes based on our strict cutoffs
icaf_male_up <- rownames(icaf_markers[icaf_markers$avg_log2FC > 0.5 & icaf_markers$p_val_adj < 0.05, ])
icaf_female_up <- rownames(icaf_markers[icaf_markers$avg_log2FC < -0.5 & icaf_markers$p_val_adj < 0.05, ])


# ============================================================
# 3. RUN DGE FOR myCAF SUBSET ONLY
# ============================================================
print("Running Male vs Female DGE for myCAFs...")

# Isolate just the myCAFs
mycaf_subset <- subset(fibroblast_subset, idents = "myCAF")

# Run the math
mycaf_markers <- FindMarkers(mycaf_subset, ident.1 = "Male", ident.2 = "Female", group.by = "Sex", 
                             logfc.threshold = 0.25, min.pct = 0.05)

# Extract the significant genes
mycaf_male_up <- rownames(mycaf_markers[mycaf_markers$avg_log2FC > 0.5 & mycaf_markers$p_val_adj < 0.05, ])
mycaf_female_up <- rownames(mycaf_markers[mycaf_markers$avg_log2FC < -0.5 & mycaf_markers$p_val_adj < 0.05, ])


# ============================================================
# 4. GENERATE THE REQUESTED VENN DIAGRAMS (MALE VS FEMALE)
# ============================================================
# Note: The intersection in these three graphs will mathematically be 0.

# 1. General Fibroblasts
venn_gen <- ggVennDiagram(list("Male UP" = gen_male_up, "Female UP" = gen_female_up), 
                          set_color = c("#0072B2", "#8B0000")) +
  scale_fill_gradient(low = "white", high = "grey80") +
  ggtitle("General Fibroblasts: Male vs Female Upregulated") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

# 2. iCAFs Only
venn_icaf <- ggVennDiagram(list("Male UP" = icaf_male_up, "Female UP" = icaf_female_up), 
                           set_color = c("#0072B2", "#8B0000")) +
  scale_fill_gradient(low = "white", high = "grey80") +
  ggtitle("iCAFs: Male vs Female Upregulated") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

# 3. myCAFs Only
venn_mycaf <- ggVennDiagram(list("Male UP" = mycaf_male_up, "Female UP" = mycaf_female_up), 
                            set_color = c("#0072B2", "#8B0000")) +
  scale_fill_gradient(low = "white", high = "grey80") +
  ggtitle("myCAFs: Male vs Female Upregulated") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

print(venn_gen)
print(venn_icaf)
print(venn_mycaf)


# ============================================================
# 5. THE PUBLICATION GRAPH: CROSS-SUBTYPE CONSERVED SIGNATURES
# ============================================================
# This compares Male genes in iCAFs vs Male genes in myCAFs. 
# The intersection here is biologically highly meaningful!

venn_conserved_male <- ggVennDiagram(
  list("iCAF Male Upregulated" = icaf_male_up, 
       "myCAF Male Upregulated" = mycaf_male_up),
  set_color = c("black", "black") # Clean black borders
) +
  # Use a blue gradient to visually denote we are looking at Male genes
  scale_fill_gradient(low = "white", high = "#0072B2") +
  ggtitle("Conserved Male Upregulated Signature Across CAFs") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

print(venn_conserved_male)




# ============================================================
# 0. LOAD LIBRARIES
# ============================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork) # For stitching figures side-by-side

# Clean up any potential NA values in the Sex column just in case
fibroblast_subset <- subset(fibroblast_subset, subset = Sex %in% c("Male", "Female"))


# ============================================================
# 1. RUN DGE FOR iCAF SUBSET ONLY
# ============================================================
print("Extracting iCAF sex-biased genes...")
icaf_subset <- subset(fibroblast_subset, idents = "iCAF")

icaf_markers <- FindMarkers(
  object = icaf_subset,
  ident.1 = "Male",
  ident.2 = "Female",
  group.by = "Sex",
  logfc.threshold = 0.25, 
  min.pct = 0.05          # Drop to 5% for snRNA-seq depth
)
icaf_markers$Gene <- rownames(icaf_markers)


# ============================================================
# 2. RUN DGE FOR myCAF SUBSET ONLY
# ============================================================
print("Extracting myCAF sex-biased genes...")
mycaf_subset <- subset(fibroblast_subset, idents = "myCAF")

mycaf_markers <- FindMarkers(
  object = mycaf_subset,
  ident.1 = "Male",
  ident.2 = "Female",
  group.by = "Sex",
  logfc.threshold = 0.25, 
  min.pct = 0.05          # Drop to 5% for snRNA-seq depth
)
mycaf_markers$Gene <- rownames(mycaf_markers)


# ============================================================
# 3. EXTRACT THE GENE LISTS FOR DATA ANALYSIS
# ============================================================
# Define publication thresholds
p_val_cutoff <- 0.05
log2fc_cutoff <- 0.5

# iCAF lists
icaf_male_up_genes   <- icaf_markers %>% filter(p_val_adj < p_val_cutoff & avg_log2FC > log2fc_cutoff) %>% pull(Gene)
icaf_female_up_genes <- icaf_markers %>% filter(p_val_adj < p_val_cutoff & avg_log2FC < -log2fc_cutoff) %>% pull(Gene)

# myCAF lists
mycaf_male_up_genes   <- mycaf_markers %>% filter(p_val_adj < p_val_cutoff & avg_log2FC > log2fc_cutoff) %>% pull(Gene)
mycaf_female_up_genes <- mycaf_markers %>% filter(p_val_adj < p_val_cutoff & avg_log2FC < -log2fc_cutoff) %>% pull(Gene)


# ============================================================
# 4. PRINT SUMMARY TO CONSOLE (Copy this data for your records)
# ============================================================
cat("\n======================================================\n")
cat("SUMMARY OF SEX-UPREGULATED GENES (Threshold log2FC > 0.5)\n")
cat("======================================================\n")
cat("iCAF Male-Upregulated Gene Count:  ", length(icaf_male_up_genes), "\n")
cat("iCAF Female-Upregulated Gene Count:", length(icaf_female_up_genes), "\n")
cat("myCAF Male-Upregulated Gene Count: ", length(mycaf_male_up_genes), "\n")
cat("myCAF Female-Upregulated Gene Count:", length(mycaf_female_up_genes), "\n")
cat("======================================================\n")


# ============================================================
# 5. GENERATE COMPLEMENTARY VOLCANO PLOTS
# ============================================================
# Helper function to generate clean, matched plots for your figure panels
make_vplot <- NA
make_vplot <- function(df_markers, title_text) {
  df_markers$Significance <- "Not Significant"
  df_markers$Significance[df_markers$p_val_adj < p_val_cutoff & df_markers$avg_log2FC > log2fc_cutoff] <- "Upregulated in Males"
  df_markers$Significance[df_markers$p_val_adj < p_val_cutoff & df_markers$avg_log2FC < -log2fc_cutoff] <- "Upregulated in Females"
  df_markers$Significance <- factor(df_markers$Significance, levels = c("Upregulated in Females", "Not Significant", "Upregulated in Males"))
  
  ggplot(df_markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
    geom_point(alpha = 0.7, size = 1.2) +
    scale_color_manual(values = c("Upregulated in Females" = "#8B0000", "Not Significant" = "grey85", "Upregulated in Males" = "#0072B2")) +
    geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", alpha = 0.4) +
    geom_hline(yintercept = -log10(p_val_cutoff), linetype = "dashed", alpha = 0.4) +
    ggtitle(title_text) + labs(x = "Log2 Fold Change", y = "-Log10(Adj P-Value)") +
    theme_classic(base_family = "sans") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), legend.position = "none")
}

# Construct the individual panels
panel_icaf  <- make_vplot(icaf_markers, "iCAF Subpopulation (M vs F)")
panel_mycaf <- make_vplot(mycaf_markers, "myCAF Subpopulation (M vs F)")

# Stitch the figures side-by-side with a collected legend at the bottom
combined_figure <- panel_icaf + panel_mycaf + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(face = "bold"))


# ============================================================
# 6. DISPLAY PLOT
# ============================================================
print(combined_figure)

# ============================================================
# GSEA: MALE VS FEMALE FIBROBLASTS IN HUMAN PDAC snRNA-seq
# ============================================================

# -------------------------------
# 0. INSTALL / LOAD LIBRARIES
# -------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages_cran <- c("dplyr", "ggplot2", "msigdbr", "tibble")
packages_bioc <- c("clusterProfiler", "enrichplot", "org.Hs.eg.db")

for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

for (pkg in packages_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}

library(Seurat)
library(dplyr)
library(ggplot2)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(tibble)

packages_cran <- c("dplyr", "ggplot2", "msigdbr", "tibble")
packages_bioc <- c("clusterProfiler", "enrichplot", "org.Hs.eg.db")

for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

for (pkg in packages_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}

# ============================================================
# 1. PREPARE FIBROBLAST OBJECT
# ============================================================

DefaultAssay(fibroblast_subset) <- "RNA"

# For Seurat v5 layered objects
fibroblast_subset <- JoinLayers(fibroblast_subset)

# Re-inject Sex metadata if needed
sex_mapping <- c(
  "Sample_01" = "Male",   "Sample_02" = "Male",   "Sample_03" = "Male",
  "Sample_04" = "Female", "Sample_05" = "Female", "Sample_06" = "Female",
  "Sample_07" = "Male",   "Sample_08" = "Female", "Sample_09" = "Male",
  "Sample_10" = "Female", "Sample_11" = "Female", "Sample_12" = "Male",
  "Sample_13" = "Female", "Sample_14" = "Male",   "Sample_15" = "Female",
  "Sample_16" = "Male",   "Sample_17" = "Female"
)

fibroblast_subset$Sex <- unname(sex_mapping[as.character(fibroblast_subset$orig.ident)])

# Remove any cells without sex annotation
fibroblast_subset <- subset(fibroblast_subset, subset = Sex %in% c("Male", "Female"))

# Check sample/cell balance
print(table(fibroblast_subset$orig.ident, fibroblast_subset$Sex))
print(table(fibroblast_subset$Sex))


# ============================================================
# 2. RUN DIFFERENTIAL EXPRESSION: MALE VS FEMALE
# ============================================================

sex_markers <- FindMarkers(
  object = fibroblast_subset,
  ident.1 = "Male",
  ident.2 = "Female",
  group.by = "Sex",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.05,
  only.pos = FALSE
)

sex_markers$Gene <- rownames(sex_markers)

# Save DGE table
write.csv(
  sex_markers,
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Male_vs_Female_DGE_for_GSEA.csv",
  row.names = FALSE
)


# ============================================================
# 3. CREATE RANKED GENE LIST FOR GSEA
# ============================================================

# For GSEA, every gene needs one ranking score.
# Positive score = male-upregulated.
# Negative score = female-upregulated.
# This ranking uses both fold-change direction and significance.

sex_markers_ranked <- sex_markers %>%
  filter(!is.na(avg_log2FC)) %>%
  filter(!is.na(p_val_adj)) %>%
  mutate(
    p_val_adj_safe = ifelse(p_val_adj == 0, min(p_val_adj[p_val_adj > 0], na.rm = TRUE), p_val_adj),
    rank_score = sign(avg_log2FC) * -log10(p_val_adj_safe)
  ) %>%
  filter(!is.na(rank_score), is.finite(rank_score)) %>%
  arrange(desc(rank_score))

# Create named ranked vector
gene_ranks <- sex_markers_ranked$rank_score
names(gene_ranks) <- sex_markers_ranked$Gene

# Remove duplicated genes if any
gene_ranks <- gene_ranks[!duplicated(names(gene_ranks))]

# Sort decreasing, required by clusterProfiler::GSEA
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

# Quick check
head(gene_ranks, 10)   # strongest male-up genes
tail(gene_ranks, 10)   # strongest female-up genes


# ============================================================
# 4. LOAD MSIGDB HALLMARK PATHWAYS
# ============================================================

hallmark_df <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

hallmark_t2g <- hallmark_df %>%
  dplyr::select(gs_name, gene_symbol)

# Optional: clean pathway names later for plotting
head(hallmark_t2g)


# ============================================================
# 5. RUN GSEA
# ============================================================

gsea_hallmark <- GSEA(
  geneList = gene_ranks,
  TERM2GENE = hallmark_t2g,
  pvalueCutoff = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 0,
  verbose = TRUE
)

gsea_results <- as.data.frame(gsea_hallmark)

# Clean pathway names
gsea_results <- gsea_results %>%
  mutate(
    pathway_clean = gsub("HALLMARK_", "", ID),
    pathway_clean = gsub("_", " ", pathway_clean),
    Direction = ifelse(NES > 0, "Enriched in Male Fibroblasts", "Enriched in Female Fibroblasts")
  ) %>%
  arrange(p.adjust)

# View top pathways
print(head(gsea_results, 20))

# Save results
write.csv(
  gsea_results,
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Male_vs_Female_Hallmark_GSEA.csv",
  row.names = FALSE
)


# ============================================================
# 6. PUBLICATION-STYLE DOTPLOT
# ============================================================

# Keep significant pathways, or top 20 if few are significant
plot_gsea <- gsea_results %>%
  filter(p.adjust < 0.25) %>%
  arrange(NES)

if (nrow(plot_gsea) < 5) {
  plot_gsea <- gsea_results %>%
    arrange(p.adjust) %>%
    head(20) %>%
    arrange(NES)
}

plot_gsea$pathway_clean <- factor(plot_gsea$pathway_clean, levels = plot_gsea$pathway_clean)

gsea_dotplot <- ggplot(plot_gsea, aes(x = NES, y = pathway_clean)) +
  geom_point(aes(size = setSize, color = p.adjust)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted P-value") +
  labs(
    title = "Hallmark GSEA: Male vs Female PDAC Fibroblasts",
    x = "Normalized Enrichment Score",
    y = "Hallmark Pathway",
    size = "Gene Set Size"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

print(gsea_dotplot)

ggsave(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Male_vs_Female_Hallmark_GSEA_Dotplot.pdf",
  plot = gsea_dotplot,
  width = 9,
  height = 7,
  dpi = 300
)


# ============================================================
# 7. BARPLOT OF TOP MALE AND FEMALE PATHWAYS
# ============================================================

top_male <- gsea_results %>%
  filter(NES > 0) %>%
  arrange(p.adjust) %>%
  head(10)

top_female <- gsea_results %>%
  filter(NES < 0) %>%
  arrange(p.adjust) %>%
  head(10)

top_combined <- bind_rows(top_female, top_male) %>%
  arrange(NES)

top_combined$pathway_clean <- factor(top_combined$pathway_clean, levels = top_combined$pathway_clean)

gsea_barplot <- ggplot(top_combined, aes(x = pathway_clean, y = NES, fill = Direction)) +
  geom_col(width = 0.75) +
  coord_flip() +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = c(
    "Enriched in Female Fibroblasts" = "#8B0000",
    "Enriched in Male Fibroblasts" = "#0072B2"
  )) +
  labs(
    title = "Top Sex-Biased Pathways in PDAC Fibroblasts",
    x = "Hallmark Pathway",
    y = "Normalized Enrichment Score",
    fill = ""
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 11, face = "bold")
  )

print(gsea_barplot)

ggsave(
  "/Users/sudharshiniram/Desktop/SRI Data Analysis/Fibroblast_Male_vs_Female_Hallmark_GSEA_Barplot.pdf",
  plot = gsea_barplot,
  width = 9,
  height = 7,
  dpi = 300
)


# ============================================================
# 8. ENRICHMENT CURVES FOR TOP PATHWAYS
# ============================================================

# Pick top male-enriched and top female-enriched pathway
top_male_pathway <- gsea_results %>%
  filter(NES > 0) %>%
  arrange(p.adjust) %>%
  slice(1) %>%
  pull(ID)

top_female_pathway <- gsea_results %>%
  filter(NES < 0) %>%
  arrange(p.adjust) %>%
  slice(1) %>%
  pull(ID)

# Plot top male pathway
if (length(top_male_pathway) > 0) {
  p_male_curve <- gseaplot2(
    gsea_hallmark,
    geneSetID = top_male_pathway,
    title = paste("Male-Enriched:", gsub("_", " ", gsub("HALLMARK_", "", top_male_pathway)))
  )
  print(p_male_curve)
}

# Plot top female pathway
if (length(top_female_pathway) > 0) {
  p_female_curve <- gseaplot2(
    gsea_hallmark,
    geneSetID = top_female_pathway,
    title = paste("Female-Enriched:", gsub("_", " ", gsub("HALLMARK_", "", top_female_pathway)))
  )
  print(p_female_curve)
}
