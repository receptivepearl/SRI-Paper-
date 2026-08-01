library(Seurat)
library(dplyr)
library(ggplot2)

# 1. Expand the memory limit
Sys.setenv('R_MAX_VSIZE'=32000000000)

# Load ONLY your final, processed object
# Note: using  'test' file  created at the bottom of your last script
Merged <- readRDS("/Users/sudharshiniram/Desktop/SRI Data Analysis/Data_Analysis_processed.rds")

gc()

# Calculate the markers with the downsampling memory-saver
Merged.markers <- FindAllMarkers(
  Merged, 
  only.pos = TRUE, 
  max.cells.per.ident = 500 
)

# Save the marker spreadsheet 
saveRDS(Merged.markers, file = '/Users/sudharshiniram/Desktop/SRI Data Analysis/Marker_Genes.rds')

# Convert the R data into CSV file
write.csv(Merged.markers, file = "/Users/sudharshiniram/Desktop/SRI Data Analysis/Marker_Genes.csv")
View(Merged.markers)

# Create  list of 44 names
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

# Assign labels to current clusters
names(my_cell_types) <- levels(Merged)
Merged <- RenameIdents(Merged, my_cell_types)

# Plot the final,labeled UMAP
final_umap <- DimPlot(Merged, reduction = "umap", label = TRUE, label.size = 3.5, repel = TRUE) + 
  ggtitle("PDAC Tumor Microenvironment by Cell Type") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#save umap as file
ggsave(filename = "/Users/sudharshiniram/Desktop/SRI Data Analysis/Final_Labeled_UMAP.png", 
       plot = final_umap, 
       width = 14,      
       height = 10,     
       dpi = 300)  

# Save  finalized object
saveRDS(Merged, file = '/Users/sudharshiniram/Desktop/SRI Data Analysis/Data_Analysis_Final_Labeled.rds')

     

