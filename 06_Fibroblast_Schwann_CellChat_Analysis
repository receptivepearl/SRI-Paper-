# ============================================================
# CELLCHAT ANALYSIS:
# Fibroblast <-> Schwann Cell Interactions
# Female vs Male PDAC snRNA-seq Samples
# ============================================================


# ============================================================
# 0. LOAD REQUIRED LIBRARIES
# ============================================================

library(CellChat)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

options(stringsAsFactors = FALSE)


# ============================================================
# 1. PREPARE SEURAT OBJECT
# ============================================================

# Your main Seurat object should already be loaded as Merged_Joined
# It should contain:
# - RNA assay
# - orig.ident
# - seurat_clusters
# - true_schwann_subset
# - Sex, or we will add Sex below

DefaultAssay(Merged_Joined) <- "RNA"

# Optional for Seurat v5 if layers are split
# Merged_Joined <- JoinLayers(Merged_Joined)


# ============================================================
# 2. ADD SEX METADATA
# ============================================================

sex_mapping <- c(
  "Sample_01" = "Male",   "Sample_02" = "Male",   "Sample_03" = "Male",
  "Sample_04" = "Female", "Sample_05" = "Female", "Sample_06" = "Female",
  "Sample_07" = "Male",   "Sample_08" = "Female", "Sample_09" = "Male",
  "Sample_10" = "Female", "Sample_11" = "Female", "Sample_12" = "Male",
  "Sample_13" = "Female", "Sample_14" = "Male",   "Sample_15" = "Female",
  "Sample_16" = "Male",   "Sample_17" = "Female"
)

Merged_Joined$Sex <- unname(
  sex_mapping[as.character(Merged_Joined$orig.ident)]
)

Merged_Joined$Sex <- factor(
  Merged_Joined$Sex,
  levels = c("Female", "Male")
)


# ============================================================
# 3. CREATE CLEAN CELLCHAT CELL TYPE LABELS
# ============================================================

# Create a fresh CellChat grouping column
Merged_Joined$cellchat_celltype <- NA

# Label fibroblasts.
# Based on your previous annotation, cluster 1 = Fibroblast.
Merged_Joined$cellchat_celltype[
  as.character(Merged_Joined$seurat_clusters) == "1"
] <- "Fibroblasts"

# Label true Schwann cells using your ERBB3+/PMP22+ Schwann metadata
Merged_Joined$cellchat_celltype[
  Merged_Joined$true_schwann_subset == "True Schwann"
] <- "Schwann cells"

# Keep only Fibroblasts and Schwann cells with clear Sex labels
target_cells <- subset(
  Merged_Joined,
  subset = cellchat_celltype %in% c("Fibroblasts", "Schwann cells") &
    Sex %in% c("Female", "Male")
)

target_cells$cellchat_celltype <- factor(
  target_cells$cellchat_celltype,
  levels = c("Fibroblasts", "Schwann cells")
)

target_cells$Sex <- factor(
  target_cells$Sex,
  levels = c("Female", "Male")
)

print("Cell counts by cell type and sex:")
print(table(target_cells$cellchat_celltype, target_cells$Sex))

print("Cell counts by patient, cell type, and sex:")
print(table(target_cells$orig.ident, target_cells$cellchat_celltype, target_cells$Sex))


# ============================================================
# 4. SPLIT INTO FEMALE AND MALE OBJECTS
# ============================================================

seurat_female <- subset(target_cells, subset = Sex == "Female")
seurat_male   <- subset(target_cells, subset = Sex == "Male")

print("Female object cell counts:")
print(table(seurat_female$cellchat_celltype))

print("Male object cell counts:")
print(table(seurat_male$cellchat_celltype))


# ============================================================
# 5. SAFETY CHECK
# ============================================================

min_cells_required <- 10

female_counts <- table(seurat_female$cellchat_celltype)
male_counts   <- table(seurat_male$cellchat_celltype)

if (any(female_counts < min_cells_required)) {
  warning("Female subset has a cell type with fewer than 10 cells. Interpret CellChat cautiously.")
}

if (any(male_counts < min_cells_required)) {
  warning("Male subset has a cell type with fewer than 10 cells. Interpret CellChat cautiously.")
}


# ============================================================
# 6. CREATE OUTPUT FOLDER
# ============================================================

cellchat_outdir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis/CellChat_Fibroblast_Schwann"
dir.create(cellchat_outdir, showWarnings = FALSE, recursive = TRUE)


# ============================================================
# 7. CELLCHAT FUNCTION
# ============================================================
# Important revisions:
# - Removed projectData(cellchat, PPI.human), since your CellChat version
#   does not have projectData().
# - Added meta$samples <- orig.ident so CellChat knows patient sample IDs.
# - Uses cellchat_celltype as the grouping column.

run_cellchat <- function(seurat_obj, dataset_name, min_cells_filter = 10) {
  
  print(paste("Running CellChat for", dataset_name, "..."))
  
  DefaultAssay(seurat_obj) <- "RNA"
  
  data.input <- GetAssayData(
    seurat_obj,
    assay = "RNA",
    layer = "data"
  )
  
  meta <- seurat_obj@meta.data
  
  meta$cellchat_celltype <- droplevels(
    factor(meta$cellchat_celltype)
  )
  
  # Fixes CellChat warning about missing samples column.
  # This keeps patient/sample information instead of assigning all cells to sample1.
  if (!"samples" %in% colnames(meta)) {
    meta$samples <- as.character(meta$orig.ident)
  }
  
  print(paste("CellChat groups for", dataset_name))
  print(table(meta$cellchat_celltype))
  
  print(paste("Samples represented in", dataset_name))
  print(table(meta$samples))
  
  cellchat <- createCellChat(
    object = data.input,
    meta = meta,
    group.by = "cellchat_celltype"
  )
  
  # Use human ligand-receptor database
  CellChatDB <- CellChatDB.human
  
  # Use full database first.
  # This includes Secreted Signaling, ECM-Receptor, and Cell-Cell Contact.
  cellchat@DB <- CellChatDB
  
  # Optional if you only want secreted signaling:
  # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  # cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # Removed:
  # cellchat <- projectData(cellchat, PPI.human)
  
  cellchat <- computeCommunProb(
    cellchat,
    type = "triMean",
    raw.use = TRUE
  )
  
  cellchat <- filterCommunication(
    cellchat,
    min.cells = min_cells_filter
  )
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  return(cellchat)
}


# ============================================================
# 8. RUN CELLCHAT FOR FEMALE AND MALE
# ============================================================

cellchat_female <- run_cellchat(
  seurat_obj = seurat_female,
  dataset_name = "Female",
  min_cells_filter = 10
)

cellchat_male <- run_cellchat(
  seurat_obj = seurat_male,
  dataset_name = "Male",
  min_cells_filter = 10
)

cellchat_female <- updateCellChat(cellchat_female)
cellchat_male   <- updateCellChat(cellchat_male)


# ============================================================
# 9. SAVE INDIVIDUAL CELLCHAT OBJECTS
# ============================================================

saveRDS(
  cellchat_female,
  file.path(cellchat_outdir, "cellchat_female_fibroblast_schwann.rds")
)

saveRDS(
  cellchat_male,
  file.path(cellchat_outdir, "cellchat_male_fibroblast_schwann.rds")
)


# ============================================================
# 10. MERGE CELLCHAT OBJECTS
# ============================================================
# Female = dataset 1
# Male   = dataset 2
# comparison = c(1, 2) should be interpreted as Male relative to Female
# in differential visualizations.

object.list <- list(
  Female = cellchat_female,
  Male = cellchat_male
)

cellchat_merged <- mergeCellChat(
  object.list,
  add.names = names(object.list)
)

saveRDS(
  cellchat_merged,
  file.path(cellchat_outdir, "cellchat_merged_female_male_fibroblast_schwann.rds")
)


# ============================================================
# 11. EXTRACT COMMUNICATION TABLES
# ============================================================

df_female <- subsetCommunication(cellchat_female)
df_male   <- subsetCommunication(cellchat_male)

df_female$Sex <- "Female"
df_male$Sex   <- "Male"

df_all <- bind_rows(df_female, df_male)

write.csv(
  df_female,
  file.path(cellchat_outdir, "Female_CellChat_All_Communications.csv"),
  row.names = FALSE
)

write.csv(
  df_male,
  file.path(cellchat_outdir, "Male_CellChat_All_Communications.csv"),
  row.names = FALSE
)

write.csv(
  df_all,
  file.path(cellchat_outdir, "Female_Male_CellChat_All_Communications.csv"),
  row.names = FALSE
)


# ============================================================
# 12. KEEP ONLY FIBROBLAST <-> SCHWANN COMMUNICATIONS
# ============================================================

fib_schwann_all <- df_all %>%
  filter(
    (source == "Fibroblasts" & target == "Schwann cells") |
      (source == "Schwann cells" & target == "Fibroblasts")
  ) %>%
  mutate(
    Direction = paste(source, "→", target),
    LR_pair = ifelse(
      "ligand" %in% colnames(.) & "receptor" %in% colnames(.),
      paste(ligand, receptor, sep = " → "),
      interaction_name
    )
  )

if (nrow(fib_schwann_all) == 0) {
  stop("No Fibroblast <-> Schwann interactions were detected. Check source/target names or CellChat results.")
}

write.csv(
  fib_schwann_all,
  file.path(cellchat_outdir, "Fibroblast_Schwann_CrossTalk_Female_vs_Male.csv"),
  row.names = FALSE
)

print("Number of Fibroblast <-> Schwann interactions by sex and direction:")
print(
  fib_schwann_all %>%
    group_by(Sex, Direction) %>%
    summarise(
      Interaction_Count = n(),
      Total_Strength = sum(prob, na.rm = TRUE),
      Mean_Strength = mean(prob, na.rm = TRUE),
      .groups = "drop"
    )
)


# ============================================================
# 13. CROSS-TALK SUMMARY TABLE
# ============================================================

cross_talk_summary <- fib_schwann_all %>%
  group_by(Sex, Direction) %>%
  summarise(
    Interaction_Count = n(),
    Total_Interaction_Strength = sum(prob, na.rm = TRUE),
    Mean_Interaction_Strength = mean(prob, na.rm = TRUE),
    .groups = "drop"
  )

print(cross_talk_summary)

write.csv(
  cross_talk_summary,
  file.path(cellchat_outdir, "Fibroblast_Schwann_CrossTalk_Summary.csv"),
  row.names = FALSE
)


# ============================================================
# 14. FIGURE A: BUILT-IN TOTAL INTERACTION COUNT
# ============================================================
# This summarizes all interactions inside the two-cell-type network:
# Fibroblast -> Fibroblast
# Fibroblast -> Schwann
# Schwann -> Fibroblast
# Schwann -> Schwann

p_count <- compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "count"
) +
  ggtitle("Number of Fibroblast/Schwann Network Interactions") +
  labs(y = "Number of inferred interactions") +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black")
  )

print(p_count)

ggsave(
  file.path(cellchat_outdir, "Figure_A_Total_Interaction_Count_Female_vs_Male.pdf"),
  plot = p_count,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(cellchat_outdir, "Figure_A_Total_Interaction_Count_Female_vs_Male.png"),
  plot = p_count,
  width = 7,
  height = 5,
  dpi = 300
)


# ============================================================
# 15. FIGURE B: BUILT-IN TOTAL INTERACTION STRENGTH
# ============================================================

p_strength <- compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "weight"
) +
  ggtitle("Total Fibroblast/Schwann Network Interaction Strength") +
  labs(y = "Interaction strength") +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black")
  )

print(p_strength)

ggsave(
  file.path(cellchat_outdir, "Figure_B_Total_Interaction_Strength_Female_vs_Male.pdf"),
  plot = p_strength,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(cellchat_outdir, "Figure_B_Total_Interaction_Strength_Female_vs_Male.png"),
  plot = p_strength,
  width = 7,
  height = 5,
  dpi = 300
)


# ============================================================
# 16. FIGURE C: CROSS-TALK ONLY COUNT BARPLOT
# ============================================================

cross_count_plot <- ggplot(
  cross_talk_summary,
  aes(x = Direction, y = Interaction_Count, fill = Sex)
) +
  geom_col(
    position = position_dodge(width = 0.75),
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_text(
    aes(label = Interaction_Count),
    position = position_dodge(width = 0.75),
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_manual(values = c(
    "Female" = "#F8766D",
    "Male" = "#00BFC4"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Fibroblast ↔ Schwann Cross-Talk Interaction Count",
    x = "Communication Direction",
    y = "Number of inferred interactions",
    fill = "Patient Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(size = 11, face = "bold", color = "black", angle = 20, hjust = 1),
    axis.text.y = element_text(size = 11, color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11)
  )

print(cross_count_plot)

ggsave(
  file.path(cellchat_outdir, "Figure_C_CrossTalk_Count_By_Direction.pdf"),
  plot = cross_count_plot,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(cellchat_outdir, "Figure_C_CrossTalk_Count_By_Direction.png"),
  plot = cross_count_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 17. FIGURE D: CROSS-TALK ONLY STRENGTH BARPLOT
# ============================================================

cross_strength_plot <- ggplot(
  cross_talk_summary,
  aes(x = Direction, y = Total_Interaction_Strength, fill = Sex)
) +
  geom_col(
    position = position_dodge(width = 0.75),
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_text(
    aes(label = round(Total_Interaction_Strength, 3)),
    position = position_dodge(width = 0.75),
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_manual(values = c(
    "Female" = "#F8766D",
    "Male" = "#00BFC4"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Fibroblast ↔ Schwann Cross-Talk Interaction Strength",
    x = "Communication Direction",
    y = "Total interaction strength",
    fill = "Patient Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(size = 11, face = "bold", color = "black", angle = 20, hjust = 1),
    axis.text.y = element_text(size = 11, color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11)
  )

print(cross_strength_plot)

ggsave(
  file.path(cellchat_outdir, "Figure_D_CrossTalk_Strength_By_Direction.pdf"),
  plot = cross_strength_plot,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(cellchat_outdir, "Figure_D_CrossTalk_Strength_By_Direction.png"),
  plot = cross_strength_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 18. BUILT-IN DIFFERENTIAL NETWORK PLOTS
# ============================================================
# These may work even if netVisual_bubble gave you an error.
# Wrapped in tryCatch so the script continues if CellChat plotting fails.

tryCatch({
  
  pdf(
    file.path(cellchat_outdir, "Figure_E_Differential_Network_Weight_Male_minus_Female.pdf"),
    width = 7,
    height = 7
  )
  
  netVisual_diffInteraction(
    cellchat_merged,
    comparison = c(1, 2),
    weight.scale = TRUE,
    measure = "weight",
    title.name = "Differential Interaction Strength: Male - Female"
  )
  
  dev.off()
  
}, error = function(e) {
  message("Differential weight network plot failed:")
  message(e$message)
  try(dev.off(), silent = TRUE)
})


tryCatch({
  
  pdf(
    file.path(cellchat_outdir, "Figure_F_Differential_Network_Count_Male_minus_Female.pdf"),
    width = 7,
    height = 7
  )
  
  netVisual_diffInteraction(
    cellchat_merged,
    comparison = c(1, 2),
    weight.scale = TRUE,
    measure = "count",
    title.name = "Differential Interaction Count: Male - Female"
  )
  
  dev.off()
  
}, error = function(e) {
  message("Differential count network plot failed:")
  message(e$message)
  try(dev.off(), silent = TRUE)
})


# ============================================================
# 19. CUSTOM BUBBLE PLOTS
# This replaces netVisual_bubble() to avoid the:
# wrong sign in 'by' argument error.
# ============================================================

top_lr_pairs <- fib_schwann_all %>%
  group_by(LR_pair) %>%
  summarise(
    Max_Prob = max(prob, na.rm = TRUE),
    Mean_Prob = mean(prob, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Max_Prob)) %>%
  slice_head(n = 30) %>%
  pull(LR_pair)

fib_schwann_top <- fib_schwann_all %>%
  filter(LR_pair %in% top_lr_pairs)

fib_schwann_top$LR_pair <- factor(
  fib_schwann_top$LR_pair,
  levels = rev(top_lr_pairs)
)

fib_schwann_top$Sex <- factor(
  fib_schwann_top$Sex,
  levels = c("Female", "Male")
)


# ============================================================
# 20. FIGURE G: CUSTOM BUBBLE PLOT, BOTH DIRECTIONS
# ============================================================

custom_bubble_both <- ggplot(
  fib_schwann_top,
  aes(x = Sex, y = LR_pair)
) +
  geom_point(
    aes(size = prob, color = -log10(pval)),
    alpha = 0.85
  ) +
  facet_wrap(~Direction, scales = "free_y") +
  scale_color_gradient(
    low = "lightgrey",
    high = "darkred",
    name = "-log10(p-value)"
  ) +
  scale_size_continuous(
    name = "Communication\nprobability",
    range = c(2, 8)
  ) +
  labs(
    title = "Fibroblast ↔ Schwann Cell Signaling: Female vs Male",
    x = "Patient Sex",
    y = "Ligand → Receptor Pair"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16,
      margin = margin(b = 15)
    ),
    axis.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold", size = 11),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

print(custom_bubble_both)

ggsave(
  file.path(cellchat_outdir, "Figure_G_Custom_Bubble_Fibroblast_Schwann_Both_Directions.pdf"),
  plot = custom_bubble_both,
  width = 11,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(cellchat_outdir, "Figure_G_Custom_Bubble_Fibroblast_Schwann_Both_Directions.png"),
  plot = custom_bubble_both,
  width = 11,
  height = 8,
  dpi = 300
)


# ============================================================
# 21. FIGURE H: CUSTOM BUBBLE, FIBROBLAST -> SCHWANN ONLY
# ============================================================

fib_to_schwann_top <- fib_schwann_top %>%
  filter(source == "Fibroblasts", target == "Schwann cells")

if (nrow(fib_to_schwann_top) > 0) {
  
  custom_bubble_fib_to_schwann <- ggplot(
    fib_to_schwann_top,
    aes(x = Sex, y = LR_pair)
  ) +
    geom_point(
      aes(size = prob, color = -log10(pval)),
      alpha = 0.85
    ) +
    scale_color_gradient(
      low = "lightgrey",
      high = "darkred",
      name = "-log10(p-value)"
    ) +
    scale_size_continuous(
      name = "Communication\nprobability",
      range = c(2, 8)
    ) +
    labs(
      title = "Fibroblast → Schwann Cell Signaling: Female vs Male",
      x = "Patient Sex",
      y = "Ligand → Receptor Pair"
    ) +
    theme_classic(base_family = "sans") +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16,
        margin = margin(b = 15)
      ),
      axis.title = element_text(face = "bold", size = 13),
      axis.text.x = element_text(size = 12, face = "bold", color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      legend.title = element_text(face = "bold")
    )
  
  print(custom_bubble_fib_to_schwann)
  
  ggsave(
    file.path(cellchat_outdir, "Figure_H_Custom_Bubble_Fibroblast_to_Schwann.pdf"),
    plot = custom_bubble_fib_to_schwann,
    width = 9,
    height = 8,
    dpi = 300
  )
  
  ggsave(
    file.path(cellchat_outdir, "Figure_H_Custom_Bubble_Fibroblast_to_Schwann.png"),
    plot = custom_bubble_fib_to_schwann,
    width = 9,
    height = 8,
    dpi = 300
  )
}


# ============================================================
# 22. FIGURE I: CUSTOM BUBBLE, SCHWANN -> FIBROBLAST ONLY
# ============================================================

schwann_to_fib_top <- fib_schwann_top %>%
  filter(source == "Schwann cells", target == "Fibroblasts")

if (nrow(schwann_to_fib_top) > 0) {
  
  custom_bubble_schwann_to_fib <- ggplot(
    schwann_to_fib_top,
    aes(x = Sex, y = LR_pair)
  ) +
    geom_point(
      aes(size = prob, color = -log10(pval)),
      alpha = 0.85
    ) +
    scale_color_gradient(
      low = "lightgrey",
      high = "darkred",
      name = "-log10(p-value)"
    ) +
    scale_size_continuous(
      name = "Communication\nprobability",
      range = c(2, 8)
    ) +
    labs(
      title = "Schwann Cell → Fibroblast Signaling: Female vs Male",
      x = "Patient Sex",
      y = "Ligand → Receptor Pair"
    ) +
    theme_classic(base_family = "sans") +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16,
        margin = margin(b = 15)
      ),
      axis.title = element_text(face = "bold", size = 13),
      axis.text.x = element_text(size = 12, face = "bold", color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      legend.title = element_text(face = "bold")
    )
  
  print(custom_bubble_schwann_to_fib)
  
  ggsave(
    file.path(cellchat_outdir, "Figure_I_Custom_Bubble_Schwann_to_Fibroblast.pdf"),
    plot = custom_bubble_schwann_to_fib,
    width = 9,
    height = 8,
    dpi = 300
  )
  
  ggsave(
    file.path(cellchat_outdir, "Figure_I_Custom_Bubble_Schwann_to_Fibroblast.png"),
    plot = custom_bubble_schwann_to_fib,
    width = 9,
    height = 8,
    dpi = 300
  )
}


# ============================================================
# 23. TOP FEMALE AND MALE CROSS-TALK INTERACTIONS
# ============================================================

top_female_interactions <- fib_schwann_all %>%
  filter(Sex == "Female") %>%
  arrange(desc(prob)) %>%
  dplyr::select(
    Sex, source, target, pathway_name, ligand, receptor,
    interaction_name, prob, pval
  ) %>%
  head(25)

top_male_interactions <- fib_schwann_all %>%
  filter(Sex == "Male") %>%
  arrange(desc(prob)) %>%
  dplyr::select(
    Sex, source, target, pathway_name, ligand, receptor,
    interaction_name, prob, pval
  ) %>%
  head(25)

print("Top Female Fibroblast <-> Schwann interactions:")
print(top_female_interactions)

print("Top Male Fibroblast <-> Schwann interactions:")
print(top_male_interactions)

write.csv(
  top_female_interactions,
  file.path(cellchat_outdir, "Top_25_Female_Fibroblast_Schwann_Interactions.csv"),
  row.names = FALSE
)

write.csv(
  top_male_interactions,
  file.path(cellchat_outdir, "Top_25_Male_Fibroblast_Schwann_Interactions.csv"),
  row.names = FALSE
)


# ============================================================
# 24. FINAL SUMMARY
# ============================================================

cat("\n======================================================\n")
cat("CELLCHAT FIBROBLAST-SCHWANN ANALYSIS COMPLETE\n")
cat("Results saved to:\n")
cat(cellchat_outdir, "\n")
cat("======================================================\n")

cat("\nIMPORTANT INTERPRETATION NOTES:\n")
cat("1. compareInteractions() summarizes the full two-cell-type network.\n")
cat("2. The custom cross-talk plots summarize only Fibroblast -> Schwann and Schwann -> Fibroblast interactions.\n")
cat("3. Communication probability is inferred from ligand-receptor expression, not direct physical contact.\n")
cat("4. Female is dataset 1 and Male is dataset 2 in the merged CellChat object.\n")
cat("======================================================\n")

# ============================================================
# CIRCULAR CELLCHAT PLOTS:
# Fibroblast <-> Schwann Cell Interactions
# ============================================================

# Make a subfolder for circle plots
circle_outdir <- file.path(cellchat_outdir, "Circle_Plots")
dir.create(circle_outdir, showWarnings = FALSE, recursive = TRUE)


# ------------------------------------------------------------
# Helper function: plot only Fibroblast <-> Schwann interactions
# Removes self-interactions:
# Fibroblast -> Fibroblast
# Schwann -> Schwann
# ------------------------------------------------------------

plot_fib_schwann_circle <- function(cellchat_obj, dataset_name, measure = "weight") {
  
  # Choose count or weight matrix
  if (measure == "weight") {
    mat <- cellchat_obj@net$weight
    plot_label <- "Interaction Strength"
  } else if (measure == "count") {
    mat <- cellchat_obj@net$count
    plot_label <- "Interaction Count"
  } else {
    stop("measure must be either 'weight' or 'count'")
  }
  
  # Keep only Fibroblasts and Schwann cells
  keep_groups <- c("Fibroblasts", "Schwann cells")
  keep_groups <- keep_groups[keep_groups %in% rownames(mat)]
  
  if (length(keep_groups) < 2) {
    stop(paste("Could not find both Fibroblasts and Schwann cells in", dataset_name))
  }
  
  mat_use <- mat[keep_groups, keep_groups, drop = FALSE]
  
  # Remove self-loops so the circle only shows Fibroblast <-> Schwann crosstalk
  diag(mat_use) <- 0
  
  # Cell counts for node sizes
  group_size <- as.numeric(table(cellchat_obj@idents)[keep_groups])
  names(group_size) <- keep_groups
  
  # Plot title
  title_use <- paste0(
    dataset_name,
    ": Fibroblast \u2194 Schwann Cell ",
    plot_label
  )
  
  # Plot
  netVisual_circle(
    mat_use,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = TRUE,
    edge.weight.max = max(mat_use),
    edge.width.max = 8,
    title.name = title_use
  )
}

###
# ============================================================
# SAVE CIRCLE PLOTS: INTERACTION STRENGTH
# ============================================================

pdf(
  file.path(circle_outdir, "Female_Fibroblast_Schwann_Circle_Strength.pdf"),
  width = 7,
  height = 7
)
plot_fib_schwann_circle(
  cellchat_obj = cellchat_female,
  dataset_name = "Female",
  measure = "weight"
)
dev.off()


pdf(
  file.path(circle_outdir, "Male_Fibroblast_Schwann_Circle_Strength.pdf"),
  width = 7,
  height = 7
)
plot_fib_schwann_circle(
  cellchat_obj = cellchat_male,
  dataset_name = "Male",
  measure = "weight"
)
dev.off()

# ============================================================
# SAVE CIRCLE PLOTS: INTERACTION COUNT
# ============================================================

pdf(
  file.path(circle_outdir, "Female_Fibroblast_Schwann_Circle_Count.pdf"),
  width = 7,
  height = 7
)
plot_fib_schwann_circle(
  cellchat_obj = cellchat_female,
  dataset_name = "Female",
  measure = "count"
)
dev.off()


pdf(
  file.path(circle_outdir, "Male_Fibroblast_Schwann_Circle_Count.pdf"),
  width = 7,
  height = 7
)
plot_fib_schwann_circle(
  cellchat_obj = cellchat_male,
  dataset_name = "Male",
  measure = "count"
)
dev.off()


# ============================================================
#SAVE PNG VERSIONS
# ============================================================

png(
  file.path(circle_outdir, "Female_Fibroblast_Schwann_Circle_Strength.png"),
  width = 1800,
  height = 1800,
  res = 300
)
plot_fib_schwann_circle(cellchat_female, "Female", measure = "weight")
dev.off()


png(
  file.path(circle_outdir, "Male_Fibroblast_Schwann_Circle_Strength.png"),
  width = 1800,
  height = 1800,
  res = 300
)
plot_fib_schwann_circle(cellchat_male, "Male", measure = "weight")
dev.off()

# Checking pathways

# ============================================================
# CHECKING PATHWAYS - PRINT ONLY, DO NOT SAVE
# ============================================================

pathway_summary <- fib_schwann_all %>%
  mutate(
    LR_pair = paste(ligand, receptor, sep = " → "),
    pathway_name = ifelse(is.na(pathway_name), "Unknown", pathway_name)
  ) %>%
  group_by(Sex, Direction, pathway_name) %>%
  summarise(
    Interaction_Count = n(),
    Total_Prob = sum(prob, na.rm = TRUE),
    Mean_Prob = mean(prob, na.rm = TRUE),
    .groups = "drop"
  )

# Check that pathway_summary has rows
print("Pathway summary:")
print(pathway_summary)

# Select top pathways
top_pathways <- pathway_summary %>%
  group_by(pathway_name) %>%
  summarise(
    Max_Total_Prob = max(Total_Prob, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Max_Total_Prob)) %>%
  slice_head(n = 12) %>%
  pull(pathway_name)

print("Top pathways:")
print(top_pathways)

# Filter data for plotting
pathway_plot_data <- pathway_summary %>%
  filter(pathway_name %in% top_pathways)

print("Number of rows being plotted:")
print(nrow(pathway_plot_data))

# Create plot object
pathway_plot <- ggplot(
  pathway_plot_data,
  aes(
    x = reorder(pathway_name, Total_Prob),
    y = Total_Prob,
    fill = Sex
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    color = "black",
    linewidth = 0.25
  ) +
  coord_flip() +
  facet_wrap(~Direction, scales = "free_y") +
  labs(
    title = "Fibroblast ↔ Schwann Signaling Pathways",
    x = "Pathway",
    y = "Total communication probability",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11)
  )

# ============================================================
# FORCE R TO PRINT TO RSTUDIO PLOTS PANE
# ============================================================

# Close any open PDF/PNG plotting devices
graphics.off()

# Test whether plotting works at all
test_plot <- ggplot(mtcars, aes(x = factor(cyl), y = mpg)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("Test Plot")

print(test_plot)

# IMPORTANT: print the plot so it appears in RStudio Plots pane
print(pathway_plot)


### Schwann - Fibroblast - iCAF only

# ============================================================
# CELLCHAT ANALYSIS:
# iCAF <-> Schwann Cell Interaction Pathways
# Female vs Male PDAC snRNA-seq Samples
# ============================================================

# ============================================================
# 0. LOAD REQUIRED LIBRARIES
# ============================================================

library(CellChat)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

options(stringsAsFactors = FALSE)

# ============================================================
# 1. SET PATHS
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

merged_path <- file.path(base_dir, "Merged_Joined.RDS")
fibroblast_path <- file.path(base_dir, "Fibroblast_Subset.RDS")

icaf_path_1 <- file.path(base_dir, "iCAF_Subset.RDS")
icaf_path_2 <- file.path(base_dir, "iCAF_subset.RDS")

cellchat_outdir <- file.path(base_dir, "CellChat_iCAF_Schwann_Pathways")

dir.create(
  cellchat_outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# 2. READ MAIN SEURAT OBJECT
# ============================================================

Merged_Joined <- readRDS(merged_path)

DefaultAssay(Merged_Joined) <- "RNA"

Merged_Joined <- tryCatch(
  JoinLayers(Merged_Joined),
  error = function(e) Merged_Joined
)

# ============================================================
# 3. READ OR CREATE iCAF SUBSET
# ============================================================
# Cluster mapping in fibroblast_subset:
# 0 = Fibroblasts
# 1 = myCAF
# 2 = iCAF
# 3 = pancreatic stellate
# 4 = apCAF

if (file.exists(icaf_path_1)) {
  
  iCAF_subset <- readRDS(icaf_path_1)
  print(paste("Loaded iCAF subset from:", icaf_path_1))
  
} else if (file.exists(icaf_path_2)) {
  
  iCAF_subset <- readRDS(icaf_path_2)
  print(paste("Loaded iCAF subset from:", icaf_path_2))
  
} else {
  
  print("iCAF subset RDS not found. Creating iCAF subset from Fibroblast_Subset.RDS...")
  
  fibroblast_subset <- readRDS(fibroblast_path)
  DefaultAssay(fibroblast_subset) <- "RNA"
  
  fibroblast_subset <- tryCatch(
    JoinLayers(fibroblast_subset),
    error = function(e) fibroblast_subset
  )
  
  iCAF_subset <- subset(
    fibroblast_subset,
    subset = seurat_clusters == 2
  )
  
  saveRDS(
    iCAF_subset,
    icaf_path_1
  )
  
  print(paste("Saved new iCAF subset to:", icaf_path_1))
}

DefaultAssay(iCAF_subset) <- "RNA"

iCAF_subset <- tryCatch(
  JoinLayers(iCAF_subset),
  error = function(e) iCAF_subset
)

print("Number of iCAF cells:")
print(ncol(iCAF_subset))

# ============================================================
# 4. ADD SEX METADATA
# ============================================================

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

Merged_Joined$Sex <- unname(
  sex_mapping[as.character(Merged_Joined$orig.ident)]
)

Merged_Joined$Sex <- factor(
  Merged_Joined$Sex,
  levels = c("Female", "Male")
)

# ============================================================
# 5. CREATE CLEAN CELLCHAT CELL TYPE LABELS
# ============================================================

Merged_Joined$cellchat_celltype <- NA_character_

# Label iCAFs using saved iCAF subset cells
icaf_cells <- intersect(
  colnames(iCAF_subset),
  colnames(Merged_Joined)
)

Merged_Joined$cellchat_celltype[
  colnames(Merged_Joined) %in% icaf_cells
] <- "iCAFs"

# Label true Schwann cells
if ("true_schwann_subset" %in% colnames(Merged_Joined@meta.data)) {
  
  Merged_Joined$cellchat_celltype[
    as.character(Merged_Joined$true_schwann_subset) %in% c(
      "True Schwann",
      "TRUE",
      "True",
      "true"
    )
  ] <- "Schwann cells"
  
} else {
  
  warning("true_schwann_subset column not found. Using seurat cluster 11 as Schwann fallback.")
  
  Merged_Joined$cellchat_celltype[
    as.character(Merged_Joined$seurat_clusters) == "11"
  ] <- "Schwann cells"
}

# Keep only iCAFs and Schwann cells with clear Sex labels
target_groups <- c(
  "iCAFs",
  "Schwann cells"
)

target_cells <- subset(
  Merged_Joined,
  subset = cellchat_celltype %in% target_groups &
    Sex %in% c("Female", "Male")
)

target_cells$cellchat_celltype <- factor(
  target_cells$cellchat_celltype,
  levels = target_groups
)

target_cells$Sex <- factor(
  target_cells$Sex,
  levels = c("Female", "Male")
)

print("Cell counts by CellChat group and sex:")
print(table(target_cells$cellchat_celltype, target_cells$Sex))

print("Cell counts by patient, CellChat group, and sex:")
print(table(target_cells$orig.ident, target_cells$cellchat_celltype, target_cells$Sex))

# ============================================================
# 6. SPLIT INTO FEMALE AND MALE OBJECTS
# ============================================================

seurat_female <- subset(
  target_cells,
  subset = Sex == "Female"
)

seurat_male <- subset(
  target_cells,
  subset = Sex == "Male"
)

print("Female object cell counts:")
print(table(seurat_female$cellchat_celltype))

print("Male object cell counts:")
print(table(seurat_male$cellchat_celltype))

# ============================================================
# 7. SAFETY CHECK
# ============================================================

min_cells_required <- 10

female_counts <- table(seurat_female$cellchat_celltype)
male_counts <- table(seurat_male$cellchat_celltype)

if (any(female_counts < min_cells_required)) {
  warning("Female subset has a cell type with fewer than 10 cells. Interpret CellChat cautiously.")
}

if (any(male_counts < min_cells_required)) {
  warning("Male subset has a cell type with fewer than 10 cells. Interpret CellChat cautiously.")
}

# ============================================================
# 8. HELPER FUNCTION TO GET RNA DATA
# ============================================================

get_rna_data <- function(seurat_obj) {
  
  data_input <- tryCatch(
    GetAssayData(
      seurat_obj,
      assay = "RNA",
      layer = "data"
    ),
    error = function(e) {
      GetAssayData(
        seurat_obj,
        assay = "RNA",
        slot = "data"
      )
    }
  )
  
  return(data_input)
}

# ============================================================
# 9. CELLCHAT FUNCTION
# ============================================================

run_cellchat <- function(seurat_obj, dataset_name, min_cells_filter = 10) {
  
  print(paste("Running CellChat for", dataset_name, "..."))
  
  DefaultAssay(seurat_obj) <- "RNA"
  
  data.input <- get_rna_data(seurat_obj)
  
  meta <- seurat_obj@meta.data
  
  meta$cellchat_celltype <- droplevels(
    factor(
      meta$cellchat_celltype,
      levels = target_groups
    )
  )
  
  # Keeps patient/sample information
  if (!"samples" %in% colnames(meta)) {
    meta$samples <- as.character(meta$orig.ident)
  }
  
  print(paste("CellChat groups for", dataset_name))
  print(table(meta$cellchat_celltype))
  
  print(paste("Samples represented in", dataset_name))
  print(table(meta$samples))
  
  cellchat <- createCellChat(
    object = data.input,
    meta = meta,
    group.by = "cellchat_celltype"
  )
  
  # Human ligand-receptor database
  CellChatDB <- CellChatDB.human
  cellchat@DB <- CellChatDB
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(
    cellchat,
    type = "triMean",
    raw.use = TRUE
  )
  
  cellchat <- filterCommunication(
    cellchat,
    min.cells = min_cells_filter
  )
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  return(cellchat)
}

# ============================================================
# 10. RUN CELLCHAT FOR FEMALE AND MALE
# ============================================================

cellchat_female <- run_cellchat(
  seurat_obj = seurat_female,
  dataset_name = "Female",
  min_cells_filter = 10
)

cellchat_male <- run_cellchat(
  seurat_obj = seurat_male,
  dataset_name = "Male",
  min_cells_filter = 10
)

cellchat_female <- updateCellChat(cellchat_female)
cellchat_male <- updateCellChat(cellchat_male)

# ============================================================
# 11. SAVE INDIVIDUAL CELLCHAT OBJECTS
# ============================================================

saveRDS(
  cellchat_female,
  file.path(cellchat_outdir, "cellchat_female_iCAF_schwann.rds")
)

saveRDS(
  cellchat_male,
  file.path(cellchat_outdir, "cellchat_male_iCAF_schwann.rds")
)

# ============================================================
# 12. MERGE CELLCHAT OBJECTS
# ============================================================

object.list <- list(
  Female = cellchat_female,
  Male = cellchat_male
)

cellchat_merged <- mergeCellChat(
  object.list,
  add.names = names(object.list)
)

saveRDS(
  cellchat_merged,
  file.path(cellchat_outdir, "cellchat_merged_female_male_iCAF_schwann.rds")
)

# ============================================================
# 13. EXTRACT COMMUNICATION TABLES
# ============================================================

df_female <- subsetCommunication(cellchat_female)
df_male <- subsetCommunication(cellchat_male)

df_female$Sex <- "Female"
df_male$Sex <- "Male"

df_all <- bind_rows(df_female, df_male)

write.csv(
  df_female,
  file.path(cellchat_outdir, "Female_CellChat_All_Communications_iCAF_Schwann.csv"),
  row.names = FALSE
)

write.csv(
  df_male,
  file.path(cellchat_outdir, "Male_CellChat_All_Communications_iCAF_Schwann.csv"),
  row.names = FALSE
)

write.csv(
  df_all,
  file.path(cellchat_outdir, "Female_Male_CellChat_All_Communications_iCAF_Schwann.csv"),
  row.names = FALSE
)

# ============================================================
# 14. KEEP ONLY iCAF <-> SCHWANN COMMUNICATIONS
# ============================================================

icaf_schwann_all <- df_all %>%
  filter(
    (source == "iCAFs" & target == "Schwann cells") |
      (source == "Schwann cells" & target == "iCAFs")
  ) %>%
  mutate(
    Direction = paste(source, "→", target),
    LR_pair = ifelse(
      "ligand" %in% colnames(.) & "receptor" %in% colnames(.),
      paste(ligand, receptor, sep = " → "),
      interaction_name
    ),
    pathway_name = ifelse(
      is.na(pathway_name),
      "Unknown",
      pathway_name
    ),
    pval_safe = ifelse(
      pval == 0,
      .Machine$double.xmin,
      pval
    ),
    neg_log10_pval = -log10(pval_safe)
  )

if (nrow(icaf_schwann_all) == 0) {
  stop("No iCAF <-> Schwann interactions were detected. Check source/target labels or CellChat results.")
}

write.csv(
  icaf_schwann_all,
  file.path(cellchat_outdir, "iCAF_Schwann_CrossTalk_Female_vs_Male.csv"),
  row.names = FALSE
)

print("Number of iCAF <-> Schwann interactions by sex and direction:")
print(
  icaf_schwann_all %>%
    group_by(Sex, Direction) %>%
    summarise(
      Interaction_Count = n(),
      Total_Strength = sum(prob, na.rm = TRUE),
      Mean_Strength = mean(prob, na.rm = TRUE),
      .groups = "drop"
    )
)

# ============================================================
# 15. PATHWAY SUMMARY TABLE
# ============================================================

pathway_summary <- icaf_schwann_all %>%
  group_by(Sex, Direction, pathway_name) %>%
  summarise(
    Interaction_Count = n(),
    Total_Prob = sum(prob, na.rm = TRUE),
    Mean_Prob = mean(prob, na.rm = TRUE),
    .groups = "drop"
  )

print("iCAF <-> Schwann pathway summary:")
print(pathway_summary)

write.csv(
  pathway_summary,
  file.path(cellchat_outdir, "iCAF_Schwann_Pathway_Summary_Female_vs_Male.csv"),
  row.names = FALSE
)

# ============================================================
# 16. SELECT TOP PATHWAYS FOR PLOTTING
# ============================================================

top_pathways <- pathway_summary %>%
  group_by(pathway_name) %>%
  summarise(
    Max_Total_Prob = max(Total_Prob, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Max_Total_Prob)) %>%
  slice_head(n = 12) %>%
  pull(pathway_name)

print("Top iCAF <-> Schwann pathways:")
print(top_pathways)

pathway_plot_data <- pathway_summary %>%
  filter(pathway_name %in% top_pathways)

print("Number of pathway rows being plotted:")
print(nrow(pathway_plot_data))

# ============================================================
# 17. PLOT iCAF <-> SCHWANN SIGNALING PATHWAYS
# ============================================================

pathway_plot <- ggplot(
  pathway_plot_data,
  aes(
    x = reorder(pathway_name, Total_Prob),
    y = Total_Prob,
    fill = Sex
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    color = "black",
    linewidth = 0.25,
    width = 0.7
  ) +
  coord_flip() +
  facet_wrap(
    ~Direction,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  labs(
    title = "iCAF ↔ Schwann Cell Signaling Pathways",
    x = "Pathway",
    y = "Total communication probability",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    strip.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.6
    ),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text = element_text(
      size = 11
    )
  )

print(pathway_plot)

# Force plot to show in RStudio Plots pane
graphics.off()
print(pathway_plot)

ggsave(
  file.path(cellchat_outdir, "iCAF_Schwann_Signaling_Pathways_Female_vs_Male.png"),
  plot = pathway_plot,
  width = 11,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(cellchat_outdir, "iCAF_Schwann_Signaling_Pathways_Female_vs_Male.pdf"),
  plot = pathway_plot,
  width = 11,
  height = 7,
  dpi = 300
)

# ============================================================
# 18. OPTIONAL: PLOT iCAF -> SCHWANN ONLY
# ============================================================

icaf_to_schwann_plot_data <- pathway_plot_data %>%
  filter(Direction == "iCAFs → Schwann cells")

if (nrow(icaf_to_schwann_plot_data) > 0) {
  
  icaf_to_schwann_plot <- ggplot(
    icaf_to_schwann_plot_data,
    aes(
      x = reorder(pathway_name, Total_Prob),
      y = Total_Prob,
      fill = Sex
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.75),
      color = "black",
      linewidth = 0.25,
      width = 0.7
    ) +
    coord_flip() +
    scale_fill_manual(
      values = c(
        "Female" = "#F8766D",
        "Male" = "#00BFC4"
      )
    ) +
    labs(
      title = "iCAF → Schwann Cell Signaling Pathways",
      x = "Pathway",
      y = "Total communication probability",
      fill = "Sex"
    ) +
    theme_classic(base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12, color = "black"),
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11)
    )
  
  print(icaf_to_schwann_plot)
  
  ggsave(
    file.path(cellchat_outdir, "iCAF_to_Schwann_Signaling_Pathways_Female_vs_Male.png"),
    plot = icaf_to_schwann_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  ggsave(
    file.path(cellchat_outdir, "iCAF_to_Schwann_Signaling_Pathways_Female_vs_Male.pdf"),
    plot = icaf_to_schwann_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# ============================================================
# 19. OPTIONAL: PLOT SCHWANN -> iCAF ONLY
# ============================================================

schwann_to_icaf_plot_data <- pathway_plot_data %>%
  filter(Direction == "Schwann cells → iCAFs")

if (nrow(schwann_to_icaf_plot_data) > 0) {
  
  schwann_to_icaf_plot <- ggplot(
    schwann_to_icaf_plot_data,
    aes(
      x = reorder(pathway_name, Total_Prob),
      y = Total_Prob,
      fill = Sex
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.75),
      color = "black",
      linewidth = 0.25,
      width = 0.7
    ) +
    coord_flip() +
    scale_fill_manual(
      values = c(
        "Female" = "#F8766D",
        "Male" = "#00BFC4"
      )
    ) +
    labs(
      title = "Schwann Cell → iCAF Signaling Pathways",
      x = "Pathway",
      y = "Total communication probability",
      fill = "Sex"
    ) +
    theme_classic(base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12, color = "black"),
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11)
    )
  
  print(schwann_to_icaf_plot)
  
  ggsave(
    file.path(cellchat_outdir, "Schwann_to_iCAF_Signaling_Pathways_Female_vs_Male.png"),
    plot = schwann_to_icaf_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  ggsave(
    file.path(cellchat_outdir, "Schwann_to_iCAF_Signaling_Pathways_Female_vs_Male.pdf"),
    plot = schwann_to_icaf_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# ============================================================
# 20. TOP FEMALE AND MALE iCAF <-> SCHWANN INTERACTIONS
# ============================================================

top_female_interactions <- icaf_schwann_all %>%
  filter(Sex == "Female") %>%
  arrange(desc(prob)) %>%
  dplyr::select(
    Sex,
    source,
    target,
    Direction,
    pathway_name,
    ligand,
    receptor,
    interaction_name,
    prob,
    pval
  ) %>%
  head(25)

top_male_interactions <- icaf_schwann_all %>%
  filter(Sex == "Male") %>%
  arrange(desc(prob)) %>%
  dplyr::select(
    Sex,
    source,
    target,
    Direction,
    pathway_name,
    ligand,
    receptor,
    interaction_name,
    prob,
    pval
  ) %>%
  head(25)

print("Top Female iCAF <-> Schwann interactions:")
print(top_female_interactions)

print("Top Male iCAF <-> Schwann interactions:")
print(top_male_interactions)

write.csv(
  top_female_interactions,
  file.path(cellchat_outdir, "Top_25_Female_iCAF_Schwann_Interactions.csv"),
  row.names = FALSE
)

write.csv(
  top_male_interactions,
  file.path(cellchat_outdir, "Top_25_Male_iCAF_Schwann_Interactions.csv"),
  row.names = FALSE
)

# ============================================================
# 21. FINAL SUMMARY
# ============================================================

cat("\n======================================================\n")
cat("CELLCHAT iCAF-SCHWANN PATHWAY ANALYSIS COMPLETE\n")
cat("Results saved to:\n")
cat(cellchat_outdir, "\n")
cat("Main pathway plot saved as:\n")
cat("iCAF_Schwann_Signaling_Pathways_Female_vs_Male.png/pdf\n")
cat("======================================================\n")

cat("\nIMPORTANT INTERPRETATION NOTES:\n")
cat("1. These are CellChat-defined signaling pathways, not GSEA pathways.\n")
cat("2. Pathway strength is based on summed CellChat communication probability.\n")
cat("3. Direction matters: iCAF -> Schwann is different from Schwann -> iCAF.\n")
cat("4. Communication probability is inferred from ligand-receptor expression, not direct physical contact.\n")
cat("======================================================\n")



### Schwann - Fibroblast - myCAF only

# ============================================================
# CELLCHAT ANALYSIS:
# myCAF <-> Schwann Cell Interaction Pathways
# Female vs Male PDAC snRNA-seq Samples
# ============================================================

# ============================================================
# 0. LOAD REQUIRED LIBRARIES
# ============================================================

library(CellChat)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

options(stringsAsFactors = FALSE)

# ============================================================
# 1. SET PATHS
# ============================================================

base_dir <- "/Users/sudharshiniram/Desktop/SRI Data Analysis"

merged_path <- file.path(base_dir, "Merged_Joined.RDS")
fibroblast_path <- file.path(base_dir, "Fibroblast_Subset.RDS")

mycaf_path_1 <- file.path(base_dir, "myCAF_Subset.RDS")
mycaf_path_2 <- file.path(base_dir, "myCAF_subset.RDS")

cellchat_outdir <- file.path(base_dir, "CellChat_myCAF_Schwann_Pathways")

dir.create(
  cellchat_outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# 2. READ MAIN SEURAT OBJECT
# ============================================================

Merged_Joined <- readRDS(merged_path)

DefaultAssay(Merged_Joined) <- "RNA"

Merged_Joined <- tryCatch(
  JoinLayers(Merged_Joined),
  error = function(e) Merged_Joined
)

# ============================================================
# 3. READ OR CREATE myCAF SUBSET
# ============================================================
# Cluster mapping in fibroblast_subset:
# 0 = Fibroblasts
# 1 = myCAF
# 2 = iCAF
# 3 = pancreatic stellate
# 4 = apCAF

if (file.exists(mycaf_path_1)) {
  
  myCAF_subset <- readRDS(mycaf_path_1)
  print(paste("Loaded myCAF subset from:", mycaf_path_1))
  
} else if (file.exists(mycaf_path_2)) {
  
  myCAF_subset <- readRDS(mycaf_path_2)
  print(paste("Loaded myCAF subset from:", mycaf_path_2))
  
} else {
  
  print("myCAF subset RDS not found. Creating myCAF subset from Fibroblast_Subset.RDS...")
  
  fibroblast_subset <- readRDS(fibroblast_path)
  DefaultAssay(fibroblast_subset) <- "RNA"
  
  fibroblast_subset <- tryCatch(
    JoinLayers(fibroblast_subset),
    error = function(e) fibroblast_subset
  )
  
  myCAF_subset <- subset(
    fibroblast_subset,
    subset = seurat_clusters == 1
  )
  
  saveRDS(
    myCAF_subset,
    mycaf_path_1
  )
  
  print(paste("Saved new myCAF subset to:", mycaf_path_1))
}

DefaultAssay(myCAF_subset) <- "RNA"

myCAF_subset <- tryCatch(
  JoinLayers(myCAF_subset),
  error = function(e) myCAF_subset
)

print("Number of myCAF cells:")
print(ncol(myCAF_subset))

# ============================================================
# 4. ADD SEX METADATA
# ============================================================

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

Merged_Joined$Sex <- unname(
  sex_mapping[as.character(Merged_Joined$orig.ident)]
)

Merged_Joined$Sex <- factor(
  Merged_Joined$Sex,
  levels = c("Female", "Male")
)

# ============================================================
# 5. CREATE CLEAN CELLCHAT CELL TYPE LABELS
# ============================================================

Merged_Joined$cellchat_celltype <- NA_character_

# Label myCAFs using saved myCAF subset cells
mycaf_cells <- intersect(
  colnames(myCAF_subset),
  colnames(Merged_Joined)
)

Merged_Joined$cellchat_celltype[
  colnames(Merged_Joined) %in% mycaf_cells
] <- "myCAFs"

# Label true Schwann cells
if ("true_schwann_subset" %in% colnames(Merged_Joined@meta.data)) {
  
  Merged_Joined$cellchat_celltype[
    as.character(Merged_Joined$true_schwann_subset) %in% c(
      "True Schwann",
      "TRUE",
      "True",
      "true"
    )
  ] <- "Schwann cells"
  
} else {
  
  warning("true_schwann_subset column not found. Using seurat cluster 11 as Schwann fallback.")
  
  Merged_Joined$cellchat_celltype[
    as.character(Merged_Joined$seurat_clusters) == "11"
  ] <- "Schwann cells"
}

# Keep only myCAFs and Schwann cells with clear Sex labels
target_groups <- c(
  "myCAFs",
  "Schwann cells"
)

target_cells <- subset(
  Merged_Joined,
  subset = cellchat_celltype %in% target_groups &
    Sex %in% c("Female", "Male")
)

target_cells$cellchat_celltype <- factor(
  target_cells$cellchat_celltype,
  levels = target_groups
)

target_cells$Sex <- factor(
  target_cells$Sex,
  levels = c("Female", "Male")
)

print("Cell counts by CellChat group and sex:")
print(table(target_cells$cellchat_celltype, target_cells$Sex))

print("Cell counts by patient, CellChat group, and sex:")
print(table(target_cells$orig.ident, target_cells$cellchat_celltype, target_cells$Sex))

# ============================================================
# 6. SPLIT INTO FEMALE AND MALE OBJECTS
# ============================================================

seurat_female <- subset(
  target_cells,
  subset = Sex == "Female"
)

seurat_male <- subset(
  target_cells,
  subset = Sex == "Male"
)

print("Female object cell counts:")
print(table(seurat_female$cellchat_celltype))

print("Male object cell counts:")
print(table(seurat_male$cellchat_celltype))

# ============================================================
# 7. SAFETY CHECK
# ============================================================

min_cells_required <- 10

female_counts <- table(seurat_female$cellchat_celltype)
male_counts <- table(seurat_male$cellchat_celltype)

if (any(female_counts < min_cells_required)) {
  warning("Female subset has a cell type with fewer than 10 cells. Interpret CellChat cautiously.")
}

if (any(male_counts < min_cells_required)) {
  warning("Male subset has a cell type with fewer than 10 cells. Interpret CellChat cautiously.")
}

# ============================================================
# 8. HELPER FUNCTION TO GET RNA DATA
# ============================================================

get_rna_data <- function(seurat_obj) {
  
  data_input <- tryCatch(
    GetAssayData(
      seurat_obj,
      assay = "RNA",
      layer = "data"
    ),
    error = function(e) {
      GetAssayData(
        seurat_obj,
        assay = "RNA",
        slot = "data"
      )
    }
  )
  
  return(data_input)
}

# ============================================================
# 9. CELLCHAT FUNCTION
# ============================================================

run_cellchat <- function(seurat_obj, dataset_name, min_cells_filter = 10) {
  
  print(paste("Running CellChat for", dataset_name, "..."))
  
  DefaultAssay(seurat_obj) <- "RNA"
  
  data.input <- get_rna_data(seurat_obj)
  
  meta <- seurat_obj@meta.data
  
  meta$cellchat_celltype <- droplevels(
    factor(
      meta$cellchat_celltype,
      levels = target_groups
    )
  )
  
  # Keeps patient/sample information
  if (!"samples" %in% colnames(meta)) {
    meta$samples <- as.character(meta$orig.ident)
  }
  
  print(paste("CellChat groups for", dataset_name))
  print(table(meta$cellchat_celltype))
  
  print(paste("Samples represented in", dataset_name))
  print(table(meta$samples))
  
  cellchat <- createCellChat(
    object = data.input,
    meta = meta,
    group.by = "cellchat_celltype"
  )
  
  # Human ligand-receptor database
  CellChatDB <- CellChatDB.human
  cellchat@DB <- CellChatDB
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(
    cellchat,
    type = "triMean",
    raw.use = TRUE
  )
  
  cellchat <- filterCommunication(
    cellchat,
    min.cells = min_cells_filter
  )
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  return(cellchat)
}

# ============================================================
# 10. RUN CELLCHAT FOR FEMALE AND MALE
# ============================================================

cellchat_female <- run_cellchat(
  seurat_obj = seurat_female,
  dataset_name = "Female",
  min_cells_filter = 10
)

cellchat_male <- run_cellchat(
  seurat_obj = seurat_male,
  dataset_name = "Male",
  min_cells_filter = 10
)

cellchat_female <- updateCellChat(cellchat_female)
cellchat_male <- updateCellChat(cellchat_male)

# ============================================================
# 11. SAVE INDIVIDUAL CELLCHAT OBJECTS
# ============================================================

saveRDS(
  cellchat_female,
  file.path(cellchat_outdir, "cellchat_female_myCAF_schwann.rds")
)

saveRDS(
  cellchat_male,
  file.path(cellchat_outdir, "cellchat_male_myCAF_schwann.rds")
)

# ============================================================
# 12. MERGE CELLCHAT OBJECTS
# ============================================================

object.list <- list(
  Female = cellchat_female,
  Male = cellchat_male
)

cellchat_merged <- mergeCellChat(
  object.list,
  add.names = names(object.list)
)

saveRDS(
  cellchat_merged,
  file.path(cellchat_outdir, "cellchat_merged_female_male_myCAF_schwann.rds")
)

# ============================================================
# 13. EXTRACT COMMUNICATION TABLES
# ============================================================

df_female <- subsetCommunication(cellchat_female)
df_male <- subsetCommunication(cellchat_male)

df_female$Sex <- "Female"
df_male$Sex <- "Male"

df_all <- bind_rows(df_female, df_male)

write.csv(
  df_female,
  file.path(cellchat_outdir, "Female_CellChat_All_Communications_myCAF_Schwann.csv"),
  row.names = FALSE
)

write.csv(
  df_male,
  file.path(cellchat_outdir, "Male_CellChat_All_Communications_myCAF_Schwann.csv"),
  row.names = FALSE
)

write.csv(
  df_all,
  file.path(cellchat_outdir, "Female_Male_CellChat_All_Communications_myCAF_Schwann.csv"),
  row.names = FALSE
)

# ============================================================
# 14. KEEP ONLY myCAF <-> SCHWANN COMMUNICATIONS
# ============================================================

mycaf_schwann_all <- df_all %>%
  filter(
    (source == "myCAFs" & target == "Schwann cells") |
      (source == "Schwann cells" & target == "myCAFs")
  ) %>%
  mutate(
    Direction = paste(source, "→", target),
    LR_pair = ifelse(
      "ligand" %in% colnames(.) & "receptor" %in% colnames(.),
      paste(ligand, receptor, sep = " → "),
      interaction_name
    ),
    pathway_name = ifelse(
      is.na(pathway_name),
      "Unknown",
      pathway_name
    ),
    pval_safe = ifelse(
      pval == 0,
      .Machine$double.xmin,
      pval
    ),
    neg_log10_pval = -log10(pval_safe)
  )

if (nrow(mycaf_schwann_all) == 0) {
  stop("No myCAF <-> Schwann interactions were detected. Check source/target labels or CellChat results.")
}

write.csv(
  mycaf_schwann_all,
  file.path(cellchat_outdir, "myCAF_Schwann_CrossTalk_Female_vs_Male.csv"),
  row.names = FALSE
)

print("Number of myCAF <-> Schwann interactions by sex and direction:")
print(
  mycaf_schwann_all %>%
    group_by(Sex, Direction) %>%
    summarise(
      Interaction_Count = n(),
      Total_Strength = sum(prob, na.rm = TRUE),
      Mean_Strength = mean(prob, na.rm = TRUE),
      .groups = "drop"
    )
)

# ============================================================
# 15. PATHWAY SUMMARY TABLE
# ============================================================

pathway_summary <- mycaf_schwann_all %>%
  group_by(Sex, Direction, pathway_name) %>%
  summarise(
    Interaction_Count = n(),
    Total_Prob = sum(prob, na.rm = TRUE),
    Mean_Prob = mean(prob, na.rm = TRUE),
    .groups = "drop"
  )

print("myCAF <-> Schwann pathway summary:")
print(pathway_summary)

write.csv(
  pathway_summary,
  file.path(cellchat_outdir, "myCAF_Schwann_Pathway_Summary_Female_vs_Male.csv"),
  row.names = FALSE
)

# ============================================================
# 16. SELECT TOP PATHWAYS FOR PLOTTING
# ============================================================

top_pathways <- pathway_summary %>%
  group_by(pathway_name) %>%
  summarise(
    Max_Total_Prob = max(Total_Prob, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Max_Total_Prob)) %>%
  slice_head(n = 12) %>%
  pull(pathway_name)

print("Top myCAF <-> Schwann pathways:")
print(top_pathways)

pathway_plot_data <- pathway_summary %>%
  filter(pathway_name %in% top_pathways)

print("Number of pathway rows being plotted:")
print(nrow(pathway_plot_data))

# ============================================================
# 17. PLOT myCAF <-> SCHWANN SIGNALING PATHWAYS
# ============================================================

pathway_plot <- ggplot(
  pathway_plot_data,
  aes(
    x = reorder(pathway_name, Total_Prob),
    y = Total_Prob,
    fill = Sex
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    color = "black",
    linewidth = 0.25,
    width = 0.7
  ) +
  coord_flip() +
  facet_wrap(
    ~Direction,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F8766D",
      "Male" = "#00BFC4"
    )
  ) +
  labs(
    title = "myCAF ↔ Schwann Cell Signaling Pathways",
    x = "Pathway",
    y = "Total communication probability",
    fill = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    strip.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.6
    ),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text = element_text(
      size = 11
    )
  )

# Force plot to show in RStudio Plots pane
graphics.off()
print(pathway_plot)

ggsave(
  file.path(cellchat_outdir, "myCAF_Schwann_Signaling_Pathways_Female_vs_Male.png"),
  plot = pathway_plot,
  width = 11,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(cellchat_outdir, "myCAF_Schwann_Signaling_Pathways_Female_vs_Male.pdf"),
  plot = pathway_plot,
  width = 11,
  height = 7,
  dpi = 300
)

# ============================================================
# 18. OPTIONAL: PLOT myCAF -> SCHWANN ONLY
# ============================================================

mycaf_to_schwann_plot_data <- pathway_plot_data %>%
  filter(Direction == "myCAFs → Schwann cells")

if (nrow(mycaf_to_schwann_plot_data) > 0) {
  
  mycaf_to_schwann_plot <- ggplot(
    mycaf_to_schwann_plot_data,
    aes(
      x = reorder(pathway_name, Total_Prob),
      y = Total_Prob,
      fill = Sex
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.75),
      color = "black",
      linewidth = 0.25,
      width = 0.7
    ) +
    coord_flip() +
    scale_fill_manual(
      values = c(
        "Female" = "#F8766D",
        "Male" = "#00BFC4"
      )
    ) +
    labs(
      title = "myCAF → Schwann Cell Signaling Pathways",
      x = "Pathway",
      y = "Total communication probability",
      fill = "Sex"
    ) +
    theme_classic(base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12, color = "black"),
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11)
    )
  
  print(mycaf_to_schwann_plot)
  
  ggsave(
    file.path(cellchat_outdir, "myCAF_to_Schwann_Signaling_Pathways_Female_vs_Male.png"),
    plot = mycaf_to_schwann_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  ggsave(
    file.path(cellchat_outdir, "myCAF_to_Schwann_Signaling_Pathways_Female_vs_Male.pdf"),
    plot = mycaf_to_schwann_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# ============================================================
# 19. OPTIONAL: PLOT SCHWANN -> myCAF ONLY
# ============================================================

schwann_to_mycaf_plot_data <- pathway_plot_data %>%
  filter(Direction == "Schwann cells → myCAFs")

if (nrow(schwann_to_mycaf_plot_data) > 0) {
  
  schwann_to_mycaf_plot <- ggplot(
    schwann_to_mycaf_plot_data,
    aes(
      x = reorder(pathway_name, Total_Prob),
      y = Total_Prob,
      fill = Sex
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.75),
      color = "black",
      linewidth = 0.25,
      width = 0.7
    ) +
    coord_flip() +
    scale_fill_manual(
      values = c(
        "Female" = "#F8766D",
        "Male" = "#00BFC4"
      )
    ) +
    labs(
      title = "Schwann Cell → myCAF Signaling Pathways",
      x = "Pathway",
      y = "Total communication probability",
      fill = "Sex"
    ) +
    theme_classic(base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12, color = "black"),
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11)
    )
  
  print(schwann_to_mycaf_plot)
  
  ggsave(
    file.path(cellchat_outdir, "Schwann_to_myCAF_Signaling_Pathways_Female_vs_Male.png"),
    plot = schwann_to_mycaf_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  ggsave(
    file.path(cellchat_outdir, "Schwann_to_myCAF_Signaling_Pathways_Female_vs_Male.pdf"),
    plot = schwann_to_mycaf_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# ============================================================
# 20. TOP FEMALE AND MALE myCAF <-> SCHWANN INTERACTIONS
# ============================================================

top_female_interactions <- mycaf_schwann_all %>%
  filter(Sex == "Female") %>%
  arrange(desc(prob)) %>%
  dplyr::select(
    Sex,
    source,
    target,
    Direction,
    pathway_name,
    ligand,
    receptor,
    interaction_name,
    prob,
    pval
  ) %>%
  head(25)

top_male_interactions <- mycaf_schwann_all %>%
  filter(Sex == "Male") %>%
  arrange(desc(prob)) %>%
  dplyr::select(
    Sex,
    source,
    target,
    Direction,
    pathway_name,
    ligand,
    receptor,
    interaction_name,
    prob,
    pval
  ) %>%
  head(25)

print("Top Female myCAF <-> Schwann interactions:")
print(top_female_interactions)

print("Top Male myCAF <-> Schwann interactions:")
print(top_male_interactions)

write.csv(
  top_female_interactions,
  file.path(cellchat_outdir, "Top_25_Female_myCAF_Schwann_Interactions.csv"),
  row.names = FALSE
)

write.csv(
  top_male_interactions,
  file.path(cellchat_outdir, "Top_25_Male_myCAF_Schwann_Interactions.csv"),
  row.names = FALSE
)

# ============================================================
# 21. FINAL SUMMARY
# ============================================================

cat("\n======================================================\n")
cat("CELLCHAT myCAF-SCHWANN PATHWAY ANALYSIS COMPLETE\n")
cat("Results saved to:\n")
cat(cellchat_outdir, "\n")
cat("Main pathway plot saved as:\n")
cat("myCAF_Schwann_Signaling_Pathways_Female_vs_Male.png/pdf\n")
cat("======================================================\n")

cat("\nIMPORTANT INTERPRETATION NOTES:\n")
cat("1. These are CellChat-defined signaling pathways, not GSEA pathways.\n")
cat("2. Pathway strength is based on summed CellChat communication probability.\n")
cat("3. Direction matters: myCAF -> Schwann is different from Schwann -> myCAF.\n")
cat("4. Communication probability is inferred from ligand-receptor expression, not direct physical contact.\n")
cat("======================================================\n")
