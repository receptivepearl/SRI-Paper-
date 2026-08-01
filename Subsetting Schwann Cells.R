# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# ============================================================
# 1. SETUP AND DATA EXTRACTION
# ============================================================
# Define your target list and the desired percentage bounds "NRXN1"
target_genes <- c("PLP1", "CDH19", "SOX10", "S100B", "GFAP", 
                  "NGFR", "PMP22", "MPZ", "GPM6B", "SCN7A", "PNLIP", "PNLIPRP1", "ERBB3")

min_pct <- 0.5
max_pct <- 6.5
expr_threshold <- 0 # Any expression > 0 counts as positive

# Verify which genes actually exist in the dataset to prevent errors
valid_genes <- intersect(target_genes, rownames(Merged_Joined))
print(paste("Valid genes found:", length(valid_genes), "out of", length(target_genes)))

# Fetch expression data for all valid genes ONCE (Highly optimized for speed)
expr_data <- FetchData(Merged_Joined, vars = c("orig.ident", valid_genes), layer = "data")

# Create a rapid logical matrix (TRUE if expression > threshold)
expr_matrix <- expr_data[, valid_genes] > expr_threshold


# ============================================================
# 2. GENERATE ALL 2- AND 3-GENE COMBINATIONS
# ============================================================
combos_2 <- combn(valid_genes, 2, simplify = FALSE)
combos_3 <- combn(valid_genes, 3, simplify = FALSE)
all_combos <- c(combos_2, combos_3)

print(paste("Testing", length(all_combos), "unique gene combinations..."))


# ============================================================
# 3. HIGH-SPEED COMBINATORIAL SCREENING
# ============================================================
successful_combos <- list()

for (i in seq_along(all_combos)) {
  current_combo <- all_combos[[i]]
  combo_name <- paste(current_combo, collapse = " + ")
  
  # A cell is positive ONLY if it expresses ALL genes in this specific combination
  # rowSums checks how many genes in the combo are TRUE. It must equal the combo length.
  is_positive <- rowSums(expr_matrix[, current_combo, drop = FALSE]) == length(current_combo)
  
  # Calculate per-sample percentages
  temp_summary <- data.frame(orig.ident = expr_data$orig.ident, Positive = is_positive) %>%
    group_by(orig.ident) %>%
    summarise(
      Total = n(),
      Pos = sum(Positive),
      Pct = (Pos / Total) * 100,
      .groups = "drop"
    )
  
  # Calculate the average percentage across all 17 samples
  mean_sample_pct <- mean(temp_summary$Pct)
  
  # If the average falls within your target 0.5% - 5.5% window, save it!
  if (mean_sample_pct >= min_pct && mean_sample_pct <= max_pct) {
    successful_combos[[combo_name]] <- data.frame(
      Combination = combo_name,
      Gene_Count = length(current_combo),
      Mean_Pct = round(mean_sample_pct, 3),
      Max_Sample_Pct = round(max(temp_summary$Pct), 3),
      Min_Sample_Pct = round(min(temp_summary$Pct), 3)
    )
  }
}


# ============================================================
# 4. FORMAT AND PRINT RESULTS
# ============================================================
# Combine all successful results into one clean table
if (length(successful_combos) > 0) {
  final_results <- bind_rows(successful_combos) %>%
    arrange(desc(Mean_Pct)) # Sort highest to lowest within your window
  
  print("======================================================")
  print(paste("SUCCESS: Found", nrow(final_results), "combinations matching your criteria!"))
  print("======================================================")
  
  # Print the top 15 best combinations to the console
  print(head(final_results, 15))
  
  # Save the full results table to your computer so you can review them all
  write.csv(final_results, "/Users/sudharshiniram/Desktop/SRI Data Analysis/Schwann_Marker_Combos.csv", row.names = FALSE)
  print("Full results saved to desktop as 'Schwann_Marker_Combos.csv'")
  
} else {
  print("======================================================")
  print("NO COMBINATIONS FOUND.")
  print("Try lowering your min_pct or reducing the expr_threshold.")
  print("======================================================")
}


# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# ============================================================
# 1. EXTRACT DATA FOR THE WINNING COMBINATION
# ============================================================
# Define the winning markers
target_genes <- c("PMP22", "SCN7A")

# Fetch data for these specific genes
expr_data <- FetchData(Merged_Joined, vars = c("orig.ident", target_genes), layer = "data")


# ============================================================
# 2. CALCULATE STRICT CO-EXPRESSION
# ============================================================
schwann_summary <- expr_data %>%
  group_by(orig.ident) %>%
  summarise(
    Total_Cells = n(),
    
    # The cell MUST have > 0 expression for BOTH genes
    Schwann_Positive_Cells = sum(
      PMP22 > 0 & 
        SCN7A > 0, 
      na.rm = TRUE
    ),
    
    Percentage_Schwann = (Schwann_Positive_Cells / Total_Cells) * 100
  ) %>%
  ungroup()


# ============================================================
# 3. RENAME SAMPLES FOR PUBLICATION GRAPH
# ============================================================
patient_mapping <- c(
  "Sample_01" = "P01", "Sample_02" = "P02", "Sample_03" = "P03",
  "Sample_04" = "P04", "Sample_05" = "P05", "Sample_06" = "P06",
  "Sample_07" = "P07", "Sample_08" = "P08", "Sample_09" = "P09",
  "Sample_10" = "P10", "Sample_11" = "P11", "Sample_12" = "P12",
  "Sample_13" = "P13", "Sample_14" = "P14", "Sample_15" = "P15",
  "Sample_16" = "P16", "Sample_17" = "P17"
)

schwann_summary$PatientID <- patient_mapping[as.character(schwann_summary$orig.ident)]
schwann_summary$PatientID <- factor(schwann_summary$PatientID, levels = sprintf("P%02d", 1:17))


# ============================================================
# 4. PRINT METRICS TO CONSOLE
# ============================================================
print("======================================================")
print("FINAL SCHWANN SUBSET: PMP22 & SCN7A Co-Expression")
print("======================================================")
print(as.data.frame(schwann_summary[, c("PatientID", "Total_Cells", "Schwann_Positive_Cells", "Percentage_Schwann")]))


# ============================================================
# 5. GENERATE PUBLICATION-READY BAR CHART
# ============================================================
schwann_graph <- ggplot(schwann_summary, aes(x = PatientID, y = Percentage_Schwann)) +
  geom_bar(stat = "identity", fill = "navy", color = "black", linewidth = 0.5) +
  
  # Ensure the y-axis starts exactly at 0
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  
  # Professional Titles
  ggtitle("Proportion of PMP22+ / SCN7A+ Schwann Cells") +
  labs(x = "Patient Sample", y = "Percentage of Total Cells (%)") +
  
  # Publication-standard formatting
  theme_classic(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 10)),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 10)),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
  )

# Print the final plot
print(schwann_graph)

