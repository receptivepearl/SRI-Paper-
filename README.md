## Methods and Computational Workflow

### Study design and data source

This repository contains the R code used to analyze single-nucleus RNA-sequencing data from 17 treatment-naïve human pancreatic ductal adenocarcinoma tumors. Each tumor sample was treated as a distinct patient-level biological sample. The cohort contained nine female and eight male patients.

The analysis examined:

1. Broad cell-type composition within the PDAC tumor microenvironment.
2. Cancer-associated fibroblast subpopulations.
3. Putative Schwann cells identified through marker co-expression.
4. Sex-associated transcriptional patterns in fibroblasts and Schwann cells.
5. Predicted fibroblast–Schwann and stromal–immune communication using CellChat.

The sample-to-sex assignments used throughout the analysis were:

| Sample    | Sex    |
| --------- | ------ |
| Sample_01 | Male   |
| Sample_02 | Male   |
| Sample_03 | Male   |
| Sample_04 | Female |
| Sample_05 | Female |
| Sample_06 | Female |
| Sample_07 | Male   |
| Sample_08 | Female |
| Sample_09 | Male   |
| Sample_10 | Female |
| Sample_11 | Female |
| Sample_12 | Male   |
| Sample_13 | Female |
| Sample_14 | Male   |
| Sample_15 | Female |
| Sample_16 | Male   |
| Sample_17 | Female |

---

## Software environment

The analysis was performed in R using packages including:

* `Seurat` for object creation, filtering, normalization, dimensionality reduction, clustering, subclustering, differential-expression testing, and visualization.
* `Harmony` through Seurat’s `HarmonyIntegration` method for sample integration.
* `Azimuth` for reference-based pancreatic cell-type annotation support.
* `CellChat` for ligand–receptor communication inference.
* `msigdbr` for retrieving Molecular Signatures Database gene sets.
* `fgsea` for ranked gene-set enrichment analysis.
* `dplyr`, `tidyr`, and `tibble` for data manipulation.
* `ggplot2`, `ggrepel`, and `patchwork` for visualization.

Because the merged object was large, the R vector-memory limit was increased to approximately 32 GB with:

```r
Sys.setenv("R_MAX_VSIZE" = 32000000000)
```

Exact package versions should be recorded with `sessionInfo()` or a package-management file when the final repository version is released.

---

## 1. Creation of sample-level Seurat objects

Raw 10x Genomics gene-expression matrices were loaded separately for each of the 17 samples using `Read10X()`.

A separate Seurat object was created for each sample with `CreateSeuratObject()`:

```r
CreateSeuratObject(
  counts = Read10X(sample_directory),
  project = sample_name,
  min.cells = 3,
  min.features = 200
)
```

The initial creation thresholds had the following effects:

* A gene had to be detected in at least three nuclei to remain in the expression matrix.
* A nucleus had to contain at least 200 detected genes to enter the initial Seurat object.

These are initial matrix-level inclusion thresholds. More restrictive nucleus-level QC was subsequently applied to the merged object.

---

## 2. Sample merging

The 17 sample-level Seurat objects were merged into one combined object using `merge()`.

Sample-specific cell prefixes, `S01` through `S17`, were added using `add.cell.ids`. This preserved the sample of origin for every nucleus and reduced the risk of duplicated barcode names across samples.

The original sample identities were retained in the `orig.ident` metadata field and were subsequently used to:

* Perform Harmony integration across samples.
* Assign patient sex.
* Calculate patient-level cell counts and proportions.
* Retain sample information during CellChat analysis.
* Generate patient-level summary tables.

The merged, unfiltered object was saved before downstream processing to preserve an intermediate checkpoint.

---

## 3. Quality control and nucleus filtering

### Mitochondrial transcript percentage

The percentage of reads assigned to mitochondrial genes was calculated with:

```r
PercentageFeatureSet(
  Merged,
  pattern = "^MT-",
  col.name = "percent.mt"
)
```

Genes beginning with `MT-` were treated as mitochondrial genes.

### Final nucleus-level filtering thresholds

Nuclei were retained only when they met all of the following criteria:

```r
nFeature_RNA > 200
nFeature_RNA < 3000
percent.mt < 5
nCount_RNA > 500
nCount_RNA < 10000
```

Therefore, the final filtering criteria were:

| QC metric                                  |   Retained range | Purpose                                                                                                                |
| ------------------------------------------ | ---------------: | ---------------------------------------------------------------------------------------------------------------------- |
| Detected genes per nucleus, `nFeature_RNA` |  >200 and <3,000 | Removed very low-complexity nuclei and unusually complex profiles that could represent multiplets or atypical droplets |
| Total RNA counts per nucleus, `nCount_RNA` | >500 and <10,000 | Removed nuclei with very low transcript recovery and nuclei with unusually large libraries                             |
| Mitochondrial RNA percentage, `percent.mt` |              <5% | Removed profiles with elevated mitochondrial signal, which may indicate damaged or low-quality nuclei                  |

The lower gene and count thresholds excluded low-information droplets. The upper gene and count thresholds reduced the contribution of unusually complex profiles that could reflect doublets, multiplets, or technical outliers.

The mitochondrial cutoff was relatively stringent. This was appropriate for a single-nucleus dataset because nuclear RNA profiles are generally expected to have lower mitochondrial transcript representation than whole-cell single-cell RNA-sequencing profiles.



## 4. Normalization and feature selection

The filtered merged object was normalized using Seurat’s global-scaling `LogNormalize` procedure:

```r
NormalizeData(
  Merged,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
```

For each nucleus, gene counts were divided by the total counts in that nucleus, multiplied by 10,000, and log-transformed.

The 2,000 most variable genes were then identified with the variance-stabilizing transformation method:

```r
FindVariableFeatures(
  Merged,
  selection.method = "vst",
  nfeatures = 2000
)
```

These highly variable genes were used to represent the major biological differences among nuclei during dimensionality reduction.

The normalized expression matrix was centered and scaled using:

```r
ScaleData(Merged)
```

No regression variables were specified in `ScaleData()`. Therefore, mitochondrial percentage, total counts, patient age, sex, and other metadata variables were not explicitly regressed out during scaling.

Sex was intentionally not used as an integration or regression variable because sex-associated biological variation was an outcome of interest.

---

## 5. Principal-component analysis and sample integration

Principal-component analysis was initially performed using the 2,000 variable genes:

```r
RunPCA(
  Merged,
  features = VariableFeatures(Merged)
)
```

Sample integration was then conducted using Harmony through Seurat’s `IntegrateLayers()` framework:

```r
IntegrateLayers(
  object = Merged,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony"
)
```

Harmony was used to reduce sample-associated technical variation while preserving shared biological structure across the 17 tumors.

The initial integrated workflow used:

* Harmony dimensions 1–20.
* A shared nearest-neighbor graph constructed from the Harmony reduction.
* Clustering resolution of 0.5.
* UMAP based on Harmony dimensions 1–20.

```r
FindNeighbors(
  Merged,
  reduction = "harmony",
  dims = 1:20
)

FindClusters(
  Merged,
  resolution = 0.5
)

RunUMAP(
  Merged,
  reduction = "harmony",
  dims = 1:20
)
```

Because the analysis used Seurat v5 layered assays, `JoinLayers()` was used before procedures that required a unified RNA-expression layer, including marker testing and downstream analyses.

---

## 6. Final broad clustering and UMAP settings

For the final broad cell-type representation, PCA was rerun with 50 principal components:

```r
RunPCA(
  Merged,
  features = VariableFeatures(Merged),
  npcs = 50
)
```

An elbow plot of the first 50 principal components was generated to inspect the amount of variation represented across components.

The final broad clustering workflow used the first 30 principal components:

```r
pc_use <- 1:30
```

A nearest-neighbor graph was constructed with `FindNeighbors()`, and broad clusters were identified with a resolution of 0.2:

```r
FindNeighbors(
  Merged,
  dims = 1:30
)

FindClusters(
  Merged,
  resolution = 0.2
)
```

The lower clustering resolution was selected to generate broad, biologically interpretable tumor-microenvironment populations rather than a large number of fine-grained subclusters.

The final broad UMAP was generated with:

```r
RunUMAP(
  Merged,
  dims = 1:30,
  n.neighbors = 30,
  min.dist = 0.3
)
```

The relevant parameters were:

| Parameter                                | Value |
| ---------------------------------------- | ----: |
| PCA components calculated                |    50 |
| PCA dimensions used for broad clustering |  1–30 |
| Clustering resolution                    |   0.2 |
| UMAP nearest neighbors                   |    30 |
| UMAP minimum distance                    |   0.3 |
| Random seed used before PCA              |  1234 |

Cluster assignments were stored in metadata before manual annotation so that the original computational cluster identities were retained.

---

## 7. Identification of cluster markers

Positive marker genes for the broad computational clusters were identified using `FindAllMarkers()`.

To reduce memory use, testing was limited to a maximum of 500 nuclei per identity:

```r
FindAllMarkers(
  object = Merged,
  only.pos = TRUE,
  max.cells.per.ident = 500
)
```

This analysis was used to support cluster annotation rather than to test sex-associated differences.

For each cluster, marker genes were ranked by average log2 fold change. The top 30 markers per cluster were exported for inspection.

Because of the `max.cells.per.ident = 500` parameter, the marker calculation used a downsampled maximum of 500 nuclei per cluster. This helped control computational memory use and prevented the largest clusters from completely dominating marker identification.

---

## 8. Broad cell-type annotation

Broad clusters were annotated by jointly evaluating:

1. Cluster-specific marker-gene results.
2. Canonical lineage-marker expression.
3. Feature plots.
4. Dot plots.
5. Cluster location and separation on the UMAP.
6. Azimuth reference-mapping predictions as supporting evidence.

Azimuth mapping was performed with the human pancreas reference:

```r
RunAzimuth(
  query = Merged,
  reference = "pancreasref"
)
```

The resulting predicted annotation and mapping score were visualized. Because the reference primarily represents pancreatic cell states and may not fully capture malignant or disease-associated populations, Azimuth predictions were used as supporting evidence rather than as the sole basis for annotation.

### Representative annotation markers

Broad populations were evaluated using marker panels including:

| Population                             | Representative markers                                                              |
| -------------------------------------- | ----------------------------------------------------------------------------------- |
| Epithelial, ductal, or PDAC-associated | `EPCAM`, `KRT19`, `KRT8`, `KRT18`, `KRT7`, `SOX9`, `MUC1`, `MMP7`, `TSPAN8`, `LCN2` |
| Acinar                                 | `PRSS1`, `PRSS2`, `CPA1`, `CPA2`, `CTRB1`, `CTRB2`, `CELA3A`, `PNLIP`               |
| Endocrine                              | `CHGA`, `CHGB`, `INS`, `GCG`, `SST`, `PPY`, `PCSK1`, `PCSK2`                        |
| Fibroblast or CAF                      | `COL1A1`, `COL1A2`, `COL3A1`, `DCN`, `LUM`, `FAP`, `PDGFRA`, `SPARC`, `FN1`         |
| Endothelial                            | `PECAM1`, `VWF`, `CDH5`, `PLVAP`, `KDR`, `FLT1`, `EMCN`                             |
| Pericyte or smooth muscle              | `RGS5`, `PDGFRB`, `CSPG4`, `MCAM`, `NOTCH3`, `ACTA2`                                |
| Myeloid or macrophage                  | `LYZ`, `LST1`, `TYROBP`, `FCER1G`, `CD68`, `C1QA`, `C1QB`, `C1QC`, `MRC1`, `CD163`  |
| T or NK cell                           | `CD3D`, `CD3E`, `CD4`, `CD8A`, `NKG7`, `GNLY`                                       |
| B cell                                 | `MS4A1`, `CD79A`, `CD19`, `CD27`                                                    |
| Schwann or peripheral glial            | `SOX10`, `S100B`, `PMP22`, `MPZ`, `ERBB3`, `PLP1`, `GFAP`                           |

Mixed-lineage marker panels were also inspected to identify clusters with possible doublet or ambiguous signals. However, a formal computational doublet-removal procedure was not run in the supplied code.

---

## 9. Fibroblast isolation and subclustering

The broad fibroblast cluster was isolated from the merged object for a dedicated fibroblast analysis.

Before subclustering, the existing PCA, Harmony, and UMAP reductions were removed so that the fibroblast population could be reanalyzed independently:

```r
fibroblast_subset@reductions$pca <- NULL
fibroblast_subset@reductions$harmony <- NULL
fibroblast_subset@reductions$umap <- NULL
```

Within the fibroblast subset:

1. The 2,000 most variable genes were recalculated.
2. Expression values were rescaled.
3. PCA was rerun.
4. A new nearest-neighbor graph was generated.
5. Fibroblasts were reclustered.
6. A fibroblast-specific UMAP was generated.

```r
FindVariableFeatures(
  fibroblast_subset,
  selection.method = "vst",
  nfeatures = 2000
)

ScaleData(fibroblast_subset)

RunPCA(
  fibroblast_subset,
  features = VariableFeatures(fibroblast_subset)
)

FindNeighbors(
  fibroblast_subset,
  reduction = "pca",
  dims = 1:15
)

FindClusters(
  fibroblast_subset,
  resolution = 0.2
)

RunUMAP(
  fibroblast_subset,
  reduction = "pca",
  dims = 1:15,
  reduction.name = "sub_umap"
)
```

The fibroblast-subclustering parameters were:

| Parameter                    |      Value |
| ---------------------------- | ---------: |
| Variable genes               |      2,000 |
| PCA dimensions used          |       1–15 |
| Clustering resolution        |        0.2 |
| Reduction used for neighbors |        PCA |
| Fibroblast UMAP name         | `sub_umap` |

Harmony was not rerun within the fibroblast subset in the supplied subclustering script. The fibroblast-specific graph and UMAP were based on the fibroblast PCA reduction.

---

## 10. Fibroblast-subtype annotation

Fibroblast subclusters were assigned using coordinated expression of canonical CAF markers rather than any single marker alone.

### Myofibroblastic CAFs

Myofibroblastic CAFs were evaluated using markers associated with contractility and extracellular-matrix production, including:

* `ACTA2`
* `TAGLN`
* `MYL9`
* `COL1A1`
* `COL1A2`
* `POSTN`

### Inflammatory CAFs

Inflammatory CAFs were evaluated using cytokine, chemokine, and fibroblast-associated markers including:

* `IL6`
* `CXCL12`
* `CXCL14`
* `LIF`
* `CFD`
* `PDGFRA`

### Antigen-presenting CAFs

Antigen-presenting CAFs were evaluated using MHC class II and antigen-presentation markers including:

* `HLA-DRA`
* `HLA-DRB1`
* `HLA-DPA1`
* `HLA-DPB1`
* `CD74`
* `CIITA`

### Putative PSC-like fibroblasts

One fibroblast cluster was annotated as a **putative pancreatic stellate cell-like fibroblast population** based on its fibroblast identity, marker-expression pattern, and enrichment for PSC-associated genes such as `PLXDC1`.

The term “putative PSC-like” is used because this population was identified computationally from transcriptomic marker expression and was not independently validated through spatial localization, protein-level assays, lineage tracing, or external reference mapping in the supplied scripts.

Feature plots and dot plots were used to compare marker expression among the computational fibroblast subclusters before labels were assigned.

---

## 11. Cell-type abundance calculations

Cell counts and proportions were calculated from Seurat metadata.

For broad cell types, the number of nuclei belonging to each annotated population was calculated separately for every patient. Proportions were calculated relative to the relevant patient-level denominator, such as:

* Total recovered nuclei from that patient for broad cell-type abundance.
* Total fibroblast nuclei from that patient for fibroblast-subtype abundance.

Sex comparisons were based on patient-level values rather than treating every nucleus as an independent sample when the scripts explicitly summarized results by `orig.ident`.

For the putative PSC-like cluster, the primary interpretable abundance measure was its percentage among fibroblast nuclei within each patient:

[
\text{PSC-like fraction} =
\frac{\text{PSC-like fibroblast nuclei}}
{\text{total fibroblast nuclei}}
\times 100
]

The numerator and denominator should both be reported in repository output tables so the compositional calculation can be reproduced.

---

## 12. Identification of the final Schwann-cell subset

Because peripheral glial cells were rare and could overlap transcriptionally with other stromal populations, the final Schwann subset was defined using co-expression of two Schwann-associated genes:

* `ERBB3`
* `PMP22`

Normalized RNA-expression values were extracted from the joined Seurat object with `FetchData()`.

A nucleus was labeled “True Schwann” when both genes had expression greater than the specified expression threshold:

```r
expr_data$ERBB3 > expr_threshold &
expr_data$PMP22 > expr_threshold
```

The resulting nuclei were stored in a metadata field called `true_schwann_subset` and isolated into a separate Seurat object.

This procedure required simultaneous expression of both genes rather than expression of either marker alone, increasing specificity for the final computational Schwann-cell definition.

The resulting subset was further inspected for expression of supporting Schwann or glial markers such as:

* `MPZ`
* `S100B`
* `SOX10`
* `PLP1`
* `GFAP`
* `NGFR`
* `POU3F1`
* `EGR2`

The final Schwann object was saved as:

```text
true_schwann_ERBB3_PMP22_subset.rds
```

### Per-patient Schwann abundance

For each patient, the scripts calculated:

* Total retained nuclei.
* Number of ERBB3-positive/PMP22-positive nuclei.
* Percentage of total nuclei classified as Schwann cells.

The percentage was calculated as:

[
\text{Schwann percentage} =
\frac{\text{ERBB3}^{+}\text{/PMP22}^{+}\text{ nuclei}}
{\text{total retained nuclei}}
\times 100
]

Because Schwann cells were rare, all Schwann-cell results should be interpreted cautiously and accompanied by per-patient nucleus counts.

---

## 13. Sex metadata assignment

Sex labels were assigned using `orig.ident` and a manually specified sample-to-sex mapping.

Only nuclei with an unambiguous female or male assignment were retained for sex-stratified analyses:

```r
subset(
  object,
  subset = Sex %in% c("Male", "Female")
)
```

Sex was stored as a factor, generally with female followed by male as the level order. In differential-expression analyses where `ident.1 = "Male"` and `ident.2 = "Female"`, positive average log2 fold changes therefore indicated higher expression in male nuclei, whereas negative values indicated higher expression in female nuclei.

---

## 14. Nucleus-level differential-expression analyses

### Fibroblast differential expression

Sex-associated fibroblast expression was explored using Seurat’s `FindMarkers()` function with the Wilcoxon rank-sum test.

The volcano-plot analysis used:

```r
FindMarkers(
  object = fibroblast_subset,
  ident.1 = "Male",
  ident.2 = "Female",
  group.by = "Sex",
  logfc.threshold = 0.25,
  min.pct = 0.10
)
```

Genes were visually classified as sex-associated in the volcano plot when they met both:

* Benjamini–Hochberg-adjusted `P < 0.05`.
* Absolute average log2 fold change greater than 0.5.

For pathway ranking, differential expression was rerun without an initial fold-change or expression-percentage threshold:

```r
FindMarkers(
  object = fibroblast_subset,
  ident.1 = "Male",
  ident.2 = "Female",
  group.by = "Sex",
  logfc.threshold = 0,
  min.pct = 0,
  only.pos = FALSE
)
```

This retained a broad ranked gene list for enrichment analysis.

### Schwann-cell differential expression

Schwann-cell differential expression was also performed using Seurat’s Wilcoxon test.

Before testing, genes were restricted to those with average normalized expression greater than 0.5 in both female and male Schwann cells:

```r
AverageExpression(
  true_schwann_subset,
  assays = "RNA",
  group.by = "Sex",
  slot = "data"
)
```

A gene passed the expression filter when:

```r
average expression in males > 0.5
average expression in females > 0.5
```

Differential expression was then calculated using:

```r
FindMarkers(
  object = true_schwann_subset,
  ident.1 = "Male",
  ident.2 = "Female",
  group.by = "Sex",
  features = genes_expr_gt_0.5,
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.05,
  only.pos = FALSE
)
```

Genes shown as sex-associated in the Schwann volcano plot met:

* Adjusted `P < 0.05`.
* Absolute average log2 fold change greater than 0.5.

### Interpretation limitation

These Seurat Wilcoxon analyses treat nuclei as observations and therefore do not, by themselves, model the patient as the biological replicate. They are best described as **nucleus-level exploratory or descriptive differential-expression analyses**.

They should not be presented as patient-level inferential evidence unless supported by a separate patient-aware analysis, such as pseudobulk edgeR, DESeq2, or limma-voom modeling.

Where the manuscript reports a separate patient-level fibroblast pseudobulk analysis, that analysis should be documented in its own README subsection using its corresponding script. It should not be implied that the Seurat `FindMarkers()` results are pseudobulk results.

---

## 15. Ranked gene-set enrichment analysis

Gene-set enrichment analysis was performed with `fgsea` using the MSigDB Hallmark collection obtained through `msigdbr`.

The Hallmark collection was loaded for Homo sapiens:

```r
msigdbr(
  species = "Homo sapiens",
  collection = "H"
)
```

For fibroblasts, genes were ranked using a signed adjusted-significance score:

[
\text{rank score} =
\operatorname{sign}(\text{average log}*2\text{FC})
\times
-\log*{10}(\text{adjusted }P)
]

This ranking incorporated both the direction of the male-versus-female fold change and the adjusted statistical evidence from the nucleus-level test.

Genes were sorted from the most strongly male-associated positive rank to the most strongly female-associated negative rank.

`fgsea()` was run with:

```r
fgsea(
  pathways = pathways_list,
  stats = gene_ranks,
  minSize = 10,
  maxSize = 500,
  eps = 0
)
```

The parameters were:

| Parameter            |           Value |
| -------------------- | --------------: |
| Gene-set collection  | MSigDB Hallmark |
| Species              |    Homo sapiens |
| Minimum pathway size |        10 genes |
| Maximum pathway size |       500 genes |
| `eps`                |               0 |
| Random seed          |              42 |

Positive normalized enrichment scores corresponded to pathways enriched toward genes with higher expression in male nuclei. Negative normalized enrichment scores corresponded to pathways enriched toward genes with higher expression in female nuclei.

The output included:

* Enrichment score.
* Normalized enrichment score.
* Nominal P value.
* Adjusted P value.
* Gene-set size.
* Leading-edge genes.

Pathways should only be described as statistically enriched when they satisfy the prespecified false-discovery-rate threshold. Pathways failing that threshold should be labeled as exploratory trends.

Because the ranking statistics were derived from nucleus-level differential-expression tests, the enrichment analysis inherits the same experimental-unit limitation and should not be represented as a patient-level pathway test.

---

## 16. CellChat analysis

CellChat was used to infer potential ligand–receptor communication among selected PDAC tumor-microenvironment populations.

Separate CellChat objects were constructed for female and male nuclei. Depending on the analysis, the selected groups included:

* Fibroblasts and Schwann cells.
* iCAFs and Schwann cells.
* myCAFs and Schwann cells.
* Multiple CAF subtypes and Schwann cells.
* Fibroblast, Schwann, macrophage, dendritic-cell, T-cell, and B-cell populations.

### Input data

The normalized RNA-expression layer was extracted from the Seurat object:

```r
GetAssayData(
  seurat_obj,
  assay = "RNA",
  layer = "data"
)
```

The metadata supplied to CellChat included:

* The CellChat group label.
* `orig.ident`.
* A `samples` column derived from `orig.ident`.

Adding `samples` preserved the patient or sample identity within the CellChat object, although the pooled female and male networks still combined nuclei from multiple tumors during network inference.

CellChat objects were created with:

```r
createCellChat(
  object = data.input,
  meta = meta,
  group.by = "cellchat_celltype"
)
```

### Ligand–receptor database

The full human CellChat database was used:

```r
CellChatDB.human
```

This included signaling interactions from categories such as:

* Secreted signaling.
* Extracellular-matrix–receptor interactions.
* Cell–cell contact.

The analysis was not restricted to secreted signaling alone.

### Communication inference

The following CellChat workflow was applied:

```r
subsetData()
identifyOverExpressedGenes()
identifyOverExpressedInteractions()
computeCommunProb(
  type = "triMean",
  raw.use = TRUE
)
filterCommunication(
  min.cells = 10
)
computeCommunProbPathway()
aggregateNet()
```

Communication probability was calculated with the trimean method. Interactions were retained only when the relevant communicating group contained at least 10 nuclei.

The minimum-cell filter helped remove inferred interactions involving extremely small groups, but it does not eliminate uncertainty caused by unequal cell recovery, unequal patient contribution, or rare populations.

`projectData()` was not used because it was unavailable in the installed CellChat version used by the scripts.

### Female–male comparisons

Female and male CellChat objects were generated separately and merged with `mergeCellChat()` for descriptive comparison.

The following network summaries were extracted:

* Number of inferred ligand–receptor interactions.
* Aggregate network interaction strength.
* Signaling pathway.
* Source population.
* Target population.
* Ligand.
* Receptor.
* Communication probability.

Directional interactions were retained separately, including:

* Fibroblast-to-Schwann signaling.
* Schwann-to-fibroblast signaling.
* iCAF-to-Schwann signaling.
* Schwann-to-iCAF signaling.
* myCAF-to-Schwann signaling.
* Schwann-to-myCAF signaling.
* Stromal-to-immune signaling.

Interaction counts and aggregate strengths were calculated by summing the corresponding CellChat network outputs.

### Interpretation of pooled networks

The pooled female and male CellChat networks are descriptive predictions. They combine nuclei from multiple patients and may be affected by:

* Unequal numbers of nuclei per sex.
* Unequal contribution of individual tumors.
* Differences in cell-type abundance.
* Rare Schwann-cell recovery.
* CellChat’s model assumptions.
* The absence of direct spatial or protein-level validation.

Pooled networks therefore indicate possible communication patterns but do not independently establish statistically significant sex differences.

---

## 17. Patient-level CellChat comparisons

For patient-aware communication summaries, separate CellChat networks were generated for individual tumors rather than treating all female or male nuclei as a single biological replicate.

For each eligible patient network, prespecified summary measures were extracted, including:

* Inferred interaction count.
* Aggregate interaction strength.

When repeated runs were available for a patient, their metrics were averaged before sex comparison.

Patient-level female and male values were compared using an exact permutation test, with the Wilcoxon rank-sum test used as a sensitivity analysis. Benjamini–Hochberg correction was applied across the prespecified network metrics.

These patient-level comparisons are distinct from pseudobulk differential-expression analysis. CellChat communication probabilities were not generated through an edgeR-style pseudobulk model.

---

## 18. Visualization and exported outputs

The repository scripts generated and exported:

* QC-filtered and integrated Seurat objects.
* Broad-cluster UMAPs.
* Marker-gene spreadsheets.
* Broad lineage-marker dot plots.
* Feature plots for canonical markers.
* Azimuth predicted annotations and mapping scores.
* Fibroblast-specific UMAPs.
* CAF-subtype marker plots.
* Patient-level cell-count and proportion tables.
* Fibroblast and Schwann volcano plots.
* Hallmark enrichment tables and plots.
* CellChat interaction tables.
* CellChat interaction-count and interaction-strength figures.
* Ligand–receptor heatmaps and pathway summaries.
* Per-patient Schwann-cell counts and percentages.

Most manuscript-ready figures were saved at 300 dots per inch using `ggsave()`.

---
