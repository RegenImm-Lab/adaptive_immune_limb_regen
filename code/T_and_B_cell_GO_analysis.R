library(ontologyIndex)
library(dplyr)
library(stringr)
library(dplyr)
library(forcats)
library(clusterProfiler)
library(edgeR)
library(tidytext)
library(Seurat)
library(patchwork)
library(hdf5r)
library(ggplot2)
library(devtools)
library(CyteTypeR)
library(viridis)
library(scCustomize)

#load in data
set.seed(1234)

setwd('/Users/nickleigh/OneDrive - Lund University/Dokument/Data/2024/single_cell/inDrops/data/h5Seurat_summed/')

files <- c("14dpa1_so_output_filtered_summed_seurat.h5", "14dpa2_so_output_filtered_summed_seurat.h5", 
           "23dpa1_so_output_filtered_summed_seurat.h5", "23dpa2_so_output_filtered_summed_seurat.h5", 
           "23dpa3_so_output_filtered_summed_seurat.h5", "23dpa4_so_output_filtered_summed_seurat.h5", 
           "23dpa5_so_output_filtered_summed_seurat.h5", "23dpa6_so_output_filtered_summed_seurat.h5", 
           "3dpa1_so_output_filtered_summed_seurat.h5", "3dpa2_so_output_filtered_summed_seurat.h5", 
           "3dpa3_so_output_filtered_summed_seurat.h5", "Contra1_so_output_filtered_seurat.h5", 
           "Contra2_so_output_filtered_seurat.h5", "Contra3_so_output_filtered_seurat.h5", 
           "Intact1_so_output_filtered_seurat.h5", "Intact2_so_output_filtered_seurat.h5")

for (file in files) {
  obj_name <- paste0("s", sub("_so_output_filtered(_summed)?_seurat\\.h5$", "", file))  
  dgCMatrix_object <- Read10X_h5(file, use.names = TRUE, unique.features = TRUE)
  seurat_object <- CreateSeuratObject(counts = dgCMatrix_object)
  assign(obj_name, seurat_object, envir = .GlobalEnv)  
}


# Batch 1: 23s14dpa1# Batch 1: 23dpa4, 5, and 6, plus 3dpa1-3
batch1 <- merge(
  x = s23dpa4,
  y = list(s23dpa5,s23dpa6 ,s3dpa1, s3dpa2, s3dpa3),
  add.cell.ids = c("23dpa4", "23dpa5", "23dpa6", "3dpa1", "3dpa2", "3dpa3"),
  project = getOption(x = "batch1", default = "indrop")
)


# Batch 2: 23dpa1-3 and 14dpa samples
batch2 <- merge(
  x = s23dpa1,
  y = list(s23dpa2, s23dpa3, s14dpa1, s14dpa2),
  add.cell.ids = c("23dpa1", "23dpa2", "23dpa3", "14dpa1", "14dpa2"),
  project = getOption(x = "batch2", default = "indrop")
)

# Batch 3: Intact and Contra
batch3 <- merge(
  x = sIntact1,
  y = list(sIntact2, sContra1, sContra2, sContra3),
  add.cell.ids = c("Intact1", "Intact2", "Contra1", "Contra2", "Contra3"),
  project = getOption(x = "batch3", default = "indrop")
)

batch1 <- JoinLayers(batch1)
#just make sure this will only return mito 
grep(paste('^MT\\.'), rownames(batch1), value = TRUE)

#add percent mito
batch1[["percent.mt"]] <- PercentageFeatureSet(batch1, pattern = "^MT\\.")


#QC plots
VlnPlot(batch1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'orig.ident')
plot1 <- FeatureScatter(batch1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(batch1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#fitler batch1
dim(batch1)
#[1] 29097 26418
batch1 <- subset(batch1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 8)
dim(batch1)
#[1] 29097 19087


#now batch 2 QC
batch2 <- JoinLayers(batch2)
batch2[["percent.mt"]] <- PercentageFeatureSet(batch2, pattern = "^MT\\.")

#QC plots
VlnPlot(batch2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'orig.ident')
plot1 <- FeatureScatter(batch2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(batch2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#fitler batch1
dim(batch2)
#[1] 28799 16039
batch2 <- subset(batch2, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 8)
dim(batch2)
#[1] 28799 11448

#batch 3 QC
batch3 <- JoinLayers(batch3)
batch3[["percent.mt"]] <- PercentageFeatureSet(batch3, pattern = "^MT\\.")


#QC plots
VlnPlot(batch3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'orig.ident')
plot1 <- FeatureScatter(batch3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(batch3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#fitler batch1
dim(batch3)
#[1] 28702 14651
batch3 <- subset(batch3, subset = nFeature_RNA > 200 & nFeature_RNA < 3750 & percent.mt < 5)
dim(batch3)
#[1] 28702 12548

#all three are cleaned now 

#merge all
#first add batch info
batch1$batch <- 'batch1'
batch2$batch <- 'batch2'
batch3$batch <- 'batch3'
merged <- merge(
  x = batch1,
  y = list(batch2, batch3))


#clean up and add sample IDs
merged$sample_id <- ''
merged$sample_id[grep("sum.23dpa4", merged$orig.ident)] <- "medium_bud_1"
merged$sample_id[grep("sum.23dpa5", merged$orig.ident)] <- "medium_bud_2"
merged$sample_id[grep("sum.S1", merged$orig.ident)] <- "medium_bud_3"
merged$sample_id[grep("sum.S2", merged$orig.ident)] <- "medium_bud_4"
merged$sample_id[grep("sum.S4", merged$orig.ident)] <- "medium_bud_5"
merged$sample_id[grep("sum.23dpa6", merged$orig.ident)] <- "medium_bud_6"
merged$sample_id[grep("sum.S3", merged$orig.ident)] <- "early_bud_1"
merged$sample_id[grep("sum.S5", merged$orig.ident)] <- "early_bud_2"
merged$sample_id[grep("sum.3dpa1", merged$orig.ident)] <- "wound_healing_1"
merged$sample_id[grep("sum.3dpa2", merged$orig.ident)] <- "wound_healing_2"
merged$sample_id[grep("sum.3dpa3", merged$orig.ident)] <- "wound_healing_3"
merged$sample_id[grep("sum.Contra1", merged$orig.ident)] <- "non-regenerating_1"
merged$sample_id[grep("sum.Contra2", merged$orig.ident)] <- "non-regenerating_2"
merged$sample_id[grep("sum.Contra3", merged$orig.ident)] <- "non-regenerating_3"
merged$sample_id[grep("sum.Intact1", merged$orig.ident)] <- "non-regenerating_4"
merged$sample_id[grep("sum.Intact2", merged$orig.ident)] <- "non-regenerating_5"

#saveRDS(merged, '20250626.cleaned.and.merged.before.PCA.rds')

merged <- JoinLayers(merged)
merged[["RNA"]] <- split(merged[["RNA"]], f = merged$batch)
merged
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged)
ElbowPlot(merged, ndims = 50)



merged <- FindNeighbors(merged, dims = 1:20, reduction = "pca")
merged <- FindClusters(merged, resolution = 1.5, cluster.name = "unintegrated_clusters")
merged <- RunUMAP(merged, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated_res1.5_20dims")

# visualize by batch and cell type annotation
DimPlot(merged, reduction = "umap.unintegrated_res1.5_20dims", group.by = c("sample_id", "unintegrated_clusters", 'batch'))

#inspect if batch correction is necessary, some admixing between batches but definitely needs correction
DimPlot(merged, reduction = "umap.unintegrated_res1.5_20dims", group.by = c('batch'), split.by = "batch")

#integration, tried various and CCA looked to be collapsing immune cells well
library(SeuratWrappers)


#takes a long time, save after and reload obj going forward
merged <- IntegrateLayers(
  object = merged, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

#saveRDS(merged,'20250814.merged.CCA.integrated.rds')
#merged <- readRDS('20250814.merged.CCA.integrated.rds')

merged <- FindNeighbors(merged, reduction = "integrated.cca", dims = 1:12)
merged <- FindClusters(merged, resolution = 0.8, cluster.name = "cca_clusters")
merged <- RunUMAP(merged, reduction = "integrated.cca", dims = 1:12, reduction.name = "umap.cca")

#assess integration
DimPlot(
  merged,
  reduction = "umap.cca",
  group.by = c("batch"), split.by = "batch",
  combine = FALSE, label.size = 2
)

#now refine clustering, this was done iteratively and these parameters were chosen
merged <- FindNeighbors(merged, dims = 1:12, k.param = 60, prune.SNN = 1 / 15, reduction =  "integrated.cca")
merged <- FindClusters(merged, graph.name = "RNA_snn", resolution = 0.5, algorithm = 1)
DimPlot(merged, reduction = "umap.cca", group.by = "RNA_snn_res.0.5", label=T) + ggtitle("louvain_0.5")

#now identify immune cell types
merged <- JoinLayers(merged)
merged@meta.data
all.markers <- FindAllMarkers(merged, only.pos = TRUE)
all.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


#automated labeling of cell types using CyteTypeR
prepped_data <- PrepareCyteTypeR(merged,
                                 all.markers,
                                 n_top_genes = 50,
                                 group_key = 'seurat_clusters',
                                 aggregate_metadata = TRUE,
                                 coordinates_key = "umap.cca")

metadata <- list(
  title = 'scRNA-seq analysis of axolotl limb regeneration, all primary cells',
  run_label = 'granular_annotation_all_cells',
  experiment_name = 'inDrops_axo',
  "DOI" = "https://doi.org/10.1038/s41467-018-07604-0",
  "GEO_Accession" = " GSE121737")

results_CyteTypeR <- CyteTypeR(obj=merged,
                              prepped_data = prepped_data, 
                              study_context = "immune and non immune cells from a time course of axolotl limb regeneration", 
                              metadata = metadata
)

#results are here: https://nygen-labs-prod--cytetype-api.modal.run/report/9681fb59-0ec1-4a96-84dd-7ad6025d0603#cluster-6


#upon inspection of results it doesnt seem the LOC annotation info was effecitivley used, so will manually add these
#after inspecting annotations
#LOC annotations are uninformative without any annotation, which clusters have a high proportion of LOC annotations in the top 50 genes
# Calculate % of LOC genes among top 40 markers per cluster
loc_summary <- all.markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%                 
  slice_head(n = 50) %>%                        
  summarise(
    total_genes = n(),
    loc_genes = sum(grepl("^LOC\\d+", gene)),   # count genes starting with LOC + numbers
    percent_loc = 100 * loc_genes / total_genes
  )

print(loc_summary)

#clusters 1, 3, 6, 7, 9, 12, all have over 50% of the top 50 genes as LOC genes. Let's add annotations
#and provide those to cytetype within the interactive html. We added annotations from the Gff file
gff_data <- read_tsv("/Users/nickleigh/Desktop/GCF_040938575.1_UKY_AmexF1_1_genomic.gff", comment = "#", col_names = FALSE)

# Filter for rows where column 3 is 'gene'
gene_rows <- gff_data %>% filter(X3 == "gene")

# Extract gene ID and description from column 9
gene_info <- gene_rows %>%
  mutate(
    gene_id = str_extract(X9, "gene=[^;]+"),
    gene_id = str_replace(gene_id, "gene=", ""),
    gene_annotation = str_extract(X9, "description=[^;]+"),
    gene_annotation = str_replace(gene_annotation, "description=", "")
  ) %>%
  dplyr::select(gene_id, gene_annotation)


# Join the annotation data to the original dataframe
annotated.markers <- all.markers %>%
  left_join(gene_info, by = c("gene" = "gene_id"))

#take first 50 of all clusters to improve annotation 

### "use these gene_annotations for the LOCs within the top 50 marker genes to refine this annotation""
annotated.markers %>%
  filter(cluster == 0) %>%
  slice_head(n = 50)
#1
annotated.markers %>%
  filter(cluster == 1) %>%
  slice_head(n = 50)
#2
annotated.markers %>%
  filter(cluster == 2) %>%
  slice_head(n = 50)
#3
annotated.markers %>%
  filter(cluster == 3) %>%
  slice_head(n = 50)
#4
annotated.markers %>%
  filter(cluster == 4) %>%
  slice_head(n = 50)
#5
annotated.markers %>%
  filter(cluster == 5) %>%
  slice_head(n = 50)
#6
annotated.markers %>%
  filter(cluster == 6) %>%
  slice_head(n = 50)
#7
annotated.markers %>%
  filter(cluster == 7) %>%
  slice_head(n = 50)

#8
annotated.markers %>%
  filter(cluster == 8) %>%
  slice_head(n = 50)
#9
annotated.markers %>%
  filter(cluster == 9) %>%
  slice_head(n = 50)
#10
annotated.markers %>%
  filter(cluster == 10) %>%
  slice_head(n = 50)
#11
annotated.markers %>%
  filter(cluster == 11) %>%
  slice_head(n = 50)
#12
annotated.markers %>%
  filter(cluster == 12) %>%
  slice_head(n = 50)
#13
annotated.markers %>%
  filter(cluster == 13) %>%
  slice_head(n = 50)
#14
annotated.markers %>%
  filter(cluster == 14) %>%
  slice_head(n = 50)
#14
annotated.markers %>%
  filter(cluster == 15) %>%
  slice_head(n = 50)

#in the html i copied in the top 50 rows with the annotation column into the 'reannotate' feature and said
#im adding annotation information to the LOC genes in the top 50 genes to further refine this annotations

#load in annotations after refinement

cytetype_annotations <- read.csv("/Users/nickleigh/Desktop/Lund/git/adaptive_immune_limb_regen/cytetype_annotations_all_cells.csv", stringsAsFactors = FALSE)

new_names <- cytetype_annotations$annotation
names(new_names) <- levels(merged)

merged <- RenameIdents(object = merged, new_names)

#now  color based on broad classsifications
clusters_broad <- Idents(merged)

# Define mapping of cluster names to classifications
classification_map <- list(
  # Immune
  "Macrophage" = "immune",
  "Antimicrobial Macrophage" = "immune",
  "CD8+ T Cell" = "immune",
  "B Cell" = "immune",
  
  # Erythroid
  "Erythroid Cell" = "erythroid",
  
  # Mesenchymal
  "Regenerative Blastema Fibroblast" = "mesenchymal",
  "Mesenchymal Progenitor Cell" = "mesenchymal",
  "Fibroblast" = "mesenchymal",
  "Schwann Cell Precursor" = "mesenchymal",
  
  # Epithelial
  "Epidermal Keratinocyte" = "epithelial",
  "Mucin-Secreting Epithelial Cell" = "epithelial",
  "Keratinocyte" = "epithelial",
  "Epithelial Cell" = "epithelial",
  "Basal Keratinocyte" = "epithelial",
  
  # Other
  "Vascular Endothelial Cell" = "endothelial" 
)

# Create a vector of classifications for each cell
classification_vector <- sapply(clusters_broad, function(x) {
  if (!is.null(classification_map[[x]])) {
    classification_map[[x]]
  } else {
    "unknown" # fallback for clusters not in the map
  }
})

# Add as metadata column
merged$classification <- classification_vector


#Figure 1A
DimPlot(
  merged,
  reduction = "umap.cca",
  group.by = 'classification',
  label = F,
  pt.size = 0.5
) +
  scale_color_viridis_d(option = "plasma") +  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    )

DimPlot(merged, label = T, reduction = 'umap.cca')

#####----------% immune at each time point-----------#########

#add time point info
merged$group <- NA
merged$group[grep("sum.23dpa4|sum.23dpa5|sum.23dpa6", merged$orig.ident)] <- "dpa23"
merged$group[grep("sum.3dpa1|sum.3dpa2|sum.3dpa3", merged$orig.ident)] <- "dpa3"
merged$group[grep("Contra|Intact", merged$orig.ident)] <- "limb"
merged$group[grep("sum.S1|sum.S2|sum.S4", merged$orig.ident)] <- "dpa23"
merged$group[grep("sum.S3|sum.S5", merged$orig.ident)] <- "dpa14"
merged$sample <- merged$orig.ident
VlnPlot(merged, features = 'PTPRC' )

md <- merged@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  as_tibble()

time_col    <- "group"       
rep_col     <- "sample_id"   
annot_col   <- "cytetype_annos" 

#define "immune"
immune_labels <- c(
  "Macrophage", "Antimicrobial Macrophage", "B Cell", "CD8+ T Cell")

#immune vs non-immune 
md <- md %>%
  mutate(is_immune = !!sym(annot_col) %in% immune_labels)

#c per replicate and time
rep_counts <- md %>%
  group_by(.data[[time_col]], .data[[rep_col]]) %>%
  summarise(
    immune_n = sum(is_immune, na.rm = TRUE),
    total_n  = n(),
    .groups = "drop"
  ) %>%
  mutate(pct_immune = 100 * immune_n / total_n)

rep_counts %>% arrange(.data[[time_col]], .data[[rep_col]]) %>% print(n = 50)


sumdf <- rep_counts %>%
  group_by(.data[[time_col]]) %>%
  summarise(
    mean_pct = mean(pct_immune),
    sd  = sd(pct_immune),
    n   = n(),
    se  = sd / sqrt(n),
    ci_low  = mean_pct - 1.96 * se,
    ci_high = mean_pct + 1.96 * se,
    .groups = "drop"
  )

#Figure 1B
ggplot() +
  geom_point(data = rep_counts,
             aes(x = .data[[time_col]], y = pct_immune, color = '#F08030'),
             position = position_jitter(width = 0.18, height = 0), alpha = 0.7, size = 4, show.legend = FALSE) +
  geom_point(data = sumdf, aes(x = .data[[time_col]], y = mean_pct), size = 3) +
  geom_errorbar(data = sumdf, aes(x = .data[[time_col]], ymin = ci_low, ymax = ci_high), width = 0.12) +
  scale_y_continuous(name = "% immune cells", limits = c(0, 100)) +
  scale_x_discrete(labels = c("limb" = "Limb", "dpa3" = "Wound\nHealing", "dpa14" = "Early\nBud", "dpa23" = "Medium\nBud")) +
  xlab("Regeneration stage") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y  = element_text(size = 14),   # y tick labels
    axis.title.y = element_text(size = 16)    # y axis title
  )



###-----subcluster immune cells--------####

#first perform less biased assessment of what clusters are immune cells
immune = merged[,Idents(merged) %in% c("Macrophage", "Antimicrobial Macrophage", "B Cell", "CD8+ T Cell")]
dim(immune)
#[1] 29653  8387

immune[["RNA"]] <- split(immune[["RNA"]], f = immune$batch)
immune = FindVariableFeatures(immune, verbose = FALSE)
immune = ScaleData(immune, assay = "RNA")
immune = RunPCA(immune, npcs = 50, verbose = F)
ElbowPlot(immune, ndims = 50)
immune <- IntegrateLayers(object = immune, 
                          method = CCAIntegration, orig.reduction = "pca", 
                          new.reduction = "integrated_immune", verbose = FALSE)

immune <- RunUMAP(immune, reduction = "integrated_immune", dims = 1:25, reduction.name = "umap_immune")

immune <- FindNeighbors(immune, reduction = "integrated_immune", dims = 1:25, k.param = 20, prune.SNN = 1 /15)

immune <- FindClusters(immune, graph.name = "RNA_snn", resolution = 0.75, algorithm = 1)

#inspect

DimPlot(immune, reduction = "umap_immune", group.by = c("RNA_snn_res.0.75", 'batch'), label=T) + ggtitle("batch")
DimPlot(immune, reduction = "umap_immune", split.by = 'batch', label=T) + ggtitle("batch") 

immune <- JoinLayers(immune)

Idents(immune) <- immune$RNA_snn_res.0.75
immune.markers <- FindAllMarkers(immune, only.pos = TRUE)


#now annotate
prepped_data_immune <- PrepareCyteTypeR(immune,
                                 immune.markers,
                                 n_top_genes = 50,
                                 group_key = 'seurat_clusters',
                                 aggregate_metadata = TRUE,
                                 coordinates_key = "umap_immune")

metadata <- list(
  title = 'My scRNA-seq analysis of axolotl limb regeneration',
  run_label = 'immune_analysis',
  experiment_name = 'inDrops_axo',
  "DOI" = "https://doi.org/10.1038/s41467-018-07604-0",
  "GEO_Accession" = " GSE121737")

results_immune <- CyteTypeR(obj=immune,
                              prepped_data = prepped_data_immune, 
                              study_context = "subsetted immune cells from axolotl limb regeneration time course, all primary cells", 
                              metadata = metadata
)
#results are here: https://nygen-labs-prod--cytetype-api.modal.run/report/4c77baa2-895a-4b9b-8c18-9c39288df9bf#cluster-21

View(results_immune@misc[["cytetype_results"]])

# Join the annotation data to the original dataframe
annotated.immune.markers <- immune.markers %>%
  left_join(gene_info, by = c("gene" = "gene_id"))

#take first 50 of all clusters to improve annotation 

### "use these gene_annotations for the LOCs within the top 50 marker genes to refine this annotation""
annotated.immune.markers %>%
  filter(cluster == 0) %>%
  slice_head(n = 50)
#1
annotated.immune.markers %>%
  filter(cluster == 1) %>%
  slice_head(n = 50)
#2
annotated.immune.markers %>%
  filter(cluster == 2) %>%
  slice_head(n = 50)
#3
annotated.immune.markers %>%
  filter(cluster == 3) %>%
  slice_head(n = 50)
#4
annotated.immune.markers %>%
  filter(cluster == 4) %>%
  slice_head(n = 50)
#5
annotated.immune.markers %>%
  filter(cluster == 5) %>%
  slice_head(n = 50)
#6
annotated.immune.markers %>%
  filter(cluster == 6) %>%
  slice_head(n = 50)
#7
annotated.immune.markers %>%
  filter(cluster == 7) %>%
  slice_head(n = 50)

#8
annotated.immune.markers %>%
  filter(cluster == 8) %>%
  slice_head(n = 50)
#9
annotated.immune.markers %>%
  filter(cluster == 9) %>%
  slice_head(n = 50)
#10
annotated.immune.markers %>%
  filter(cluster == 10) %>%
  slice_head(n = 50)
#11
annotated.immune.markers %>%
  filter(cluster == 11) %>%
  slice_head(n = 50)
#12
annotated.immune.markers %>%
  filter(cluster == 12) %>%
  slice_head(n = 50)
#13
annotated.immune.markers %>%
  filter(cluster == 13) %>%
  slice_head(n = 50)
#14
annotated.immune.markers %>%
  filter(cluster == 14) %>%
  slice_head(n = 50)
#15
annotated.immune.markers %>%
  filter(cluster == 15) %>%
  slice_head(n = 50)
#16
annotated.immune.markers %>%
  filter(cluster == 16) %>%
  slice_head(n = 50)
#17
annotated.immune.markers %>%
  filter(cluster == 17) %>%
  slice_head(n = 50)
#18
annotated.immune.markers %>%
  filter(cluster == 18) %>%
  slice_head(n = 50)
#19
annotated.immune.markers %>%
  filter(cluster == 19) %>%
  slice_head(n = 50)
#20
annotated.immune.markers %>%
  filter(cluster == 20) %>%
  slice_head(n = 50)


#in the html i copied in the top 50 rows with the annotation column into the 'reannotate' feature and said
#im adding annotation information to the LOC genes in the top 50 genes to further refine this annotations
#13 still seems wierd with the annotation so will provide top 100
#13
annotated.immune.markers %>%
  filter(cluster == 13) %>%
  slice_head(n = 100)
#didn't change anything
#see final info here: https://nygen-labs-prod--cytetype-api.modal.run/report/4c77baa2-895a-4b9b-8c18-9c39288df9bf#cluster-21

#manually changed one of the macrophage labels to phagocytic as this was a sub annotation and differentiates it from the other macrophage cluster
cytetype_annotations_immune <- read.csv("/Users/nickleigh/Desktop/Lund/git/adaptive_immune_limb_regen/cytetype_annotations_immune.csv", stringsAsFactors = FALSE)
new_names_immune <- cytetype_annotations_immune$annotation
names(new_names_immune) <- levels(immune)

immune <- RenameIdents(object = immune, new_names_immune)
#remove non immune cells from umap

immune_filtered <- subset(immune, idents = c('Keratinocyte', 'Fibroblast'), invert = TRUE)

#supplemental Figure 
annotated.immune.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

#heatmap will be better if LOC is replaced with informative gene_annotation informatoin
genes <- top10$gene

#anywith with LOC turns sto gene_annotation instead
labels <- ifelse(
  grepl("^LOC", top10$gene),       
  top10$gene_annotation,           
  top10$gene                       
)

# Name the vector for scale_y_discrete()
names(labels) <- genes

heatmap <- DoHeatmap(immune_filtered, features = genes) + NoLegend()

#now swap on LOCs where we can
heat <- heatmap + scale_y_discrete(labels = labels)

ggsave("../../../../../../../../../../Desktop/papers/2025/Rag1/SupplementalFigureheatmap.pdf", plot = heat, width = 35, height = 20, units = "in", limitsize = F) 


#Figure 1C
DimPlot(
  immune_filtered,
  reduction = "umap_immune",
  label = T,
  pt.size = 0.5
) +
  scale_color_viridis_d(option = "viridis") +  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
  ) 




#now get dynamics across regeneration
immune_md <- immune_filtered@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  as_tibble()

#set order of samples
immune_md[[time_col]] <- factor(immune_md[[time_col]], levels = c("limb", "dpa3", "dpa14", "dpa23"))

#total cell counts per time point and rep
totals <- immune_md %>%
  group_by(.data[[time_col]], .data[[rep_col]]) %>%
  summarise(total_n = n(), .groups = "drop")
colnames(immune_md)

#valuea of interest
time_col  <- "group"
rep_col   <- "sample_id"
annot_col <- "cytetype_annos"

#total per time poitn and rep
totals <- immune_md %>%
  group_by(.data[[time_col]], .data[[rep_col]]) %>%
  summarise(total_n = n(), .groups = "drop")

#cluster  per stage, replicate, and cluster
cluster_counts <- immune_md %>%
  group_by(.data[[time_col]], .data[[rep_col]], .data[[annot_col]]) %>%
  summarise(cluster_n = n(), .groups = "drop") %>%
  complete(
    !!sym(time_col),
    !!sym(rep_col),
    !!sym(annot_col),
    fill = list(cluster_n = 0)
  ) %>%
  left_join(totals, by = c(time_col, rep_col)) %>%
  mutate(
    pct_cluster = ifelse(total_n > 0, 100 * cluster_n / total_n, 0)
  )

#summarize across replicates
sumdf_cluster <- cluster_counts %>%
  filter(total_n > 0) %>%  # exclude fake combos
  group_by(.data[[time_col]], .data[[annot_col]]) %>%
  summarise(
    mean_pct = mean(pct_cluster, na.rm = TRUE),
    sd       = sd(pct_cluster, na.rm = TRUE),
    n        = dplyr::n(),  # now reflects actual replicates
    se       = sd / sqrt(n),
    ci_low   = mean_pct - 1.96 * se,
    ci_high  = mean_pct + 1.96 * se,
    .groups  = "drop"
  )


lvl <- levels(Idents(immune_filtered))

pal <- setNames(viridisLite::viridis(length(lvl), option = "D"), lvl)


# 1) Build a named palette once (names = cell types, values = hex colors)
#    You can generate from current levels or hand-pick colors.
library(viridisLite)

annot_col <- "cytetype_annos"
celltypes <- levels(factor(cluster_counts[[annot_col]]))
pal <- setNames(viridisLite::viridis(length(celltypes), option = "D"), celltypes)

# 2) Plot with color mapped to cell type and manual scale
dynamics <- ggplot() +
  # replicate-level points colored by cell type
  geom_point(
    data = cluster_counts,
    aes(x = .data[[time_col]], y = pct_cluster, color = .data[[annot_col]]),
    position = position_jitter(width = 0.18, height = 0),
    alpha = 0.7, size = 3
  ) +
  # mean per stage colored by cell type
  geom_point(
    data = sumdf_cluster,
    aes(x = .data[[time_col]], y = mean_pct, color = .data[[annot_col]]),
    size = 3
  ) +
  #add CIs
  geom_errorbar(
    data = sumdf_cluster,
    aes(x = .data[[time_col]], ymin = ci_low, ymax = ci_high),
    width = 0.12, color = "black"
  ) +
  scale_color_manual(values = pal) +                       
  scale_y_continuous(name = "% of Immune cells", limits = c(0, 100)) +
  scale_x_discrete(labels = c("limb"  = "Limb",
                              "dpa3"  = "Wound\nHealing",
                              "dpa14" = "Early\nBud",
                              "dpa23" = "Medium\nBud")) +
  xlab("Regeneration stage") +
  facet_wrap(vars(.data[[annot_col]]), nrow = 5) +
  theme_classic(base_size = 12) +
  theme(
    strip.text  = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14), 
    strip.background = ggplot2::element_rect(fill = "grey85", color = NA),
    legend.position = "none")
  


ggsave("../../../../../../../../../../Desktop/papers/2025/Rag1/SupplmentalFigure2.pdf", plot = dynamics, width = 10, height =10, units = "in") 
dev.off()



#saveRDS(immune_filtered, '20251114.immune_only.rds')

##----------GO Analysis of abundant cells-----------------##

#readRDS('20251114.immune_only.rds') -> immune_filtered

background_genes <- rownames(immune)[rowSums(GetAssayData(immune, assay = "RNA", layer = "counts")) > 0]
eggnog_mapping <- read.delim("/Users/nickleigh/Library/CloudStorage/OneDrive-LundUniversity/Dokument/Data/2025/axolotl/axolotl_eggnog_union_per_gene.tsv", header = T)

term2gene <- eggnog_mapping %>%
  filter(!is.na(GOs_union)) %>%
  separate_rows(GOs_union, sep = ",") %>%
  transmute(GO = str_trim(GOs_union), gene_id = gene_id) %>%
  distinct()


go_names <- AnnotationDbi::select(
  GO.db,
  keys = unique(term2gene$GO),
  columns = c("TERM", "DEFINITION", "ONTOLOGY"),
  keytype = "GOID"
)

term2name <- go_names %>%
  transmute(GO = GOID, Term = TERM) %>%
  distinct()



# Function to run DE (+ optional ORA) for a given cluster and timepoint
run_enrichment <- function(immune_obj, 
                           cluster_id, 
                           timepoint, 
                           sel.clust, 
                           term2gene, 
                           term2name, 
                           do_enrichment = TRUE) {
  # ---- Basic checks ----
  if (!sel.clust %in% colnames(immune_obj@meta.data)) {
    stop("Column '", sel.clust, "' not found in metadata.")
  }
  if (!"group" %in% colnames(immune_obj@meta.data)) {
    stop("Column 'group' not found in metadata. Please add or set identities accordingly.")
  }
  if (!"orig.ident" %in% colnames(immune_obj@meta.data)) {
    stop("Column 'orig.ident' not found in metadata.")
  }
  
  # ---- Subset the cluster ----
  this_cluster <- subset(
    immune_obj, 
    cells = colnames(immune_obj)[immune_obj@meta.data[[sel.clust]] == cluster_id]
  )
  this_cluster <- SetIdent(this_cluster, value = "group")
  
  # ---- Subset for limb vs timepoint ----
  cluster_sub <- subset(this_cluster, group %in% c("limb", timepoint))
  if (ncol(cluster_sub) == 0) {
    message("No cells found for cluster ", cluster_id, " at ", timepoint)
    return(list(edgeR_res = data.frame(), go_res = data.frame()))
  }
  
  # ---- Pseudobulk ----
  # Aggregate by sample (orig.ident), not by condition
  pseudobulk <- AggregateExpression(
    cluster_sub, 
    group.by = "orig.ident", 
    assays = "RNA"
  )$RNA
  
  # ---- Build sample labels -> map to conditions (limb vs timepoint) ----
  bulk.labels <- gsub("^sum\\.", "", colnames(pseudobulk))
  
  # Match your naming patterns; keep as close to your original rules as possible
  bulk.labels <- dplyr::case_when(
    grepl("Contra|Intact", bulk.labels, ignore.case = TRUE) ~ "limb",
    grepl("(^|_)3dpa",    bulk.labels, ignore.case = TRUE) ~ "dpa3",
    grepl("(^|_)23dpa",   bulk.labels, ignore.case = TRUE) ~ "dpa23",
    grepl("(^|_)S",       bulk.labels, ignore.case = TRUE) & timepoint == "dpa23" ~ "dpa23",
    TRUE ~ bulk.labels
  )
  
  bulk.labels <- factor(bulk.labels)
  
  # Keep only the two groups we care about (limb vs given timepoint)
  keep_samples <- bulk.labels %in% c("limb", timepoint)
  if (!any(keep_samples)) {
    message("No matching samples for cluster ", cluster_id, " at ", timepoint)
    return(list(edgeR_res = data.frame(), go_res = data.frame()))
  }
  pseudobulk <- pseudobulk[, keep_samples, drop = FALSE]
  bulk.labels <- droplevels(bulk.labels[keep_samples])
  
  # Ensure limb exists as reference
  if (!"limb" %in% levels(bulk.labels)) {
    message("Reference 'limb' not present after filtering for cluster ", cluster_id, " at ", timepoint)
    return(list(edgeR_res = data.frame(), go_res = data.frame()))
  }
  bulk.labels <- relevel(bulk.labels, ref = "limb")
  
  # ---- edgeR pipeline ----
  dge <- DGEList(counts = pseudobulk, group = bulk.labels)
  keep <- filterByExpr(dge, group = bulk.labels)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  if (nrow(dge) == 0) {
    message("No genes passed filtering for cluster ", cluster_id, " at ", timepoint)
    return(list(edgeR_res = data.frame(), go_res = data.frame()))
  }
  dge <- calcNormFactors(dge)
  
  design <- model.matrix(~ bulk.labels)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  
  # coef=2 corresponds to the non-reference level (timepoint vs limb)
  qlf <- glmQLFTest(fit, coef = 2)
  edgeR_tab <- topTags(qlf, n = Inf)$table
  edgeR_tab$gene <- rownames(edgeR_tab)
  
  # ---- If enrichment is not requested, return DE only ----
  if (!do_enrichment) {
    return(list(edgeR_res = edgeR_tab, go_res = data.frame()))
  }
  
  # ---- Prepare up/down gene lists for ORA ----
  up_genes <- edgeR_tab %>%
    dplyr::filter(logFC >= 1, FDR <= 0.05) %>%
    dplyr::pull(gene)
  
  down_genes <- edgeR_tab %>%
    dplyr::filter(logFC <= -1, FDR <= 0.05) %>%
    dplyr::pull(gene)
  
  gene_clusters <- list(Up = up_genes, Down = down_genes)
  
  universe_genes <- rownames(dge)  # background based on filtered genes
  
  # ---- Enrichment (ORA) ----
  # If both lists are empty, skip compareCluster
  if ((length(up_genes) == 0) && (length(down_genes) == 0)) {
    message("No DE genes pass thresholds; skipping enrichment for cluster ", cluster_id, " at ", timepoint)
    return(list(edgeR_res = edgeR_tab, go_res = data.frame()))
  }
  
  cc_res <- tryCatch({
    compareCluster(
      geneClusters = gene_clusters,
      fun = "enricher",
      universe = universe_genes,
      TERM2GENE = term2gene,
      TERM2NAME = term2name
    )
  }, error = function(e) {
    message("compareCluster failed: ", conditionMessage(e))
    return(NULL)
  })
  
  if (is.null(cc_res)) {
    return(list(edgeR_res = edgeR_tab, go_res = data.frame()))
  }
  
  cc_df <- cc_res@compareClusterResult
  if (is.null(cc_df) || nrow(cc_df) == 0) {
    message("No enriched terms for cluster ", cluster_id, " at ", timepoint)
    return(list(edgeR_res = edgeR_tab, go_res = data.frame()))
  }
  
  sig_terms <- cc_df %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::mutate(cluster = paste0("cluster", cluster_id),
                  timepoint = timepoint)
  
  return(list(edgeR_res = edgeR_tab, go_res = sig_terms))
}


clusters   <- c("T Cell", "B Cell", "Phagocytic Macrophage", "Neutrophil-Like Granulocyte", "Macrophage")
timepoints <- c("dpa3", "dpa23")
Idents(immune_filtered) -> immune_filtered$cytetype_annos
sel.clust <- "cytetype_annos"

edger_results <- list()

for (cl in clusters) {
  for (tp in timepoints) {
    message("Processing DE-only for ", cl, " at ", tp)
    res <- run_enrichment(
      immune_obj  = immune_filtered,
      cluster_id  = cl,
      timepoint   = tp,
      sel.clust   = sel.clust,
      term2gene   = term2gene,
      term2name   = term2name,
      do_enrichment = FALSE  # <--- key
    )
    edger_results[[paste(cl, tp, sep = "_")]] <- res$edgeR_res
  }
}

#annotate LOC genesin output

annotated_edgeR_results <- lapply(
  edger_results,
  function(df) {
    if (is.null(df) || nrow(df) == 0) return(df)
    dplyr::left_join(df, gene_info, by = c("gene" = "gene_id"))
  }
)


#dotplot of differentiallye expressed genes for each cell type
#filter for seurat obj with only cells we did DE on
immune_de <- subset(immune_filtered,  idents = c("T Cell", "B Cell", "Phagocytic Macrophage", "Neutrophil-Like Granulocyte", "Macrophage"))

table(immune_filtered$cond_rep)

immune_filtered$cond_rep <- interaction(immune_filtered$group,
                                        immune_filtered$sample_id,
                                  sep = "_")

#set order of samples for the graphs
#check out genes with biggest difference
samples <- unique(immune_filtered$sample_id)
samples
nonreg  <- samples[grepl("^non-regenerating", samples)]
wh      <- samples[grepl("^wound_healing", samples)]
eb      <- samples[grepl("^early_bud", samples)]
mb      <- samples[grepl("^medium_bud", samples)]
final_order <- c(nonreg, wh, eb, mb)
immune_filtered$sample_id <- factor(
  immune_filtered$sample_id,
  levels = final_order
)

cells_t <- WhichCells(immune_filtered, idents = "T Cell")

#downregulated first, upregulated second
#class II histocompatibility antigen%2C M beta 1 chain = LOC138522041
#ETS1 Th1 https://pmc.ncbi.nlm.nih.gov/articles/PMC2213045/
#STAT4 th1 https://pmc.ncbi.nlm.nih.gov/articles/PMC2768040/
#IFNGR1 th1
#'IL18RAP'
t_cell_genes <- c('GZMA', 'STAT4', 'ETS1','IFNGR1', 'CEBPD','IRF1' ,'LOC138522041', 'LGALS8')

t_cell_p <- dittoDotPlot(
  immune_filtered,
  t_cell_genes,
  group.by = "cond_rep",
  split.by = "cytetype_annos",
  cells.use = cells_t,
  assay = "RNA",
  min.color = "#800080",      # color-blind friendly blue for low expression
  mid.color = "grey90",       # grey for zero/midpoint
  max.color = "#FF8C00",      # purple for high expression
  legend.show = TRUE
) +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = NULL, color = "Avg.Exp", size = "% cells") +
  ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey85", color = NA))


t_cell_p


cells_b <- WhichCells(immune_filtered, idents = "B Cell")
#downregualted first
#LOC138522002 TLR13 recognizes bacterial RNA (this is down)
#ZNF593 https://www.nature.com/articles/s41418-025-01508-5 (suppresses cGAS) so should be more responsible 
#IFNGR1 importnat for GC formation https://rupress.org/jem/article/213/5/715/42108/IFN-receptor-and-STAT1-signaling-in-B-cells-are
#IKZF3 https://www.cell.com/immunity/fulltext/S1074-7613(00)80637-8?cc=y (plays role in B cell function, down)
#TRAF3IP3 (MZ b cell surival https://pubmed.ncbi.nlm.nih.gov/26011558/, up)
#BTG1 proliferative suppressor
#GALECTIN-1 Breg (https://pmc.ncbi.nlm.nih.gov/articles/PMC5807431/, down so less regulaotry)
#cd209 LOC138447708 b cell surival https://academic.oup.com/jimmunol/article/202/1_Supplement/123.5/7957367?login=true 
#LOC138580979 beta-2 micoglobulin, so apc is down? 
#narrative: seems to me tthat genes associated with B cell surival are up and a bit of a mixed bag with fxn

b_cell_genes <- c('IFNGR1','LGALS1','ZNF593', 'LOC138580979','ETS1', 'BTG1', 'TRAF3IP3', 'LOC138447708', 'LBH')
b_cell_p <- dittoDotPlot(
  immune_filtered,
  b_cell_genes,                               # features
    group.by = "cond_rep",
    split.by = "cytetype_annos",
    cells.use = cells_b,
    assay = "RNA",
    min.color = "#800080",      # color-blind friendly blue for low expression
    mid.color = "grey90",       # grey for zero/midpoint
    max.color = "#FF8C00",      # purple for high expression
    legend.show = TRUE
  ) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = NULL, color = "Avg.Exp", size = "% cells") +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey85", color = NA))
b_cell_p


table(immune_filtered$cytetype_annos, immune_filtered$group)
VlnPlot(immune_filtered, features = 'FOXP3')


#now get enriched GO terms
all_results <- list()

for (cl in clusters) {
  for (tp in timepoints) {
    message("Processing GO enrichment for ", cl, " at ", tp)
    res <- run_enrichment(
      immune_obj  = immune_filtered,
      cluster_id  = cl,
      timepoint   = tp,
      sel.clust   = sel.clust,
      term2gene   = term2gene,
      term2name   = term2name,
      do_enrichment = TRUE   # <--- key
    )
    # Only store GO results; skip if empty
    if (!is.null(res$go_res) && nrow(res$go_res) > 0) {
      all_results[[paste(cl, tp, sep = "_")]] <- res$go_res
    }
  }
}

# Combine results
combined_df <- bind_rows(all_results)

# ---- Plotting ----
combined_df <- combined_df %>%
  mutate(
    timepoint = factor(timepoint, levels = c("dpa3", "dpa23")),
    comparison = paste(cluster, timepoint, sep = "_"),
    significance = -log10(p.adjust)
  )


# Define mapping
cluster_labels <- c(
 "clusterT Cell" = "T Cell",
 "clusterB Cell" = "B Cell",
 "clusterPhagocytic Macrophage" = "Phagocytic Macrophage",
 "clusterMacrophage" = "Macrophage",
 'clusterNeutrophil-Like Granulocyte' = 'Neutrophil-Like Granulocyte')


# Apply mapping
combined_df <- combined_df %>%
 mutate(
   cluster = recode(cluster, !!!cluster_labels),
   cluster = factor(cluster, levels = cluster_labels)
 )


library(GO.db)
library(rrvgo)
library(GSEABase)

install.packages("../../../../Genome/axolotl/org.Amexicanum.eg.db/", repos = NULL, type = "source")
library(org.Amexicanum.eg.db)

simMatrix <- calculateSimMatrix(
  combined_df$ID,
  orgdb = "org.Amexicanum.eg.db",
  ont = "BP",
  keytype = "GID"
)

scores <- setNames(-log10(combined_df$qvalue), combined_df$ID)

reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Amexicanum.eg.db", 
                                keytype = 'GID')

#remove terms not in reduced terms
df_reduced <- combined_df %>%
  filter(ID %in% reducedTerms$go)


df_plot <- df_reduced %>%
  filter(qvalue < 0.005, RichFactor > 0.5, FoldEnrichment > 5) %>%
  left_join(
    reducedTerms %>%
      as.data.frame() %>%
      dplyr::select(term, cluster), 
    by = c("Description" = "term")
  )

#clean up some of rhe column names for better organization
df_plot <- df_plot %>%
  dplyr::rename(
    cluster   = any_of("cluster.y"),   # semantic similarity cluster
    Direction = any_of("Cluster"),     # up/down direction
    CellType  = any_of("cluster.x"), 
    go = any_of("ID")
  )

#filter strigent to get most enriched terms
#add semantic clustering information
df_plot <- df_plot %>%
  left_join(
    reducedTerms %>% dplyr::select(go, cluster),
    by = c("go", "cluster")
  ) 


df_plot$timepoint <- factor(df_plot$timepoint, levels = c("dpa3", "dpa23"))

#group based on semanitc clustering
#f_plot <- df_plot %>%
# arrange(cluster) %>%
# mutate(
#   cluster = factor(cluster, levels = unique(cluster)),
#   Description = factor(Description, levels = unique(Description))
# )

#give useful names to the clusters
cluster_labels <- reducedTerms %>%
  group_by(cluster) %>%
  summarize(cluster_label = unique(parentTerm))

df_plot$cluster <- as.integer(as.character(df_plot$cluster))


df_plot <- df_plot %>%
  left_join(cluster_labels, by = "cluster")


#group in a coherent fashion based on semantic similarity
custom_order <-c('maturation of SSU-rRNA','antibacterial humoral response', 'T cell activation', 'negative regulation of leukocyte mediated cytotoxicity'
                 ,'negative regulation of natural killer cell mediated cytotoxicity', 'regulation of natural killer cell mediated cytotoxicity','antigen processing and presentation of exogenous peptide antigen via MHC class I, TAP-dependent',
                  'negative regulation of mononuclear cell proliferation', 
                 'regulation of leukocyte cell-cell adhesion','production of molecular mediator of immune response',
                 'positive regulation of cytokine production','cellular response to type II interferon')

df_plot$cluster_label <- factor(df_plot$cluster_label, levels = custom_order)


df_plot <- df_plot %>%
  arrange(cluster_label, desc(RichFactor))

df_plot$Description <- factor(df_plot$Description, levels = unique(df_plot$Description))


#make names shorthand to condense figure


# Define abbreviation rules
rules <- list(
  "\\bpositive\\b" = "pos",
  "\\bnegative\\b" = "neg",
  "\\bimmunity\\b" = "imm",
  "\\bimmune\\b" = "imm",
  "\\bregulation\\b" = "reg",
  "natural\\s+killer" = "NK",
  "antigen" = "Ag",
  "presentation" = "pres",
  "interferon" = "IFN",
  "MHC\\s+class\\s+I" = "MHCI"
)

# Apply replacements
shorten <- function(text) {
  for (pattern in names(rules)) {
    text <- gsub(pattern, rules[[pattern]], text, ignore.case = TRUE)
  }
  text
}
short_names <- sapply(df_plot$Description, shorten)
df_plot$short_description <- shorten(df_plot$Description)

df_plot$Description <- factor(
  df_plot$short_description,
  levels = df_plot %>% 
    arrange(cluster_label, desc(RichFactor)) %>% 
    pull(short_description) %>% 
    unique() %>%
    rev()
)


#figure 1E
GO_plot <- ggplot(df_plot, aes(x = timepoint, y = Description, color = Direction, size = RichFactor)) +
  geom_point(alpha = 0.9, na.rm = TRUE) +
  # Keep legend entries only for Up/Down (no NA swatch)
  scale_color_manual(values = c("Up" = "#8B0000", "Down" = "#00008B"), na.translate = FALSE) +
  facet_wrap(~ CellType, ncol = 4, scales = "fixed") +  # B (left) | T (right)
  scale_y_discrete(drop = FALSE) +                      
  scale_x_discrete(labels = c("dpa3" = "Wound\nHealing", "dpa23" = "Medium\nBud")) +
  labs(
    x = "Timepoint",
    y = "Enriched Gene Sets",
    color = "Direction",         
    size  = "Rich Factor",
    title = "GO Term Enrichment for immune cells across regeneration") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 7),
    panel.grid.major.y = element_line(color = "gray90"),
    strip.background = element_rect(fill = "grey85", color = NA)
  )

GO_plot


# Save to PDF with same size
ggsave("../../../../../../../../../../Desktop/papers/2025/Rag1/Figure1E.pdf", plot = GO_plot, width = 15, height = 6, units = "in") 
dev.off()



ggplot(df_reduced, aes(x = timepoint, y = Description, color = Cluster, size = RichFactor)) +
  geom_point(alpha = 0.9, na.rm = TRUE) +
  # Keep legend entries only for Up/Down (no NA swatch)
  scale_color_manual(values = c("Up" = "#8B0000", "Down" = "#00008B"), na.translate = FALSE) +
  facet_wrap(~ cluster, ncol = 2, scales = "fixed") +  # B (left) | T (right)
  scale_y_discrete(drop = FALSE) +                      
  scale_x_discrete(labels = c("dpa3" = "Wound Healing", "dpa23" = "Medium Bud")) +
  labs(
    x = "Timepoint",
    y = "Enriched Gene Sets",
    color = "Direction",         
    size  = "Rich Factor",
    title = "GO Term Enrichment for B and T cells across regeneration") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 7),
    panel.grid.major.y = element_line(color = "gray90")
  )











