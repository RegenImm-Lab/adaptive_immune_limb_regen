library(ontologyIndex)
library(dplyr)
library(stringr)
library(dplyr)
library(ggplot2)
library(forcats)
library(clusterProfiler)
library(edgeR)
library(tidytext)

readRDS('20250811.immune_only.rds') -> immune

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
  columns = c("TERM"),
  keytype = "GOID"
)

term2name <- go_names %>%
  transmute(GO = GOID, Term = TERM) %>%
  distinct()



# Function to run DE + ORA for a given cluster and timepoint
run_enrichment <- function(immune_obj, cluster_id, timepoint, sel.clust, term2gene, term2name) {
  # Subset the cluster
  this_cluster <- subset(immune_obj, cells = colnames(immune_obj)[immune_obj@meta.data[, sel.clust] == cluster_id])
  this_cluster <- SetIdent(this_cluster, value = "group")
  
  # Subset for limb vs timepoint
  cluster_sub <- subset(this_cluster, group %in% c("limb", timepoint))
  
  # Pseudobulk
  pseudobulk <- AggregateExpression(cluster_sub, group.by = "orig.ident", assays = "RNA")$RNA
  # Define minimum number of cells required in each condition
  min_cells <- 10
  
  # Clean sample labels
  bulk.labels <- gsub("^sum\\.", "", colnames(pseudobulk))
  
  # Map sample names to experimental groups based on patterns
  bulk.labels <- case_when(
    grepl("Contra|Intact", bulk.labels, ignore.case = TRUE) ~ "limb",
    grepl("(^|_)3dpa", bulk.labels, ignore.case = TRUE) ~ "dpa3",
    grepl("(^|_)23dpa", bulk.labels, ignore.case = TRUE) ~ "dpa23",
    grepl("(^|_)S", bulk.labels, ignore.case = TRUE) & timepoint == "dpa23" ~ "dpa23",
    TRUE ~ bulk.labels
  )
  
  bulk.labels <- as.factor(bulk.labels)
  print(table(bulk.labels))
  
  # Relevel to make limb the reference
  bulk.labels <- relevel(bulk.labels, ref = "limb")
  
  # edgeR pipeline
  dge <- DGEList(counts = pseudobulk, group = factor(bulk.labels))
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  design <- model.matrix(~ bulk.labels)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = 2)
  res <- topTags(qlf, n = Inf)$table
  res$gene <- rownames(res)
  
  # Get up/downregulated genes
  up_genes <- res %>% filter(logFC >= 1 & FDR <= 0.05) %>% pull(gene)
  down_genes <- res %>% filter(logFC <= -1 & FDR <= 0.05) %>% pull(gene)
  
  gene_clusters <- list(Up = up_genes, Down = down_genes)
  
  # Enrichment
  cc_res <- compareCluster(
    geneClusters = gene_clusters,
    fun = "enricher",
    TERM2GENE = term2gene,
    TERM2NAME = term2name
  )
  
  sig_terms <- cc_res@compareClusterResult %>%
    filter(p.adjust < 0.05) %>%
    mutate(cluster = paste0("cluster", cluster_id),
           timepoint = timepoint)
  
  return(sig_terms)
}

# ---- Run for multiple clusters and timepoints ----
clusters <- c("0", "1", "2", "7")
timepoints <- c("dpa3", "dpa23")
sel.clust <- "RNA_snn_res.0.75"

all_results <- list()
for (cl in clusters) {
  for (tp in timepoints) {
    message("Processing cluster ", cl, " at ", tp)
    res <- run_enrichment(immune, cl, tp, sel.clust, term2gene, term2name)
    if (nrow(res) > 0) {
      all_results[[paste(cl, tp, sep = "_")]] <- res
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
  "cluster0" = "Myeloid enriched WH",
  "cluster1" = "T cells",
  "cluster2" = "Myeloid enriched blastema",
  "cluster7" = "B cells"
)

# Apply mapping
combined_df <- combined_df %>%
  mutate(
    cluster = recode(cluster, !!!cluster_labels),
    # Keep the same order but now with new names
    cluster = factor(cluster, levels = cluster_labels)
  )


#get just t and be cells
tb_df <- combined_df %>%
  mutate(
    cluster = recode(as.character(cluster), !!!cluster_labels, .default = as.character(cluster)),
    timepoint = factor(timepoint, levels = c("dpa3", "dpa23")),
    significance = -log10(p.adjust)
  ) %>%
  filter(cluster %in% c("T cells", "B cells"))


top50_bp_df <- tb_df %>%
  group_by(comparison, timepoint, Cluster) %>%
  arrange(desc(RichFactor)) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 50) %>%
  mutate(Description = forcats::fct_reorder(Description, RichFactor)) %>%
  ungroup()

keep_terms <- c(
  "regulation of memory T cell activation",
  "positive regulation of memory T cell activation",
  "antigen processing and presentation of exogenous peptide antigen via MHC class I, TAP-independent",
  "antigen processing and presentation of endogenous peptide antigen via MHC class I via ER pathway, TAP-independent",
  "immune response-inhibiting signal transduction",
  "immune response-inhibiting cell surface receptor signaling pathway",
  "regulation of TRAIL production",
  "positive regulation of TRAIL production",
  "antigen processing and presentation of peptide antigen via MHC class Ib",
  "positive regulation of CD8-positive, alpha-beta T cell proliferation",
  "antigen processing and presentation via MHC class Ib",
  "negative regulation of dendritic cell differentiation",
  "antigen processing and presentation of endogenous peptide antigen via MHC class I via ER pathway, TAP-dependent",
  "antigen processing and presentation of endogenous peptide antigen via MHC class Ib via ER pathway",
  "antigen processing and presentation of endogenous peptide antigen via MHC class Ib via ER pathway, TAP-dependent",
  "interleukin-17 production",
  "antigen processing and presentation of endogenous peptide antigen via MHC class Ib",
  "antigen processing and presentation of endogenous peptide antigen via MHC class I via ER pathway",
  "positive regulation of T cell tolerance induction",
  "regulation of T cell anergy",
  "regulation of lymphocyte anergy",
  "regulation of CD8-positive, alpha-beta T cell proliferation",
  "antigen processing and presentation of endogenous peptide antigen",
  "antigen processing and presentation of endogenous peptide antigen via MHC class I",
  "regulation of dendritic cell differentiation",
  "positive regulation of tolerance induction",
  "regulation of T cell tolerance induction",
  "antigen processing and presentation of endogenous antigen",
  "positive regulation of CD8-positive, alpha-beta T cell activation",
  "positive regulation of T cell mediated cytotoxicity",
  "cellular response to interleukin-2",
  "positive regulation of regulatory T cell differentiation",
  "response to interleukin-2",
  "antigen processing and presentation of exogenous peptide antigen via MHC class I, TAP-dependent",
  "CD8-positive, alpha-beta T cell activation",
  "regulation of CD8-positive, alpha-beta T cell activation",
  "positive regulation of T-helper 1 cell differentiation",
  "regulation of T cell mediated cytotoxicity",
  "regulation of tolerance induction",
  "positive regulation of T cell cytokine production",
  "positive regulation of natural killer cell mediated immunity",
  "antigen processing and presentation of exogenous peptide antigen via MHC class I",
  "positive regulation of interleukin-13 production",
  "cytokine production involved in immune response",
  "positive regulation of leukocyte mediated cytotoxicity",
  "DNA ligation involved in DNA repair",
  "positive regulation of interferon-beta production",
  "monocyte chemotaxis",
  "mononuclear cell migration",
  "regulation of stem cell division",
  "positive regulation of stem cell population maintenance",
  "V(D)J recombination",
  "mammary gland epithelial cell differentiation",
  "apoptotic DNA fragmentation",
  "positive regulation of metaphase/anaphase transition of cell cycle",
  "positive regulation of chromosome separation",
  "positive regulation of meiotic nuclear division",
  "cellular component disassembly involved in execution phase of apoptosis",
  "inflammatory response to antigenic stimulus",
  "positive regulation of interferon-beta production",
  "negative regulation of stem cell differentiation",
  "somatic diversification of immune receptors via germline recombination within a single locus",
  "negative regulation of extrinsic apoptotic signaling pathway via death domain receptors",
  "regulation of stem cell population maintenance",
  "somatic cell DNA recombination",
  "somatic diversification of immune receptors"
)


keep_terms <- keep_terms[!duplicated(keep_terms)]


df_filtered <- top50_bp_df %>%
  arrange(RichFactor) %>%
  group_by(comparison, timepoint, Cluster) %>%
  filter(Description %in% keep_terms) %>%
  mutate(
    Description = factor(Description, levels = keep_terms), 
    Description = forcats::fct_drop(Description)             
  )


# Define manual groups for related terms
go_groups <- list(
  "Antigen Presentation" = c("antigen processing and presentation of endogenous peptide antigen",
                             "antigen processing and presentation via MHC class Ib",
                             "antigen processing and presentation of exogenous peptide antigen via MHC class I, TAP-dependent",
                             "antigen processing and presentation of endogenous peptide antigen via MHC class I",
                             "antigen processing and presentation of exogenous peptide antigen via MHC class I, TAP-independent",
                             "antigen processing and presentation of endogenous peptide antigen via MHC class I via ER pathway, TAP-dependent",
                             "antigen processing and presentation of endogenous peptide antigen via MHC class Ib via ER pathway, TAP-dependent"),
  
  "T Cell Regulation" = c("positive regulation of T cell mediated cytotoxicity",
                          "regulation of T cell tolerance induction",
                          "positive regulation of T cell tolerance induction",
                          "regulation of T cell anergy",
                          "regulation of CD8-positive, alpha-beta T cell activation",
                          "positive regulation of CD8-positive, alpha-beta T cell proliferation"),
  
  "Stem Cell Processes" = c("regulation of stem cell division",
                            "positive regulation of stem cell population maintenance",
                            "negative regulation of stem cell differentiation"),
  
  "DNA Repair & Apoptosis" = c("DNA ligation involved in DNA repair",
                               "cellular component disassembly involved in execution phase of apoptosis",
                               "apoptotic DNA fragmentation"),
  
  "Cytokine & Interleukin Response" = c("interleukin-17 production",
                                        "response to interleukin-2",
                                        "positive regulation of interleukin-13 production",
                                        "positive regulation of interferon-beta production"),
  
  "Other Immune" = c("immune response-inhibiting cell surface receptor signaling pathway",
                     "positive regulation of natural killer cell mediated immunity",
                     "monocyte chemotaxis",
                     "mononuclear cell migration")
)

# Make a lookup table
term_to_group <- enframe(go_groups, name = "Group", value = "Terms") %>%
  unnest_longer(Terms)

# Add to your dataframe
df_grouped <- df_filtered %>%
  left_join(term_to_group, by = c("Description" = "Terms")) %>%
  mutate(Group = ifelse(is.na(Group), "Other / Unclassified", Group))

df_grouped <- df_grouped %>%
  arrange(Group, Description) %>%
  mutate(Description = fct_inorder(Description))

ggplot(df_grouped, aes(x = timepoint, y = Description, color = Cluster, size = RichFactor)) +
  geom_point(alpha = 0.9, na.rm = TRUE) +
  # Keep legend entries only for Up/Down (no NA swatch)
  scale_color_manual(values = c("Up" = "#8B0000", "Down" = "#00008B"), na.translate = FALSE) +
  facet_wrap(~ cluster, ncol = 2, scales = "fixed") +  # B (left) | T (right)
  scale_y_discrete(drop = FALSE) +                      # show ALL rows in both panels
  labs(
    x = "Timepoint",
    y = "Enriched Gene Sets",
    color = "Direction",         # legend title for Up/Down
    size  = "Rich Factor",
    title = "GO Term Enrichment for B and T Cells",
    subtitle = "Selected GO Terms"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 7),
    panel.grid.major.y = element_line(color = "gray90")
  )














