library(tidyverse)
library(readr)
library(stringr)

# --- Paths ---
eggnog_path <- "/Users/nickleigh/Library/CloudStorage/OneDrive-LundUniversity/Dokument/Data/2025/axolotl/axolotl_eggnog.emapper.annotations"
gff_path    <- "/Users/nickleigh/Desktop/GCF_040938575.1_UKY_AmexF1_1_genomic.gff"

# --- Read EggNOG emapper annotations (skip commented header lines) ---
eggnog <- read.delim(
  eggnog_path,
  header = TRUE,
  comment.char = "#",
  sep = "\t",
  quote = "",
  check.names = FALSE
)

stopifnot(all(c("query", "evalue", "score") %in% names(eggnog)))
if (!"GOs" %in% names(eggnog)) {
  warning("Column 'GOs' not found in EggNOG file; downstream GO union will be skipped.")
  eggnog$GOs <- NA_character_
}

# --- Read GFF ---
gff <- read_tsv(
  gff_path,
  comment = "#",
  col_names = FALSE,
  show_col_types = FALSE
)
colnames(gff)[1:9] <- c("seqid","source","type","start","end","score","strand","phase","attributes")

# Helper: extract key=value from GFF attributes
extract_attr <- function(attr, key) {
  m <- str_match(attr, paste0("(^|;)", key, "=([^;]+)"))
  ifelse(is.na(m[,3]), NA_character_, m[,3])
}

# Helper: normalize protein IDs (strip version numbers etc.)
#normalize_protein <- function(x) sub("\\.\\d+$", "", x)

# --- Gene table ---
genes <- gff %>%
  filter(type == "gene") %>%
  transmute(
    gene_id = coalesce(
      extract_attr(attributes, "gene"),
      extract_attr(attributes, "Name"),
      extract_attr(attributes, "ID")
    )
  ) %>%
  distinct()

# --- mRNA table ---
mrnas <- gff %>%
  filter(type %in% c("mRNA","transcript")) %>%
  transmute(
    transcript_id = coalesce(
      extract_attr(attributes, "ID"),
      extract_attr(attributes, "transcript_id")
    ),
    gene_parent = coalesce(
      extract_attr(attributes, "Parent"),
      extract_attr(attributes, "gene")
    )
  ) %>%
  filter(!is.na(transcript_id)) %>%
  distinct()

# --- CDS table ---
cds <- gff %>%
  filter(type == "CDS") %>%
  transmute(
    transcript_id = coalesce(
      extract_attr(attributes, "Parent"),
      extract_attr(attributes, "transcript_id")
    ),
    protein_id = coalesce(
      extract_attr(attributes, "protein_id"),
      extract_attr(attributes, "proteinId"),
      extract_attr(attributes, "ID")
    )
  ) %>%
  filter(!is.na(protein_id)) %>%
  distinct()

# --- Protein -> gene mapping ---
prot_tx <- cds %>%
  inner_join(mrnas, by = "transcript_id") %>%
  transmute(
    protein_id,
    gene_id = gene_parent
  ) %>%
  distinct() %>%
  mutate(
    protein_id_norm = normalize_protein(protein_id),
    gene_id = sub("^gene-", "", gene_id)   # <--- drop 'gene-' prefix
  )

# --- Join EggNOG with mapping ---
eggnog_joined <- eggnog %>%
  mutate(protein_id = normalize_protein(query)) %>%
  left_join(prot_tx, by = c("protein_id" = "protein_id_norm"))

match_rate <- mean(!is.na(eggnog_joined$gene_id))
message(sprintf("Initial join match rate (query -> protein_id): %.1f%%", 100 * match_rate))

n_missing <- sum(is.na(eggnog_joined$gene_id))
if (n_missing > 0) {
  message(sprintf("Warning: %d EggNOG rows could not be mapped to a gene via GFF.", n_missing))
}

#----Union GO terms per gene ---
union_na_safe <- function(x) {
  x <- x[!is.na(x) & x != "-" & x != ""]
  if (length(x) == 0) return(NA_character_)
  u <- unique(unlist(strsplit(x, ",")))
  paste(u, collapse = ",")
}

eggnog_union_per_gene <- eggnog_joined %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id) %>%
  summarise(
    n_isoforms = n(),
    .best_row = which.max(suppressWarnings(as.numeric(score))),
    query_best = query[.best_row],
    score_best = suppressWarnings(max(as.numeric(score), na.rm = TRUE)),
    evalue_best = suppressWarnings(min(as.numeric(evalue), na.rm = TRUE)),
    GOs_union = union_na_safe(GOs)
  ) %>%
  ungroup()

# --- Summary ---
message(sprintf("EggNOG rows: %d", nrow(eggnog)))
message(sprintf("Mapped to genes: %d", sum(!is.na(eggnog_joined$gene_id))))
message(sprintf("Distinct genes (mapped): %d", n_distinct(eggnog_joined$gene_id[!is.na(eggnog_joined$gene_id)])))
message(sprintf("Union-per-gene rows: %d", nrow(eggnog_union_per_gene)))

save_path <- "/Users/nickleigh/Library/CloudStorage/OneDrive-LundUniversity/Dokument/Data/2025/axolotl/axolotl_eggnog_union_per_gene.tsv"

write_tsv(eggnog_union_per_gene, file = save_path)

