# Install Bioconductor packages if needed 
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DECIPHER", quietly = TRUE)) BiocManager::install("DECIPHER")
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")

# Install CRAN packages if missing
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")

# Load libraries
library(DECIPHER)
library(Biostrings)
library(ggplot2)
library(reshape2)

# Detect gap characters used by DECIPHER
is_gap <- function(x) {
  x %in% c("-", ".", "–")  # hyphen, dot, en dash
}

# -----------------------
# Example sequences: replace with your own
# -----------------------
seqs <- DNAStringSet(c(
  Clone1 = "CAAGATCCTATAGGGCAAGCAGTGGGTATCAACGCAGAGTACATGGGGACCACAGTCACAATACACACGAGGGGCAGCGGCTGTGAGCCTATCTCAGCAAAGGGAAGTAGGATGGCGCAAATCATCGTGTTTGTGGTTTTCATGGTCTTTCCTGCCTGGGTTCTCTCTGAGGTGTCCCTGGTGGAGTCCGGTGTAGGGATGCTGAAGCCCTCGGAGAGCCTGAGGCTGAGCTGTAAAGTTACTGGACACACCATCACTACAAGCCACGACTGGCACTGGATCCGCCAACCTCACGGGAAAGCCCTGGAATGGATCGGGCTAATAAGCTGGAATACAGTTACCGTGTATTACGCACCTTCTCTTCAAAGCAGATCAAAAATCACGATGGACCCCTCCAAAAATGAGTATTCTCTACAGATGACAGACATGAGGGCTGAAGACGAGGGCAAGTACTTCTGTGCAAGAGGCACAGCCTGGGAGAGGAAGGAGGAGCTGCAACAAGAATCTGCTCGGAGAGCCGGGGCAGCCTACTTCGATTACTGGGGACAAGGGACCTTTGTCACCGTCACCACAGAGCATATTTCCAAGCCTACGGTTTATCCGTTGGTAAATTGTTGTGGGAGCCAGTCTGGCACCAATGTGGCAATCGGATGCCTAGTGACTGACTATCTCCCTGAGCCGGTGGAGGTGACTTGGAGCAACGGCAACCCTATCTCAAATGAGGTGAAAATCTACTCTTCATTTGTGCTGAAAAGCACTCAAAAATACATCCGAAGCAGCCACCTGTATGTCCCAACATCAGACTGGGAGAGCAAGAAGTACACCTGTACGGTGAAACACCAATCCACAAGCCCTGAAGTCAGCCATACCGTGTCAAGAAGCTTGGCGTAATCATG",
  Clone2 = "CCAGATCCTATAGGGGCAAGCAGTGGGTATCAACGCAGAGTACATGGGGACCACAGATCACACAGCACAGGAGGGGCAGCTGTGAATCTGTGCTGTTAACTTCAGGGCACTATGGCGCATCTCATTGTGTTTATCGTCTTCTTGGTCTTTCCTTCCTACGTCCATCCGCAGATCACGCTGCTGGAATCTGGTCCAGGGATGCTGAAGCCAGCGGAGACCTTCACGCTGAGCTGTAAAGTTACCGGATTCACAATCAGTAGCAATTATTGGTGGGCTTGGATCCGTCAGCTCCCTGGAAAGGGCCTAGAGTGGATAGGGCAAATCGAACATAATGATGCTAGAAAATTTTATTCTGATTTGCTGAAAAGTCGAACGTCTATTACGAAGGACAACGCCAGGAATGAGTATTACCTGCAACTAAAGGATGCGAGAGCCGAAGATTCGGGCGTGTATTACTGCGGGAGAGATGGGAGGGGAGGTGGAGCTTTTGATATTTGGGGGCCAGGGACCATGGTGACTGTCACTTCAGAGCATATTTCCAAGCCTACGGTTTATCCGTTGGTAAATTGTTGTGGGAGCCAGTCTGGCACCAATGTGGCAATCGGATGCCTAGTGACTGACTATCTCCCTGAGCCGGTGGAGGTGACTTGGAGCAACGGCAACCCTATCTCAAATGAGGTGAAAATCTACTCTTCATTTGTGCTGAAAAGCACTCAAAAATACATCCGAAGCAGCCACCTGTATGTCCCAACATCAGACTGGGAGAGCAAGAAGTACACCTGTACGGTGAAACACCAATCCACAAGCCCTGAAGTCAGCCATACCGTGTCAAGAAGCTTGGCGTAATCATGGTCA",
  Clone3 = "CATGATCTATAAGGGCAAGCAGTGGGTATCAACGCAGAGTACATGGGGACCACAGATCACACAGCACAGGAGGGGCAGCTGTGAATCTGTGCTGTTAACTTCAGGGCACTATGGCGCATCTCATTGTGTTTATCGTCTTCTTGGTCTTTCCTTCCTACGTCCATCCGCAGATCACGCTGCTGGAATCTGGTCCAGGGATGCTGAAGCCAGCGGAGACCTTCACGCTGAGCTGTAAAGTTACCGGATTCACAATCAGTAGCAATTATTGGTGGGCTTGGATCCGTCAGCTCCCTGGAAAGGGCCTAGAGTGGATAGGGCAAATCGAACATAATGATGCTAGAAAATTTTATTCTGATTTGCTGAAAAGTCGAACGTCTATTACGAAGGACAACGCCAGGAATGAGTATTACCTGCAACTAAAGGATGCGAGAGCCGAAGATTCGGGCGTGTATTACTGCGGGAGAGATGGGGGCAGCATAGGTAGCAGCTTCGATGTGTGGGGGAAAGGAGTCCAGGTCACCGTTACTTCAGAGCATATTTCCAAGCCTACGGTTTATCCGTTGGTAAATTGTTGTGGGAGCCAGTCTGGCACCAATGTGGCAATCGGATGCCTAGTGACTGACTATCTCCCTGAGCCGGTGGAGGTGACTTGGAGCAACGGCAACCCTATCTCAAATGAGGTGAAAATCTACTCTTCATTTGTGCTGAAAAGCACTCAAAAATACATCCGAAGCAGCCACCTGTATGTCCCAACATCAGACTGGGAGAGCAAGAAGTACACCTGTACGGTGAAACACCAATCCACAAGCCCTGAAGTCAGCCATACCGTGTCAAGAAGCTTGGCGTAATCATG",
  Clone4 = "AATGCTCCTATAGGGCAAGCAGTGGGTATCAACGCAGTAGTACATGGGACAGGAGGCGGATTATTCTTACCGCACAGAAGCAGGCAGAATGTTATCCACTGGTTTTAAGAGAAATATGTCCCCGTCTTTACAAATCTCACTTTTGCTGACTGTATTGTCCTGTGTCCAGTCACAGATAACTCTGACCCAATCGGGGTCAGAAATCAGGAAACCAGGGGAGTCTGTGAAGCTGAAGTGTCTTGTGAGTGGATTTAACATAAACAGCTACTGGATGAACTGGATTAGGCAGGCTCCGGGGCGAGGTCTGGAATGGGTAGCGAGATATAATTCTGGGAGCAGCCCCCCACATTATTCCTCCGATGCAGTAAAGGGCAGATTCACCGCCTCCACAGACAGTTCCTCTCTCTATCTGCAAATGAACAACCTGAAGACGGAAGACACTGGCGTTTATTACTGTGCACGGGACCACTACAGCAAGAACAATGCTTTCTTTGATTACTGGGGTCAGGGGACCCTAGTCACAGTCACACAAGAGCATATTTCCAAGCCTACGGTTTATCCGTTGGTAAATTGTTGTGGGAGCCAGTCTGGCACCAATGTGGCAATCGGATGCCTAGTGACTGACTATCTCCCTGAGCCGGTGGAGGTGACTTGGAGCAACGGCAACCCTATCTCAGATGAGGTGAAAATCTACTCTTCATTTGTGCTGAAAAGCACTCAAAAATACATCCGAAGCAGCCACCTGTATGTCCCAACATCAGACTGGGAGAGCAAGAAGTACACCTGTACGGTGAAACACCAATCCACAAGCCCTGAAGTCAGCCATACCGTGTCAAGAAGCTTGGCGTAATCATG"
))

# -----------------------
# Align sequences
# -----------------------
aligned <- AlignSeqs(seqs, verbose = FALSE)

# Convert to character matrix
aln_char <- as.character(aligned)
seq_names <- names(aligned)
aln_list <- strsplit(aln_char, split = "")
aln_mat <- do.call(rbind, aln_list)
rownames(aln_mat) <- seq_names

n_seqs <- nrow(aln_mat)
n_cols <- ncol(aln_mat)

# -----------------------
# Column conservation accounting for gaps
# -----------------------
cat_level <- character(n_cols)
for (j in seq_len(n_cols)) {
  column <- aln_mat[, j]
  gap_mask <- is_gap(column)
  
  if (all(gap_mask)) {
    cat_level[j] <- "gap"
  } else if (any(gap_mask)) {
    # Any gap in the column -> not conserved
    cat_level[j] <- "none"
  } else {
    tab <- table(column)
    top_count <- max(tab)
    if (top_count == n_seqs) {
      cat_level[j] <- "identical"
    } else if (top_count / n_seqs >= 0.5) {
      cat_level[j] <- "partial"
    } else {
      cat_level[j] <- "none"
    }
  }
}

# -----------------------
# Prepare data frame for plotting
# -----------------------
df <- as.data.frame(aln_mat, stringsAsFactors = FALSE)
df$Sequence <- rownames(aln_mat)
df_m <- melt(df, id.vars = "Sequence", variable.name = "V", value.name = "Base")
df_m$Position <- as.integer(sub("^V", "", df_m$V))

# Split alignment into 3 paragraphs
# Split alignment into 4 paragraphs
n_panels <- 6
L <- max(df_m$Position, na.rm = TRUE)
chunk_size <- ceiling(L / n_panels)
df_m$Paragraph <- ceiling(df_m$Position / chunk_size)
df_m$Paragraph <- factor(df_m$Paragraph, levels = 1:n_panels, labels = paste("Paragraph", 1:n_panels))

# Assign column-level conservation
df_m$Conservation <- cat_level[df_m$Position]

# Replace all gap-like characters with "-"
df_m$Base[is_gap(df_m$Base)] <- "-"

# Individual gaps override conservation
df_m$ConservationPlot <- df_m$Conservation
df_m$ConservationPlot[is_gap(df_m$Base)] <- "gap"
df_m$ConservationPlot <- factor(df_m$ConservationPlot, levels = c("identical","partial","none","gap"))

# Top-down sequence order
df_m$Sequence <- factor(df_m$Sequence, levels = rev(seq_names))

# -----------------------
# Color palette
# -----------------------
palette_vals <- c(
  "identical" = "goldenrod3",
  "partial"   = "gray",
  "none"      = "steelblue3",
  "gap"       = "white"
)

# -----------------------
# Plot
# -----------------------
p <- ggplot(df_m, aes(x = Position, y = Sequence, fill = ConservationPlot)) +
  geom_tile(color = "black", linewidth = 0.2) +
  geom_text(aes(label = Base), size = 3) +
  scale_fill_manual(values = palette_vals, na.value = "white") +
  facet_wrap(~ Paragraph, ncol = 1, scales = "free_x") +
  theme_minimal() +
  labs(title = "Multiple Sequence Alignment — Conservation",
       x = "Alignment position",
       y = "",
       fill = "Conservation") +
  theme(
    text = element_text(family = "Arial"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )

print(p)

# -----------------------
# Summary stats
# -----------------------
tbl <- table(cat_level)
cat("Conservation summary (columns):\n")
print(tbl)

