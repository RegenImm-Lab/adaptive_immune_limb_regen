#make annotation db for axolotl with eggnog annotations
library(AnnotationForge)
library(dplyr)
library(tidyr)

#make "GID","GO","EVIDENCE"
eggnog_mapping <- read.delim("/Users/nickleigh/Library/CloudStorage/OneDrive-LundUniversity/Dokument/Data/2025/axolotl/axolotl_eggnog_union_per_gene.tsv", header = T)

nGO <- eggnog_mapping %>%
  filter(!is.na(GOs_union)) %>%
  separate_rows(GOs_union, sep = ",") %>%
  transmute(
    GO        = str_trim(GOs_union),
    GID = query_best,
  ) %>%
  distinct()

#format how the package wants
nGO <- nGO %>%
  transmute(
    GID = GID,          
    GO = GO,           
    EVIDENCE = "IEA"    
  )


#make "GID","SYMBOL","GENENAME"
gff <- read.delim("/Users/nickleigh/Desktop/GCF_040938575.1_UKY_AmexF1_1_genomic.gff",
                  comment.char = "#",
                  header = FALSE,
                  sep = "\t",
                  quote = "",
                  fill = TRUE,
                  stringsAsFactors = FALSE)
colnames(gff) <- c("seqid","source","type","start","end","score","strand","phase","attributes")

prot_id <- str_extract(gff$attributes, "protein_id=[^;]+") |>
  str_replace("protein_id=", "")
gene_symbol <- str_extract(gff$attributes, "gene=[^;]+") |>
  str_replace("gene=", "")
gene_name <- str_extract(gff$attributes, "product=[^;]+") |>
  str_replace("product=", "")

# make "GID","SYMBOL","GENENAME"
nSym <- data.frame(
  GID = prot_id,
  SYMBOL = gene_symbol,
  GENENAME = gene_name,
  stringsAsFactors = FALSE
)

#now GID and CHR
nChr <- gff %>%
  transmute(
    GID = str_extract(gff$attributes, "protein_id=[^;]+") |>
      str_replace("protein_id=", ""),
    CHROMOSOME = seqid,         
    
  )

#prepare for makeorgpacckage
#remove useless rows and duplicates

nSym <- nSym %>%
  filter(!is.na(GID)) %>%                   
  distinct(GID, .keep_all = TRUE)  

#now nChr
nChr <- nChr %>%
  filter(!is.na(GID)) %>%                   
  distinct(GID, .keep_all = TRUE) 


makeOrgPackage(gene_info=nSym, chromosome=nChr, go=nGO,
               version="0.1",
               maintainer="Nick Leigh <nicholas.leigh@med.lu.se>",
               author="Nick Leigh <nicholas.leigh@med.lu.se>",
               outputDir = "/Users/nickleigh/Library/CloudStorage/OneDrive-LundUniversity/Dokument/Data/2024/Genome/axolotl",
               tax_id="8296",
               genus="Ambystoma",
               species="mexicanum",
               goTable="go")



