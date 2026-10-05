# Load  packages
if(!requireNamespace("httr")) install.packages("httr")
if(!requireNamespace("jsonlite")) install.packages("jsonlite")
if(!requireNamespace("ggplot2")) install.packages("ggplot2")
if(!requireNamespace("dplyr")) install.packages("dplyr")

library(httr)
library(jsonlite)
library(ggplot2)
library(dplyr)

# ---------------------------
# 1. Function to retrieve UniProt domain features
# ---------------------------
get_uniprot_features <- function(uniprot_id){
  url <- paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id, "?format=json")
  resp <- GET(url)
  
  if (resp$status_code != 200)
    stop("UniProt request failed for ", uniprot_id, " (HTTP ", resp$status_code, ")")
  
  txt <- content(resp, as = "text", encoding = "UTF-8")
  json <- fromJSON(txt, simplifyDataFrame = FALSE)
  
  if (is.null(json$features)) return(NULL)
  
  df <- do.call(rbind, lapply(json$features, function(f) {
    start <- f$location$start$value
    end   <- f$location$end$value
    
    if (!is.null(start) && !is.null(end)) {
      data.frame(
        protein = uniprot_id,
        type = f$type,
        description = ifelse(is.null(f$description), f$type, f$description),
        start = as.numeric(start),
        end = as.numeric(end),
        stringsAsFactors = FALSE
      )
    } else NULL
  }))
  
  return(df)
}

# ---------------------------
# 2. UniProt IDs for species
# ---------------------------
protein_ids <- c(
  "A0A8J0QTA5",   # Xenopus tropicalis
  "A0AAV7U3D9",   # Pleurodeles waltl
  "A7KS68",       # Axolotl
  "P15919",       # Mouse
  "P15918"        # Human
)

species_names <- c("Xenopus", "Pleurodeles", "Axolotl", "Mouse", "Human")

# ---------------------------
# 3. Download features
# ---------------------------
all_features <- lapply(protein_ids, get_uniprot_features)
all_features <- do.call(rbind, all_features)

# Label species
all_features$protein_label <- factor(
  all_features$protein,
  levels = rev(protein_ids),
  labels = rev(species_names)
)

# ---------------------------
# 4. Collapse rare descriptions into "Other"
# ---------------------------
desc_counts <- table(all_features$description)
rare <- names(desc_counts[desc_counts < 4])

all_features$description_collapsed <- ifelse(
  all_features$description %in% rare, "Other", all_features$description
)

# ---------------------------
# 5. Keep ONLY domains shared by all species
# ---------------------------
shared_domains <- all_features %>%
  group_by(description_collapsed) %>%
  summarise(n_species = n_distinct(protein_label)) %>%
  filter(n_species == length(species_names)) %>%
  pull(description_collapsed)

all_features_common <- all_features %>%
  filter(description_collapsed %in% shared_domains)

# ---------------------------
# 6. Colors for domain descriptions (COLOR-SAFE, NO RED/GREEN)
# ---------------------------
unique_desc <- sort(unique(all_features_common$description_collapsed))

# A palette with **no red, no green**
safe_palette <- c(
  "steelblue1",  # blue
  "purple3",  # purple
  "brown3",  # brown
  "gray50",  # grey
  "#bcbd22",  # mustard yellow (not green)
  "#17becf"   # teal-blue
)

desc_colors <- setNames(
  safe_palette[seq_along(unique_desc)],
  unique_desc
)

# ---------------------------
# 7. Determine AA length for tick marks
# ---------------------------
max_len <- max(all_features_common$end)

# ---------------------------
# 8. Plot
# ---------------------------
ggplot(
  all_features_common,
  aes(
    x = start, xend = end,
    y = protein_label, yend = protein_label,
    color = description_collapsed
  )
) +
  
  # Backbone line starting at 0
  geom_segment(
    data = all_features_common %>%
      group_by(protein_label) %>%
      summarise(max_x = max(end)),
    aes(x = 0, xend = max_x, y = protein_label, yend = protein_label),
    inherit.aes = FALSE,
    size = 1.4,
    color = "black"
  ) +
  
  # Domain blocks
  geom_segment(size = 6) +
  
  # Domain colors
  scale_color_manual(values = desc_colors) +
  
  # X-axis ticks every 100 AA
  scale_x_continuous(
    breaks = seq(0, max_len + 50, by = 100),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  theme_classic() +
  labs(
    title = "Rag1 protein domains shared across species",
    x = "Amino Acid Position",
    y = "Rag1 protein",
    color = "Domain Description"
  ) +
  theme(
    axis.text.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )




