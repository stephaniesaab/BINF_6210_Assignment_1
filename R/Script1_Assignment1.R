# Assignment 1: BINF6210: Software for Bioinformatics
# Date: October 13, 2025
# Author: Stephanie Saab
#Secondary Contributor: Iroayo Toki
#Reviewed by: Kexin Gong and Eman Tahir
# Taxonomy of interest: The Canidae
# Research Questions: How similar are Canidae BIN compositions across continental regions? Do domestic dogs (Canis familiaris) show a wider geographical diversity patterns than wild canidae species?

# Install and Load Packages ====
# Loop checks if the packages are already installed, if not it installs them, if they are it loads them.

required_packages <- c("tidyverse", "dplyr", "vegan", "countrycode", "viridis", "ggplot2", "igraph")
for (package in required_packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package)
  } else {
    library(package, character.only = TRUE)
  }
}

# Importing Data ===============
# Ensure working directory is set to source file location (".../BINF_6210_Assignment_1/R")
raw_data <- read_tsv("../data/result.tsv")

# Cleaning Data =================
# Want to collect data on the country of sample collection, BIN diversity, and species identity. Remove NA values and only get valid country IDs.

# Finding which columns contain "familiaris"
# Result: found that one entry lists "Canis familiaris" as species instead of subspecies.
cols_familiaris <- sapply(raw_data, function(x) any(grepl("familiaris", x, ignore.case = TRUE)))

# Finding if there are cases where country/ocean is NA but country_iso or province/state is not NA
# Result: 0, so only need to filter for cases where country/ocean is not NA
df_country_NA <- raw_data %>%
  select(`country/ocean`, country_iso, `province/state`) %>%
  filter(is.na(`country/ocean`) & (!is.na(country_iso) | !is.na(`province/state`)))
nrow(df_country_NA) #Returns 0

# Want to filter for valid country/ocean, and species or subspecies data, change the one row with species "Canis familiaris" to have subspecies be "Canis lupus" to be consistent in species naming
cleaned_data_bins <- raw_data %>%
  select(bin_uri, species, subspecies, `country/ocean`) %>%
  filter(!is.na(bin_uri) & !is.na(`country/ocean`) & (!is.na(species) | !(is.na(subspecies)))) %>%
  filter(!`country/ocean` %in% c("Exception - Culture", "Exception - Zoological Park", "Unrecoverable")) %>%
  mutate(
    subspecies = ifelse(
      test = species == "Canis familiaris",
      yes = "Canis lupus familiaris",
      no = subspecies
    ),
    species = ifelse(
      test = species == "Canis familiaris",
      yes = "Canis lupus",
      no = species
    )
  )

# Checkpoint for Cleaning data
# There are some cases in which subspecies is NA, that is okay for the data analysis
nrow(cleaned_data_bins) # Should be 1338
colnames(cleaned_data_bins) # Should be "bin_uri", "species", "subspecies", "country/ocean"
any(is.na(cleaned_data_bins$bin_uri) & is.na(cleaned_data_bins$`country/ocean`) &
  is.na(cleaned_data_bins$species)) # Should be FALSE


# Uniqueness of Data
# Check if the size of the dataset is acceptable for analysis (>10 records of each variable)

length(unique(cleaned_data_bins$species)) # 25 unique species in the sample
length(unique(cleaned_data_bins$bin_uri)) # 24 unique BIN IDs
length(unique(cleaned_data_bins$`country/ocean`)) # 67 unique country/oceans

# Data Analysis =====================

## BIN Richness per continent ====
# Building a presence-absence table to identify the diversity of the sample

#Converting geographic data into continents by creating a new column with the corresponding continents of collection using countrycode package, using destination as "region" instead of "continent" to separate the Americas
cleaned_data_bins$Continent <- countrycode(
  cleaned_data_bins$`country/ocean`,
  origin = "country.name.en",
  destination = "region"
)

#Set Reunion to Sub-Saharan Africa
cleaned_data_bins <- cleaned_data_bins %>%
  mutate(Continent = ifelse(
    test = `country/ocean` == "Reunion",
    yes = "Sub-Saharan Africa",
    no = Continent
  ))

# Create a data frame with columns as unique continents, rows as unique BIN IDs, values as 0 (absence) or 1 (presence)
df_presence_absence <- as.data.frame(cleaned_data_bins %>%
  select(bin_uri, Continent) %>%
  mutate(presence = 1) %>% # Make a new column with the presence value of the species
  distinct(bin_uri, Continent, .keep_all = TRUE) %>% # Remove duplicate species
  pivot_wider(
    names_from = Continent,
    values_from = presence, # New column we made is presence, 1 by default if present
    values_fill = 0
  ))

# Setting the BIN IDs as row names (Original code)
#rownames(df_presence_absence) <- df_presence_absence$bin_uri
#df_presence_absence <- df_presence_absence[, !(names(df_presence_absence) %in% "bin_uri")]

#edit 1: Shorten rowname conversion 
df_presence_absence <- df_presence_absence %>% column_to_rownames("bin_uri")


# Checkpoint for presence-absence table
nrow(df_presence_absence) # Should be 24 as there are 24 unique BINs
ncol(df_presence_absence) # Should be 7

### Richness per continent ========

# Calculate richness per continent (column)
continent_richness <- colSums(df_presence_absence)

# Create df to plot richness (in Visualization section)
df_richness <- data.frame(
  Continent = names(continent_richness),
  Richness = as.numeric(continent_richness)
)

# Arrange continents in descending order of richness, convert into factor to stay ordered in plot
df_richness <- df_richness %>%
  arrange(desc(Richness)) %>%
  mutate(Continent = factor(Continent, levels = Continent))

# Checkpoint for richness
max(df_richness$Richness) # Should be 10 (East Asia & Pacific)
min(df_richness$Richness) # Should be 3 (South Asia)

### Dissimilarity between continents ====
# Shows how similar or distinct the Canidae BIN compositions are between continental regions, based on shared and unique BINs

# Get distance matrix with pairwise dissimilarities between continents (lower numbers = more similarity)
jaccard_dist_cont <- vegdist(t(df_presence_absence), method = "jaccard")

# Make it a matrix so it's got complete data values
jaccard_matrix_continent <- as.matrix(jaccard_dist_cont)

# Make a tidy and complete set of the pairwise comparisons, to later plot a heatmap (Visualization section)
df_heatmap_cont <- as.data.frame(as.table(jaccard_matrix_continent)) %>%
  rename("Continent1" = Var1, "Continent2" = Var2, "Dissimilarity" = Freq)

# Checkpoint for dissimilarity analysis
max(jaccard_matrix_continent) # Should give 0.9285714
min(jaccard_matrix_continent) # Should give 0
nrow(df_heatmap_cont) # Should give 49 (7 pairwise comparisons)

## Canidae familiaris distribution ====
# Compare the number of BINs or average range size for C. familiaris vs. other species.


# Create the data frame that will be used to make the edges and nodes: Get a count for the occurrences for each taxon-continent pair to be able to get the highest weight nodes later and create a unified column with the taxonomic unit
df_edges <- cleaned_data_bins %>%
  mutate(
    taxon = ifelse(
      test = !is.na(subspecies),
      yes = subspecies,
      no = species
    )
  ) %>%
  select(taxon, Continent) %>%
  group_by(taxon, Continent) %>%
  summarise(weight = n(), .groups = "drop")

# Checkpoint for creating df_edges
max(df_edges$weight) # Should return 338
min(df_edges$weight) # should return 1
nrow(df_edges) # Should return 55

### Build bipartite network ====
# Build a graph that will be plotted in Visuazliations section, it creates edges with the first two columns (taxon and continents)
ls_g_ntwk <- graph_from_data_frame(df_edges)

# Assign the weight attribute based on the record counts
E(ls_g_ntwk)$weight <- df_edges$weight

# Distinguish the nodes as either taxonomic unit (species or subspecies) or continents, based on the name attribute of the nodes (TRUE = continent, FALSE = taxon)
V(ls_g_ntwk)$type <- V(ls_g_ntwk)$name %in% df_edges$Continent

# Assign the nodes a colour based on their identity, differentiate Canis lupus familiaris
V(ls_g_ntwk)$color <- ifelse(
  test = V(ls_g_ntwk)$type == TRUE,
  yes = "gold",
  no = "blue"
)
V(ls_g_ntwk)$color[V(ls_g_ntwk)$name == "Canis lupus familiaris"] <- "red"

# Filter for strong edges and non-isolated vertices to make the plot readable
weight_threshold <- 10 # Arbitrary threshold set
ls_g_high <- delete_edges(ls_g_ntwk, E(ls_g_ntwk)[weight <= weight_threshold])
ls_g_high <- delete_vertices(ls_g_high, which(degree(ls_g_high) == 0))


# Checkpoint for edge and vertex counts in graph
ecount(ls_g_ntwk) # Should return 55 edges (all edges)
ecount(ls_g_high) # Should return 18 (only strong edges)
vcount(ls_g_ntwk) # Should return 42
vcount(ls_g_high) # Should return 17 (only connected vertices)

# Visualizations =========================

## Bar plot for BIN Richness per continent ====

ggplot(df_richness, aes(x = Continent, y = Richness)) +
  geom_col() +
  theme_classic() +
  labs(
    title = "Canidae BIN Richness per Continental Region",
    x = "Continental Region",
    y = "BIN Richness"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1.05),
    legend.position = "none"
  )

# Save plot
ggsave("../figures/plot_continent_richenss.png", width = 6, height = 4, dpi = 300)

## Heat map for Jaccard disimilarities ====
ggplot(df_heatmap_cont, aes(x = Continent1, y = Continent2, fill = Dissimilarity)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_classic() +
  labs(
    title = "Jaccard Dissimilarity Score Heatmap",
    x = "Continental Region",
    y = "Continental Region",
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1.05)
  )

# Save heatmap
ggsave("../figures/heatmap_continent_jaccard.png", width = 6, height = 4, dpi = 300)


## Bipartite Network for Canis familiaris and species range ====
# Ensure it's reproducible with set.seed function, compute the positions of each node in a matrix of 2 columns (species on one side and continents on other). Scale the edge thickness by weight, distinguish shape by identity. Save manually if needed.
set.seed(42)
plot(ls_g_high,
  layout = layout_as_bipartite(ls_g_high)[, c(2, 1)],
  vertex.shape = ifelse(
    test = V(ls_g_high)$type == "TRUE",
    yes = "square",
    no = "circle"
  ),
  vertex.size = 7,
  vertex.label.color = "black",
  vertex.label.cex = 1.0,
  edge.color = "gray70",
  main = "Network of Domestic vs. Wild Canidae Species \n in Continental regions"
)

#edit 2 Create Abundance dataframe and calculate Shannon's diversity index(H')
# It quantifies community diversity by incorporating both richness (number of BINs present) and evenness (how evenly individuals are distributed across BINs). Higher H' values indicate a community with more BIN types and more balanced abundances.
#Create dataframe with BIN Richness(columns) by continent(rows) 
df_abundance <- as.data.frame(cleaned_data_bins) %>% 
  group_by(Continent, bin_uri) %>% 
  count() %>% 
  pivot_wider(
    names_from = bin_uri, 
    values_from = n, 
    values_fill = 0 )
#Make continent row names
df_abundance <- df_abundance %>% column_to_rownames("Continent")
#add Shannons index 

df_abundance <- df_abundance %>% mutate(Shannon = vegan::diversity(df_abundance, index ="shannon"))
df
#Highest BIN diversity in the Middle east and   North Africa
print(df_abundance$Shannon)

#Edit 3 : Bar plot for shannons index 
df_abundance$Continent <- rownames(df_abundance)

ggplot(df_abundance, aes(x = Continent, y = Shannon, fill = Continent)) +
  geom_col() +
  labs(title = "Shannon Diversity Index by Continent",
       y = "Shannon Index (H')",
       x = "") +
  theme_minimal()

ggsave("../figures/Shannon_diversity_barplot.png", width = 13.5, height = 8, dpi = 300)

