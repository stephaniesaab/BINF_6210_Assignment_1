#Assignment 1: BINF6210: Software for Bioinformatics
#Date: October 10, 2025
#Author: Stephanie Saab
#Taxonomy of interest: The Canidae 
#Research Question: How similar are Canidae BIN compositions across continents?
#Bonus Discussion Question: Do domestic dogs show different barcode diversity patterns than wild canidae?


#References:
#Help Package for countrycode package: https://cran.r-project.org/web/packages/countrycode/refman/countrycode.html


# Downloading Packages if not installed: Loop checks if the packages are already installed, if not it installs them, if they are it loads them.

required_packages <- c('tidyverse', 'dplyr', 'vegan', 'countrycode', 'viridis', 'ggplot2', 'igraph')
for (package in required_packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package)
  }
  else {
    library(package, character.only = TRUE)
  }
}

# Importing Data ===============
#Ensure working directory is set to source file location
raw_data <- read_tsv("../data/result.tsv")

# Cleaning Data =================

#Want to collect data on only the country collected (columns: geoid, country/ocean, country.iso, province/state) and BIN diversity (columns: bin_uri) and extra data on taxonomy (kingdom, phylum, class, order, family, genus, species). Remove NA values in geographical or BIN data, check for duplicates.
#Using the OR to ensure that we're getting all cases with at least one geographic identifier (country/ocean, country_iso, or province/state)
#It was observed that for cases where country/ocean was NA, the country_iso and province/state were also NA, so only need to filter for NA in country/ocean column
#Want to filter for valid country IDs (some invalid ones were found to be Exception - Culture, Exception - Zoological Park, Unrecoverable")
#Filtering for NA in the country/ocean column and bin_uri columns as they are most important
#Results have 1343 rows (nrow(cleaned_data_bins) = 1343)

cleaned_data_bins <- raw_data %>%
  select(bin_uri, species, `country/ocean`) %>% 
  filter(!is.na(bin_uri) & !is.na(`country/ocean`)) %>% 
  filter(!`country/ocean`%in% c("Exception - Culture", "Exception - Zoological Park", "Unrecoverable")) 
nrow(cleaned_data_bins) #Should be 1343

# Summarizing Data ==================
#Summarize and describe data (counts, means, etc.)

# Data Analysis =====================

## Building a presence/absence table ====
#A table to identify the diversity of the sample, in which rows are regions (continents), columns are BIN IDs, values are 1 for presence and 0 for absence

### Converting geographic data into continents ====
#Creating a new column with the corresponding continents
#Using code of country name in English
#Destination is the continent as defined in the World Bank Development Indicators
#This initially didn't match some values as there are some values which are not country codes: "Exception - Culture, Exception - Zoological Park, Unrecoverable"
#Want to filter these out as we want to identify biodiversity in the native habitats
cleaned_data_bins$Continent <- countrycode(
  cleaned_data_bins$`country/ocean`,
  origin = "country.name.en",
  destination = "continent")

#Removing NAs from species column for species analysis in part 3
#Should give table of 1334 observations

cleaned_data_species <- cleaned_data_bins %>% 
  filter(!is.na(species)) %>% 
  select(species, Continent)

### Creating presence-absence table ====
#data frame with columns as unique continents, rows as unique BIN IDs, values as 0 (absence) or 1 (presence)

df_presence_absence <- as.data.frame(cleaned_data_bins %>% 
  select(bin_uri, Continent) %>% #This table only needs the BIN IDs and continents
  mutate(presence = 1) %>% #Make a new column with the presence value of the species
  distinct(bin_uri, Continent, .keep_all = TRUE) %>% #Remove duplicate species 
  pivot_wider(
    names_from = Continent, 
    values_from = presence, #New column we made is presence, 1 by default if present
    values_fill = 0
  )
)

#Setting the BIN IDs as row names
rownames(df_presence_absence) <- df_presence_absence$bin_uri
df_presence_absence <- df_presence_absence[, !(names(df_presence_absence) %in% "bin_uri")]
 

#Save the tsv (uncomment to run again)
#write_tsv(presence_absence_table, file = "../figures/Presence_absence_table.tsv")

### Basic data summary: Richness per continent ========


#Calculate richness per continent (column)
continent_richness <- colSums(df_presence_absence)

#Create df to plot
df_richness <- data.frame(
  Continent = names(continent_richness),
  Richness = as.numeric(continent_richness)
)

#Arrange continents in descending order of richness, convert into factor to stay ordered in plot
df_richness <- df_richness %>% 
  arrange(desc(Richness)) %>% 
  mutate(Continent = factor(Continent, levels = Continent))


### Disimiarlity between continents ====
#Asks - how many species do two continents share, and how many are unique to each? -> overlap and uniqueness of species between regions

#Gets a distance matrix showing pairwise dissimilarities between continents (lower numbers = more similarity), transpose to compare continents
jaccard_dist_cont <- vegdist(t(df_presence_absence), method = "jaccard")

#Make it a matrix so it's got complete data value, then as a table then dataframe to get a tidy and complete set of pairwise comparisons
#Highest disimilarity is between Africa and Americas
#Lowest disimilarity is between Europe and Asia 
jaccard_matrix_continent <- as.matrix(jaccard_dist_cont)
df_heatmap_cont <- as.data.frame(as.table(jaccard_matrix_continent)) %>% 
  rename("Continent1" = Var1, "Continent2" = Var2, "Disimilarity" = Freq)

### Looking at Canidae familiaris distribution and average range size ====

#Compare the number of BINs or average range size for C. familiaris vs. other species.Expect dogs to have global distribution, while wild species are continent-restricted.


#Making bipartite network, with nodes representing species and continents, and edges representing presence

df_species_ntwk <- cleaned_data_species %>% 
  distinct() #Remove duplicate connections

g_species_ntwk <- graph_from_data_frame(df_species_ntwk)




#Setting vertices type based on continent or species
V(g_species_ntwk)$type <- ifelse(
  V(g_species_ntwk)$name %in% df_species_ntwk$Continent, "continent", "species"
  )

#Set vertex colour based on type
V(g_species_ntwk)$color <- ifelse(
  V(g_species_ntwk)$type == "continent", "gold", "skyblue")
  )

#Make Canis familiaris stand out
V(g_species_ntwk)$color[V(g_species_ntwk)$name == "Canis  familiaris"] <- "red"






v(g)$type 
#Create bipartite graph
species_bip_graph <- make_bipartite_graph(edges = df_species_ntwk)


# 
# 
# 
# nmds <- metaMDS(t(df_presence_absence), distance = "jaccard")
# plot(nmds)
# 
# plot(hclust(jaccard_dist_continent))
# 
# 
# 
# df_bin_table <- df_presence_absence %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   rownames_to_column(var = "Continent")
# 
# #Renaming the columns with the first row
# bin_table <- bin_table %>% 
#   rename_with(function(x) bin_table[1, ])
# 
# #Add 
# bin_richness <- lapply(bin_table, sum)
# 
# 
# bin_richness <- vegdist(bin_table, method = "jaccard")


# 
# 
# continents <- c(unique(cleaned_data_bins$Continent))
# species <- c(unique(cleaned_data_bins$bin_uri))
# cleaned_data_bins_tib <- tibble(
#   bin_uri = species,
#   Continent = continents
# )
# df_pres_abs <- data.frame(matrix(NA, nrow = length(df_rows), ncol = length(df_cols)))
# colnames(df_pres_abs) <- df_cols
# row.names(df_pres_abs) <- df_rows
# presence_absence_table <- df_pres_abs %>% 
#   mutate(present = 1) %>% #Make a new column with values of presence
#   pivot_wider(
#     id_cols = df_cols,
#     names_from = df_rows,
#     values_from = present,
#     values_file = 0
#   )
# 
#   

  


#dist_matrix 
# Visualizations =========================

## Bar plot for BIN Richness per continent ====

ggplot(df_richness, aes(x = Continent, y = Richness))+
  geom_col(fill = "violet") + 
  theme_classic() + 
  labs(title = "Canidae BIN Richness per Continent",
       x = "Continent", 
       y = "BIN Richness") +
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none")

#Save plot
ggsave("../figures/plot_continent_richenss.png", width = 6, height = 4, dpi = 300)
 
## Heat map for Jaccard disimilarities ====
ggplot(df_heatmap_cont, aes(x = Continent1, y = Continent2, fill = Disimilarity))+
  geom_tile() +
  scale_fill_viridis_c() + 
  theme_classic()+
  labs(title = "Jaccard Disimilarity Score",
       x = "Continent",
       y = "Continent",)

#Save heatmap
ggsave("../figures/heatmap_continent_jaccard.png", width = 6, height = 4, dpi = 300)


## Bipartite Network for Canis familiaris and species range ====
plot(g_species_ntwk,
     layout = layout_with_fr(g_species_ntwk),
     vertex.size = 6,
     vertex.label.cex = 0.7,
     vertex.label.color = 'black', 
     edge.color = 'gray70',
     main = 'Network of Canidae Species & Continents range')

# TODO ====
#ASK if we can bring in other packages in R
#ASK if we can use loops (to add the install packages loop, so they'll install if they aren't loaded)
#Ask if that's enough?
#Add the testing for if there's country/ocean that's NA but country_iso that isn't?
#Ask about the GEOID, unrecoverable, exception - Culture and Exception - Zoo
#Add data summary? --> can add other stuff but put a note saying it was inconclusive
#Try to make presence absence table with different package or in a pipe?
#We include summary / exploration of data?
#Maybe take out the bar plot if they're only marking 3?


