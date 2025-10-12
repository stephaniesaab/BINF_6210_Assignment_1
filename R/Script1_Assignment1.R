#Assignment 1: BINF6210: Software for Bioinformatics
#Date: October 13, 2025
#Author: Stephanie Saab
#Taxonomy of interest: The Canidae 
#Research Question: How similar are Canidae BIN compositions across continents?
#Second Hypothesis: Do domestic dogs (Canis familiaris) show a wider geographical barcode diversity patterns than wild canidae?


#References:
#Help Package for countrycode package: https://cran.r-project.org/web/packages/countrycode/refman/countrycode.html


# Downloading and loading required packages: Loop checks if the packages are already installed, if not it installs them, if they are it loads them.

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
#Ensure working directory is set to source file location (".../BINF_6210_Assignment_1/R")
raw_data <- read_tsv("../data/result.tsv")

# Cleaning Data =================

#Want to collect data on only the country of sample collection and BIN diversity and extra some potentially useful data on taxonomy (kingdom, phylum, class, order, family, genus, species). Remove NA values and only get valid country IDs.

#Finding which columns contain "familiaris"
#Result: found that some list "familiaris" under species and some under subspecies, so include both in data, one row contains "Canis familiaris" as species instead of subspecies
cols_familiaris <- sapply(raw_data, function(x) any(grepl("familiaris", x, ignore.case = TRUE)))

#Finding if there are cases where country/ocean is NA but country_iso or province/state is not NA
#Result: 0, so only need to filter for cases where country/ocean is not NA
df_country_NA <- raw_data %>% 
  select(`country/ocean`, country_iso, `province/state`) %>% 
  filter(is.na(`country/ocean`) & (!is.na(country_iso) | !is.na(`province/state`)))
nrow(df_country_NA)

#Want to filter for valid country/ocean, and species or subspecies data, change the one row with species "Canis familiaris" to have subspecies be "Canis lupus" to be consistent in species naming
cleaned_data_bins <- raw_data %>%
  select(bin_uri, species, subspecies, `country/ocean`) %>% 
  filter(!is.na(bin_uri) & !is.na(`country/ocean`) & (!is.na(species) | !(is.na(subspecies)))) %>% 
  filter(!`country/ocean`%in% c("Exception - Culture", "Exception - Zoological Park", "Unrecoverable")) %>% 
   mutate(
     subspecies = ifelse(
       test = species == "Canis familiaris", 
       yes = "Canis lupus familiaris", 
       no = subspecies),
     species = ifelse(
       test = species == "Canis familiaris",
       yes = "Canis lupus", 
       no = species)
     )

#Checkpoint for Cleaning data
#There are some cases in which subspecies is NA, that is okay for the data analysis
nrow(cleaned_data_bins) #Should be 1338
colnames(cleaned_data_bins) #Should be "bin_uri", "species", "subspecies", "country/ocean"
any(is.na(cleaned_data_bins$bin_uri) & is.na(cleaned_data_bins$`country/ocean`)
    & is.na(cleaned_data_bins$species)) #Should be FALSE


#Uniqueness of Data
#Check if the size of the dataset is acceptable for analysis (>10 records of each variable)

length(unique(cleaned_data_bins$species)) #26 unique species in the sample
length(unique(cleaned_data_bins$bin_uri)) #24 unique BIN IDs
length(unique(cleaned_data_bins$`country/ocean`)) #67 unique country/oceans

# Data Analysis =====================

## BIN Richness per continent ====
#Building a presence-absence table to identify the diversity of the sample

#Converting geographic data into continents by creating a new column with the corresponding continents of collection using countrycode package, using destination as "region" instead of "continent" to separate the Americas
cleaned_data_bins$Continent <- countrycode(
  cleaned_data_bins$`country/ocean`,
  origin = "country.name.en",
  destination = "region")

#Got NA continent for two samples in which country/ocean = "Reunion", manually set them to Africa
cleaned_data_bins <- cleaned_data_bins %>% 
  mutate(Continent = ifelse(
    test = `country/ocean` == "Reunion",
    yes = "Sub-Saharan Africa", 
    no = Continent))

#Create a data frame with columns as unique continents, rows as unique BIN IDs, values as 0 (absence) or 1 (presence)
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

#Checkpoint for presence-absence table
nrow(df_presence_absence) #Should be 24 as there are 24 unique BINs
ncol(df_presence_absence) #Should be 7

### Richness per continent ========

#Calculate richness per continent (column)
continent_richness <- colSums(df_presence_absence)

#Create df to plot richness (in Visualization section)
df_richness <- data.frame(
  Continent = names(continent_richness),
  Richness = as.numeric(continent_richness)
)

#Arrange continents in descending order of richness, convert into factor to stay ordered in plot
df_richness <- df_richness %>% 
  arrange(desc(Richness)) %>% 
  mutate(Continent = factor(Continent, levels = Continent))

#Checkpoint for richness
max(df_richness$Richness) #Should be 10 (East Asia & Pacific)
min(df_richness$Richness) #Should be 3 (South Asia)

### Dissimilarity between continents ====
#Finding how many species do two continents share, and how many are unique to each

#Get a distance matrix showing pairwise dissimilarities between continents (lower numbers = more similarity), transpose df_presence_absence to compare continents
jaccard_dist_cont <- vegdist(t(df_presence_absence), method = "jaccard")

#Make it a matrix so it's got complete data values
jaccard_matrix_continent <- as.matrix(jaccard_dist_cont)

#Make a table then dataframe of the matrix to get a tidy and complete set of the pairwise comparisons, to later plot a heatmap (Visualization section)
df_heatmap_cont <- as.data.frame(as.table(jaccard_matrix_continent)) %>% 
  rename("Continent1" = Var1, "Continent2" = Var2, "Disimilarity" = Freq)

#Checkpoint for dissimilarity
max(jaccard_matrix_continent) #Should give 0.9285714
min(jaccard_matrix_continent) #Should give 0
nrow(df_heatmap_cont) #Should give 49 (7 pairwise comparisons)

## Canidae familiaris distribution ====
#Compare the number of BINs or average range size for C. familiaris vs. other species. It is expect dogs to have global distribution, while wild species are continent-restricted.


#Create the data frame that will be used to make the edges and nodes: Get a count for the occurrences for each taxon-continent pair to be able to get the highest weight nodes later and create a unified column with the taxonomic unit
df_edges <- cleaned_data_bins %>% 
  mutate(
    taxon = ifelse(
      test = !is.na(subspecies), 
      yes = subspecies, 
      no = species)
  ) %>% 
  select(taxon, Continent) %>% 
  group_by(taxon, Continent) %>% 
  summarise(weight = n(), .groups = "drop")

#Checkpoint for creating df_edges
max(df_edges$weight) #Should return 338
min(df_edges$weight) #should return 1
nrow(df_edges) #Should return 55

### Build bipartite network ====
#Build a graph that will be plotted in Visuazliations section, it creates edges with the first two columns (taxon and continents)
ls_g_ntwk <- graph_from_data_frame(df_edges)

#Assign the weight attribute based on the counts done previously
E(ls_g_ntwk)$weight <- df_edges$weight

#Distinguish the nodes as either taxonomic unit (species or subspecies) or continents, based on the name attribute of the nodes, set type as logical variable (TRUE = continent, FALSE = taxon)
V(ls_g_ntwk)$type <- V(ls_g_ntwk)$name %in% df_edges$Continent

#Assign the nodes a colour based on their identity
V(ls_g_ntwk)$color <- ifelse(
  test = V(ls_g_ntwk)$type == TRUE, 
  yes = "gold", 
  no = "blue"
)

#Assign the Canis lupus familiaris taxon a different colour to highlight it
V(ls_g_ntwk)$color[V(ls_g_ntwk)$name == "Canis lupus familiaris"] <- "red"

#Filter for strong edges only to make the plot readable
weight_threshold <- 10 #Arbitrary threshold set to remove weak connections
ls_g_high <- delete_edges(ls_g_ntwk, E(ls_g_ntwk)[weight <= weight_threshold])

#Remove any isolated vertices to clean up the plot
ls_g_high <- delete_vertices(ls_g_high, which(degree(ls_g_high) == 0))


#Checkpoint for edge and vertex counts in graph
ecount(ls_g_ntwk) #Should return 55 edges (all edges)
ecount(ls_g_high) #Should return 18 (only strong edges)
vcount(ls_g_ntwk) #Should return 42
vcount(ls_g_high) #Should return 17 (only connected vertices)

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
#Ensure it's reproducible with set.seed function, compute the positions of each node in a matrix of 2 columns, swap the positions of those columns to flip the axes (species on one side and continents on other). Use the weights to scale the edge thickness to visualize stronger occurrences
set.seed(42)
plot(ls_g_high,
     layout = layout_as_bipartite(ls_g_high)[, c(2, 1)],
     vertex.shape = ifelse(
       test = V(ls_g_high)$type == "continent",
       yes = "square", 
       no = "circle"
     ),
     vertex.size = 6,
     vertex.label.cex = 0.7,
     vertex.label.color = 'black', 
     edge.color = 'gray70',
     main = 'Network of Domestic vs. Wild Canidae Species in Continents'
     )

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




# 
# #Making bipartite network, with nodes representing species and continents, and edges representing presence
# df_species_ntwk <- cleaned_data_bins %>% 
#   distinct() #Remove duplicate connections
# 
# 
# 
# 
# 
# g_species_ntwk <- graph_from_data_frame(df_species_ntwk)
# 
# #Setting vertices type based on continent or species
# V(g_species_ntwk)$type <- ifelse(
#   V(g_species_ntwk)$name %in% df_species_ntwk$Continent, "continent", "species"
#   )
# 
# #Set vertex colour based on type
# V(g_species_ntwk)$color <- ifelse(
#   V(g_species_ntwk)$type == "continent", "gold", "skyblue")
#   )
# 
# #Make Canis familiaris stand out
# V(g_species_ntwk)$color[V(g_species_ntwk)$name == "Canis lupus familiaris"] <- "red"
# 
# 
# 
# 

# 
# v(g)$type 
# #Create bipartite graph
# species_bip_graph <- make_bipartite_graph(edges = df_species_ntwk)
# 
# 
# # 
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


