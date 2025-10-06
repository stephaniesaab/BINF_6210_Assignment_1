#Assignment 1: BINF6210: Software for Bioinformatics
#Date: October 10, 2025
#Author: Stephanie Saab
#Taxonomy of interest: The Canidae 
#Research Question: How similar are Canidae BIN compositions across continents?
#Bonus Discussion Question: Do domestic dogs show different barcode diversity patterns than wild canidae?


#References:
#Help Package for countrycode package: https://cran.r-project.org/web/packages/countrycode/refman/countrycode.html


# Downloading Packages if not installed:

#Add for loop here ---> 
install.packages('countrycode')

# Loading Libraries ===============
library(tidyverse)
library(dplyr)
library(vegan)
library(countrycode)

# Importing Data ===============

raw_data <- read_tsv("data/result.tsv")

# Cleaning Data =================

#Want to collect data on only the country collected (columns: geoid, country/ocean, country.iso, province/state) and BIN diversity (columns: bin_uri) and extra data on taxonomy (kingdom, phylum, class, order, family, genus, species). Remove NA values in geographical or BIN data, check for duplicates.
#Using the OR to ensure that we're getting all cases with at least one geographic identifier (country/ocean, country_iso, or province/state)
#It was observed that for cases where country/ocean was NA, the country_iso and province/state were also NA, so only need to filter for NA in country/ocean column
#Want to filter for valid country IDs (some invalid ones were found to be Exception - Culture, Exception - Zoological Park, Unrecoverable")
#Results in a table with 1343 rows (nrow(cleaned_data) = 1343)

cleaned_data <- raw_data %>%
  select(bin_uri, species, `country/ocean`) %>% 
  filter(!is.na(bin_uri) & !is.na(`country/ocean`)) %>% 
  filter(!`country/ocean`%in% c("Exception - Culture", "Exception - Zoological Park", "Unrecoverable")) 
  
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
#Want to filter these out as we want to identify biodiversity in the native habitats, not 
cleaned_data$Continent <- countrycode(
  cleaned_data$`country/ocean`,
  origin = "country.name.en",
  destination = "continent")

# Creating data frame with columns as unique continents, rows as unique BIN IDs, values as 0 (absence) or 1 (presence)

presence_absence_table <- cleaned_data %>% 
  select(bin_uri, Continent) %>% #This table only needs the BIN IDs and continents
  mutate(presence = 1) %>% #Make a new column with the presence value of the species
  distinct(bin_uri, Continent, .keep_all = TRUE) %>% #Remove duplicate species 
  pivot_wider(
    names_from = Continent, 
    values_from = presence, #New column we made is presence, 1 by default if present
    values_fill = 0
  )
#https://www.youtube.com/watch?v=XUQl31aLMhI
#Calculate similarity, 
dist_matrix <- vegdist(bin_table, method = "jaccard")


# 
# 
# continents <- c(unique(cleaned_data$Continent))
# species <- c(unique(cleaned_data$bin_uri))
# cleaned_data_tib <- tibble(
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

# TODO ====
#ASK if we can bring in other packages in R
#ASK if we can use loops (to add the install packages loop, so they'll install if they aren't loaded)
#Ask if that's enough?
#Add the testing for if there's country/ocean that's NA but country_iso that isn't?
#Ask about the GEOID, unrecoverable, exception - Culture and Exception - Zoo
#Add data summary?
#Try to make presence absence table with different package or in a pipe?


