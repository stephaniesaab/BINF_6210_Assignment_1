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
cleaned_data <- raw_data %>% 
  select(bin_uri, family, genus, species, `country/ocean`, country_iso, `province/state`) %>% 
  filter(!is.na(bin_uri) & (!is.na(`country/ocean`) | !is.na(country_iso) | !is.na(`province/state`))) 

# Data Analysis =====================

## Building a presence/absence table ====
#A table to identify the diversity of the sample, in which rows are regions (continents), columns are BIN IDs, values are 1 for presence and 0 for absence

### Converting geographic data into continents





# TODO ====
#ASK if we can bring in other packages in R
#ASK if we can use loops (to add the install packages loop, so they'll install if they aren't loaded)
#Ask if that's enough?



