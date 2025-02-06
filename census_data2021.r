# Load required libraries
library(readr)
library(dplyr)
library(tidyverse)
library(sf)
library(tmap)

# Load Data ----------------------------------------------------
setwd("/lustre06/project/6050814/chanalex/msc_thesis/Data/StatsCan")

# 2021 census division shapefile
censusDivisions_sf <- st_read("census_divisions2021.shp")

# 2021 census race data
census_data <- read_delim("CensusProfile2021-VisibleMinority_CensusDivisions.csv",
                          delim = ",")

# Filtering Data -----------------------------------------------
# Remove total counts
census_dataFiltered <- census_data %>%
  dplyr::filter(!RACE %in% c("Total - Visible minority", "Total—Visible minority", "Total visible minority population")))

# Aggregate racial categories to reduce memory requirements
# Combine "Chinese", "Korean", "Japanese" racialized individuals into "East Asian"
east_asian <- census_dataFiltered %>%
  dplyr::filter(RACE %in% c("Chinese", "Korean", "Japanese")) %>%
  group_by(REF_DATE, GEO, DGUID) %>%
  summarise(VALUE = sum(VALUE), .groups = "drop") %>%
  mutate(RACE = "East Asian")

# Combine "Filipino" and "Southeast Asian" racialized individuals
southeast_asian <- census_dataFiltered %>%
  dplyr::filter(RACE %in% c("Filipino", "Southeast Asian")) %>%
  group_by(REF_DATE, GEO, DGUID) %>%
  summarise(VALUE = sum(VALUE), .groups = "drop") %>%
  mutate(RACE = "Southeast Asian")

# Remove individual rows that were aggregated
census_dataClean <- census_dataFiltered %>%
  dplyr::filter(!RACE %in% c("Chinese", "Korean", "Japanese", "Filipino", "Southeast Asian"))

# Combine the aggregated racial categories with the cleaned data
census_dataFinal <- bind_rows(census_dataClean, east_asian, southeast_asian)

# Verify filtering
levels(census_data$RACE)
levels(census_dataFinal$RACE)

# Join the census race data with spatial data by DGUID (i.e., Dissemination Geography Unique Identifier)
censusRace_sf <- left_join(census_dataFinal, censusDivisions_sf, by = "DGUID")

# Generate dot density points ----------------------------------
# Function to generate random points inside a polygon (i.e., the census division)
# One dot represents 10 people (dot_ratio = 10).
generate_dots <- function(polygon, count, dot_ratio = 10){
  num_dots <- round(as.numeric(count) / dot_ratio)
  if(num_dots > 0){
    st_sample(polygon, size = num_dots, type = "random")
  } 
  else{
    NULL
  }
}

# Loop over each row of the joined data to generate points 
# for each race in each census division
dots_list <- list()
for(i in seq_len(nrow(censusRace_sf))){
  row_data <- censusRace_sf[i, ]
  pts      <- generate_dots(row_data$geometry, row_data$VALUE, dot_ratio = 10)
  
  if(!is.null(pts)){
    pts_sf <- st_sf(
      RACE     = row_data$RACE,
      GEO      = row_data$GEO,
      DGUID    = row_data$DGUID,
      geometry = pts
    )
    dots_list[[length(dots_list) + 1]] <- pts_sf
  }
}

# Combine all dot sf objects into one
dots_all <- do.call(rbind, dots_list)

# Create the dot density map -----------------------------------
# Following: https://walker-data.com/census-r/mapping-census-data-with-r.html
dot_map <- tm_shape(cd_sf) +
  tm_polygons(fill = "white", col = "gray", fill_alpha = 0.5) +
  tm_shape(dots_all) +
  tm_dots(col = "RACE", palette = "brewer.set1", size = 0.005) +
  tm_title("Racial Diversity in 2021 Census Divisions") +
  tm_layout(legend.outside = TRUE)

# Save the map
tmap_save(dot_map, filename = "dot_density_map.png")