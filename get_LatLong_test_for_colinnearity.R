##### Set your Working Directory ############################
setwd("~/Desktop/Projects/Valerianaceae_biogeography")

#############################################################
####### you'll need to install there packages ###############
####### with install.packages() command #####################
#############################################################
library(tidyverse)
library(rgbif)
library(dplyr)
library(tibble)
library(geodata)
library(terra)
library(corrplot)
library(usdm)
library(ENMTools)
library(caret)
library(Hmisc)


###  Species List
species_list <- c("Valeriana hebecarpa","Valeriana leucocarpa","Valeriana grandifolia","Valeriana virescens","Valeriana laxiflora","Valeriana lepidota","Valeriana bracteosa","Valeriana bridgesii","Valeriana papilla","Valeriana vaga","Valeriana verticillata","Valeriana polemoniifolia","Valeriana boelckei","Valeriana carnosa","Valeriana fonkii","Valeriana clarionifolia","Plectritis macrocera","Plectritis congesta","Valeriana moyanoi","Valeriana chilensis","Valeriana chaerophylloides","Valeriana philippiana","Valeriana polystachya","Valeriana graciliceps","Valeriana radicalis","Valeriana adscendens","Valeriana secunda","Valeriana plantaginea","Valeriana pilosa","Valeriana tatamana","Valeriana hirtella","Valeriana microphylla","Valeriana connata","Valeriana lapathifolia","Valeriana nivalis","Valeriana arborea","Valeriana tomentosa","Valeriana hornschuchiana","Valeriana stricta","Valeriana crispa","Valeriana clematitis","Valeriana naidae","Valeriana pyramidalis","Valeriana effusa","Valeriana interrupta","Valeriana gallinae","Valeriana rzedowskiorum","Valeriana densiflora","Valeriana sorbifolia","Valeriana selerorum","Valeriana mexicana","Valeriana urticifolia","Valeriana candolleana","Valeriana scandens","Valeriana apiifolia","Valeriana tanacetifolia","Valeriana bryophila","Valeriana robertianifolia","Valeriana palmeri","Valeriana albonervata","Valeriana ceratophylla","Valeriana texana","Valeriana edulis","Valeriana prionophylla","Valeriana stuckertii","Valeriana polybotrya","Valeriana macrorhiza","Valeriana niphobia","Valeriana aretioides","Valeriana bracteata","Valeriana rumicoides","Valeriana coarctata","Valeriana stenoptera","Valeriana officinalis","Valeriana minutiflora","Valeriana kawakamii", "Valeriana flaccidissima","Valeriana trichostoma","Valeriana sichuanica","Valeriana fauriei","Valeriana collina","Valeriana tuberosa","Valeriana asarifolia","Valeriana phu","Valeriana occidentalis","Valeriana californica","Valeriana sitchensis","Valeriana pauciflora","Valeriana acutiloba","Valeriana dioica","Valeriana arizonica","Centranthus trinervis","Centranthus sieberi","Centranthus ruber","Centranthus lecoqii","Centranthus angustifolius","Centranthus longiflorus","Centranthus macrosiphon","Centranthus calcitrapae","Valeriana tripteris","Valeriana supina","Valeriana montana","Valeriana pyrenaica","Valerianella stenocarpa","Valerianella texana","Valerianella radiata","Valerianella amarella","Valerianella dentata","Valerianella eriocarpa","Valerianella vesicaria","Valerianella coronata","Valerianella discoidea","Valerianella carinata","Valerianella locusta","Valerianella pumila","Patrinia saniculifolia","Fedia cornucopiae","Valeriana celtica","Valeriana saxatilis","Valeriana elongata","Nardostachys jatamansi","Patrinia gibbosa","Patrinia villosa","Patrinia triloba","Patrinia scabiosifolia","Triplostegia glandulifera")

# Create an empty list to store results
gbif_data_list <- list()

# Loop through the species list to fetch data

for (spp in species_list) {
  
  cat("Downloading data for:", spp, "\n")
  
  data <- try(
    occ_data(
      scientificName = spp,
      hasCoordinate = TRUE,
      limit = 10000
    ),
    silent = TRUE
  )
  
  if (inherits(data, "try-error")) {
    cat("GBIF error for:", spp, "\n")
    next
  }
  
  if (!is.null(data$data) && nrow(data$data) > 0) {
    
    coords <- data$data %>%
      dplyr::select(decimalLongitude, decimalLatitude) %>%
      dplyr::filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
      dplyr::mutate(species = spp)
    
    gbif_data_list[[spp]] <- coords
    
  } else {
    cat("No coordinate data found for:", spp, "\n")
  }
  
  Sys.sleep(0.2)  # avoid GBIF throttling
}


# Combine all species data frames into one table
final_occurrence_table <- bind_rows(gbif_data_list)

# View the head of the final table
head(final_occurrence_table)

# View the column names to confirm structure
names(final_occurrence_table)

write.csv(final_occurrence_table, "~/Desktop/Projects/Valerianaceae_biogeography/my_coordinates.csv", row.names = FALSE)

#final_occurrence_table <- read.csv("~/Desktop/Projects/Valerianaceae_biogeography/my_coordinates.csv")
###########################################################################
############ Test for co-linearity ########################################
############ Ignore this for now !!!#######################################
###########################################################################


# 1. Download WorldClim bioclimatic variables (bio) at 10-minute resolution
# Options for res: 0.5, 2.5, 5, 10

clim_data <- geodata::worldclim_global(global = TRUE, var = "bio", res = 10)

# 2. Prepare your coordinate table (ensure columns are lon, lat)
# Example: my_coords <- data.frame(lon = c(-74.0, -0.1), lat = c(40.7, 51.5))
points <- vect(final_occurrence_table, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")

# 3. Extract climate values for your points
extracted_values <- raster::extract(clim_data, points)

# 4. Create new data frame with occurrence data, species, and extracted environmental data
new_df <- merge(final_occurrence_table, extracted_values)

# Remove the 'ID' column generated by extract() and handle NAs

cor_matrix <- cor(extracted_values[, -1], use = "complete.obs")
print(cor_matrix)

corrplot(cor_matrix, method = "circle")


##### caret Package############
#install.packages("caret")

findCorrelation(cor_matrix, cutoff = 75)

# Install and load the Hmisc package if you haven't already
#install.packages("Hmisc")

rcorr(as.matrix(cor_matrix), type = "pearson")

vif(extracted_values)

# Automatically excludes variables with VIF > 10
vif_result <- vifstep(extracted_values, th = 10)
print(vif_result)

# Get the names of the remaining 'independent' variables
keep_vars <- vif_result@results$Variables


