install.packages("geodata")   # if not already installed

install.packages("FactoMineR")
install.packages("factoextra")

library(tidyverse)
library(rgbif)
library(terra)
library(dplyr)
library(ggplot2)
library(FactoMineR)
library(factoextra)


species_list <- c(
  "Valeriana officinalis",
  "Valeriana sitchensis",
  "Valeriana edulis"
)


bioclim <- geodata::worldclim_global(var = "bio", res = 2.5, path = ".")

all_species_data <- list()

for (spp in species_list) {
  
  cat("Downloading:", spp, "\n")
  
  occ <- occ_data(
    scientificName = spp,
    hasCoordinate = TRUE,
    limit = 5000
  )
  
  if (!is.null(occ$data) && nrow(occ$data) > 0) {
    
    coords <- occ$data %>%
      dplyr::select(decimalLongitude, decimalLatitude) %>%
      na.omit()
    
    # Convert to spatial points
    pts <- vect(coords, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
    
    # Extract BIOCLIM
    clim_vals <- terra::extract(bioclim, pts)
    
    species_df <- cbind(
      species = spp,
      clim_vals
    )
    
    all_species_data[[spp]] <- species_df
  }
}

#Step 4: Combine + Clean Data
climate_data <- bind_rows(all_species_data)

# Remove ID column from terra extract
climate_data <- climate_data %>% dplyr::select(-ID)

# Remove rows with NA
climate_data <- na.omit(climate_data)

#Step 5: Prepare PCA Matrix
#Separate species label from climate variables:

species_labels <- climate_data$species
clim_matrix <- climate_data %>% dplyr::select(-species)

# Scale variables (VERY IMPORTANT)
clim_scaled <- scale(clim_matrix)

#Step 6: Run PCA

pca_res <- prcomp(clim_scaled, center = TRUE, scale. = TRUE)
summary(pca_res)

#Step 7: Visualize PCA

pca_scores <- as.data.frame(pca_res$x)
pca_scores$species <- species_labels

ggplot(pca_scores, aes(PC1, PC2, color = species)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "PCA of BIOCLIM Variables")

fviz_pca_ind(pca_res,
             habillage = species_labels,
             addEllipses = TRUE,
             ellipse.level = 0.95)
#Optional: Variable Contributions

fviz_pca_var(pca_res)



