install.packages("rgbif", "stringi", "geodata")
library(rgbif)
library(dplyr)
library(terra)
library(stringi)
library(geodata)
library(future)
library(future.apply)

#### SET WORKING DIRECTORY##############################

wd <- ("~/Desktop/Projects/Valerianaceae_biogeography")
setwd(wd)

#### ReWrite for multiple species list
wc <- geodata::worldclim_global(global = TRUE, var = "bio", res = 10)
elev <- geodata::worldclim_global(global = TRUE, var = "elev", res = 10)

##########################################################################################
#####################Test for Co-linearity of bioclim variables ##########################
##########################################################################################
#install.packages("raster")
#install.packages("sp")
#install.packages("usdm")
#install.packages("virtualspecies")
#if (!requireNamespace("ENMTools", quietly = TRUE)) install.packages("ENMTools")

# Load libraries
library(ENMTools)
library(raster)
library(sp)
library(virtualspecies)
library(usdm)


##### variables to keep###################################################################
#### from get_LatLong_test_for_colinearity.R scipt #######################################
##########################################################################################
#bio_1 Annual Mean Temp 
#bio_2 Mean Diurnal Range
#bio_3 Isothermality
#bio_4 Temp Seasonality
#bio_8 Mean Temp. Wettest Quarter
#bio_9 Mean Temp. Driest Quarter
#bio_13 Precipitation of Wettest Month 
#bio_14 Precipitation of Driest Month 
#bio_15 Precipitation Seasonality 
#bio_19 Precipitation of Coldest Quarter
##########################################################################################
############ Redo All of These!!!! with ABOVE variables. #################################
##########################################################################################

bio1 <- wc[[1]]     # Mean Annual Temperature (MAT)
bio2 <- wc[[2]]     # Mean Diurnal Range (MDT)
bio3 <- wc[[3]]     # Isothermality (ISO)
bio4 <- wc[[4]]     # Temperature seasonality (TS)
bio5 <- wc[[5]]     # Max temp of warmest month (MaxTWM)
bio6 <- wc[[6]]     # Min temp of coolest month (MinTCM) !!!!! Breaks everything
bio7 <- wc[[7]]     # Temp Annual Range (TAR)
bio8 <- wc[[8]]     # KEEP: 
bio9 <- wc[[9]]     # KEEP:
bio10 <- wc[[10]]   #
bio11 <- wc[[11]]   #
bio12 <- wc[[12]]   # Mean Annual Precipitation
bio13 <- wc[[13]]  # Prec of Wettest Month
bio14 <- wc[[14]]  # KEEP:
bio15 <- wc[[15]]  # KEEP: Prec Seasonality
bio16 <- wc[[16]]  #
bio17 <- wc[[17]]  # 
bio18 <- wc[[18]]  # KEEP:
bio19 <- wc[[19]]  # KEEP:   
elev  <- elev[[1]]  # Elevation

get_species_environment <- function(species, bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9, bio10, bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19, elev) {
  
  message("Processing: ", species)
  
  occ <- try(occ_search(
    scientificName = species,
    hasCoordinate = TRUE,
    limit = 10000
  ), silent = TRUE)
  
  if(inherits(occ, "try-error") || nrow(occ$data) == 0) {
    return(data.frame(
      species = species,
      n_records = 0,
      mean_bio1 = NA,
      mean_bio2 = NA,
      mean_bio3 = NA,
      mean_bio4 = NA,
      mean_bio8 = NA,
      mean_bio9 = NA,
      mean_bio13 = NA,
      mean_bio14 = NA,
      mean_bio15 = NA,
      mean_bio19 = NA,
      mean_elevation = NA
    ))
  }
  
  dat <- occ$data %>%
    filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
    distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)
  
  ### Include non-correlated varables
  if(nrow(dat) == 0) {
    return(data.frame(
      species = species,
      n_records = 0,
      mean_bio1 = NA,
      mean_bio2 = NA,
      mean_bio3 = NA,
      mean_bio4 = NA,
      mean_bio8 = NA,
      mean_bio9 = NA,
      mean_bio13 = NA,
      mean_bio14 = NA,
      mean_bio15 = NA,
      mean_bio19 = NA,
      mean_elevation = NA
    ))
  }
  
  # Convert to spatial points
  pts <- vect(dat[, c("decimalLongitude", "decimalLatitude")],
              geom = c("decimalLongitude", "decimalLatitude"),
              crs = "EPSG:4326")
  
  # Take log10 of bioclim values, except for elevation

  # Extract values
  dat$BIO1      <- terra::extract(bio1, pts)[,2]
  dat$BIO2      <- terra::extract(bio2, pts)[,2]
  dat$BIO3     <- terra::extract(bio3, pts)[,2]
  dat$BIO4     <- terra::extract(bio4, pts)[,2]
  dat$BIO8      <- terra::extract(bio8, pts)[,2]
  dat$BIO9      <- terra::extract(bio9, pts)[,2]
  dat$BIO13   <- terra::extract(bio13, pts)[,2]
  dat$BIO14     <- terra::extract(bio14, pts)[,2]
  dat$BIO15     <- terra::extract(bio15, pts)[,2]
  dat$BIO19    <- terra::extract(bio19, pts)[,2]
  dat$Elev    <- terra::extract(elev, pts)[,2]
  
  data.frame(
    species        = species,
    n_records      = nrow(dat),
    mean_bio1       = mean(dat$BIO1,  na.rm = TRUE), #bio1
    mean_bio2       = mean(dat$BIO2,  na.rm = TRUE), #bio1
    mean_bio3       = mean(dat$BIO3,  na.rm = TRUE), #bio2
    mean_bio4       = mean(dat$BIO4,  na.rm = TRUE), #bio2
    mean_bio8       = mean(dat$BIO8,  na.rm = TRUE), #bio3
    mean_bio9       = mean(dat$BIO9,   na.rm = TRUE), #bio4
    mean_bio13      = mean(dat$BIO13,   na.rm = TRUE), #bio5
    mean_bio14       = mean(dat$BIO14,   na.rm = TRUE), #bio5
    mean_bio15       = mean(dat$BIO15,  na.rm = TRUE), #bio12
    mean_bio19       = mean(dat$BIO19,  na.rm = TRUE), #bio12
    mean_elevation = mean(dat$Elev, na.rm = TRUE)  #elevation
  )
}



### Test Species List

species_list <- c("Valeriana hebecarpa","Valeriana leucocarpa","Valeriana grandifolia","Valeriana virescens","Valeriana laxiflora","Valeriana lepidota","Valeriana bracteosa","Valeriana bridgesii","Valeriana papilla","Valeriana vaga","Valeriana verticillata","Valeriana polemoniifolia","Valeriana boelckei","Valeriana carnosa","Valeriana fonkii","Valeriana clarionifolia","Plectritis macrocera","Plectritis congesta","Valeriana moyanoi","Valeriana chilensis","Valeriana chaerophylloides","Valeriana philippiana","Valeriana polystachya","Valeriana graciliceps","Valeriana radicalis","Valeriana adscendens","Valeriana secunda","Valeriana plantaginea","Valeriana pilosa","Valeriana tatamana","Valeriana hirtella","Valeriana microphylla","Valeriana connata","Valeriana lapathifolia","Valeriana nivalis","Valeriana arborea","Valeriana tomentosa","Valeriana hornschuchiana","Valeriana stricta","Valeriana crispa","Valeriana clematitis","Valeriana naidae","Valeriana pyramidalis","Valeriana effusa","Valeriana interrupta","Valeriana gallinae","Valeriana rzedowskiorum","Valeriana densiflora","Valeriana sorbifolia","Valeriana selerorum","Valeriana mexicana","Valeriana urticifolia","Valeriana candolleana","Valeriana scandens","Valeriana apiifolia","Valeriana tanacetifolia","Valeriana bryophila","Valeriana robertianifolia","Valeriana palmeri","Valeriana albonervata","Valeriana ceratophylla","Valeriana texana","Valeriana edulis","Valeriana prionophylla","Valeriana stuckertii","Valeriana polybotrya","Valeriana macrorhiza","Valeriana niphobia","Valeriana aretioides","Valeriana bracteata","Valeriana rumicoides","Valeriana coarctata","Valeriana stenoptera","Valeriana officinalis","Valeriana minutiflora","Valeriana kawakamii", "Valeriana flaccidissima","Valeriana trichostoma","Valeriana sichuanica","Valeriana fauriei","Valeriana collina","Valeriana tuberosa","Valeriana asarifolia","Valeriana phu","Valeriana occidentalis","Valeriana californica","Valeriana sitchensis","Valeriana pauciflora","Valeriana acutiloba","Valeriana dioica","Valeriana arizonica","Centranthus trinervis","Centranthus sieberi","Centranthus ruber","Centranthus lecoqii","Centranthus angustifolius","Centranthus longiflorus","Centranthus macrosiphon","Centranthus calcitrapae","Valeriana tripteris","Valeriana supina","Valeriana montana","Valeriana pyrenaica","Valerianella stenocarpa","Valerianella texana","Valerianella radiata","Valerianella amarella","Valerianella dentata","Valerianella eriocarpa","Valerianella vesicaria","Valerianella coronata","Valerianella discoidea","Valerianella carinata","Valerianella locusta","Valerianella pumila","Patrinia saniculifolia","Fedia cornucopiae","Valeriana celtica","Valeriana saxatilis","Valeriana elongata","Nardostachys jatamansi","Patrinia gibbosa","Patrinia villosa","Patrinia triloba","Patrinia scabiosifolia","Triplostegia glandulifera")


### Read CSV file of species -HEADER is species
###species_list <- read.csv("valeriana_species_list.csv")

### multiprocesors
plan(multisession)


results <- do.call(
  rbind,
  lapply(species_list, get_species_environment,
         bio1 = bio1,
         bio2 = bio2,
         bio3 = bio3,
         bio4 = bio4,
         bio8 = bio8,
         bio9 = bio9,
         bio13 = bio13,
         bio14 = bio14,
         bio15 = bio15,
         bio19 = bio19,
         elev = elev)
)

results

###write results to file
write.csv(results, "~/Desktop/Projects/Valerianaceae_biogeography/species_environmental_traits2.csv", row.names = FALSE)


















