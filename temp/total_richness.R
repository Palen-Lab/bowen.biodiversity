# Bowen Island Boundary
area <- sf::st_read("RawData/bowen_boundary/Bowen_boundary.shp")

# Available SDMs
SDM_dir <- "RawData/2024_11_28_SDM_Bowen_Island/"
SDM_filepaths <- list.files(SDM_dir, recursive = T, full.names = T) 
SDM_filepaths <- SDM_filepaths[grepl("_NAs_", SDM_filepaths)]

# Create species list from available SDMs
species_list <- SDM_filepaths %>%
  basename() %>%
  stringr::str_extract("^\\S+\\.[a-zA-Z]+_") %>%
  stringr::str_remove("_") %>%
  stringr::str_replace("\\.", " ")

type_list <- SDM_filepaths %>%
  dirname() %>%
  str_remove("_mask")


# iNat observations from rinat fork for Palen-Lab
# SDM maximum values anywhere on Bowen Island, generated from Zonation Collaboration 2014 SDMs 
bowen_SDM_iNat <- SDM_iNat_congruence(SDM_dir = SDM_dir, 
                                      species_list = species_list,
                                      area = area,
                                      mask_NA = T)
bowen_SDM_iNat$taxon_group <- type_list

# iNat query can take a while
# write.csv(bowen_SDM_iNat, "ProcessedData/2024_11_29_SDM_CS_all_species")
bowen_SDM_iNat <- read.csv("ProcessedData/2024_11_29_SDM_CS_all_species")

# Mammal presences from Pottinger Gaherty Environmental Consultants, Ltd.(January 2005) 
bowen_mammals <- readxl::read_xlsx("RawData/Appendix 2 - Expected Occurrence of Wildlife Species In The Cape Roger Curtis Area.xlsx", 
                                   sheet = "mammals") %>%
  filter(`Documented on Bowen` == "y" | `Documented on Bowen` == "?" | `Documented on Bowen` == "One or more spp" )

# number of eBird observations by bird species on Bowen Island
bowen_birds <- read.csv("RawData/final_bird.csv") %>%
  select(c("scientific_name", "n_ebird"))

bowen_SDM_iNat_ebird <- merge(bowen_SDM_iNat, 
                              bowen_birds, 
                              by.x = "species", 
                              by.y = "scientific_name",
                              all.x = TRUE) %>%
  rename(n_obs_inat = inat_n_obs,
         n_obs_ebird = n_ebird) %>%
  select(
    species,
    sdm_max_value,
    taxon_group,
    n_obs_inat,
    n_obs_ebird
  )

include_threshold <- 0.5 # threshold for including species as present 
bowen_SDM_obs_focal_df <- bowen_SDM_iNat_ebird %>%
  mutate(include = case_when(
    (n_obs_inat > 0 | n_obs_ebird > 0 ) & sdm_max_value > include_threshold ~ "y",
    species %in% bowen_mammals$`Scientific name` ~ "y"
  )
  )

# Create list of species names to include
bowen_focal_list <- bowen_SDM_obs_focal_df %>%
  filter(include == "y") %>%
  select(species) 
bowen_focal_list <- bowen_focal_list$species %>%
  unlist() %>%
  unname() 

# Get Bowen Island SDMs filepaths for the focal species
# Available SDMs
# SDM_dir <- "RawData/2024_11_28_SDM_Bowen_Island/"
# SDM_filepaths <- list.files(SDM_dir, recursive = T, full.names = T) 
# SDM_filepaths <- SDM_filepaths[grepl("_NAs_", SDM_filepaths)]
# SDM_species <- SDM_filepaths %>%
#   basename() %>%
#   stringr::str_extract("^\\S+\\.[a-zA-Z]+_") %>%
#   stringr::str_remove("_") %>%
#   stringr::str_replace("\\.", " ")

SDM_df <- data.frame(filepath = SDM_filepaths,
                     filename = basename(SDM_filepaths),
                     species = species_list)

SDM_obs_df <- merge(SDM_df, bowen_SDM_obs_focal_df, by = "species")
write.csv(SDM_obs_df, "ProcessedData/2024_12_10_SDM_CS_all_species.csv")

SDM_obs_df_bowen <- SDM_obs_df %>%
  filter(species %in% bowen_focal_list)

# Prepare SDM stack
SDM_stack <- list()
the_crs <- "+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
bowen_shoreline <- terra::vect("RawData/shoreline_dem_smoothed2/shoreline_dem_smoothed2.shp") %>%
  terra::project(the_crs)
weights <- matrix(c(0, 0, 1, 1, 1, 0, 0,
                    0, 1, 1, 2, 1, 1, 0,
                    1, 1, 3, 3, 3, 1, 1,
                    1, 2, 3, 5, 3, 2, 1,
                    1, 1, 3, 3, 3, 1, 1,
                    0, 1, 1, 2, 1, 1, 0,
                    0, 0, 1, 1, 1, 0, 0
), nrow=7)
for(i in 1:length(SDM_obs_df_bowen$filepath)) {
  SDM <- terra::rast(SDM_obs_df_bowen$filepath[i]) %>%
    terra::project(y = the_crs) %>%
    terra::focal(w = weights, 
                 fun = "mean",
                 na.policy = "only") %>%
    terra::mask(bowen_shoreline)
  SDM_stack[[i]] <- SDM
}
SDM_stack <- SDM_stack %>%
  terra::rast()

SDM_stack %>% sum(na.rm = T) %>% plot()

source("R/utils-sdm.R")
total_richness_0.5 <- sdm_to_species_richness(SDM_stack = SDM_stack,
                                          presence_threshold = 0.5)
writeRaster(total_richness_0.5, 
            "ProcessedData/2024_12_10_species_richness_0.5.tif",
            overwrite = TRUE)

total_richness_0.7 <- sdm_to_species_richness(SDM_stack = SDM_stack,
                                          presence_threshold = 0.7)
writeRaster(total_richness_0.7, 
            "ProcessedData/2024_12_10_species_richness_0.7.tif",
            overwrite = TRUE)

source("R/plots.R")
total_richness_0.5_plot <- bowen_map(raster_layer = total_richness_0.5,
                                     title = "Total Species Richness",
                                     subtitle = "Based on >0.5 probability of occurrence in Species Distribution Models. Includes Birds, Small Mammals, Amphibians, Reptiles, Fish.",
                                     caption = "Created December 10 2024",
                                     legend_label = "Number of Species")
ggplot2::ggsave(
  "Figures/2024_12_10_total_richness_0.5.png", 
  total_richness_0.5_plot, 
  device = ragg::agg_png, 
  width = 9, height = 12, units = "in", res = 300
)

total_richness_0.7_plot <- bowen_map(raster_layer = total_richness_0.7,
                                     title = "Total Species Richness",
                                     subtitle = "Based on >0.7 probability of occurrence in Species Distribution Models. Includes Birds, Small Mammals, Amphibians, Reptiles, Fish.",
                                     caption = "Created December 10 2024",
                                     legend_label = "Number of Species")
ggplot2::ggsave(
  "Figures/2024_12_10_total_richness_0.7.png", 
  total_richness_0.7_plot, 
  device = ragg::agg_png, 
  width = 9, height = 12, units = "in", res = 300
)

