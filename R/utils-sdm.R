#' Create species richness map from SDM rasters
#'
#' @param sdm_stack raster stack, to be combined into a richness map 
#' @param threshold numeric, value from 0 to 1 corresponding to probability 
#' threshold in the cells of the SDM to include into species richness map
#'
#' @return SpatRaster, with species present in each cell
#' @export
#'
sdm_to_species_richness <- function(SDM_stack, 
                                    presence_threshold = 0.7
) {
  presence_stack <- sapp(SDM_stack, fun = function(x) {
    x[x > presence_threshold] <- 1
    x[x <= presence_threshold] <- 0
    return(x)
  })
  
  total_richness <- presence_stack %>% sum(na.rm = T) 
}

# Defining CRS to make sure all rasters have the same
#the_crs <- "+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

# Bowen Shoreline polygon for masking the rasters 
# bowen_shoreline <- terra::vect("RawData/shoreline_dem_smoothed2/shoreline_dem_smoothed2.shp") %>%
#   terra::project(the_crs)
# bowen_boundary <- terra::vect("RawData/bowen_boundary/Bowen_boundary.shp")

#' Prepare rasters to match Bowen Island
#'
#' @param input_rast 
#'
#' @return
#' @export
#'
#' @examples
raster_prep_bowen <- function(input_rast,
                              bowen_boundary,
                              bowen_shoreline,
                              crs) {
  # Weights for terra::focal() function that provides moving window average
  ## - Smoothing algorithm essentially
  ## - Use this to fill some empty cells on Bowen Island with weighted mean of the 
  ##   neighbouring cells. 
  weights <- matrix(c(0, 0, 1, 1, 1, 0, 0,
                      0, 1, 1, 2, 1, 1, 0,
                      1, 1, 3, 3, 3, 1, 1,
                      1, 2, 3, 5, 3, 2, 1,
                      1, 1, 3, 3, 3, 1, 1,
                      0, 1, 1, 2, 1, 1, 0,
                      0, 0, 1, 1, 1, 0, 0
  ), nrow=7)
  
  output_rast <- input_rast %>%
    # Figured out why the SDMs were cropped incorrectly
    # The default for cropping a raster by a polygon is to snap = "near"
    # This means that the polygon extent is snapped to the closest raster
    # We need snap = "out" to make sure the full extent of the polygon is covered by the output raster
    terra::crop(bowen_boundary, snap = "out") %>%
    terra::project(y = crs) %>%
    # Fills the NA values missing in south part of Bowen Island
    terra::focal(w = weights, 
                 fun = "mean",
                 na.policy = "only") %>%
    # Downsamples from 400 m to 100 m resolution
    terra::disagg(fact = 4) %>%
    # Smooths the raster
    terra::focal(w = weights,
                 fun = "mean",
                 na.rm = T,
                 na.policy = "omit") %>%
    terra::mask(bowen_shoreline)
}