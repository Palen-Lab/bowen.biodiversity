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
                                    presence_threshold = 0.7,
                                    crs = "+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
                                    ) {

  presence_stack <- sapp(SDM_stack, fun = function(x) {
    x[x > presence_threshold] <- 1
    x[x <= presence_threshold] <- 0
    return(x)
  })
  
  total_richness <- presence_stack %>% sum(na.rm = T) 
}
