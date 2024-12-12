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
