#' @keywords internal
#' @importFrom dplyr filter select mutate lag lead
#' @importFrom R6 R6Class
#' @importFrom utils head tail
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

# TODO document
# Currente limitations
# No forages/trees
# Only sims one season at a time.
utils::globalVariables(c(
  "%>%", ".", "AL", "BL", "LAI", "LAIMX", "NDS", "SRATE", "VPDF", "WSFDS", "WSFL", "WSSG", "average_temp",
  "crop", "doy", "frBLS", "frEMR", "fraction_usable_water_top", "irrigation_mm", "rain_mm",
  "rain_sum", "sforc", "sliding_average_temp", "snow_melt", "snow_mm", "sowing_day", "srad",
  "t_max", "t_min", "values", "water", "year", "sat_drain_factor", "saturation_hd", "slope"
))
