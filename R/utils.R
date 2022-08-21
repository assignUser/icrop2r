#' Calculate solar radiation from sunshine hours
#'
#' @param lat Latitude in degrees. South latitudes should be entered
#'   as negative.
#' @param
#' @return
#' @details
#' @examples
#' @exported
#' @md
srad_from_sunh <- function() {

}
pf <- parent.frame
create_sim_env <- function(crop, soil, weather, management = NULL, .env = parent.frame()) {
  stopifnot(is.list(crop), is.list(soil), is.data.frame(weather))
  sim_env <- list2env(crop, parent = .env)
  sim_env <- list2env(soil, envir = sim_env)
  sim_env$data <- weather

  with(sim_env, {
    # TODO add management
    MAI <- MAI1 <- 0.9
    # water
    hydration_depth <- min(hydration_depth, soil_depth)
    maximum_transpirable_water <- 0 * extractable_water
    maximum_usable_water_top <- depth_top_layer * extractable_water
    maximum_usable_water_hd <- hydration_depth * extractable_water
    initial_usable_water_hd <- hydration_depth * extractable_water * MAI
    initial_usable_water_top <- depth_top_layer * extractable_water * MAI1
    initial_water_below_roots <- soil_depth * extractable_water * MAI
    lower_limit_hd <- hydration_depth * lower_limit
    lower_limit_top <- soil_depth * lower_limit # TODO move to soil func
    drained_upper_limit_hd <- hydration_depth * drained_upper_limit
    saturation_hd <- hydration_depth * saturation
    initial_water_hd <- lower_limit_hd + initial_usable_water_hd

    # Dry matter production
    CO2 <- 385 # TODO location inputs
    c3c4 <- `C3/C4`
    rue <- function(co2, c3c4) 1 * (1 + c3c4 * (log(co2 / 330) / log(10)))
    RUE_385 <- rue(385, c3c4)
    RUE_CO2 <- rue(co2, c3c4)
    CO2RUE <- RUE_CO2 / RUE_385
    # TODO var naming scheme CONSTs
    transpiration_efficiency <- TEC * CO2RUE
    # Phenology
    DAP <- NDS <- CTU <- DAYT <- SRAINT <- STMINT <-
      STMAXT <- SSRADT <- SUMETT <- DAY3 <- SRAIN3 <- STMIN3 <-
      STMAX3 <- SSRAD3 <- SUMET3 <- DAY2 <- SRAIN2 <- STMIN2 <-
      STMAX2 <- SSRAD2 <- SUMET2 <- 0
    # water stress factor development seneacence
    WSFDS <- 1

    # Crop LAI
    PART1 <- log((1 / y1 - 1) / (1 / x1))
    PART2 <- log((1 / y2 - 1) / (1 / x2))
    BL <- (PART2 - PART1) / (x1 - x2)
    AL <- PART1 + BL * x1
    # TODO create CONSTs and setup dataframe with all cols
  })

  sim_env
}

ensure_var <- function(env, var, value, force = FALSE) {
  if (!exists(var, where = env) || force) {
    assign(var, value, envir = env)
  }
}
