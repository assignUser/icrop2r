
# function generator with internal state
crop_LAI_gen <- function(sim_env) {
  sub_env <- new.env(parent = sim_env)
  # daily INcrease in LAI
  sub_env$GLAI <- 0
  # daily decrease leaf area increase
  sub_env$DLAI <- 0

  crop_lai <- function() {
    if (NDS >= frEMR && NDS < frBLS) {
      # between emergence and effective seed growth (bsg/bls)
      LAI2 <- NDS / (NDS + exp(AL - BL * NDS)) * LAIMX
      GLAI <- (LAI2 - LAI1) * WSFL
      LAI1 <- LAI2
      BLSLAI <- LAI # save value of LAI at BLS (there is abetter wy)
    } else if (NDS >= frBLS) {
      LAI2 <- BLSLAI * ((1.000001 - NDS) / (1 - frBLS))^SRATE
      DLAI <- (LAI - LAI2) * WSFDS
    }
    # frost n heat
  }
  environment(crop_lai) <- sub_env
  return(crop_lai)
}


total_transpiration <- function(t_min, t_max, daily_dry_matter,
                                vapour_pressure_coefficient = VPDF,
                                transpiration_efficiency = transpiration_efficiency) {
  stopifnot(length(t_min) == length(t_max))
  vapor_pressure <- function(temp) 0.6108 * exp(17.27 * temp / (temp + 237.3))

  vapor_pressure_deficit <- vapour_pressure_coefficient * (vapor_pressure(t_max) - vapor_pressure(t_min))

  transpiration <- daily_dry_matter * vapor_pressure_deficit / transpiration_efficiency
  transpiration <- max(transpiration, 0)
}

top_transpiration <- function(transpiration,
                              root_depth,
                              fraction_usable_water_top,
                              depth_top_layer = get(depth_top_layer, pf()),
                              transpiration_threshold = get(WSSG, pf())) {
  stopifnot(length(transpiration) == length(fraction_usable_water_top))
  top_share <- ifelse(fraction_usable_water_top > transpiration_threshold,
    1,
    fraction_usable_water_top / transpiration_threshold
  )

  ifelse(root_depth <= depth_top_layer,
    transpiration,
    transpiration * top_share
  )
}
