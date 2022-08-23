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
