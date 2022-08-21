# DUL = Feldkapazität: GPV - LK
# EXTR = nFK: FK - TOT/LL
# LL = lower limit -> tot wasser: DUL - EXTR
# SAT = GPV

#' convert gesättigte Wasserleitfähigkeit to CN runoff curve number
#' @param kf gesättigte Wasserleitfähigkeit in cm/d
#'
kf2CN <- function(kf) {
  # convert to mm/h
  kf <- kf * 10 / 24

  # based on the SCS hydrological soil groups
  dplyr::case_when(
    kf > 7.6 ~ 77, # A
    kf > 3.8 & kf <= 7.6 ~ 86, # B
    kf > 1.3 & kf <= 3.8 ~ 91, # C
    kf <= 1.3 ~ 94 # D
  )
}

#' Compute iCrop2 compatible soil values from german soil info
#'
#' Units mm/mm so divide by 100 for mm/dm
#' @param drain_factor Percentage of water over ATSW that drains per day. If not set
#'   a crude linear model will be used. Generally related to soil composition:
#'   more clay = smaller df (~0.2 for predominantly clay), more sand = larger df (~0.75
#'   for predom. sand). See Sinclair & Soltani (2012) p. 177, table 14.2
#' @param we Effektive Wurzelungstiefe in cm.
german_soil <- function(name, GPV, LK, nFK, kf, we,
                        depth_top_layer = 200,
                        albedo = 0.1,
                        drain_factor = NULL,
                        sat_drain_factor = 0.5,
                        slope = 10,
                        hydration_depth = 600) {
  # this is only a rough approximation
  cn <- kf2CN(kf)
  if (is.null(drain_factor)) {
    drain_factor <- round(cn * -0.0291 + 2.6139, 2)
  }
  stopifnot(GPV >= LK + nFK, is.character(name))
  list(
    name = name,
    soil_depth = we * 10, # maximum depth for roots to grow to in mm
    depth_top_layer = depth_top_layer,
    albedo = albedo,
    CN = cn,
    drain_factor = drain_factor,
    saturation = GPV,
    drained_upper_limit = GPV - LK,
    extractable_water = nFK,
    lower_limit = GPV - LK - nFK,
    sat_drain_factor = sat_drain_factor,
    slope = slope,
    hydration_depth = hydration_depth
  )
}

soil_water <- function(sim_env) { # NO
  with(sim_env, {
    data <- data %>% mutate(
      irrigation_mm = 0, # TODO make conditional/setup structure in the beginning
      et_LAI = 0, # TODO make conditional/setup structure in the beginning
      days_since_wetting = calculate_days_since_wetting(rain_mm, irrigation_mm),
      potential_et = calculate_potential_et(t_min, t_max, srad, albedo, et_LAI),
      # TODO water
      soil_evaporation = calculate_soil_evaporation(
        potential_et,
        days_since_wetting,
        fraction_transpirable_water, current_usable_water_top,
        et_LAI
      )
    )
  })
}


drain <- function(current, maximum, drain_factor = get("drain_factor", pf())) {
  ifelse(current > maximum,
    (current - maximum) * drain_factor,
    0
  )
  # water_below_roots <- max(water_below_roots + drain_transpirable - root_growth_water(), 0)
}


calculate_fallow_water <- function(sim_env, hydration_depth = 600, MAI = 0.9, root_depth = 0) {
  ensure_var(sim_env, "hydration_depth", hydration_depth)
  ensure_var(sim_env, "MAI", MAI)

  with(sim_env, {
    # data <- data %>% mutate(
    #   current_usable_water_top = initial_usable_water_top,
    #   current_usable_water_hd = initial_usable_water_hd,
    #   total_water_hd = initial_water_hd,
    #   current_water_below_roots = initial_water_below_roots,
    #   current_transpirable_water = 0,
    #   irrigation_mm = 0
    # )
    days <- nrow(data)
    current_usable_water_top <- current_usable_water_hd <- double(days)
    fraction_usable_water_top <- fraction_usable_water_hd <- double(days)
    total_water_top <- total_water_hd <- double(days)
    current_top <- initial_usable_water_top
    current_hd <- initial_usable_water_hd
    drain_top <- drain_hd <- double(days)
    transpiration <- transpiration_top <- 0
    evaporation <- double(days)

    for (day in seq_len(days)) {
      total_top <- lower_limit_top + current_top
      total_hd <- lower_limit_hd + current_hd
      drain_top[[day]] <- drain(current_top, maximum_usable_water_top)
      drain_hd[[day]] <- drain(current_hd, maximum_usable_water_hd)
      # drain_transp
      rain <- data$rain_mm[[day]]
      irrigation <- data$irrigation_mm[[day]]
      s_runoff <- surface_runoff(current_hd, rain, )
      runoff <- s_runoff + depth_runoff(total_hd, drain_hd[[day]], s_runoff)
      evaporation[[day]] <- calculate_soil_evaporation(
        data$potential_et[[day]], data$days_since_wetting[[day]],
        0, current_top, 0, 0
      )
      fixed_change <- rain + irrigation - runoff - evaporation[[day]]

      current_top <- current_top - drain_top[[day]] - transpiration_top + fixed_change
      current_top <- max(current_top, 0)
      current_usable_water_top[[day]] <- current_top
      fraction_usable_water_top[[day]] <- current_top / maximum_usable_water_top
      total_water_top[[day]] <- lower_limit_top + current_top

      current_hd <- current_hd - drain_hd[[day]] - transpiration + fixed_change
      current_hd <- max(current_hd, 0)
      fraction_usable_water_hd[[day]] <- current_hd / maximum_usable_water_hd
      total_water_hd[[day]] <- lower_limit_hd + current_hd
    }

    data <- data %>% mutate(
      current_usable_water_top,
      fraction_usable_water_top,
      total_water_top,
      current_usable_water_hd,
      fraction_usable_water_hd,
      total_water_hd
    )
  })

  # maximum_transpirable_water = root_depth * soil$extractable_water,
  # current_transpirable_water = root_depth * soil$extractable_water
}

root_growth_water <- function(daily_temp_unit, # TODO
                              potential_root_growth = crop$GRTDP,
                              frBRG = crop$frBRG,
                              frTRG = crop$frTRG,
                              MEED = crop$MEED,
                              NDS,
                              dry_matter_production,
                              root_depth,
                              extractable_water = extractable_water,
                              soil_depth = soil_depth,
                              water_below_roots = water$water_below_roots) {
  actual_root_growth <- potential_root_growth * daily_temp_unit
  if (NDS < frBRG ||
    NDS > frTRG ||
    dry_matter_production == 0 ||
    root_depth >= soil_depth ||
    root_depth >= MEED ||
    water_below_roots == 0) {
    actual_root_growth <- 0
  }
  root_depth <- root_depth + actual_root_growth
  root_water <- max(actual_root_growth * extractable_water, water_below_roots)
}

surface_runoff <- function(current_usable_water_hd, rain_mm,
                           CN = get("CN", pf()),
                           maximum_usable_water_hd = get("maximum_usable_water_hd", pf()),
                           et_LAI = get("et_LAI", pf()),
                           slope = get("slope", pf()), rain_fed = TRUE) {
  if (!rain_fed) {
    return(0)
  }
  if (missing(et_LAI)) et_LAI <- 0
  KET <- 0.5
  runoff <- 0
  # Surface runoff. Only for rain fed -> assumption irrigation is well managed so no "waste"
  # TODO origin of formulae?
  CN2 <- CN * exp(0.00673 * (100 - CN))
  CNS <- 0.333 * (CN2 - CN) * (1 - 2 * exp(-13.86 * slope)) + CN
  cover <- (1 - exp(-KET * et_LAI)) * 100
  CNC <- CNS - cover * 0.25
  if ((CNS - CNC) > 20) CNS <- CNS - 20
  CN_ <- CNC - (20 * (100 - CNC)) / (100 - CNC + exp(2.533 - 0.0636 * (100 - CNC)))
  S_max <- 254 * (100 / CN_ - 1)
  S <- S_max * (1 - current_usable_water_hd / (1.12 * maximum_usable_water_hd))

  if (rain_mm > 0.2 * S) {
    runoff <- ((rain_mm - 0.2 * S)^2) / (rain_mm + 0.8 * S)
  }
  runoff
}

depth_runoff <- function(total_water_hd, drain_hd, surface_runoff = 0, saturation_hd = get("saturation_hd", pf()),
                         min_water_height = 0, sat_drain_factor = get("sat_drain_factor", pf()),
                         slope = get("slope", pf()),
                         KET = 0.5) {
  # From saturated soil under rain fed and irrigated land (excl. rice)
  runoff_hd <- 0
  if ((total_water_hd - drain_hd - surface_runoff) > saturation_hd && min_water_height == 0) {
    runoff_hd <- max((total_water_hd - saturation_hd - drain_hd - runoff) * sat_drain_factor, 0)
  }

  runoff_hd
}

calculate_potential_et <- function(t_min, t_max, srad, albedo, et_LAI, chlorophyll_albedo = 0.23, KET = 0.5) {
  weighted_temp <- 0.6 * t_max + 0.4 * t_min
  total_albedo <- chlorophyll_albedo * (1 - exp(-KET * et_LAI)) + albedo * exp(-KET * et_LAI)
  equilibrium_evaporation <- srad * (0.004876 - 0.004374 * total_albedo) * (weighted_temp + 29)

  et_mod <- dplyr::case_when(
    t_max > 34 ~ ((t_max - 34) * 0.05 + 1.1), # increased to account for advection
    t_max < 5 ~ 0.01 * exp(0.18 * (t_max + 20)), # reduced due to stoma closure and frozen ground
    TRUE ~ 1.1 # + 10% to account for effect of unsaturated air
  )

  equilibrium_evaporation * et_mod
}

calculate_soil_evaporation <- function(potential_et,
                                       days_since_wetting,
                                       fraction_transpirable_water,
                                       current_usable_water_top,
                                       et_LAI, minimum_evaporation = 1.5) {
  stopifnot(length(potential_et) == length(days_since_wetting))
  KET <- 0.5
  potential_evaporation <- potential_et * exp(-KET * et_LAI)
  # even with full crop cover there is a minimum evaporation from soil happening
  potential_evaporation <- ifelse(
    potential_et > minimum_evaporation & potential_evaporation < minimum_evaporation,
    minimum_evaporation,
    potential_evaporation
  )

  # enter stage I when soil not completely wet
  ifelse(
    days_since_wetting > 1 |
      fraction_transpirable_water < 0.5 |
      current_usable_water_top <= 1,
    potential_evaporation * ((days_since_wetting + 1)^0.5 - days_since_wetting^0.5),
    potential_evaporation
  )
}

calculate_days_since_wetting <- function(rain_mm, irrigation_mm, WET_MM = 10) {
  ifelse(rain_mm + irrigation_mm > WET_MM, 1, 0) %>%
    rle() %>%
    {
      purrr::map2(
        .$values, .$lengths,
        ~ dplyr::case_when(.x == 1 ~ rep(1, .y), .x == 0 ~ (1:.y) + 1)
      )
    } %>%
    unlist()
}
