phenology <- function(crop, weather, ini) {
  if (ini == 0) {
    DAP <- NDS <- CTU <- DAYT <- SRAINT <- STMINT <-
      STMAXT <- SSRADT <- SUMETT <- DAY3 <- SRAIN3 <- STMIN3 <-
      STMAX3 <- SSRAD3 <- SUMET3 <- DAY2 <- SRAIN2 <- STMIN2 <-
      STMAX2 <- SSRAD2 <- SUMET2 <- 0
    # water stress factor development seneacence
    water_stress_development <- WSFDS <- 1
    ini <- 1
  }

  daily_temp_unit <- calculate_daily_temp_unit(average_temp, temp_min, temp_max, ideal_temp_min, ideal_temp_max)

  if (NDS > scenario$crop$fremr) {
    daily_temp_unit <- daily_temp_unit * water_stress_development
  }
  # cummulative temp unit
  CTU <- CTU + daily_temp_unit
  # normalized development score
  NDS <- CTU / scenario$crop$tuhar
  # days already processed?
  DAP <- DAP + 1

  # Updating day counters
  if (NDS < scenario$crop$frbsg) days2BSG <- DAP + 1
  if (NDS < scenario$crop$frtsg) days2TSG <- DAP + 1
  if (NDS < 1) days2HAR <- DAP + 1

  # Stop if plants are fully grown or maximum time reached
  if (NDS >= 1 || doy == scenario$management$latest_harvest_doy) MAT <- 1

  # the following bit calculates sums for certain periods during the growth cycle
  # only used to printout in summary -> dplyr
}

crop_LAI <- function() {
  if (iniLAI) {
    PART1 <- log((1 / y1LAI - 1) / (1 / x1NDS))
    PART2 <- log((1 / y2LAI - 1) / (1 / x2NDS))
    BL <- (PART2 - PART1) / (x1NDS - x2NDS)
    AL <- PART1 + BL * x1NDS
    LAI1 <- LAI2 <- LAI <- MXXLAI <- 0
    iniLAI <- FALSE
  }

  # daily INcrease in LAI
  GLAI <- 0
  # daily decrease leaf area increase
  DLAI <- 0

  if (NDS >= fremr && NDS < frbls) {
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
                              depth_top_layer = depth_top_layer,
                              transpiration_threshold = WSSG) {
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
