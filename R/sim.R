scenario <- list(
  name = character(),
  location = list(
    name = character(),
    latitude = double(),
    vapour_pressure_deficit_coefficent = 0.75,
    weather = list(
      year = 1984,
      day_of_year = 1,
      solar_radiation = double(), # in MJ/m^2/day
      temp_max = double(), # celsius obv.
      temp_min = double(),
      rain_mm = double()
    ),
    temp_mod = 0,
    rain__mod = 1,
    co2_ppm = 385
  ),
  management = list(
    name = character(),
    find_sowing_date = factor(1:9,
      labels = c(
        "Fixed sowing date",
        "Sow in the 5th day of a 5-day rainfree period",
        "Sow in the 5th day  of a 5-day rainfree period + average temp > sow_temp.",
        "Sow in the 5th day  of a 5-day rainfree period + average temp < sow_temp.",
        "Sow when FTSW1 > SowWat.",
        "Sow when FTSW1 < SowWat.",
        "Sow when cumulative rainfall over a 5-day period > SowWat.",
        "Sow when cumulative rainfall over a 5-day period > SowWat and average temp < sow_temp.",
        "Prediction of bud burst in tree crops."
      )
    ),
    start_year = integer(),
    sim_n_years = integer(),
    start_doy = integer(),
    sowing_doy = integer(), # fixed sowing date or start of sowing window
    sowing_window_duration = integer(), # days
    rain_free_preiod = 5, # fixed don't change ?
    sowing_temp = integer(),
    sow_water = double(), # fraction transporable soil water
    water = factor(0:3, labels = c(
      "potential prodution",
      "automated irrigation",
      "rain feed",
      "fixed irrigation"
    )),
    irrigation_level = 0.5,
    top_soil_MAI = 0.9, # moisture availability index
    soil_MAI = 0.9,
    latest_harvest_doy = integer(),
    forages_clipping_n = integer(),
    min_water_mm_rice = integer(),
    max_water_mm_rice = integer()
  ),
  soil = list(
    name = character(),
    soil_depth = integer(), # mm
    depth_top_layer = integer(), # mm
    soil_albedo = double(),
    curve_number = integer(), # https://en.wikipedia.org/wiki/Runoff_curve_number
    drainage_factor = double(),
    drained_upper_limit = double(), # mm mm^-1 volumetric soil water content at dul
    extractable_water = double(), # mm mm^-1  Volumetric soil water content available for extraction by crop roots
    liquid_limit = double(), # https://en.wikipedia.org/wiki/Atterberg_limits
    surface_drainage_factor = 0.5, # surface run off?
    slope = 10, # ?
    ec = NULL # ?
  ),
  crop = list(
    tbd = integer(), # Base temp for development
    tp1d = integer(), # Lower optimum temperature for development
    tp2d = integer(), # Upper optimum temperature for development
    tcd = integer(), # Ceiling temperature for development
    forcereg = integer(), # TODO
    tuhar = integer(),
    fremr = double(),
    frbsg = double(),
    frtsg = double(),
    frpm = double(),
    x1 = double(),
    y1 = double(),
    x2 = double(),
    y2 = double(),
    laimx = double(),
    frbls = double(),
    srate = integer(),
    frzth = integer(), # oC
    frzldr = double(), # m2/m2/oC
    heatth = integer(), # oC
    HeatTH = double(),
    HtLDR = double(),
    TBRUE = double(),
    TP1RUE = double(),
    TP2RUE = double(),
    TCRUE = double(),
    KPAR = double(),
    IRUE = double(),
    "C3/C4" = double(),
    HImax = double(),
    FRTRL = double(),
    GCC = double(),
    frBRG = double(),
    frTRG = double(),
    iDEPORT = double(),
    MEED = double(),
    TEC = double(),
    WSSG = double(),
    WSSL = double(),
    WSSD = double(),
    "MC%" = double(),
    SaltTH = double(),
    SaltSlope = double()
  )
)

main_loop <- function(scenarios) {
  for (scenario in scenarios) {
    for (year in seq_len(scenario$management$sim_n_years)) {
      if (start_doy == sowing_doy) start_doy <- start_doy - 1
      MAT <- iniPheno <- iniLAI <- iniSW <- SNOW <- 0

      find_sow_date()
      while (MAT != 1) {
        update_weather()
        phenology(iniPheneo)
        crop_LAI()
      }
    }
  }
}



#' Get sowing window
#'
#' @param year Year sowing window starts in.
#' @param sowing_doy Start of sowing window
#' @param duration Duration of sowing window in days.
#' @return Year and DOY columns of the sowing window.
#' @noRd
get_sowing_window <- function(year, sowing_doy, duration) {
  start_date <- lubridate::make_date(year) %>% lubridate::`yday<-`(sowing_doy)
  end_date <- start_date + (duration - 1)
  window <- start_date:end_date %>% lubridate::as_date()

  dplyr::tibble(year = lubridate::year(window), doy = lubridate::yday(window))
}

n_day_sum <- function(n, values) {
  sapply(1:n - 1, function(x) lag(values, x)) %>% rowSums()
}

nth_rain_free <- function(weather, n) {
  enc <- rle(weather$RAIN)

  ok <- values == 0 & lengths >= 5
  ends <- cumsum(lengths)
  starts <- ends - lengths + 1
  data.frame(starts, ends)[ok, ]
}

find_sow_date <- function(data,
                          sowing_doy,
                          start_year,
                          sowing_window_duration,
                          sowing_temp,
                          sowing_water,
                          find_sowing_date = 1) {
  # instead of having a long if else chain use function objects with fixed api
  # find simulation start day
  # -- just use dplyr instead of looping through everything
  stopifnot(start_year %in% data$year, is.data.frame(data))

  # get from input
  dur <- 5 # consecutive days conditions have to be fulfilled for

  # sowing_temp <- 15
  # sowing_doy <- 100
  # sowing_window_duration <- 120
  # start_year <- 1992
  sowing_window <- get_sowing_window(start_year, sowing_doy, sowing_window_duration)
  sowing_data <- data %>%
    filter(year %in% sowing_window$year & doy %in% sowing_window$doy) %>%
    mutate(
      sliding_average_temp = n_day_sum(5, average_temp) / dur,
      rain_sum = n_day_sum(dur, rain_mm)
    )

  # if fixed day nothing to do
  if (find_sowing_date == 1) {
    sowing_data <- data %>%
      filter(year == start_year & doy == sowing_doy) %>%
      mutate(sowing_day = TRUE)
  }

  if (find_sowing_date %in% c(2, 3, 4)) {
    if (find_sowing_date == 2) {
      # find a 5 day rain free period starting from sowing_doy in start_year
      sowing_data <- sowing_data %>%
        mutate(sowing_day = rain_sum == 0)
    }

    if (find_sowing_date == 3) {
      # find a 5 day rain free period starting from sowing_doy in start_year
      # and average temp in 5 day period > sowing_temp
      sowing_data <- sowing_data %>% mutate(sowing_day = rain_sum == 0 & sliding_average_temp > sowing_temp)
    }
    # TODO this excludes days where sliding_average_temp = sowing_temp
    if (find_sowing_date == 4) {
      # find a 5 day rain free period starting from sowing_doy in start_year
      # and average temp in 5 day period < sowing_temp
      sowing_data <- sowing_data %>% mutate(sowing_day = rain_sum == 0 & sliding_average_temp < sowing_temp)
    }
  }

  if (find_sowing_date == 5) {
    # sow when FTSW1 >= sow_water
    # do run soil_water()/ water should be >= 1
    sowing_data <- sowing_data %>% mutate(sowing_day = fraction_usable_water_top >= sowing_water)
  }

  if (find_sowing_date == 6) {
    # sow when FTSW1 <= sow_water
    # do run soil_water()/ water should be >= 1
    sowing_data <- sowing_data %>% mutate(sowing_day = fraction_usable_water_top <= sowing_water)
  }

  if (find_sowing_date == 7) {
    # sow when cumsum rainfall > sow_water
    sowing_data <- sowing_data %>% mutate(
      sowing_day = rain_sum > sowing_water
    )
  }

  if (find_sowing_date == 8) {
    # sow when cumsum rainfall > sow_water and sliding_average_temp < sow_temp
    sowing_data <- sowing_data %>% mutate(
      sowing_day = rain_sum > sowing_water & sliding_average_temp < sowing_temp
    )
  }

  if (find_sowing_date == 9) {
    forcreq <- 100
    temp_min <- 10
    # trees always start at 1.1.
    # bud burst is based on temp accumulation from 1.1
    # temp_mod applied
    sowing_data <- data %>%
      filter(year == start_year) %>%
      mutate(
        sforc = cumsum(pmax(average_temp - temp_min, 0)), # using pmax is important here
        sowing_day = sforc > forcreq
      )
  }

  sowing_data <- sowing_data %>%
    filter(sowing_day == TRUE) %>%
    head(1)
  if (nrow(sowing_data) == 0) stop("Could not find sowing date.")
  stopifnot(nrow(sowing_data) == 1)
  sowing_data
}

#' Update Weather
#'
#' Updates the weather with modifiers and calculates snow and snow melt.
#' @param sim_env Simulation environment created with [create_sim_env()].
#' @param temp_mod Additive modifier for temperature.
#' @param rain_mod Multiplicative modifier for precipitation.
#' @return Updated weather data.
#' @export
prepare_weather <- function(weather, crop, albedo, rain_mod = 1, temp_mod = 0) {
    calc_snow <- function(snow_mm, rain_mm, t_max) {
      if (t_max <= 1) {
        snow_mm <- snow_mm + rain_mm
      }
      return(snow_mm)
    }

    calc_snow_melt <- function(snow_mm, snow_melt, t_max) {
      if (t_max >= 1) {
        snow_melt <- min(snow_melt, snow_mm)
        snow_mm <- max(snow_mm - snow_melt, 0)
      }
      return(snow_mm)
    }

    weather %>%
      mutate(
        t_min = t_min + temp_mod,
        t_max = t_max + temp_mod,
        rain_mm = rain_mm * rain_mod,
        average_temp = (t_max + t_min) / 2,
        daily_temp_unit = calculate_daily_temp_unit(
          average_temp,
          crop$TBD, crop$TCD,
          crop$TP1D, crop$TP2D
        ),
        # TODO where does the formula & 0.4 magic number come from
        snow_mm = purrr::accumulate2(rain_mm, t_max, calc_snow, .init = 0) %>%
          unlist() %>%
          tail(-1),
        snow_melt = ifelse(t_max > 1, pmin(t_max + rain_mm * 0.4, snow_mm), 0),
      ) %>%
      mutate(
        rain_mm = ifelse(t_max <= 1, 0, rain_mm + snow_melt),
        snow_mm = ifelse(t_max > 1,
          snow_mm - snow_melt, # TODO test not negative
          snow_mm
        ),
        irrigation_mm = 0,
        days_since_wetting = calculate_days_since_wetting(rain_mm, irrigation_mm),
        potential_et = calculate_potential_et(t_min, t_max, srad, albedo, 0)
      )
}

# use with fresh sim_env to pre-calc all values possible & necessary to find sow date
calculate_fallow_data <- function(sim_env) {
  prepare_weather(sim_env)
  calculate_fallow_water(sim_env)
  invisible(sim_env)
}

# Ablauf
# load weather
# load crops
# load/gen soil
# calc_fallow
# find_sowing_date
# simulation
#   phenology
#   lai
#   dry matter
#   water -> no water, rain fed, irrigated (auto/fixed)
