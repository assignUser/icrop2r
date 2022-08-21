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

weather <- read_csv_arrow("C:\\Users\\jacob\\Downloads\\SSM-iCrop2\\Weather\\SEM_Garmsar_40758.csv", skip = 9)

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

  tibble(Year = year(window), DOY = yday(window))
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
find_sow_date <- function(find_sowing_Date = 1, weather) {
  # instead of having a long if else chain use function objects with fixed api
  # find simulation start day
  # -- just use dplyr instead of looping through everything

  # if fixed day nothing to do
  if (find_sowing_date == 1) {
    return()
  }

  # get from input
  dur <- 5 # consecutive days conditions have to be fulfilled for
  temp_mod <- 0
  rain_mod <- 1
  sowing_temp <- 15
  sowing_doy <- 100
  sowing_window_duration <- 120
  start_year <- 1992
  sowing_window <- get_sowing_window(start_year, sowing_doy, sowing_window_duration)
  sowing_data <- weather %>%
    filter(Year %in% sowing_window$Year, DOY %in% sowing_window$DOY) %>%
    mutate(
      sliding_average_temp = n_day_sum(dur, average_temp) / dur,
      rain_sum = n_day_sum(dur, RAIN)
    )

  if (find_sowing_date %in% c(2, 3, 4)) {
    if (find_sowing_date == 2) {
      # find a 5 day rain free period starting from sowing_day in start_year
      sowing_day <- sowing_data %>%
        mutate(sowing_day = rain_sum == 0)
    }

    if (find_sowing_date == 3) {
      # find a 5 day rain free period starting from sowing_day in start_year
      # and average temp in 5 day period > sowing_temp
      sowing_day <- sowing_data %>% mutate(sowing_day = rain_sum == 0 & sliding_average_temp > sowing_temp)
    }
    # TODO this excludes days where sliding_average_temp = sowing_temp
    if (find_sowing_date == 4) {
      # find a 5 day rain free period starting from sowing_day in start_year
      # and average temp in 5 day period < sowing_temp
      sowing_day <- sowing_data %>% mutate(sowing_day = rain_sum == 0 & sliding_average_temp < sowing_temp)
    }
  }

  if (scenario$management$find_sowing_date == 5) {
    # sow when FTSW1 >= sow_water
    # do run soil_water()/ water should be >= 1
  }
  if (scenario$management$find_sowing_date == 6) {
    # sow when FTSW1 <= sow_water
    # do run soil_water()/ water should be >= 1
  }
  if (scenario$management$find_sowing_date == 7) {
    # sow when cumsum rainfall > sow_water
    sowing_day <- sowing_Data %>% mutate(
      sowing_day = rain_sum > sowing_water
    )
  }
  if (scenario$management$find_sowing_date == 8) {
    # sow when cumsum rainfall > sow_water and sliding_average_temp < sow_temp
    sowing_day <- sowing_Data %>% mutate(
      sowing_day = rain_sum > sowing_water & sliding_average_temp < sowing_temp
    )
  }

  if (find_sowing_date == 9) {
    forcreq <- 100
    temp_min <- 10
    # trees always start at 1.1.
    # bud burst is based on temp accumulation from 1.1
    # temp_mod applied
    sowing_day <- weather %>%
      filter(Year == start_year) %>%
      mutate(
        sforc = cumsum(pmax(average_temp - temp_min, 0)), # using pmax is important here
        sowing_day = sforc > forcreq
      )
  }

  sowing_day <- sowing_day %>%
    filter(sowing_day) %>%
    head(1)
  if (nrow(sowing_day) == 0) stop("Could not find sowing date.")
  stopifnot(nrow(sowing_day) == 1)
  sowing_day
}

#' Update Weather
#'
#' Updates the weather with modifiers and calculates snow and snow melt.
#' @param sim_env Simulation environment created with [create_sim_env()].
#' @param temp_mod Additive modifier for temperature.
#' @param rain_mod Multiplicative modifier for precipitation.
#' @return Updated weather data.
#' @export
update_weather <- function(sim_env, temp_mod = 0, rain_mod = 1) {
  ensure_var(sim_env, "temp_mod", temp_mod, !missing(temp_mod))
  ensure_var(sim_env, "rain_mod", rain_mod, !missing(rain_mod))
  with(sim_env, {
    

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

    data <- data %>%
      mutate(
        t_min = t_min + temp_mod,
        t_max = t_max + temp_mod,
        rain_mm = rain_mm * rain_mod,
        average_temp = (t_max + t_min) / 2,
        daily_temp_unit = calculate_daily_temp_unit(
          average_temp,
          temp_min, temp_max,
          ideal_temp_min, ideal_temp_max
        ),
        # TODO where does the formula & 0.4 magic number come from
        snow_melt = ifelse(t_max > 1, t_max + rain_mm * 0.4, 0),
        snow_mm = purrr::accumulate2(rain_mm, t_max, calc_snow, .init = 0) %>%
          unlist() %>%
          tail(-1)
      ) %>%
      mutate(
        rain_mm = ifelse(t_max <= 1, 0, rain_mm),
        snow_mm = purrr::accumulate2(snow_melt, t_max, calc_snow_melt, .init = 0) %>%
          unlist() %>%
          tail(-1),
        irrigation_mm = 0,
        days_since_wetting = calculate_days_since_wetting(rain_mm, irrigation_mm)
      )
  })
}

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
