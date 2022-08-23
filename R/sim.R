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
                          find_sowing_date = 0) {
  # instead of having a long if else chain use function objects with fixed api
  # find simulation start day
  # -- just use dplyr instead of looping through everything
  stopifnot(start_year %in% data$year, is.data.frame(data))

  # get from input
  dur <- 5 # consecutive days conditions have to be fulfilled for
  
  # if fixed day nothing to do
  if (find_sowing_date == 0) {
    sowing_data <- data %>%
      filter(year == start_year & doy == sowing_doy) %>%
      mutate(sowing_day = TRUE)
    return(sowing_data)
  }

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

  if (find_sowing_date %in% c(1, 2, 3)) {
    if (find_sowing_date == 1) {
      # find a 5 day rain free period starting from sowing_doy in start_year
      sowing_data <- sowing_data %>%
        mutate(sowing_day = rain_sum == 0)
    }

    if (find_sowing_date == 2) {
      # find a 5 day rain free period starting from sowing_doy in start_year
      # and average temp in 5 day period > sowing_temp
      sowing_data <- sowing_data %>% mutate(sowing_day = rain_sum == 0 & sliding_average_temp > sowing_temp)
    }
    # TODO this excludes days where sliding_average_temp = sowing_temp
    if (find_sowing_date == 3) {
      # find a 5 day rain free period starting from sowing_doy in start_year
      # and average temp in 5 day period < sowing_temp
      sowing_data <- sowing_data %>% mutate(sowing_day = rain_sum == 0 & sliding_average_temp < sowing_temp)
    }
  }

  if (find_sowing_date == 4) {
    # sow when FTSW1 >= sow_water
    # do run soil_water()/ water should be >= 1
    sowing_data <- sowing_data %>% mutate(sowing_day = fraction_usable_water_top >= sowing_water)
  }

  if (find_sowing_date == 5) {
    # sow when FTSW1 <= sow_water
    # do run soil_water()/ water should be >= 1
    sowing_data <- sowing_data %>% mutate(sowing_day = fraction_usable_water_top <= sowing_water)
  }

  if (find_sowing_date == 6) {
    # sow when cumsum rainfall > sow_water
    sowing_data <- sowing_data %>% mutate(
      sowing_day = rain_sum > sowing_water
    )
  }

  if (find_sowing_date == 7) {
    # sow when cumsum rainfall > sow_water and sliding_average_temp < sow_temp
    sowing_data <- sowing_data %>% mutate(
      sowing_day = rain_sum > sowing_water & sliding_average_temp < sowing_temp
    )
  }

  if (find_sowing_date == 91) {
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
prepare_weather <- function(weather, crop, albedo, rain_mod = 1, temp_mod = 0) {
    calc_snow <- function(snow_mm, rain_mm, t_max) {
      if (t_max <= 1) {
        snow_mm <- snow_mm + rain_mm
      } else {
        snow_melt <- t_max + rain_mm * 0.4
        snow_melt <- min(snow_melt, snow_mm)
        snow_mm <- snow_mm - snow_melt
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
       snow_melt = snow_mm - lead(snow_mm, default = 0),
      ) %>%
      mutate(
        rain_mm = ifelse(t_max <= 1, 0, rain_mm + snow_melt),
        # snow_mm = ifelse(t_max > 1,
        #   snow_mm - snow_melt, # TODO test not negative
        #   snow_mm
        # ),
        irrigation_mm = 0,
        days_since_wetting = calculate_days_since_wetting(rain_mm, irrigation_mm),
        potential_et = calculate_potential_et(t_min, t_max, srad, albedo, 0)
      )
}
