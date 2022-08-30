#' Calculate Daily Temperature Unit
#'
#' The daily temperature unit is used to modify growth speed and efficiency.
#' All temperature are in Celsius.
#' @param average_temp Average daily temperature or vector of temperatures.
#' @param temp_min Minimum viable growing temperature.
#' @param temp_max Maximum viable growing temperature.
#' @param ideal_temp_min Lower bound of ideal growing window.
#' @param ideal_temp_max Upper bound of ideal growing window.
#' @return Daily temperature unit.
calculate_daily_temp_unit <- function(average_temp,
                                   temp_min, temp_max,
                                   ideal_temp_min, ideal_temp_max) {
  stopifnot(
    is.numeric(average_temp),
    is.numeric(temp_min), length(temp_min) == 1,
    is.numeric(temp_max), length(temp_max) == 1,
    is.numeric(ideal_temp_min), length(ideal_temp_min) == 1,
    is.numeric(ideal_temp_max), length(ideal_temp_max) == 1
  )

  temp_factor <- dplyr::case_when(
    average_temp <= temp_min | average_temp >= temp_max ~ 0,
    average_temp > temp_min & average_temp < ideal_temp_min ~ (average_temp - temp_min) / (ideal_temp_min - temp_min),
    average_temp > ideal_temp_max & average_temp < temp_max ~ (temp_max - average_temp) / (temp_max - ideal_temp_max),
    average_temp >= ideal_temp_min & average_temp <= ideal_temp_max ~ 1
  )

  (ideal_temp_min - temp_min) * temp_factor
}


calculate_rue_factor <- function(average_temp, crop) {
    dplyr::case_when(
    average_temp <= crop$TBRUE | average_temp >= crop$TCRUE ~ 0,
    average_temp > crop$TBRUE & average_temp < crop$TP1RUE ~ (average_temp - crop$TBRUE) / (crop$TP1RUE - crop$TBRUE),
    average_temp > crop$TP2RUE & average_temp < crop$TCRUE ~ (crop$TCRUE - average_temp) / (crop$TCRUE - crop$TP2RUE),
    average_temp >= crop$TP1RUE & average_temp <= crop$TP2RUE ~ 1
  )
}
