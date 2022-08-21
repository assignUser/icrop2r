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

create_sim_env <- function(crop, soil, weather, management = NULL, .env = parent.frame()) {
  stopifnot(is.list(crop), is.list(soil), is.data.frame(weather))
  sim_env <- list2env(crop, parent = .env)
  sim_env <- list2env(soil, envir = sim_env)
  sim_env$data <- weather
  sim_env
}

ensure_var <- function(env, var, value, force = FALSE) {
  if (!exists(var, where = env) || force) {
    assign(var, value, envir = env)
  }
}
