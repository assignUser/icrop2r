pf <- function(...) {
  parent.frame(...)
}

`%||%` <- function(x, y) {
  if (purrr::is_null(x)) {
    y
  } else {
    x
  }
}

n_day_sum <- function(n, values) {
  sapply(1:n - 1, function(x) lag(values, x)) %>% rowSums()
}
