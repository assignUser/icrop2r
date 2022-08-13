# Find sowing date
test_that("sowing window intra-year is correct", {
  duration <- 150
  start_year <- 1992
  sowing_doy <- 120

  window <- get_sowing_window(start_year, sowing_doy, duration)
  expect_length(window$DOY, duration)
  expect_equal(window$DOY[[1]], sowing_doy)
  expect_equal(window$DOY[[duration]], sowing_doy + (duration - 1))
  expect_equal(window$DOY[[1]], sowing_doy)
})

test_that("sowing window inter-year is correct", {
  duration <- 150
  start_year <- 1992 # leap year
  sowing_doy <- 300

  window <- get_sowing_window(start_year, sowing_doy, duration)
  expect_length(window$DOY, duration)
  expect_equal(window$DOY[[1]], sowing_doy)
  expect_equal(window$DOY[[duration]], (sowing_doy + (duration - 1)) - 366)
})
