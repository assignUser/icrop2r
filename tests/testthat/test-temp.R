test_that("DTU calc works", {
  avg <- c(15.23, -15.23, 0, 17)

  expect_equal(calculate_daily_temp_unit(avg, 10, 20, 15, 16), c(5, 0, 0, 3.75))

})
