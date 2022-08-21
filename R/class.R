Simulation <- R6Class("Simulation", # nolint
  public = list(
    crop = NULL,
    location = NULL,
    management = NULL,
    soil = NULL,
    weather = NULL,
    initialize = function(crop, location, management, soil, weather, hydration_depth = 600) {
      self$crop <- crop
      self$location <- location
      self$management <- management
      self$soil <- soil
      if (!"hydration_depth" %in% names(self$soil)) {
        self$soil$hydration_depth <- min(hydration_depth, self$soil$soil_depth)
      }
      self$weather <- weather
      private$pre_calc_data()
    },
    get_data = function() private$data,
    reset = function(sure = FALSE) {
      if (!sure) {
        stop(paste0(
          "This will delete all simulated data!\\n",
          "If you are sure set the arg `sure = TRUE`."
        ))
      }
      private$pre_calc_data()
    }
  ),
  private = list(
    data = NULL,
    consts = NULL,
    current_step = NA,
    sowing_window = NA,
    pre_calc_data = function() {
      private$data <- prepare_weather(
        self$weather,
        self$crop,
        self$soil$albedo,
        self$location$pchng,
        self$location$tchng
      )

      soil_consts <- with(self$soil, list(
        # maximum_transpirable_water = 0 * extractable_water, TODO move else where
        maximum_usable_water_top = depth_top_layer * extractable_water,
        maximum_usable_water_hd = hydration_depth * extractable_water,
        initial_usable_water_hd = hydration_depth * extractable_water * self$management$MAI,
        initial_usable_water_top = depth_top_layer * extractable_water * self$management$MAI1,
        initial_water_below_roots = soil_depth * extractable_water * self$management$MAI,
        lower_limit_hd = hydration_depth * lower_limit,
        lower_limit_top = soil_depth * lower_limit,
        drained_upper_limit_hd = hydration_depth * drained_upper_limit,
        saturation_hd = hydration_depth * saturation,
        initial_water_hd = lower_limit_hd + initial_usable_water_hd,
      ))

      rue <- function(co2, c3c4) 1 * (1 + c3c4 * (log(co2 / 330) / log(10)))
      dry_matter_consts <- with(self$crop, {
        c3c4 <- `C3/C4`
        list(
          RUE_385 = rue(385, c3c4),
          RUE_CO2 = rue(self$location$co2, c3c4),
          CO2RUE = RUE_CO2 / RUE_385,
          transpiration_efficiency = self$crop$TEC * CO2RUE
        )
      })
      crop_consts <- with(self$crop, {
        list(
          PART1 = log((1 / y1 - 1) / (1 / x1)),
          PART2 = log((1 / y2 - 1) / (1 / x2)),
          BL = (PART2 - PART1) / (x1 - x2),
          AL = PART1 + BL * x1
        )
      })
      private$consts <- c(soil_consts, dry_matter_consts, crop_consts)
      # TODo at fallow water
    },

    # TODO validators
  )
)
