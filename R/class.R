#' @export
Simulation <- R6Class("Simulation", # nolint
  public = list(
    crop = NULL,
    location = NULL,
    management = NULL,
    soil = NULL,
    weather = NULL,
    fd = NULL,
    result = NA,
    initialize = function(crop, location, management, soil, weather, hydration_depth = 600) {
      self$crop <- crop
      self$location <- location
      self$management <- management
      self$soil <- soil
      if (!"hydration_depth" %in% names(self$soil)) {
        self$soil$hydration_depth <- min(hydration_depth, self$soil$soil_depth)
      }
      self$weather <- weather
      self$prepare_weather()
      private$pre_calc_data()
    },
    prepare_weather = function() {
      private$data <- prepare_weather(
        self$weather,
        self$crop,
        self$soil$albedo,
        self$location$pchng,
        self$location$tchng
      )
    },
    get_data = function() private$data,
    get_state = function() private$state,
    get_fallow_data = function() private$fallow_data,
    get_sowing_date = function() private$sowing_day,
    run_simulation = function() {
      private$state <- private$get_initial_state()
      state <- private$state
      self$result <- private$data %>%
        filter(year == private$sowing_day$year &
          doy >= private$sowing_day$doy &
          doy <= self$management$StopDoy) %>%
        dplyr::rowwise() %>%
        dplyr::group_map(.f = function(day, key) {
          state$day <- day
          state$WSFDS <- 1
          state$WSFL <- 1
          private$phenology_step(state)
          private$LAI_step(state)
          private$dry_matter_step(state)
          if (self$management$water %in% c(1, 2)) {
            # TODO implement checks on water selection
            private$water_step(state)
          }

          with(state, data.frame(
            irrigation_no = IRGNO,
            current_usable_water_top = current_top,
            fraction_usable_water_top = current_hd,
            total_water_top = total_top,
            current_usable_water_hd = current_hd,
            fraction_usable_water_hd = fraction_hd,
            total_water_hd = total_hd,
            water_below_roots = water_below_roots,
            current_transpirable_water = current_transp,
            maximum_transpirable_water = maximum_transp,
            drain = drain_transp,
            surface_runoff = s_runoff,
            total_runoff = runoff,
            soil_evaporation = evaporation,
            total_transpiration = transpiration,
            top_layer_transpiration = transpiration_top,
            NDS = NDS,
            DAP = DAP,
            HI = HI,
            root_depth = DEPORT,
            LAI = LAI,
            RUE = RUE,
            WVEG = WVEG,
            WTOP = WTOP,
            WGRN = WGRN,
            MAT = MAT
          )) %>% dplyr::bind_cols(day, .)
        }) %>%
        dplyr::bind_rows()
    },
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
    fallow_data = NULL,
    consts = NULL,
    current_step = NA,
    state = NA,
    sowing_day = NULL,
    get_initial_state = function() {
      state <- new.env()
      list2env(private$consts, envir = state)

      list2env(self$crop, envir = state)
      list2env(self$location, envir = state)
      list2env(self$management, envir = state)
      list2env(self$soil, envir = state)

      with(state, {
        DAP <- 0
        NDS <- 0
        CTU <- 0
        MAT <- 0
        current_top <- initial_usable_water_top
        current_hd <- initial_usable_water_hd

        total_top <- lower_limit_top + current_top
        total_hd <- lower_limit_hd + current_hd
        fraction_top <- current_top / maximum_usable_water_top
        fraction_hd <- current_hd / maximum_usable_water_hd

        maximum_transp <- iDEPORT * extractable_water
        current_transp <- maximum_transp * MAI
        fraction_transp <- current_transp / maximum_transp
        water_below_roots <- soil_depth * extractable_water * MAI - current_transp
        DEPORT <- iDEPORT
        IRGNO <- 0
        EWAT <- 0
        # LAI
        BLSLAI <- LAI1 <- LAI2 <- LAI <- MXXLAI <- 0
        ETLAIMN <- 0
        # Dry matter
        seed_index <- 0
        IRUE <- IRUE * CO2RUE
        SUMFINT <- SUMIPAR <- WGRN <- DDMP <- 0
        WVEG <- 1
        WTOP <- WVEG + WGRN
        WSFG <- 1
      })

      state
    },
    phenology_step = function(state) {
      with(state, {
        # influenced by available water, potentially stops de
        if (NDS > frEMR) {
          day$daily_temp_unit <- day$daily_temp_unit * WSFDS
        }

        CTU <- CTU + day$daily_temp_unit
        # normalized development score
        NDS <- CTU / tuHAR
        # days already processed?
        DAP <- DAP + 1

        # Stop if plants are fully grown or maximum time reached
        if (NDS >= 1 || day$doy == StopDoy) MAT <- 1 # probably unnecessary
      })
    },
    LAI_step = function(state) {
      with(state, {
        GLAI <- DLAI <- 0

        if (NDS >= frEMR && NDS < frBLS) {
          # between emergence and effective seed growth (bsg/bls)
          LAI2 <- NDS / (NDS + exp(AL - BL * NDS)) * LAIMX
          GLAI <- (LAI2 - LAI1) * WSFL
          LAI1 <- LAI2
          BLSLAI <- LAI # save value of LAI at BLS
        } else if (NDS >= frBLS) {
          LAI2 <- BLSLAI * ((1.000001 - NDS) / (1 - frBLS))^SRATE
          DLAI <- (LAI - LAI2) * WSFDS
        }

        # frost & heat
        DLAIF <- 0
        if (NDS > frEMR && day$t_min < FrzTh) {
          frstf <- abs(day$t_min - FrzTh) * FrzLDR
          if (frstf < 0) frstf <- 0
          if (frstf > 1) frstf <- 1
          DLAIF <- LAI * frstf
        }
        DLAI <- max(DLAI, DLAIF)

        DLAIH <- DLAI
        if (NDS > frEMR && day$t_max > HeatTH) {
          heatf <- max(1 + (day$t_max - HeatTH) * HtLDR, 1)
          DLAIH <- DLAI * heatf
        }
        DLAI <- max(DLAI, DLAIH)

        LAI <- max(LAI + GLAI - DLAI, 0)
      })
    },
    dry_matter_step = function(state) {
      with(state, {
        RUE <- IRUE * calculate_rue_factor(day$average_temp, self$crop) * WSFG
        if (NDS < frEMR || NDS > frPM) RUE <- 0

        FINT <- 1 - exp(-KPAR * LAI)

        DDMP <- day$srad * 0.48 * FINT * RUE

        HI <- WGRN / WTOP
        TRANSL <- SGR <- 0

        if (NDS < frBSG) {
          TRLDM <- WTOP * FRTRL
        } else if (NDS >= frBSG && NDS <= frTSG) {
          DHI <- PDHI * day$daily_temp_unit
          SGR <- DHI * (WTOP + DDMP) + DDMP * HI

          if (HI >= HImax) SGR <- 0

          seed_index <- (SGR / GCC)

          if (seed_index > DDMP) {
            TRANSL <- min(seed_index - DDMP, TRLDM)
          } else if (seed_index <= DDMP) {
            TRANSL <- 0
          }

          TRLDM <- TRLDM - TRANSL

          if (SGR > (DDMP + TRANSL) * GCC) SGR <- (DDMP + TRANSL) * GCC
        }

        WGRN <- WGRN + SGR
        WVEG <- WVEG + DDMP - seed_index
        WTOP <- WVEG + WGRN
      })
    },
    water_step = function(state) {
      with(state, {
        if (NDS <= frBLS) {
          et_LAI <- LAI
        } else {
          et_LAI <- BLSLAI
        }
        et_LAI <- min(et_LAI, ETLAIMN)

        # Drain
        drain_top <- drain(current_top, maximum_usable_water_top, drain_factor)
        drain_hd <- drain(current_hd, maximum_usable_water_hd, drain_factor)
        drain_transp <- drain(current_transp, maximum_transp, drain_factor)
        water_below_roots <- max(water_below_roots + drain_transp - EWAT, 0)

        # Automatic Irrigation
        IRGW <- 0
        if (water == 1 && fraction_transp <= IRGLVL && NDS > 0 && NDS < (0.975 * frPM)) {
          IRGW <- maximum_transp - current_transp
          IRGNO <- IRGNO + 1
        }
        day$irrigation_mm <- IRGW
        # TODO Irrigation by fixed days
        
        # EWAT - Water exploitation by root growth
        GRTD <- GRTDP * day$daily_temp_unit
        if (NDS < frBRG ||
          NDS > frTRG ||
          DDMP == 0 ||
          DEPORT >= soil_depth ||
          DEPORT >= MEED ||
          water_below_roots == 0) {
          GRTD <- 0
        }
        DEPORT <- DEPORT + GRTD
        EWAT <- min(GRTD * extractable_water, water_below_roots)

        # Runoff
        s_runoff <- surface_runoff(current_hd, day$rain_mm, rain_fed = water == 2)
        runoff <- s_runoff + depth_runoff(total_hd, drain_hd, s_runoff, saturation_hd, 0, sat_drain_factor, slope)
        evaporation <- calculate_soil_evaporation(
          day$potential_et, day$days_since_wetting,
          fraction_transp, current_top, et_LAI
        )

        # Transpiration
        transpiration <- total_transpiration(day$t_min, day$t_max, DDMP, VPDF, transpiration_efficiency)
        transpiration_top <- top_transpiration(transpiration, DEPORT, fraction_top, depth_top_layer, WSSG)

        # Updating
        fixed_change <- day$rain_mm + day$irrigation_mm - runoff - evaporation

        current_top <- current_top - drain_top - transpiration_top + fixed_change
        current_top <- max(current_top, 0)
        fraction_top <- current_top / maximum_usable_water_top
        total_top <- lower_limit_top + current_top

        current_transp <- current_transp + EWAT - drain_transp - transpiration + fixed_change
        current_transp <- max(current_transp, 0)
        maximum_transp <- DEPORT * extractable_water
        fraction_transp <- current_transp / maximum_transp
        WATRL <- DEPORT * lower_limit + current_transp
        WSATRL <- DEPORT * saturation

        current_hd <- current_hd - drain_hd - transpiration + fixed_change
        current_hd <- max(current_hd, 0)
        fraction_hd <- current_hd / maximum_usable_water_hd
        total_hd <- lower_limit_hd + current_hd

        # Water stress-factors
        WSFL <- min(fraction_transp / WSSL, 1)
        WSFG <- min(fraction_transp / WSSG, 1)
        WSFDS <- (1 - WSFG) * WSSD + 1

        # currently no rice so no minWH
        if (WATRL > (0.99 * WSATRL) && TRUE) {
          WSFN <- WSFG <- WSFL <- 0
        }
      })
    },
    set_consts = function() {
      soil_consts <- with(self$soil, {
        lower_limit_hd <- hydration_depth * lower_limit
        lower_limit_top <- soil_depth * lower_limit
        initial_usable_water_hd <- hydration_depth * extractable_water * self$management$MAI
        initial_usable_water_top <- depth_top_layer * extractable_water * self$management$MAI1
        list(
          maximum_usable_water_top = depth_top_layer * extractable_water,
          maximum_usable_water_hd = hydration_depth * extractable_water,
          initial_usable_water_hd = initial_usable_water_hd,
          initial_usable_water_top = initial_usable_water_top,
          initial_water_below_roots = soil_depth * extractable_water * self$management$MAI,
          initial_water_top = lower_limit_top + initial_usable_water_top,
          lower_limit_hd = lower_limit_hd,
          lower_limit_top = lower_limit_top,
          drained_upper_limit_hd = hydration_depth * drained_upper_limit,
          saturation_hd = hydration_depth * saturation,
          initial_water_hd = lower_limit_hd + initial_usable_water_hd
        )
      })

      dry_matter_consts <- with(self$crop, {
        rue <- function(co2, c3c4) 1 * (1 + c3c4 * (log(co2 / 330) / log(10)))
        c3c4 <- `C3/C4`
        RUE_385 <- rue(385, c3c4)
        RUE_CO2 <- rue(self$location$CO2, c3c4)
        CO2RUE <- RUE_CO2 / RUE_385
        transpiration_efficiency <- TEC * CO2RUE

        list(
          GRTDP = (MEED - iDEPORT) / ((frTRG - frBRG) * tuHAR), # mm per oC
          PDHI = HImax / (tuHAR * (frTSG - frBSG)),
          RUE_385 = RUE_385,
          RUE_CO2 = RUE_CO2,
          CO2RUE = CO2RUE,
          transpiration_efficiency = transpiration_efficiency
        )
      })

      crop_consts <- with(self$crop, {
        PART1 <- log((1 / y1 - 1) / (1 / x1))
        PART2 <- log((1 / y2 - 1) / (1 / x2))
        BL <- (PART2 - PART1) / (x1 - x2)
        AL <- PART1 + BL * x1
        list(
          PART1 = PART1,
          PART2 = PART2,
          BL = BL,
          AL = AL
        )
      })

      private$consts <- c(soil_consts, dry_matter_consts, crop_consts)
    },
    pre_calc_data = function() {
      private$set_consts()
      fallow_state <- private$get_initial_state()

        private$fallow_data <- private$data %>%
          dplyr::rowwise() %>%
          dplyr::group_map(.f = function(day, key) {
            fallow_state$day <- day
            private$water_step(fallow_state)
            with(fallow_state, data.frame(
              current_usable_water_top = current_top,
              fraction_usable_water_top = current_hd,
              total_water_top = total_top,
              current_usable_water_hd = current_hd,
              fraction_usable_water_hd = fraction_hd,
              total_water_hd = total_hd
            ))
          }) %>%
          dplyr::bind_rows() %>%
          dplyr::bind_cols(private$data, .)

      start_doy <- self$management$Fpdoy %||% self$management$SimDoy # TODO improve input names -> via function?
      private$sowing_day <- find_sow_date(
        private$fallow_data,
        start_doy,
        self$management$Fyear,
        self$management$SearchDur,
        self$management$SowTemp,
        self$management$SowWat,
        self$management$FixFind
      )
    }
    # TODO validators
  )
)
