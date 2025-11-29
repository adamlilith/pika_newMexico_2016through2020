
say('#######################################################################')
say('### analyze patterns of occupancy for multi-year data (Erik Beever) ###')
say('#######################################################################')

	# model terms
	terms <- sort(c('subLethalCold', 'subLethalHeat', 'acuteHeat', 'acuteCold', 'summerRespite', 'chronicCold', 'chronicHeat', 'coldNoSnow', 'meltRefreeze'))
	terms_univariate <- sort(c('subLethalCold', 'subLethalHeat', 'acuteHeat', 'acuteCold', 'summerRespite', 'chronicCold', 'chronicHeat', 'coldNoSnow', 'meltRefreeze'))

	terms <- c(terms, 'meanDistToClosest4Patches', 'numHomeRanges', 'occy_prev')
	terms_univariate <- c(terms_univariate, 'meanDistToClosest4Patches', 'numHomeRanges', 'occy_prev')

	# Bandolier sites with microclimate over > 1 year
	band <- fread('./Data/Microclimate/Multi-Year 2025-10-23 Erik Beever & Dylan Ryals/band_clim_14.csv')
	band_sites <- tolower(band$patch)
	band_sites[band_sites %in% c('band101a', 'band101b')] <- 'band101'

	# add isolation and home range size of patches
	load('./Data/05 New Mexico Pika - Added PRISM Cell Number & Cell-Based Weight.rda')
	all_sites <- as.data.table(pika)
	all_sites_sites <- tolower(all_sites$polygonName)
	all_sites_sites <- gsub(all_sites_sites, pattern = ' ', replacement = '')

	matches <- match(band_sites, all_sites_sites)
	band[ , meanDistToClosest4Patches := all_sites$meanDistToClosest4Patches[matches]]
	band[ , numHomeRanges := all_sites$numHomeRanges[matches]]

	# amalgamate BAND101A and BAND101B
	years <- sort(unique(band$year))
	for (year in years) {
	
		idx <- which(band$year == year & band$patch %in% c('BAND101A', 'BAND101B'))

		if (length(idx) > 1) {
		
			this <- band[idx]
			this_new <- this[1]

			this_new$patch <- 'BAND101'
			this_new$subLethalCold <- mean(this$subLethalCold, na.rm = TRUE)
			this_new$meltRefreeze <- mean(this$meltRefreeze, na.rm = TRUE)
			this_new$chronicCold <- mean(this$chronicCold, na.rm = TRUE)
			this_new$coldNoSnow <- mean(this$coldNoSnow, na.rm = TRUE)
			this_new$acuteCold <- mean(this$acuteCold, na.rm = TRUE)
			this_new$acuteHeat <- mean(this$acuteHeat, na.rm = TRUE)
			this_new$subLethalHeat <- mean(this$subLethalHeat, na.rm = TRUE)
			this_new$chronicHeat <- mean(this$chronicHeat, na.rm = TRUE)
			this_new$summerRespite <- mean(this$summerRespite, na.rm = TRUE)

			this_new$occy <- max(this$occy, na.rm = TRUE)

			band <- band[-idx]
			band <- rbind(band, this_new)

		}

	
	}

	band$patch[band$patch == 'BAND101A'] <- 'BAND101'

	# center and scale selected columns in 'band'
	cols_to_scale <- c('subLethalCold', 'meltRefreeze', 'chronicCold', 'coldNoSnow', 'acuteCold', 'acuteHeat', 'subLethalHeat', 'chronicHeat', 'summerRespite', 'meanDistToClosest4Patches', 'numHomeRanges')

	band[ , occy := as.numeric(occy)]

	# # remove rows with NA in any of the selected columns
	# band <- band[complete.cases(band[, ..cols_to_scale])]

	# dummy variable for AR1
	band[ , year_recode := year - min(year) + 1]
	band[ , year_recode := as.integer(year_recode)]

	### tally patch by year
	#######################

	# table with 0 = recorded patch status in that year, NA = did not
	patches <- sort(unique(band$patch))
	years <- sort(unique(band$year))

	tallies <- data.table(patch = patches)
	for (year in years) {
	
		tallies[ , this_year := NA_integer_]
		for (patch in patches) {
		
			n <- sum(band$year == year & band$patch == patch)
			if (n > 0) tallies$this_year[tallies$patch == patch] <- 0
		
		}
	
		setnames(tallies, old = 'this_year', new = paste0('year_', year))
		# colnames(tallies)[ncol(tallies)] <- paste0('year_', year)
		# n <- ncol(tallies)

	}

	# re-code non-NA values to 0/1 if detection or not
	for (year in years) {
		idx <- which(!is.na(tallies[[paste0('year_', year)]]))
		for (i in idx) {

			patch_name <- tallies$patch[i]
			occy_val <- band$occy[band$patch == patch_name & band$year == year]
			if (length(occy_val) > 0 && any(occy_val == 1, na.rm = TRUE)) {
				tallies[[paste0('year_', year)]][i] <- 1
			}
		
		}
	}

	year_cols <- names(tallies)[names(tallies) != 'patch']
	tallies[ , n_years_obs := rowSums(!is.na(.SD)), .SDcols = year_cols]

	# # remove any patches that have just one year of observation
	# if (any(tallies$n_years_obs <= 1)) {
	
	# 	bads <- tallies$patch[tallies$n_years_obs <= 1]
	# 	band <- band[patch %notin% bads]
	# 	tallies <- tallies[patch %notin% bads]
	
	# }

	# add predictors
	patches <- tallies$patch
	for (i in seq_along(terms_univariate)) {
	
		term <- terms_univariate[i]
		tallies[ , DUMMY := NA_real_]
		
		for (this_patch in patches) {
			
			y <- band[band$patch == this_patch, ..term]
			y <- y[[1]]
			y <- mean(y)
			tallies$DUMMY[tallies$patch == this_patch] <- y

		}
		
		setnames(tallies, 'DUMMY', term)

	}

	tallies <- tallies[order(meanDistToClosest4Patches)]
	fwrite(tallies, './Figures & Tables/Microclimate Analysis - Multiple Years/Microclimate Analysis - Multiple Years - Tallies of Detections by Site and Year.csv')

	# model formulae
	model_forms <- c(
		'1',
		# '1 + subLethalCold + subLethalHeat',
		'subLethalHeat',
		# 'acuteHeat + acuteCold',
		'acuteHeat',
		# 'acuteCold + acuteHeat + acuteCold:acuteHeat',
		'subLethalHeat + summerRespite',
		'acuteHeat + summerRespite',
		# 'acuteHeat + subLethalCold',
		'summerRespite',
		'chronicCold',
		'acuteCold',
		'subLethalCold',
		'chronicHeat + subLethalHeat',
		'chronicHeat',
		'coldNoSnow',
		'coldNoSnow + subLethalCold',
		'coldNoSnow + acuteCold',
		'coldNoSnow + chronicCold',
		# 'coldNoSnow + summerRespite',
		# 'coldNoSnow + chronicHeat',
		# 'coldNoSnow + subLethalCold + chronicHeat',
		# 'coldNoSnow + subLethalCold + summerRespite',
		# 'subLethalHeat + subLethalCold + subLethalHeat:subLethalCold',
		# 'meltRefreeze + acuteHeat',
		# 'chronicHeat + subLethalCold'
		'meltRefreeze',
		# 'meltRefreeze + summerRespite',
		# 'meltRefreeze + chronicHeat',
		'meltRefreeze + acuteCold',
		'meltRefreeze + subLethalCold'
	)

	model_forms_season <- c(
		'1' = NA,
		# '1 + subLethalCold + subLethalHeat',
		'subLethalHeat' = 'summer',
		# 'acuteHeat + acuteCold',
		'acuteHeat' = 'summer',
		# 'acuteCold + acuteHeat + acuteCold:acuteHeat',
		'subLethalHeat + summerRespite' = 'summer',
		'acuteHeat + summerRespite' = 'summer',
		# 'acuteHeat + subLethalCold',
		'summerRespite' = 'summer',
		'chronicCold' = 'winter',
		'acuteCold' = 'winter',
		'subLethalCold' = 'winter',
		'chronicHeat + subLethalHeat' = 'summer',
		'chronicHeat' = 'summer',
		'coldNoSnow' = 'winter',
		'coldNoSnow + subLethalCold' = 'winter',
		'coldNoSnow + acuteCold' = 'winter',
		'coldNoSnow + chronicCold' = 'winter',
		# 'coldNoSnow + summerRespite',
		# 'coldNoSnow + chronicHeat',
		# 'coldNoSnow + subLethalCold + chronicHeat',
		# 'coldNoSnow + subLethalCold + summerRespite',
		# 'subLethalHeat + subLethalCold + subLethalHeat:subLethalCold',
		# 'meltRefreeze + acuteHeat',
		# 'chronicHeat + subLethalCold'
		'meltRefreeze' = 'winter',
		# 'meltRefreeze + summerRespite',
		# 'meltRefreeze + chronicHeat',
		'meltRefreeze + acuteCold' = 'winter',
		'meltRefreeze + subLethalCold' = 'winter'
	)

	control <- glmmTMBControl(optCtrl = list(iter.max = 2e4, maxfun = 1e5))
	# control <- glmmTMBControl(optCtrl = list(iter.max = 2e4), optimizer = optim, optArgs = list(method = 'BFGS'))
	# control <- glmmTMBControl(optCtrl = list(iter.max = 2e4), optimizer = optim, optArgs = list(method = 'L-BFGS-B'))
	# control <- glmmTMBControl(optCtrl = list(iter.max = 2e4), optimizer = 'bobyqa')

	# temporal autocorrelation formula... site-specific
	ar_form <- ~ ar1(year_recode + 0 | patch)

	# center/scale
	band[ , (cols_to_scale) := lapply(.SD, scale), .SDcols = cols_to_scale]

	# null model for Nagelkerke's R2
	null_model <- glmmTMB(occy ~ 1, family = binomial, data = band, control = control)
	intercept_only_model_ll <- as.numeric(logLik(null_model))
	intercept_only_model_lik <- exp(intercept_only_model_ll)

		# function to remember multi-year model
		remember_multi_year_model <- function(model, form, results, intercept_only_model_lik, site_intercept, tac, season) {

			# if (model$fit$convergence != 0) stop('Convergence failure on model ', form, '. Value ', model$fit$convergence, '.')
			if (model$fit$convergence != 0) return(results)

			# remember model
			n <- nrow(band)
			ll <- as.numeric(logLik(model))
			lik_mean <- ll / n
			lik <- exp(ll)

			nr2 <- nagelR2(likeNull = intercept_only_model_lik, likeFull = lik, n = n)

			# aicc <- AICc(model)

			# get deviance from fitted glmmTMB model
			deviance <- -2 * as.numeric(logLik(model))

			this_results <- data.table(
				model = form,
				season = season,
				site_intercept = site_intercept,
				# tac = tac,
				n_obs = n,
				log_like = ll,
				lik_mean = lik_mean,
				pseudo_r2 = nr2,
				deviance = deviance
			)

			# remember coefficients			
			for (count_term in seq_along(terms)) {
				this_results[ , DUMMY := NA_real_]
				colnames(this_results)[ncol(this_results)] <- terms[count_term]
			}

			coeffs <- fixef(model)$cond
			names(coeffs) <- gsub(names(coeffs), pattern = 'cond\\(\\(', replacement = '')
			names(coeffs) <- gsub(names(coeffs), pattern = 'cond\\(', replacement = '')
			names(coeffs) <- gsub(names(coeffs), pattern = '\\)', replacement = '')

			coeffs <- coeffs[names(coeffs) != 'Int']

			for (count_coeff in seq_along(coeffs)) {
				
				col <- which(names(this_results) == names(coeffs)[count_coeff])
				this_results[1, col] <- coeffs[count_coeff]

			}

			results <- rbind(results, this_results)
			results

		}

	results <- data.table()
	n_forms <- length(model_forms)
	for (previous_year_effect in c(FALSE, TRUE)) {
	
		for (i in 1:n_forms) {
		
			this_form <- model_forms[i]
			if (previous_year_effect) this_form <- paste0(this_form, ' + occy_prev')

			this_season <- model_forms_season[i]
			say(this_form)

			this_form <- paste('occy ~', this_form)
			this_form_iso <- paste(this_form, ' + meanDistToClosest4Patches')
			this_form_hr <- paste(this_form, ' + numHomeRanges')
			this_form_iso_hr <- paste(this_form, ' + meanDistToClosest4Patches + numHomeRanges')

			### models: no TAC and no site effects
			######################################

			site_intercept <- FALSE
			tac <- FALSE

			# model with no TAC, no site effect, no HR, no isolation
			form <- this_form
			model <- glmmTMB(as.formula(form),
				family = binomial,
				data = band,
				control = control
			)

			results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# model with no TAC, no site effect, no HR, YES isolation
			form <- this_form_iso
			model <- glmmTMB(as.formula(form),
				family = binomial,
				data = band,
				control = control
			)

			results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# model with no TAC, no site effect, YES HR, NO isolation
			form <- this_form_hr
			model <- glmmTMB(as.formula(form),
				family = binomial,
				data = band,
				control = control
			)

			results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# model with no TAC, no site effect, YES HR, YES isolation
			form <- this_form_iso_hr
			model <- glmmTMB(as.formula(form),
				family = binomial,
				data = band,
				control = control
			)

			results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			### models: no TAC and YES site effects
			#######################################

			site_intercept <- TRUE
			tac <- FALSE

			# model with no TAC, YES site effect, no HR, no isolation
			form <- this_form
			form <- paste(form, '+ (1 | patch)')
			model <- glmmTMB(as.formula(form),
				family = binomial,
				data = band,
				control = control
			)

			results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# model with no TAC, YES site effect, no HR, YES isolation
			form <- this_form_iso
			form <- paste(form, '+ (1 | patch)')
			model <- glmmTMB(as.formula(form),
				family = binomial,
				data = band,
				control = control
			)

			results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# model with no TAC, YES site effect, YES HR, NO isolation
			form <- this_form_hr
			form <- paste(form, '+ (1 | patch)')
			model <- glmmTMB(as.formula(form),
				family = binomial,
				data = band,
				control = control
			)

			results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# model with no TAC, YES site effect, YES HR, YES isolation
			form <- this_form_iso_hr
			form <- paste(form, '+ (1 | patch)')
			model <- glmmTMB(as.formula(form),
				family = binomial,
				data = band,
				control = control
			)

			results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# ### models: YES TAC and NO site effects
			# #######################################

			# site_intercept <- FALSE
			# tac <- TRUE

			# # model with YES TAC, no site effect, no HR, no isolation
			# form <- this_form
			# model <- glmmTMB(as.formula(form),
			# 	family = binomial,
			# 	data = band,
			# 	control = control,
			# 	dispformula = ar_form
			# )

			# results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# # model with YES TAC, no site effect, no HR, YES isolation
			# form <- this_form_iso
			# model <- glmmTMB(as.formula(form),
			# 	family = binomial,
			# 	data = band,
			# 	control = control,
			# 	dispformula = ar_form
			# )

			# results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# # model with YES TAC, no site effect, YES HR, NO isolation
			# form <- this_form_hr
			# model <- glmmTMB(as.formula(form),
			# 	family = binomial,
			# 	data = band,
			# 	control = control,
			# 	dispformula = ar_form
			# )

			# results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# # model with YES TAC, no site effect, YES HR, YES isolation
			# form <- this_form_iso_hr
			# model <- glmmTMB(as.formula(form),
			# 	family = binomial,
			# 	data = band,
			# 	control = control,
			# 	dispformula = ar_form
			# )

			# results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# ### models: YES TAC and YES site effects
			# ########################################

			# site_intercept <- TRUE
			# tac <- TRUE

			# # model with YES TAC, YES site effect, no HR, no isolation
			# form <- this_form
			# form <- paste(form, '+ (1 | patch)')
			# model <- glmmTMB(as.formula(form),
			# 	family = binomial,
			# 	data = band,
			# 	control = control,
			# 	dispformula = ar_form
			# )

			# results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# # model with YES TAC, YES site effect, no HR, YES isolation
			# form <- this_form_iso
			# form <- paste(form, '+ (1 | patch)')
			# model <- glmmTMB(as.formula(form),
			# 	family = binomial,
			# 	data = band,
			# 	control = control,
			# 	dispformula = ar_form
			# )

			# results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# # model with YES TAC, YES site effect, YES HR, NO isolation
			# form <- this_form_hr
			# form <- paste(form, '+ (1 | patch)')
			# model <- glmmTMB(as.formula(form),
			# 	family = binomial,
			# 	data = band,
			# 	control = control,
			# 	dispformula = ar_form
			# )

			# results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# # model with YES TAC, YES site effect, YES HR, YES isolation
			# form <- this_form_iso_hr
			# form <- paste(form, '+ (1 | patch)')
			# model <- glmmTMB(as.formula(form),
			# 	family = binomial,
			# 	data = band,
			# 	control = control,
			# 	dispformula = ar_form
			# )

			# results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

		} # next model form

	} # include previous year site status as predictor?

	results <- results[order(lik_mean, decreasing = TRUE)]
	# results[ , delta_aicc := aicc - min(aicc)]
	# results[ , weight := exp(-0.5 * delta_aicc) / sum(exp(-0.5 * delta_aicc))]

	fwrite(results, './Figures & Tables/Microclimate Analysis - Multiple Years/Microclimate Analysis - Multiple Years - Models.csv')

	### variable importance
	#######################

	var_imp <- data.table(
		coefficient = terms,
		n_models = NA_integer_,
		n_models_positive = NA_integer_,
		n_models_negative = NA_integer_,
		mean_rank = NA_real_
	)

	results$rank <- 1:nrow(results)
	for (i in seq_along(terms)) {
	
		term <- terms[i]
		col <- which(names(results) == term)
		idx <- which(complete.cases(results[ , ..col]))

		if (length(idx) > 0) {

			this_results <- results[idx]
			coeff <- this_results[[term]]

			var_imp$n_models[i] <- nrow(this_results)
			var_imp$n_models_positive[i] <- sum(coeff > 0)
			var_imp$n_models_negative[i] <- sum(coeff < 0)
			var_imp$mean_rank[i] <- mean(this_results$rank)

		}
		
	}

	var_imp <- var_imp[order(mean_rank, decreasing = FALSE)]

	fwrite(var_imp, './Figures & Tables/Microclimate Analysis - Multiple Years/Microclimate Analysis - Multiple Years - Variable Importance.csv')
 
# 	# fwrite(band, './Figures & Tables/Microclimate Analysis - Multiple Years/Microclimate Analysis - Multiple Years - Data Centered & Scaled.csv')

