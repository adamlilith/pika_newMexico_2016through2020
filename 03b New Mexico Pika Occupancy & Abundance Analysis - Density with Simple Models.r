### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('C:/Ecology/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/03b New Mexico Pika Occupancy & Abundance Analysis - Density with Simple Models.r')
### source('E:/Adam/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/03b New Mexico Pika Occupancy & Abundance Analysis - Density with Simple Models.r')
###
### CONTENTS ###
### setup ###
### climate/isolation DENSITY models ###
### compile table of predictor importance for climate/isolation DENSITY models ###
### climate/isolation/ecology DENSITY models ###
### double-checking number of models each variable should be in ###

#############
### setup ###
#############

	rm(list=ls())

	drive <- 'C:/Ecology/'
	# drive <- 'E:/Adam/'

	source('C:/Ecology/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r')
	# source('E:/Adam/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r')

say('########################################')
say('### climate/isolation DENSITY models ###')
say('########################################')

	load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')
	
	pika$meanDistToClosest4Patches <- log10(pika$meanDistToClosest4Patches)
	
	pika <- pika[!is.na(pika$latestDensity), ]
	pika$region <- as.factor(pika$region)
	vars <- getVars('density')
	vars <- c(vars, 'meanDistToClosest4Patches')
	pika[ , vars] <- scale(pika[ , vars])

	### intercept-only models
	#########################

		model1 <- glm(latestDensity ~ 1, data=pika, family=Gamma(link='log'))
		model2 <- glm(latestDensity ~ 1 + region, data=pika, family=Gamma(link='log'))
		
		model3 <- glm(latestDensity ~ meanDistToClosest4Patches, data=pika, family=Gamma(link='log'))
		model4 <- glm(latestDensity ~ 1 + region + meanDistToClosest4Patches, data=pika, family=Gamma(link='log'))
		
		region <- c(FALSE, TRUE, FALSE, TRUE)
		
		aicc1 <- AICc(model1)
		aicc2 <- AICc(model2)
		aicc3 <- AICc(model3)
		aicc4 <- AICc(model4)
		aicc <- c(aicc1, aicc2, aicc3, aicc4)
		
		ll1 <- logLik(model1)
		ll2 <- logLik(model2)
		ll3 <- logLik(model3)
		ll4 <- logLik(model4)
		llNull <- ll1
		
		pseudoR2_1 <- nagelR2(llNull, ll1, n=nrow(pika))
		pseudoR2_2 <- nagelR2(llNull, ll2, n=nrow(pika))
		pseudoR2_3 <- nagelR2(llNull, ll3, n=nrow(pika))
		pseudoR2_4 <- nagelR2(llNull, ll4, n=nrow(pika))
		pseudoR2 <- c(pseudoR2_1, pseudoR2_2, pseudoR2_3, pseudoR2_4)
		
		coef3 <- coefficients(model3)
		coef4 <- coefficients(model4)
		
		isolationCoef3 <- coef3[['meanDistToClosest4Patches']]
		isolationCoef4 <- coef4[['meanDistToClosest4Patches']]
		
		isolationCoef <- c(NA, NA, isolationCoef3, isolationCoef4)
		
		results <- data.frame(
			model = '(Intercept)',
			term1 = NA,
			term2 = NA,
			term3 = NA,
			isolationCoef = isolationCoef,
			region = region,
			aicc = aicc,
			pseudoR2 = pseudoR2
		)

	### climate models
	##################
	
		formulae <- getFormulaeDens()

		for (formula in formulae) {
	
			say(formula)

			form1 <- as.formula(paste0('latestDensity ~ 1 + ', formula))
			form2 <- as.formula(paste0('latestDensity ~ 1 + ', formula, ' + region'))
			form3 <- as.formula(paste0('latestDensity ~ 1 + ', formula, ' + meanDistToClosest4Patches'))
			form4 <- as.formula(paste0('latestDensity ~ 1 + ', formula, ' + region + meanDistToClosest4Patches'))
	
			model1 <- glm(form1, data=pika, family=Gamma(link='log'))
			model2 <- glm(form2, data=pika, family=Gamma(link='log'))
			model3 <- glm(form3, data=pika, family=Gamma(link='log'))
			model4 <- glm(form4, data=pika, family=Gamma(link='log'))
			
			terms <- extractTerms(model1, model2, model3, model4)
			term1 <- terms$term1
			term2 <- terms$term2
			term3 <- terms$term3
			
			region <- c(FALSE, TRUE, FALSE, TRUE)
			
			aicc1 <- AICc(model1)
			aicc2 <- AICc(model2)
			aicc3 <- AICc(model3)
			aicc4 <- AICc(model4)
			aicc <- c(aicc1, aicc2, aicc3, aicc4)
			
			ll1 <- logLik(model1)
			ll2 <- logLik(model2)
			ll3 <- logLik(model3)
			ll4 <- logLik(model4)
			
			pseudoR2_1 <- nagelR2(llNull, ll1, n=nrow(pika))
			pseudoR2_2 <- nagelR2(llNull, ll2, n=nrow(pika))
			pseudoR2_3 <- nagelR2(llNull, ll3, n=nrow(pika))
			pseudoR2_4 <- nagelR2(llNull, ll4, n=nrow(pika))
			pseudoR2 <- c(pseudoR2_1, pseudoR2_2, pseudoR2_3, pseudoR2_4)
			
			coef3 <- coefficients(model3)
			coef4 <- coefficients(model4)
			
			isolationCoef3 <- coef3[['meanDistToClosest4Patches']]
			isolationCoef4 <- coef4[['meanDistToClosest4Patches']]
			
			isolationCoef <- c(NA, NA, isolationCoef3, isolationCoef4)
			
			results <- rbind(
				results,
				data.frame(
					model = formula,
					term1 = term1,
					term2 = term2,
					term3 = term3,
					isolationCoef = isolationCoef,
					region = region,
					aicc = aicc,
					pseudoR2 = pseudoR2
				)
			)
			
		} # next formula

		### reports
		###########
		
			results$deltaAicc <- results$aicc - min(results$aicc)
			w <- exp(-0.5 * results$deltaAicc)
			results$weight <- w / sum(w)

			results <- results[order(results$weight, decreasing=TRUE), ]
			rownames(results) <- NULL

			file <- paste0('./Figures & Tables/Density - Simple Models/Density - Simple GLMs.csv')
			write.csv(results, file, row.names=FALSE)
	
say('##################################################################################')
say('### compile table of predictor importance for climate/isolation DENSITY models ###')
say('##################################################################################')
	
	# rank variables by mean AICc weight

	# get models
	file <- paste0('./Figures & Tables/Density - Simple Models/Density - Simple GLMs.csv')
	models <- read.csv(file)
	
	# get variables
	vars <- getVars('density')
	vars <- c(vars, 'meanDistToClosest4Patches')

	imp <- data.frame()
	for (var in vars) {
	
		index <- if (var != 'meanDistToClosest4Patches') {
			which(grepl(models$model, pattern=var))
		} else {
			which(!is.na(models$isolationCoef))
		}
		n <- length(index)
		sumWeight <- sum(models$weight[index])
		meanWeight <- sumWeight / n
		
		imp <- rbind(
			imp,
			data.frame(
				variable = var,
				niceVar = makeNiceVars(var, 'density'),
				numModels = n,
				sumWeight = sumWeight,
				meanWeight = meanWeight
			)
		)
		
	}
	
	imp <- imp[order(imp$meanWeight, decreasing=TRUE), ]

	write.csv(imp, paste0('./Figures & Tables/Density - Simple Models/Density - Simple GLMs - Var Import.csv'), row.names=FALSE)
	
say('################################################')
say('### climate/isolation/ecology DENSITY models ###')
say('################################################')

	### collate data
	################

		### data
		load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')
		pika <- pika[!is.na(pika$latestDensity), ]

		### collate "latest" version of each predictor
		ecoVars <- c('perimBurnedMostRecent_perc', 'latestGrazing', 'latestGrassForb') # NB grass, forbs, and grass/form are highly correlated

		# pika$latestGrazing <- pika$latestGrass <- pika$latestForb <- pika$latestGrassForb <- NA
		pika$latestGrazing <- pika$latestGrassForb <- NA
		
		for (i in 1:nrow(pika)) {
		
			pika$latestGrazing[i] <- pika[i, paste0('grazing', pika$latestDensSurveyYear[i])]
			# pika$latestGrass[i] <- pika[i, paste0('grass', pika$latestDensSurveyYear[i], '_perc')]
			# pika$latestForb[i] <- pika[i, paste0('forb', pika$latestDensSurveyYear[i], '_perc')]
			pika$latestGrassForb[i] <- pika[i, paste0('grassForb', pika$latestDensSurveyYear[i], '_perc')]
		
		
		}
		
		# remove any records with NAs in predictors
		nas <- which(
			# is.na(pika$latestGrazing) | is.na(pika$perimBurnedMostRecent_perc) | is.na(pika$latestGrass) | is.na(pika$latestForb) | is.na(pika$latestGrassForb)
			is.na(pika$latestGrazing) | is.na(pika$perimBurnedMostRecent_perc) | is.na(pika$latestGrassForb)
		)

		if (length(nas) > 0) pika <- pika[-nas, ]

		pika[ , ecoVars[ecoVars != 'latestGrazing']] <- scale(pika[ , ecoVars[ecoVars != 'latestGrazing']])
		pika$meanDistToClosest4Patches <- log10(pika$meanDistToClosest4Patches)
		pika$meanDistToClosest4Patches <- scale(pika$meanDistToClosest4Patches)
		
		pika$region <- as.factor(pika$region)
		
	### models
	##########
	
		allClimModels <- read.csv('./Figures & Tables/Density - Simple Models/Density - Simple GLMs.csv')
		topClimModels <- allClimModels[allClimModels$deltaAicc < maxDeltaAic_density, ]

		report <- data.frame()
		
		nullModel <- glm(latestDensity ~ 1, data=pika, family=Gamma(link='log'))
		nullModelRegion <- glm(latestDensity ~ region, data=pika, family=Gamma(link='log'))
		nullModelIsolation <- glm(latestDensity ~ meanDistToClosest4Patches, data=pika, family=Gamma(link='log'))
		nullModelRegionIsolation <- glm(latestDensity ~ region + meanDistToClosest4Patches, data=pika, family=Gamma(link='log'))
		
		llNull <- logLik(nullModel)
		llNullRegion <- logLik(nullModelRegion)
		llNullIsolation <- logLik(nullModelIsolation)
		llNullRegionIsolation <- logLik(nullModelRegionIsolation)

		ecoTermsGrid <- expand.grid(a = c(TRUE, FALSE), b = c(TRUE, FALSE), c = c(TRUE, FALSE))
		names(ecoTermsGrid) <- ecoVars
		
		for (countClimModel in 1:nrow(topClimModels)) {

			# climate formula
			climForm <- topClimModels$model[countClimModel]
			climForm <- paste0('latestDensity ~ ', climForm)
			hasRegion <- topClimModels$region[countClimModel]
			if (hasRegion) climForm <- paste(climForm, ' + region')
			hasIsolation <- !is.na(topClimModels$isolationCoef[countClimModel])
			if (hasIsolation) climForm <- paste0(climForm, ' + meanDistToClosest4Patches')
				
			for (countEcoModel in 1:nrow(ecoTermsGrid)) {

				climEcoForm <- climForm

				# eco-variable formula
				ecoTermsInModel <- unlist(ecoTermsGrid[countEcoModel, ])
				ecoTermsInModel <- ecoVars[ecoTermsInModel]
				ecoTermsInModel <- paste(ecoTermsInModel, collapse = ' + ')
				if (ecoTermsInModel != '') climEcoForm <- paste(climEcoForm, '+', ecoTermsInModel)

				# model
				model <- glm(climEcoForm, data=pika, family=Gamma(link='log'))
				terms <- extractTerms(model)

				term1 <- terms$term1
				term2 <- if (is.null(terms$term2)) { c(NA, NA) } else { terms$term2 }
				term3 <- if (is.null(terms$term3)) { c(NA, NA) } else { terms$term3 }
				term4 <- if (is.null(terms$term4)) { c(NA, NA) } else { terms$term4 }
				term5 <- if (is.null(terms$term5)) { c(NA, NA) } else { terms$term5 }
				term6 <- if (is.null(terms$term6)) { c(NA, NA) } else { terms$term6 }
				term7 <- if (is.null(terms$term7)) { c(NA, NA) } else { terms$term7 }
				term8 <- if (is.null(terms$term8)) { c(NA, NA) } else { terms$term8 }
				term9 <- if (is.null(terms$term9)) { c(NA, NA) } else { terms$term9 }
				
				aicc <- AICc(model)
				
				ll <- logLik(model)
				pseudoR2 <- nagelR2(llNull, ll, n=nrow(pika))
				
				# report
				report <- rbind(
					report,
					data.frame(
						model = climEcoForm,
						term1 = term1,
						term2 = term2,
						term3 = term3,
						term4 = term4,
						term5 = term5,
						term6 = term6,
						term7 = term7,
						term8 = term8,
						term9 = term9,
						region = hasRegion,
						hasIsolation = hasIsolation,
						aiccClim = topClimModels$aicc[countClimModel],
						aiccClimEco = aicc,
						pseudoR2Clim = pseudoR2,
						pseudoR2ClimEco = topClimModels$pseudoR2[countClimModel]
					)
				)
				
			} # next eco-term model
		
		} # next climate model
		
		# intercept-only model
		aiccClim <- allClimModels$aicc[allClimModels$model == '(Intercept)' & is.na(allClimModels$isolation) & !allClimModels$region]
		
		report <- rbind(
			report,
			data.frame(
				model = 'latestDensity ~ 1',
				term1 = NA,
				term2 = NA,
				term3 = NA,
				term4 = NA,
				term5 = NA,
				term6 = NA,
				term7 = NA,
				term8 = NA,
				term9 = NA,
				region = FALSE,
				hasIsolation = FALSE,
				aiccClim = aiccClim,
				aiccClimEco = AICc(nullModel),
				pseudoR2Clim = 0,
				pseudoR2ClimEco = NA
			)
		)

		# region-only model
		aiccClim <- allClimModels$aicc[allClimModels$model == '(Intercept)' & is.na(allClimModels$isolation) & allClimModels$region]

		report <- rbind(
			report,
			data.frame(
				model = 'latestDensity ~ 1',
				term1 = NA,
				term2 = NA,
				term3 = NA,
				term4 = NA,
				term5 = NA,
				term6 = NA,
				term7 = NA,
				term8 = NA,
				term9 = NA,
				region = TRUE,
				hasIsolation = FALSE,
				aiccClim = aiccClim,
				aiccClimEco = AICc(nullModelRegion),
				pseudoR2Clim = nagelR2(llNull, llNullRegion, n=nrow(pika)),
				pseudoR2ClimEco = NA
			)
		)

		# isolation-only model
		aiccClim <- allClimModels$aicc[allClimModels$model == '(Intercept)' & !is.na(allClimModels$isolation) & !allClimModels$region]

		report <- rbind(
			report,
			data.frame(
				model = 'latestDensity ~ 1',
				term1 = NA,
				term2 = NA,
				term3 = NA,
				term4 = NA,
				term5 = NA,
				term6 = NA,
				term7 = NA,
				term8 = NA,
				term9 = NA,
				region = FALSE,
				hasIsolation = TRUE,
				aiccClim = aiccClim,
				aiccClimEco = AICc(nullModelIsolation),
				pseudoR2Clim = nagelR2(llNull, llNullIsolation, n=nrow(pika)),
				pseudoR2ClimEco = NA
			)
		)

		# isolation/region-only model
		aiccClim <- allClimModels$aicc[allClimModels$model == '(Intercept)' & !is.na(allClimModels$isolation) & !allClimModels$region]

		report <- rbind(
			report,
			data.frame(
				model = 'latestDensity ~ 1',
				term1 = NA,
				term2 = NA,
				term3 = NA,
				term4 = NA,
				term5 = NA,
				term6 = NA,
				term7 = NA,
				term8 = NA,
				term9 = NA,
				region = TRUE,
				hasIsolation = TRUE,
				aiccClim = aiccClim,
				aiccClimEco = AICc(nullModelRegionIsolation),
				pseudoR2Clim = nagelR2(llNull, llNullRegionIsolation, n=nrow(pika)),
				pseudoR2ClimEco = NA
			)
		)

	### report
	##########
	
		report$deltaAiccClim <- report$aiccClim - min(report$aiccClim)
		w <- exp(-0.5 * report$deltaAiccClim)
		report$weightClim <- w / sum(w)

		report$deltaAiccClimEco <- report$aiccClimEco - min(report$aiccClimEco)
		w <- exp(-0.5 * report$deltaAiccClimEco)
		report$weightClimEco <- w / sum(w)

		report <- report[order(report$aiccClimEco), ]

		rownames(report) <- NULL

		file <- paste0('./Figures & Tables/Density - Simple Models/Top Climate & Isolation Density Models with Ecological Variables.csv')
		write.csv(report, file, row.names=FALSE)
	
say('###################################################################')
say('### double-checking number of models each variable should be in ###')
say('###################################################################')
	
	vars <- getVars('density')
	formulae <- getFormulaeDens()
	
	counts <- data.frame()
	for (var in vars) {
	
		ins <- sum(grepl(var, formulae))
		counts <- rbind(
			counts,
			data.frame(
				var = var,
				n = ins
			)
		)
	
	}
	
	write.csv(counts, './Figures & Tables/Density - Simple Models/Number of Base Models with Each Variable.csv', row.names=FALSE)
	
say('DONE!!!', level=1, deco='%')
