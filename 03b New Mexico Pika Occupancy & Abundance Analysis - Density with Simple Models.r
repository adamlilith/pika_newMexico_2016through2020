### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### Analysis of sites where non-zero density of pikas was recorded. Constructs and compares multiple simple models of density and estimates of variable importance.
###
### source('C:/Kaji/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/03b New Mexico Pika Occupancy & Abundance Analysis - Density with Simple Models.r')
###
### CONTENTS ###
### setup ###
### climate/biogeography DENSITY models ###
### compile table of predictor importance for climate/biogeography DENSITY models ###
### climate/biogeography/management DENSITY models ###
### double-checking number of models each variable should be in ###
### make response plots for top climate/biogeography/management DENSITY models ###

#############
### setup ###
#############

	rm(list=ls())

	drive <- 'C:/Kaji/'

	source(paste0(drive, '/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r'))

	# use "extended" list of climate models provided by Erik 2025-09
	# affects all analyses herein
	extended <- TRUE # original plus extended set
	# extended <- FALSE # original set

	analysisDate <- Sys.time()
	analysisDate <- substr(analysisDate, 1, 10)

	if (extended) {
		say('Using EXTENDED set of climate models!!!', level = 2)
		dirCreate('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models ', analysisDate)
	} else {
		say('Using ORIGINAL set of climate models!!!', level = 2)
		dirCreate('./Figures & Tables/Density - Simple Models with Original Set of Climate Models ', analysisDate)
	}

say('###########################################')
say('### climate/biogeography DENSITY models ###')
say('###########################################')

	# Models of pika density at sites where pika were detected. Models use only climate and biogeographic (isolation, patch size) variables.

	for (focalRegion in c('all', 'southwest', 'non-southwest')) {

		say(focalRegion, level = 2)

		load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')

		if (focalRegion == 'southwest') pika <- pika[pika$region == 'southwest', ]
		if (focalRegion == 'non-southwest') pika <- pika[pika$region != 'southwest', ]

		pika$meanDistToClosest4Patches <- log10(pika$meanDistToClosest4Patches)
		
		pika <- pika[!is.na(pika$latestDensity), ]
		pika$region <- as.factor(pika$region)
		vars <- getVars('density')
		vars <- c(vars, 'meanDistToClosest4Patches')
		pika[ , vars] <- scale(pika[ , vars])

		### intercept-only models
		#########################

			model1 <- glm(latestDensity ~ 1, data=pika, family=Gamma(link='log'))
			if (focalRegion != 'southwest') model2 <- glm(latestDensity ~ 1 + region, data=pika, family=Gamma(link='log'))
			model3 <- glm(latestDensity ~ meanDistToClosest4Patches, data=pika, family=Gamma(link='log'))
			if (focalRegion != 'southwest') model4 <- glm(latestDensity ~ 1 + region + meanDistToClosest4Patches, data=pika, family=Gamma(link='log'))
			
			region <- if (focalRegion != 'southwest') {
				c(FALSE, TRUE, FALSE, TRUE)
			} else {
				c(FALSE, TRUE, FALSE, NA)
			}
			
			aicc1 <- AICc(model1)
			aicc2 <- if (focalRegion != 'southwest') { aicc2 <- AICc(model2) } else { NA }
			aicc3 <- AICc(model3)
			aicc4 <- if (focalRegion != 'southwest') { AICc(model4) } else { NA }
			aicc <- c(aicc1, aicc2, aicc3, aicc4)
			
			ll1 <- logLik(model1)
			ll2 <- if (focalRegion != 'southwest') { logLik(model2) } else { NA }
			ll3 <- logLik(model3)
			ll4 <- if (focalRegion != 'southwest') { logLik(model4) } else { NA }	
			llNull <- ll1

			# pseudoR2_1 <- nagelR2(llNull, ll1, n=nrow(pika))
			# pseudoR2_2 <- if (focalRegion != 'southwest') { nagelR2(llNull, ll2, n=nrow(pika)) } else { NA }
			# pseudoR2_3 <- nagelR2(llNull, ll3, n=nrow(pika))
			# pseudoR2_4 <- if (focalRegion != 'southwest') { nagelR2(llNull, ll4, n=nrow(pika)) } else { NA }
			# pseudoR2 <- c(pseudoR2_1, pseudoR2_2, pseudoR2_3, pseudoR2_4)
			
			modelNull <- model1
			pseudoR2_1 <- 1 - deviance(model1) / deviance(modelNull)
			pseudoR2_2 <- if (focalRegion != 'southwest') { 1 - deviance(model2) / deviance(modelNull) } else { NA }
			pseudoR2_3 <- 1 - deviance(model3) / deviance(modelNull)
			pseudoR2_4 <- if (focalRegion != 'southwest') { 1 - deviance(model4) / deviance(modelNull) } else { NA }
			pseudoR2 <- c(pseudoR2_1, pseudoR2_2, pseudoR2_3, pseudoR2_4)
			
			coef3 <- coefficients(model3)
			coef4 <- if (focalRegion != 'southwest') { coefficients(model4) } else { NA }
			
			isolationCoef3 <- coef3[['meanDistToClosest4Patches']]
			isolationCoef4 <- if (focalRegion != 'southwest') { coef4[['meanDistToClosest4Patches']] } else { NA }
			
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
			if (extended) {
			
				formulaeExtended <- getFormulaeDensExtended()
				formulae <- sort(unique(c(formulae, formulaeExtended)))

			}

			for (formula in formulae) {
		
				say(formula)

				form1 <- as.formula(paste0('latestDensity ~ 1 + ', formula))
				if (focalRegion != 'southwest') form2 <- as.formula(paste0('latestDensity ~ 1 + ', formula, ' + region'))
				form3 <- as.formula(paste0('latestDensity ~ 1 + ', formula, ' + meanDistToClosest4Patches'))
				if (focalRegion != 'southwest') form4 <- as.formula(paste0('latestDensity ~ 1 + ', formula, ' + region + meanDistToClosest4Patches'))
		
				model1 <- glm(form1, data=pika, family=Gamma(link='log'))
				if (focalRegion != 'southwest') model2 <- glm(form2, data=pika, family=Gamma(link='log'))
				model3 <- glm(form3, data=pika, family=Gamma(link='log'))
				if (focalRegion != 'southwest') model4 <- glm(form4, data=pika, family=Gamma(link='log'))
				
				if (focalRegion != 'southwest') {
					terms <- extractTerms(model1, model2, model3, model4)
				} else {
					terms <- extractTerms(model1, model3)
				}
				term1 <- terms$term1
				term2 <- terms$term2
				term3 <- terms$term3
				
				hasRegion <- if (focalRegion != 'southwest') { TRUE } else { NA }
				region <- c(FALSE, hasRegion, FALSE, hasRegion)
				
				aicc1 <- AICc(model1)
				aicc2 <- if (focalRegion != 'southwest') { AICc(model2) } else { NA }
				aicc3 <- AICc(model3)
				aicc4 <- if (focalRegion != 'southwest') { AICc(model4) } else { NA }
				aicc <- c(aicc1, aicc2, aicc3, aicc4)
				
				ll1 <- logLik(model1)
				ll2 <- if (focalRegion != 'southwest') { logLik(model2) } else { NA }
				ll3 <- logLik(model3)
				ll4 <- if (focalRegion != 'southwest') { logLik(model4) } else { NA }
				
				# pseudoR2_1 <- nagelR2(llNull, ll1, n=nrow(pika))
				# pseudoR2_2 <- if (focalRegion != 'southwest') { nagelR2(llNull, ll2, n=nrow(pika)) } else { NA }
				# pseudoR2_3 <- nagelR2(llNull, ll3, n=nrow(pika))
				# pseudoR2_4 <- if (focalRegion != 'southwest') { nagelR2(llNull, ll4, n=nrow(pika)) } else { NA }
				# pseudoR2 <- c(pseudoR2_1, pseudoR2_2, pseudoR2_3, pseudoR2_4)
				
				pseudoR2_1 <- 1 - deviance(model1) / deviance(modelNull)
				pseudoR2_2 <- if (focalRegion != 'southwest') { 1 - deviance(model2) / deviance(modelNull) } else { NA }
				pseudoR2_3 <- 1 - deviance(model3) / deviance(modelNull)
				pseudoR2_4 <- if (focalRegion != 'southwest') { 1 - deviance(model4) / deviance(modelNull) } else { NA }
				pseudoR2 <- c(pseudoR2_1, pseudoR2_2, pseudoR2_3, pseudoR2_4)

				coef3 <- coefficients(model3)
				coef4 <- if (focalRegion != 'southwest') { coefficients(model4) } else { NA }
				
				isolationCoef3 <- coef3[['meanDistToClosest4Patches']]
				isolationCoef4 <- if (focalRegion != 'southwest') { coef4[['meanDistToClosest4Patches']] } else { NA }
				
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
	
				if (focalRegion == 'southwest') results <- results[!is.na(results$aicc), ]
				results$deltaAicc <- results$aicc - min(results$aicc)
				w <- exp(-0.5 * results$deltaAicc)
				results$weight <- w / sum(w)

				results <- results[order(results$weight, decreasing=TRUE), ]
				rownames(results) <- NULL

				file <- if (extended) {
					paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models ', analysisDate, '/Density - Models without Management Variables ', capIt(focalRegion), ' Region ', analysisDate, '.csv')
				} else {
					paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models ', analysisDate, '/Density - Models without Management Variables ', capIt(focalRegion), ' Region ', analysisDate, '.csv')
				}
				write.csv(results, file, row.names=FALSE)

	}

say('#####################################################################################')
say('### compile table of predictor importance for climate/biogeography DENSITY models ###')
say('#####################################################################################')
	
	# Assess importance of variables in models of pika density at sites where pika were detected. Models use only climate and biogeographic (isolation, patch size) variables.

	for (focalRegion in c('all', 'southwest', 'non-southwest')) {

		# get models
		file <- if (extended) {
			paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models ', analysisDate, '/Density - Models without Management Variables ', capIt(focalRegion), ' Region ', analysisDate, '.csv')
		} else {
			paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models ', analysisDate, '/Density - Models without Management Variables ', capIt(focalRegion), ' Region ', analysisDate, '.csv')
		}
		models <- read.csv(file)
		
		# get variables
		vars <- getVars('density')
		vars <- c(vars, 'meanDistToClosest4Patches')

		imp <- data.frame()
		for (var in vars) {
		
			index <- if (var != 'meanDistToClosest4Patches') {
				which(grepl(pattern=var, x=models$model))
			} else {
				which(!is.na(models$isolationCoef))
			}
			n <- length(index)
			sumWeight <- sum(models$weight[index])
			meanWeight <- sumWeight / n

			# number of pos/necg cases
			forms <- models[index, c('model', 'term1', 'term2', 'term3'), drop = FALSE]
			negatives <- positives <- 0
			
			if (var != 'meanDistToClosest4Patches') {
				
				for (i in 1:nrow(forms)) {
							
					modelExploded <- strsplit(forms$model[i], split = ' \\+ ')[[1]]
					match <- which(modelExploded == var)
					term <- forms[i, paste0('term', match), drop = TRUE]
					if (substr(term, 1, 1) == '-') {
						negatives <- negatives + 1
					} else {
						positives <- positives + 1
					}
					
				}
				
			} else if (var == 'meanDistToClosest4Patches') {
				
				negatives <- sum(substr(models$isolationCoef, 1, 1) == '-', na.rm = TRUE)
				positives <- sum(substr(models$isolationCoef, 1, 1) != '-' & substr(models$isolationCoef, 1, 2) != 'NA', na.rm = TRUE)
				
			}

			imp <- rbind(
				imp,
				data.frame(
					variable = var,
					niceVar = makeNiceVars(var, 'density'),
					numModels = n,
					sumWeight = sumWeight,
					meanWeight = meanWeight,
					negatives = negatives,
					positives = positives
				)
			)
			
		}
		
		imp <- imp[order(imp$meanWeight, decreasing=TRUE), ]

		file <- if (extended) {
			paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models ', analysisDate, '/Density - Variable Importance in Models without Management Variables ', capIt(focalRegion), ' Region ', analysisDate, '.csv')
		} else {
			paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models ', analysisDate, '/Density - Variable Importance in Models without Management Variables ', capIt(focalRegion), ' Region ', analysisDate, '.csv')
		}
		write.csv(imp, file, row.names=FALSE)

	}

say('######################################################')
say('### climate/biogeography/management DENSITY models ###')
say('######################################################')

	# Models of pika density at sites where pika were detected. Models use climate, biogeographic (isolation, patch size), and management (perimeter burned, grazing, grass/forb) variables.

	### function to create a "long" table with variable, coefficients, and AICc

		tallyVariableImp <- function(varImport, model) {
		
			coeffs <- coefficients(model)
			
			thisOut <- data.frame(var = names(coeffs), coefficient = as.numeric(coeffs))
			thisOut$AICc <- AICc(model)
			
			varImport <- rbind(varImport, thisOut)
		
		}

	### analyze all regions together and just Jemez alone
	#####################################################

	for (focalRegion in c('all', 'southwest', 'non-southwest')) {

	  	say(focalRegion, level = 2)

		### collate data
		################

			### data
			load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')
			pika <- pika[!is.na(pika$latestDensity), ]
	  
			# Jemez
			if (focalRegion == 'southwest') pika <- pika[pika$region == 'southwest', ]
			if (focalRegion == 'non-southwest') pika <- pika[pika$region != 'southwest', ]

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
		
			# use climate/biogeographic models as basis
			file <- if (extended) {
				paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models ', analysisDate, '/Density - Models without Management Variables ', capIt(focalRegion), ' Region ', analysisDate, '.csv')
			} else {
				paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models ', analysisDate, '/Density - Models without Management Variables ', capIt(focalRegion), ' Region ', analysisDate, '.csv')
			}
			allClimModels <- read.csv(file)
			topClimModels <- allClimModels[allClimModels$deltaAicc < maxDeltaAic_density, ]

			report <- data.frame()
			varImport <- data.frame()
			
			nullModel <- glm(latestDensity ~ 1, data=pika, family=Gamma(link='log'))
			if (focalRegion %in% c('all', 'non-southwest')) {
			  	
			  	nullModelRegion <- glm(latestDensity ~ region, data=pika, family=Gamma(link='log'))
				
			 	nullModelRegionIsolation <- glm(latestDensity ~ region + 
				meanDistToClosest4Patches, data=pika, family=Gamma(link='log'))
			  
			}
	  
			nullModelIsolation <- glm(latestDensity ~ meanDistToClosest4Patches, data=pika, family=Gamma(link='log'))
					
			varImport <- tallyVariableImp(varImport, model = nullModelIsolation)
			if (focalRegion %in% c('all', 'non-southwest')) varImport <- tallyVariableImp(varImport, model = nullModelRegionIsolation)
			
			llNull <- logLik(nullModel)
			if (focalRegion %in% c('all', 'non-southwest')) {
			  	llNullRegion <- logLik(nullModelRegion)
			    llNullRegionIsolation <- logLik(nullModelRegionIsolation)
			}
			llNullIsolation <- logLik(nullModelIsolation)
			
			ecoTermsGrid <- expand.grid(a = c(TRUE, FALSE), b = c(TRUE, FALSE), c = c(TRUE, FALSE))
			names(ecoTermsGrid) <- ecoVars
			
			for (countClimModel in 1:nrow(topClimModels)) {

				# climate formula
				climForm <- topClimModels$model[countClimModel]
				climForm <- paste0('latestDensity ~ ', climForm)
				hasRegion <- focalRegion %in% c('all', 'non-southwest') & topClimModels$region[countClimModel]
				if (hasRegion & focalRegion %in% c('all', 'non-southwest')) climForm <- paste(climForm, ' + region')
				hasIsolation <- !is.na(topClimModels$isolationCoef[countClimModel])
				if (hasIsolation) climForm <- paste0(climForm, ' + meanDistToClosest4Patches')
					
				for (countEcoModel in 1:nrow(ecoTermsGrid)) {

					climMgmtForm <- climForm

					# eco-variable formula
					ecoTermsInModel <- unlist(ecoTermsGrid[countEcoModel, ])
					ecoTermsInModel <- ecoVars[ecoTermsInModel]
					ecoTermsInModel <- paste(ecoTermsInModel, collapse = ' + ')
					if (ecoTermsInModel != '') climMgmtForm <- paste(climMgmtForm, '+', ecoTermsInModel)

					# model
					model <- glm(climMgmtForm, data=pika, family=Gamma(link='log'))
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
					# pseudoR2 <- nagelR2(llNull, ll, n=nrow(pika))				
					pseudoR2 <- 1 - model$deviance / nullModel$deviance
					
					# report
					varImport <- tallyVariableImp(varImport, model = model)
					
					report <- rbind(
						report,
						data.frame(
							model = climMgmtForm,
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
							aiccClimMgmt = aicc,
							pseudoR2Clim = pseudoR2,
							pseudoR2climMgmt = topClimModels$pseudoR2[countClimModel]
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
					aiccClimMgmt = AICc(nullModel),
					pseudoR2Clim = 0,
					pseudoR2climMgmt = NA
				)
			)

	  	if (focalRegion %in% c('all', 'non-southwest')) {
			
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
						aiccClimMgmt = AICc(nullModelRegion),
						# pseudoR2Clim = nagelR2(llNull, llNullRegion, n=nrow(pika)),
						pseudoR2Clim = 1 - nullModelRegion$deviance / nullModel$deviance,
						pseudoR2climMgmt = NA
					)
				)
			
			}

			# isolation-only model
			aiccClim <- allClimModels$aicc[allClimModels$model == '(Intercept)' & !is.na(allClimModels$isolation) & !allClimModels$region]

			report <- rbind(
				report,
				data.frame(
					model = 'latestDensity ~ meanDistToClosest4Patches',
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
					aiccClimMgmt = AICc(nullModelIsolation),
					pseudoR2Clim = 1 - nullModelIsolation$deviance / nullModel$deviance,
					pseudoR2climMgmt = NA
				)
			)
	
	  	if (focalRegion %in% c('all', 'non-southwest')) {
		
			# isolation/region-only model
			aiccClim <- allClimModels$aicc[allClimModels$model == '(Intercept)' & !is.na(allClimModels$isolation) & allClimModels$region]

			report <- rbind(
				report,
				data.frame(
					model = 'latestDensity ~ meanDistToClosest4Patches',
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
					aiccClimMgmt = AICc(nullModelRegionIsolation),
					pseudoR2Clim = 1 - nullModelRegionIsolation$deviance / nullModel$deviance,
					pseudoR2climMgmt = NA
				)
			)
		
		}

		### report
		##########

			# sample size
	  		n <- nrow(pika)
	  
			file <- if (extended) {
				paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models ', analysisDate, '/Density - Models with Management Variables ', capIt(focalRegion), ' Region ', analysisDate, ' Sample Size.txt')
			} else {
				paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models ', analysisDate, '/Density - Models with Management Variables - ', capIt(focalRegion), ' Region ', analysisDate, ' Sample Size.txt')
			}

			sink(file, split = TRUE)
			say('Number of sites with management variables in ', focalRegion, ' region: ', n)
			sink()

			# model performance
			report$deltaAiccClim <- report$aiccClim - min(report$aiccClim)
			w <- exp(-0.5 * report$deltaAiccClim)
			report$weightClim <- w / sum(w)

			report$deltaAiccClimMgmt <- report$aiccClimMgmt - min(report$aiccClimMgmt)
			w <- exp(-0.5 * report$deltaAiccClimMgmt)
			report$weightClimMgmt <- w / sum(w)

			report <- report[order(report$aiccClimMgmt), ]

			rownames(report) <- NULL

			file <- if (extended) {
				paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models ', analysisDate, '/Density - Models with Management Variables ', capIt(focalRegion), ' Region ', analysisDate, '.csv')
			} else {
				paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models ', analysisDate, '/Density - Models with Management Variables - ', capIt(focalRegion), ' Region ', analysisDate, '.csv')
			}
			write.csv(report, file, row.names=FALSE)
			
			# variable importance
			varImport$deltaAicc <- varImport$AICc - min(varImport$AICc)
			w <- exp(-0.5 * varImport$deltaAicc)
			varImport$aiccWeight <- w / sum(w)
			
			vars <- unique(varImport$var)
			vars <- vars[vars %notin% c('(Intercept)', 'regionsouthwest', 'regionnorthwest', 'regionsoutheast')]
		
			varImportOverall <- data.frame()
			for (var in vars) {
			
				subVarImport <- varImport[varImport$var == var, ]
			
				varImportOverall <- rbind(
					varImportOverall,
					data.frame(
						var = var,
						sumAiccWeight = sum(subVarImport$aiccWeight),
						meanAiccWeight = sum(subVarImport$aiccWeight) / nrow(subVarImport),
						avgCoeff = sum(subVarImport$coefficient * subVarImport$aiccWeight) / sum(subVarImport$aiccWeight),
						nModels = nrow(subVarImport),
						nPos = sum(subVarImport$coefficient > 0),
						nNeg = sum(subVarImport$coefficient < 0)
					)
				)
			
			}
		
			file <- if (extended) {
				paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models ', analysisDate, '/Density - Variable Importance in Models with Management Variables - ', capIt(focalRegion), ' Region ', analysisDate, '.csv')
			} else {
				paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models ', analysisDate, '/Density - Variable Importance in Models with Management Variables - ', capIt(focalRegion), ' Region ', analysisDate, '.csv')
			}
			write.csv(varImportOverall, file, row.names=FALSE)

	} # next region

say('###################################################################')
say('### double-checking number of models each variable should be in ###')
say('###################################################################')
	
	vars <- getVars('density')
	formulae <- getFormulaeDens()
	if (extended) {
		formulaeExtended <- getFormulaeDensExtended()
		formulae <- sort(unique(c(formulae, formulaeExtended)))
	}
	
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
	file <- if (extended) {
		paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models ', analysisDate, '/Number of Base Models with Each Variable ', analysisDate, '.csv')
	} else {
		paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models ', analysisDate, '/Number of Base Models with Each Variable ', analysisDate, '.csv')
	}
	write.csv(counts, file, row.names=FALSE)

# say('##################################################################################')
# say('### make response plots for top climate/biogeography/management DENSITY models ###')
# say('##################################################################################')

# 	### models
# 	##########

# 		# all regions
# 		file <- if (extended) {
# 			paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models ', analysisDate, '/Density - Models with Management Variables All Region ', analysisDate, '.csv')
# 		} else {
# 			paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models ', analysisDate, '/Density - Models with Management Variables All Region ', analysisDate, '.csv')
# 		}
# 		allModelsReport <- read.csv(file)
# 		allModelsReport <- allModelsReport[order(allModelsReport$aiccClimMgmt), ]
# 		topAllModel <- allModelsReport[1, ]

# 		# Jemez region
# 		file <- if (extended) {
# 			paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models ', analysisDate, '/Density - Models with Management Variables Southwest Region ', analysisDate, '.csv')
# 		} else {
# 			paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models ', analysisDate, '/Density - Models with Management Variables Southwest Region ', analysisDate, '.csv')
# 		}
# 		southwestModelsReport <- read.csv(file)
# 		southwestModelsReport <- southwestModelsReport[order(southwestModelsReport$aiccClimMgmt), ]
# 		topSouthwestModel <- southwestModelsReport[1, ]

# 		# all non-Jemez regions
# 		file <- if (extended) {
# 			paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models ', analysisDate, '/Density - Models with Management Variables Non-southwest Region ', analysisDate, '.csv')
# 		} else {
# 			paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models ', analysisDate, '/Density - Models with Management Variables Non-southwest Region ', analysisDate, '.csv')
# 		}
# 		nonSouthwestModelsReport <- read.csv(file)
# 		nonSouthwestModelsReport <- nonSouthwestModelsReport[order(nonSouthwestModelsReport$aiccClimMgmt), ]
# 		topNonSouthwestModel <- nonSouthwestModelsReport[1, ]

# 	### pika data
# 	#############
		
# 		load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')

# 		pika <- pika[!is.na(pika$latestDensity), ]
# 		pika$region <- as.factor(pika$region)
# 		vars <- getVars('density')
# 		scaledVars <- scale(pika[ , vars])
# 		centers <- attr(scaledVars, 'scaled:center')
# 		scales <- attr(scaledVars, 'scaled:scale')
# 		pika[ , vars] <- scaledVars

# 		### collate "latest" version of each predictor
# 		ecoVars <- c('perimBurnedMostRecent_perc', 'latestGrazing', 'latestGrassForb') # NB grass, forbs, and grass/form are highly correlated

# 		pika$latestGrazing <- pika$latestGrassForb <- NA
		
# 		for (i in 1:nrow(pika)) {
		
# 			pika$latestGrazing[i] <- pika[i, paste0('grazing', pika$latestDensSurveyYear[i])]
# 			pika$latestGrassForb[i] <- pika[i, paste0('grassForb', pika$latestDensSurveyYear[i], '_perc')]
		
		
# 		}
		
# 		# remove any records with NAs in predictors
# 		nas <- which(
# 			is.na(pika$latestGrazing) | is.na(pika$perimBurnedMostRecent_perc) | is.na(pika$latestGrassForb)
# 		)

# 		if (length(nas) > 0) pika <- pika[-nas, ]

# 		ecoScales <- scale(pika[ , ecoVars[ecoVars != 'latestGrazing']])
# 		centers <- c(centers, attr(ecoScales, 'scaled:center'))
# 		scales <- c(scales, attr(ecoScales, 'scaled:scale'))
# 		pika[ , ecoVars[ecoVars != 'latestGrazing']] <- ecoScales
# 		pika$meanDistToClosest4Patches <- log10(pika$meanDistToClosest4Patches)
# 		centers <- c(centers, 'meanDistToClosest4Patches' = mean(pika$meanDistToClosest4Patches))
# 		scales <- c(scales, 'meanDistToClosest4Patches' = sd(pika$meanDistToClosest4Patches))
# 		pika$meanDistToClosest4Patches <- as.numeric(scale(pika$meanDistToClosest4Patches))
		
# 		pika$region <- as.factor(pika$region)

# 		allModel <- glm(as.formula(topAllModel$model), data = pika, family = Gamma(link='log'))
# 		southwestModel <- glm(as.formula(topSouthwestModel$model), data = pika, family = Gamma(link='log'))
# 		nonSouthwestModel <- glm(as.formula(topNonSouthwestModel$model), data = pika, family = Gamma(link='log'))

# 	### plot each term
# 	##################
		
# 		topAllVars <- attr(terms(allModel), 'term.labels')
# 		topSouthwestVars <- attr(terms(southwestModel), 'term.labels')
# 		topNonSouthwestVars <- attr(terms(nonSouthwestModel), 'term.labels')
# 		topVars <- unique(c(topAllVars, topSouthwestVars, topNonSouthwestVars))

# 		if (any(grepl(topVars, pattern = ':'))) topVars <- topVars[!grepl(topVars, pattern = ':')]

# 		for (var in c('perimBurnedMostRecent_perc', 'latestGrazing', 'meanDistToClosest4Patches')) {

# 			if (any(topVars == var)) {
# 				index <- which(topVars == var)
# 				topVars <- if (index == 1) {
# 					c(topVars[2:length(topVars)], topVars[index])
# 				} else {
# 					c(topVars[1:(index - 1)], topVars[(index + 1):length(topVars)], topVars[index])
# 				}
# 			}

# 		}

# 		baseData <- pika

# 		ymax <- -Inf
# 		plots <- list()
# 		for (i in seq_along(topVars)) {

# 			var <- topVars[i]

# 			if (var != 'latestGrazing') {
				
# 				plotData <- data.frame(
# 					var = seq(min(pika[ , var]), max(pika[ , var]), length.out = 100)
# 				)

# 			} else {

# 				plotData <- data.frame(
# 					var = c(0, 1)
# 				)

# 			}

# 			names(plotData) <- var
# 			nonFocalVars <- topVars[topVars != var]
# 			for (j in seq_along(nonFocalVars)) {
# 				plotData$DUMMY <- 0
# 				names(plotData)[ncol(plotData)] <- nonFocalVars[j]
# 			}

# 			# all-regions model

# 				preds <- predict(allModel, newdata = plotData, type = 'link', se.fit = TRUE)
# 				plotDataAll <- plotData
# 				plotDataAll$predMean <- preds$fit
# 				plotDataAll$predMeanPlus2SE <- preds$fit + 2 * preds$se.fit
# 				plotDataAll$predMeanMinus2SE <- preds$fit - 2 * preds$se.fit

# 				plotDataAll$predMean <- exp(plotDataAll$predMean)
# 				plotDataAll$predMeanPlus2SE <- exp(plotDataAll$predMeanPlus2SE)
# 				plotDataAll$predMeanMinus2SE <- exp(plotDataAll$predMeanMinus2SE)
			
# 			# non-southwest model

# 				preds <- predict(nonSouthwestModel, newdata = plotData, type = 'link', se.fit = TRUE)
# 				plotDataNonSouthwest <- plotData
# 				plotDataNonSouthwest$predMean <- preds$fit
# 				plotDataNonSouthwest$predMeanPlus2SE <- preds$fit + 2 * preds$se.fit
# 				plotDataNonSouthwest$predMeanMinus2SE <- preds$fit - 2 * preds$se.fit

# 				plotDataNonSouthwest$predMean <- exp(plotDataNonSouthwest$predMean)
# 				plotDataNonSouthwest$predMeanPlus2SE <- exp(plotDataNonSouthwest$predMeanPlus2SE)
# 				plotDataNonSouthwest$predMeanMinus2SE <- exp(plotDataNonSouthwest$predMeanMinus2SE)
			
# 			# southwest model

# 				preds <- predict(southwestModel, newdata = plotData, type = 'link', se.fit = TRUE)
# 				plotDataSouthwest <- plotData
# 				plotDataSouthwest$predMean <- preds$fit
# 				plotDataSouthwest$predMeanPlus2SE <- preds$fit + 2 * preds$se.fit
# 				plotDataSouthwest$predMeanMinus2SE <- preds$fit - 2 * preds$se.fit

# 				plotDataSouthwest$predMean <- exp(plotDataSouthwest$predMean)
# 				plotDataSouthwest$predMeanPlus2SE <- exp(plotDataSouthwest$predMeanPlus2SE)
# 				plotDataSouthwest$predMeanMinus2SE <- exp(plotDataSouthwest$predMeanMinus2SE)
			
# 			# unscale for plotting
# 			if (var != 'latestGrazing') {
# 				plotDataAll$varUnscaled <- plotData[ , var] * scales[var] + centers[var]
# 				plotDataNonSouthwest$varUnscaled <- plotData[ , var] * scales[var] + centers[var]
# 				plotDataSouthwest$varUnscaled <- plotData[ , var] * scales[var] + centers[var]
# 			} else {
# 				plotDataAll$varUnscaled <- plotData[ , var]
# 				plotDataNonSouthwest$varUnscaled <- plotData[ , var]
# 				plotDataSouthwest$varUnscaled <- plotData[ , var]
# 			}

# 			plotDataAllSE <- data.frame(
# 				var = c(plotDataAll$varUnscaled, rev(plotDataAll$varUnscaled)),
# 				se = c(plotDataAll$predMeanPlus2SE, rev(plotDataAll$predMeanMinus2SE))
# 			)

# 			plotDataNonSouthwestSE <- data.frame(
# 				var = c(plotDataNonSouthwest$varUnscaled, rev(plotDataNonSouthwest$varUnscaled)),
# 				se = c(plotDataNonSouthwest$predMeanPlus2SE, rev(plotDataNonSouthwest$predMeanMinus2SE))
# 			)

# 			plotDataSouthwestSE <- data.frame(
# 				var = c(plotDataSouthwest$varUnscaled, rev(plotDataSouthwest$varUnscaled)),
# 				se = c(plotDataSouthwest$predMeanPlus2SE, rev(plotDataSouthwest$predMeanMinus2SE))
# 			)

# 			if (var %in% c('latestGrazing', 'perimBurnedMostRecent_perc', 'meanDistToClosest4Patches')) {
# 				if (var == 'latestGrazing') {
# 					xlab <- 'Grazing'
# 					title <- 'Grazing'
# 				} else if (var == 'perimBurnedMostRecent_perc') {
# 					xlab <- 'Perimeter burned (%)'
# 					title <- 'Perimeter burned'
# 				} else if (var == 'meanDistToClosest4Patches') {
# 					xlab <- 'Isolation (log10, m)'
# 					title <- 'Isolation'
# 				}
	
# 			} else {
# 				xlab <- makeNiceVars(var, 'density')
# 				xlab <- paste0(toupper(substr(xlab, 1, 1)), substr(xlab, 2, nchar(xlab)))
# 				title <- xlab
# 				if (xlab == 'Chronic heat (0 yr)') {
# 					xlab <- 'Chronic heat (0 yr, °C)'
# 				} else if (xlab == 'GS precipitation (0 yr)') {
# 					xlab <- 'GS precipitation (0 yr, mm)'
# 				}

# 			}
			
# 			if (var != 'latestGrazing') {

# 				plots[[i]] <- ggplot() +
# 					geom_polygon(data = plotDataAllSE, mapping = aes(x = var, y = se), color = NA, fill = 'gray', alpha = 0.4) +
# 					geom_line(data = plotDataAll, mapping = aes(x = varUnscaled, y = predMean), linewidth = 1.6) +

# 					geom_polygon(data = plotDataNonSouthwestSE, mapping = aes(x = var, y = se), color = NA, fill = '#1b9e77', alpha = 0.4) +
# 					geom_line(data = plotDataNonSouthwest, mapping = aes(x = varUnscaled, y = predMean), color = '#1b9e77', linewidth = 1.6) +

# 					geom_polygon(data = plotDataSouthwestSE, mapping = aes(x = var, y = se), color = NA, fill = '#d95f02', alpha = 0.4) +
# 					geom_line(data = plotDataSouthwest, mapping = aes(x = varUnscaled, y = predMean), color = '#d95f02', linewidth = 1.6) +

# 					labs(x = xlab, y = 'Predicted density') +
# 					ggtitle(paste0(letters[i], ') ', title)) +
# 					theme_bw()

# 			} else {

# 				plotDataAll$grazing <- plotDataAll$varUnscaled == 1
# 				plotDataAllSE$grazing <- plotDataAllSE$var == 1

# 				plotDataNonSouthwest$grazing <- plotDataNonSouthwest$varUnscaled == 1
# 				plotDataNonSouthwestSE$grazing <- plotDataNonSouthwestSE$var == 1

# 				plotDataSouthwest$grazing <- plotDataSouthwest$varUnscaled == 1
# 				plotDataSouthwestSE$grazing <- plotDataSouthwestSE$var == 1

# 				nudge <- 0.2
# 				size <- 2.2
# 				pch <- 3
# 				linewidth <- 1

# 				plots[[i]] <- ggplot() +
					
# 					geom_point(data = plotDataAll, mapping = aes(x = grazing, y = predMean), position = position_nudge(x = -1 * nudge), size = size, pch = pch) +
# 					geom_line(data = plotDataAllSE[!plotDataAllSE$grazing, ], mapping = aes(x = grazing, y = se), position = position_nudge(x = -1 * nudge), linewidth = linewidth) +
# 					geom_line(data = plotDataAllSE[plotDataAllSE$grazing, ], mapping = aes(x = grazing, y = se), position = position_nudge(x = -1 * nudge), linewidth = linewidth) +

# 					geom_point(data = plotDataNonSouthwest, mapping = aes(x = grazing, y = predMean), color = '#1b9e77', size = size, pch = pch) +
# 					geom_line(data = plotDataNonSouthwestSE[!plotDataNonSouthwestSE$grazing, ], mapping = aes(x = grazing, y = se), color = '#1b9e77', linewidth = linewidth) +
# 					geom_line(data = plotDataNonSouthwestSE[plotDataNonSouthwestSE$grazing, ], mapping = aes(x = grazing, y = se), color = '#1b9e77', linewidth = linewidth) +
					
# 					geom_point(data = plotDataSouthwest, mapping = aes(x = grazing, y = predMean), color = '#d95f02', position = position_nudge(x = nudge), size = size, pch = pch) +
# 					geom_line(data = plotDataSouthwestSE[!plotDataSouthwestSE$grazing, ], mapping = aes(x = grazing, y = se), color = '#d95f02', position = position_nudge(x = nudge), linewidth = linewidth) +
# 					geom_line(data = plotDataSouthwestSE[plotDataSouthwestSE$grazing, ], mapping = aes(x = grazing, y = se), color = '#d95f02', position = position_nudge(x = nudge), linewidth = linewidth) +
					
# 					labs(x = xlab, y = 'Predicted density') +
# 					ggtitle(paste0(letters[i], ') ', title)) +
# 					theme_bw()

# 			}

# 			ymax <- max(ymax, quantile(plotDataAllSE$se, 0.75), quantile(plotDataNonSouthwestSE$se, 0.75), quantile(plotDataSouthwestSE$se, 0.75))

# 		}

# 	for (i in seq_along(plots)) {
# 		plots[[i]] <- plots[[i]] +
# 			coord_cartesian(ylim = c(0, ymax)) +
# 			theme(
# 				axis.title = element_text(size = 12)
# 			)
# 	}
	
# 	combo <- plot_grid(plotlist = plots, nrow = 2, align = 'hv')

# 	file <- if (extended) {
# 		paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models ', analysisDate, '/Density - Response Curves for Top Model ', analysisDate, '.png')
# 	} else {
# 		paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models ', analysisDate, '/Density - Response Curves for Top Model ', analysisDate, '.png')
# 	}
# 	ggsave(combo, filename = file, width = 10, height = 7, bg = 'white', dpi = 600)


say('DONE!!!', level=1, deco='%')
