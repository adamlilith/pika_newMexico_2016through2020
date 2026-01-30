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

	if (extended) {
		say('Using EXTENDED set of climate models!!!', level = 2)
		dirCreate('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models')
	} else {
		say('Using ORIGINAL set of climate models!!!', level = 2)
		dirCreate('./Figures & Tables/Density - Simple Models with Original Set of Climate Models')
	}

say('###########################################')
say('### climate/biogeography DENSITY models ###')
say('###########################################')

	# Models of pika density at sites where pika were detected. Models use only climate and biogeographic (isolation, patch size) variables.

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
		if (extended) {
		
			formulaeExtended <- getFormulaeDensExtended()
			formulae <- sort(unique(c(formulae, formulaeExtended)))

		}

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

			file <- if (extended) {
				paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models/Density - Models without Management Variables.csv')
			} else {
				paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models/Density - Models without Management Variables.csv')
			}
			write.csv(results, file, row.names=FALSE)
	
say('#####################################################################################')
say('### compile table of predictor importance for climate/biogeography DENSITY models ###')
say('#####################################################################################')
	
	# Assess importance of variables in models of pika density at sites where pika were detected. Models use only climate and biogeographic (isolation, patch size) variables.

	# get models
	file <- if (extended) {
		paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models/Density - Models without Management Variables.csv')
	} else {
		paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models/Density - Models without Management Variables.csv')
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
		paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models/Density - Variable Importance in Models without Management Variables.csv')
	} else {
		paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models/Density - Variable Importance in Models without Management Variables.csv')
	}
	write.csv(imp, file, row.names=FALSE)
	
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
	
		# use climate/biogeographic models as basis
		file <- if (extended) {
			'./Figures & Tables/Density - Simple Models with Extended Set of Climate Models/Density - Models without Management Variables.csv'
		} else {
			'./Figures & Tables/Density - Simple Models with Original Set of Climate Models/Density - Models without Management Variables.csv'
		}
		allClimModels <- read.csv(file)
		topClimModels <- allClimModels[allClimModels$deltaAicc < maxDeltaAic_density, ]

		report <- data.frame()
		varImport <- data.frame()
		
		nullModel <- glm(latestDensity ~ 1, data=pika, family=Gamma(link='log'))
		nullModelRegion <- glm(latestDensity ~ region, data=pika, family=Gamma(link='log'))
		nullModelIsolation <- glm(latestDensity ~ meanDistToClosest4Patches, data=pika, family=Gamma(link='log'))
		nullModelRegionIsolation <- glm(latestDensity ~ region + meanDistToClosest4Patches, data=pika, family=Gamma(link='log'))
		
		varImport <- tallyVariableImp(varImport, model = nullModelIsolation)
		varImport <- tallyVariableImp(varImport, model = nullModelRegionIsolation)
		
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
				pseudoR2 <- nagelR2(llNull, ll, n=nrow(pika))				
				
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
				pseudoR2Clim = nagelR2(llNull, llNullRegion, n=nrow(pika)),
				pseudoR2climMgmt = NA
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
				aiccClimMgmt = AICc(nullModelIsolation),
				pseudoR2Clim = nagelR2(llNull, llNullIsolation, n=nrow(pika)),
				pseudoR2climMgmt = NA
			)
		)

		# isolation/region-only model
		aiccClim <- allClimModels$aicc[allClimModels$model == '(Intercept)' & !is.na(allClimModels$isolation) & allClimModels$region]

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
				aiccClimMgmt = AICc(nullModelRegionIsolation),
				pseudoR2Clim = nagelR2(llNull, llNullRegionIsolation, n=nrow(pika)),
				pseudoR2climMgmt = NA
			)
		)

	### report
	##########
	
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
			paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models/Density - Models with Management Variables.csv')
		} else {
			paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models/Density - Models with Management Variables.csv')
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
			paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models/Density - Variable Importance in Models with Management Variables.csv')
		} else {
			paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models/Density - Variable Importance in Models with Management Variables.csv')
		}
		write.csv(varImportOverall, file, row.names=FALSE)
			
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
		'./Figures & Tables/Density - Simple Models with Extended Set of Climate Models/Number of Base Models with Each Variable.csv'
	} else {
		'./Figures & Tables/Density - Simple Models with Original Set of Climate Models/Number of Base Models with Each Variable.csv'
	}
	write.csv(counts, file, row.names=FALSE)
	
say('DONE!!!', level=1, deco='%')
