### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### Analysis of sites where non-zero density of pikas was recorded. Constructs and compares multiple simple models of density and estimates of variable importance. Exactly the same as the other "03b" script, except that it uses a kriging term to "absorb" spatial variation not explained by the predictor(s) in each model. The output folder and outputs have "Accounting for SAC" appended to their names. The kriging is implemented using gam(y ~ x + s(long, lat, bs = 'gp'))
###
### source('C:/Kaji/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/03b New Mexico Pika Occupancy & Abundance Analysis - Density with Simple Models - Accounting for Spatial Autocorrelation.r')
###
### CONTENTS ###
### setup ###
### climate/biogeography DENSITY models ###
### compile table of predictor importance for climate/biogeography DENSITY models ###

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
		dirCreate('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models ', analysisDate, ' - Accounting for SAC')
	} else {
		say('Using ORIGINAL set of climate models!!!', level = 2)
		dirCreate('./Figures & Tables/Density - Simple Models with Original Set of Climate Models ', analysisDate, ' - Accounting for SAC')
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

			model1 <- gam(latestDensity ~ 1 + s(longitude, latitude, bs = 'gp'), data=pika, family=Gamma(link='log'))
			if (focalRegion != 'southwest') model2 <- gam(latestDensity ~ 1 + region + s(longitude, latitude, bs = 'gp'), data=pika, family=Gamma(link='log'))
			model3 <- gam(latestDensity ~ meanDistToClosest4Patches + s(longitude, latitude, bs = 'gp'), data=pika, family=Gamma(link='log'))
			if (focalRegion != 'southwest') model4 <- gam(latestDensity ~ 1 + region + meanDistToClosest4Patches + s(longitude, latitude, bs = 'gp'), data=pika, family=Gamma(link='log'))
			
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

				form1 <- as.formula(paste0('latestDensity ~ 1 + ', formula, ' + s(longitude, latitude, bs = "gp")'))
				if (focalRegion != 'southwest') form2 <- as.formula(paste0('latestDensity ~ 1 + ', formula, ' + region', ' + s(longitude, latitude, bs = "gp")'))
				form3 <- as.formula(paste0('latestDensity ~ 1 + ', formula, ' + meanDistToClosest4Patches', ' + s(longitude, latitude, bs = "gp")'))
				if (focalRegion != 'southwest') form4 <- as.formula(paste0('latestDensity ~ 1 + ', formula, ' + region + meanDistToClosest4Patches', ' + s(longitude, latitude, bs = "gp")'))
		
				model1 <- gam(form1, data=pika, family=Gamma(link='log'))
				if (focalRegion != 'southwest') model2 <- gam(form2, data=pika, family=Gamma(link='log'))
				model3 <- gam(form3, data=pika, family=Gamma(link='log'))
				if (focalRegion != 'southwest') model4 <- gam(form4, data=pika, family=Gamma(link='log'))
				
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
					paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models ', analysisDate, ' - Accounting for SAC/Density - Models without Management Variables ', capIt(focalRegion), ' Region ', analysisDate, ' - Accounting for SAC.csv')
				} else {
					paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models ', analysisDate, ' - Accounting for SAC/Density - Models without Management Variables ', capIt(focalRegion), ' Region ', analysisDate, ' - Accounting for SAC.csv')
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
			paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models ', analysisDate, ' - Accounting for SAC/Density - Models without Management Variables ', capIt(focalRegion), ' Region ', analysisDate, ' - Accounting for SAC.csv')
		} else {
			paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models ', analysisDate, ' - Accounting for SAC/Density - Models without Management Variables ', capIt(focalRegion), ' Region ', analysisDate, ' - Accounting for SAC.csv')
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
			paste0('./Figures & Tables/Density - Simple Models with Extended Set of Climate Models ', analysisDate, ' - Accounting for SAC/Density - Variable Importance in Models without Management Variables ', capIt(focalRegion), ' Region ', analysisDate, ' - Accounting for SAC.csv')
		} else {
			paste0('./Figures & Tables/Density - Simple Models with Original Set of Climate Models ', analysisDate, ' - Accounting for SAC/Density - Variable Importance in Models without Management Variables ', capIt(focalRegion), ' Region ', analysisDate, ' - Accounting for SAC.csv')
		}
		write.csv(imp, file, row.names=FALSE)

	}

say('DONE!!!', level=1, deco='%')
