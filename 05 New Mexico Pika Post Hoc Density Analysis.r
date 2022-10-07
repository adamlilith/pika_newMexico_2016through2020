### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('C:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/05 New Mexico Pika Post Hoc Density Analysis.r')
### source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/05 New Mexico Pika Post Hoc Density Analysis.r')
###
### CONTENTS ###
### setup ###
### calculate distance to nearest patch(es) ###
### post hoc density analysis ###

#############
### setup ###
#############

	# drive <- 'C:'
	# drive <- 'D:'
	drive <- 'E:'

	source(paste0(drive, '/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r'))

say('#################################')
say('### post hoc density analysis ###')
say('#################################')

	### collate data
	################

		### data
		load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')
		pika <- pika[!is.na(pika$latestDensity), ]

		### collate "latest" version of each predictor
		pika$latestGrazing <- pika$latestGrass <- pika$latestForb <- pika$latestGrassForb <- NA
		
		for (i in 1:nrow(pika)) {
		
			pika$latestGrazing[i] <- pika[i, paste0('grazing', pika$latestDensSurveyYear[i])]
			pika$latestGrass[i] <- pika[i, paste0('grass', pika$latestDensSurveyYear[i], '_perc')]
			pika$latestForb[i] <- pika[i, paste0('forb', pika$latestDensSurveyYear[i], '_perc')]
			pika$latestGrassForb[i] <- pika[i, paste0('grassForb', pika$latestDensSurveyYear[i], '_perc')]
		
		
		}
		
		# remove any records with NAs in predictors
		nas <- which(
			is.na(pika$latestGrazing) | is.na(pika$perimBurnedMostRecent_perc) | is.na(pika$latestGrass) | is.na(pika$latestForb) | is.na(pika$latestGrassForb) | is.na(pika$latestDensSurveyYear)
		)

		if (length(nas) > 0) pika <- pika[-nas, ]

		# distances to nearest patches
		pika$logMeanDistToClosest3Patches_m <- log10(rowMeans(pika[ , paste0('distClosestPatch_patch', 1:3, '_m')]))
		pika$logMeanDistToClosest4Patches_m <- log10(rowMeans(pika[ , paste0('distClosestPatch_patch', 1:4, '_m')]))
		
	### models
	##########
	
		allClimModels <- read.csv('./Figures & Tables/Density - Simple Models/Density - Simple GLMs.csv')
		topClimModels <- allClimModels[allClimModels$deltaAicc < maxDeltaAic_density, ]

		report <- data.frame()
		
		nullModel <- glm(latestDensity ~ 1, data=pika, family=Gamma(link='log'))
		nullModelRegion <- glm(latestDensity ~ region, data=pika, family=Gamma(link='log'))
		llNull <- logLik(nullModel)
		llNullRegion <- logLik(nullModelRegion)

		ecoTerms <- c('latestGrazing', 'perimBurnedMostRecent_perc', 'latestGrassForb', 'latestDensSurveyYear', 'logMeanDistToClosest3Patches_m', 'logMeanDistToClosest4Patches_m')
		ecoTermsGrid <- expand.grid(a = c(TRUE, FALSE), b = c(TRUE, FALSE), c = c(TRUE, FALSE), d = c(TRUE, FALSE), e = c(TRUE, FALSE), f= c(TRUE, FALSE))
		names(ecoTermsGrid) <- ecoTerms
		
		ecoTermsGrid <- ecoTermsGrid[!(ecoTermsGrid$logMeanDistToClosest3Patches_m & ecoTermsGrid$logMeanDistToClosest4Patches_m), ]

		for (countClimModel in 1:nrow(topClimModels)) {

			# climate formula
			climForm <- topClimModels$model[countClimModel]
			climForm <- paste0('latestDensity ~ ', climForm)
			hasRegion <- topClimModels$region[countClimModel]
			if (hasRegion) climForm <- paste(climForm, ' + region')
				
			for (countEcoModel in 1:nrow(ecoTermsGrid)) {

				form <- climForm

				# eco-variable formula
				ecoTermsInModel <- unlist(ecoTermsGrid[countEcoModel, ])
				ecoTermsInModel <- ecoTerms[ecoTermsInModel]
				ecoTermsInModel <- paste(ecoTermsInModel, collapse = ' + ')
				if (ecoTermsInModel != '') form <- paste(form, '+', ecoTermsInModel)

				formText <- form
				form <- as.formula(form)

				# center and scale
				terms <- attr(terms(form), 'term.labels')
				thisPika <- pika
				thisPika <- thisPika[ , c('latestDensity', 'region', terms)]
				for (term in terms) {
					if (!(term %in% c('latestGrazing', 'region'))) thisPika[ , term] <- scale(thisPika[ , term, drop=TRUE])
				}
				
				if (hasRegion) thisPika$region <- as.factor(thisPika$region)

				# model
				model <- glm(form, data=thisPika, family=Gamma(link='log'))
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
				pseudoR2 <- nagelR2(llNull, ll, n=nrow(thisPika))

				# report
				report <- rbind(
					report,
					data.frame(
						model = formText,
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
						aicc = aicc,
						pseudoR2 = pseudoR2
					)
				)
				
			} # next eco-term model
		
		} # next climate model
		
		# intercept-only model
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
				aicc = AICc(nullModel),
				pseudoR2 = 0
			)
		)

		# region-only model
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
				aicc = AICc(nullModelRegion),
				pseudoR2 = nagelR2(llNull, llNullRegion, n=nrow(thisPika))
			)
		)

	### report
	##########
	
		report$deltaAicc <- report$aicc - min(report$aicc)
		w <- exp(-0.5 * report$deltaAicc)
		report$weight <- w / sum(w)
		report <- report[order(report$aicc), ]
		rownames(report) <- NULL

		dirCreate('./Figures & Tables/Density - Post Hoc Models')
		file <- paste0('./Figures & Tables/Density - Post Hoc Models/Post Hoc Density Models.csv')
		write.csv(report, file, row.names=FALSE)
		
say('DONE!!!', level=1, deco='%')
