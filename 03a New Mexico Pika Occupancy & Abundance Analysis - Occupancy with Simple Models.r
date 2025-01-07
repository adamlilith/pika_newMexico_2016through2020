### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('C:/Ecology/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/03a New Mexico Pika Occupancy & Abundance Analysis - Occupancy with Simple Models.r')
### source('C:/Subarashi/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/03a New Mexico Pika Occupancy & Abundance Analysis - Occupancy with Simple Models.r')
###
### CONTENTS ###
### setup ###
###
### ORDINAL simple OCCUPANCY analysis ###
### report relative odds of class change across regions across all ORDINAL OCCUPANCY models ###
### summarize support for 7- and 10-yr windows ###
### compile table of predictor weights for ORDINAL simple OCCUPANCY analysis: all-data models ###
###
### BINARY simple OCCUPANCY analysis ###
### report relative odds of class change across regions across all BINARY OCCUPANCY models ###
### compile table of predictor weights for BINARY simple OCCUPANCY analysis ###
### double-checking number of models each variable should be in ###
###
### BINARY and ORDINAL analyses ###
### make maps of predicted probabilities of each class from best ordinal and binary models ###
### make maps of future predicted probabilities from binary models and extract predictions at sites ###
### cross-validation of models ###

#############
### setup ###
#############

	rm(list=ls())

	drive <- 'C:/Ecology/'
	# drive <- 'C:/Subarashi/'

	source(paste0(drive, '/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r'))

# say('#########################################')
# say('### ORDINAL simple OCCUPANCY analysis ###')
# say('#########################################')

# 	load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')
# 	pika$latestOccStatus <- factor(pika$latestOccStatus, levels=c('0 never', '1 old', '2 occupied'), ordered=TRUE)
	
# 	pika$region <- as.factor(pika$region)
	
# 	pika$meanDistToClosest4Patches <- log10(pika$meanDistToClosest4Patches)
# 	pika$numHomeRangesScaled <- scale(log10(pika$numHomeRanges + 1))
	
# 	vars <- getVars('occupancy')
# 	vars <- c(vars, 'meanDistToClosest4Patches')
# 	pika[ , vars] <- scale(pika[ , vars])

# 	# models
# 	formulae <- getFormulaeOcc()

# 	# coefficients and model AICc's
# 	vars <- getVars('occupancy')
# 	vars <- c(vars, 'meanDistToClosest4Patches')
# 	accumCoeffs <- list()
# 	for (i in seq_along(vars)) accumCoeffs[[i]] <- numeric()
# 	names(accumCoeffs) <- vars
# 	aiccs <- accumCoeffs

# 	### using all sites (no cross-validation)
# 	#########################################
	
# 		### null models
# 		###############

# 			model1 <- polr(latestOccStatus ~ 1, data=pika, Hess=TRUE)
# 			model2 <- polr(latestOccStatus ~ numHomeRangesScaled, data=pika, Hess=TRUE)
# 			model3 <- polr(latestOccStatus ~ region, data=pika, Hess=TRUE)
# 			model4 <- polr(latestOccStatus ~ numHomeRangesScaled + region, data=pika, Hess=TRUE)

# 			model5 <- polr(latestOccStatus ~ meanDistToClosest4Patches, data=pika, Hess=TRUE)
# 			model6 <- polr(latestOccStatus ~ numHomeRangesScaled + meanDistToClosest4Patches, data=pika, Hess=TRUE)
# 			model7 <- polr(latestOccStatus ~ region + meanDistToClosest4Patches, data=pika, Hess=TRUE)
# 			model8 <- polr(latestOccStatus ~ numHomeRangesScaled + region + meanDistToClosest4Patches, data=pika, Hess=TRUE)

# 			aicc1 <- AICc(model1)
# 			aicc2 <- AICc(model2)
# 			aicc3 <- AICc(model3)
# 			aicc4 <- AICc(model4)
# 			aicc5 <- AICc(model5)
# 			aicc6 <- AICc(model6)
# 			aicc7 <- AICc(model7)
# 			aicc8 <- AICc(model8)

# 			coeffs <- c(NA, NA, NA, NA, NA, NA, NA, NA)
# 			numHomeRanges <- c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)
# 			region <- c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)
# 			isolation <- c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
			
# 			aicc <- c(aicc1, aicc2, aicc3, aicc4, aicc5, aicc6, aicc7, aicc8)

# 			like1 <- logLik(model1)
# 			like2 <- logLik(model2)
# 			like3 <- logLik(model3)
# 			like4 <- logLik(model4)
# 			like5 <- logLik(model5)
# 			like6 <- logLik(model6)
# 			like7 <- logLik(model7)
# 			like8 <- logLik(model8)
			
# 			likeNull <- like1
			
# 			n <- nrow(pika)
# 			pseudoR2_1 <- nagelR2(likeNull, like1, n)
# 			pseudoR2_2 <- nagelR2(likeNull, like2, n)
# 			pseudoR2_3 <- nagelR2(likeNull, like3, n)
# 			pseudoR2_4 <- nagelR2(likeNull, like4, n)
# 			pseudoR2_5 <- nagelR2(likeNull, like5, n)
# 			pseudoR2_6 <- nagelR2(likeNull, like6, n)
# 			pseudoR2_7 <- nagelR2(likeNull, like7, n)
# 			pseudoR2_8 <- nagelR2(likeNull, like8, n)
			
# 			numHomRangesCoef2 <- coefficients(model2)['numHomeRangesScaled']
# 			numHomRangesCoef4 <- coefficients(model4)['numHomeRangesScaled']
# 			numHomRangesCoef6 <- coefficients(model6)['numHomeRangesScaled']
# 			numHomRangesCoef8 <- coefficients(model8)['numHomeRangesScaled']
			
# 			isolationCoef5 <- coefficients(model5)['meanDistToClosest4Patches']
# 			isolationCoef6 <- coefficients(model6)['meanDistToClosest4Patches']
# 			isolationCoef7 <- coefficients(model7)['meanDistToClosest4Patches']
# 			isolationCoef8 <- coefficients(model8)['meanDistToClosest4Patches']

# 			pseudoR2 <- c(pseudoR2_1, pseudoR2_2, pseudoR2_3, pseudoR2_4, pseudoR2_5, pseudoR2_6, pseudoR2_7, pseudoR2_8)

# 			coef1 <- coefficients(model1)
# 			coef2 <- coefficients(model2)
# 			coef3 <- coefficients(model3)
# 			coef4 <- coefficients(model4)
# 			coef5 <- coefficients(model5)
# 			coef6 <- coefficients(model6)
# 			coef7 <- coefficients(model7)
# 			coef8 <- coefficients(model8)

# 			nw3 <- coef3['regionnorthwest']
# 			se3 <- coef3['regionsoutheast']
# 			sw3 <- coef3['regionsouthwest']

# 			nw4 <- coef4['regionnorthwest']
# 			se4 <- coef4['regionsoutheast']
# 			sw4 <- coef4['regionsouthwest']

# 			nw7 <- coef7['regionnorthwest']
# 			se7 <- coef7['regionsoutheast']
# 			sw7 <- coef7['regionsouthwest']

# 			nw8 <- coef8['regionnorthwest']
# 			se8 <- coef8['regionsoutheast']
# 			sw8 <- coef8['regionsouthwest']

# 			results <- data.frame(
# 				model = '(Intercept)',
# 				term1 = NA,
# 				term2 = NA,
# 				term3 = NA,
# 				term4 = NA,
# 				period = NA,
# 				numHomeRanges = numHomeRanges,
# 				homeRangeCoeff = c(NA, numHomRangesCoef2, NA, numHomRangesCoef4, NA, numHomRangesCoef6, NA, numHomRangesCoef8),
# 				isolationCoeff = c(NA, NA, NA, NA, isolationCoef5, isolationCoef6, isolationCoef7, isolationCoef8),
# 				region = region,
# 				aicc = aicc,
# 				pseudoR2 = pseudoR2,
# 				nw = c(NA, NA, nw3, nw4, NA, NA, nw7, nw8),
# 				se = c(NA, NA, se3, se4, NA, NA, se7, se8),
# 				sw = c(NA, NA, sw3, sw4, NA, NA, sw7, sw8)
# 			)
			
# 			# remember coefficients to do AICc-based coefficient weighting
# 			for (k in 5:8) {
			
# 				coefk <- get(paste0('coef', k))
# 				aicck <- get(paste0('aicc', k))

# 				for (i in seq_along(accumCoeffs)) {
# 					for (j in seq_along(coefk)) {
# 						if (names(coefk)[j] == names(accumCoeffs)[i]) {
# 							accumCoeffs[[i]] <- c(accumCoeffs[[i]], coefk[j])
# 							aiccs[[i]] <- c(aiccs[[i]], aicck)
# 						}
# 					}
				
# 				}
				
# 			}

# 		### by climate variable
# 		#######################
		
# 		for (formula in formulae) {
			
# 			say(formula)
			
# 			form1 <- as.formula(paste0('latestOccStatus ~ 1 + ', formula))
# 			form2 <- as.formula(paste0('latestOccStatus ~ 1 + ', formula, ' + numHomeRangesScaled'))
# 			form3 <- as.formula(paste0('latestOccStatus ~ 1 + ', formula, ' + region'))
# 			form4 <- as.formula(paste0('latestOccStatus ~ 1 + ', formula, ' + numHomeRangesScaled + region'))
			
# 			form5 <- as.formula(paste0('latestOccStatus ~ 1 + ', formula, ' + meanDistToClosest4Patches'))
# 			form6 <- as.formula(paste0('latestOccStatus ~ 1 + ', formula, ' + numHomeRangesScaled + meanDistToClosest4Patches'))
# 			form7 <- as.formula(paste0('latestOccStatus ~ 1 + ', formula, ' + region + meanDistToClosest4Patches'))
# 			form8 <- as.formula(paste0('latestOccStatus ~ 1 + ', formula, ' + numHomeRangesScaled + region + meanDistToClosest4Patches'))
			
# 			model1 <- polr(form1, data=pika, Hess=TRUE)
# 			model2 <- polr(form2, data=pika, Hess=TRUE)
# 			model3 <- polr(form3, data=pika, Hess=TRUE)
# 			model4 <- polr(form4, data=pika, Hess=TRUE)
# 			model5 <- polr(form5, data=pika, Hess=TRUE)
# 			model6 <- polr(form6, data=pika, Hess=TRUE)
# 			model7 <- polr(form7, data=pika, Hess=TRUE)
# 			model8 <- polr(form8, data=pika, Hess=TRUE)

# 			aicc1 <- AICc(model1)
# 			aicc2 <- AICc(model2)
# 			aicc3 <- AICc(model3)
# 			aicc4 <- AICc(model4)
# 			aicc5 <- AICc(model5)
# 			aicc6 <- AICc(model6)
# 			aicc7 <- AICc(model7)
# 			aicc8 <- AICc(model8)
	
# 			terms <- extractTerms(model1, model2, model3, model4, model5, model6, model7, model8)
# 			term1 <- terms$term1
# 			term2 <- terms$term2
# 			term3 <- terms$term3
# 			term4 <- terms$term4
# 			term5 <- terms$term5
# 			term6 <- terms$term6
# 			term7 <- terms$term7
# 			term8 <- terms$term8
			
# 			numHomeRanges <- c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)
# 			region <- c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)
# 			isolation <- c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
			
# 			aicc <- c(aicc1, aicc2, aicc3, aicc4, aicc5, aicc6, aicc7, aicc8)
			
# 			like1 <- logLik(model1)
# 			like2 <- logLik(model2)
# 			like3 <- logLik(model3)
# 			like4 <- logLik(model4)
# 			like5 <- logLik(model5)
# 			like6 <- logLik(model6)
# 			like7 <- logLik(model7)
# 			like8 <- logLik(model8)
			
# 			pseudoR2_1 <- nagelR2(likeNull, like1, n)
# 			pseudoR2_2 <- nagelR2(likeNull, like2, n)
# 			pseudoR2_3 <- nagelR2(likeNull, like3, n)
# 			pseudoR2_4 <- nagelR2(likeNull, like4, n)
# 			pseudoR2_5 <- nagelR2(likeNull, like5, n)
# 			pseudoR2_6 <- nagelR2(likeNull, like6, n)
# 			pseudoR2_7 <- nagelR2(likeNull, like7, n)
# 			pseudoR2_8 <- nagelR2(likeNull, like8, n)
			
# 			numHomRangesCoef2 <- coefficients(model2)['numHomeRangesScaled']
# 			numHomRangesCoef4 <- coefficients(model4)['numHomeRangesScaled']

# 			numHomRangesCoef6 <- coefficients(model6)['numHomeRangesScaled']
# 			numHomRangesCoef8 <- coefficients(model8)['numHomeRangesScaled']

# 			isolationCoef5 <- coefficients(model5)['meanDistToClosest4Patches']
# 			isolationCoef6 <- coefficients(model6)['meanDistToClosest4Patches']
# 			isolationCoef7 <- coefficients(model7)['meanDistToClosest4Patches']
# 			isolationCoef8 <- coefficients(model8)['meanDistToClosest4Patches']
			
# 			pseudoR2 <- c(pseudoR2_1, pseudoR2_2, pseudoR2_3, pseudoR2_4, pseudoR2_5, pseudoR2_6, pseudoR2_7, pseudoR2_8)

# 			coef1 <- coefficients(model1)
# 			coef2 <- coefficients(model2)
# 			coef3 <- coefficients(model3)
# 			coef4 <- coefficients(model4)
# 			coef5 <- coefficients(model5)
# 			coef6 <- coefficients(model6)
# 			coef7 <- coefficients(model7)
# 			coef7 <- coefficients(model8)

# 			nw3 <- coef3['regionnorthwest']
# 			se3 <- coef3['regionsoutheast']
# 			sw3 <- coef3['regionsouthwest']
	
# 			nw4 <- coef4['regionnorthwest']
# 			se4 <- coef4['regionsoutheast']
# 			sw4 <- coef4['regionsouthwest']

# 			nw7 <- coef7['regionnorthwest']
# 			se7 <- coef7['regionsoutheast']
# 			sw7 <- coef7['regionsouthwest']

# 			nw8 <- coef8['regionnorthwest']
# 			se8 <- coef8['regionsoutheast']
# 			sw8 <- coef8['regionsouthwest']
			
# 			period <- if (grepl(formula, pattern=paste0(7, 'yrWindow'))) {
# 				'7-yr'
# 			} else {
# 				'10-yr'
# 			}
			
# 			# remember
# 			results <- rbind(
# 				results,
# 				data.frame(
# 					model = formula,
# 					term1 = term1,
# 					term2 = term2,
# 					term3 = term3,
# 					term4 = term4,
# 					period = period,
# 					numHomeRanges = numHomeRanges,
# 					homeRangeCoeff = c(NA, numHomRangesCoef2, NA, numHomRangesCoef4, NA, numHomRangesCoef6, NA, numHomRangesCoef8),
# 					isolationCoeff = c(NA, NA, NA, NA, isolationCoef5, isolationCoef6, isolationCoef7, isolationCoef8),
# 					region = region,
# 					aicc = aicc,
# 					pseudoR2 = pseudoR2,
# 					nw = c(NA, NA, nw3, nw4, NA, NA, nw7, nw8),
# 					se = c(NA, NA, se3, se4, NA, NA, se7, se8),
# 					sw = c(NA, NA, sw3, sw4, NA, NA, sw7, sw8)
# 				)
# 			)
			
# 			# remember coefficients to do AICc-based coefficient weighting
# 			for (k in 1:8) {
			
# 				coefk <- get(paste0('coef', k))
# 				aicck <- get(paste0('aicc', k))

# 				for (i in seq_along(accumCoeffs)) {
# 					for (j in seq_along(coefk)) {
# 						if (names(coefk)[j] == names(accumCoeffs)[i]) {
# 							accumCoeffs[[i]] <- c(accumCoeffs[[i]], coefk[j])
# 							aiccs[[i]] <- c(aiccs[[i]], aicck)
# 						}
# 					}
				
# 				}
				
# 			}

# 		} # next variable

# 		### reports
# 		###########
		
# 			for (occWindow in c(occWindows_y, NA)) {
			
# 				if (is.na(occWindow)) {
# 					thisResults <- results
# 					nice <- paste0(paste(occWindows_y, collapse=' & '), '-yr Windows')
# 				} else {

# 					thisResults <- rbind(
# 						results[grepl(results$model, pattern=paste0(occWindow, 'yrWindow')), ],
# 						results[results$model == '(Intercept)', ]
# 					)
						
# 					nice <- paste0(occWindow, '-yr Window')
# 				}
				
# 				thisResults$deltaAicc <- thisResults$aicc - min(thisResults$aicc)
# 				w <- exp(-0.5 * thisResults$deltaAicc)
# 				thisResults$weight <- w / sum(w)

# 				thisResults <- thisResults[order(thisResults$weight, decreasing=TRUE), ]
# 				rownames(thisResults) <- NULL

# 				file <- paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Ordinal Models Using All Data - ', nice, '.csv')
# 				write.csv(thisResults, file, row.names=FALSE)
				
# 			} # next window

# 	##############################
# 	### summarize coefficients ###
# 	##############################
	
# 	### AICc-weighted
	
# 	minAicc <- min(results$aicc)
# 	deltaAicc <- results$aicc - minAicc
# 	w <- exp(-0.5 * deltaAicc)
# 	wSum <- sum(w)

# 	for (i in seq_along(accumCoeffs)) {
# 		aiccs[[i]] <- aiccs[[i]] - minAicc
# 		aiccs[[i]] <- exp(-0.5 * aiccs[[i]])
# 		aiccs[[i]] <- aiccs[[i]] / wSum
# 		accumCoeffs[[i]] <- accumCoeffs[[i]] * aiccs[[i]] / sum(aiccs[[i]])
# 	}
	
# 	accumCoeffs <- lapply(accumCoeffs, sum)

# 	sink('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Ordinal Models Using All Data - 7- & 10-yr Window Coefficient Summary.txt')
		
# 		say('AICc-weighted coefficient values for ORDINAL OCCUPANCY models', post=2)
# 		print(accumCoeffs)
		
# 	sink()

# say('###############################################################################################')
# say('### report relative odds of class change across regions across all ORDINAL OCCUPANCY models ###')
# say('###############################################################################################')

# 	file <- paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Ordinal Models Using All Data - 7 & 10-yr Windows.csv')
# 	results <- read.csv(file)
# 	results <- results[results$region, ]
# 	deltaAicc <- results$aicc - min(results$aicc)
# 	w <- exp(-0.5 * deltaAicc)
# 	w <- w / sum(w)
# 	nw <- sum(results$nw * w, na.rm=TRUE)
# 	se <- sum(results$se * w, na.rm=TRUE)
# 	sw <- sum(results$sw * w, na.rm=TRUE)
	
# 	sink('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Ordinal Models Using All Data Odds of Switching between Classes.txt')
		
# 		say('Ordinal odds ratios (arithmetic scale):')
		
# 		say('Ordinal (arithmetic) odds of NW region relative to NE region: ', exp(nw))
# 		say('Ordinal (arithmetic) odds of SE region relative to NE region: ', exp(se))
# 		say('Ordinal (arithmetic) odds of SW region relative to NE region: ', exp(sw))
		
# 	sink()

# say('##################################################')
# say('### summarize support for 7- and 10-yr windows ###')
# say('##################################################')

# 	file <- paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Ordinal Models Using All Data - 7 & 10-yr Windows.csv')
# 	results <- read.csv(file)
# 	results <- results[results$region, ]
# 	deltaAicc <- results$aicc - min(results$aicc)
# 	w <- exp(-0.5 * deltaAicc)
# 	w <- w / sum(w)

# 	say('Relative support for 7- vs 10-yr windows for climatic variables:')
# 	say(date(), post=1)
	
# 	aicc7yr <- sum(w[results$period == '7-yr'], na.rm=TRUE)
# 	aicc10yr <- sum(w[results$period == '10-yr'], na.rm=TRUE)
# 	say('Sum of AICc across all models using 7-yr period:  ', aicc7yr)
# 	say('Sum of AICc across all models using 10-yr period: ', aicc10yr)

# say('#################################################################################################')
# say('### compile table of predictor weights for ORDINAL simple OCCUPANCY analysis: all-data models ###')
# say('#################################################################################################')

	# # rank variables by mean AICc weight
	# for (occWindow in c(occWindows_y, NA)) {
	
		# nice <- if (is.na(occWindow)) {
			# paste0(paste(occWindows_y, collapse=' & '), '-yr Windows')
		# } else {
			# paste0(occWindow, '-yr Window')
		# }
	
		# # get variables
		# vars <- getVars('occupancy')
		# if (!is.na(occWindow)) vars <- vars[grepl(vars, pattern=paste0(occWindow, 'yrWindow'))]
		# vars <- c(vars, 'meanDistToClosest4Patches')

		# # get models
		# file <- paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Ordinal Models Using All Data - ', nice, '.csv')
		# models <- read.csv(file)
		
		# imp <- data.frame()
		# for (var in vars) {
		
			# index <- if (var != 'meanDistToClosest4Patches') {
				# which(grepl(models$model, pattern=var))
			# } else {
				# which(!is.na(models$isolationCoeff))
			# }

			# n <- length(index)
			# sumWeight <- sum(models$weight[index])
			# meanWeight <- sumWeight / n

			# # tally how many times this variable has + or - coefficient
			# coeffs <- models[index, c('term1', 'term2', 'term3', 'term4')]
			# if (var != 'meanDistToClosest4Patches') {

				# numPos <- numNeg <- 0
				# for (countIndex in seq_along(index)) {
					
					# thisIndex <- index[countIndex]
					# terms <- strsplit(models$model[thisIndex], split = ' \\+ ')[[1]]
					# whichOne <- which(grepl(terms, pattern = var))
					
					# # do not count if variable is part of an interaction term
					# if (!grepl(terms[whichOne], pattern = '\\*')) {
						
						# coeff <- coeffs[countIndex, whichOne]
						# coeff <- strsplit(coeff, split = ' ')[[1]]
						# coeff <- coeff[1]
						# coeff <- as.numeric(coeff)
						
						# if (coeff < 0) { numNeg <- numNeg + 1 } else { numPos <- numPos + 1 }
						
					# }
						
				# }

			# } else {
				# numNeg <- sum(models[index, 'isolationCoeff'] < 0)
				# numPos <- sum(models[index, 'isolationCoeff'] > 0)
			# }

			# imp <- rbind(
				# imp,
				# data.frame(
					# variable = var,
					# niceVar = makeNiceVars(var, 'occupancy'),
					# numModels = n,
					# sumWeight = sumWeight,
					# meanWeight = meanWeight,
					# numPos = numPos,
					# numNeg = numNeg
				# )
			# )
			
		# }
		
		# imp <- imp[order(imp$meanWeight, decreasing=TRUE), ]

		# imp$niceVar[imp$variable == 'meanDistToClosest4Patches'] <- 'isolation'

		# write.csv(imp, paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Ordinal Models Using All Data - ', nice, ' - Var Import.csv'), row.names=FALSE)
	
	# } # next occupancy window

# say('########################################')
# say('### BINARY simple OCCUPANCY analysis ###')
# say('########################################')

# 	load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')
# 	pika$region <- as.factor(pika$region)
	
# 	pika$meanDistToClosest4Patches <- log10(pika$meanDistToClosest4Patches)
# 	pika$numHomeRangesScaled <- scale(log10(pika$numHomeRanges + 1))
	
# 	vars <- getVars('occupancy')
# 	vars <- c(vars, 'meanDistToClosest4Patches')
# 	pika[ , vars] <- scale(pika[ , vars])

# 	# models
# 	formulae <- getFormulaeOcc()

# 	# coefficients and model AICc's
# 	accumCoeffs <- list()
# 	for (i in seq_along(vars)) accumCoeffs[[i]] <- numeric()
# 	names(accumCoeffs) <- vars
# 	aiccs <- accumCoeffs

# 	### null models
# 	###############

# 		model1 <- glm(presAbs ~ 1, data=pika, family=binomial)
# 		model2 <- glm(presAbs ~ numHomeRangesScaled, data=pika, family=binomial)
# 		model3 <- glm(presAbs ~ region, data=pika, family=binomial)
# 		model4 <- glm(presAbs ~ numHomeRangesScaled + region, data=pika, family=binomial)

# 		model5 <- glm(presAbs ~ meanDistToClosest4Patches, data=pika, family=binomial)
# 		model6 <- glm(presAbs ~ numHomeRangesScaled + meanDistToClosest4Patches, data=pika, family=binomial)
# 		model7 <- glm(presAbs ~ region + meanDistToClosest4Patches, data=pika, family=binomial)
# 		model8 <- glm(presAbs ~ numHomeRangesScaled + region + meanDistToClosest4Patches, data=pika, family=binomial)

# 		aicc1 <- AICc(model1)
# 		aicc2 <- AICc(model2)
# 		aicc3 <- AICc(model3)
# 		aicc4 <- AICc(model4)
# 		aicc5 <- AICc(model5)
# 		aicc6 <- AICc(model6)
# 		aicc7 <- AICc(model7)
# 		aicc8 <- AICc(model8)

# 		numHomeRanges <- c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)
# 		region <- c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)
# 		isolation <- c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
		
# 		aicc <- c(aicc1, aicc2, aicc3, aicc4, aicc5, aicc6, aicc7, aicc8)
		
# 		like1 <- logLik(model1)
# 		like2 <- logLik(model2)
# 		like3 <- logLik(model3)
# 		like4 <- logLik(model4)
# 		like5 <- logLik(model5)
# 		like6 <- logLik(model6)
# 		like7 <- logLik(model7)
# 		like8 <- logLik(model8)
		
# 		likeNull <- like1
		
# 		numHomRangesCoef2 <- coefficients(model2)['numHomeRangesScaled']
# 		numHomRangesCoef4 <- coefficients(model4)['numHomeRangesScaled']

# 		numHomRangesCoef6 <- coefficients(model6)['numHomeRangesScaled']
# 		numHomRangesCoef8 <- coefficients(model8)['numHomeRangesScaled']

# 		isolationCoef5 <- coefficients(model5)['meanDistToClosest4Patches']
# 		isolationCoef6 <- coefficients(model6)['meanDistToClosest4Patches']
# 		isolationCoef7 <- coefficients(model7)['meanDistToClosest4Patches']
# 		isolationCoef8 <- coefficients(model8)['meanDistToClosest4Patches']

# 		n <- nrow(pika)
# 		pseudoR2_1 <- nagelR2(likeNull, like1, n)
# 		pseudoR2_2 <- nagelR2(likeNull, like2, n)
# 		pseudoR2_3 <- nagelR2(likeNull, like3, n)
# 		pseudoR2_4 <- nagelR2(likeNull, like4, n)
# 		pseudoR2_5 <- nagelR2(likeNull, like5, n)
# 		pseudoR2_6 <- nagelR2(likeNull, like6, n)
# 		pseudoR2_7 <- nagelR2(likeNull, like7, n)
# 		pseudoR2_8 <- nagelR2(likeNull, like8, n)
		
# 		pseudoR2 <- c(pseudoR2_1, pseudoR2_2, pseudoR2_3, pseudoR2_4, pseudoR2_5, pseudoR2_6, pseudoR2_7, pseudoR2_8)

# 		coef1 <- coefficients(model1)
# 		coef2 <- coefficients(model2)
# 		coef3 <- coefficients(model3)
# 		coef4 <- coefficients(model4)
# 		coef5 <- coefficients(model5)
# 		coef6 <- coefficients(model6)
# 		coef7 <- coefficients(model7)
# 		coef8 <- coefficients(model8)

# 		nw3 <- coef3['regionnorthwest']
# 		se3 <- coef3['regionsoutheast']
# 		sw3 <- coef3['regionsouthwest']

# 		nw4 <- coef4['regionnorthwest']
# 		se4 <- coef4['regionsoutheast']
# 		sw4 <- coef4['regionsouthwest']

# 		nw7 <- coef7['regionnorthwest']
# 		se7 <- coef7['regionsoutheast']
# 		sw7 <- coef7['regionsouthwest']

# 		nw8 <- coef8['regionnorthwest']
# 		se8 <- coef8['regionsoutheast']
# 		sw8 <- coef8['regionsouthwest']

# 		results <- data.frame(
# 			model = '(Intercept)',
# 			term1 = NA,
# 			term2 = NA,
# 			term3 = NA,
# 			term4 = NA,
# 			numHomeRanges = numHomeRanges,
# 			homeRangeCoeff = c(NA, numHomRangesCoef2, NA, numHomRangesCoef4, NA, numHomRangesCoef6, NA, numHomRangesCoef8),
# 			isolationCoeff = c(NA, NA, NA, NA, isolationCoef5, isolationCoef6, isolationCoef7, isolationCoef8),
# 			region = region,
# 			aicc = aicc,
# 			pseudoR2 = pseudoR2,
# 			nw = c(NA, NA, nw3, nw4, NA, NA, nw7, nw8),
# 			se = c(NA, NA, se3, se4, NA, NA, se7, se8),
# 			sw = c(NA, NA, sw3, sw4, NA, NA, sw7, sw8)
# 		)

# 		# remember coefficients to do AICc-based coefficient weighting
# 		for (k in 5:8) {
		
# 			coefk <- get(paste0('coef', k))
# 			aicck <- get(paste0('aicc', k))

# 			for (i in seq_along(accumCoeffs)) {
# 				for (j in seq_along(coefk)) {
# 					if (names(coefk)[j] == names(accumCoeffs)[i]) {
# 						accumCoeffs[[i]] <- c(accumCoeffs[[i]], coefk[j])
# 						aiccs[[i]] <- c(aiccs[[i]], aicck)
# 					}
# 				}
			
# 			}
			
# 		}

# 	### by climate variable
# 	#######################

# 		for (formula in formulae) {
			
# 			say(formula)
			
# 			form1 <- as.formula(paste0('presAbs ~ 1 + ', formula))
# 			form2 <- as.formula(paste0('presAbs ~ 1 + ', formula, ' + numHomeRangesScaled'))
# 			form3 <- as.formula(paste0('presAbs ~ 1 + ', formula, ' + region'))
# 			form4 <- as.formula(paste0('presAbs ~ 1 + ', formula, ' + numHomeRangesScaled + region'))

# 			form5 <- as.formula(paste0('presAbs ~ 1 + ', formula, ' + meanDistToClosest4Patches'))
# 			form6 <- as.formula(paste0('presAbs ~ 1 + ', formula, ' + numHomeRangesScaled + meanDistToClosest4Patches'))
# 			form7 <- as.formula(paste0('presAbs ~ 1 + ', formula, ' + region + meanDistToClosest4Patches'))
# 			form8 <- as.formula(paste0('presAbs ~ 1 + ', formula, ' + numHomeRangesScaled + region + meanDistToClosest4Patches'))

# 			model1 <- glm(form1, data=pika, family=binomial)
# 			model2 <- glm(form2, data=pika, family=binomial)
# 			model3 <- glm(form3, data=pika, family=binomial)
# 			model4 <- glm(form4, data=pika, family=binomial)
# 			model5 <- glm(form5, data=pika, family=binomial)
# 			model6 <- glm(form6, data=pika, family=binomial)
# 			model7 <- glm(form7, data=pika, family=binomial)
# 			model8 <- glm(form8, data=pika, family=binomial)
			
# 			aicc1 <- AICc(model1)
# 			aicc2 <- AICc(model2)
# 			aicc3 <- AICc(model3)
# 			aicc4 <- AICc(model4)
# 			aicc5 <- AICc(model5)
# 			aicc6 <- AICc(model6)
# 			aicc7 <- AICc(model7)
# 			aicc8 <- AICc(model8)
			
# 			terms <- extractTerms(model1, model2, model3, model4, model5, model6, model7, model8)
# 			term1 <- terms$term1
# 			term2 <- terms$term2
# 			term3 <- terms$term3
# 			term4 <- terms$term4
			
# 			numHomeRanges <- c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)
# 			region <- c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)
# 			isolation <- c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
			
# 			aicc <- c(aicc1, aicc2, aicc3, aicc4, aicc5, aicc6, aicc7, aicc8)
			
# 			like1 <- logLik(model1)
# 			like2 <- logLik(model2)
# 			like3 <- logLik(model3)
# 			like4 <- logLik(model4)
# 			like5 <- logLik(model5)
# 			like6 <- logLik(model6)
# 			like7 <- logLik(model7)
# 			like8 <- logLik(model8)
			
# 			numHomRangesCoef2 <- coefficients(model2)['numHomeRangesScaled']
# 			numHomRangesCoef4 <- coefficients(model4)['numHomeRangesScaled']
# 			numHomRangesCoef6 <- coefficients(model6)['numHomeRangesScaled']
# 			numHomRangesCoef8 <- coefficients(model8)['numHomeRangesScaled']
			
# 			isolationCoef5 <- coefficients(model5)['meanDistToClosest4Patches']
# 			isolationCoef6 <- coefficients(model6)['meanDistToClosest4Patches']
# 			isolationCoef7 <- coefficients(model7)['meanDistToClosest4Patches']
# 			isolationCoef8 <- coefficients(model8)['meanDistToClosest4Patches']
			
# 			pseudoR2_1 <- nagelR2(likeNull, like1, n)
# 			pseudoR2_2 <- nagelR2(likeNull, like2, n)
# 			pseudoR2_3 <- nagelR2(likeNull, like3, n)
# 			pseudoR2_4 <- nagelR2(likeNull, like4, n)
# 			pseudoR2_5 <- nagelR2(likeNull, like5, n)
# 			pseudoR2_6 <- nagelR2(likeNull, like6, n)
# 			pseudoR2_7 <- nagelR2(likeNull, like7, n)
# 			pseudoR2_8 <- nagelR2(likeNull, like8, n)
			
# 			pseudoR2 <- c(pseudoR2_1, pseudoR2_2, pseudoR2_3, pseudoR2_4, pseudoR2_5, pseudoR2_6, pseudoR2_7, pseudoR2_8)

# 			coef1 <- coefficients(model1)
# 			coef2 <- coefficients(model2)
# 			coef3 <- coefficients(model3)
# 			coef4 <- coefficients(model4)
# 			coef5 <- coefficients(model5)
# 			coef6 <- coefficients(model6)
# 			coef7 <- coefficients(model7)
# 			coef8 <- coefficients(model8)

# 			nw3 <- coef3['regionnorthwest']
# 			se3 <- coef3['regionsoutheast']
# 			sw3 <- coef3['regionsouthwest']

# 			nw4 <- coef4['regionnorthwest']
# 			se4 <- coef4['regionsoutheast']
# 			sw4 <- coef4['regionsouthwest']

# 			nw7 <- coef7['regionnorthwest']
# 			se7 <- coef7['regionsoutheast']
# 			sw7 <- coef7['regionsouthwest']

# 			nw8 <- coef8['regionnorthwest']
# 			se8 <- coef8['regionsoutheast']
# 			sw8 <- coef8['regionsouthwest']

# 			results <- rbind(
# 				results,
# 				data.frame(
# 					model = formula,
# 					term1 = term1,
# 					term2 = term2,
# 					term3 = term3,
# 					term4 = term4,
# 					numHomeRanges = numHomeRanges,
# 					homeRangeCoeff = c(NA, numHomRangesCoef2, NA, numHomRangesCoef4, NA, numHomRangesCoef6, NA, numHomRangesCoef8),
# 					isolationCoeff = c(NA, NA, NA, NA, isolationCoef5, isolationCoef6, isolationCoef7, isolationCoef8),
# 					region = region,
# 					aicc = aicc,
# 					pseudoR2 = pseudoR2,
# 					nw = c(NA, NA, nw3, nw4, NA, NA, nw7, nw8),
# 					se = c(NA, NA, se3, se4, NA, NA, se7, se8),
# 					sw = c(NA, NA, sw3, sw4, NA, NA, sw7, sw8)
# 				)
# 			)
				
# 			# remember coefficients to do AICc-based coefficient weighting
# 			for (k in 1:8) {
			
# 				coefk <- get(paste0('coef', k))
# 				aicck <- get(paste0('aicc', k))

# 				for (i in seq_along(accumCoeffs)) {
# 					for (j in seq_along(coefk)) {
# 						if (names(coefk)[j] == names(accumCoeffs)[i]) {
# 							accumCoeffs[[i]] <- c(accumCoeffs[[i]], coefk[j])
# 							aiccs[[i]] <- c(aiccs[[i]], aicck)
# 						}
# 					}
				
# 				}
				
# 			}

# 		} # next formula

# 	### reports
# 	###########
	
# 		for (occWindow in c(occWindows_y, NA)) {
		
# 			if (is.na(occWindow)) {
# 				thisResults <- results
# 				nice <- paste0(paste(occWindows_y, collapse=' & '), '-yr Windows')
# 			} else {
# 				thisResults <- rbind(
# 					results[grepl(results$model, pattern=paste0(occWindow, 'yrWindow')), ],
# 					results[results$model == '(Intercept)', ]
# 				)
# 				nice <- paste0(occWindow, '-yr Window')
# 			}
			
# 			thisResults$deltaAicc <- thisResults$aicc - min(thisResults$aicc)
# 			w <- exp(-0.5 * thisResults$deltaAicc)
# 			thisResults$weight <- w / sum(w)

# 			thisResults <- thisResults[order(thisResults$weight, decreasing=TRUE), ]
# 			rownames(thisResults) <- NULL

# 			file <- paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Binary Models Using All Data - ', nice, '.csv')
# 			write.csv(thisResults, file, row.names=FALSE)
			
# 		} # next window

# 	### summarize coefficients
# 	### AICc-weighted

# 	minAicc <- min(results$aicc)
# 	deltaAicc <- results$aicc - minAicc
# 	w <- exp(-0.5 * deltaAicc)
# 	wSum <- sum(w)

# 	for (i in seq_along(accumCoeffs)) {
# 		aiccs[[i]] <- aiccs[[i]] - minAicc
# 		aiccs[[i]] <- exp(-0.5 * aiccs[[i]])
# 		aiccs[[i]] <- aiccs[[i]] / wSum
# 		accumCoeffs[[i]] <- accumCoeffs[[i]] * aiccs[[i]] / sum(aiccs[[i]])
# 	}
	
# 	accumCoeffs <- lapply(accumCoeffs, sum)

# 	sink('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Binary Models Using All Data - 7- & 10-yr Window Coefficient Summary.txt')
		
# 		say('AICc-weighted coefficient values for BINARY OCCUPANCY models', post=2)
# 		print(accumCoeffs)
		
# 	sink()


# say('##############################################################################################')
# say('### report relative odds of class change across regions across all BINARY OCCUPANCY models ###')
# say('##############################################################################################')

# 	file <- paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Binary Models Using All Data - 7 & 10-yr Windows.csv')
# 	results <- read.csv(file)
# 	results <- results[results$region, ]
# 	deltaAicc <- results$aicc - min(results$aicc)
# 	w <- exp(-0.5 * deltaAicc)
# 	w <- w / sum(w)
# 	nw <- sum(results$nw * w, na.rm=TRUE)
# 	se <- sum(results$se * w, na.rm=TRUE)
# 	sw <- sum(results$sw * w, na.rm=TRUE)
		
# 	sink('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Binary Models Using All Data Odds of Switching between Classes.txt')

# 		say('Binary odds ratios (arithmetic scale):')
			
# 		say('Binary odds of NW region relative to NE region: ', exp(nw))
# 		say('Binary odds of SE region relative to NE region: ', exp(se))
# 		say('Binary odds of SW region relative to NE region: ', exp(sw))
		
# 	sink()

# say('###############################################################################')
# say('### compile table of predictor weights for BINARY simple OCCUPANCY analysis ###')
# say('###############################################################################')

	# # rank variables by mean AICc weight
	# for (occWindow in c(occWindows_y, NA)) {
	
		# nice <- if (is.na(occWindow)) {
			# paste0(paste(occWindows_y, collapse=' & '), '-yr Windows')
		# } else {
			# paste0(occWindow, '-yr Window')
		# }
	
		# # get variables
		# vars <- getVars('occupancy')
		# if (!is.na(occWindow)) vars <- vars[grepl(vars, pattern=paste0(occWindow, 'yrWindow'))]
		# vars <- c(vars, 'meanDistToClosest4Patches')

		# # get models
		# file <- paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Binary Models Using All Data - ', nice, '.csv')
		# models <- read.csv(file)
		
		# imp <- data.frame()
		# for (var in vars) {
		
			# index <- if (var != 'meanDistToClosest4Patches') {
				# which(grepl(models$model, pattern=var))
			# } else {
				# which(!is.na(models$isolationCoef))
			# }
			# n <- length(index)
			# sumWeight <- sum(models$weight[index])
			# meanWeight <- sumWeight / n
			
			# # tally how many times this variable has + or - coefficient
			# coeffs <- models[index, c('term1', 'term2', 'term3', 'term4')]
			# if (var != 'meanDistToClosest4Patches') {

				# numPos <- numNeg <- 0
				# for (countIndex in seq_along(index)) {
					
					# thisIndex <- index[countIndex]
					# terms <- strsplit(models$model[thisIndex], split = ' \\+ ')[[1]]
					# whichOne <- which(grepl(terms, pattern = var))
					
					# # do not count if variable is part of an interaction term
					# if (!grepl(terms[whichOne], pattern = '\\*')) {
						
						# coeff <- coeffs[countIndex, whichOne]
						# coeff <- strsplit(coeff, split = ' ')[[1]]
						# coeff <- coeff[1]
						# coeff <- as.numeric(coeff)
						
						# if (coeff < 0) { numNeg <- numNeg + 1 } else { numPos <- numPos + 1 }
						
					# }
						
				# }

			# } else {
				# numNeg <- sum(models[index, 'isolationCoeff'] < 0)
				# numPos <- sum(models[index, 'isolationCoeff'] > 0)
			# }
			
			# imp <- rbind(
				# imp,
				# data.frame(
					# variable = var,
					# niceVar = makeNiceVars(var, 'occupancy'),
					# numModels = n,
					# sumWeight = sumWeight,
					# meanWeight = meanWeight,
					# numPos = numPos,
					# numNeg = numNeg
				# )
			# )
			
		# }
		
		# imp <- imp[order(imp$meanWeight, decreasing=TRUE), ]
		# imp$niceVar[imp$variable == 'meanDistToClosest4Patches'] <- 'isolation'

		# write.csv(imp, paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Binary Models Using All Data - ', nice, ' - Var Import.csv'), row.names=FALSE)
	
	# } # next occupancy window

# say('###################################################################')
# say('### double-checking number of models each variable should be in ###')
# say('###################################################################')
	
# 	vars <- getVars('occupancy')
# 	formulae <- getFormulaeOcc()
	
# 	counts <- data.frame()
# 	for (var in vars) {
	
# 		ins <- sum(grepl(var, formulae))
# 		counts <- rbind(
# 			counts,
# 			data.frame(
# 				var = var,
# 				n = ins
# 			)
# 		)
	
# 	}
	
# 	write.csv(counts, './Figures & Tables/Occupancy - Simple Models/Number of Base Models with Each Variable.csv', row.names=FALSE)

# say('##############################################################################################')
# say('### make maps of predicted probabilities of each class from best ordinal and binary models ###')
# say('##############################################################################################')

	# ### user-defined
	# ################
	
		# # Dates across which to calculate climate variables. The *months* will be used, so, for example, 2019-08-01 means to use August of 2019 (so is the same as 2019-08-31).
		# beginEndDates <- c('2009-09-01', '2019-08-01')
	
		# # assume all regions are this region when making predictions
		# region <- 'northwest'

		# # home folder in which PRISM rasters are stored... should end in "/an81"
		# prDir <- 'E:/Ecology/PRISM/working/an81' # HAL9000
		# # prDir <- 'F:/Ecology/PRISM/working/an81' # GRN
	
		# # formula of most-supported models obtained manually from analyses above
		# bestOrdinalForm <- latestOccStatus ~ occVar_chronicCold_C_10yrWindow + occVar_gsPpt_mm_10yrWindow + numHomeRangesScaled + region + meanDistToClosest4Patches

		# bestBinaryForm <- presAbs ~ occVar_chronicHeat_C_10yrWindow + occVar_monsoonPpt_mm_10yrWindow + numHomeRangesScaled + region + meanDistToClosest4Patches

	# ### pika data
	# #############
		
		# load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')

		# pika$latestOccStatus <- factor(pika$latestOccStatus, levels=c('0 never', '1 old', '2 occupied'), ordered=TRUE)
		# pika$region <- as.factor(pika$region)
		
		# pika$meanDistToClosest4Patches <- log10(pika$meanDistToClosest4Patches)
		# pika$numHomeRangesScaled <- as.numeric(scale(log10(pika$numHomeRanges + 1)))
		
		# vars <- getVars('occupancy')
		# vars <- c(vars, 'meanDistToClosest4Patches')
		# varsScaled <- scale(pika[ , vars])
		# varsCenter <- attr(varsScaled, 'scaled:center')
		# varsScale <- attr(varsScaled, 'scaled:scale')
		# pika[ , vars] <- varsScaled

	# ### models
	# ##########
	
		# bestOrdinalModel <- polr(bestOrdinalForm, data=pika, Hess=TRUE)
		# bestBinaryModel <- glm(bestBinaryForm, data=pika, family=binomial)

	# ### plot extent
	# ###############
		
		# pikaVectUnproj <- vect(as.matrix(pika[ , ll]), 'points', crs=getCRS('wgs84'))
		# pikaBuffUnprojSm <- terra::buffer(pikaVectUnproj, width = 40000)

	# ### create predictor rasters
	# ############################
	
		# devtools::load_all(paste0(drive, '/R/airUpThere'))
		
		# # chronic cold
		# x <- prStack(
			# prDir = prDir,
			# vars = 'tmean',
			# dates = beginEndDates,
			# span = TRUE,
			# by = 'month',
			# res = 800,
			# rastSuffix = 'tif'
		# )
		
		# names <- names(x)
		# months <- c(11, 12, 1, 2, 3)
		# months <- prefix(months, 2)
		# inFocalTimePeriod <- which(substr(names, 24, 25) %in% months)
		# x <- x[[inFocalTimePeriod]]
		
		# x <- crop(x, pikaBuffUnprojSm)
		# occVar_chronicCold_C_10yrWindow <- mean(x)
		# names(occVar_chronicCold_C_10yrWindow) <- 'occVar_chronicCold_C_10yrWindow'
		
		# # chronic heat
		# x <- prStack(
			# prDir = prDir,
			# vars = 'tmean',
			# dates = beginEndDates,
			# span = TRUE,
			# by = 'month',
			# res = 800,
			# rastSuffix = 'tif'
		# )
		
		# names <- names(x)
		# months <- 6:9
		# months <- prefix(months, 2)
		# inFocalTimePeriod <- which(substr(names, 24, 25) %in% months)
		# x <- x[[inFocalTimePeriod]]
		
		# x <- crop(x, pikaBuffUnprojSm)
		# occVar_chronicHeat_C_10yrWindow <- mean(x)
		# names(occVar_chronicHeat_C_10yrWindow) <- 'occVar_chronicHeat_C_10yrWindow'
		
		# # GS precipitation
		# x <- prStack(
			# prDir = prDir,
			# vars = 'ppt',
			# dates = beginEndDates,
			# span = TRUE,
			# by = 'month',
			# res = 800,
			# rastSuffix = 'tif'
		# )
		
		# names <- names(x)
		# months <- 5:9
		# months <- prefix(months, 2)
		# inFocalTimePeriod <- which(substr(names, 22, 23) %in% months)
		# x <- x[[inFocalTimePeriod]]

		# x <- crop(x, pikaBuffUnprojSm)
		# occVar_gsPpt_mm_10yrWindow <- sum(x) / 10
		# names(occVar_gsPpt_mm_10yrWindow) <- 'occVar_gsPpt_mm_10yrWindow'
				
		# # monsoon precipitation
		# x <- prStack(
			# prDir = prDir,
			# vars = 'ppt',
			# dates = beginEndDates,
			# span = TRUE,
			# by = 'month',
			# res = 800,
			# rastSuffix = 'tif'
		# )
		
		# names <- names(x)
		# months <- 6:8
		# months <- prefix(months, 2)
		# inFocalTimePeriod <- which(substr(names, 22, 23) %in% months)
		# x <- x[[inFocalTimePeriod]]

		# x <- crop(x, pikaBuffUnprojSm)
		# occVar_monsoonPpt_mm_10yrWindow <- sum(x) / 10
		# names(occVar_monsoonPpt_mm_10yrWindow) <- 'occVar_monsoonPpt_mm_10yrWindow'
		
		# env <- c(occVar_chronicCold_C_10yrWindow, occVar_chronicHeat_C_10yrWindow, occVar_gsPpt_mm_10yrWindow, occVar_monsoonPpt_mm_10yrWindow)
		# envDF <- as.data.frame(env, cells = TRUE)

		# # scale
		# vars <- names(envDF)
		# vars <- vars[vars != 'cell']
		# for (var in vars) {
			# mu <- varsCenter[[var]]
			# sigma <- varsScale[[var]]
			# envDF[ , var] <- (envDF[ , var] - mu) / sigma
		# }		

		# envDF$numHomeRangesScaled <- median(pika$numHomeRangesScaled)
		# envDF$meanDistToClosest4Patches <- median(pika$meanDistToClosest4Patches)
		# envDF$region <- region

	# ### make predictions
	# ####################

		# predOrdinal <- predict(bestOrdinalModel, envDF, type = 'p')
		# predBinary <- predict(bestBinaryModel, envDF, type = 'response')

	# ### assign predicted values to rasters
	# ######################################

		# template <- rast(paste0(drive, '/Research Data/PRISM/PRISM_us_dem_800m.tif'))
		# template <- crop(template, pikaBuffUnprojSm)

		# ordinalRast_noEvidence <- setValueByCell(template, val = predOrdinal[ , '0 never'], cell = envDF$cell)
		# ordinalRast_oldEvidence <- setValueByCell(template, val = predOrdinal[ , '1 old'], cell = envDF$cell)
		# ordinalRast_occupied <- setValueByCell(template, val = predOrdinal[ , '2 occupied'], cell = envDF$cell)

		# binaryRast <- template
		# binaryRast[] <- predBinary

		# ordinalRast <- c(ordinalRast_noEvidence, ordinalRast_oldEvidence, ordinalRast_occupied)
		# names(ordinalRast) <- c('no evidence', 'previously occupied', 'occupied')

		# binaryRast <- project(binaryRast, getCRS('NA Albers'))
		# ordinalRast <- project(ordinalRast, getCRS('NA Albers'))

	# ### hillshade
	# #############

		# elev <- rast('./Data/elev_fine_m.tif')
		# elev <- crop(elev, pikaBuffUnprojSm)
		
		# slope <- terrain(elev, 'slope', unit = 'radians')
		# aspect <- terrain(elev, 'aspect', unit = 'radians')

		# hs <- shade(slope, aspect, direction = 45)
		# hs <- project(hs, getCRS('NA Albers'))

	# ### plot
	# ########

		# pikaVect <- vect(pika, geom = c('longitude', 'latitude'), crs = getCRS('WGS84'))
		# pikaVect <- project(pikaVect, getCRS('NA Albers'))
		# pikaVect$presAbs <- as.character(pikaVect$presAbs)

		# pikaCols <- rep(NA_character_, nrow(pikaVect))
		# pikaCols[pikaVect$latestOccStatus == '0 never'] <- 'red'
		# pikaCols[pikaVect$latestOccStatus == '1 old'] <- 'orange'
		# pikaCols[pikaVect$latestOccStatus == '2 occupied'] <- 'green'

		# usa <- vect(paste0(drive, './Research Data/GADM/Version 4.1/High Res North America Level 2 sans Great Lakes SpatVector WGS84.gpkg'))
		# focusCounties <- crop(usa, ext(pikaBuffUnprojSm))
		# focusCounties <- project(focusCounties, getCRS('NA Albers'))

		# binaryMap <- ggplot() +
			# layer_spatial(hs, dpi = 600) +
			# layer_spatial(binaryRast, dpi = 600, alpha = 0.7) +
			# scale_fill_continuous(
				# name = 'Prediction',
				# na.value = NA,
				# low = 'red',
				# high = 'green3',
				# limits = c(0, 1)
			# ) +
			# layer_spatial(focusCounties, fill = NA) +
			# ggtitle('(a) Binary model predictions') +
			# theme_bw() +
			# theme(
				# panel.border = element_blank(),
				# plot.title = element_text(size = 10),
				# legend.title = element_text(size = 9),
				# legend.text = element_text(size = 8),
				# axis.text = element_text(size = 5),
				# axis.ticks = element_blank()
			# )		
		
		# pikaVect <- pikaVect[order(pikaVect$presAbs)]
		# presAbsSiteMap <- ggplot() +
			# layer_spatial(hs, dpi = 600) +
			# scale_fill_gradient(
				# low = 'gray0',
				# high = 'gray80',
				# guide = 'none',
				# na.value = NA
			# ) +
			# layer_spatial(focusCounties, fill = NA) +
			# layer_spatial(pikaVect, aes(color = presAbs, shape = presAbs), size = 0.8) +
			# scale_shape_manual(
				# name = 'Site\nstatus',
				# labels = c('Absent', 'Present'),
				# values = c('0' = 4, '1' = 1)
			# ) +
			# scale_color_manual(
				# name = 'Site\nstatus',
				# labels = c('Absent', 'Present'),
				# values = c('0' = 'red', '1' = 'green3')
			# ) +
			# ggtitle('(b) Binary site status') +
			# theme_bw() +
			# theme(
				# panel.border = element_blank(),
				# plot.title = element_text(size = 10),
				# axis.text = element_text(size = 5),
				# axis.ticks = element_blank(),
				# legend.title = element_text(size = 8),
				# legend.text = element_text(size = 7)
			# )		

		# ordinalMaps <- ordinalSiteMaps <- list()
		# low <- NA
		# for (i in 1:3) {
			
			# if (i == 1) {
				# x <- ordinalRast_noEvidence
				# title <- 'No evidence'
				# values = c('0 never' = 'red', '1 old' = 'gray30', '2 occupied' = 'gray30')
				# high <- 'red'
				# pikaVect <- pikaVect[order(pikaVect$latestOccStatus, decreasing = TRUE)]
			# } else if (i == 2) {
				# x <- ordinalRast_oldEvidence
				# title <- 'Previously occupied'
				# values = c('0 never' = 'gray30', '1 old' = 'orange', '2 occupied' = 'gray30')
				# high <- 'orange'
				# pikaVect <- rbind(pikaVect[pikaVect$latestOccStatus != '1 old'], pikaVect[pikaVect$latestOccStatus == '1 old'])
			# } else {
				# x <- ordinalRast_occupied
				# title <- 'Occupied'
				# values = c('0 never' = 'gray30', '1 old' = 'gray30', '2 occupied' = 'green3')
				# high <- 'green3'
				# pikaVect <- pikaVect[order(pikaVect$latestOccStatus)]
			# }
			
			# predictionTitle <- paste0('(', letters[2 * i + 1], ') Ordinal model predictions:\n      ', title)
			# siteTitle <- paste0('(', letters[2 * i + 2], ') Ordinal site status:\n      ', title, ' highlighted')
				
			# ordinalMaps[[i]] <- ggplot() +
				# layer_spatial(hs, dpi = 600) +
				# layer_spatial(x, dpi = 600, alpha = 0.7) +
				# scale_fill_continuous(
					# name = 'Prediction',
					# na.value = NA,
					# low = low,
					# high = high,
					# limits = c(0, 1)
				# ) +
				# layer_spatial(focusCounties, fill = NA) +
				# # layer_spatial(pikaVect, pch = 1, aes(color = latestOccStatus)) +
				# # scale_color_manual(
					# # name = 'Site\nstatus',
					# # labels = c('No\n  evidence', 'Previously\n  occupied', 'Occupied'),
					# # values = c('0 never' = 'red', '1 old' = 'orange', '2 occupied' = 'green3')
				# # ) +
				# ggtitle(predictionTitle)
				
			# ordinalSiteMaps[[i]] <- ggplot() +
				# layer_spatial(hs, dpi = 600) +
				# scale_fill_gradient(
					# low = 'gray0',
					# high = 'gray80',
					# guide = 'none',
					# na.value = NA
				# ) +
				# layer_spatial(focusCounties, fill = NA) +
				# layer_spatial(pikaVect, aes(color = latestOccStatus, shape = latestOccStatus), size = 0.8) +
				# scale_shape_manual(
					# name = 'Site\nstatus',
					# labels = c('No\n  evidence', 'Previously\n    occupied', 'Occupied'),
					# values = c('0 never' = 4, '1 old' = 2, '2 occupied' = 1)
				# ) +
				# scale_color_manual(
					# name = 'Site\nstatus',
					# labels = c('No\n  evidence', 'Previously\n    occupied', 'Occupied'),
					# values = values
				# ) +
				# ggtitle(siteTitle) +
				# theme_bw() +
				# theme(
					# panel.border = element_blank(),
					# plot.title = element_text(size = 10),
					# axis.text = element_text(size = 5),
					# axis.ticks = element_blank()
				# )

				# ordinalMaps[[i]] <- ordinalMaps[[i]] +
					# theme_bw() +
					# theme(
						# panel.border = element_blank(),
						# plot.title = element_text(size = 10),
						# legend.title = element_text(size = 9),
						# legend.text = element_text(size = 8),
						# axis.text = element_text(size = 5),
						# axis.ticks = element_blank()
					# )			

				# ordinalSiteMaps[[i]] <- ordinalSiteMaps[[i]] +
					# theme_bw() +
					# theme(
						# panel.border = element_blank(),
						# plot.title = element_text(size = 10),
						# axis.text = element_text(size = 5),
						# axis.ticks = element_blank(),
						# legend.title = element_text(size = 8),
						# legend.text = element_text(size = 7)
					# )			

		# }
		
		# mapsArranged <- list(
			# binaryMap, presAbsSiteMap,
			# ordinalMaps[[1]], ordinalSiteMaps[[1]],
			# ordinalMaps[[2]], ordinalSiteMaps[[2]],
			# ordinalMaps[[3]], ordinalSiteMaps[[3]]
		# )
		
		# maps <- plot_grid(plotlist = mapsArranged, ncol = 2)

		# region <- capIt(region)
		# ggsave(maps, file = paste0('./Figures & Tables/Predictions Assuming ', region, '.png'), dpi = 600, width = 7, height = 11, bg = 'white')

say('#######################################################################################################')
say('### make maps of future predicted probabilities from binary models and extract predictions at sites ###')
say('#######################################################################################################')

	### user-defined
	################
	
		# are we using the 7- or 10-yr window models?
		# window_y <- 7
		window_y <- 10
	
		# use all models with deltaAICc < this amount
		# NB if we use 4 for the 10-yr window, we need to somehow calculate number of days with heat >18 C
		deltaAICcThreshold <- 2
		
		# SSPs for future scenarios
		ssps <- c(245, 370)
		
		# year ranges for future scenarios
		periods <- c('2041_2070', '2071_2100')
		
		# where are ClimateNA rasters stored?
		climateNAPath <- paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest')
		
		# NB Future climate rasters are from ClimateNA. To define variables in the top models, I had to manually look at the top models and and select the variables. This chunk won't necessarily run if the predictors in the top models change.

	### models
	##########
	
		models <- read.csv(paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Binary Models Using All Data - ', window_y, '-yr Window.csv'))
		models <- models[order(models$deltaAicc), ]
		nBestModels <- sum(models$deltaAicc < deltaAICcThreshold)

	### pika data
	#############

		load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')
		pika$region <- as.factor(pika$region)
		
		scaled <- scale(log10(pika$meanDistToClosest4Patches))
		isolationCenter <- attr(scaled, 'scaled:center')
		isolationScale <- attr(scaled, 'scaled:scale')
		pika$meanDistToClosest4Patches <- as.numeric(scaled)
		
		scaled <- scale(log10(pika$numHomeRanges + 1))
		numHomeRangesCenter <- attr(scaled, 'scaled:center')
		numHomeRangesScale <- attr(scaled, 'scaled:scale')
		pika$numHomeRangesScaled <- as.numeric(scaled)
		
		vars <- getVars('occupancy')
		vars <- c(vars, 'meanDistToClosest4Patches')
		scaled <- scale(pika[ , vars])
		centers <- attr(scaled, 'scaled:center')
		scales <- attr(scaled, 'scaled:scale')
		pika[ , vars] <- scaled

	### elevation
	#############
	
		elevation <- rast(paste0(climateNAPath, '/elevation.tif'))

	### study (plot) region
	#######################
	
		pikaSpatial <- vect(pika, geom = c('longitude', 'latitude'), crs = getCRS('WGS84'))
		pikaSpatial <- project(pikaSpatial, elevation)
		studyRegionFocus <- buffer(pikaSpatial, 40000)
		studyRegionExtended <- buffer(pikaSpatial, 50000)

		studyRegionFocus <- ext(studyRegionFocus)
		studyRegionExtended <- ext(studyRegionExtended)
		
		studyRegionFocus <- as.polygons(studyRegionFocus, crs = elevation)
		studyRegionExtended <- as.polygons(studyRegionExtended, crs = elevation)

	### calculate best models
	#########################

		trainedModels <- list()
		for (i in 1:nBestModels) {

			# calibrate model
			form <- paste0('presAbs ~ ', models$model[i])
			if (models$region[i]) form <- paste0(form, ' + region')
			if (models$numHomeRanges[i]) form <- paste0(form, ' + numHomeRangesScaled')
			if (!is.na(models$isolationCoef[i])) form <- paste0(form, ' + meanDistToClosest4Patches')
			
			form <- as.formula(form)
			trainedModels[[i]] <- glm(form, data = pika, family = binomial)
		
		}
		
	### predict to future rasters
	#############################

	# holds rasters of future predictions... will have 4 elements, each of which is a raster stack of 4 rasters, one per region
	futPredictions <- list()
	for (ssp in ssps) {
		for (period in periods) {	

			# predictor rasters
			futPredictors <- list()

			### compile futPredictors from ClimateNA
			#####################################
			
				### chronicHeat_C: Mean of summer average temperature (tmean), June-September
				#############################################################################
				
				months <- 6:9
				name <- paste0('occVar_chronicHeat_C_', window_y, 'yrWindow')
				cnaVar <- 'Tave'
				
				x <- rast(paste0(climateNAPath, '/ensemble_8GCMs_ssp', ssp, '_', period, '_monthly/ensemble_8GCMs_ssp', ssp, '_', period, '_', cnaVar, prefix(months, 2), '.tif'))
				x <- crop(x, studyRegionExtended)
				x <- mean(x)
				x <- scale(x, center = centers[[name]], scale = scales[[name]])
				names(x) <- name
				futPredictors[[length(futPredictors) + 1]] <- x
								
				### monsoonPpt_mm: Total monsoon precipitation (ppt), mid-June-August
				#####################################################################
				
				months1 <- 6
				months2 <- 7:8
				name <- paste0('occVar_monsoonPpt_mm_', window_y, 'yrWindow')
				cnaVar <- 'PPT'

				# use half of June ppt
				x1 <- rast(paste0(climateNAPath, '/ensemble_8GCMs_ssp', ssp, '_', period, '_monthly/ensemble_8GCMs_ssp', ssp, '_', period, '_', cnaVar, prefix(months1, 2), '.tif'))
				x2 <- rast(paste0(climateNAPath, '/ensemble_8GCMs_ssp', ssp, '_', period, '_monthly/ensemble_8GCMs_ssp', ssp, '_', period, '_', cnaVar, prefix(months2, 2), '.tif'))

				x1 <- crop(x2, studyRegionExtended)
				x2 <- crop(x2, studyRegionExtended)

				x1 <- 0.5 * x1
				x <- c(x1, x2)

				x <- sum(x)
				x <- scale(x, center = centers[[name]], scale = scales[[name]])
				names(x) <- name
				futPredictors[[length(futPredictors) + 1]] <- x
								
				### gsPpt_mm: Total growing-season precipitation (ppt), May-September
				#####################################################################

				months <- 5:9
				name <- paste0('occVar_gsPpt_mm_', window_y, 'yrWindow')
				cnaVar <- 'PPT'
				
				x <- rast(paste0(climateNAPath, '/ensemble_8GCMs_ssp', ssp, '_', period, '_monthly/ensemble_8GCMs_ssp', ssp, '_', period, '_', cnaVar, prefix(months, 2), '.tif'))
				x <- crop(x, studyRegionExtended)
				x <- sum(x)
				x <- scale(x, center = centers[[name]], scale = scales[[name]])
				names(x) <- name
				futPredictors[[length(futPredictors) + 1]] <- x
				
			### non-climatic future
			#######################
				
				# isolation
				x <- futPredictors[[1]]
				x[] <- isolationCenter
				names(x) <- 'meanDistToClosest4Patches'
				futPredictors[[length(futPredictors) + 1]] <- x
				
				# number of home ranges
				x <- futPredictors[[1]]
				x[] <- numHomeRangesCenter
				names(x) <- 'numHomeRangesScaled'
				futPredictors[[length(futPredictors) + 1]] <- x
				
				futPredictors <- do.call(c, futPredictors)
		
			### make predictions... cycle through each region
			#################################################
			
				regions <- c('northwest', 'southwest', 'northeast', 'southeast')
				regionPredictions <- list()
				for (region in regions) {

					x <- futPredictors[[1]]
					x[] <- 1
					names(x) <- paste0('region', region)
					x <- as.factor(x)
					levs <- data.frame(value = 1, region = region)
					levels(x) <- levs
					thesePredictors <- c(futPredictors, x)

					thesePredictions <- list()
					weights <- numeric()
					for (i in 1:nBestModels) {
					
						weight <- models$weight[i]
						weights[i] <- weight

						thesePredictions[[i]] <- predict(thesePredictors, trainedModels[[i]], type = 'response')
						thesePredictions[[i]] <- thesePredictions[[i]] * weight
					
					}
					
					thesePredictions <- do.call(c, thesePredictions)
					thesePredictions <- sum(thesePredictions) / sum(weights)
				
					regionPredictions[[length(regionPredictions) + 1]] <- thesePredictions
					names(regionPredictions)[[length(regionPredictions)]] <- region
				
				}

				regionPredictions <- c(regionPredictions[[1]], regionPredictions[[2]], regionPredictions[[3]], regionPredictions[[4]])
				names(regionPredictions) <- regions
		
				futPredictions[[length(futPredictions) + 1]] <- regionPredictions
				names(futPredictions)[length(futPredictions)] <- paste(ssp, period)
		
		} # next period
		
	} # next SSP

	### predict to present-day conditions using PRISM
	#################################################
	
		### climate predictors
		######################
	
			devtools::load_all(paste0(drive, '/R/airUpThere'))
			
			studyRegionExtendedNAD83 <- project(studyRegionExtended, getCRS('NAD83'))
			
			# home folder in which PRISM rasters are stored... should end in "/an81"
			prDir <- 'E:/Ecology/PRISM/working/an81' # HAL9000
			# prDir <- 'F:/Ecology/PRISM/working/an81' # GRN

			beginEndDates <- c('2009-09-01', '2019-08-01')
			
			predictors <- list()
			
			# chronic heat
			name <- paste0('occVar_chronicHeat_C_', window_y, 'yrWindow')
			x <- prStack(
				prDir = prDir,
				vars = 'tmean',
				dates = beginEndDates,
				span = TRUE,
				by = 'month',
				res = 800,
				rastSuffix = 'tif'
			)
			
			names <- names(x)
			months <- 6:9
			months <- prefix(months, 2)
			inFocalTimePeriod <- which(substr(names, 24, 25) %in% months)
			x <- x[[inFocalTimePeriod]]
			
			x <- crop(x, studyRegionExtendedNAD83)
			x <- mean(x)
			x <- scale(x, center = centers[[name]], scale = scales[[name]])
			names(x) <- name
			predictors[[length(predictors) + 1]] <- x
			
			# GS precipitation
			name <- paste0('occVar_gsPpt_mm_', window_y, 'yrWindow')
			x <- prStack(
				prDir = prDir,
				vars = 'ppt',
				dates = beginEndDates,
				span = TRUE,
				by = 'month',
				res = 800,
				rastSuffix = 'tif'
			)
			
			names <- names(x)
			months <- 5:9
			months <- prefix(months, 2)
			inFocalTimePeriod <- which(substr(names, 22, 23) %in% months)
			x <- x[[inFocalTimePeriod]]

			x <- crop(x, studyRegionExtendedNAD83)
			x <- sum(x) / 10
			x <- scale(x, center = centers[[name]], scale = scales[[name]])
			names(x) <- name
			predictors[[length(predictors) + 1]] <- x
					
			# monsoon precipitation
			name <- paste0('occVar_monsoonPpt_mm_', window_y, 'yrWindow')
			x <- prStack(
				prDir = prDir,
				vars = 'ppt',
				dates = beginEndDates,
				span = TRUE,
				by = 'month',
				res = 800,
				rastSuffix = 'tif'
			)
			
			names <- names(x)
			months <- 6:8
			months <- prefix(months, 2)
			inFocalTimePeriod <- which(substr(names, 22, 23) %in% months)
			x <- x[[inFocalTimePeriod]]

			x <- crop(x, studyRegionExtendedNAD83)
			x <- sum(x) / 10
			x <- scale(x, center = centers[[name]], scale = scales[[name]])
			names(x) <- name
			predictors[[length(predictors) + 1]] <- x
		
		### non-climatic predictors
		###########################
			
			# isolation
			x <- predictors[[1]]
			x[] <- isolationCenter
			names(x) <- 'meanDistToClosest4Patches'
			predictors[[length(predictors) + 1]] <- x
			
			# number of home ranges
			x <- predictors[[1]]
			x[] <- numHomeRangesCenter
			names(x) <- 'numHomeRangesScaled'
			predictors[[length(predictors) + 1]] <- x
			
		### make predictions... cycle through each region
		#################################################

			sqPredictors <- do.call(c, predictors)
		
			regions <- c('northwest', 'southwest', 'northeast', 'southeast')
			regionPredictions <- list()
			for (region in regions) {

				x <- sqPredictors[[1]]
				x[] <- 1
				names(x) <- paste0('region', region)
				x <- as.factor(x)
				levs <- data.frame(value = 1, region = region)
				levels(x) <- levs
				thesePredictors <- c(sqPredictors, x)

				thesePredictions <- list()
				weights <- numeric()
				for (i in 1:nBestModels) {
				
					weight <- models$weight[i]
					weights[i] <- weight

					thesePredictions[[i]] <- predict(thesePredictors, trainedModels[[i]], type = 'response')
					thesePredictions[[i]] <- thesePredictions[[i]] * weight
				
				}
				
				thesePredictions <- do.call(c, thesePredictions)
				thesePredictions <- sum(thesePredictions) / sum(weights)
			
				regionPredictions[[length(regionPredictions) + 1]] <- thesePredictions
				names(regionPredictions)[[length(regionPredictions)]] <- region
			
			}

			sqPredictions <- c(regionPredictions[[1]], regionPredictions[[2]], regionPredictions[[3]], regionPredictions[[4]])
			names(sqPredictions) <- regions

			sqPredictions <- project(sqPredictions, futPredictions[[1]])

		### crop rasters to respective regions
		######################################
		
			### make polygon with 4 boxes, one per region
			extent <- ext(studyRegionExtended)
			extent <- as.vector(extent)
			
			eastWest <- -534756.3
			northSouth <- -947687.6
			
			sw <- c(extent[1], eastWest, extent[3], northSouth)
			nw <- c(extent[1], eastWest, northSouth, extent[4])
			ne <- c(eastWest, extent[2], northSouth, extent[4])
			se <- c(eastWest, extent[2], extent[3], northSouth)

			sw <- ext(sw)
			nw <- ext(nw)
			ne <- ext(ne)
			se <- ext(se)
			
			sw <- as.polygons(sw, crs = elevation)
			nw <- as.polygons(nw, crs = elevation)
			ne <- as.polygons(ne, crs = elevation)
			se <- as.polygons(se, crs = elevation)

			boundaries <- list(northwest = nw, northeast = ne, southwest = sw, southeast = se)

			# crop rasters to respective region
			for (region in regions) {
				sqPredictions[[region]] <- mask(sqPredictions[[region]], boundaries[[region]])
			}
			
			for (ssp in ssps) {
				for (period in periods) {
					for (region in regions) {
						futPredictions[[paste(ssp, period)]][[region]] <- mask(futPredictions[[paste(ssp, period)]][[region]], boundaries[[region]])
					}
				}
			}

		### change in predictions
		#########################
		
			deltaPredictions <- list()
			i <- 1
			for (ssp in ssps) {
				for (period in periods) {
					deltaPredictions[[i]] <- futPredictions[[i]] - sqPredictions
					i <- i + 1
				}
			}
			names(deltaPredictions) <- names(futPredictions)
		
		### extract predictions to sampled sites in present and future
		##############################################################

			pikaPoints <- vect(pika, geom = ll, crs = getCRS('NAD83'))
			presPoints <- pikaPoints[pikaPoints$latestOccStatus == '2 occupied']
			presPoints <- project(presPoints, elevation)

			# extract for present
			presPoints$predictionsPresent <- NA_real_
			for (region in regions) {
				
				index <- which(presPoints$region == region)
				presThisRegion <- presPoints[index]
				preds <- extract(sqPredictions, presThisRegion)[ , region, drop = TRUE]
				presPoints$predictionsPresent[index] <- preds

			}


			# extract for future
			for (ssp in ssps) {
				for (period in periods) {

					fut <- paste(ssp, gsub(period, pattern = '-', replacement = '_'))
					presPoints$DUMMY <- NA_real_

					for (region in regions) {

						index <- which(presPoints$region == region)
						presThisRegion <- presPoints[index]
						preds <- extract(futPredictions[[fut]], presThisRegion)[ , region, drop = TRUE]
						presPoints$DUMMY[index] <- preds

					}
					names(presPoints)[ncol(presPoints)] <- paste0('ssp', ssp, '_', gsub(period, pattern = '-', replacement = '_'))
					i <- i + 1
				}
			}

		### summaries of change in probability of presence through time
		###############################################################

		sink('./Figures & Tables/Change in Suitability through Time from Binary Models.txt', split = TRUE)
		say('SUMMARY OF PROBABILITIES OF PRESENCE AT PRESENCE SITES THROUGH TIME FROM BINARY MODELS')
		say(date())

		n <- nrow(presPoints)
		say('Total number of presence sites: ', n, pre = 1)

		say('Predictions to present:', pre = 1)
		print(summary(presPoints$predictionsPresent))

		say('Predictions to SSP 245, 2041-2070:', pre = 1)
		print(summary(presPoints$ssp245_2041_2070))

		say('Predictions to SSP 245, 2071-2100:', pre = 1)
		print(summary(presPoints$ssp245_2071_2100))

		say('Predictions to SSP 370, 2041-2070:', pre = 1)
		print(summary(presPoints$ssp370_2041_2070))

		say('Predictions to SSP 370, 2071-2100:', pre = 1)
		print(summary(presPoints$ssp370_2071_2100))

		say('CHANGE in predictions from present to SSP 245, 2041-2070 (future - present) & proportion of sites with losses:', pre = 1)
		delta <- presPoints$ssp245_2041_2070 - presPoints$predictionsPresent
		print(summary(delta))
		print(sum(delta < 0) / n)

		say('CHANGE in predictions from present to SSP 245, 2071-2100 (future - present) & proportion of sites with losses:', pre = 1)
		delta <- presPoints$ssp245_2071_2100 - presPoints$predictionsPresent
		print(summary(delta))
		print(sum(delta < 0) / n)

		say('CHANGE in predictions from present to SSP 370, 2041-2070 (future - present) & proportion of sites with losses:', pre = 1)
		delta <- presPoints$ssp370_2041_2070 - presPoints$predictionsPresent
		print(summary(delta))
		print(sum(delta < 0) / n)

		say('CHANGE in predictions from present to SSP 370, 2071-2100 (future - present) & proportion of sites with losses:', pre = 1)
		delta <- presPoints$ssp370_2071_2100 - presPoints$predictionsPresent
		print(summary(delta))
		print(sum(delta < 0) / n)

		for (region in regions) {
		
			say(region, level = 2)

			index <- which(presPoints$region == region)
			n <- length(index)

			say('CHANGE in predictions from present to SSP 245, 2041-2070 (future - present) & proportion of sites with losses:', pre = 1)
			delta <- presPoints$ssp245_2041_2070[index] - presPoints$predictionsPresent[index]
			print(summary(delta))
			print(sum(delta < 0) / n)

			say('CHANGE in predictions from present to SSP 245, 2071-2100 (future - present) & proportion of sites with losses:', pre = 1)
			delta <- presPoints$ssp245_2071_2100[index] - presPoints$predictionsPresent[index]
			print(summary(delta))
			print(sum(delta < 0) / n)

			say('CHANGE in predictions from present to SSP 370, 2041-2070 (future - present) & proportion of sites with losses:', pre = 1)
			delta <- presPoints$ssp370_2041_2070[index] - presPoints$predictionsPresent[index]
			print(summary(delta))
			print(sum(delta < 0) / n)

			say('CHANGE in predictions from present to SSP 370, 2071-2100 (future - present) & proportion of sites with losses:', pre = 1)
			delta <- presPoints$ssp370_2071_2100[index] - presPoints$predictionsPresent[index]
			print(summary(delta))
			print(sum(delta < 0) / n)
		
		}

		sink()

		# # ### graph and summaries of changes in suitability across pika-occupied sites
		# # ############################################################################

		# # 	cols <- c('predictionsPresent', 'ssp245_2041_2070', 'ssp245_2071_2100', 'ssp370_2041_2070', 'ssp370_2071_2100')

		# # 	sqStartEndYears <- beginEndDates
		# # 	sqStartEndYears <- substr(sqStartEndYears, 1, 4)
		# # 	sqStartEndYears <- paste(sqStartEndYears, collapse = '-')

		# # 	presPointsLong <- data.frame(
		# # 		site = presPoints$polygonName,
		# # 		region = presPoints$region,
		# # 		period = sqStartEndYears,
		# # 		prediction = presPointsWide$predictionsPresent
		# # 	)

		# # 	for (ssp in ssps) {
		# # 		for (period in periods) {

		# # 			presPointsLong <- rbind(
		# # 				presPointsLong,
		# # 				data.frame(
		# # 					site = presPoints$polygonName,
		# # 					region = presPoints$region,
		# # 					period = paste0('SSP', ssp, ' ', gsub(period, pattern = '_', replacement = '-')),
		# # 					prediction = presPointsWide[ , paste0('ssp', ssp, '_', period)]
		# # 				)
		# # 			)

		# # 		}
		# # 	}

		# # 	# plot boxplot for each region/time period and join values for same site across boxes
		# # 	predByRegion <- list()
		# # 	for (region in regions) {

		# # 		data <- presPointsLong[presPointsLong$region == region, ]

		# # 		col <- if (presPoints$region[i] == 'northwest') {
		# # 			'#1b9e77'
		# # 		} else if (presPoints$region[i] == 'southwest') {
		# # 			'#d95f02'
		# # 		} else if (presPoints$region[i] == 'northeast') {
		# # 			'#7570b3'
		# # 		} else if (presPoints$region[i] == 'southeast') {
		# # 			'#e7298a'
		# # 		}

		# # 		main <- capIt(region)
		# # 		main <- as.character(main)

		# # 		predByRegion[[length(predByRegion) + 1]] <-
		# # 			ggplot(data, aes(x = period, y = prediction)) +
		# # 			geom_boxplot(fill = col) +
		# # 			geom_point() +
		# # 			geom_line(aes(group = site)) +
		# # 			ylab('Probability of occurrence') + ylim(0, 1) +
		# # 			ggtitle(main)

		# # 	}


		### hillshade map
		#################

			elev <- rast('./Data/elev_fine_m.tif')
			elev <- crop(elev, studyRegionExtendedNAD83)
			
			slope <- terrain(elev, 'slope', unit = 'radians')
			aspect <- terrain(elev, 'aspect', unit = 'radians')

			hs <- shade(slope, aspect, direction = 45)
			hs <- project(hs, elevation)
			hs <- trim(hs)
			
		### contours
		############
		
			contours <- as.contour(elev, levels = 2308)
			contours <- project(contours, elevation)
			
		### counties
		############
		
			# New Mexico		
			usa1 <- vect(paste0(drive, './Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector WGS84.gpkg'))
			nm <- usa1[usa1$NAME_1 == 'New Mexico']
			nm <- project(nm, elevation)

			# counties
			usa2 <- vect(paste0(drive, './Research Data/GADM/Version 4.1/High Res North America Level 2 sans Great Lakes SpatVector WGS84.gpkg'))

			focusCounties <- crop(usa2, ext(studyRegionExtendedNAD83))
			focusCounties <- project(focusCounties, elevation)

		### map of present suitability
		##############################
		
			dirCreate('./Figures & Tables/Maps of Predictions')

			maxcell <- 1E5
			# maxcell <- Inf

			period <- substr(beginEndDates, 1, 4)
			period <- paste0(period[1], '-', period[2])
			
			extent <- as.vector(ext(studyRegionFocus))
			xlim <- extent[1:2]
			ylim <- extent[3:4]

			png(paste0('./Figures & Tables/Maps of Predictions/Binary Occurrence ', period, '.png'), res = 600, width = 1200, height = 1300)

				main <- paste0('(a) Suitability ', period)
				# col <- paste0('gray', 0:50)
				col <- paste0('gray', 0:100)
				plot(hs, col = col, legend = FALSE, axes = FALSE, ext = studyRegionFocus, box = FALSE, maxcell = maxcell, mar = c(0.4, 0.1, 0.1, 0.1))
				
				col <- colorRampPalette(c('gray90', 'green2'))(10)
				
				plot(sqPredictions[['northwest']], col = col, alpha = 0.5, legend = FALSE, axes = FALSE, ext = studyRegionFocus, box = FALSE, maxcell = maxcell, range = c(0, 1), add = TRUE)
				plot(sqPredictions[['southwest']], col = col, alpha = 0.5, legend = FALSE, axes = FALSE, ext = studyRegionFocus, box = FALSE, maxcell = maxcell, range = c(0, 1), add = TRUE)
				plot(sqPredictions[['southeast']], col = col, alpha = 0.5, legend = FALSE, axes = FALSE, ext = studyRegionFocus, box = FALSE, maxcell = maxcell, range = c(0, 1), add = TRUE)
				plot(sqPredictions[['northeast']], col = col, alpha = 0.5, legend = FALSE, axes = FALSE, ext = studyRegionFocus, box = FALSE, maxcell = maxcell, range = c(0, 1), add = TRUE)
				
				plot(focusCounties, lwd = 0.3, add = TRUE)
				plot(studyRegionFocus, add = TRUE, lwd = 0.7)
				plot(nm, lwd = 0.7, add = TRUE)
				
				plot(contours, lwd = 0.4, add = TRUE)

				legendGrad(
					x = 'bottom',
					inset = -0.02,
					vert = FALSE,
					width = 1,
					height = 0.06,
					adjX = c(0.25, 1),
					col = c('gray90', 'green2'),
					title = 'Occurrence\nProbability',
					titleCex = 0.45,
					titleAdj = c(0.12, 0.4),
					labels = '',
					boxBorder = NA,
					lwd = 0.8,
					xpd = NA
				)
				
				text(-592000, -1065937, labels = 0, cex = 0.4, xpd = NA)
				text(-515000, -1065937, labels = 0.5, cex = 0.4, xpd = NA)
				text(-430000, -1065937, labels = 1, cex = 0.4, xpd = NA)
				labelFig(label = main, adj = c(-0.05, -0.02), cex = 0.4)
				
			dev.off()
		
		### map of change in future suitability
		#######################################

			# color ranges
			ranges <- c(Inf, -Inf)
			for (ssp in ssps) {
				for (period in periods) {
					mm <- minmax(deltaPredictions[[paste(ssp, period)]])
					ranges[1] <- min(ranges[1], mm[1, ])
					ranges[2] <- max(ranges[2], mm[2, ])
				}
			}

			colNeg <- colorRampPalette(c('red', 'gray90'))
			# colPos <- colorRampPalette(c('gray90', 'green3'))
			# colNeg <- colorRampPalette(c('#d73027', 'gray90'))
			colPos <- colorRampPalette(c('gray90', '#4575b4'))
			negRanges <- abs(100 * round(ranges[1], 2))
			posRanges <- 100 * round(ranges[2], 2)
			
			maxRanges <- max(negRanges, posRanges)
			colNeg <- colNeg(maxRanges)
			colPos <- colPos(maxRanges)
			
			colNeg <- colNeg[1:negRanges]
			colPos <- colPos[1:posRanges]

			deltaCol <- c(colNeg, colPos)

			png(paste0('./Figures & Tables/Maps of Predictions/Binary Occurrence Futures.png'), res = 600, width = 2 * 1200, height = 2 * 1300)
			
				par(mfrow = c(2, 2), mar = c(0, 0, 0, 0))

				i <- 1
				for (ssp in ssps) {
				
					for (period in periods) {

						nicePeriod <- if (period == '2041_2070') { '2050s' } else { '2080s' }
						
						extent <- as.vector(ext(studyRegionFocus))
						xlim <- extent[1:2]
						ylim <- extent[3:4]

						if (i %in% 1:2) { mar <- c(0, 0.1, 0.1, 0.1) } else { mar <- c(0.1, 0.1, 0, 0.1) }

						main <- paste0('(', letters[i + 1], ') SSP ', ssp, ': ', nicePeriod)
						# col <- paste0('gray', 0:50)
						col <- paste0('gray', 0:100)
						plot(hs, col = col, legend = FALSE, axes = FALSE, ext = studyRegionFocus, box = FALSE, maxcell = maxcell, mar = c(0.4, 0.1, 0.1, 0.1))
						# plot(sqPredictions[['northwest']], legend = FALSE, axes = FALSE, ext = studyRegionFocus, box = FALSE, maxcell = maxcell, mar = mar)
						
						for (region in regions) {
						
							x <- deltaPredictions[[paste(ssp, period)]][[region]]
							
							plot(x, col = deltaCol, alpha = 0.65, legend = FALSE, axes = FALSE, ext = studyRegionFocus, box = FALSE, maxcell = maxcell, range = ranges, add = TRUE)
							
						}
						
						plot(focusCounties, lwd = 0.3, add = TRUE)
						plot(studyRegionFocus, add = TRUE, lwd = 0.7)
						plot(nm, lwd = 0.7, add = TRUE)

						plot(contours, lwd = 0.4, add = TRUE)

						if (i == 4) {

							legendGrad(
								x = 'bottom',
								inset = -0.02,
								vert = FALSE,
								width = 1,
								height = 0.06,
								adjX = c(-0.8, 1),
								adjY = c(-0.1, 1),
								col = deltaCol,
								title = 'Change',
								titleCex = 0.7,
								titleAdj = c(-0.92, 0.45),
								labels = '',
								boxBorder = NA,
								lwd = 0.8,
								xpd = NA
							)
							
							niceRanges <- roundTo(ranges, 0.01)
							
							y <- -1065937
							text(-825000, y, labels = niceRanges[1], cex = 0.63, xpd = NA)
							text(-575000, y, labels = 0, cex = 0.63, xpd = NA)
							text(-440000, y, labels = paste0('+', niceRanges[2]), cex = 0.63, xpd = NA)
							
						}
						
						labelFig(label = main, adj = c(-0.05, -0.015), cex = 0.75)
						
						i <- i + 1
						
					} # next period
				
				} # next ssp

			dev.off()

# say('##################################')
# say('### cross-validation of models ###')
# say('##################################')

	# ### Divide data by region/ordinal site status then evaluate binary and ordinal models.

	# ### user-defined
	# ################
	
		# # seed
		# set.seed(1)
	
		# # proportion of sites used for validation
		# propTest <- 1 / 3 # NB smallest number of cases x region is 9
		
		# # number of k-folds
		# nFolds <- 100
	
		# # formula of most-supported models obtained manually from analyses above
		# bestOrdinalForm <- latestOccStatus ~ occVar_chronicCold_C_10yrWindow + occVar_gsPpt_mm_10yrWindow + numHomeRangesScaled + region + meanDistToClosest4Patches
		# bestBinaryForm <- presAbs ~ occVar_chronicHeat_C_10yrWindow + occVar_monsoonPpt_mm_10yrWindow + numHomeRangesScaled + region + meanDistToClosest4Patches

	# ### pika data
	# #############
		
		# load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')

		# pika$latestOccStatus <- factor(pika$latestOccStatus, levels=c('0 never', '1 old', '2 occupied'), ordered=TRUE)
		# pika$region <- as.factor(pika$region)
		
		# pika$meanDistToClosest4Patches <- log10(pika$meanDistToClosest4Patches)
		# pika$numHomeRangesScaled <- as.numeric(scale(log10(pika$numHomeRanges + 1)))
		
		# vars <- getVars('occupancy')
		# vars <- c(vars, 'meanDistToClosest4Patches')
		# varsScaled <- scale(pika[ , vars])
		# varsCenter <- attr(varsScaled, 'scaled:center')
		# varsScale <- attr(varsScaled, 'scaled:scale')
		# pika[ , vars] <- varsScaled

	# ### cross-validation
	# ####################
	
		# xvalid <- data.table()
		
		# for (fold in 1:nFolds) {

			# # compile data
			# trainData <- testData <- data.table()
			# for (region in c('northwest', 'southwest', 'northeast', 'southeast')) {
			
				# for (status in c('0 never', '1 old', '2 occupied')) {
				
					# thisData <- pika[pika$region == region & pika$latestOccStatus == status, ]
					# nTrain <- round((1 - propTest) * nrow(thisData))
					# trainIndex <- sample(nrow(thisData), nTrain)
					# thisTrainData <- thisData[trainIndex, ]
					# thisTestData <- thisData[1:nrow(thisData) %notin% trainIndex, ]
					
					# trainData <- rbind(trainData, thisTrainData)
					# testData <- rbind(testData, thisTestData)
					
				# }
				
			# }
		
			# # model and predict
			# ordinalModel <- polr(bestOrdinalForm, data = trainData, Hess = TRUE)
			# binaryModel <- glm(bestBinaryForm, data = trainData, family = binomial)

			# predOrdinal <- predict(ordinalModel, testData, type = 'p')
			# predBinary <- predict(binaryModel, testData, type = 'response')
	
			# # binary AUC
			# predBinaryPres <- predBinary[testData$presAbs == 1]
			# predBinaryAbs <- predBinary[testData$presAbs == 0]
			# aucBinary <- evalAUC(predBinaryPres, predBinaryAbs)
			# cbiBinary <- evalContBoyce(predBinaryPres, predBinaryAbs)

			# # ordinal AUC by site status
			# lower <- '0 never'
			# upper <- '1 old'
			# predOrdinalLower <- predOrdinal[testData$latestOccStatus == lower, lower]
			# predOrdinalUpper <- predOrdinal[testData$latestOccStatus == upper, upper]
			# aucOrdinalNeverVsOld <- evalAUC(predOrdinalUpper, predOrdinalLower)
			# cbiOrdinalNeverVsOld <- evalContBoyce(predOrdinalUpper, predOrdinalLower)

			# lower <- '1 old'
			# upper <- '2 occupied'
			# predOrdinalLower <- predOrdinal[testData$latestOccStatus == lower, lower]
			# predOrdinalUpper <- predOrdinal[testData$latestOccStatus == upper, upper]
			# aucOrdinalOldVsOcc <- evalAUC(predOrdinalUpper, predOrdinalLower)
			# cbiOrdinalOldVsOcc <- evalContBoyce(predOrdinalUpper, predOrdinalLower)

			# lower <- '0 never'
			# upper <- '2 occupied'
			# predOrdinalLower <- predOrdinal[testData$latestOccStatus == lower, lower]
			# predOrdinalUpper <- predOrdinal[testData$latestOccStatus == upper, upper]
			# aucOrdinalNeverVsOcc <- evalAUC(predOrdinalUpper, predOrdinalLower)
			# cbiOrdinalNeverVsOcc <- evalContBoyce(predOrdinalUpper, predOrdinalLower)

			# lower <- c('0 never', '1 old')
			# upper <- '2 occupied'
			# predOrdinalLower <- predOrdinal[testData$latestOccStatus %in% lower, colnames(predOrdinal) %in% lower]
			# predOrdinalUpper <- predOrdinal[testData$latestOccStatus == upper, upper]
			# aucOrdinalNeverAndOldVsOcc <- evalAUC(predOrdinalUpper, predOrdinalLower)
			# cbiOrdinalNeverAndOldVsOcc <- evalContBoyce(predOrdinalUpper, predOrdinalLower)
		
			# xvalid <- rbind(
				# xvalid,
				# data.table(

					# fold = fold,
					
					# aucBinary = aucBinary,
					# aucOrdinalNeverAndOldVsOcc = aucOrdinalNeverAndOldVsOcc,
					
					# cbiBinary = cbiBinary,
					# cbiOrdinalNeverAndOldVsOcc = cbiOrdinalNeverAndOldVsOcc,
					
					# aucOrdinalNeverVsOld = aucOrdinalNeverVsOld,
					# cbiOrdinalNeverVsOld = cbiOrdinalNeverVsOld,
					
					# aucOrdinalOldVsOcc = aucOrdinalOldVsOcc,
					# cbiOrdinalOldVsOcc = cbiOrdinalOldVsOcc,

					# aucOrdinalNeverVsOld = aucOrdinalNeverVsOld,
					# cbiOrdinalNeverVsOld = cbiOrdinalNeverVsOld
					
				# )
			# )

		# }
		
	# write.csv(xvalid, './Figures & Tables/Occupancy - Simple Models/Occupancy - Cross-validation.csv', row.names = FALSE)

	# sink('./Figures & Tables/Occupancy - Simple Models/Occupancy - Cross-validation Summary.txt', split = TRUE)
	# say('Summary of cross-validation results')
	# say(date(), post = 1)
	# print(colMeans(xvalid))
	# sink()

say('DONE!!!', level=1, deco='%')
