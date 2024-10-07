### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('C:/Ecology/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/03a New Mexico Pika Occupancy & Abundance Analysis - Occupancy with Simple Models.r')
### source('E:/Adam/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/03a New Mexico Pika Occupancy & Abundance Analysis - Occupancy with Simple Models.r')
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
### cross-validation of models ###

#############
### setup ###
#############

	rm(list=ls())

	drive <- 'C:/Ecology/'
	# drive <- 'E:/Adam/'

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

# 	# rank variables by mean AICc weight

# 	for (occWindow in c(occWindows_y, NA)) {
	
# 		nice <- if (is.na(occWindow)) {
# 			paste0(paste(occWindows_y, collapse=' & '), '-yr Windows')
# 		} else {
# 			paste0(occWindow, '-yr Window')
# 		}
	
# 		# get variables
# 		vars <- getVars('occupancy')
# 		if (!is.na(occWindow)) vars <- vars[grepl(vars, pattern=paste0(occWindow, 'yrWindow'))]
# 		vars <- c(vars, 'meanDistToClosest4Patches')

# 		# get models
# 		file <- paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Ordinal Models Using All Data - ', nice, '.csv')
# 		models <- read.csv(file)
		
# 		imp <- data.frame()
# 		for (var in vars) {
		
# 			index <- if (var != 'meanDistToClosest4Patches') {
# 				which(grepl(models$model, pattern=var))
# 			} else {
# 				which(!is.na(models$isolationCoeff))
# 			}
			
# 			n <- length(index)
# 			sumWeight <- sum(models$weight[index])
# 			meanWeight <- sumWeight / n

# 			imp <- rbind(
# 				imp,
# 				data.frame(
# 					variable = var,
# 					niceVar = makeNiceVars(var, 'occupancy'),
# 					numModels = n,
# 					sumWeight = sumWeight,
# 					meanWeight = meanWeight
# 				)
# 			)
			
# 		}
		
# 		imp <- imp[order(imp$meanWeight, decreasing=TRUE), ]

# 		imp$niceVar[imp$variable == 'meanDistToClosest4Patches'] <- 'isolation'

# 		write.csv(imp, paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Ordinal Models Using All Data - ', nice, ' - Var Import.csv'), row.names=FALSE)
	
# 	} # next occupancy window

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

# 	# rank variables by mean AICc weight
# 	for (occWindow in c(occWindows_y, NA)) {
	
# 		nice <- if (is.na(occWindow)) {
# 			paste0(paste(occWindows_y, collapse=' & '), '-yr Windows')
# 		} else {
# 			paste0(occWindow, '-yr Window')
# 		}
	
# 		# get variables
# 		vars <- getVars('occupancy')
# 		if (!is.na(occWindow)) vars <- vars[grepl(vars, pattern=paste0(occWindow, 'yrWindow'))]
# 		vars <- c(vars, 'meanDistToClosest4Patches')

# 		# get models
# 		file <- paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Binary Models Using All Data - ', nice, '.csv')
# 		models <- read.csv(file)
		
# 		imp <- data.frame()
# 		for (var in vars) {
		
# 			index <- if (var != 'meanDistToClosest4Patches') {
# 				which(grepl(models$model, pattern=var))
# 			} else {
# 				which(!is.na(models$isolationCoef))
# 			}
# 			n <- length(index)
# 			sumWeight <- sum(models$weight[index])
# 			meanWeight <- sumWeight / n
			
# 			imp <- rbind(
# 				imp,
# 				data.frame(
# 					variable = var,
# 					niceVar = makeNiceVars(var, 'occupancy'),
# 					numModels = n,
# 					sumWeight = sumWeight,
# 					meanWeight = meanWeight
# 				)
# 			)
			
# 		}
		
# 		imp <- imp[order(imp$meanWeight, decreasing=TRUE), ]
# 		imp$niceVar[imp$variable == 'meanDistToClosest4Patches'] <- 'isolation'

# 		write.csv(imp, paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Binary Models Using All Data - ', nice, ' - Var Import.csv'), row.names=FALSE)
	
# 	} # next occupancy window

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

say('##############################################################################################')
say('### make maps of predicted probabilities of each class from best ordinal and binary models ###')
say('##############################################################################################')

	### user-defined
	################
	
		# Dates across which to calculate climate variables. The *months* will be used, so, for example, 2019-08-01 means to use August of 2019 (so is the same as 2019-08-31).
		beginEndDates <- c('2009-09-01', '2019-08-01')
	
		# assume all regions are this region when making predictions
		region <- 'northwest'

		# home folder in which PRISM rasters are stored... should end in "/an81"
		prDir <- 'E:/Ecology/PRISM/working/an81' # HAL9000
		# prDir <- 'F:/Ecology/PRISM/working/an81' # GRN
	
		# formula of most-supported models obtained manually from analyses above
		bestOrdinalForm <- latestOccStatus ~ occVar_chronicCold_C_10yrWindow + occVar_gsPpt_mm_10yrWindow + numHomeRangesScaled + region + meanDistToClosest4Patches

		bestBinaryForm <- presAbs ~ occVar_chronicHeat_C_10yrWindow + occVar_monsoonPpt_mm_10yrWindow + numHomeRangesScaled + region + meanDistToClosest4Patches

	### pika data
	#############
		
		load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')

		pika$latestOccStatus <- factor(pika$latestOccStatus, levels=c('0 never', '1 old', '2 occupied'), ordered=TRUE)
		pika$region <- as.factor(pika$region)
		
		pika$meanDistToClosest4Patches <- log10(pika$meanDistToClosest4Patches)
		pika$numHomeRangesScaled <- as.numeric(scale(log10(pika$numHomeRanges + 1)))
		
		vars <- getVars('occupancy')
		vars <- c(vars, 'meanDistToClosest4Patches')
		varsScaled <- scale(pika[ , vars])
		varsCenter <- attr(varsScaled, 'scaled:center')
		varsScale <- attr(varsScaled, 'scaled:scale')
		pika[ , vars] <- varsScaled

	### models
	##########
	
		bestOrdinalModel <- polr(bestOrdinalForm, data=pika, Hess=TRUE)
		bestBinaryModel <- glm(bestBinaryForm, data=pika, family=binomial)

	### plot extent
	###############
		
		pikaVectUnproj <- vect(as.matrix(pika[ , ll]), 'points', crs=getCRS('wgs84'))
		pikaBuffUnprojSm <- terra::buffer(pikaVectUnproj, width = 40000)

	### create predictor rasters
	############################
	
		devtools::load_all(paste0(drive, '/R/airUpThere'))
		
		# chronic cold
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
		months <- c(11, 12, 1, 2, 3)
		months <- prefix(months, 2)
		inFocalTimePeriod <- which(substr(names, 24, 25) %in% months)
		x <- x[[inFocalTimePeriod]]
		
		x <- crop(x, pikaBuffUnprojSm)
		occVar_chronicCold_C_10yrWindow <- mean(x)
		names(occVar_chronicCold_C_10yrWindow) <- 'occVar_chronicCold_C_10yrWindow'
		
		# chronic heat
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
		
		x <- crop(x, pikaBuffUnprojSm)
		occVar_chronicHeat_C_10yrWindow <- mean(x)
		names(occVar_chronicHeat_C_10yrWindow) <- 'occVar_chronicHeat_C_10yrWindow'
		
		# GS precipitation
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

		x <- crop(x, pikaBuffUnprojSm)
		occVar_gsPpt_mm_10yrWindow <- sum(x) / 10
		names(occVar_gsPpt_mm_10yrWindow) <- 'occVar_gsPpt_mm_10yrWindow'
				
		# monsoon precipitation
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

		x <- crop(x, pikaBuffUnprojSm)
		occVar_monsoonPpt_mm_10yrWindow <- sum(x) / 10
		names(occVar_monsoonPpt_mm_10yrWindow) <- 'occVar_monsoonPpt_mm_10yrWindow'
		
		env <- c(occVar_chronicCold_C_10yrWindow, occVar_chronicHeat_C_10yrWindow, occVar_gsPpt_mm_10yrWindow, occVar_monsoonPpt_mm_10yrWindow)
		envDF <- as.data.frame(env, cells = TRUE)

		# scale
		vars <- names(envDF)
		vars <- vars[vars != 'cell']
		for (var in vars) {
			mu <- varsCenter[[var]]
			sigma <- varsScale[[var]]
			envDF[ , var] <- (envDF[ , var] - mu) / sigma
		}		

		envDF$numHomeRangesScaled <- median(pika$numHomeRangesScaled)
		envDF$meanDistToClosest4Patches <- median(pika$meanDistToClosest4Patches)
		envDF$region <- region

	### make predictions
	####################

		predOrdinal <- predict(bestOrdinalModel, envDF, type = 'p')
		predBinary <- predict(bestBinaryModel, envDF, type = 'response')

	### assign predicted values to rasters
	######################################

		template <- rast(paste0(drive, '/Research Data/PRISM/PRISM_us_dem_800m.tif'))
		template <- crop(template, pikaBuffUnprojSm)

		ordinalRast_noEvidence <- setValueByCell(template, val = predOrdinal[ , '0 never'], cell = envDF$cell)
		ordinalRast_oldEvidence <- setValueByCell(template, val = predOrdinal[ , '1 old'], cell = envDF$cell)
		ordinalRast_occupied <- setValueByCell(template, val = predOrdinal[ , '2 occupied'], cell = envDF$cell)

		binaryRast <- template
		binaryRast[] <- predBinary

		ordinalRast <- c(ordinalRast_noEvidence, ordinalRast_oldEvidence, ordinalRast_occupied)
		names(ordinalRast) <- c('no evidence', 'previously occupied', 'occupied')

		binaryRast <- project(binaryRast, getCRS('NA Albers'))
		ordinalRast <- project(ordinalRast, getCRS('NA Albers'))

	### hillshade
	#############

		elev <- rast('./Data/elev_fine_m.tif')
		elev <- crop(elev, pikaBuffUnprojSm)
		
		slope <- terrain(elev, 'slope', unit = 'radians')
		aspect <- terrain(elev, 'aspect', unit = 'radians')

		hs <- shade(slope, aspect, direction = 45)
		hs <- project(hs, getCRS('NA Albers'))

	### plot
	########

		pikaVect <- vect(pika, geom = c('longitude', 'latitude'), crs = getCRS('WGS84'))
		pikaVect <- project(pikaVect, getCRS('NA Albers'))
		pikaVect$presAbs <- as.character(pikaVect$presAbs)

		pikaCols <- rep(NA_character_, nrow(pikaVect))
		pikaCols[pikaVect$latestOccStatus == '0 never'] <- 'red'
		pikaCols[pikaVect$latestOccStatus == '1 old'] <- 'orange'
		pikaCols[pikaVect$latestOccStatus == '2 occupied'] <- 'green'

		usa <- vect(paste0(drive, './Research Data/GADM/Version 4.1/High Res North America Level 2 sans Great Lakes SpatVector WGS84.gpkg'))
		focusCounties <- crop(usa, ext(pikaBuffUnprojSm))
		focusCounties <- project(focusCounties, getCRS('NA Albers'))

		binaryMap <- ggplot() +
			layer_spatial(hs, dpi = 600) +
			layer_spatial(binaryRast, dpi = 600, alpha = 0.7) +
			scale_fill_continuous(
				name = 'Prediction',
				na.value = NA,
				low = 'red',
				high = 'green3',
				limits = c(0, 1)
			) +
			layer_spatial(focusCounties, fill = NA) +
			ggtitle('(a) Binary model predictions') +
			theme_bw() +
			theme(
				panel.border = element_blank(),
				plot.title = element_text(size = 10),
				legend.title = element_text(size = 9),
				legend.text = element_text(size = 8),
				axis.text = element_text(size = 5),
				axis.ticks = element_blank()
			)		
		
		pikaVect <- pikaVect[order(pikaVect$presAbs)]
		presAbsSiteMap <- ggplot() +
			layer_spatial(hs, dpi = 600) +
			scale_fill_gradient(
				low = 'gray0',
				high = 'gray80',
				guide = 'none',
				na.value = NA
			) +
			layer_spatial(focusCounties, fill = NA) +
			layer_spatial(pikaVect, aes(color = presAbs, shape = presAbs), size = 0.8) +
			scale_shape_manual(
				name = 'Site\nstatus',
				labels = c('Absent', 'Present'),
				values = c('0' = 4, '1' = 1)
			) +
			scale_color_manual(
				name = 'Site\nstatus',
				labels = c('Absent', 'Present'),
				values = c('0' = 'red', '1' = 'green3')
			) +
			ggtitle('(b) Binary site status') +
			theme_bw() +
			theme(
				panel.border = element_blank(),
				plot.title = element_text(size = 10),
				axis.text = element_text(size = 5),
				axis.ticks = element_blank(),
				legend.title = element_text(size = 8),
				legend.text = element_text(size = 7)
			)		

		ordinalMaps <- ordinalSiteMaps <- list()
		low <- NA
		for (i in 1:3) {
			
			if (i == 1) {
				x <- ordinalRast_noEvidence
				title <- 'No evidence'
				values = c('0 never' = 'red', '1 old' = 'gray30', '2 occupied' = 'gray30')
				high <- 'red'
				pikaVect <- pikaVect[order(pikaVect$latestOccStatus, decreasing = TRUE)]
			} else if (i == 2) {
				x <- ordinalRast_oldEvidence
				title <- 'Previously occupied'
				values = c('0 never' = 'gray30', '1 old' = 'orange', '2 occupied' = 'gray30')
				high <- 'orange'
				pikaVect <- rbind(pikaVect[pikaVect$latestOccStatus != '1 old'], pikaVect[pikaVect$latestOccStatus == '1 old'])
			} else {
				x <- ordinalRast_occupied
				title <- 'Occupied'
				values = c('0 never' = 'gray30', '1 old' = 'gray30', '2 occupied' = 'green3')
				high <- 'green3'
				pikaVect <- pikaVect[order(pikaVect$latestOccStatus)]
			}
			
			predictionTitle <- paste0('(', letters[2 * i + 1], ') Ordinal model predictions:\n      ', title)
			siteTitle <- paste0('(', letters[2 * i + 2], ') Ordinal site status:\n      ', title, ' highlighted')
				
			ordinalMaps[[i]] <- ggplot() +
				layer_spatial(hs, dpi = 600) +
				layer_spatial(x, dpi = 600, alpha = 0.7) +
				scale_fill_continuous(
					name = 'Prediction',
					na.value = NA,
					low = low,
					high = high,
					limits = c(0, 1)
				) +
				layer_spatial(focusCounties, fill = NA) +
				# layer_spatial(pikaVect, pch = 1, aes(color = latestOccStatus)) +
				# scale_color_manual(
					# name = 'Site\nstatus',
					# labels = c('No\n  evidence', 'Previously\n  occupied', 'Occupied'),
					# values = c('0 never' = 'red', '1 old' = 'orange', '2 occupied' = 'green3')
				# ) +
				ggtitle(predictionTitle)
				
			ordinalSiteMaps[[i]] <- ggplot() +
				layer_spatial(hs, dpi = 600) +
				scale_fill_gradient(
					low = 'gray0',
					high = 'gray80',
					guide = 'none',
					na.value = NA
				) +
				layer_spatial(focusCounties, fill = NA) +
				layer_spatial(pikaVect, aes(color = latestOccStatus, shape = latestOccStatus), size = 0.8) +
				scale_shape_manual(
					name = 'Site\nstatus',
					labels = c('No\n  evidence', 'Previously\n    occupied', 'Occupied'),
					values = c('0 never' = 4, '1 old' = 2, '2 occupied' = 1)
				) +
				scale_color_manual(
					name = 'Site\nstatus',
					labels = c('No\n  evidence', 'Previously\n    occupied', 'Occupied'),
					values = values
				) +
				ggtitle(siteTitle) +
				theme_bw() +
				theme(
					panel.border = element_blank(),
					plot.title = element_text(size = 10),
					axis.text = element_text(size = 5),
					axis.ticks = element_blank()
				)

				ordinalMaps[[i]] <- ordinalMaps[[i]] +
					theme_bw() +
					theme(
						panel.border = element_blank(),
						plot.title = element_text(size = 10),
						legend.title = element_text(size = 9),
						legend.text = element_text(size = 8),
						axis.text = element_text(size = 5),
						axis.ticks = element_blank()
					)			

				ordinalSiteMaps[[i]] <- ordinalSiteMaps[[i]] +
					theme_bw() +
					theme(
						panel.border = element_blank(),
						plot.title = element_text(size = 10),
						axis.text = element_text(size = 5),
						axis.ticks = element_blank(),
						legend.title = element_text(size = 8),
						legend.text = element_text(size = 7)
					)			

		}
		
		mapsArranged <- list(
			binaryMap, presAbsSiteMap,
			ordinalMaps[[1]], ordinalSiteMaps[[1]],
			ordinalMaps[[2]], ordinalSiteMaps[[2]],
			ordinalMaps[[3]], ordinalSiteMaps[[3]]
		)
		
		maps <- plot_grid(plotlist = mapsArranged, ncol = 2)

		region <- capIt(region)
		ggsave(maps, file = paste0('./Figures & Tables/Predictions Assuming ', region, '.png'), dpi = 600, width = 7, height = 11, bg = 'white')

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
