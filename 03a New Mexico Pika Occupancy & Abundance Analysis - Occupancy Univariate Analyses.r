### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/Code/03a New Mexico Pika Occupancy & Abundance Analysis - Occupancy Univariate Analyses.r')
###
### CONTENTS ###
### setup ###
### ORDINAL univariate OCCUPANCY analysis: all-data models ###
### BINARY univariate OCCUPANCY analysis: all-data models ###
### ORDINAL univariate OCCUPANCY analysis: cross-validated models ###
### BINARY univariate OCCUPANCY analysis: cross-validated models ###
### ORDINAL univariate OCCUPANCY analysis: cross-validated models summary ###
### BINARY univariate OCCUPANCY analysis: cross-validated models summary ###
### ORDINAL multivariate OCCUPANCY analysis: all-data models ###

#############
### setup ###
#############

	# drive <- 'C:'
	drive <- 'D:'
	# drive <- 'E:'

	source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/Code/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r')

say('##############################################################')
say('### ORDINAL univariate OCCUPANCY analysis: all-data models ###')
say('##############################################################')

	load('./Data/03 New Mexico Pika - Assigned Folds.rda')
	pika$latestOccStatus <- factor(pika$latestOccStatus, levels=c('0 never', '1 old', '2 occupied'), ordered=TRUE)
	pika$region <- as.factor(pika$region)
	
	pika$numHomeRangesScaled <- scale(pika$numHomeRanges)

	vars <- names(pika)[grepl(names(pika), pattern='occVar_')]

	# predictors
	file <- './Predictors/Predictor Variables 2021-04-28 Edits by Adam.xlsx'
	univariatePredTable <- read_excel(file, sheet='Predictor Variables', skip=1)

	### using all sites (no cross-validation)
	#########################################
	
		### null models
		###############

			model1 <- polr(latestOccStatus ~ 1, data=pika, Hess=TRUE)
			model2 <- polr(latestOccStatus ~ numHomeRangesScaled, data=pika, Hess=TRUE)
			model3 <- polr(latestOccStatus ~ region, data=pika, Hess=TRUE)
			model4 <- polr(latestOccStatus ~ numHomeRangesScaled + region, data=pika, Hess=TRUE)

			aicc1 <- AICc(model1)
			aicc2 <- AICc(model2)
			aicc3 <- AICc(model3)
			aicc4 <- AICc(model4)

			coeffs <- c(NA, NA, NA, NA)
			numHomeRanges <- c(FALSE, TRUE, FALSE, TRUE)
			region <- c(FALSE, FALSE, TRUE, TRUE)
			
			aicc <- c(aicc1, aicc2, aicc3, aicc4)

			like1 <- logLik(model1)
			like2 <- logLik(model2)
			like3 <- logLik(model3)
			like4 <- logLik(model4)
			
			likeNull <- like1
			
			n <- nrow(pika)
			pseudoR2_1 <- nagelR2(likeNull, like1, n)
			pseudoR2_2 <- nagelR2(likeNull, like2, n)
			pseudoR2_3 <- nagelR2(likeNull, like3, n)
			pseudoR2_4 <- nagelR2(likeNull, like4, n)
			
			pseudoR2 <- c(pseudoR2_1, pseudoR2_2, pseudoR2_3, pseudoR2_4)

			results <- data.frame(
				var = '(Intercept)',
				coeff = coeffs,
				numHomeRanges = numHomeRanges,
				region = region,
				aicc = aicc,
				pseudoR2 = pseudoR2
			)
			

		### by climate variable
		#######################
		
		for (var in vars) {
			
			x <- scale(pika[ , var])
			
			model1 <- polr(latestOccStatus ~ x, data=pika, Hess=TRUE)
			model2 <- polr(latestOccStatus ~ x + numHomeRangesScaled, data=pika, Hess=TRUE)
			model3 <- polr(latestOccStatus ~ x + region, data=pika, Hess=TRUE)
			model4 <- polr(latestOccStatus ~ x + numHomeRangesScaled + region, data=pika, Hess=TRUE)
			
			aicc1 <- AICc(model1)
			aicc2 <- AICc(model2)
			aicc3 <- AICc(model3)
			aicc4 <- AICc(model4)
			
			coeffs <- c(
				coefficients(model1)[['x']],
				coefficients(model2)[['x']],
				coefficients(model3)[['x']],
				coefficients(model4)[['x']]
			)
			
			numHomeRanges <- c(FALSE, TRUE, FALSE, TRUE)
			region <- c(FALSE, FALSE, TRUE, TRUE)
			
			aicc <- c(aicc1, aicc2, aicc3, aicc4)
			
			like1 <- logLik(model1)
			like2 <- logLik(model2)
			like3 <- logLik(model3)
			like4 <- logLik(model4)
			
			pseudoR2_1 <- nagelR2(likeNull, like1, n)
			pseudoR2_2 <- nagelR2(likeNull, like2, n)
			pseudoR2_3 <- nagelR2(likeNull, like3, n)
			pseudoR2_4 <- nagelR2(likeNull, like4, n)
			
			pseudoR2 <- c(pseudoR2_1, pseudoR2_2, pseudoR2_3, pseudoR2_4)
			
			results <- rbind(
				results,
				data.frame(
					var = var,
					coeff = coeffs,
					numHomeRanges = numHomeRanges,
					region = region,
					aicc = aicc,
					pseudoR2 = pseudoR2
				)
			)

		} # next variable

		### reports
		###########
		
			for (occWindow in c(occWindows_y, NA)) {
			
				if (is.na(occWindow)) {
					thisResults <- results
					nice <- paste0(paste(occWindows_y, collapse=' & '), '-yr Windows')
				} else {

					thisResults <- rbind(
						results[grepl(results$var, pattern=paste0(occWindow, 'yrWindow')), ],
						results[results$var == '(Intercept)', ]
					)
						
					nice <- paste0(occWindow, '-yr Window')
				}
				
				thisResults$deltaAicc <- thisResults$aicc - min(thisResults$aicc)
				w <- exp(-0.5 * thisResults$deltaAicc)
				thisResults$weight <- w / sum(w)

				thisResults <- thisResults[order(thisResults$weight, decreasing=TRUE), ]
				rownames(thisResults) <- NULL

				file <- paste0('./Figures & Tables/Occupancy - Univariate/Occupancy - Univariate Ordinal Models Using All Data - ', nice, '.csv')
				write.csv(thisResults, file, row.names=FALSE)
				
			} # next window
				
		### model diagnostics
		#####################

say('##############################################################')
say('### BINARY univariate OCCUPANCY analysis: all-data models ###')
say('##############################################################')

	load('./Data/03 New Mexico Pika - Assigned Folds.rda')
	pika$region <- as.factor(pika$region)
	
	pika$numHomeRangesScaled <- scale(pika$numHomeRanges)

	vars <- names(pika)[grepl(names(pika), pattern='occVar_')]

	# predictors
	file <- './Predictors/Predictor Variables 2021-04-28 Edits by Adam.xlsx'
	univariatePredTable <- read_excel(file, sheet='Predictor Variables', skip=1)
	
	### null models
	###############

		model1 <- glm(presAbs ~ 1, data=pika, family=binomial)
		model2 <- glm(presAbs ~ numHomeRangesScaled, data=pika, family=binomial)
		model3 <- glm(presAbs ~ region, data=pika, family=binomial)
		model4 <- glm(presAbs ~ numHomeRangesScaled + region, data=pika, family=binomial)

		aicc1 <- AICc(model1)
		aicc2 <- AICc(model2)
		aicc3 <- AICc(model3)
		aicc4 <- AICc(model4)

		coeffs <- c(NA, NA, NA, NA)
		numHomeRanges <- c(FALSE, TRUE, FALSE, TRUE)
		region <- c(FALSE, FALSE, TRUE, TRUE)
		
		aicc <- c(aicc1, aicc2, aicc3, aicc4)
		
		like1 <- logLik(model1)
		like2 <- logLik(model2)
		like3 <- logLik(model3)
		like4 <- logLik(model4)
		
		likeNull <- like1
		
		n <- nrow(pika)
		pseudoR2_1 <- nagelR2(likeNull, like1, n)
		pseudoR2_2 <- nagelR2(likeNull, like2, n)
		pseudoR2_3 <- nagelR2(likeNull, like3, n)
		pseudoR2_4 <- nagelR2(likeNull, like4, n)
		
		pseudoR2 <- c(pseudoR2_1, pseudoR2_2, pseudoR2_3, pseudoR2_4)

		results <- data.frame(
			var = '(Intercept)',
			coeff = coeffs,
			numHomeRanges = numHomeRanges,
			region = region,
			aicc = aicc,
			pseudoR2 = pseudoR2
		)

	### by climate variable
	#######################

		for (var in vars) {
			
			x <- scale(pika[ , var])
			
			model1 <- glm(presAbs ~ x, data=pika, family=binomial)
			model2 <- glm(presAbs ~ x + numHomeRangesScaled, data=pika, family=binomial)
			model3 <- glm(presAbs ~ x + region, data=pika, family=binomial)
			model4 <- glm(presAbs ~ x + numHomeRangesScaled + region, data=pika, family=binomial)
			
			aicc1 <- AICc(model1)
			aicc2 <- AICc(model2)
			aicc3 <- AICc(model3)
			aicc4 <- AICc(model4)
			
			coeffs <- c(
				coefficients(model1)[['x']],
				coefficients(model2)[['x']],
				coefficients(model3)[['x']],
				coefficients(model4)[['x']]
			)
			
			numHomeRanges <- c(FALSE, TRUE, FALSE, TRUE)
			region <- c(FALSE, FALSE, TRUE, TRUE)
			
			aicc <- c(aicc1, aicc2, aicc3, aicc4)
			
			like1 <- logLik(model1)
			like2 <- logLik(model2)
			like3 <- logLik(model3)
			like4 <- logLik(model4)
			
			pseudoR2_1 <- nagelR2(likeNull, like1, n)
			pseudoR2_2 <- nagelR2(likeNull, like2, n)
			pseudoR2_3 <- nagelR2(likeNull, like3, n)
			pseudoR2_4 <- nagelR2(likeNull, like4, n)
			
			pseudoR2 <- c(pseudoR2_1, pseudoR2_2, pseudoR2_3, pseudoR2_4)

			results <- rbind(
				results,
				data.frame(
					var = var,
					coeff = coeffs,
					numHomeRanges = numHomeRanges,
					region = region,
					aicc = aicc,
					pseudoR2 = pseudoR2
				)
			)
			
		} # next variable

	### reports
	###########
	
		for (occWindow in c(occWindows_y, NA)) {
		
			if (is.na(occWindow)) {
				thisResults <- results
				nice <- paste0(paste(occWindows_y, collapse=' & '), '-yr Windows')
			} else {
				thisResults <- results[grepl(results$var, pattern=paste0(occWindow, 'yrWindow')), ]
				nice <- paste0(occWindow, '-yr Window')
			}
			
			thisResults$deltaAicc <- thisResults$aicc - min(thisResults$aicc)
			w <- exp(-0.5 * thisResults$deltaAicc)
			thisResults$weight <- w / sum(w)

			thisResults <- thisResults[order(thisResults$weight, decreasing=TRUE), ]
			rownames(thisResults) <- NULL

			file <- paste0('./Figures & Tables/Occupancy - Univariate/Occupancy - Univariate Binary Models Using All Data - ', nice, '.csv')
			write.csv(thisResults, file, row.names=FALSE)
			
		} # next window
			
	### model diagnostics
	#####################

say('#####################################################################')
say('### ORDINAL univariate OCCUPANCY analysis: cross-validated models ###')
say('#####################################################################')

	load('./Data/03 New Mexico Pika - Assigned Folds.rda')
	pika$latestOccStatus <- factor(pika$latestOccStatus, levels=c('0 never', '1 old', '2 occupied'), ordered=TRUE)
	pika$region <- as.factor(pika$region)
	pika$numHomeRangesScaled <- scale(pika$numHomeRanges)

	vars <- names(pika)[grepl(names(pika), pattern='occVar_')]

	folds <- expand.grid(nwFold=1:2, swFold=1:2, neFold=1:2, seFold=1:2)

	results <- data.frame()

	### intercept-only models
	#########################
	
		for (k in 1:nrow(folds)) {

			trainIndex <- which(
				(pika$region == 'nw' & pika$fold == folds$nwFold[k]) |
				(pika$region == 'sw' & pika$fold == folds$swFold[k]) |
				(pika$region == 'ne' & pika$fold == folds$neFold[k]) |
				(pika$region == 'se' & pika$fold == folds$seFold[k])
			)
			
			testIndex <- which(
				(pika$region == 'nw' & pika$fold != folds$nwFold[k]) |
				(pika$region == 'sw' & pika$fold != folds$swFold[k]) |
				(pika$region == 'ne' & pika$fold != folds$neFold[k]) |
				(pika$region == 'se' & pika$fold != folds$seFold[k])
			)
			
			trainData <- pika[trainIndex, ]
			testData <- pika[testIndex, ]

			# train models

			model1 <- polr(latestOccStatus ~ 1, data=trainData, Hess=TRUE)
			model2 <- polr(latestOccStatus ~ 1 + numHomeRangesScaled, data=trainData, Hess=TRUE)
			model3 <- polr(latestOccStatus ~ 1 + region, data=trainData, Hess=TRUE)
			model4 <- polr(latestOccStatus ~ 1 + numHomeRangesScaled + region, data=trainData, Hess=TRUE)
			
			# test models

			pred1 <- predict(model1, testData, type='probs')
			pred2 <- predict(model2, testData, type='probs')
			pred3 <- predict(model3, testData, type='probs')
			pred4 <- predict(model4, testData, type='probs')
			
			pred1_never <- pred1[testData$latestOccStatus == '0 never', 1]
			pred1_old <- pred1[testData$latestOccStatus == '1 old', 2]
			pred1_occ <- pred1[testData$latestOccStatus == '2 occupied', 3]
			auc1 <- aucMultiWeighted(pred1_never, pred1_old, pred1_occ)
			
			pred2_never <- pred2[testData$latestOccStatus == '0 never', 1]
			pred2_old <- pred2[testData$latestOccStatus == '1 old', 2]
			pred2_occ <- pred2[testData$latestOccStatus == '2 occupied', 3]
			auc2 <- aucMultiWeighted(pred2_never, pred2_old, pred2_occ)
			
			pred3_never <- pred3[testData$latestOccStatus == '0 never', 1]
			pred3_old <- pred3[testData$latestOccStatus == '1 old', 2]
			pred3_occ <- pred3[testData$latestOccStatus == '2 occupied', 3]
			auc3 <- aucMultiWeighted(pred3_never, pred3_old, pred3_occ)
			
			pred4_never <- pred4[testData$latestOccStatus == '0 never', 1]
			pred4_old <- pred4[testData$latestOccStatus == '1 old', 2]
			pred4_occ <- pred4[testData$latestOccStatus == '2 occupied', 3]
			auc4 <- aucMultiWeighted(pred4_never, pred4_old, pred4_occ)
			
			aucs <- c(auc1[['multivariate']], auc2[['multivariate']], auc3[['multivariate']], auc4[['multivariate']])
						
			results <- rbind(
				results,
				data.frame(
					var = '(Intercept)',
					timeFrame = NA,
					fold = k,
					numHomeRanges = c(FALSE, TRUE, FALSE, TRUE),
					region = c(FALSE, FALSE, TRUE, TRUE),
					auc = aucs
				)
			)
			
		} # next fold
	
	### by climate variable
	#######################
	
	for (var in vars) {
		
		xScaled <- scale(pika[ , var])

		for (k in 1:nrow(folds)) {

			trainIndex <- which(
				(pika$region == 'nw' & pika$fold == folds$nwFold[k]) |
				(pika$region == 'sw' & pika$fold == folds$swFold[k]) |
				(pika$region == 'ne' & pika$fold == folds$neFold[k]) |
				(pika$region == 'se' & pika$fold == folds$seFold[k])
			)
			
			testIndex <- which(
				(pika$region == 'nw' & pika$fold != folds$nwFold[k]) |
				(pika$region == 'sw' & pika$fold != folds$swFold[k]) |
				(pika$region == 'ne' & pika$fold != folds$neFold[k]) |
				(pika$region == 'se' & pika$fold != folds$seFold[k])
			)
			
			trainData <- pika[trainIndex, ]
			testData <- pika[testIndex, ]

			# train models
			trainData$x <- xScaled[trainIndex]

			model1 <- polr(latestOccStatus ~ x, data=trainData, Hess=TRUE)
			model2 <- polr(latestOccStatus ~ x + numHomeRangesScaled, data=trainData, Hess=TRUE)
			model3 <- polr(latestOccStatus ~ x + region, data=trainData, Hess=TRUE)
			model4 <- polr(latestOccStatus ~ x + numHomeRangesScaled + region, data=trainData, Hess=TRUE)
			
			# test models
			testData$x <- xScaled[testIndex]

			pred1 <- predict(model1, testData, type='probs')
			pred2 <- predict(model2, testData, type='probs')
			pred3 <- predict(model3, testData, type='probs')
			pred4 <- predict(model4, testData, type='probs')
			
			pred1_never <- pred1[testData$latestOccStatus == '0 never', 1]
			pred1_old <- pred1[testData$latestOccStatus == '1 old', 2]
			pred1_occ <- pred1[testData$latestOccStatus == '2 occupied', 3]
			auc1 <- aucMultiWeighted(pred1_never, pred1_old, pred1_occ)
			
			pred2_never <- pred2[testData$latestOccStatus == '0 never', 1]
			pred2_old <- pred2[testData$latestOccStatus == '1 old', 2]
			pred2_occ <- pred2[testData$latestOccStatus == '2 occupied', 3]
			auc2 <- aucMultiWeighted(pred2_never, pred2_old, pred2_occ)
			
			pred3_never <- pred3[testData$latestOccStatus == '0 never', 1]
			pred3_old <- pred3[testData$latestOccStatus == '1 old', 2]
			pred3_occ <- pred3[testData$latestOccStatus == '2 occupied', 3]
			auc3 <- aucMultiWeighted(pred3_never, pred3_old, pred3_occ)
			
			pred4_never <- pred4[testData$latestOccStatus == '0 never', 1]
			pred4_old <- pred4[testData$latestOccStatus == '1 old', 2]
			pred4_occ <- pred4[testData$latestOccStatus == '2 occupied', 3]
			auc4 <- aucMultiWeighted(pred4_never, pred4_old, pred4_occ)
			
			aucs <- c(auc1[['multivariate']], auc2[['multivariate']], auc3[['multivariate']], auc4[['multivariate']])
			
			timeFrame <- substr(var, nchar(var) - 9, nchar(var) - 8)
			if (substr(timeFrame, 1, 1) == '_') timeFrame <- substr(timeFrame, 2, 2)
			timeFrame <- as.integer(timeFrame)
			
			results <- rbind(
				results,
				data.frame(
					var = var,
					timeFrame = timeFrame,
					fold = k,
					numHomeRanges = c(FALSE, TRUE, FALSE, TRUE),
					region = c(FALSE, FALSE, TRUE, TRUE),
					auc = c(auc1[['multivariate']], auc2[['multivariate']], auc3[['multivariate']], auc4[['multivariate']])
				)
			)
			
		} # next fold
			
	} # next variable
	
	write.csv(results, file='./Figures & Tables/Occupancy - Univariate/Occupancy - Univariate Ordinal Models Cross-validation Results.csv', row.names=FALSE)

say('####################################################################')
say('### BINARY univariate OCCUPANCY analysis: cross-validated models ###')
say('####################################################################')

	load('./Data/03 New Mexico Pika - Assigned Folds.rda')
	pika$region <- as.factor(pika$region)
	
	pika$numHomeRangesScaled <- scale(pika$numHomeRanges)

	vars <- names(pika)[grepl(names(pika), pattern='occVar_')]

	folds <- expand.grid(nwFold=1:2, swFold=1:2, neFold=1:2, seFold=1:2)

	results <- data.frame()
	
	### intercept-only models
	#########################
	
		for (k in 1:nrow(folds)) {

			trainIndex <- which(
				(pika$region == 'nw' & pika$fold == folds$nwFold[k]) |
				(pika$region == 'sw' & pika$fold == folds$swFold[k]) |
				(pika$region == 'ne' & pika$fold == folds$neFold[k]) |
				(pika$region == 'se' & pika$fold == folds$seFold[k])
			)
			
			testIndex <- which(
				(pika$region == 'nw' & pika$fold != folds$nwFold[k]) |
				(pika$region == 'sw' & pika$fold != folds$swFold[k]) |
				(pika$region == 'ne' & pika$fold != folds$neFold[k]) |
				(pika$region == 'se' & pika$fold != folds$seFold[k])
			)
			
			trainData <- pika[trainIndex, ]
			testData <- pika[testIndex, ]

			# train models

			model1 <- glm(presAbs ~ 1, data=trainData, family=binomial)
			model2 <- glm(presAbs ~ 1 + numHomeRangesScaled, data=trainData, family=binomial)
			model3 <- glm(presAbs ~ 1 + region, data=trainData, family=binomial)
			model4 <- glm(presAbs ~ 1 + numHomeRangesScaled + region, data=trainData, family=binomial)
			
			# test models

			pred1 <- predict(model1, testData, type='response')
			pred2 <- predict(model2, testData, type='response')
			pred3 <- predict(model3, testData, type='response')
			pred4 <- predict(model4, testData, type='response')
			
			predPres1 <- pred1[testData$presAbs == 1]
			predPres2 <- pred2[testData$presAbs == 1]
			predPres3 <- pred3[testData$presAbs == 1]
			predPres4 <- pred4[testData$presAbs == 1]
			
			predAbs1 <- pred1[testData$presAbs == 0]
			predAbs2 <- pred2[testData$presAbs == 0]
			predAbs3 <- pred3[testData$presAbs == 0]
			predAbs4 <- pred4[testData$presAbs == 0]
			
			auc1 <- aucWeighted(predPres1, predAbs1)
			auc2 <- aucWeighted(predPres2, predAbs2)
			auc3 <- aucWeighted(predPres3, predAbs3)
			auc4 <- aucWeighted(predPres4, predAbs4)
	
			aucs <- c(auc1, auc2, auc3, auc4)
						
			results <- rbind(
				results,
				data.frame(
					var = '(Intercept)',
					timeFrame = NA,
					fold = k,
					numHomeRanges = c(FALSE, TRUE, FALSE, TRUE),
					region = c(FALSE, FALSE, TRUE, TRUE),
					auc = aucs
				)
			)
			
		} # next fold

	### by climate variable
	#######################

	for (var in vars) {
		
		xScaled <- scale(pika[ , var])

		for (k in 1:nrow(folds)) {

			trainIndex <- which(
				(pika$region == 'nw' & pika$fold == folds$nwFold[k]) |
				(pika$region == 'sw' & pika$fold == folds$swFold[k]) |
				(pika$region == 'ne' & pika$fold == folds$neFold[k]) |
				(pika$region == 'se' & pika$fold == folds$seFold[k])
			)
			
			testIndex <- which(
				(pika$region == 'nw' & pika$fold != folds$nwFold[k]) |
				(pika$region == 'sw' & pika$fold != folds$swFold[k]) |
				(pika$region == 'ne' & pika$fold != folds$neFold[k]) |
				(pika$region == 'se' & pika$fold != folds$seFold[k])
			)
			
			trainData <- pika[trainIndex, ]
			testData <- pika[testIndex, ]

			# train models
			x <- xScaled[trainIndex]

			model1 <- glm(presAbs ~ x, data=trainData, family=binomial)
			model2 <- glm(presAbs ~ x + numHomeRangesScaled, data=trainData, family=binomial)
			model3 <- glm(presAbs ~ x + region, data=trainData, family=binomial)
			model4 <- glm(presAbs ~ x + numHomeRangesScaled + region, data=trainData, family=binomial)
			
			# test models
			testData$x <- xScaled[testIndex]

			pred1 <- predict(model1, testData, type='response')
			pred2 <- predict(model2, testData, type='response')
			pred3 <- predict(model3, testData, type='response')
			pred4 <- predict(model4, testData, type='response')
			
			predPres1 <- pred1[testData$presAbs == 1]
			predPres2 <- pred2[testData$presAbs == 1]
			predPres3 <- pred3[testData$presAbs == 1]
			predPres4 <- pred4[testData$presAbs == 1]
			
			predAbs1 <- pred1[testData$presAbs == 0]
			predAbs2 <- pred2[testData$presAbs == 0]
			predAbs3 <- pred3[testData$presAbs == 0]
			predAbs4 <- pred4[testData$presAbs == 0]
			
			auc1 <- aucWeighted(predPres1, predAbs1)
			auc2 <- aucWeighted(predPres2, predAbs2)
			auc3 <- aucWeighted(predPres3, predAbs3)
			auc4 <- aucWeighted(predPres4, predAbs4)
			
			aucs <- c(auc1, auc2, auc3, auc4)
			
			timeFrame <- substr(var, nchar(var) - 9, nchar(var) - 8)
			if (substr(timeFrame, 1, 1) == '_') timeFrame <- substr(timeFrame, 2, 2)
			timeFrame <- as.integer(timeFrame)
			
			results <- rbind(
				results,
				data.frame(
					var = var,
					timeFrame = timeFrame,
					fold = k,
					numHomeRanges = c(FALSE, TRUE, FALSE, TRUE),
					region = c(FALSE, FALSE, TRUE, TRUE),
					auc = aucs
				)
			)
			
		} # next fold
			
	} # next variable
	
	write.csv(results, file='./Figures & Tables/Occupancy - Univariate/Occupancy - Univariate Binary Models Cross-validation Results.csv', row.names=FALSE)

say('#############################################################################')
say('### ORDINAL univariate OCCUPANCY analysis: cross-validated models summary ###')
say('#############################################################################')
	
	results <- read.csv('./Figures & Tables/Occupancy - Univariate/Occupancy - Univariate Ordinal Models Cross-validation Results.csv')

	# intercept-only models
	nulls <- results[results$var == '(Intercept)', ]
	nulls <- aggregate(nulls, by=list(nulls$numHomeRanges, nulls$region), FUN=mean)
	nulls$var <- nulls$timeFrame <- nulls$fold <- nulls$numHomeRanges <- nulls$region <- NULL
	names(nulls)[1:2] <- c('numHomeRanges', 'region')
	n <- nrow(nulls)

	first <- data.frame(
		var = rep('(Intercept)', n),
		timeFrame = rep(NA, n)
	)
	
	last <- data.frame(varNice = rep('(Intercept)', n))

	nulls <- cbind(first, nulls, last)

	# climate models
	clims <- results[results$var != '(Intercept)', ]
	clims <- aggregate(clims, by=list(clims$var, clims$timeFrame, clims$numHomeRanges, clims$region), FUN=mean)
	clims$var <- clims$fold <- clims$timeFrame <- clims$numHomeRanges <- clims$region <- NULL
	names(clims)[1:4] <- c('var', 'timeFrame', 'numHomeRanges', 'region')
	clims$varNice <- makeNiceVars(clims$var, occOrDens='occ')
	
	results <- rbind(nulls, clims)
	results <- results[order(results$auc, decreasing=TRUE), ]
	results$varNice <- factor(results$varNice, levels=rev(unique(results$varNice)))

	ylim <- c(min(0.5, min(results$auc)), 1)

	p <- ggplot(data=results, aes(x=varNice, y=auc, col=region, pch=numHomeRanges)) +
		geom_point(size=2) +
		labs(shape='Number of Home\nRanges as Covariate', color='Region Covariate') +
		xlab(NULL) + ylab('Multivariate AUC') +
		# ylim(ylim[1], ylim[2]) +
		theme(axis.text.y=element_text(size=12)) +
		coord_flip()
	
	pdf('./Figures & Tables/Occupancy - Univariate/Occupancy - Univariate Ordinal Models Cross-validation Results.pdf', width=8, height=8)
		print(p)
	dev.off()
	
say('############################################################################')
say('### BINARY univariate OCCUPANCY analysis: cross-validated models summary ###')
say('############################################################################')
	
	results <- read.csv('./Figures & Tables/Occupancy - Univariate/Occupancy - Univariate Binary Models Cross-validation Results.csv')
	
	# intercept-only models
	nulls <- results[results$var == '(Intercept)', ]
	nulls <- aggregate(nulls, by=list(nulls$numHomeRanges, nulls$region), FUN=mean)
	nulls$var <- nulls$timeFrame <- nulls$fold <- nulls$numHomeRanges <- nulls$region <- NULL
	names(nulls)[1:2] <- c('numHomeRanges', 'region')
	n <- nrow(nulls)

	first <- data.frame(
		var = rep('(Intercept)', n),
		timeFrame = rep(NA, n)
	)
	
	last <- data.frame(varNice = rep('(Intercept)', n))

	nulls <- cbind(first, nulls, last)

	# climate models
	clims <- results[results$var != '(Intercept)', ]
	clims <- aggregate(clims, by=list(clims$var, clims$timeFrame, clims$numHomeRanges, clims$region), FUN=mean)
	clims$var <- clims$fold <- clims$timeFrame <- clims$numHomeRanges <- clims$region <- NULL
	names(clims)[1:4] <- c('var', 'timeFrame', 'numHomeRanges', 'region')
	clims$varNice <- makeNiceVars(clims$var, occOrDens='occ')
	
	results <- rbind(nulls, clims)
	results <- results[order(results$auc, decreasing=TRUE), ]
	results$varNice <- factor(results$varNice, levels=rev(unique(results$varNice)))

	ylim <- c(min(0.5, min(results$auc)), 1)

	p <- ggplot(data=results, aes(x=varNice, y=auc, col=region, pch=numHomeRanges)) +
		geom_point(size=2) +
		labs(shape='Number of Home\nRanges as Covariate', color='Region Covariate') +
		xlab(NULL) + ylab('AUC') +
		# ylim(ylim[1], ylim[2]) +
		theme(axis.text.y=element_text(size=12)) +
		coord_flip()
	
	pdf('./Figures & Tables/Occupancy - Univariate/Occupancy - Univariate Binary Models Cross-validation Results.pdf', width=8, height=8)
		print(p)
	dev.off()
	

say('DONE!!!', level=1, deco='%')
