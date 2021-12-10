### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('E:/Ecology/Drive/Research Active/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/03a New Mexico Pika Occupancy & Abundance Analysis - Occupancy with Simple Models.r')
###
### CONTENTS ###
### setup ###
### ORDINAL simple OCCUPANCY analysis: all-data models ###
### BINARY simple OCCUPANCY analysis: all-data models ###
### ORDINAL simple OCCUPANCY analysis: cross-validated models ###
### BINARY simple OCCUPANCY analysis: cross-validated models ###
### ORDINAL simple OCCUPANCY analysis: cross-validated models summary ###
### BINARY simple OCCUPANCY analysis: cross-validated models summary ###

#############
### setup ###
#############

	# drive <- 'C:'
	# drive <- 'D:'
	drive <- 'E:'

	source(paste0(drive, '/Ecology/Drive/Research Active/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r'))

# say('##########################################################')
# say('### ORDINAL simple OCCUPANCY analysis: all-data models ###')
# say('##########################################################')

	# load('./Data/03 New Mexico Pika - Assigned Folds.rda')
	# pika$latestOccStatus <- factor(pika$latestOccStatus, levels=c('0 never', '1 old', '2 occupied'), ordered=TRUE)
	# pika$region <- as.factor(pika$region)
	
	# pika$numHomeRangesScaled <- scale(pika$numHomeRanges)
	
	# vars <- getVars('occupancy')
	# pika[ , vars] <- scale(pika[ , vars])

	# # models
	# formulae <- getFormulae('occupancy')

	# ### using all sites (no cross-validation)
	# #########################################
	
		# ### null models
		# ###############

			# model1 <- polr(latestOccStatus ~ 1, data=pika, Hess=TRUE)
			# model2 <- polr(latestOccStatus ~ numHomeRangesScaled, data=pika, Hess=TRUE)
			# model3 <- polr(latestOccStatus ~ region, data=pika, Hess=TRUE)
			# model4 <- polr(latestOccStatus ~ numHomeRangesScaled + region, data=pika, Hess=TRUE)

			# aicc1 <- AICc(model1)
			# aicc2 <- AICc(model2)
			# aicc3 <- AICc(model3)
			# aicc4 <- AICc(model4)

			# coeffs <- c(NA, NA, NA, NA)
			# numHomeRanges <- c(FALSE, TRUE, FALSE, TRUE)
			# region <- c(FALSE, FALSE, TRUE, TRUE)
			
			# aicc <- c(aicc1, aicc2, aicc3, aicc4)

			# like1 <- logLik(model1)
			# like2 <- logLik(model2)
			# like3 <- logLik(model3)
			# like4 <- logLik(model4)
			
			# likeNull <- like1
			
			# n <- nrow(pika)
			# pseudoR2_1 <- nagelR2(likeNull, like1, n)
			# pseudoR2_2 <- nagelR2(likeNull, like2, n)
			# pseudoR2_3 <- nagelR2(likeNull, like3, n)
			# pseudoR2_4 <- nagelR2(likeNull, like4, n)
			
			# pseudoR2 <- c(pseudoR2_1, pseudoR2_2, pseudoR2_3, pseudoR2_4)

			# results <- data.frame(
				# model = '(Intercept)',
				# term1 = NA,
				# term2 = NA,
				# term3 = NA,
				# term4 = NA,
				# numHomeRanges = numHomeRanges,
				# region = region,
				# aicc = aicc,
				# pseudoR2 = pseudoR2
			# )
			

		# ### by climate variable
		# #######################
		
		# for (formula in formulae) {
			
			# say(formula)
			
			# form1 <- as.formula(paste0('latestOccStatus ~ 1 + ', formula))
			# form2 <- as.formula(paste0('latestOccStatus ~ 1 + ', formula, ' + numHomeRangesScaled'))
			# form3 <- as.formula(paste0('latestOccStatus ~ 1 + ', formula, ' + region'))
			# form4 <- as.formula(paste0('latestOccStatus ~ 1 + ', formula, ' + numHomeRangesScaled + region'))
			
			# model1 <- polr(form1, data=pika, Hess=TRUE)
			# model2 <- polr(form2, data=pika, Hess=TRUE)
			# model3 <- polr(form3, data=pika, Hess=TRUE)
			# model4 <- polr(form4, data=pika, Hess=TRUE)
			
			# aicc1 <- AICc(model1)
			# aicc2 <- AICc(model2)
			# aicc3 <- AICc(model3)
			# aicc4 <- AICc(model4)
	
			# terms <- extractTerms(model1, model2, model3, model4)
			# term1 <- terms$term1
			# term2 <- terms$term2
			# term3 <- terms$term3
			# term4 <- terms$term4
			
			# numHomeRanges <- c(FALSE, TRUE, FALSE, TRUE)
			# region <- c(FALSE, FALSE, TRUE, TRUE)
			
			# aicc <- c(aicc1, aicc2, aicc3, aicc4)
			
			# like1 <- logLik(model1)
			# like2 <- logLik(model2)
			# like3 <- logLik(model3)
			# like4 <- logLik(model4)
			
			# pseudoR2_1 <- nagelR2(likeNull, like1, n)
			# pseudoR2_2 <- nagelR2(likeNull, like2, n)
			# pseudoR2_3 <- nagelR2(likeNull, like3, n)
			# pseudoR2_4 <- nagelR2(likeNull, like4, n)
			
			# pseudoR2 <- c(pseudoR2_1, pseudoR2_2, pseudoR2_3, pseudoR2_4)
			
			# results <- rbind(
				# results,
				# data.frame(
					# model = formula,
					# term1 = term1,
					# term2 = term2,
					# term3 = term3,
					# term4 = term4,
					# numHomeRanges = numHomeRanges,
					# region = region,
					# aicc = aicc,
					# pseudoR2 = pseudoR2
				# )
			# )

		# } # next variable

		# ### reports
		# ###########
		
			# for (occWindow in c(occWindows_y, NA)) {
			
				# if (is.na(occWindow)) {
					# thisResults <- results
					# nice <- paste0(paste(occWindows_y, collapse=' & '), '-yr Windows')
				# } else {

					# thisResults <- rbind(
						# results[grepl(results$model, pattern=paste0(occWindow, 'yrWindow')), ],
						# results[results$model == '(Intercept)', ]
					# )
						
					# nice <- paste0(occWindow, '-yr Window')
				# }
				
				# thisResults$deltaAicc <- thisResults$aicc - min(thisResults$aicc)
				# w <- exp(-0.5 * thisResults$deltaAicc)
				# thisResults$weight <- w / sum(w)

				# thisResults <- thisResults[order(thisResults$weight, decreasing=TRUE), ]
				# rownames(thisResults) <- NULL

				# file <- paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Ordinal Models Using All Data - ', nice, '.csv')
				# write.csv(thisResults, file, row.names=FALSE)
				
			# } # next window

# say('#########################################################')
# say('### BINARY simple OCCUPANCY analysis: all-data models ###')
# say('#########################################################')

	# load('./Data/03 New Mexico Pika - Assigned Folds.rda')
	# pika$region <- as.factor(pika$region)
	
	# pika$numHomeRangesScaled <- scale(pika$numHomeRanges)

	# vars <- getVars('occupancy')
	# pika[ , vars] <- scale(pika[ , vars])

	# # models
	# formulae <- getFormulae('occupancy')

	# ### null models
	# ###############

		# model1 <- glm(presAbs ~ 1, data=pika, family=binomial)
		# model2 <- glm(presAbs ~ numHomeRangesScaled, data=pika, family=binomial)
		# model3 <- glm(presAbs ~ region, data=pika, family=binomial)
		# model4 <- glm(presAbs ~ numHomeRangesScaled + region, data=pika, family=binomial)

		# aicc1 <- AICc(model1)
		# aicc2 <- AICc(model2)
		# aicc3 <- AICc(model3)
		# aicc4 <- AICc(model4)

		# numHomeRanges <- c(FALSE, TRUE, FALSE, TRUE)
		# region <- c(FALSE, FALSE, TRUE, TRUE)
		
		# aicc <- c(aicc1, aicc2, aicc3, aicc4)
		
		# like1 <- logLik(model1)
		# like2 <- logLik(model2)
		# like3 <- logLik(model3)
		# like4 <- logLik(model4)
		
		# likeNull <- like1
		
		# n <- nrow(pika)
		# pseudoR2_1 <- nagelR2(likeNull, like1, n)
		# pseudoR2_2 <- nagelR2(likeNull, like2, n)
		# pseudoR2_3 <- nagelR2(likeNull, like3, n)
		# pseudoR2_4 <- nagelR2(likeNull, like4, n)
		
		# pseudoR2 <- c(pseudoR2_1, pseudoR2_2, pseudoR2_3, pseudoR2_4)

		# results <- data.frame(
			# model = '(Intercept)',
			# term1 = NA,
			# term2 = NA,
			# term3 = NA,
			# term4 = NA,
			# numHomeRanges = numHomeRanges,
			# region = region,
			# aicc = aicc,
			# pseudoR2 = pseudoR2
		# )

	# ### by climate variable
	# #######################

		# for (formula in formulae) {
			
			# say(formula)
			
			# form1 <- as.formula(paste0('presAbs ~ 1 + ', formula))
			# form2 <- as.formula(paste0('presAbs ~ 1 + ', formula, ' + numHomeRangesScaled'))
			# form3 <- as.formula(paste0('presAbs ~ 1 + ', formula, ' + region'))
			# form4 <- as.formula(paste0('presAbs ~ 1 + ', formula, ' + numHomeRangesScaled + region'))

			# model1 <- glm(form1, data=pika, family=binomial)
			# model2 <- glm(form2, data=pika, family=binomial)
			# model3 <- glm(form3, data=pika, family=binomial)
			# model4 <- glm(form4, data=pika, family=binomial)
			
			# aicc1 <- AICc(model1)
			# aicc2 <- AICc(model2)
			# aicc3 <- AICc(model3)
			# aicc4 <- AICc(model4)
			
			# terms <- extractTerms(model1, model2, model3, model4)
			# term1 <- terms$term1
			# term2 <- terms$term2
			# term3 <- terms$term3
			# term4 <- terms$term4
			
			# numHomeRanges <- c(FALSE, TRUE, FALSE, TRUE)
			# region <- c(FALSE, FALSE, TRUE, TRUE)
			
			# aicc <- c(aicc1, aicc2, aicc3, aicc4)
			
			# like1 <- logLik(model1)
			# like2 <- logLik(model2)
			# like3 <- logLik(model3)
			# like4 <- logLik(model4)
			
			# pseudoR2_1 <- nagelR2(likeNull, like1, n)
			# pseudoR2_2 <- nagelR2(likeNull, like2, n)
			# pseudoR2_3 <- nagelR2(likeNull, like3, n)
			# pseudoR2_4 <- nagelR2(likeNull, like4, n)
			
			# pseudoR2 <- c(pseudoR2_1, pseudoR2_2, pseudoR2_3, pseudoR2_4)

			# results <- rbind(
				# results,
				# data.frame(
					# model = formula,
					# term1 = term1,
					# term2 = term2,
					# term3 = term3,
					# term4 = term4,
					# numHomeRanges = numHomeRanges,
					# region = region,
					# aicc = aicc,
					# pseudoR2 = pseudoR2
				# )
			# )
			
		# } # next variable

	# ### reports
	# ###########
	
		# for (occWindow in c(occWindows_y, NA)) {
		
			# if (is.na(occWindow)) {
				# thisResults <- results
				# nice <- paste0(paste(occWindows_y, collapse=' & '), '-yr Windows')
			# } else {
				# thisResults <- rbind(
					# results[grepl(results$model, pattern=paste0(occWindow, 'yrWindow')), ],
					# results[results$model == '(Intercept)', ]
				# )
				# nice <- paste0(occWindow, '-yr Window')
			# }
			
			# thisResults$deltaAicc <- thisResults$aicc - min(thisResults$aicc)
			# w <- exp(-0.5 * thisResults$deltaAicc)
			# thisResults$weight <- w / sum(w)

			# thisResults <- thisResults[order(thisResults$weight, decreasing=TRUE), ]
			# rownames(thisResults) <- NULL

			# file <- paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Binary Models Using All Data - ', nice, '.csv')
			# write.csv(thisResults, file, row.names=FALSE)
			
		# } # next window

# say('#################################################################')
# say('### ORDINAL simple OCCUPANCY analysis: cross-validated models ###')
# say('#################################################################')

	# load('./Data/03 New Mexico Pika - Assigned Folds.rda')
	# pika$latestOccStatus <- factor(pika$latestOccStatus, levels=c('0 never', '1 old', '2 occupied'), ordered=TRUE)
	# pika$region <- as.factor(pika$region)
	# pika$numHomeRangesScaled <- scale(pika$numHomeRanges)

	# vars <- getVars('occupancy')
	# pika[ , vars] <- scale(pika[ , vars])

	# # models
	# formulae <- getFormulae('occupancy')

	# folds <- expand.grid(nwFold=1:2, swFold=1:2, neFold=1:2, seFold=1:2)

	# results <- data.frame()

	# ### intercept-only models
	# #########################
	
		# for (k in 1:nrow(folds)) {

			# trainIndex <- which(
				# (pika$region == 'northwest' & pika$fold == folds$nwFold[k]) |
				# (pika$region == 'southwest' & pika$fold == folds$swFold[k]) |
				# (pika$region == 'northeast' & pika$fold == folds$neFold[k]) |
				# (pika$region == 'southeast' & pika$fold == folds$seFold[k])
			# )
			
			# testIndex <- which(
				# (pika$region == 'northwest' & pika$fold != folds$nwFold[k]) |
				# (pika$region == 'southwest' & pika$fold != folds$swFold[k]) |
				# (pika$region == 'northeast' & pika$fold != folds$neFold[k]) |
				# (pika$region == 'southeast' & pika$fold != folds$seFold[k])
			# )
			
			# trainData <- pika[trainIndex, ]
			# testData <- pika[testIndex, ]

			# # train models

			# model1 <- polr(latestOccStatus ~ 1, data=trainData, Hess=TRUE)
			# model2 <- polr(latestOccStatus ~ 1 + numHomeRangesScaled, data=trainData, Hess=TRUE)
			# model3 <- polr(latestOccStatus ~ 1 + region, data=trainData, Hess=TRUE)
			# model4 <- polr(latestOccStatus ~ 1 + numHomeRangesScaled + region, data=trainData, Hess=TRUE)
			
			# # test models

			# pred1 <- predict(model1, testData, type='probs')
			# pred2 <- predict(model2, testData, type='probs')
			# pred3 <- predict(model3, testData, type='probs')
			# pred4 <- predict(model4, testData, type='probs')
			
			# pred1_never <- pred1[testData$latestOccStatus == '0 never', 1]
			# pred1_old <- pred1[testData$latestOccStatus == '1 old', 2]
			# pred1_occ <- pred1[testData$latestOccStatus == '2 occupied', 3]
			# auc1 <- aucMultiWeighted(pred1_never, pred1_old, pred1_occ)
			
			# pred2_never <- pred2[testData$latestOccStatus == '0 never', 1]
			# pred2_old <- pred2[testData$latestOccStatus == '1 old', 2]
			# pred2_occ <- pred2[testData$latestOccStatus == '2 occupied', 3]
			# auc2 <- aucMultiWeighted(pred2_never, pred2_old, pred2_occ)
			
			# pred3_never <- pred3[testData$latestOccStatus == '0 never', 1]
			# pred3_old <- pred3[testData$latestOccStatus == '1 old', 2]
			# pred3_occ <- pred3[testData$latestOccStatus == '2 occupied', 3]
			# auc3 <- aucMultiWeighted(pred3_never, pred3_old, pred3_occ)
			
			# pred4_never <- pred4[testData$latestOccStatus == '0 never', 1]
			# pred4_old <- pred4[testData$latestOccStatus == '1 old', 2]
			# pred4_occ <- pred4[testData$latestOccStatus == '2 occupied', 3]
			# auc4 <- aucMultiWeighted(pred4_never, pred4_old, pred4_occ)
			
			# aucs <- c(auc1[['multivariate']], auc2[['multivariate']], auc3[['multivariate']], auc4[['multivariate']])
						
			# results <- rbind(
				# results,
				# data.frame(
					# model = '(Intercept)',
					# timeFrame = NA,
					# fold = k,
					# numHomeRanges = c(FALSE, TRUE, FALSE, TRUE),
					# region = c(FALSE, FALSE, TRUE, TRUE),
					# auc = aucs
				# )
			# )
			
		# } # next fold
	
	# ### by climate variable
	# #######################
	
	# for (formula in formulae) {
		
		# say(formula)
		
		# for (k in 1:nrow(folds)) {

			# trainIndex <- which(
				# (pika$region == 'northwest' & pika$fold == folds$nwFold[k]) |
				# (pika$region == 'southwest' & pika$fold == folds$swFold[k]) |
				# (pika$region == 'northeast' & pika$fold == folds$neFold[k]) |
				# (pika$region == 'southeast' & pika$fold == folds$seFold[k])
			# )
			
			# testIndex <- which(
				# (pika$region == 'northwest' & pika$fold != folds$nwFold[k]) |
				# (pika$region == 'southwest' & pika$fold != folds$swFold[k]) |
				# (pika$region == 'northeast' & pika$fold != folds$neFold[k]) |
				# (pika$region == 'southeast' & pika$fold != folds$seFold[k])
			# )
			
			# trainData <- pika[trainIndex, ]
			# testData <- pika[testIndex, ]

			# form1 <- as.formula(paste0('latestOccStatus ~ 1 + ', formula))
			# form2 <- as.formula(paste0('latestOccStatus ~ 1 + ', formula, ' + numHomeRangesScaled'))
			# form3 <- as.formula(paste0('latestOccStatus ~ 1 + ', formula, ' + region'))
			# form4 <- as.formula(paste0('latestOccStatus ~ 1 + ', formula, ' + numHomeRangesScaled + region'))

			# # train models
			# model1 <- polr(form1, data=trainData, Hess=TRUE)
			# model2 <- polr(form2, data=trainData, Hess=TRUE)
			# model3 <- polr(form3, data=trainData, Hess=TRUE)
			# model4 <- polr(form4, data=trainData, Hess=TRUE)
			
			# # test models
			# pred1 <- predict(model1, testData, type='probs')
			# pred2 <- predict(model2, testData, type='probs')
			# pred3 <- predict(model3, testData, type='probs')
			# pred4 <- predict(model4, testData, type='probs')
			
			# pred1_never <- pred1[testData$latestOccStatus == '0 never', 1]
			# pred1_old <- pred1[testData$latestOccStatus == '1 old', 2]
			# pred1_occ <- pred1[testData$latestOccStatus == '2 occupied', 3]
			# auc1 <- aucMultiWeighted(pred1_never, pred1_old, pred1_occ)
			
			# pred2_never <- pred2[testData$latestOccStatus == '0 never', 1]
			# pred2_old <- pred2[testData$latestOccStatus == '1 old', 2]
			# pred2_occ <- pred2[testData$latestOccStatus == '2 occupied', 3]
			# auc2 <- aucMultiWeighted(pred2_never, pred2_old, pred2_occ)
			
			# pred3_never <- pred3[testData$latestOccStatus == '0 never', 1]
			# pred3_old <- pred3[testData$latestOccStatus == '1 old', 2]
			# pred3_occ <- pred3[testData$latestOccStatus == '2 occupied', 3]
			# auc3 <- aucMultiWeighted(pred3_never, pred3_old, pred3_occ)
			
			# pred4_never <- pred4[testData$latestOccStatus == '0 never', 1]
			# pred4_old <- pred4[testData$latestOccStatus == '1 old', 2]
			# pred4_occ <- pred4[testData$latestOccStatus == '2 occupied', 3]
			# auc4 <- aucMultiWeighted(pred4_never, pred4_old, pred4_occ)
			
			# aucs <- c(auc1[['multivariate']], auc2[['multivariate']], auc3[['multivariate']], auc4[['multivariate']])
			
			# timeFrame <- substr(formula, nchar(formula) - 9, nchar(formula) - 8)
			# if (substr(timeFrame, 1, 1) == '_') timeFrame <- substr(timeFrame, 2, 2)
			# timeFrame <- as.integer(timeFrame)
			
			# results <- rbind(
				# results,
				# data.frame(
					# model = formula,
					# timeFrame = timeFrame,
					# fold = k,
					# numHomeRanges = c(FALSE, TRUE, FALSE, TRUE),
					# region = c(FALSE, FALSE, TRUE, TRUE),
					# auc = c(auc1[['multivariate']], auc2[['multivariate']], auc3[['multivariate']], auc4[['multivariate']])
				# )
			# )
			
		# } # next fold
			
	# } # next variable
	
	# results <- results[order(results$auc, decreasing=TRUE), ]
	
	# write.csv(results, file='./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Ordinal Models Cross-validation Results.csv', row.names=FALSE)

# say('################################################################')
# say('### BINARY simple OCCUPANCY analysis: cross-validated models ###')
# say('################################################################')

	# load('./Data/03 New Mexico Pika - Assigned Folds.rda')
	# pika$region <- as.factor(pika$region)
	
	# pika$numHomeRangesScaled <- scale(pika$numHomeRanges)

	# vars <- getVars('occupancy')
	# pika[ , vars] <- scale(pika[ , vars])

	# # models
	# formulae <- getFormulae('occupancy')

	# folds <- expand.grid(nwFold=1:2, swFold=1:2, neFold=1:2, seFold=1:2)

	# results <- data.frame()
	
	# ### intercept-only models
	# #########################
	
		# for (k in 1:nrow(folds)) {

			# trainIndex <- which(
				# (pika$region == 'northwest' & pika$fold == folds$nwFold[k]) |
				# (pika$region == 'southwest' & pika$fold == folds$swFold[k]) |
				# (pika$region == 'northeast' & pika$fold == folds$neFold[k]) |
				# (pika$region == 'southeast' & pika$fold == folds$seFold[k])
			# )
			
			# testIndex <- which(
				# (pika$region == 'northwest' & pika$fold != folds$nwFold[k]) |
				# (pika$region == 'southwest' & pika$fold != folds$swFold[k]) |
				# (pika$region == 'northeast' & pika$fold != folds$neFold[k]) |
				# (pika$region == 'southeast' & pika$fold != folds$seFold[k])
			# )
			
			# trainData <- pika[trainIndex, ]
			# testData <- pika[testIndex, ]

			# # train models

			# model1 <- glm(presAbs ~ 1, data=trainData, family=binomial)
			# model2 <- glm(presAbs ~ 1 + numHomeRangesScaled, data=trainData, family=binomial)
			# model3 <- glm(presAbs ~ 1 + region, data=trainData, family=binomial)
			# model4 <- glm(presAbs ~ 1 + numHomeRangesScaled + region, data=trainData, family=binomial)
			
			# # test models

			# pred1 <- predict(model1, testData, type='response')
			# pred2 <- predict(model2, testData, type='response')
			# pred3 <- predict(model3, testData, type='response')
			# pred4 <- predict(model4, testData, type='response')
			
			# predPres1 <- pred1[testData$presAbs == 1]
			# predPres2 <- pred2[testData$presAbs == 1]
			# predPres3 <- pred3[testData$presAbs == 1]
			# predPres4 <- pred4[testData$presAbs == 1]
			
			# predAbs1 <- pred1[testData$presAbs == 0]
			# predAbs2 <- pred2[testData$presAbs == 0]
			# predAbs3 <- pred3[testData$presAbs == 0]
			# predAbs4 <- pred4[testData$presAbs == 0]
			
			# auc1 <- aucWeighted(predPres1, predAbs1)
			# auc2 <- aucWeighted(predPres2, predAbs2)
			# auc3 <- aucWeighted(predPres3, predAbs3)
			# auc4 <- aucWeighted(predPres4, predAbs4)
	
			# aucs <- c(auc1, auc2, auc3, auc4)
						
			# results <- rbind(
				# results,
				# data.frame(
					# model = '(Intercept)',
					# timeFrame = NA,
					# fold = k,
					# numHomeRanges = c(FALSE, TRUE, FALSE, TRUE),
					# region = c(FALSE, FALSE, TRUE, TRUE),
					# auc = aucs
				# )
			# )
			
		# } # next fold

	# ### by climate variable
	# #######################

	# for (formula in formulae) {
		
		# say(formula)

		# for (k in 1:nrow(folds)) {

			# trainIndex <- which(
				# (pika$region == 'northwest' & pika$fold == folds$nwFold[k]) |
				# (pika$region == 'southwest' & pika$fold == folds$swFold[k]) |
				# (pika$region == 'northeast' & pika$fold == folds$neFold[k]) |
				# (pika$region == 'southeast' & pika$fold == folds$seFold[k])
			# )
			
			# testIndex <- which(
				# (pika$region == 'northwest' & pika$fold != folds$nwFold[k]) |
				# (pika$region == 'southwest' & pika$fold != folds$swFold[k]) |
				# (pika$region == 'northeast' & pika$fold != folds$neFold[k]) |
				# (pika$region == 'southeast' & pika$fold != folds$seFold[k])
			# )
			
			# trainData <- pika[trainIndex, ]
			# testData <- pika[testIndex, ]

			# # train models
			# form1 <- as.formula(paste0('presAbs ~ 1 + ', formula))
			# form2 <- as.formula(paste0('presAbs ~ 1 + ', formula, ' + numHomeRangesScaled'))
			# form3 <- as.formula(paste0('presAbs ~ 1 + ', formula, ' + region'))
			# form4 <- as.formula(paste0('presAbs ~ 1 + ', formula, ' + numHomeRangesScaled + region'))

			# model1 <- glm(form1, data=trainData, family=binomial)
			# model2 <- glm(form2, data=trainData, family=binomial)
			# model3 <- glm(form3, data=trainData, family=binomial)
			# model4 <- glm(form4, data=trainData, family=binomial)
			
			# # test models
			# pred1 <- predict(model1, testData, type='response')
			# pred2 <- predict(model2, testData, type='response')
			# pred3 <- predict(model3, testData, type='response')
			# pred4 <- predict(model4, testData, type='response')
			
			# predPres1 <- pred1[testData$presAbs == 1]
			# predPres2 <- pred2[testData$presAbs == 1]
			# predPres3 <- pred3[testData$presAbs == 1]
			# predPres4 <- pred4[testData$presAbs == 1]
			
			# predAbs1 <- pred1[testData$presAbs == 0]
			# predAbs2 <- pred2[testData$presAbs == 0]
			# predAbs3 <- pred3[testData$presAbs == 0]
			# predAbs4 <- pred4[testData$presAbs == 0]
			
			# auc1 <- aucWeighted(predPres1, predAbs1)
			# auc2 <- aucWeighted(predPres2, predAbs2)
			# auc3 <- aucWeighted(predPres3, predAbs3)
			# auc4 <- aucWeighted(predPres4, predAbs4)
			
			# aucs <- c(auc1, auc2, auc3, auc4)
			
			# timeFrame <- substr(formula, nchar(formula) - 9, nchar(formula) - 8)
			# if (substr(timeFrame, 1, 1) == '_') timeFrame <- substr(timeFrame, 2, 2)
			# timeFrame <- as.integer(timeFrame)
			
			# results <- rbind(
				# results,
				# data.frame(
					# model = formula,
					# timeFrame = timeFrame,
					# fold = k,
					# numHomeRanges = c(FALSE, TRUE, FALSE, TRUE),
					# region = c(FALSE, FALSE, TRUE, TRUE),
					# auc = aucs
				# )
			# )
			
		# } # next fold
			
	# } # next variable
	
	# results <- results[order(results$auc, decreasing=TRUE), ]
	
	# write.csv(results, file='./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Binary Models Cross-validation Results.csv', row.names=FALSE)

# say('#########################################################################')
# say('### ORDINAL simple OCCUPANCY analysis: cross-validated models summary ###')
# say('#########################################################################')
	
	# results <- read.csv('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Ordinal Models Cross-validation Results.csv')

	# # intercept-only models
	# nulls <- results[results$model == '(Intercept)', ]
	# nulls <- aggregate(nulls, by=list(nulls$numHomeRanges, nulls$region), FUN=mean)
	# nulls$model <- nulls$timeFrame <- nulls$fold <- nulls$numHomeRanges <- nulls$region <- NULL
	# names(nulls)[1:2] <- c('numHomeRanges', 'region')
	# n <- nrow(nulls)

	# first <- data.frame(
		# model = rep('(Intercept)', n),
		# timeFrame = rep(NA, n)
	# )
	
	# last <- data.frame(modelNice = rep('(Intercept)', n))

	# nulls <- cbind(first, nulls, last)

	# # climate models
	# clims <- results[results$model != '(Intercept)', ]
	# clims <- aggregate(clims, by=list(clims$model, clims$timeFrame, clims$numHomeRanges, clims$region), FUN=mean)
	# clims$model <- clims$fold <- clims$timeFrame <- clims$numHomeRanges <- clims$region <- NULL
	# names(clims)[1:4] <- c('model', 'timeFrame', 'numHomeRanges', 'region')
	# clims$modelNice <- niceFormulae(clims$model, occOrDens='occupancy')
	
	# results <- rbind(nulls, clims)
	# results <- results[order(results$auc, decreasing=TRUE), ]
	# results$modelNice <- factor(results$modelNice, levels=unique(results$modelNice))
	# # results$modelNice <- factor(results$modelNice)

	# ylim <- c(min(0.5, min(results$auc)), 1)

	# p <- ggplot(data=results, aes(x=modelNice, y=auc, col=region, pch=numHomeRanges)) +
		# geom_point(size=2) +
		# labs(shape='Number of Home\nRanges as Covariate', color='Region Covariate') +
		# xlab(NULL) + ylab('Multivariate AUC') +
		# # ylim(ylim[1], ylim[2]) +
		# theme(
			# axis.text.y=element_text(size=6),
			# legend.position=c(0.85, 0.9)
		# ) +
		# coord_flip()
	
	# pdf('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Ordinal Models Cross-validation Results.pdf', width=8, height=10.5)
		# print(p)
	# dev.off()

# say('########################################################################')
# say('### BINARY simple OCCUPANCY analysis: cross-validated models summary ###')
# say('########################################################################')
	
	# results <- read.csv('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Binary Models Cross-validation Results.csv')
	
	# # intercept-only models
	# nulls <- results[results$model == '(Intercept)', ]
	# nulls <- aggregate(nulls, by=list(nulls$numHomeRanges, nulls$region), FUN=mean)
	# nulls$model <- nulls$timeFrame <- nulls$fold <- nulls$numHomeRanges <- nulls$region <- NULL
	# names(nulls)[1:2] <- c('numHomeRanges', 'region')
	# n <- nrow(nulls)

	# first <- data.frame(
		# model = rep('(Intercept)', n),
		# timeFrame = rep(NA, n)
	# )
	
	# last <- data.frame(modelNice = rep('(Intercept)', n))

	# nulls <- cbind(first, nulls, last)

	# # climate models
	# clims <- results[results$model != '(Intercept)', ]
	# clims <- aggregate(clims, by=list(clims$model, clims$timeFrame, clims$numHomeRanges, clims$region), FUN=mean)
	# clims$model <- clims$fold <- clims$timeFrame <- clims$numHomeRanges <- clims$region <- NULL
	# names(clims)[1:4] <- c('model', 'timeFrame', 'numHomeRanges', 'region')
	# clims$modelNice <- niceFormulae(clims$model, occOrDens='occupancy')
	
	# results <- rbind(nulls, clims)
	# results <- results[order(results$auc, decreasing=TRUE), ]
	# results$modelNice <- factor(results$modelNice, levels=unique(results$modelNice))

	# ylim <- c(min(0.5, min(results$auc)), 1)

	# p <- ggplot(data=results, aes(x=modelNice, y=auc, col=region, pch=numHomeRanges)) +
		# geom_point(size=2) +
		# labs(shape='Number of Home\nRanges as Covariate', color='Region Covariate') +
		# xlab(NULL) + ylab('AUC') +
		# # ylim(ylim[1], ylim[2]) +
		# theme(
			# axis.text.y=element_text(size=8),
			# legend.title=element_text(size=6),
			# legend.text=element_text(size=6),
			# legend.position=c(-1.1, 0.9),
			# plot.margin = unit(c(0.5, 0.5, 0.5, 2), 'cm')
		# ) +
		# coord_flip()
	
	# pdf('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Binary Models Cross-validation Results.pdf', width=8, height=10.5)
		# print(p)
	# dev.off()
	

say('DONE!!!', level=1, deco='%')
