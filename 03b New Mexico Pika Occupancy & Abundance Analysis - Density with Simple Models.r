### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/03b New Mexico Pika Occupancy & Abundance Analysis - Density with Simple Models.r')
###
### CONTENTS ###
### setup ###
### univariate DENSITY models: all-data models ###
### univariate DENSITY analysis: cross-validated models ###
### univariate DENSITY analysis: cross-validated models summary ###

### post hoc analysis of DENSITY using all data and distance to nearest patches ###



#############
### setup ###
#############

	source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r')

# say('##################################################')
# say('### univariate DENSITY models: all-data models ###')
# say('##################################################')

	# load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')
	
	# pika <- pika[!is.na(pika$latestDensity), ]
	# pika$region <- as.factor(pika$region)
	# vars <- getVars('density')
	# pika[ , vars] <- scale(pika[ , vars])

	# ### intercept-only models
	# #########################

		# model1 <- glm(latestDensity ~ 1, data=pika, family=Gamma(link='log'))
		# model2 <- glm(latestDensity ~ 1 + region, data=pika, family=Gamma(link='log'))
		
		# region <- c(FALSE, TRUE)
		
		# aicc1 <- AICc(model1)
		# aicc2 <- AICc(model2)
		# aicc <- c(aicc1, aicc2)
		
		# ll1 <- logLik(model1)
		# ll2 <- logLik(model2)
		# llNull <- ll1
		
		# pseudoR2_1 <- nagelR2(llNull, ll1, n=nrow(pika))
		# pseudoR2_2 <- nagelR2(llNull, ll2, n=nrow(pika))
		# pseudoR2 <- c(pseudoR2_1, pseudoR2_2)
		
		# results <- data.frame(
			# model = '(Intercept)',
			# term1 = NA,
			# term2 = NA,
			# term3 = NA,
			# term4 = NA,
			# region = region,
			# aicc = aicc,
			# pseudoR2 = pseudoR2
		# )

	# ### climate models
	# ##################
	
		# formulae <- getFormulaeDens()

		# for (formula in formulae) {
	
			# say(formula)

			# form1 <- as.formula(paste0('latestDensity ~ 1 + ', formula))
			# form2 <- as.formula(paste0('latestDensity ~ 1 + ', formula, ' + region'))
	
			# model1 <- glm(form1, data=pika, family=Gamma(link='log'))
			# model2 <- glm(form2, data=pika, family=Gamma(link='log'))
			
			# terms <- extractTerms(model1, model2)
			# term1 <- terms$term1
			# term2 <- terms$term2
			# term3 <- terms$term3
			# term4 <- terms$term4
			
			# region <- c(FALSE, TRUE)
			
			# aicc1 <- AICc(model1)
			# aicc2 <- AICc(model2)
			# aicc <- c(aicc1, aicc2)
			
			# ll1 <- logLik(model1)
			# ll2 <- logLik(model2)
			
			# pseudoR2_1 <- nagelR2(llNull, ll1, n=nrow(pika))
			# pseudoR2_2 <- nagelR2(llNull, ll2, n=nrow(pika))
			# pseudoR2 <- c(pseudoR2_1, pseudoR2_2)
			
			# results <- rbind(
				# results,
				# data.frame(
					# model = formula,
					# term1 = term1,
					# term2 = term2,
					# term3 = term3,
					# term4 = term4,
					# region = region,
					# aicc = aicc,
					# pseudoR2 = pseudoR2
				# )
			# )
			
		# } # next formula

		# ### reports
		# ###########
		
			# results$deltaAicc <- results$aicc - min(results$aicc)
			# w <- exp(-0.5 * results$deltaAicc)
			# results$weight <- w / sum(w)

			# results <- results[order(results$weight, decreasing=TRUE), ]
			# rownames(results) <- NULL

			# file <- paste0('./Figures & Tables/Density - Simple Models/Density - Simple GLMs.csv')
			# write.csv(results, file, row.names=FALSE)
	
# say('############################################################################################')
# say('### compile table of predictor importance for univariate DENSITY models: all-data models ###')
# say('############################################################################################')
	

	# # rank variables by mean AICc weight

	# # get models
	# file <- paste0('./Figures & Tables/Density - Simple Models/Density - Simple GLMs.csv')
	# models <- read.csv(file)
	
	# # get variables
	# vars <- getVars('density')

	# imp <- data.frame()
	# for (var in vars) {
	
		# index <- which(grepl(models$model, pattern=var))
		# n <- length(index)
		# sumWeight <- sum(models$weight[index])
		# meanWeight <- sumWeight / n
		
		# imp <- rbind(
			# imp,
			# data.frame(
				# variable = var,
				# niceVar = makeNiceVars(var, 'density'),
				# numModels = n,
				# sumWeight = sumWeight,
				# meanWeight = meanWeight
			# )
		# )
		
	# }
	
	# imp <- imp[order(imp$meanWeight, decreasing=TRUE), ]

	# write.csv(imp, paste0('./Figures & Tables/Density - Simple Models/Density - Simple GLMs - Var Import.csv'), row.names=FALSE)
	
# say('###########################################################')
# say('### univariate DENSITY analysis: cross-validated models ###')
# say('###########################################################')

	# load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')
	
	# pika <- pika[!is.na(pika$latestDensity), ]
	# pika$region <- as.factor(pika$region)
	# vars <- getVars('density')
	# pika[ , vars] <- scale(pika[ , vars])

	# folds <- expand.grid(nwFold=1:2, swFold=1:2, neFold=1:2, seFold=1:2)
	# kFoldResults <- data.frame()

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

			# model1 <- glm(latestDensity ~ 1, data=trainData, family=Gamma(link='log'))
			# model2 <- glm(latestDensity ~ 1 + region, data=trainData, family=Gamma(link='log'))
			
			# # test models

			# pred1 <- predict(model1, testData, type='response')
			# pred2 <- predict(model2, testData, type='response')
			
			# rmse1 <- sqrt(mean((testData$latestDensity - pred1)^2))
			# rmse2 <- sqrt(mean((testData$latestDensity - pred2)^2))
			
			# meanAbsPercError1 <- mean(abs(pred1 - testData$latestDensity) / testData$latestDensity)
			# meanAbsPercError2 <- mean(abs(pred2 - testData$latestDensity) / testData$latestDensity)
			# meanAbsPercError <- c(meanAbsPercError1, meanAbsPercError2)
			
			# kFoldResults <- rbind(
				# kFoldResults,
				# data.frame(
					# model = '(Intercept)',
					# fold = k,
					# region = c(FALSE, TRUE),
					# rmse = c(rmse1, rmse2),
					# meanAbsPercError = meanAbsPercError
				# )
			# )
			
		# } # next fold

	# ### climate models
	# ##################
	
	# formulae <- getFormulaeDens()
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
			
			# form1 <- as.formula(paste0('latestDensity ~ 1 + ', formula))
			# form2 <- as.formula(paste0('latestDensity ~ 1 + ', formula, ' + region'))
	
			# start <- c(1, rep(0, length(attr(terms(form1), 'term.labels'))))
			# model1 <- glm(form1, data=trainData, family=Gamma(link='log'), start=start)
			# coeffs <- coefficients(model1)
			# start <- c(1, rep(0, length(attr(terms(form2), 'term.labels')) - 1), 0, 0, 0)
			# model2 <- glm(form2, data=trainData, family=Gamma(link='log'), start=start)
			
			# # test models
			# pred1 <- predict(model1, testData, type='response')
			# pred2 <- predict(model2, testData, type='response')
			
			# rmse1 <- sqrt(mean((testData$latestDensity - pred1)^2))
			# rmse2 <- sqrt(mean((testData$latestDensity - pred2)^2))
			
			# meanAbsPercError1 <- mean(abs(pred1 - testData$latestDensity) / testData$latestDensity)
			# meanAbsPercError2 <- mean(abs(pred2 - testData$latestDensity) / testData$latestDensity)
			# meanAbsPercError <- c(meanAbsPercError1, meanAbsPercError2)
			
			# kFoldResults <- rbind(
				# kFoldResults,
				# data.frame(
					# model = formula,
					# fold = k,
					# region = c(FALSE, TRUE),
					# rmse = c(rmse1, rmse2),
					# meanAbsPercError = meanAbsPercError
				# )
			# )
			
		# } # next fold
			
	# } # next model
	
	# write.csv(kFoldResults, file='./Figures & Tables/Density - Simple Models/Density - Simple GLMs Cross-validation Results.csv', row.names=FALSE)

# say('###################################################################')
# say('### univariate DENSITY analysis: cross-validated models summary ###')
# say('###################################################################')
	
	# results <- read.csv('./Figures & Tables/Density - Simple Models/Density - Simple GLMs Cross-validation Results.csv')
	
	# nulls <- results[results$model == '(Intercept)', ]
	# nulls <- aggregate(nulls, by=list(nulls$region), FUN=mean)
	# nulls$region <- nulls$fold <- NULL
	# names(nulls)[1] <- 'region'
	
	# nulls <- cbind(nulls, data.frame(modelNice = rep('(Intercept)', nrow(nulls))))
	
	# clim <- results[results$model != '(Intercept)', ]
	# clim <- aggregate(clim, by=list(clim$model, clim$region), FUN=mean)
	# clim$model <- clim$fold <- clim$region <- NULL
	# names(clim)[1:2] <- c('model', 'region')
	# clim$modelNice <- niceFormulae(clim$model, occOrDens='density')

	# results <- rbind(nulls, clim)
	
	# results <- results[order(results$meanAbsPercError, decreasing=FALSE), ]
	# results$modelNice <- factor(results$modelNice, levels=unique(results$modelNice))

	# p <- ggplot(data=results, aes(x=modelNice, y=meanAbsPercError, pch=region)) +
		# geom_point(size=2) +
		# labs(shape='Region\nCovariate') +
		# xlab(NULL) + ylab('Mean absolute percent error') +
		# theme(
			# axis.text.y=element_text(size=10)
		# ) +
		# # scale_y_continuous(trans='log10') +
		# scale_y_continuous(labels=scales::percent) +
		# coord_flip()
	
	# pdf('./Figures & Tables/Density - Simple Models/Density - Simple GLMs Cross-validation Results.pdf', width=8, height=10.5)
		# print(p)
	# dev.off()


say('###################################################################################')
say('### post hoc analysis of DENSITY using all data and distance to nearest patches ###')
say('###################################################################################')

	### prepare climate data
	########################
	
		load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')
		pika <- pika[!is.na(pika$latestDensity), ]
		pika$region <- as.factor(pika$region)
		vars <- getVars('density')
		pika[ , vars] <- scale(pika[ , vars])
		
	### generate and evaluate models
	################################
		
		# Assuming no models are intercept-only models!

		# intercept-only null model for pseudo-R2
		modelNull <- glm(latestDensity ~ 1, data=pika, family=Gamma(link='log'))
		logLikeNull <- logLik(modelNull)

		results <- data.frame()

		climModels <- read.csv(paste0('./Figures & Tables/Density - Simple Models/Density - Simple GLMs.csv'))
		
		climModels <- climModels[order(climModels$deltaAicc), ]
		climModels <- climModels[climModels$deltaAicc <= maxDeltaAic_density, , drop=FALSE]
		
		for (countClimModel in 1:nrow(climModels)) {
		
			form <- climModels$model[countClimModel]
			thisForm <- paste('latestDensity ~', form)
			
			region <- climModels$region[countClimModel]
			if (region) thisForm <- paste(thisForm, '+ region')
			
			thisForm_isolation <- paste(thisForm, '+ meanDistToClosestPatchesScaled')
			
			# evaluate n closest patches
			for (thisNumClosestPatches in c(3, 4)) {
			
				thisPika <- pika
				thisPika$meanDistToClosestPatches_m <-
					rowMeans(thisPika[ , paste0('distClosestPatch_patch', 1:thisNumClosestPatches, '_m')])
				thisPika$meanDistToClosestPatchesScaled <- scale(thisPika$meanDistToClosestPatches_m)
			
				model_noIsolation <- glm(thisForm, data=thisPika, family=Gamma(link='log'))
				model_isolation <- glm(thisForm_isolation, data=thisPika, family=Gamma(link='log'))

				aicc_noIsolation <- AICc(model_noIsolation)
				aicc_isolation <- AICc(model_isolation)
				
				logLikeModel_noIsolation <- logLik(model_noIsolation)
				pseudoR2_noIsolation <- nagelR2(logLikeNull, logLikeModel_noIsolation, nrow(thisPika))

				logLikeModel_isolation <- logLik(model_isolation)
				pseudoR2_isolation <- nagelR2(logLikeNull, logLikeModel_isolation, nrow(thisPika))

				isolationCoeff <- coeffs(model_isolation)['meanDistToClosestPatchesScaled']

				results <- rbind(
					results,
					data.frame(
						thisNumClosestPatches = thisNumClosestPatches,
						climateModel = form,
						region = region,
						isolationCoeff = isolationCoeff,
						pseudoR2_noIsolation = pseudoR2_noIsolation,
						pseudoR2_isolation = pseudoR2_isolation,
						aicc_noIsolation = aicc_noIsolation,
						aicc_isolation = aicc_isolation
					)
				
				)
			
			} # next number of closest patches
		
		} # next climate model
		
		results$deltaAic_isolationNoIsolation <- results$aicc_noIsolation - results$aicc_isolation
		results$deltaPseudoR2_isolationNoIsolation <- results$pseudoR2_isolation - results$pseudoR2_noIsolation
		
		rownames(results) <- NULL
		write.csv(results, paste0('./Figures & Tables/Density - Simple Models/Density - Simple GLMs - Isolation.csv'), row.names=FALSE)

	
say('DONE!!!', level=1, deco='%')
