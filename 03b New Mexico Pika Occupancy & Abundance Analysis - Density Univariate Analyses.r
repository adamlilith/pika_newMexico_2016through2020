### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/03b New Mexico Pika Occupancy & Abundance Analysis - Density Univariate Analyses.r')
###
### CONTENTS ###
### setup ###
### univariate DENSITY models: all-data models ###
### univariate DENSITY analysis: cross-validated models ###
### univariate DENSITY analysis: cross-validated models summary ###

#############
### setup ###
#############

	source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r')

say('##################################################')
say('### univariate DENSITY models: all-data models ###')
say('##################################################')

	load('./Data/03 New Mexico Pika - Assigned Folds.rda')
	
	pika <- pika[!is.na(pika$latestDensity), ]
	pika$region <- as.factor(pika$region)
	vars <- names(pika)[grepl(names(pika), pattern='densVar_')]

	# predictors
	
	### intercept-only models
	#########################

		model1 <- glm(latestDensity ~ 1, data=pika, family=Gamma(link='log'))
		model2 <- glm(latestDensity ~ 1 + region, data=pika, family=Gamma(link='log'))
		
		coeffs <- c(NA, NA)
		
		region <- c(FALSE, TRUE)
		
		aicc1 <- AICc(model1)
		aicc2 <- AICc(model2)
		aicc <- c(aicc1, aicc2)
		
		ll1 <- logLik(model1)
		ll2 <- logLik(model2)
		llNull <- ll1
		
		pseudoR2_1 <- nagelR2(llNull, ll1, n=nrow(pika))
		pseudoR2_2 <- nagelR2(llNull, ll2, n=nrow(pika))
		pseudoR2 <- c(pseudoR2_1, pseudoR2_2)
		
		results <- data.frame(
			var = '(Intercept)',
			coeff = coeffs,
			region = region,
			aicc = aicc,
			pseudoR2 = pseudoR2
		)

	### climate models
	##################
	
		for (var in vars) {
			
if (var != 'densVar_subLethalHeat_d_0yrPrior') {
			
			x <- scale(pika[ , var])
			
			model1 <- glm(latestDensity ~ x, data=pika, family=Gamma(link='log'))
			model2 <- glm(latestDensity ~ x + region, data=pika, family=Gamma(link='log'))
			
			coeffs <- c(
				coefficients(model1)[['x']],
				coefficients(model2)[['x']]
			)
			
			region <- c(FALSE, TRUE)
			
			aicc1 <- AICc(model1)
			aicc2 <- AICc(model2)
			aicc <- c(aicc1, aicc2)
			
			ll1 <- logLik(model1)
			ll2 <- logLik(model2)
			
			pseudoR2_1 <- nagelR2(llNull, ll1, n=nrow(pika))
			pseudoR2_2 <- nagelR2(llNull, ll2, n=nrow(pika))
			pseudoR2 <- c(pseudoR2_1, pseudoR2_2)
			
			results <- rbind(
				results,
				data.frame(
					var = var,
					coeff = coeffs,
					region = region,
					aicc = aicc,
					pseudoR2 = pseudoR2
				)
			)
			
			
}
			
		} # next variable

		### reports
		###########
		
			for (densWindow in c(densWindows_y, NA)) {

				if (is.na(densWindow)) {
					densWindow <- results
					nice <- 'Both Years Prior'
				} else {
			
					thisResults <- rbind(
						results[grepl(results$var, pattern=paste0(densWindow, 'yrPrior')), ],
						results[results$var=='(Intercept)', ]
					)
					nice <- paste0(densWindow, ' yr Prior')
					
				}
				
				thisResults$deltaAicc <- thisResults$aicc - min(thisResults$aicc)
				w <- exp(-0.5 * thisResults$deltaAicc)
				thisResults$weight <- w / sum(w)

				thisResults <- thisResults[order(thisResults$weight, decreasing=TRUE), ]
				rownames(thisResults) <- NULL

				file <- paste0('./Figures & Tables/Density - Univariate/Density - Univariate GLMs Using All Data - ', nice, '.csv')
				write.csv(thisResults, file, row.names=FALSE)
			
			} # next window
				
		### model diagnostics
		#####################
	
say('###########################################################')
say('### univariate DENSITY analysis: cross-validated models ###')
say('###########################################################')

	load('./Data/03 New Mexico Pika - Assigned Folds.rda')
	pika <- pika[!is.na(pika$latestDensity), ]
	
	vars <- names(pika)[grepl(names(pika), pattern='densVar_')]

	pika$region <- as.factor(pika$region)

	folds <- expand.grid(nwFold=1:2, swFold=1:2, neFold=1:2, seFold=1:2)
	kFoldResults <- data.frame()

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

			model1 <- glm(latestDensity ~ 1, data=trainData, family=Gamma(link='log'))
			model2 <- glm(latestDensity ~ 1 + region, data=trainData, family=Gamma(link='log'))
			
			# test models

			pred1 <- predict(model1, testData, type='response')
			pred2 <- predict(model2, testData, type='response')
			
			rmse1 <- sqrt(mean((testData$latestDensity - pred1)^2))
			rmse2 <- sqrt(mean((testData$latestDensity - pred2)^2))
			
			meanAbsPercError1 <- mean(abs(pred1 - testData$latestDensity) / testData$latestDensity)
			meanAbsPercError2 <- mean(abs(pred2 - testData$latestDensity) / testData$latestDensity)
			meanAbsPercError <- c(meanAbsPercError1, meanAbsPercError2)
			
			kFoldResults <- rbind(
				kFoldResults,
				data.frame(
					var = '(Intercept)',
					timeFrame = NA,
					fold = k,
					region = c(FALSE, TRUE),
					rmse = c(rmse1, rmse2),
					meanAbsPercError = meanAbsPercError
				)
			)
			
		} # next fold

	### climate models
	##################
	
	for (var in vars) {
		
if (var != 'densVar_subLethalHeat_d_0yrPrior') {
		
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

			model1 <- glm(latestDensity ~ x, data=trainData, family=Gamma(link='log'), start=c(1, 0))
			coeffs <- coefficients(model1)
			
			model2 <- glm(latestDensity ~ x + region, data=trainData, family=Gamma(link='log'), start=c(coeffs, 0, 0, 0))
			
			# test models
			testData$x <- xScaled[testIndex]

			pred1 <- predict(model1, testData, type='response')
			pred2 <- predict(model2, testData, type='response')
			
			rmse1 <- sqrt(mean((testData$latestDensity - pred1)^2))
			rmse2 <- sqrt(mean((testData$latestDensity - pred2)^2))
			
			meanAbsPercError1 <- mean(abs(pred1 - testData$latestDensity) / testData$latestDensity)
			meanAbsPercError2 <- mean(abs(pred2 - testData$latestDensity) / testData$latestDensity)
			meanAbsPercError <- c(meanAbsPercError1, meanAbsPercError2)
			
			timeFrame <- substr(var, nchar(var) - 7, nchar(var) - 7)
			timeFrame <- as.integer(timeFrame)
			
			kFoldResults <- rbind(
				kFoldResults,
				data.frame(
					var = var,
					timeFrame = timeFrame,
					fold = k,
					region = c(FALSE, TRUE),
					rmse = c(rmse1, rmse2),
					meanAbsPercError = meanAbsPercError
				)
			)
			
		} # next fold
			
}			
			
	} # next variable
	
	write.csv(kFoldResults, file='./Figures & Tables/Density - Univariate/Density - Univariate GLMs Cross-validation Results.csv', row.names=FALSE)

say('###################################################################')
say('### univariate DENSITY analysis: cross-validated models summary ###')
say('###################################################################')
	
	results <- read.csv('./Figures & Tables/Density - Univariate/Density - Univariate GLMs Cross-validation Results.csv')
	
	nulls <- results[results$var == '(Intercept)', ]
	nulls <- aggregate(nulls, by=list(nulls$region), FUN=mean)
	nulls$region <- nulls$fold <- NULL
	names(nulls)[1] <- 'region'
	
	nulls <- cbind(nulls, data.frame(varNice = rep('(Intercept)', nrow(nulls))))
	
	clim <- results[results$var != '(Intercept)', ]
	clim <- aggregate(clim, by=list(clim$var, clim$timeFrame, clim$region), FUN=mean)
	clim$var <- clim$fold <- clim$timeFrame <- clim$region <- NULL
	names(clim)[1:3] <- c('var', 'timeFrame', 'region')
	clim$varNice <- makeNiceVars(clim$var, occOrDens='dens')

	results <- rbind(nulls, clim)
	
	results <- results[order(results$meanAbsPercError, decreasing=FALSE), ]
	results$varNice <- factor(results$varNice, levels=rev(unique(results$varNice)))

	p <- ggplot(data=results, aes(x=varNice, y=meanAbsPercError, pch=region)) +
		geom_point(size=2) +
		labs(shape='Region\nCovariate') +
		xlab(NULL) + ylab('Mean absolute percent error') +
		theme(axis.text.y=element_text(size=12)) +
		# scale_y_continuous(trans='log10') +
		scale_y_continuous(labels=scales::percent) +
		coord_flip()
	
	pdf('./Figures & Tables/Density - Univariate/Density - Univariate GLMs Cross-validation Results.pdf', width=8, height=8)
		print(p)
	dev.off()
	
say('DONE!!!', level=1, deco='%')
