### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('E:/Ecology/Drive/Research Active/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/03b New Mexico Pika Occupancy & Abundance Analysis - Density with Simple Models.r')
###
### CONTENTS ###
### setup ###
### univariate DENSITY models: all-data models ###
### univariate DENSITY analysis: cross-validated models ###
### univariate DENSITY analysis: cross-validated models summary ###

#############
### setup ###
#############

	source('E:/Ecology/Drive/Research Active/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r')

say('##################################################')
say('### univariate DENSITY models: all-data models ###')
say('##################################################')

	load('./Data/03 New Mexico Pika - Assigned Folds.rda')
	
	pika <- pika[!is.na(pika$latestDensity), ]
	pika$region <- as.factor(pika$region)
	vars <- getVars('density')
	pika[ , vars] <- scale(pika[ , vars])

	### intercept-only models
	#########################

		model1 <- glm(latestDensity ~ 1, data=pika, family=Gamma(link='log'))
		model2 <- glm(latestDensity ~ 1 + region, data=pika, family=Gamma(link='log'))
		
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
			model = '(Intercept)',
			term1 = NA,
			term2 = NA,
			term3 = NA,
			term4 = NA,
			region = region,
			aicc = aicc,
			pseudoR2 = pseudoR2
		)

	### climate models
	##################
	
		formulae <- getFormulaeDens()

		for (formula in formulae) {
	
			say(formula)

			form1 <- as.formula(paste0('latestDensity ~ 1 + ', formula))
			form2 <- as.formula(paste0('latestDensity ~ 1 + ', formula, ' + region'))
	
			model1 <- glm(form1, data=pika, family=Gamma(link='log'))
			model2 <- glm(form2, data=pika, family=Gamma(link='log'))
			
			terms <- extractTerms(model1, model2)
			term1 <- terms$term1
			term2 <- terms$term2
			term3 <- terms$term3
			term4 <- terms$term4
			
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
					model = formula,
					term1 = term1,
					term2 = term2,
					term3 = term3,
					term4 = term4,
					region = region,
					aicc = aicc,
					pseudoR2 = pseudoR2
				)
			)
			
		} # next formula

		### reports
		###########
		
			for (densWindow in c(densWindows_y, NA)) {

				if (is.na(densWindow)) {
					thisResults <- results
					nice <- 'Both Years Prior'
				} else {
			
					thisResults <- rbind(
						results[grepl(results$model, pattern=paste0(densWindow, 'yrPrior')), ],
						results[results$var=='(Intercept)', ]
					)
					nice <- paste0(densWindow, ' yr Prior')
					
				}
				
				thisResults$deltaAicc <- thisResults$aicc - min(thisResults$aicc)
				w <- exp(-0.5 * thisResults$deltaAicc)
				thisResults$weight <- w / sum(w)

				thisResults <- thisResults[order(thisResults$weight, decreasing=TRUE), ]
				rownames(thisResults) <- NULL

				file <- paste0('./Figures & Tables/Density - Simple Models/Density - Simple GLMs Using All Data - ', nice, '.csv')
				write.csv(thisResults, file, row.names=FALSE)
			
			} # next window
				
		### model diagnostics
		#####################
	
say('###########################################################')
say('### univariate DENSITY analysis: cross-validated models ###')
say('###########################################################')

	load('./Data/03 New Mexico Pika - Assigned Folds.rda')
	
	pika <- pika[!is.na(pika$latestDensity), ]
	pika$region <- as.factor(pika$region)
	vars <- getVars('density')
	pika[ , vars] <- scale(pika[ , vars])

	folds <- expand.grid(nwFold=1:2, swFold=1:2, neFold=1:2, seFold=1:2)
	kFoldResults <- data.frame()

	### intercept-only models
	#########################

		for (k in 1:nrow(folds)) {

			trainIndex <- which(
				(pika$region == 'northwest' & pika$fold == folds$nwFold[k]) |
				(pika$region == 'southwest' & pika$fold == folds$swFold[k]) |
				(pika$region == 'northeast' & pika$fold == folds$neFold[k]) |
				(pika$region == 'southeast' & pika$fold == folds$seFold[k])
			)
			
			testIndex <- which(
				(pika$region == 'northwest' & pika$fold != folds$nwFold[k]) |
				(pika$region == 'southwest' & pika$fold != folds$swFold[k]) |
				(pika$region == 'northeast' & pika$fold != folds$neFold[k]) |
				(pika$region == 'southeast' & pika$fold != folds$seFold[k])
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
					model = '(Intercept)',
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
	
	formulae <- getFormulaeDens()
	for (formula in formulae) {

		say(formula)

		for (k in 1:nrow(folds)) {

			trainIndex <- which(
				(pika$region == 'northwest' & pika$fold == folds$nwFold[k]) |
				(pika$region == 'southwest' & pika$fold == folds$swFold[k]) |
				(pika$region == 'northeast' & pika$fold == folds$neFold[k]) |
				(pika$region == 'southeast' & pika$fold == folds$seFold[k])
			)
			
			testIndex <- which(
				(pika$region == 'northwest' & pika$fold != folds$nwFold[k]) |
				(pika$region == 'southwest' & pika$fold != folds$swFold[k]) |
				(pika$region == 'northeast' & pika$fold != folds$neFold[k]) |
				(pika$region == 'southeast' & pika$fold != folds$seFold[k])
			)
			
			trainData <- pika[trainIndex, ]
			testData <- pika[testIndex, ]

			# train models
			
			form1 <- as.formula(paste0('latestDensity ~ 1 + ', formula))
			form2 <- as.formula(paste0('latestDensity ~ 1 + ', formula, ' + region'))
	
			start <- c(1, rep(0, length(attr(terms(form1), 'term.labels'))))
			model1 <- glm(form1, data=trainData, family=Gamma(link='log'), start=start)
			coeffs <- coefficients(model1)
			start <- c(1, rep(0, length(attr(terms(form2), 'term.labels')) - 1), 0, 0, 0)
			model2 <- glm(form2, data=trainData, family=Gamma(link='log'), start=start)
			
			# test models
			pred1 <- predict(model1, testData, type='response')
			pred2 <- predict(model2, testData, type='response')
			
			rmse1 <- sqrt(mean((testData$latestDensity - pred1)^2))
			rmse2 <- sqrt(mean((testData$latestDensity - pred2)^2))
			
			meanAbsPercError1 <- mean(abs(pred1 - testData$latestDensity) / testData$latestDensity)
			meanAbsPercError2 <- mean(abs(pred2 - testData$latestDensity) / testData$latestDensity)
			meanAbsPercError <- c(meanAbsPercError1, meanAbsPercError2)
			
			for (densWindow_y in densWindows_y) if (grepl(formula, pattern=paste0(densWindow_y, 'yrPrior'))) timeFrame <- densWindow_y
			
			kFoldResults <- rbind(
				kFoldResults,
				data.frame(
					model = formula,
					timeFrame = timeFrame,
					fold = k,
					region = c(FALSE, TRUE),
					rmse = c(rmse1, rmse2),
					meanAbsPercError = meanAbsPercError
				)
			)
			
		} # next fold
			
	} # next model
	
	write.csv(kFoldResults, file='./Figures & Tables/Density - Simple Models/Density - Simple GLMs Cross-validation Results.csv', row.names=FALSE)

say('###################################################################')
say('### univariate DENSITY analysis: cross-validated models summary ###')
say('###################################################################')
	
	results <- read.csv('./Figures & Tables/Density - Simple Models/Density - Simple GLMs Cross-validation Results.csv')
	
	nulls <- results[results$model == '(Intercept)', ]
	nulls <- aggregate(nulls, by=list(nulls$region), FUN=mean)
	nulls$region <- nulls$fold <- NULL
	names(nulls)[1] <- 'region'
	
	nulls <- cbind(nulls, data.frame(modelNice = rep('(Intercept)', nrow(nulls))))
	
	clim <- results[results$model != '(Intercept)', ]
	clim <- aggregate(clim, by=list(clim$model, clim$timeFrame, clim$region), FUN=mean)
	clim$model <- clim$fold <- clim$timeFrame <- clim$region <- NULL
	names(clim)[1:3] <- c('model', 'timeFrame', 'region')
	clim$modelNice <- niceFormulae(clim$model, occOrDens='density')

	results <- rbind(nulls, clim)
	
	results <- results[order(results$meanAbsPercError, decreasing=FALSE), ]
	results$modelNice <- factor(results$modelNice, levels=unique(results$modelNice))

	p <- ggplot(data=results, aes(x=modelNice, y=meanAbsPercError, pch=region)) +
		geom_point(size=2) +
		labs(shape='Region\nCovariate') +
		xlab(NULL) + ylab('Mean absolute percent error') +
		theme(
			axis.text.y=element_text(size=10)
		) +
		# scale_y_continuous(trans='log10') +
		scale_y_continuous(labels=scales::percent) +
		coord_flip()
	
	pdf('./Figures & Tables/Density - Simple Models/Density - Simple GLMs Cross-validation Results.pdf', width=8, height=10.5)
		print(p)
	dev.off()
	
say('DONE!!!', level=1, deco='%')
