### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/03c New Mexico Pika Occupancy & Abundance Analysis - Occupancy Multivariate Analyses.r')
###
### CONTENTS ###
### setup ###
### ORDINAL multivariate OCCUPANCY analysis: all-data models ###
### ORDINAL multivariate OCCUPANCY analysis: all-data models summary ###
### BINARY multivariate OCCUPANCY analysis: all-data models ###
### BINARY multivariate OCCUPANCY analysis: all-data models summary ###

#############
### setup ###
#############

	source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r')

	### generalization
	alphas <- seq(0, 1, by=0.05)
	# alphas <- seq(0, 1, by=1)
	
say('################################################################')
say('### ORDINAL multivariate OCCUPANCY analysis: all-data models ###')
say('################################################################')

	# generalization
	lambdaMinRatio <- 0.001
	maxiterOut <- 10000

	load('./Data/03 New Mexico Pika - Assigned Folds.rda')
	pika$latestOccStatus <- factor(pika$latestOccStatus, levels=c('0 never', '1 old', '2 occupied'), ordered=TRUE)
	pika$region <- as.factor(pika$region)
	
	pika$numHomeRangesScaled <- scale(pika$numHomeRanges)

	# predictors
	file <- './Predictors/Predictor Variables 2021-04-28 Edits by Adam.xlsx'
	univariatePredTable <- read_excel(file, sheet='Predictor Variables', skip=1)

	y <- pika$latestOccStatus
	
	results <- data.frame()

	### model
	#########

	for (occWindow in c(occWindows_y, NA)) {
	# for (occWindow in c(occWindows_y)) {

		vars <- names(pika)[grepl(names(pika), pattern='occVar_')]
		if (!is.na(occWindow)) vars <- vars[grepl(vars, pattern=paste0('_', occWindow, 'yrWindow'))]

		x <- scale(pika[ , vars])
		colnames(x) <- vars

		## climate only
		
			say('climate model ', occWindow, '-yr window')

			modelNum <- 1
			numHomeRanges <- FALSE
			region <- FALSE

			# implement elastic net for ordinal model
			model <- optOrdinal(x=x, y=y, alphas=alphas, lambdaMinRatio=lambdaMinRatio, maxiterOut=maxiterOut)
			
			this <- which.min(model$aic)
			aic <- model$aic[this]
			coeffVal <- model$coefs[this, vars]
			alpha <- model$args$alpha
			
			results <- rbind(
				results,
				data.frame(
					modelNum=modelNum,
					occWindow=occWindow,
					aic=aic,
					alpha=alpha,
					numHomeRanges=numHomeRanges,
					region=region,
					coeff=vars,
					coeffVal=coeffVal
				)
			)
			
		### climate + number of home ranges
		
			say('climate + HR model ', occWindow, '-yr window')

			modelNum <- 2
			numHomeRanges <- TRUE
			region <- FALSE

			hr <- matrix(pika$numHomeRangesScaled, ncol=1)
			colnames(hr) <- 'numHomeRanges'
			xx <- cbind(x, hr)

			# implement elastic net for ordinal model
			model <- optOrdinal(x=xx, y=y, alphas=alphas, lambdaMinRatio=lambdaMinRatio, maxiterOut=maxiterOut)
			
			this <- which.min(model$aic)
			aic <- model$aic[this]
			coeffVal <- model$coefs[this, vars]
			alpha <- model$args$alpha
			
			results <- rbind(
				results,
				data.frame(
					modelNum=modelNum,
					occWindow=occWindow,
					aic=aic,
					alpha=alpha,
					numHomeRanges=numHomeRanges,
					region=region,
					coeff=vars,
					coeffVal=coeffVal
				)
			)
			
		### climate + region
		
			say('climate + region model ', occWindow, '-yr window')

			modelNum <- 3
			numHomeRanges <- FALSE
			region <- TRUE

			regionMat <- matrix(0, nrow=nrow(x), ncol=4)
			colnames(regionMat) <- levels(pika$region)
			for (i in 1:nrow(regionMat)) {
				regionMat[i, colnames(regionMat) == pika$region[i]] <- 1
			}
			xx <- cbind(x, regionMat)

			# implement elastic net for ordinal model
			model <- optOrdinal(x=xx, y=y, alphas=alphas, lambdaMinRatio=lambdaMinRatio, maxiterOut=maxiterOut)

			this <- which.min(model$aic)
			aic <- model$aic[this]
			coeffVal <- model$coefs[this, vars]
			alpha <- model$args$alpha
			
			results <- rbind(
				results,
				data.frame(
					modelNum=modelNum,
					occWindow=occWindow,
					aic=aic,
					alpha=alpha,
					numHomeRanges=numHomeRanges,
					region=region,
					coeff=vars,
					coeffVal=coeffVal
				)
			)
			
		### climate + home ranges + region

			say('climate + HR + region model ', occWindow, '-yr window')

			modelNum <- 4
			numHomeRanges <- TRUE
			region <- TRUE

			xx <- cbind(x, hr, regionMat)

			# implement elastic net for ordinal model
			model <- optOrdinal(x=xx, y=y, alphas=alphas, lambdaMinRatio=lambdaMinRatio, maxiterOut=maxiterOut)
			
			this <- which.min(model$aic)
			aic <- model$aic[this]
			coeffVal <- model$coefs[this, vars]
			alpha <- model$args$alpha
			
			results <- rbind(
				results,
				data.frame(
					modelNum=modelNum,
					occWindow=occWindow,
					aic=aic,
					alpha=alpha,
					numHomeRanges=numHomeRanges,
					region=region,
					coeff=vars,
					coeffVal=coeffVal
				)
			)
		
	}
	
	rownames(results) <- NULL
	
	write.csv(results, './Figures & Tables/Occupancy - Multivariate/Occupancy - Multivariate Ordinal Models Using All Data.csv', row.names=FALSE)
say('###############################################################')
say('### BINARY multivariate OCCUPANCY analysis: all-data models ###')
say('###############################################################')

	load('./Data/03 New Mexico Pika - Assigned Folds.rda')
	pika$presAbs <- NA
	pika$presAbs <- ifelse(pika$latestOccStatus %in% c('0 never', '1 old'), 0, 1)
	pika$presAbs <- as.numeric(pika$presAbs)
	pika$region <- as.factor(pika$region)
	pika$numHomeRangesScaled <- scale(pika$numHomeRanges)
	
	y <- pika$presAbs

	# predictors
	file <- './Predictors/Predictor Variables 2021-04-28 Edits by Adam.xlsx'
	univariatePredTable <- read_excel(file, sheet='Predictor Variables', skip=1)

	results <- data.frame()

	### model
	#########

	for (occWindow in c(occWindows_y, NA)) {
	# for (occWindow in c(occWindows_y)) {

		vars <- names(pika)[grepl(names(pika), pattern='occVar_')]
		if (!is.na(occWindow)) vars <- vars[grepl(vars, pattern=paste0('_', occWindow, 'yrWindow'))]

		x <- scale(pika[ , vars])
		colnames(x) <- vars

		## climate only
		
			say('climate model ', occWindow, '-yr window')

			modelNum <- 1
			numHomeRanges <- FALSE
			region <- FALSE

			# implement elastic net for binary model
			thisResults <- optBinary(x, y, alphas=alphas)
			
			results <- rbind(
				results,
				cbind(
					data.frame(
						modelNum=modelNum,
						occWindow=occWindow,
						numHomeRanges=numHomeRanges,
						region=region
					),
					thisResults
				)
			)
			
		### climate + number of home ranges
		
			say('climate + HR model ', occWindow, '-yr window')

			modelNum <- 2
			numHomeRanges <- TRUE
			region <- FALSE

			hr <- matrix(pika$numHomeRangesScaled, ncol=1)
			colnames(hr) <- 'numHomeRanges'
			xx <- cbind(x, hr)

			# implement elastic net for binary model
			thisResults <- optBinary(xx, y, alphas=alphas)
			
			results <- rbind(
				results,
				cbind(
					data.frame(
						modelNum=modelNum,
						occWindow=occWindow,
						numHomeRanges=numHomeRanges,
						region=region
					),
					thisResults
				)
			)
			
		### climate + region
		
			say('climate + region model ', occWindow, '-yr window')

			modelNum <- 3
			numHomeRanges <- FALSE
			region <- TRUE

			regionMat <- matrix(0, nrow=nrow(x), ncol=4)
			colnames(regionMat) <- levels(pika$region)
			for (i in 1:nrow(regionMat)) {
				regionMat[i, colnames(regionMat) == pika$region[i]] <- 1
			}
			xx <- cbind(x, regionMat)

			# implement elastic net for binary model
			thisResults <- optBinary(xx, y, alphas=alphas)
			
			results <- rbind(
				results,
				cbind(
					data.frame(
						modelNum=modelNum,
						occWindow=occWindow,
						numHomeRanges=numHomeRanges,
						region=region
					),
					thisResults
				)
			)
			
		### climate + home ranges + region

			say('climate + HR + region model ', occWindow, '-yr window')

			modelNum <- 4
			numHomeRanges <- TRUE
			region <- TRUE

			xx <- cbind(x, hr, regionMat)

			# implement elastic net for binary model
			thisResults <- optBinary(xx, y, alphas=alphas)
			
			results <- rbind(
				results,
				cbind(
					data.frame(
						modelNum=modelNum,
						occWindow=occWindow,
						numHomeRanges=numHomeRanges,
						region=region
					),
					thisResults
				)
			)
		
	}
	
	rownames(results) <- NULL
	
	write.csv(results, './Figures & Tables/Occupancy - Multivariate/Occupancy - Multivariate Binary Models Using All Data.csv', row.names=FALSE)
	
	
say('########################################################################')
say('### ORDINAL multivariate OCCUPANCY analysis: all-data models summary ###')
say('########################################################################')

	# predictors
	file <- './Predictors/Predictor Variables 2021-04-28 Edits by Adam.xlsx'
	univariatePredTable <- read_excel(file, sheet='Predictor Variables', skip=1)

	results <- read.csv('./Figures & Tables/Occupancy - Multivariate/Occupancy - Multivariate Ordinal Models Using All Data.csv')
	
	minAic <- min(results$aic)
	best <- results[which(results$aic == minAic), ]
	
	best <- best[order(abs(best$coeffVal)), ]

	if (any(best$coeff == 'nw')) best <- best[-which(best$coeff == 'nw'), ] 
	if (any(best$coeff == 'sw')) best <- best[-which(best$coeff == 'sw'), ] 
	if (any(best$coeff == 'ne')) best <- best[-which(best$coeff == 'ne'), ] 
	if (any(best$coeff == 'se')) best <- best[-which(best$coeff == 'se'), ] 
	if (any(best$coeff == 'numHomeRanges')) best <- best[-which(best$coeff == 'numHomeRanges'), ] 
	
	incTime <- (nrow(best) > length(univariatePredTable$var[univariatePredTable$useOccAbund]))

	best$varNice <- makeNiceVars(best$coeff, occOrDens='occ', incTime=TRUE)
	best$varNice <- factor(best$varNice, levels=best$varNice)
	
	coeffs <- ggplot(best, aes(x=coeffVal, y=varNice)) +
		geom_point() +
		xlab('Coefficient value') + ylab(NULL) +
		theme(axis.text = element_text(size=7), axis.title = element_text(size = 8))
		
	pdf('./Figures & Tables/Occupancy - Multivariate/Occupancy Multivariate Ordinal Model Coefficients.pdf', width=3.25, height=4.5)
		print(coeffs)
	dev.off()

say('#######################################################################')
say('### BINARY multivariate OCCUPANCY analysis: all-data models summary ###')
say('#######################################################################')

	# predictors
	file <- './Predictors/Predictor Variables 2021-04-28 Edits by Adam.xlsx'
	univariatePredTable <- read_excel(file, sheet='Predictor Variables', skip=1)

	results <- read.csv('./Figures & Tables/Occupancy - Multivariate/Occupancy - Multivariate Binary Models Using All Data.csv')
	
	minMse <- min(results$mse)
	best <- results[which(results$mse == minMse), ]
	
	best <- best[order(abs(best$coeffVal)), ]

	if (any(best$coeff == 'nw')) best <- best[-which(best$coeff == 'nw'), ] 
	if (any(best$coeff == 'sw')) best <- best[-which(best$coeff == 'sw'), ] 
	if (any(best$coeff == 'ne')) best <- best[-which(best$coeff == 'ne'), ] 
	if (any(best$coeff == 'se')) best <- best[-which(best$coeff == 'se'), ] 
	if (any(best$coeff == 'numHomeRanges')) best <- best[-which(best$coeff == 'numHomeRanges'), ] 
	
	incTime <- (nrow(best) > length(univariatePredTable$var[univariatePredTable$useOccAbund]))

	best$varNice <- makeNiceVars(best$coeff, occOrDens='occ', incTime=TRUE)
	best$varNice <- factor(best$varNice, levels=best$varNice)
	
	coeffs <- ggplot(best, aes(x=coeffVal, y=varNice)) +
		geom_point() +
		xlab('Coefficient value') + ylab(NULL) +
		theme(axis.text = element_text(size=7), axis.title = element_text(size = 8))

	pdf('./Figures & Tables/Occupancy - Multivariate/Occupancy Multivariate Binary Model Coefficients.pdf', width=3.25, height=4.5)
		print(coeffs)
	dev.off()

say('DONE!!!', level=1, deco='%')
