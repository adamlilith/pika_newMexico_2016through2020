### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/Code/03d New Mexico Pika Occupancy & Abundance Analysis - Density Multivariate Analyses.r')
###
### CONTENTS ###
### setup ###
### multivariate DENSITY analysis: all-data models ###
### multivariate DENSITY analysis: all-data models summary ###

#############
### setup ###
#############

	source('/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/Code/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r')

	### generalization
	alphas <- seq(0, 1, by=0.05)
	# alphas <- seq(0, 1, by=1)
	
say('######################################################')
say('### multivariate DENSITY analysis: all-data models ###')
say('######################################################')

	# generalization
	lambdaMinRatio <- 0.001
	maxiterOut <- 10000

	load('./Data/03 New Mexico Pika - Assigned Folds.rda')
	pika$region <- as.factor(pika$region)
	pika <- pika[!is.na(pika$latestDensity), ]
	pika$numHomeRangesScaled <- scale(pika$numHomeRanges)

	# predictors
	file <- './Predictors/Predictor Variables 2021-04-28 Edits by Adam.xlsx'
	univariatePredTable <- read_excel(file, sheet='Predictor Variables', skip=1)

	y <- pika$latestDensity
	
	results <- data.frame()

	### model
	#########

	for (densWindow in c(densWindows_y, NA)) {

		vars <- names(pika)[grepl(names(pika), pattern='densVar_')]
		if (!is.na(densWindow)) vars <- vars[grepl(vars, pattern=paste0('_', densWindow, 'yrPrior'))]

if (any(vars=='densVar_subLethalHeat_d_0yrPrior')) vars <- vars[-which(vars=='densVar_subLethalHeat_d_0yrPrior')]

		x <- scale(pika[ , vars])
		colnames(x) <- vars

		## climate only
		
			say('climate model ', densWindow, '-yr window')

			modelNum <- 1
			region <- FALSE

			# implement elastic net for density model
			model <- optDensity(x=x, y=y, alphas=alphas)
			
			results <- rbind(
				results,
				cbind(
					data.frame(
						modelNum=modelNum,
						densWindow=densWindow,
						region=region
					),
					model
				)
			)
			
		### climate + region
		
			say('climate + region model ', densWindow, '-yr window')

			modelNum <- 2
			region <- TRUE

			regionMat <- matrix(0, nrow=nrow(x), ncol=4)
			colnames(regionMat) <- levels(pika$region)
			for (i in 1:nrow(regionMat)) {
				regionMat[i, colnames(regionMat) == pika$region[i]] <- 1
			}
			xx <- cbind(x, regionMat)

			# implement elastic net for density model
			model <- optDensity(x=xx, y=y, alphas=alphas)
			
			results <- rbind(
				results,
				cbind(
					data.frame(
						modelNum=modelNum,
						densWindow=densWindow,
						region=region
					),
					model
				)
			)
			
		
	}
	
	rownames(results) <- NULL
	
	write.csv(results, './Figures & Tables/Density - Multivariate/Density - Multivariate Models Using All Data.csv', row.names=FALSE)
	
say('##############################################################')
say('### multivariate DENSITY analysis: all-data models summary ###')
say('##############################################################')

	# predictors
	file <- './Predictors/Predictor Variables 2021-04-28 Edits by Adam.xlsx'
	univariatePredTable <- read_excel(file, sheet='Predictor Variables', skip=1)

	results <- read.csv('./Figures & Tables/Density - Multivariate/Density - Multivariate Models Using All Data.csv')
	
	best <- min(results$mse)
	best <- results[which(results$mse == best), ]
	
	best <- best[order(abs(best$coeffVal)), ]

	if (any(best$coeff == 'nw')) best <- best[-which(best$coeff == 'nw'), ] 
	if (any(best$coeff == 'sw')) best <- best[-which(best$coeff == 'sw'), ] 
	if (any(best$coeff == 'ne')) best <- best[-which(best$coeff == 'ne'), ] 
	if (any(best$coeff == 'se')) best <- best[-which(best$coeff == 'se'), ] 
	
	best$varNice <- makeNiceVars(best$coeff, occOrDens='dens', incTime=TRUE)
	best$varNice <- factor(best$varNice, levels=best$varNice)
	
	coeffs <- ggplot(best, aes(x=coeffVal, y=varNice)) +
		geom_point() +
		xlab('Coefficient value') + ylab(NULL) +
		theme(axis.text = element_text(size=7), axis.title = element_text(size = 8))


	pdf('./Figures & Tables/Density - Multivariate/Density Multivariate Model Coefficients.pdf', width=3.25, height=4.5)
		print(coeffs)
	dev.off()

say('DONE!!!', level=1, deco='%')
