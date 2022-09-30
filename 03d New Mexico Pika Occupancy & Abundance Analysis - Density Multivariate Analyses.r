### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/03d New Mexico Pika Occupancy & Abundance Analysis - Density Multivariate Analyses.r')
###
### CONTENTS ###
### setup ###
### multivariate DENSITY analysis: all-data models ###
### multivariate DENSITY analysis: all-data models summary ###

#############
### setup ###
#############

	drive <- 'E:'
	source(paste0('/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r'))

	### generalization
	alphas <- seq(0, 1, by=0.05)
	
say('######################################################')
say('### multivariate DENSITY analysis: all-data models ###')
say('######################################################')

	# generalization
	lambdaMinRatio <- 0.001
	maxiterOut <- 10000

	load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')
	pika$region <- as.factor(pika$region)
	pika <- pika[!is.na(pika$latestDensity), ]
	pika$numHomeRangesScaled <- scale(pika$numHomeRanges)

	y <- pika$latestDensity
	
	vars <- getVars('density')

	### climate model without region
	###############################

		say('sans region')
		
		x <- scale(pika[ , vars])
		colnames(x) <- vars

		modelNum <- 1
		region <- FALSE
		model <- optDensity(x=x, y=y, alphas=alphas)

		results <- cbind(
			data.frame(
				modelNum=modelNum,
				region=region
			),
			model
		)

	### climate + region
	####################
		
		say('with region')
		
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
					region=region
				),
				model
			)
		)
	
	rownames(results) <- NULL
	
	write.csv(results, './Figures & Tables/Density - Multivariate/Density - Multivariate Models Using All Data.csv', row.names=FALSE)
	
say('##############################################################')
say('### multivariate DENSITY analysis: all-data models summary ###')
say('##############################################################')

	# predictors

	results <- read.csv('./Figures & Tables/Density - Multivariate/Density - Multivariate Models Using All Data.csv')
	
	best <- min(results$mse)
	best <- results[which(results$mse == best), ]
	
	best <- best[order(abs(best$coeffVal)), ]

	if (any(best$coeff == 'northwest')) best <- best[-which(best$coeff == 'northwest'), ] 
	if (any(best$coeff == 'southwest')) best <- best[-which(best$coeff == 'southwest'), ] 
	if (any(best$coeff == 'northeast')) best <- best[-which(best$coeff == 'northeast'), ] 
	if (any(best$coeff == 'southeast')) best <- best[-which(best$coeff == 'southeast'), ] 
	
	best$varNice <- makeNiceVars(best$coeff, occOrDens='density', incTime=TRUE)
	best$varNice <- factor(best$varNice, levels=best$varNice)
	
	coeffs <- ggplot(best, aes(x=coeffVal, y=varNice)) +
		geom_point() +
		xlab('Coefficient value') + ylab(NULL) +
		theme(axis.text = element_text(size=7), axis.title = element_text(size = 8))

	pdf('./Figures & Tables/Density - Multivariate/Density Multivariate Model Coefficients.pdf', width=3.25, height=4.5)
		print(coeffs)
	dev.off()

say('DONE!!!', level=1, deco='%')
