### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/04 New Mexico Pika Occupancy & Abundance Analysis - Grand Summary.r')
###
### CONTENTS ###
### setup ###
### summarize occupancy ###
### summarize density ###

#############
### setup ###
#############

	drive <- 'E:'
	source(paste0(drive, '/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r'))

	# generalization
	minAicWeight <- 0.05 # minimum AIC weight of a model for it to be included in the "best" set
	topModels <- 5 # top number of models/coefficients to get
	minAbsCoeffValDens <- 0 # minimum absolute value of a coefficient in a multivariate density model for it to be "important" in a multivariate model

say('###########################')
say('### summarize occupancy ###')
say('###########################')

	vars <- getVars('occupancy')

	# score models by number of times they appear in "important" models
	imp <- data.frame(var = vars)
	imp$binaryAic <- NA
	imp$binaryCrossValid <- NA
	imp$binaryCrossValid <- NA
	imp$binaryMultivar <- NA
	imp$ordinalAic <- NA
	imp$ordinalCrossValid <- NA
	imp$ordinalMultivar <- NA
	
	### binary AIC
	##############
	
		results <- read.csv('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Binary Models Using All Data - 10-yr Window.csv')
		
		results <- results[results$weight >= minAicWeight, , drop=FALSE]
	
		for (i in 1:nrow(results)) {
			for (j in 1:nrow(imp)) {
				if (grepl(results$model[i], pattern=imp$var[j])) imp$binaryAic[j] <- 1
			}
		}
		
	### binary cross-validation
	###########################
	
		results <- read.csv('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Binary Models Cross-validation Results.csv')
		results$timeFrame[is.na(results$timeFrame)] <- -1
		
		results <- aggregate(results, by=list(results$model, results$timeFrame, results$numHomeRanges, results$region), FUN=mean)
		results$model <- results$timeFrame <- results$numHomeRanges <- results$region <- NULL
		names(results)[1:4] <- c('model', 'timeFrame', 'numHomeRanges', 'region')
		
		results <- results[order(results$auc, decreasing=TRUE), ]
		results <- results[1:topModels, , drop=FALSE]
	
		for (i in 1:nrow(results)) {
			for (j in 1:nrow(imp)) {
				if (grepl(results$model[i], pattern=imp$var[j])) imp$binaryCrossValid[j] <- 1
			}
		}
		
	### ordinal AIC
	###############
	
		results <- read.csv('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Ordinal Models Using All Data - 7 & 10-yr Windows.csv')
		
		results <- results[results$weight >= minAicWeight, , drop=FALSE]
	
		for (i in 1:nrow(results)) {
			for (j in 1:nrow(imp)) {
				if (grepl(results$model[i], pattern=imp$var[j])) imp$ordinalAic[j] <- 1
			}
		}
		
	### ordinal cross-validation
	############################
	
		results <- read.csv('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Ordinal Models Cross-validation Results.csv')
		results$timeFrame[is.na(results$timeFrame)] <- -1
		
		results <- aggregate(results, by=list(results$model, results$timeFrame, results$numHomeRanges, results$region), FUN=mean)
		results$model <- results$timeFrame <- results$numHomeRanges <- results$region <- NULL
		names(results)[1:4] <- c('model', 'timeFrame', 'numHomeRanges', 'region')
		
		results <- results[order(results$auc, decreasing=TRUE), ]
		results <- results[1:topModels, , drop=FALSE]
	
		for (i in 1:nrow(results)) {
			for (j in 1:nrow(imp)) {
				if (grepl(results$model[i], pattern=imp$var[j])) imp$ordinalCrossValid[j] <- 1
			}
		}
		
	### binary multivariate
	########################
	
		results <- read.csv('./Figures & Tables/Occupancy - Multivariate/Occupancy - Multivariate Binary Models Using All Data.csv')
		
		results <- results[which(results$mse == min(results$mse)), ]
		results <- results[order(abs(results$coeffVal), decreasing=TRUE), ]
		results <- results[1:topModels, ]
	
		for (i in 1:nrow(results)) {
			for (j in 1:nrow(imp)) {
				if (results$coeff[i] == imp$var[j]) imp$binaryMultivar[j] <- 1
			}
		}
		
	### ordinal multivariate
	########################
	
		results <- read.csv('./Figures & Tables/Occupancy - Multivariate/Occupancy - Multivariate Ordinal Models Using All Data.csv')
		
		results <- results[which(results$aic == min(results$aic)), ]
		results <- results[1:topModels, ]
	
		for (i in 1:nrow(results)) {
			for (j in 1:nrow(imp)) {
				if (results$coeff[i] == imp$var[j]) imp$ordinalMultivar[j] <- 1
			}
		}
		
	imp$var <- gsub(imp$var, pattern='occVar_', replacement='')
	imp$sum <- rowSums(imp[ , 2:ncol(imp)], na.rm=TRUE)
	imp <- imp[order(imp$sum, decreasing=TRUE), ]
	write.csv(imp, './Figures & Tables/Summary of Occupancy Results.csv', row.names=FALSE)
	
say('#########################')
say('### summarize density ###')
say('#########################')

	vars <- getVars('density')

	# score models by number of times they appear in "important" models
	imp <- data.frame(var = vars)
	imp$simpleAic <- NA
	imp$crossValid <- NA
	imp$multivar <- NA
	
	### binary AIC
	##############
	
		results <- read.csv('./Figures & Tables/Density - Simple Models/Density - Simple GLMs.csv')
		
		results <- results[results$weight >= minAicWeight, , drop=FALSE]
	
		for (i in 1:nrow(results)) {
			for (j in 1:nrow(imp)) {
				if (grepl(results$model[i], pattern=imp$var[j])) imp$simpleAic[j] <- 1
			}
		}
		
	### simple cross-validation
	###########################
	
		results <- read.csv('./Figures & Tables/Density - Simple Models/Density - Simple GLMs Cross-validation Results.csv')
		
		results <- aggregate(results, by=list(results$model, results$region), FUN=mean)
		results$model <- results$numHomeRanges <- results$region <- NULL
		names(results)[1:3] <- c('model', 'numHomeRanges', 'region')
		
		results <- results[order(results$meanAbsPercError, decreasing=FALSE), ]
		results <- results[1:topModels, , drop=FALSE]
	
		for (i in 1:nrow(results)) {
			for (j in 1:nrow(imp)) {
				if (grepl(results$model[i], pattern=imp$var[j])) imp$crossValid[j] <- 1
			}
		}
		
	### multivariate
	################
	
		results <- read.csv('./Figures & Tables/Density - Multivariate/Density - Multivariate Models Using All Data.csv')
		
		results <- results[which(results$mse == min(results$mse)), ]
		results <- results[which(abs(results$coeffVal) > minAbsCoeffValDens), ]
	
		for (i in 1:nrow(results)) {
			for (j in 1:nrow(imp)) {
				if (results$coeff[i] == imp$var[j]) imp$multivar[j] <- 1
			}
		}
		
	imp$var <- gsub(imp$var, pattern='densVar_', replacement='')
	imp$sum <- rowSums(imp[ , 2:ncol(imp)], na.rm=TRUE)
	imp <- imp[order(imp$sum, decreasing=TRUE), ]
	write.csv(imp, './Figures & Tables/Summary of Density Results.csv', row.names=FALSE)
	
say('DONE!!!', level=1, deco='%')
