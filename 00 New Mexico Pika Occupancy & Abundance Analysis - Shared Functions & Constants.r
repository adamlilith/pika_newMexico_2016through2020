### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('C:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r')
### source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r')
###
### CONTENTS ###
### setup ###
### convert variable name to nice name ###
### convert name of month to month number ###
### sum of days meeting a particular condition ###
### implement elastic net for ordinal model ###
### implement elastic net for binary model ###
### implement elastic net for density model ###
### extract PRISM values and calculate derived variable for OCCUPANCY predictors ###
### extract PRISM values and calculate derived variable for DENSITY predictors ###
### generate variables ###
### generate models: occupancy ###
### generate models: density ###
### make formulae with nice names ###
### extract coefficients and terms for a set of 4 models ###

#############
### setup ###
#############

	cat(date(), '\n'); flush.console()
	# rm(list=ls())
	gc()
	options(stringsAsFactors=FALSE)
	
	drive <- 'C:/Ecology/'
	# drive <- 'E:/Adam/'
	
	setwd(paste0(drive, '/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)'))

	# library(rainbow) # www.github.com/adamlilith/rainbow
	library(legendary) # www.github.com/adamlilith/legendary
	# # library(enmSdm) # www.github.com/adamlilith/enmSdm
	library(enmSdmX) # www.github.com/adamlilith/enmSdmX
	
	library(cowplot)
	library(data.table)
	library(geodata)
	library(ggcorrplot)
	library(ggnewscale)
	library(ggplot2)
	library(ggspatial)
	library(glmnet)
	library(lubridate)
	library(MASS)
	library(MuMIn)
	library(ordinalNet)
	library(omnibus)
	library(patchwork)
	library(readxl)
	library(raster)
	library(sf)
	library(statisfactory)
	library(terra)
	# library(tidyverse)
	
	# # only needed for extracting from PRISM
	# ff <- listFiles(paste0(drive, '/ecology/Drive/R/airUpThere/R'))
	# for (f in ff) source(f)
	
	ll <- c('longitude', 'latitude')
	surveyYears <- 2016:2020
	occWindows_y <- c(7, 10) # time span for occupancy variables (calculate up to these years prior to survey year)
	densWindows_y <- c(0, 1) # time offset for density variables (calculate these years prior to survey year)
	regions <- c('northwest', 'southwest', 'northeast', 'southeast')
	occStatuses <- c('0 never', '1 old', '2 occupied')

	# predictor table
	file <- './Predictors/Univariate Predictor Variables 2021-12-10 Added Flag for Density Predictors with 0- or 1-yr Lag.xlsx'
	predTable <- read_excel(file, sheet='Predictor Variables', skip=1)
	predTable <- predTable[(predTable$var %notin% c('subLethalHeat22deg_d', 'subLethalHeat20deg_d')), ]

	dirCreate('./Figures & Tables')
	dirCreate('./Figures & Tables/Occupancy - Simple Models')
	dirCreate('./Figures & Tables/Occupancy - Simple Models - Accounting for Spatial Redundancies')
	dirCreate('./Figures & Tables/Density - Simple Models')

	# colors
	occCol <- 'chartreuse' # occupied sites
	oldCol <- 'darkgoldenrod3' # old occupancy
	neverCol <- 'firebrick3' # never occupied
	
	# for post hoc analyses including isolation, what cutoff of delta AICc to use for selecting top "climate-only" models?
	maxDeltaAic_occupancy <- 10
	maxDeltaAic_density <- 3

##########################################
### convert variable name to nice name ###
##########################################

	makeNiceVars <- function(vars, occOrDens, incTime=TRUE, wrapTime=FALSE) {
	
		# vars 		vector of variable names
		# occOrDens 'occupancy' or 'density'
		# incTime	include time frame in parentheses?
		# wrapTime  insert "\n" before time?
	
		origVars <- vars
	
		# time frame
		if (occOrDens == 'occupancy') {
		
			timeFrame <- substr(vars, nchar(vars) - 9, nchar(vars) - 8)
			for (i in seq_along(timeFrame)) if (substr(timeFrame[i], 1, 1) == '_') timeFrame[i] <- substr(timeFrame[i], 2, 2)

			for (occWindow in occWindows_y) {
				vars <- gsub(vars, pattern=paste0('_', occWindow, 'yrWindow'), replacement='')
			}
		
		} else if (occOrDens == 'density') {
		
			timeFrame <- substr(vars, nchar(vars) - 7, nchar(vars) - 7)

			for (densWindow in densWindows_y) {
				vars <- gsub(vars, pattern=paste0('_', densWindow, 'yrPrior'), replacement='')
			}
		
		}

		# variables
		vars <- gsub(vars, pattern='occVar_', replacement='')
		vars <- gsub(vars, pattern='densVar_', replacement='')

		vars <- predTable$varNice[match(vars, predTable$var)]
		
		wrap <- ifelse(wrapTime, '\n', ' ')
		if (incTime & occOrDens == 'occupancy') vars <- paste0(vars, wrap, '(', timeFrame, '-yr)')
		if (incTime & occOrDens == 'density') vars <- paste0(vars, wrap, '(', timeFrame, ' yr)')
		
		if (any(origVars == 'meanDistToClosest4Patches')) vars[origVars == 'meanDistToClosest4Patches'] <- 'isolation'
		
		vars
	
	}

#############################################
### convert name of month to month number ###
#############################################

	monthNameToNum <- function(x) {
	
		# x name of month
	
		mos <- c('january', 'february', 'march', 'april', 'may', 'june', 'july', 'august', 'september', 'october', 'november', 'december')
		which(tolower(x) == mos)
		
	}

##################################################
### sum of days meeting a particular condition ###
##################################################

	#' Sum of days of values meeting a particular condition
	#'
	#' Sum of days of values meeting a particular condition.
	#'
	#' @param x Vector of numeric, character, or other values.
	#' @param fx A function that returns \code{TRUE}, \code{FALSE}, or (optionally) \code{NA}. The function must use \code{x} as its first argument. For example, \code{function(x) x == 0} is allowable, but \code{function(y) y == 0} is not. Values that count as \code{TRUE} will be counted toward a run.
	#' @param args A \emph{list} object with additional arguments to supply to the function \code{fx}.
	#' @param failIfAllNA If \code{TRUE}, fail if all values are \code{NA} after being evaluated by \code{fx}.
	#'
	#' @return Lengths of successive runs of elements that meet the criterion. A single value of 0 indicates no conditions meet the criterion.
	#' @examples
	#'
	#' x <- c(1, 4, 0, 0, 0, 2, 0, 10)
	#' fx <- function(x) x == 0
	#' totalDays(x, fx)
	#' 
	#' fx <- function(x) x > 0
	#' totalDays(x, fx)
	#'  
	#' fx <- function(x) x > 0 & x < 5
	#' totalDays(x, fx)
	#' 
	#' x <- c(1, 4, 0, 0, 0, 2, 0, 10)
	#' fx <- function(x, th) x == th
	#' totalDays(x, fx, args=list(th=0))
	#' 
	#' # "count" NA as an observation 
	#' x <- c(1, 4, 0, 0, 0, NA, 0, 10)
	#' fx <- function(x, th) ifelse(is.na(x), FALSE, x == th)
	#' totalDays(x, fx, args=list(th=0))
	#'  
	#' # include NAs as part of a run
	#' x <- c(1, 4, 0, 0, 0, NA, 0, 10)
	#' fx <- function(x, th) ifelse(is.na(x), TRUE, x == th)
	#' totalDays(x, fx, args=list(th=0))
	#'  
	#' @export
	totalDays <- function(x, fx, args=NULL, failIfAllNA = FALSE) {

		theseArgs <- c(list(x=x), args)
		y <- do.call(fx, args=theseArgs)

		if (all(is.na(y))) {
			if (failIfAllNA) {
				stop('All evaluated values are NA.')
			} else {
				out <- 0
			}
		} else {
			out <- sum(y, na.rm=TRUE)
		}
		out

	}

###############################################
### implement elastic net for ordinal model ###
###############################################

	optOrdinal <- function(x, y, alphas=seq(0, 1, by=0.05), lambdaMinRatio=0.001, includeLambda0=TRUE, maxiterOut=1000) {
	
		models <- list()
		summaries <- data.frame()
		for (i in seq_along(alphas)) {
		
			alpha <- alphas[i]
			models[[i]] <- ordinalNet(x, y, alpha=alpha, lambdaMinRatio=lambdaMinRatio, includeLambda0=includeLambda0, maxiterOut=includeLambda0, reverse=TRUE)
			# models[[i]] <- ordinalNetCV(x, y, alpha=alpha, lambdaMinRatio=lambdaMinRatio, includeLambda0=includeLambda0, maxiterOut=includeLambda0, tuneMethod='aic', printProgress=FALSE)
			summ <- summary(models[[i]])
			summ <- cbind(data.frame(i=rep(i, nrow(summ)), alpha=rep(alpha, nrow(summ))), summ)
			summaries <- rbind(summaries, summ)
		
		}
		
		bestIndex <- which.min(summaries$aic)
		bestModel <- summaries$i[bestIndex]
		model <- models[[bestModel]]
		model
	
	}
		
##############################################
### implement elastic net for binary model ###
##############################################

	optBinary <- function(x, y, alphas=seq(0, 1, by=0.05)) {
	
		models <- list()
		summaries <- data.frame()
		for (i in seq_along(alphas)) {
		
			alpha <- alphas[i]
			models[[i]] <- cv.glmnet(x=x, y=y, family='binomial', alpha=alpha, standardize=FALSE)
			best <- models[[i]]$index['min', 'Lambda']
			
			mse <- models[[i]]$cvm[best]
			summ <- models[[i]]$glmnet.fit$beta[ , best, drop=FALSE]
			summ <- as.matrix(summ)
			colnames(summ) <- 'coeffVal'
			n <- nrow(summ)
			summ <- cbind(data.frame(alpha=rep(alpha, n), mse=rep(mse, n), coeff=rownames(summ)), summ)
			summaries <- rbind(summaries, summ)
		
		}
		
		bestIndex <- which(summaries$mse == min(summaries$mse))
		out <- summaries[bestIndex, ]
		out
	
	}
		
###############################################
### implement elastic net for density model ###
###############################################

	optDensity <- function(x, y, alphas=seq(0, 1, by=0.05)) {
	
		# y <- log10(y)
	
		models <- list()
		summaries <- data.frame()
		for (i in seq_along(alphas)) {
		
			alpha <- alphas[i]
			models[[i]] <- cv.glmnet(x=x, y=y, family=Gamma(link='log'), alpha=alpha, standardize=FALSE)
			best <- models[[i]]$index['min', 'Lambda']
			
			mse <- models[[i]]$cvm[best]
			summ <- models[[i]]$glmnet.fit$beta[ , best, drop=FALSE]
			summ <- as.matrix(summ)
			colnames(summ) <- 'coeffVal'
			n <- nrow(summ)
			summ <- cbind(data.frame(alpha=rep(alpha, n), mse=rep(mse, n), coeff=rownames(summ)), summ)
			summaries <- rbind(summaries, summ)
		
		}
		
		bestIndex <- which(summaries$mse == min(summaries$mse))
		out <- summaries[bestIndex, ]
		out
	
	}
		
####################################################################################
### extract PRISM values and calculate derived variable for OCCUPANCY predictors ###
####################################################################################

	extractAndCalcVarForOcc <- function(
		thisPred,					# pika predictor name
		prVar,						# relevant PRISM variable
		predTable,					# predictor table
		priorYears,					# series of years prior to sample year (eg, 7:0)
		pika,						# pika data frame with all data
		summaryFx,					# function to summarize values within a year
		...,						# named arguments to pass to summaryFx
		dateField = 'latestOccSurveyYear', # name of field with survey date
		fail = TRUE 				# stop if any value of summaryFx is infinite or NA,
	) {
		
		say('occupancy predictor: ', thisPred, ' up to ', max(priorYears), ' year(s) prior', level=2)

		# start/end months/days
		startMo <- predTable$startMonthOccDens[predTable$var == thisPred]
		startDay <- predTable$startDayOccDens[predTable$var == thisPred]

		endMo <- predTable$endMonthOccDens[predTable$var == thisPred]
		endDay <- predTable$endDayOccDens[predTable$var == thisPred]

		### extract for each period
		if (exists('subSummary')) { rm(subSummary) }
		
		for (yearDelta in priorYears) {

			## get start/end dates for this period
			# for periods beginning in previous year
			thisStartYear <- if (monthAsNum(startMo, FALSE) > monthAsNum(endMo, FALSE)) {
				as.numeric(substr(as.character(pika[ , dateField]), 1, 4)) - yearDelta - 1
			} else {
				as.numeric(substr(as.character(pika[ , dateField]), 1, 4)) - yearDelta
			}
		
			startEndDates <- data.frame(
				startYear = thisStartYear,
				startMo = monthAsNum(startMo, names=FALSE),
				startDay = startDay,

				endYear = as.numeric(substr(as.character(pika[ , dateField]), 1, 4)) - yearDelta,
				endMo = monthAsNum(endMo, names=FALSE),
				endDay = endDay
			)
			
			startDates <- paste0(startEndDates$startYear, '-', prefix(startEndDates$startMo, 2), '-', prefix(startEndDates$startDay, 2))
			endDates <- paste0(startEndDates$endYear, '-', prefix(startEndDates$endMo, 2), '-', prefix(startEndDates$endDay, 2))
			
			startDates <- as.Date(startDates)
			endDates <- as.Date(endDates)
			startEndDatesTogether <- data.frame(startDate=startDates, endDate=endDates)
			
			# extract
			pikaWithPeriod <- cbind(pika, startEndDatesTogether)
			pikaWithPeriod <- vect(pikaWithPeriod, geom = c('longitude', 'latitude'), crs = getCRS('NAD83'))
			
			ext <- prExtractAbsoluteDaily(
				x=pikaWithPeriod,
				startDate='startDate',
				endDate='endDate',
				prDir=prDir,
				vars=prVar,
				res=30,
				rastSuffix = 'tif',
				removeNAs = FALSE,
				verbose = TRUE
			)

			thisSumm <- apply(ext, MARGIN=1, FUN=summaryFx, ...)
			# thisSumm <- apply(ext, MARGIN=1, FUN=summaryFx)

			if (fail) if (any(is.infinite(thisSumm)) | any(is.na(thisSumm))) stop('Infinity! NA-NA-NA-NA~NA!')
			
			subSummary <- if (exists('subSummary')) {
				cbind(subSummary, cbind(thisSumm))
			} else {
				cbind(thisSumm)
			}
			
			# if (yearDelta < max(priorYears)) say('')
			
		}

		colnames(subSummary) <- paste0('delta', priorYears)
		masterSumm <- apply(subSummary, 1, mean)

		pika$DUMMY <- masterSumm
		names(pika)[ncol(pika)] <- paste0('occVar_', thisPred, '_', max(priorYears),'yrWindow')

		pika
		
	}

##################################################################################
### extract PRISM values and calculate derived variable for DENSITY predictors ###
##################################################################################

	extractAndCalcVarForDens <- function(
		thisPred,					# pika predictor name
		prVar,						# relevant PRISM variable
		predTable,					# predictor table
		pika,						# pika data frame with all data
		summaryFx,					# function to summarize values within a year
		...,						# named arguments to pass to summaryFx
		fail = TRUE 				# stop if any value of summaryFx is infinite or NA,
	) {
		
		say('density predictor: ', thisPred, level=2)

		### indices of density surveys
		###############################
		
			densSurveyIndices <- which(!is.na(pika$latestDensSurveyYear))

		### YEAR PRIOR
		##############

			yearDelta <- 1 # 1 year prior
		
			# start/end months/days
			startMo <- predTable$startMonthOccDens[predTable$var == thisPred]
			startDay <- predTable$startDayOccDens[predTable$var == thisPred]
		
			endMo <- predTable$endMonthOccDens[predTable$var == thisPred]
			endDay <- predTable$endDayOccDens[predTable$var == thisPred]
		
			## get start/end dates for this period
			# for periods beginning in previous year
			thisStartYear <- if (monthAsNum(startMo, FALSE) > monthAsNum(endMo, FALSE)) {
				pika$latestDensSurveyYear[densSurveyIndices] - yearDelta - 1
			} else {
				pika$latestDensSurveyYear[densSurveyIndices] - yearDelta
			}
		
			startEndDates <- data.frame(
				startYear = thisStartYear,
				startMo = monthAsNum(startMo, names=FALSE),
				startDay = startDay,

				endYear = pika$latestDensSurveyYear[densSurveyIndices] - yearDelta,
				endMo = monthAsNum(endMo, names=FALSE),
				endDay = endDay
			)
			
			startDates <- paste0(startEndDates$startYear, '-', prefix(startEndDates$startMo, 2), '-', prefix(startEndDates$startDay, 2))
			endDates <- paste0(startEndDates$endYear, '-', prefix(startEndDates$endMo, 2), '-', prefix(startEndDates$endDay, 2))
			startEndDatesTogether <- data.frame(startDate=startDates, endDate=endDates)
			
			# extract
			pikaWithPeriod <- cbind(pika[densSurveyIndices, ll], startEndDatesTogether)
			pikaWithPeriod <- vect(pikaWithPeriod, geom = c('longitude', 'latitude'), crs = getCRS('NAD83'))
			
			ext <- prExtractAbsoluteDaily(
				x=pikaWithPeriod,
				startDate='startDate',
				endDate='endDate',
				prDir=prDir,
				vars=prVar,
				res=30,
				rastSuffix = 'tif',
				removeNAs = FALSE,
				verbose = TRUE
			)
			thisSumm <- apply(ext, MARGIN=1, FUN=summaryFx, ...)

			if (fail) if (any(is.infinite(thisSumm)) | any(is.na(thisSumm))) stop('Infinity! NA-NA-NA-NA~NA!')
				
			pika$DUMMY <- NA
			pika$DUMMY[densSurveyIndices] <- thisSumm
			names(pika)[ncol(pika)] <- paste0('densVar_', thisPred, '_', yearDelta,'yrPrior')

		### SAME YEAR
		#############
		
			if (predTable$useAbundSameYear[predTable$var == thisPred]) {

				say('')
				yearDelta <- 0

				# start/end months/days
				startMo <- predTable$startMonthDens[predTable$var == thisPred]
				startDay <- predTable$startDayDens[predTable$var == thisPred]
			
				endMo <- predTable$endMonthDens[predTable$var == thisPred]
				endDay <- predTable$endDayDens[predTable$var == thisPred]
			
				## get start/end dates for this period
				# for periods beginning in previous year
				thisStartYear <- if (monthAsNum(startMo, FALSE) > monthAsNum(endMo, FALSE)) {
					pika$latestDensSurveyYear[densSurveyIndices] - yearDelta - 1
				} else {
					pika$latestDensSurveyYear[densSurveyIndices] - yearDelta
				}
			
				startEndDates <- data.frame(
					startYear = thisStartYear,
					startMo = monthAsNum(startMo, names=FALSE),
					startDay = startDay,

					endYear = pika$latestDensSurveyYear[densSurveyIndices] - yearDelta,
					endMo = monthAsNum(endMo, names=FALSE),
					endDay = endDay
				)
				
				startDates <- paste0(startEndDates$startYear, '-', prefix(startEndDates$startMo, 2), '-', prefix(startEndDates$startDay, 2))
				endDates <- paste0(startEndDates$endYear, '-', prefix(startEndDates$endMo, 2), '-', prefix(startEndDates$endDay, 2))
				startEndDatesTogether <- data.frame(startDate=startDates, endDate=endDates)
				
				# extract
				pikaWithPeriod <- cbind(pika[densSurveyIndices, ll], startEndDatesTogether)
				pikaWithPeriod <- vect(pikaWithPeriod, geom = c('longitude', 'latitude'), crs = getCRS('NAD83'))
				
				ext <- prExtractAbsoluteDaily(
					x=pikaWithPeriod,
					startDate='startDate',
					endDate='endDate',
					prDir=prDir,
					vars=prVar,
					res=30,
					rastSuffix = 'tif',
					removeNAs = FALSE,
					verbose = TRUE
				)

				thisSumm <- apply(ext, MARGIN=1, FUN=summaryFx, ...)

				if (fail) if (any(is.infinite(thisSumm)) | any(is.na(thisSumm))) stop('Infinity! NA-NA-NA-NA~NA!')
				
				pika$DUMMY <- NA
				pika$DUMMY[densSurveyIndices] <- thisSumm
				names(pika)[ncol(pika)] <- paste0('densVar_', thisPred, '_', yearDelta,'yrPrior')
				
			} # if wanting this year's environmental values

		pika
		
	}

##########################
### generate variables ###
##########################

	### return a list of variables
	getVars <- function(occOrDens) {
		
		# occOrDens	'occupancy' or 'density'

		univars <- predTable$var[predTable$useOccAbund]
		
		if (occOrDens == 'occupancy') {
			univars <- rep(univars, each=2)
			univars <- paste0(univars, '_', occWindows_y, 'yrWindow')
		} else if (occOrDens == 'density') {
			newVars <- character()
			for (i in 1:nrow(predTable)) {
				if (predTable$useOccAbund[i]) {
					yearsPrior <- if (predTable$densityVarYearsPrior[i] == 'zero') {
						0
					} else if (predTable$densityVarYearsPrior[i] == 'one') {
						1
					} else if (predTable$densityVarYearsPrior[i] == 'both') {
						0:1
					}
					thisVar <- paste0(predTable$var[i], '_', yearsPrior, 'yrPrior')
					newVars <- c(newVars, thisVar)
				}
			}
			univars <- newVars
		}

		univars <- if (occOrDens == 'occupancy') {
			paste0('occVar_', univars)
		} else if (occOrDens == 'density') {
			paste0('densVar_', univars)
		}
		
		# univars <- if (occOrDens == 'occupancy') {
			# c(univars, 'numHomeRangesScaled', 'meanDistToClosest4Patches')
		# } else if (occOrDens == 'density') {
			# c(univars, 'meanDistToClosest4Patches')
		# }
		univars
		
	}

##################################
### generate models: occupancy ###
##################################

	getFormulaeOcc <- function() {
	
		### univariate models
		univars <- getVars('occupancy')
		
		### bivariate
		bivars <- read_xlsx(paste0('./Predictors/Bivariate Model Suite__3June2021 (002) [2021-06-07b occupancy].xlsx'))
		bivars <- as.data.frame(bivars)
		
		# remove redundant models
		for (i in 1:nrow(bivars)) {
			term1 <- bivars$term1[i]
			term2 <- bivars$term2[i]
			if (term1 > term2) {
				bivars$term1[i] <- term2
				bivars$term2[i] <- term1
			}
		}
		
		dups <- duplicated(bivars[ , c('term1', 'term2')])
		bivars <- bivars[!dups, ]
		term1 <- bivars$term1
		term2 <- bivars$term2

		term1 <- paste0('occVar_', term1)
		term2 <- paste0('occVar_', term2)
		
		term1 <- paste0(rep(term1, each=2), '_', occWindows_y, 'yrWindow')
		term2 <- paste0(rep(term2, each=2), '_', occWindows_y, 'yrWindow')
			
		bivars <- paste0(term1, ' + ', term2)
		
		### trivariate
		trivars <- read_xlsx(paste0('./Predictors/Trivariate 3-predictor models (including some interactions) [2021-06-07b occupancy].xlsx'))
		trivars <- as.data.frame(trivars)
		
		# remove redundant models
		for (i in 1:nrow(trivars)) {
			term1 <- trivars$term1[i]
			term2 <- trivars$term2[i]
			if (term1 > term2) {
				trivars$term1[i] <- term2
				trivars$term2[i] <- term1
			}
		}
		
		dups <- duplicated(trivars[ , c('term1', 'term2', 'term3')])
		trivars <- trivars[!dups, ]
		term1 <- trivars$term1
		term2 <- trivars$term2
		term3 <- trivars$term3

		term1 <- paste0('occVar_', term1)
		term2 <- paste0('occVar_', term2)
		term3 <- paste0('occVar_', term3)

		term1 <- paste0(rep(term1, each=2), '_', occWindows_y, 'yrWindow')
		term2 <- paste0(rep(term2, each=2), '_', occWindows_y, 'yrWindow')
		term3 <- paste0(rep(term3, each=2), '_', occWindows_y, 'yrWindow')

		trivars <- paste0(term1, ' ', rep(trivars$fxBetweenTerms1and2, each=2), ' ', term2, ' + ', term3)

		models <- c(univars, bivars, trivars)
		
		# models <- c(
			# models,
			# paste0(models, ' + numHomeRangesScaled')
		# )
		
		# models <- c(
			# models,
			# paste0(models, ' + meanDistToClosest4Patches')
		# )
		
		models
		
	}

################################
### generate models: density ###
################################

	getFormulaeDens <- function() {
	
		models <- read_xlsx(paste0('./Predictors/!Models for Density 2021-12-09.xlsx'))
		models <- as.data.frame(models)
	
		### remove redundant bivariate models
		bivars <- models[models$modelOrder == 'bivariate', ]
		
		for (i in 1:nrow(bivars)) {
			term1 <- bivars$term1[i]
			term2 <- bivars$term2[i]
			if (term1 > term2) {
			
				bivars$term1[i] <- term2
				bivars$term2[i] <- term1
				
				yearsPriorTerm1 <- bivars$yearsPriorTerm1[i]
				yearsPriorTerm2 <- bivars$yearsPriorTerm2[i]
				bivars$yearsPriorTerm1[i] <- yearsPriorTerm2
				bivars$yearsPriorTerm2[i] <- yearsPriorTerm1
				
			}
		}
		
		dups <- duplicated(bivars[ , c('term1', 'term2')])
		bivars <- bivars[!dups, ]

		### remove redundant trivariate models
		trivars <- models[models$modelOrder == 'trivariate', ]
		
		for (i in 1:nrow(trivars)) {
			term1 <- trivars$term1[i]
			term2 <- trivars$term2[i]
			if (term1 > term2) {
				trivars$term1[i] <- term2
				trivars$term2[i] <- term1

				yearsPriorTerm1 <- trivars$yearsPriorTerm1[i]
				yearsPriorTerm2 <- trivars$yearsPriorTerm2[i]
				trivars$yearsPriorTerm1[i] <- yearsPriorTerm2
				trivars$yearsPriorTerm2[i] <- yearsPriorTerm1
			}
		}
		
		dups <- duplicated(trivars[ , c('term1', 'term2', 'term3')])
		trivars <- trivars[!dups, ]

		### add years-prior windows: bivariate
		bivarModels <- character()
		for (i in 1:nrow(bivars)) {
		
			term1 <- if (bivars$yearsPriorTerm1[i] == 'zero') {
				paste0('densVar_', bivars$term1[i], '_0yrPrior')
			} else if (bivars$yearsPriorTerm1[i] == 'one') {
				paste0('densVar_', bivars$term1[i], '_1yrPrior')
			} else if (bivars$yearsPriorTerm1[i] == 'both') {
				c(
					paste0('densVar_', bivars$term1[i], '_0yrPrior'),
					paste0('densVar_', bivars$term1[i], '_1yrPrior')
				)
			}
			
			term2 <- if (bivars$yearsPriorTerm2[i] == 'zero') {
				paste0('densVar_', bivars$term2[i], '_0yrPrior')
			} else if (bivars$yearsPriorTerm2[i] == 'one') {
				paste0('densVar_', bivars$term2[i], '_1yrPrior')
			} else if (bivars$yearsPriorTerm2[i] == 'both') {
				c(
					paste0('densVar_', bivars$term2[i], '_0yrPrior'),
					paste0('densVar_', bivars$term2[i], '_1yrPrior')
				)
			}
			
			theseModels <- if (length(term1 < 2) | length(term2) < 2) {
				paste(term1, term2, sep=paste0(' ', bivars$fxBetweenTerm1AndTerm2[i], ' '))
			} else {
				c(
					paste(term1[1], term2[1], sep=paste0(' ', bivars$fxBetweenTerm1AndTerm2[i], ' ')),
					paste(term1[2], term2[1], sep=paste0(' ', bivars$fxBetweenTerm1AndTerm2[i], ' ')),
					paste(term1[1], term2[2], sep=paste0(' ', bivars$fxBetweenTerm1AndTerm2[i], ' ')),
					paste(term1[2], term2[2], sep=paste0(' ', bivars$fxBetweenTerm1AndTerm2[i], ' '))
				)
			}
			
			bivarModels <- c(bivarModels, theseModels)
			
		} # next bivariate model form
		
		### add years-prior windows: trivariate
		trivarModels <- character()
		for (i in 1:nrow(trivars)) {
		
			term1 <- if (trivars$yearsPriorTerm1[i] == 'zero') {
				paste0('densVar_', trivars$term1[i], '_0yrPrior')
			} else if (trivars$yearsPriorTerm1[i] == 'one') {
				paste0('densVar_', trivars$term1[i], '_1yrPrior')
			} else if (trivars$yearsPriorTerm1[i] == 'both') {
				c(
					paste0('densVar_', trivars$term1[i], '_0yrPrior'),
					paste0('densVar_', trivars$term1[i], '_1yrPrior')
				)
			}
			
			term2 <- if (trivars$yearsPriorTerm2[i] == 'zero') {
				paste0('densVar_', trivars$term2[i], '_0yrPrior')
			} else if (trivars$yearsPriorTerm2[i] == 'one') {
				paste0('densVar_', trivars$term2[i], '_1yrPrior')
			} else if (trivars$yearsPriorTerm2[i] == 'both') {
				c(
					paste0('densVar_', trivars$term2[i], '_0yrPrior'),
					paste0('densVar_', trivars$term2[i], '_1yrPrior')
				)
			}
			
			term3 <- if (trivars$yearsPriorTerm3[i] == 'zero') {
				paste0('densVar_', trivars$term3[i], '_0yrPrior')
			} else if (trivars$yearsPriorTerm3[i] == 'one') {
				paste0('densVar_', trivars$term3[i], '_1yrPrior')
			} else if (trivars$yearsPriorTerm3[i] == 'both') {
				c(
					paste0('densVar_', trivars$term3[i], '_0yrPrior'),
					paste0('densVar_', trivars$term3[i], '_1yrPrior')
				)
			}
			
			theseModels <- if (length(term1 < 2) | length(term2) < 2) {
				paste(term1, term2, sep=paste0(' ', trivars$fxBetweenTerm1AndTerm2[i], ' '))
			} else {
				c(
					paste(term1[1], term2[1], sep=paste0(' ', trivars$fxBetweenTerm1AndTerm2[i], ' ')),
					paste(term1[2], term2[1], sep=paste0(' ', trivars$fxBetweenTerm1AndTerm2[i], ' ')),
					paste(term1[1], term2[2], sep=paste0(' ', trivars$fxBetweenTerm1AndTerm2[i], ' ')),
					paste(term1[2], term2[2], sep=paste0(' ', trivars$fxBetweenTerm1AndTerm2[i], ' '))
				)
			}
			
			theseModels <- paste(theseModels, rep(term3, each=length(theseModels)), sep=' + ')
			trivarModels <- c(trivarModels, theseModels)
			
		} # next trivariate model form

		models <- c(bivarModels, trivarModels)
		models
		
	}

#####################################
### make formulae with nice names ###
#####################################

	niceFormulae <- function(x, occOrDens, around1='', around2='') {
	
		# x			formulae (as characters)
		# occOrDens	'occupancy' or 'density'
		# around1 	symbols to put around each term
		# around2 	symbols to put around each term
		
		timeFrame <- rep(NA, length(x))
		if (occOrDens == 'occupancy') {
			for (occWindow_y in occWindows_y) {
				foundIt <- grepl(x, pattern=paste0(occWindow_y, 'yrWindow'))
				if (any(foundIt)) timeFrame[foundIt] <- occWindow_y
			}
		} else if (occOrDens == 'density') {
			for (densWindow_y in densWindows_y) {
				foundIt <- grepl(x, pattern=paste0(densWindow_y, 'yrPrior'))
				if (any(foundIt)) timeFrame[foundIt] <- densWindow_y
			}
		}

		for (i in seq_along(x)) x[i] <- gsub(x[i], pattern=paste0('_', timeFrame[i], 'yrWindow'), replacement='')
		for (i in seq_along(x)) x[i] <- gsub(x[i], pattern=paste0('_', timeFrame[i], 'yrPrior'), replacement='')
		timeFrame <- paste0(' (', timeFrame, '-yr)')
		
		x <- gsub(x, pattern='occVar_', replacement='')
		x <- gsub(x, pattern='densVar_', replacement='')
			
		for (i in 1:nrow(predTable)) {
		
			pattern <- predTable$var[i]
			replaceWith <- predTable$varNice[i]
			replaceWith <- paste0(around1, replaceWith, around2)
			
			x <- gsub(x, pattern=pattern, replacement=replaceWith)
		
		}
		
		x <- paste0(x, timeFrame)
		x
		
	}

##########################################################
### extract coefficients and terms for a set of models ###
##########################################################
extractTerms <- function(...) {
	
	# ...	objects of class "glm"
	
	models <- list(...)
	n <- length(models)
	
	term1 <- term2 <- term3 <- term4 <- term5 <- term6 <- term7 <- term8 <- term9 <- rep(NA, n)
	for (countModel in 1:n) {
	
		coefs <- coefficients(models[[countModel]])
		if (any(names(coefs) == '(Intercept)')) coefs <- coefs[names(coefs) != '(Intercept)']
		if (any(names(coefs) == 'numHomeRangesScaled')) coefs <- coefs[names(coefs) != 'numHomeRangesScaled']
		if (any(names(coefs) == 'meanDistToClosest4Patches')) coefs <- coefs[names(coefs) != 'meanDistToClosest4Patches']
		if (any(names(coefs) == 'regionnorthwest')) coefs <- coefs[names(coefs) != 'regionnorthwest']
		if (any(names(coefs) == 'regionsoutheast')) coefs <- coefs[names(coefs) != 'regionsoutheast']
		if (any(names(coefs) == 'regionsouthwest')) coefs <- coefs[names(coefs) != 'regionsouthwest']

		term <- coefs[1]
		term1[countModel] <- paste(sprintf('%.2f', round(term, 2)), names(term))

		if (length(coefs) > 1) {
			term <- coefs[2]
			term2[countModel] <- paste(sprintf('%.2f', round(term, 2)), names(term))
		}
	
		if (length(coefs) > 2) {
			term <- coefs[3]
			term3[countModel] <- paste(sprintf('%.2f', round(term, 2)), names(term))
		}
	
		if (length(coefs) > 3) {
			term <- coefs[4]
			term4[countModel] <- paste(sprintf('%.2f', round(term, 2)), names(term))
		}
	
		if (length(coefs) > 4) {
			term <- coefs[5]
			term5[countModel] <- paste(sprintf('%.2f', round(term, 2)), names(term))
		}
	
		if (length(coefs) > 5) {
			term <- coefs[6]
			term6[countModel] <- paste(sprintf('%.2f', round(term, 2)), names(term))
		}
	
		if (length(coefs) > 6) {
			term <- coefs[7]
			term7[countModel] <- paste(sprintf('%.2f', round(term, 2)), names(term))
		}
	
		if (length(coefs) > 7) {
			term <- coefs[8]
			term8[countModel] <- paste(sprintf('%.2f', round(term, 2)), names(term))
		}
	
		if (length(coefs) > 8) {
			term <- coefs[9]
			term9[countModel] <- paste(sprintf('%.2f', round(term, 2)), names(term))
		}
	
	}
	
	occOrDens <- if (grepl(term1[1], pattern='occVar_')) {
		'occupancy'
	} else {
		'density'
	}
	
	term1 <- niceFormulae(term1, occOrDens=occOrDens)
	term2 <- niceFormulae(term2, occOrDens=occOrDens)
	term3 <- niceFormulae(term3, occOrDens=occOrDens)
	term4 <- niceFormulae(term4, occOrDens=occOrDens)

	out <- list(term1, term2, term3, term4, term5, term6, term7, term8, term9)
	names(out) <- paste0('term', 1:9)
	out
	
}

say('Loaded shared functions and constants for pika New Mexico analysis.', level=1)
