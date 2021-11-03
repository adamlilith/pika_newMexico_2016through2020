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

#############
### setup ###
#############

	cat(date(), '\n'); flush.console()
	memory.limit(memory.limit() * 2^29)
	rm(list=ls())
	gc()
	options(stringsAsFactors=FALSE)
	
	# drive <- 'C:'
	# drive <- 'D:'
	drive <- 'E:'
	
	setwd(paste0(drive, '/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)'))

	.libPaths(paste0(drive, '/Ecology/Drive/R/libraries'))

	# library(rainbow) # www.github.com/adamlilith/rainbow
	library(legendary) # www.github.com/adamlilith/legendary
	library(omnibus) # www.github.com/adamlilith/omnibus
	library(enmSdm) # www.github.com/adamlilith/enmSdm
	library(statisfactory) # www.github.com/adamlilith/statisfactory
	
	library(cowplot)
	# library(elsa)
	library(ggplot2)
	library(ggcorrplot)
	library(glmnet)
	library(lubridate)
	library(MASS)
	library(MuMIn)
	library(ordinalNet)
	library(readxl)
	library(raster)
	library(terra)
	library(tidyverse)
	
	ll <- c('longitude', 'latitude')
	surveyYears <- 2016:2020
	occWindows_y <- c(7, 10) # time span for occupancy variables (calculate up to these years prior to survey year)
	densWindows_y <- c(0, 1) # time offset for density variables (calculate these years prior to survey year)

	# predictor table
	file <- './Predictors/Univariate Predictor Variables 2021-06-09 Modified Sub-lethal Heat.xlsx'
	predTable <- read_excel(file, sheet='Predictor Variables', skip=1)

	dirCreate('./Figures & Tables')
	dirCreate('./Figures & Tables/Occupancy - Univariate')
	dirCreate('./Figures & Tables/Occupancy - Multivariate')
	dirCreate('./Figures & Tables/Density - Univariate')
	dirCreate('./Figures & Tables/Density - Multivariate')

	occCol <- 'chartreuse'
	oldCol <- 'darkgoldenrod3'
	neverCol <- 'firebrick3'

##########################################
### convert variable name to nice name ###
##########################################

	makeNiceVars <- function(vars, occOrDens, incTime=TRUE, wrapTime=FALSE) {
	
		# vars 		vector of variable names
		# occOrDens 'occ' or 'dens'
		# incTime	include time frame in parentheses?
		# wrapTime  insert "\n" before time?
	
		file <- './Predictors/Univariate Predictor Variables 2021-06-09 Modified Sub-lethal Heat.xlsx'
		predTable <- read_excel(file, sheet='Predictor Variables', skip=1)
		
		# time frame
		if (occOrDens == 'occ') {
		
			timeFrame <- substr(vars, nchar(vars) - 9, nchar(vars) - 8)
			for (i in seq_along(timeFrame)) if (substr(timeFrame[i], 1, 1) == '_') timeFrame[i] <- substr(timeFrame[i], 2, 2)

			for (occWindow in occWindows_y) {
				vars <- gsub(vars, pattern=paste0('_', occWindow, 'yrWindow'), replacement='')
			}
		
		} else if (occOrDens == 'dens') {
		
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
		if (incTime & occOrDens == 'occ') vars <- paste0(vars, wrap, '(', timeFrame, '-yr window)')
		if (incTime & occOrDens == 'dens') vars <- paste0(vars, wrap, '(', timeFrame, ' yr prior)')
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

#' Maximum number of continuous "runs" of values meeting a particular condition
#'
#' Consider an ordered set of values, say {1, 4, 0, 0, 0, 2, 0, 10}. We can ask, what is the number of times in which zeroes appear successively? In this example, we have one set of three continuous zeros, and one set of a single zero. So the number of runs with 0 is 2, and the maximum run length is 3. This function calculates the number of runs based on a certain condition for defining the run. The condition is stated as a function that returns a logical value. Continuing the example, \code{function(x) x == 0}.
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
#' maxRuns(x, fx)
#' 
#' fx <- function(x) x > 0
#' maxRuns(x, fx)
#'  
#' fx <- function(x) x > 0 & x < 5
#' maxRuns(x, fx)
#' 
#' x <- c(1, 4, 0, 0, 0, 2, 0, 10)
#' fx <- function(x, th) x == th
#' maxRuns(x, fx, args=list(th=0))
#' 
#' # "count" NA as an observation 
#' x <- c(1, 4, 0, 0, 0, NA, 0, 10)
#' fx <- function(x, th) ifelse(is.na(x), FALSE, x == th)
#' maxRuns(x, fx, args=list(th=0))
#'  
#' # include NAs as part of a run
#' x <- c(1, 4, 0, 0, 0, NA, 0, 10)
#' fx <- function(x, th) ifelse(is.na(x), TRUE, x == th)
#' maxRuns(x, fx, args=list(th=0))
#'  
#' @export
maxRuns <- function(x, fx, args=NULL, failIfAllNA = FALSE) {

	theseArgs <- c(list(x=x), args)
	y <- do.call(fx, args=theseArgs)

	if (all(is.na(y))) {
		if (failIfAllNA) {
			stop('All evaluated values are NA.')
		} else {
			out <- 0
		}
	} else {
		
		rl <- rle(y)
		whichMeetCriteria <- which(!is.na(rl$values) & rl$values)
		out <- if (length(whichMeetCriteria) == 0) {
			0
		} else {
			rl$lengths[whichMeetCriteria]
		}
	}
	
	out <- max(out)
	out

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
			
			ext <- prExtractAbsoluteDaily(
				x=pikaWithPeriod,
				startDate='startDate',
				endDate='endDate',
				prDir=prDir,
				vars=prVar,
				res=30,
				rastSuffix = 'tif',
				longLat = ll,
				verbose = TRUE,
				method='bilinear'
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
			
			ext <- prExtractAbsoluteDaily(
				x=pikaWithPeriod,
				startDate='startDate',
				endDate='endDate',
				prDir=prDir,
				vars=prVar,
				res=30,
				rastSuffix = 'tif',
				longLat = ll,
				verbose = FALSE,
				method='bilinear'
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
				
				ext <- prExtractAbsoluteDaily(
					x=pikaWithPeriod,
					startDate='startDate',
					endDate='endDate',
					prDir=prDir,
					vars=prVar,
					res=30,
					rastSuffix = 'tif',
					longLat = ll,
					verbose = FALSE,
					method='bilinear'
				)

				thisSumm <- apply(ext, MARGIN=1, FUN=summaryFx, ...)

				if (fail) if (any(is.infinite(thisSumm)) | any(is.na(thisSumm))) stop('Infinity! NA-NA-NA-NA~NA!')
				
				pika$DUMMY <- NA
				pika$DUMMY[densSurveyIndices] <- thisSumm
				names(pika)[ncol(pika)] <- paste0('densVar_', thisPred, '_', yearDelta,'yrPrior')
				
			} # if wanting this year's environmental values

		pika
		
	}

