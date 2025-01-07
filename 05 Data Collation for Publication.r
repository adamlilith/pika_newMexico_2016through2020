### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('C:/Ecology/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/05 Data Collation for Publication.r')
### source('C:/Subarashi/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/05 Data Collation for Publication.r')
###
### CONTENTS ###
### setup ###
### collate pika data for Beever et al. 202X Ordinal modeling. Ecography ###

#############
### setup ###
#############

	rm(list=ls())

	drive <- 'C:/Ecology/'
	# drive <- 'C:/Subarashi/'
	
	library(omnibus)

	source(paste0(drive, '/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r'))

############################################################################
### collate pika data for Beever et al. 202X Ordinal modeling. Ecography ###
############################################################################

	load('./Data/05 New Mexico Pika - Added PRISM Cell Number & Cell-Based Weight.rda')
	
	pika$surveyDate2016 <- as.Date(pika$surveyDate2016)
	pika$surveyDate2017 <- as.Date(pika$surveyDate2017)
	pika$surveyDate2018 <- as.Date(pika$surveyDate2018)
	pika$surveyDate2019 <- as.Date(pika$surveyDate2019)
	pika$surveyDate2020 <- as.Date(pika$surveyDate2020)

	fields <- c('surveyDate2016', 'surveyDate2017', 'surveyDate2018', 'surveyDate2019', 'surveyDate2020')

	pika$lastSurveyDate <- apply(pika[ , fields], 1, max, na.rm = TRUE)
	
	fields <- c('polygonName', 'longitude', 'latitude', 'latestOccSurveyYear', 'lastSurveyDate', 'latestOccStatus', 'presAbs', 'region', 'elevation_m', 'numHomeRanges', 'meanDistToClosest4Patches')
	
	clim_vars <- names(pika)[grepl(names(pika), pattern = 'occVar_') & grepl(names(pika), pattern = '_10yrWindow')]
	
	pikaCleaned <- pika[ , c(fields, clim_vars)]
	pikaCleaned <- renameCol(pikaCleaned, 'polygonName', 'siteName')
	pikaCleaned <- renameCol(pikaCleaned, 'latestOccSurveyYear', 'surveyYear')
	pikaCleaned <- renameCol(pikaCleaned, 'latestOccStatus', 'ordinalStatus')
	pikaCleaned <- renameCol(pikaCleaned, 'presAbs', 'binaryStatus')
	pikaCleaned <- renameCol(pikaCleaned, 'meanDistToClosest4Patches', 'meanDistToClosest4Patches_m')

	names(pikaCleaned) <- gsub(names(pikaCleaned), pattern = 'occVar_', replacement = '')

	write.csv(pikaCleaned, './Figures & Tables/Pika Site Data Cleaned for Ordinal Modeling.csv', row.names = FALSE)


say('DONE!!!', level=1, deco='%')
