### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('C:/Ecology/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/03c Univariate Occupancy and Density Models.r')
### source('E:/Adam/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/03c Univariate Occupancy and Density Models.r')
###
### CONTENTS ###
### setup ###
### univariate occupancy models with standardized coefficients ###
### univariate density models with standardized coefficients ###

#############
### setup ###
#############

	rm(list=ls())

	drive <- 'C:/Ecology/'
	# drive <- 'E:/Adam/'

	source(paste0(drive, '/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r'))

	library(jtools)
	library(huxtable)

say('##################################################################')
say('### univariate OCCUPANCY models with standardized coefficients ###')
say('##################################################################')

	load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')
	pika$region <- as.factor(pika$region)
	
	vars <- getVars('occupancy')
	pika[ , vars] <- scale(pika[ , vars])

	models <- list()
	for (var in vars) {
	
		form <- paste0('presAbs ~ 1 +', var)
		models[[length(models) + 1]] <- glm(form, data=pika, family=binomial)
	
	}

	table <- export_summs(models, error_format = "[{conf.low}, {conf.high}]")
	dirCreate('./Figures & Tables/Standardized Coefficient Tables for Univariate Models')
	saveRDS(table, './Figures & Tables/Standardized Coefficient Tables for Univariate Models/Occupancy.rds')

say('################################################################')
say('### univariate DENSITY models with standardized coefficients ###')
say('################################################################')

	load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')
	pika <- pika[!is.na(pika$latestDensity), ]
	pika$region <- as.factor(pika$region)
	
	vars <- getVars('density')
	pika[ , vars] <- scale(pika[ , vars])

	models <- list()
	for (var in vars) {
	
		form <- paste0('latestDensity ~ 1 +', var)
		models[[length(models) + 1]] <- glm(form, data=pika, family=Gamma(link='log'))
	
	}

	table <- export_summs(models, error_format = "[{conf.low}, {conf.high}]")
	dirCreate('./Figures & Tables/Standardized Coefficient Tables for Univariate Models')
	saveRDS(table, './Figures & Tables/Standardized Coefficient Tables for Univariate Models/Density.rds')

say('DONE!!!', level=1, deco='%')
