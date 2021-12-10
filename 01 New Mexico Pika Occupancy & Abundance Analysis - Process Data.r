### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('C:/Ecology/Drive/Research Active/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/01 New Mexico Pika Occupancy & Abundance Analysis - Process Data.r')
### source('E:/Ecology/Drive/Research Active/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/01 New Mexico Pika Occupancy & Abundance Analysis - Process Data.r')
###
### CONTENTS ###
### process and clean pika data ###
### extract environmental data and calculate predictors ###
### distributions of predictors ###
### define regions and folds ###

#############
### setup ###
#############

	source('E:/Ecology/Drive/Research Active/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r')

# say('###################################')
# say('### process and clean pika data ###')
# say('###################################')

	# # reformat Excel data into R
	# file <- paste0('./Data/00 NM Database 2021 Master [Downloaded 2021-04-23 ORIGINAL].xlsx')
	# pika <- read_excel(file, sheet='All Patches', skip=4, col_types='text', na=c('N/A', 'NA', ''))

	# pika <- as.data.frame(pika)

	# ### nice names
	# ##############

		# names(pika)[names(pika) == 'Polygon Name'] <- 'polygonName'
		# names(pika)[names(pika) == 'Latitude'] <- 'latitude'
		# names(pika)[names(pika) == 'Longitude'] <- 'longitude'
		# names(pika)[names(pika) == 'Elevation (m)'] <- 'elevation_m'
		# names(pika)[names(pika) == 'Elevation (ft)'] <- 'elevation_ft'
		# names(pika)[names(pika) == 'Mean Aspect'] <- 'aspectMean_deg'
		# names(pika)[names(pika) == 'Mean Slope'] <- 'slopeMean_deg'
		# names(pika)[names(pika) == 'Size (#HR)'] <- 'numHomeRanges'
		# names(pika)[names(pika) == 'Color Talus (most recent survey):'] <- 'origOccCodeMostRecent'
		# names(pika)[names(pika) == '2016 Survey Results (Color)'] <- 'origOccCode2016'
		# names(pika)[names(pika) == '2017 Survey Results (Color)'] <- 'origOccCode2017'
		# names(pika)[names(pika) == '2018 Survey Results (Color)'] <- 'origOccCode2018'
		# names(pika)[names(pika) == '2019 Survey Results (Color)'] <- 'origOccCode2019'
		# names(pika)[names(pika) == '2020 Suvey Results (Color)'] <- 'origOccCode2020'
		# names(pika)[names(pika) == 'Change in Occupancy'] <- 'changeInOcc'
		# names(pika)[names(pika) == 'Patch Notes'] <- 'patchNotes'
		# pika$`...17` <- NULL
		# names(pika)[names(pika) == '2016 Survey Results (Density)'] <- 'origDensity2016'
		# names(pika)[names(pika) == '2017 Survey Results (Density)'] <- 'origDensity2017'
		# names(pika)[names(pika) == '2018 Survey Results (Density)'] <- 'origDensity2018'
		# names(pika)[names(pika) == '2019 Survey Results (Density)'] <- 'origDensity2019'
		# names(pika)[names(pika) == '2020 Survey Results (Density)'] <- 'origDensity2020'
		# pika$`...23` <- NULL
		# names(pika)[names(pika) == 'Grazing Status (Most Recent Survey; Yes=1, No=0)'] <- 'grazingMostRecent'
		# names(pika)[names(pika) == '2016 Grazing Status (Yes=1, No=0)'] <- 'grazing2016'
		# names(pika)[names(pika) == '2017 Grazing Status (Yes=1, No=0)'] <- 'grazing2017'
		# names(pika)[names(pika) == '2018 Grazing Status (Yes=1, No=0)'] <- 'grazing2018'
		# names(pika)[names(pika) == '2019 Grazing Status (Yes=1, No=0)'] <- 'grazing2019'
		# names(pika)[names(pika) == '2020 Grazing Status (Yes=1, No=0)'] <- 'grazing2020'
		# pika$`...30` <- NULL
		# names(pika)[names(pika) == 'Percent Patch Perimeter Burned       (Most Recent Survey)'] <- 'perimBurnedMostRecent_perc'
		# names(pika)[names(pika) == '2016 Percent Patch Perimeter Burned'] <- 'perimBurned2016_perc'
		# names(pika)[names(pika) == '2017 Percent Patch Perimeter Burned'] <- 'perimBurned2017_perc'
		# names(pika)[names(pika) == '2018 Percent Patch Perimeter Burned'] <- 'perimBurned2018_perc'
		# names(pika)[names(pika) == '2019 Percent Patch Perimeter Burned'] <- 'perimBurned2019_perc'
		# names(pika)[names(pika) == '2020 Percent Patch Perimeter Burned'] <- 'perimBurned2020_perc'
		# names(pika)[names(pika) == 'Year of Most Recent Fire (Prior to 2017)'] <- 'mostRecentFireYear'
		# names(pika)[names(pika) == 'Most Recent Fire Name (Prior to 2017)'] <- 'mostRecentFireName'
		# names(pika)[names(pika) == '# Acres Buned in Most Recent Fire (Prior to 2017)'] <- 'mostRecentFireArea_ac'
		# pika$`...40` <- NULL
		# names(pika)[names(pika) == '2016 Patch-Wide Mean Grass Cover, all spp. (%)'] <- 'grass2016_perc'
		# names(pika)[names(pika) == '2016 Patch-Wide Mean Forb Cover, all spp. (%)'] <- 'forb2016_perc'
		# names(pika)[names(pika) == '2016 Patch-Wide Mean Total G+F, all spp. (%)'] <- 'grassForb2016_perc'
		# names(pika)[names(pika) == '2017 Patch-Wide Mean Grass Cover, all spp. (%)'] <- 'grass2017_perc'
		# names(pika)[names(pika) == '2017 Patch-Wide Mean Forb Cover, all spp. (%)'] <- 'forb2017_perc'
		# names(pika)[names(pika) == '2017 Patch-Wide Mean Total G+F, all spp. (%)'] <- 'grassForb2017_perc'
		# names(pika)[names(pika) == '2018 Patch-Wide Mean Grass Cover, all spp. (%)'] <- 'grass2018_perc'
		# names(pika)[names(pika) == '2018 Patch-Wide Mean Forb Cover, all spp. (%)'] <- 'forb2018_perc'
		# names(pika)[names(pika) == '2018 Patch-Wide Mean Total G+F, all spp. (%)'] <- 'grassForb2018_perc'
		# names(pika)[names(pika) == '2019 Patch-Wide Mean Grass Cover, all spp. (%)'] <- 'grass2019_perc'
		# names(pika)[names(pika) == '2019 Patch-Wide Mean Forb Cover, all spp. (%)'] <- 'forb2019_perc'
		# names(pika)[names(pika) == '2019 Patch-Wide Mean Total G+F, all spp. (%)'] <- 'grassForb2019_perc'
		# names(pika)[names(pika) == '2020 Patch-Wide Mean Grass Cover, all spp. (%)'] <- 'grass2020_perc'
		# names(pika)[names(pika) == '2020 Patch-Wide Mean Forb Cover, all spp. (%)'] <- 'forb2020_perc'
		# names(pika)[names(pika) == '2020 Patch-Wide Mean Total G+F, all spp. (%)'] <- 'grassForb2020_perc'
		# pika$`...56` <- NULL
		# names(pika)[names(pika) == '2016 Survey Date'] <- 'surveyDate2016'
		# names(pika)[names(pika) == '2017 Survey Date'] <- 'surveyDate2017'
		# names(pika)[names(pika) == '2018 Survey Date'] <- 'surveyDate2018'
		# names(pika)[names(pika) == '2019 Survey Date'] <- 'surveyDate2019'
		# names(pika)[names(pika) == '2020 Survey Date'] <- 'surveyDate2020'
		# pika$`...62` <- NULL
		# names(pika)[names(pika) == 'Area'] <- 'areaName'
		# names(pika)[names(pika) == 'Best Access Notes'] <- 'bestAccessNotes'
		# names(pika)[names(pika) == 'Added By'] <- 'addedBy'
		# names(pika)[names(pika) == 'Historical Site'] <- 'historicalSite'
		# pika$`...67` <- NULL

	# ### remove all "problem sites"
	# ##############################
	
		# probsStart <- which(pika$polygonName == 'PROBLEM SITES')
		# pika <- pika[1:(probsStart - 1), ]
		
		# empty <- which(is.na(pika$polygonName))
		# if (length(empty) > 0) pika <- pika[-empty, ]

		# if (any(pika$polygonName == 'Misattributed Survey: Not Pecos5')) pika <- pika[-which(pika$polygonName == 'Misattributed Survey: Not Pecos5'), ]
		
	# ### checks
	# ##########
	
		# say('Number of unique polygon names is same as number of rows: ', length(unique(pika$polygonName)) == nrow(pika))
		# say('Number of records with NA for longitude ........ ', sum(is.na(pika$longitude)))
		# say('Number of records with NA for latitude ......... ', sum(is.na(pika$latitude)))
		# say('Number of records with NA for elevation_m ...... ', sum(is.na(pika$origOccCodeMostRecent)))
		# # pika$longitude
		# # pika$latitude
		# # pika$elevation_m
		# # pika$origOccCodeMostRecent
		# # pika$surveyDate2016
		# # pika$surveyDate2017
		# # pika$surveyDate2018
		# # pika$surveyDate2019
		# # pika$surveyDate2020
		
	# ### column-specific changes and data corrections
	# ################################################
	
		# pika$longitude <- as.numeric(pika$longitude)
		# pika$latitude <- as.numeric(pika$latitude)
		# pika$elevation_m <- as.numeric(pika$elevation_m)

		# pika$surveyDate2016 <- as.Date(as.numeric(pika$surveyDate2016), origin ='1900-01-01') # warnings are expected
		# pika$surveyDate2016[pika$polygonName == 'PecosBaldyLake'] <- as.Date('2016-09-10') # original: "7/6/16 and 9/10/16"

		# pika$surveyDate2017 <- as.Date(as.numeric(pika$surveyDate2017), origin ='1900-01-01') # warnings are expected
		# pika$surveyDate2017[pika$polygonName == 'Barb1'] <- as.Date('2017-10-21') # original: "6/14/2017 and 10/21/17"
		# pika$surveyDate2017[pika$polygonName == 'BARB1.5'] <- as.Date('2017-10-21') # original: "6/14/2017 and 10/21/17"
		# pika$surveyDate2017[pika$polygonName == 'Madrid 2'] <- as.Date('2017-06-10') # original: "6/9/2017 - 6/10/17"
		# pika$surveyDate2017[pika$polygonName == 'Redondo103'] <- as.Date('2017-06-07') # original: "6/7/2016" but in wrong  column for it to be in 2016
		
		# pika$surveyDate2018 <- as.Date(as.numeric(pika$surveyDate2018), origin ='1900-01-01') # warnings are expected
		# pika$surveyDate2018[pika$polygonName == 'BAND9'] <- as.Date('2018-07-13') # original: "7/4/2018 and 7/13/2018"
		
		# pika$surveyDate2019 <- as.Date(as.numeric(pika$surveyDate2019), origin ='1900-01-01') # warnings are expected
		# pika$surveyDate2020 <- as.Date(as.numeric(pika$surveyDate2020), origin ='1900-01-01') # warnings are expected
		
		# pika$origDensity2017[pika$polygonCode == 'Flats3'] <- 0 # was "R"
		# pika$origDensity2017[pika$polygonCode == 'Flats4'] <- 0 # was "O"
		
	# ### clean OCCUPANCY codes
	# #########################
		
		# for (year in surveyYears) {

			# site <- pika$polygonName
			
			# x <- pika[ , paste0('origOccCode', year), drop=TRUE]
			# x <- tolower(x)

			# pika$DUMMY <- ifelse(x == 'no survey', NA,
				# ifelse(x == 'r', '0 never',
				# ifelse(x == 'o', '1 old',
				# ifelse(x == 'g', '2 occupied',
				# ifelse(site=='BANDnew' & year == 2016 & x == 'r/o', 'o', # Peter Billman 2021-04-19 says to use "o"
				# ifelse(site=='Redondo20' & year == 2017 & x == 'r/g', 'r/g', # using survey from 2018 anyway
				# '!!!!! INSPECT ME !!!!!'
			# ))))))
			
			# if (any(pika$DUMMY %==na% '!!!!! INSPECT ME !!!!!')) stop(paste0('Indeterminate occupancy code for ', year, '.'))

			# names(pika)[names(pika) == 'DUMMY'] <- paste0('occCode', year)
		# }

	# ### clean DENSITY
	# #################
		
		# for (year in surveyYears) {

			# site <- pika$polygonName
			
			# x <- pika[ , paste0('origDensity', year), drop=TRUE]
			# x <- tolower(x)
			
			# x <- ifelse(x == 'no survey', NA, x)
			# x <- as.numeric(x)
			
			# pika$DUMMY <- ifelse(x == 0, NA, x)
			# names(pika)[names(pika) == 'DUMMY'] <- paste0('density', year)
		
		# }

	# ### date/year of latest OCCUPANCY survey
	# ########################################
		
		# pika$latestOccSurveyYear <- pika$latestOccSurveyDate <- NA
		# for (i in 1:nrow(pika)) {
		
			# sampleDates <- pika[i, paste0('surveyDate', surveyYears), drop=TRUE]
			# notNas <- which(!is.na(sampleDates))
			
			# # some surveys don't have dates of survey, so using the year of latest survey with a site-status code
			# if (length(notNas) == 0) {
			
				# x <- unlist(pika[i, paste0('occCode', surveyYears)])
				# latestSurvey <- tail(surveyYears[!is.na(x)], 1)
				# pika$latestOccSurveyDate[i] <- as.character(latestSurvey)
				# pika$latestOccSurveyYear[i] <- latestSurvey
			
			# } else {
				
				# pika$latestOccSurveyDate[i] <- as.character(sampleDates[[tail(notNas, 1)]])
				# pika$latestOccSurveyYear[i] <- year(sampleDates[[tail(notNas, 1)]])
				
			# }
			
		# }
		
		# nas <- which(is.na(pika$latestOccSurveyDate))
		# say('There are ', length(nas), ' remaining sites surveyed for OCCUPANCY that still do not have a year of survey assigned to them.')

		# ## assign latest occupancy status
		# pika$latestOccStatus <- NA
		# for (i in 1:nrow(pika)) {
			# pika$latestOccStatus[i] <- pika[i, paste0('occCode', pika$latestOccSurveyYear[i])]
		# }

	# ### date/year of latest DENSITY survey
	# ######################################
		
		# pika$latestDensSurveyYear <- pika$latestDensSurveyDate <- NA
		# for (i in 1:nrow(pika)) {
		
			# sampleDates <- pika[i, paste0('surveyDate', surveyYears), drop=TRUE]
			# notNas <- which(!is.na(unlist(sampleDates)))
			
			# # some surveys don't have dates of survey, so using the year of latest survey with a site-status code
			# if (length(notNas) == 0) {
			
				# x <- unlist(pika[i, paste0('density', surveyYears)])
				
				# if (!all(is.na(x))) {
					
					# latestSurvey <- tail(surveyYears[!is.na(x)], 1)
					# pika$latestDensSurveyDate[i] <- as.character(latestSurvey)
					# pika$latestDensSurveyYear[i] <- latestSurvey
					
				# }
				
			# } else {

				# latestDate <- as.character(sampleDates[[tail(notNas, 1)]])
				# latestYear <- year(latestDate)
				
				# if (!is.na(pika[i, paste0('density', latestYear)])) {
					
					# pika$latestDensSurveyDate[i] <- latestDate
					# pika$latestDensSurveyYear[i] <- latestYear
					
				# }
					
			# }
			
		# }
		
		
		# nas <- which(is.na(pika$latestDensSurveyYear))
		# say('There are ', length(nas), ' remaining sites surveyed for DENSITY that still do not have a year of survey assigned to them.')

		# ## assign latest occupancy status
		# pika$latestDensity <- NA
		# for (i in 1:nrow(pika)) {
			# if (!is.na(pika$latestDensSurveyYear[i])) {
				# pika$latestDensity[i] <- pika[i, paste0('density', pika$latestDensSurveyYear[i])]
			# }
		# }

	# ### save
	# write.csv(pika, './Data/01 New Mexico Pika - Cleaned.csv', row.names=FALSE)

# say('###########################################################')
# say('### extract environmental data and calculate predictors ###')
# say('###########################################################')

	# say('This step extracts environmental data from PRISM and calculated the derived climate variables as described in the Excel document created by Erik and Maria.', breaks=80)

	# ### generalization
	# # fail <- FALSE # if TRUE then fail if predictor has any infinite or NA values!
	# fail <- TRUE # if TRUE then fail if predictor has any infinite or NA values!
	
	# say('Fail on NA/infinite value: ', fail)

	# # PRISM daily base directory
	# # prDir <- 'F:/ecology/Climate/PRISM/acquired_2020/an81'
	# prDir <- 'I:/Ecology/Climate/PRISM/working/an81'
	
	# # load airUpThere functions for extracting climate data
	# fxs <- listFiles(paste0(drive, '/Ecology/Drive/R/airUpThere/R'), pattern='.r')
	# for (fx in fxs) source(fx)
	
	# fxs <- listFiles(paste0(drive, '/Ecology/Drive/R/airUpThere/data'), pattern='.rda')
	# for (fx in fxs) load(fx)
	
	# ### data
	# pika <- read.csv('./Data/01 New Mexico Pika - Cleaned.csv')

	# ### extract PRISM for OCCCUPANCY variables
	# ##########################################

		# for (mostPriorYear in occWindows_y) {

			# ### chronicCold_C: Average winter temperature (tmean), November-March
			# #####################################################################
				
				# # generalization
				# thisPred <- 'chronicCold_C' # name of this variable
				# prVar <- 'tmean' # variable to extract from PRISM/TerraClimate
				# priorYears <- mostPriorYear:0 # years along window to extract
				
				# # function to summarize values within a year
				# summaryFx <- mean # may need to put arguments in the actual function call below!

				# pika <- extractAndCalcVarForOcc(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# priorYears=priorYears,
					# pika=pika,
					# args=args,
					# fail=fail,
					# summaryFx=summaryFx,
					# na.rm=TRUE
				# )	

			# ### subLethalCold_d: Length of longest run with winter minimum temperature (tmin) <5°C, September-March
			# #######################################################################################################
				
				# # generalization
				# thisPred <- 'subLethalCold_d' # name of this variable
				# prVar <- 'tmin' # variable to extract from PRISM/TerraClimate
				# priorYears <- mostPriorYear:0 # years along window to extract
				
				# # function to summarize values within a year
				# summaryFx <- maxRuns # may need to put arguments in the actual function call below!
				
				# fx <- function(x, th=th) { x < th }
				# args <- list(th = 5)
				
				# pika <- extractAndCalcVarForOcc(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# priorYears=priorYears,
					# pika=pika,
					# summaryFx=summaryFx,
					# fail=fail,
					# failIfAllNA=TRUE,
					# fx=fx,
					# args=args
				# )	

			# ### acuteCold_d: Number of winter days with minimum temperature (tmin) <-10°C, September-March
			# ##############################################################################################
				
				# # generalization
				# thisPred <- 'acuteCold_d' # name of this variable
				# prVar <- 'tmin' # variable to extract from PRISM/TerraClimate
				# priorYears <- mostPriorYear:0 # years along window to extract
				
				# # function to summarize values within a year
				# summaryFx <- totalDays # may need to put arguments in the actual function call below!
				
				# fx <- function(x, th=th) { x < th }
				# args <- list(th = -10)
				
				# pika <- extractAndCalcVarForOcc(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# priorYears=priorYears,
					# pika=pika,
					# summaryFx=summaryFx,
					# fail=fail,
					# failIfAllNA=TRUE,
					# fx=fx,
					# args=args
				# )	

			# ### acuteHeat_d: Number of days with maximum temperature ≥26°C, January-June
			# ############################################################################
				
				# # generalization
				# thisPred <- 'acuteHeat_d' # name of this variable
				# prVar <- 'tmax' # variable to extract from PRISM/TerraClimate
				# priorYears <- mostPriorYear:1 # years along window to extract
				
				# # function to summarize values within a year
				# summaryFx <- totalDays # may need to put arguments in the actual function call below!
				
				# fx <- function(x, th=th) { x >= th }
				# args <- list(th = 26)
				
				# pika <- extractAndCalcVarForOcc(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# priorYears=priorYears,
					# pika=pika,
					# summaryFx=summaryFx,
					# fail=fail,
					# failIfAllNA=TRUE,
					# fx=fx,
					# args=args
				# )	

			# ### subLethalHeat22deg_d: Longest run of days with mean (tmean) temperature ≥22°C, January-June
			# ###############################################################################################
				
				# # generalization
				# thisPred <- 'subLethalHeat22deg_d' # name of this variable
				# prVar <- 'tmean' # variable to extract from PRISM/TerraClimate
				# priorYears <- mostPriorYear:1 # years along window to extract
				
				# # function to summarize values within a year
				# summaryFx <- maxRuns # may need to put arguments in the actual function call below!
				
				# fx <- function(x, th=th) { x >= th }
				# args <- list(th = 22)
				
				# pika <- extractAndCalcVarForOcc(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# priorYears=priorYears,
					# pika=pika,
					# summaryFx=summaryFx,
					# fail=fail,
					# failIfAllNA=TRUE,
					# fx=fx,
					# args=args
				# )	

			# ### subLethalHeat20deg_d: Longest run of days with mean (tmean) temperature ≥20°C, January-June
			# ###############################################################################################
				
				# # generalization
				# thisPred <- 'subLethalHeat20deg_d' # name of this variable
				# prVar <- 'tmean' # variable to extract from PRISM/TerraClimate
				# priorYears <- mostPriorYear:1 # years along window to extract
				
				# # function to summarize values within a year
				# summaryFx <- maxRuns # may need to put arguments in the actual function call below!
				
				# fx <- function(x, th=th) { x >= th }
				# args <- list(th = 20)
				
				# pika <- extractAndCalcVarForOcc(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# priorYears=priorYears,
					# pika=pika,
					# summaryFx=summaryFx,
					# fail=fail,
					# failIfAllNA=TRUE,
					# fx=fx,
					# args=args
				# )	

			# ### subLethalHeat18deg_d: Longest run of days with mean (tmean) temperature ≥18°C, January-June
			# ###############################################################################################
				
				# # generalization
				# thisPred <- 'subLethalHeat18deg_d' # name of this variable
				# prVar <- 'tmean' # variable to extract from PRISM/TerraClimate
				# priorYears <- mostPriorYear:1 # years along window to extract
				
				# # function to summarize values within a year
				# summaryFx <- maxRuns # may need to put arguments in the actual function call below!
				
				# fx <- function(x, th=th) { x >= th }
				# args <- list(th = 18)
				
				# pika <- extractAndCalcVarForOcc(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# priorYears=priorYears,
					# pika=pika,
					# summaryFx=summaryFx,
					# fail=fail,
					# failIfAllNA=TRUE,
					# fx=fx,
					# args=args
				# )	

			# ### chronicHeat_C: Mean of summer average temperature (tmean), June-September
			# #############################################################################
				
				# # generalization
				# thisPred <- 'chronicHeat_C' # name of this variable
				# prVar <- 'tmean' # variable to extract from PRISM/TerraClimate
				# priorYears <- mostPriorYear:1 # years along window to extract
				
				# # function to summarize values within a year
				# summaryFx <- mean # may need to put arguments in the actual function call below!
				
				# pika <- extractAndCalcVarForOcc(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# priorYears=priorYears,
					# pika=pika,
					# summaryFx=summaryFx,
					# fail=fail,
					# na.rm=TRUE
				# )	

			# ### summerRespiteHeat_C: Mean of summer minimum temperature (tmin), June-September
			# ###################################################################################
				
				# # generalization
				# thisPred <- 'summerRespiteHeat_C' # name of this variable
				# prVar <- 'tmin' # variable to extract from PRISM/TerraClimate
				# priorYears <- mostPriorYear:1 # years along window to extract
				
				# # function to summarize values within a year
				# summaryFx <- mean # may need to put arguments in the actual function call below!
				
				# pika <- extractAndCalcVarForOcc(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# priorYears=priorYears,
					# pika=pika,
					# summaryFx=summaryFx,
					# fail=fail,
					# na.rm=TRUE
				# )	

			# ### chronicMoistStress_hPa: Mean of daily summer minimum vapor pressure deficit (vpdmin), June-September
			# ########################################################################################################
				
				# # generalization
				# thisPred <- 'chronicMoistStress_hPa' # name of this variable
				# prVar <- 'vpdmin' # variable to extract from PRISM/TerraClimate
				# priorYears <- mostPriorYear:1 # years along window to extract
				
				# # function to summarize values within a year
				# summaryFx <- mean # may need to put arguments in the actual function call below!
				
				# pika <- extractAndCalcVarForOcc(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# priorYears=priorYears,
					# pika=pika,
					# summaryFx=summaryFx,
					# fail=fail,
					# na.rm=TRUE
				# )	

			# ### gsVpdMin_hPa: Mean of daily summer minimum vapor pressure deficit (vpdmin), May-September
			# #############################################################################################
				
				# # generalization
				# thisPred <- 'gsVpdMin_hPa' # name of this variable
				# prVar <- 'vpdmin' # variable to extract from PRISM/TerraClimate
				# priorYears <- mostPriorYear:1 # years along window to extract
				
				# # function to summarize values within a year
				# summaryFx <- mean # may need to put arguments in the actual function call below!
				
				# pika <- extractAndCalcVarForOcc(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# priorYears=priorYears,
					# pika=pika,
					# summaryFx=summaryFx,
					# fail=fail,
					# na.rm=TRUE
				# )	

			# ### peakMoistStress_hPa: Mean of summer maximum vapor pressure deficit (vpdmax), May-September
			# ##############################################################################################
				
				# # generalization
				# thisPred <- 'peakMoistStress_hPa' # name of this variable
				# prVar <- 'vpdmax' # variable to extract from PRISM/TerraClimate
				# priorYears <- mostPriorYear:1 # years along window to extract
				
				# # function to summarize values within a year
				# summaryFx <- mean # may need to put arguments in the actual function call below!
				
				# pika <- extractAndCalcVarForOcc(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# priorYears=priorYears,
					# pika=pika,
					# summaryFx=summaryFx,
					# fail=fail,
					# na.rm=TRUE
				# )	

			# ### gsPpt_mm: Total growing-season precipitation (ppt), May-September
			# #####################################################################
				
				# # generalization
				# thisPred <- 'gsPpt_mm' # name of this variable
				# prVar <- 'ppt' # variable to extract from PRISM/TerraClimate
				# priorYears <- mostPriorYear:1 # years along window to extract
				
				# # function to summarize values within a year
				# summaryFx <- sum # may need to put arguments in the actual function call below!
				
				# pika <- extractAndCalcVarForOcc(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# priorYears=priorYears,
					# pika=pika,
					# summaryFx=summaryFx,
					# fail=fail,
					# na.rm=TRUE
				# )	

			# ### monsoonPpt_mm: Total monsoon precipitation (ppt), min-June-August
			# #####################################################################
				
				# # generalization
				# thisPred <- 'monsoonPpt_mm' # name of this variable
				# prVar <- 'ppt' # variable to extract from PRISM/TerraClimate
				# priorYears <- mostPriorYear:1 # years along window to extract
				
				# # function to summarize values within a year
				# summaryFx <- sum # may need to put arguments in the actual function call below!
				
				# pika <- extractAndCalcVarForOcc(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# priorYears=priorYears,
					# pika=pika,
					# summaryFx=summaryFx,
					# fail=fail,
					# na.rm=TRUE
				# )	

			# ### winterSnow_mm: Total precipitation, December-March
			# ######################################################
				
				# # generalization
				# thisPred <- 'winterSnow_mm' # name of this variable
				# prVar <- 'ppt' # variable to extract from PRISM/TerraClimate
				# priorYears <- mostPriorYear:0 # years along window to extract
				
				# # function to summarize values within a year
				# summaryFx <- sum # may need to put arguments in the actual function call below!
				
				# pika <- extractAndCalcVarForOcc(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# priorYears=priorYears,
					# pika=pika,
					# summaryFx=summaryFx,
					# fail=fail,
					# na.rm=TRUE
				# )

			# ### annualPpt_mm: Total annual precipitation, June-May
			# ######################################################
				
				# # generalization
				# thisPred <- 'annualPpt_mm' # name of this variable
				# prVar <- 'ppt' # variable to extract from PRISM/TerraClimate
				# priorYears <- mostPriorYear:1 # years along window to extract
				
				# # function to summarize values within a year
				# summaryFx <- sum # may need to put arguments in the actual function call below!
				
				# pika <- extractAndCalcVarForOcc(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# priorYears=priorYears,
					# pika=pika,
					# summaryFx=summaryFx,
					# fail=fail,
					# na.rm=TRUE
				# )	

			# ### springPpt_mm: Total spring precipitation (ppt), mid-March-May
			# #################################################################
				
				# # generalization
				# thisPred <- 'springPpt_mm' # name of this variable
				# prVar <- 'ppt' # variable to extract from PRISM/TerraClimate
				# priorYears <- mostPriorYear:0 # years along window to extract
				
				# # function to summarize values within a year
				# summaryFx <- sum # may need to put arguments in the actual function call below!
				
				# pika <- extractAndCalcVarForOcc(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# priorYears=priorYears,
					# pika=pika,
					# summaryFx=summaryFx,
					# fail=fail,
					# na.rm=TRUE
				# )
				
		# } # next most prior year

	# ### extract PRISM for DENSITY variables
	# #######################################

			# ### chronicCold_C: Average winter temperature (tmean), November-March
			# #####################################################################
				
				# # generalization
				# thisPred <- 'chronicCold_C' # name of this variable
				# prVar <- 'tmean' # variable to extract from PRISM/TerraClimate
				
				# # function to summarize values within a year
				# summaryFx <- mean # may need to put arguments in the actual function call below!

				# pika <- extractAndCalcVarForDens(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# pika=pika,
					# summaryFx=summaryFx,
					# fail=fail,
					# failIfAllNA=TRUE,
					# na.rm=TRUE
				# )	

			# ### subLethalCold_d: Length of longest run with winter minimum temperature (tmin) <5°C, September-March
			# #######################################################################################################
				
				# # generalization
				# thisPred <- 'subLethalCold_d' # name of this variable
				# prVar <- 'tmin' # variable to extract from PRISM/TerraClimate
				
				# # function to summarize values within a year
				# summaryFx <- maxRuns # may need to put arguments in the actual function call below!
				
				# fx <- function(x, th=th) { x < th }
				# args <- list(th = 5)
				
				# pika <- extractAndCalcVarForDens(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# pika=pika,
					# summaryFx=summaryFx,
					# fail=fail,
					# failIfAllNA=TRUE,
					# fx=fx,
					# args=args
				# )	

			# ### acuteCold_d: Number of winter days with minimum temperature (tmin) <-10°C, September-March
			# ##############################################################################################
				
				# # generalization
				# thisPred <- 'acuteCold_d' # name of this variable
				# prVar <- 'tmin' # variable to extract from PRISM/TerraClimate
				
				# # function to summarize values within a year
				# summaryFx <- totalDays # may need to put arguments in the actual function call below!
				
				# fx <- function(x, th=th) { x < th }
				# args <- list(th = -10)
				
				# pika <- extractAndCalcVarForDens(
					# thisPred=thisPred,
					# prVar=prVar,
					# predTable=predTable,
					# pika=pika,
					# summaryFx=summaryFx,
					# fail=fail,
					# failIfAllNA=TRUE,
					# fx=fx,
					# args=args
				# )	

		# ### acuteHeat_d: Number of days with maximum temperature ≥26°C, January-June
		# ############################################################################
			
			# # generalization
			# thisPred <- 'acuteHeat_d' # name of this variable
			# prVar <- 'tmax' # variable to extract from PRISM/TerraClimate
			
			# # function to summarize values within a year
			# summaryFx <- totalDays # may need to put arguments in the actual function call below!
			
			# fx <- function(x, th=th) { x >= th }
			# args <- list(th = 26)
			
			# pika <- extractAndCalcVarForDens(
				# thisPred=thisPred,
				# prVar=prVar,
				# predTable=predTable,
				# pika=pika,
				# summaryFx=summaryFx,
				# fail=fail,
				# failIfAllNA=TRUE,
				# fx=fx,
				# args=args
			# )	

		# ### subLethalHeat22deg_d: Longest run of days with mean (tmean) temperature ≥22°C, January-June
		# ###############################################################################################
			
			# # generalization
			# thisPred <- 'subLethalHeat22deg_d' # name of this variable
			# prVar <- 'tmean' # variable to extract from PRISM/TerraClimate
			
			# # function to summarize values within a year
			# summaryFx <- maxRuns # may need to put arguments in the actual function call below!
			
			# fx <- function(x, th=th) { x >= th }
			# args <- list(th = 22)
			
			# pika <- extractAndCalcVarForDens(
				# thisPred=thisPred,
				# prVar=prVar,
				# predTable=predTable,
				# pika=pika,
				# summaryFx=summaryFx,
				# fail=fail,
				# failIfAllNA=TRUE,
				# fx=fx,
				# args=args
			# )	

		# ### subLethalHeat20deg_d: Longest run of days with mean (tmean) temperature ≥20°C, January-June
		# ###############################################################################################
			
			# # generalization
			# thisPred <- 'subLethalHeat20deg_d' # name of this variable
			# prVar <- 'tmean' # variable to extract from PRISM/TerraClimate
			
			# # function to summarize values within a year
			# summaryFx <- maxRuns # may need to put arguments in the actual function call below!
			
			# fx <- function(x, th=th) { x >= th }
			# args <- list(th = 20)
			
			# pika <- extractAndCalcVarForDens(
				# thisPred=thisPred,
				# prVar=prVar,
				# predTable=predTable,
				# pika=pika,
				# summaryFx=summaryFx,
				# fail=fail,
				# failIfAllNA=TRUE,
				# fx=fx,
				# args=args
			# )	

		# ### subLethalHeat18deg_d: Longest run of days with mean (tmean) temperature ≥18°C, January-June
		# ###############################################################################################
			
			# # generalization
			# thisPred <- 'subLethalHeat18deg_d' # name of this variable
			# prVar <- 'tmean' # variable to extract from PRISM/TerraClimate
			
			# # function to summarize values within a year
			# summaryFx <- maxRuns # may need to put arguments in the actual function call below!
			
			# fx <- function(x, th=th) { x >= th }
			# args <- list(th = 18)
			
			# pika <- extractAndCalcVarForDens(
				# thisPred=thisPred,
				# prVar=prVar,
				# predTable=predTable,
				# pika=pika,
				# summaryFx=summaryFx,
				# fail=fail,
				# failIfAllNA=TRUE,
				# fx=fx,
				# args=args
			# )	

		# ### chronicHeat_C: Mean of summer average temperature (tmean), June-September
		# #############################################################################
			
			# # generalization
			# thisPred <- 'chronicHeat_C' # name of this variable
			# prVar <- 'tmean' # variable to extract from PRISM/TerraClimate
			
			# # function to summarize values within a year
			# summaryFx <- mean # may need to put arguments in the actual function call below!
			
			# pika <- extractAndCalcVarForDens(
				# thisPred=thisPred,
				# prVar=prVar,
				# predTable=predTable,
				# pika=pika,
				# summaryFx=summaryFx,
				# fail=fail,
				# failIfAllNA=TRUE,
				# na.rm=TRUE
			# )	

		# ### summerRespiteHeat_C: Mean of summer minimum temperature (tmin), June-September
		# ###################################################################################
			
			# # generalization
			# thisPred <- 'summerRespiteHeat_C' # name of this variable
			# prVar <- 'tmin' # variable to extract from PRISM/TerraClimate
			
			# # function to summarize values within a year
			# summaryFx <- mean # may need to put arguments in the actual function call below!
			
			# pika <- extractAndCalcVarForDens(
				# thisPred=thisPred,
				# prVar=prVar,
				# predTable=predTable,
				# pika=pika,
				# summaryFx=summaryFx,
				# fail=fail,
				# failIfAllNA=TRUE,
				# na.rm=TRUE
			# )	

		# ### chronicMoistStress_hPa: Mean of daily summer minimum vapor pressure deficit (vpdmin), June-September
		# ########################################################################################################
			
			# # generalization
			# thisPred <- 'chronicMoistStress_hPa' # name of this variable
			# prVar <- 'vpdmin' # variable to extract from PRISM/TerraClimate
			
			# # function to summarize values within a year
			# summaryFx <- mean # may need to put arguments in the actual function call below!
			
			# pika <- extractAndCalcVarForDens(
				# thisPred=thisPred,
				# prVar=prVar,
				# predTable=predTable,
				# pika=pika,
				# summaryFx=summaryFx,
				# fail=fail,
				# failIfAllNA=TRUE,
				# na.rm=TRUE
			# )	

		# ### gsVpdMin_hPa: Mean of daily summer minimum vapor pressure deficit (vpdmin), May-September
		# #############################################################################################
			
			# # generalization
			# thisPred <- 'gsVpdMin_hPa' # name of this variable
			# prVar <- 'vpdmin' # variable to extract from PRISM/TerraClimate
			
			# # function to summarize values within a year
			# summaryFx <- mean # may need to put arguments in the actual function call below!
			
			# pika <- extractAndCalcVarForDens(
				# thisPred=thisPred,
				# prVar=prVar,
				# predTable=predTable,
				# pika=pika,
				# summaryFx=summaryFx,
				# fail=fail,
				# failIfAllNA=TRUE,
				# na.rm=TRUE
			# )	

		# ### peakMoistStress_hPa: Mean of summer maximum vapor pressure deficit (vpdmax), May-September
		# ##############################################################################################
			
			# # generalization
			# thisPred <- 'peakMoistStress_hPa' # name of this variable
			# prVar <- 'vpdmax' # variable to extract from PRISM/TerraClimate
			
			# # function to summarize values within a year
			# summaryFx <- mean # may need to put arguments in the actual function call below!
			
			# pika <- extractAndCalcVarForDens(
				# thisPred=thisPred,
				# prVar=prVar,
				# predTable=predTable,
				# pika=pika,
				# summaryFx=summaryFx,
				# fail=fail,
				# failIfAllNA=TRUE,
				# na.rm=TRUE
			# )	

		# ### gsPpt_mm: Total growing-season precipitation (ppt), May-September
		# #####################################################################
			
			# # generalization
			# thisPred <- 'gsPpt_mm' # name of this variable
			# prVar <- 'ppt' # variable to extract from PRISM/TerraClimate
			
			# # function to summarize values within a year
			# summaryFx <- sum # may need to put arguments in the actual function call below!
			
			# pika <- extractAndCalcVarForDens(
				# thisPred=thisPred,
				# prVar=prVar,
				# predTable=predTable,
				# pika=pika,
				# summaryFx=summaryFx,
				# fail=fail,
				# failIfAllNA=TRUE,
				# na.rm=TRUE
			# )	

		# ### monsoonPpt_mm: Total monsoon precipitation (ppt), min-June-August
		# #####################################################################
			
			# # generalization
			# thisPred <- 'monsoonPpt_mm' # name of this variable
			# prVar <- 'ppt' # variable to extract from PRISM/TerraClimate
			
			# # function to summarize values within a year
			# summaryFx <- sum # may need to put arguments in the actual function call below!
			
			# pika <- extractAndCalcVarForDens(
				# thisPred=thisPred,
				# prVar=prVar,
				# predTable=predTable,
				# pika=pika,
				# summaryFx=summaryFx,
				# fail=fail,
				# na.rm=TRUE
			# )	
		
		# ### winterSnow_mm: Total precipitation, December-March
		# ######################################################
			
			# # generalization
			# thisPred <- 'winterSnow_mm' # name of this variable
			# prVar <- 'ppt' # variable to extract from PRISM/TerraClimate
			
			# # function to summarize values within a year
			# summaryFx <- sum # may need to put arguments in the actual function call below!
			
			# pika <- extractAndCalcVarForDens(
				# thisPred=thisPred,
				# prVar=prVar,
				# predTable=predTable,
				# pika=pika,
				# summaryFx=summaryFx,
				# fail=fail,
				# na.rm=TRUE
			# )

		# ### annualPpt_mm: Total annual precipitation, June-May
		# ######################################################
			
			# # generalization
			# thisPred <- 'annualPpt_mm' # name of this variable
			# prVar <- 'ppt' # variable to extract from PRISM/TerraClimate
			
			# # function to summarize values within a year
			# summaryFx <- sum # may need to put arguments in the actual function call below!
			
			# pika <- extractAndCalcVarForDens(
				# thisPred=thisPred,
				# prVar=prVar,
				# predTable=predTable,
				# pika=pika,
				# summaryFx=summaryFx,
				# fail=fail,
				# na.rm=TRUE
			# )	

		# ### springPpt_mm: Total spring precipitation (ppt), mid-March-May
		# #################################################################
			
			# # generalization
			# thisPred <- 'springPpt_mm' # name of this variable
			# prVar <- 'ppt' # variable to extract from PRISM/TerraClimate
			
			# # function to summarize values within a year
			# summaryFx <- sum # may need to put arguments in the actual function call below!
			
			# pika <- extractAndCalcVarForDens(
				# thisPred=thisPred,
				# prVar=prVar,
				# predTable=predTable,
				# pika=pika,
				# summaryFx=summaryFx,
				# fail=fail,
				# failIfAllNA=TRUE,
				# na.rm=TRUE
			# )	

	# save(pika, file='./Data/02 New Mexico Pika - Environmental Values Extracted and Calculated.rda')
	
# say('###################################')
# say('### distributions of predictors ###')
# say('###################################')

	# ### generalization
	# thold <- 0.7 # too much correlation!

	# ### data
	# load('./Data/02 New Mexico Pika - Environmental Values Extracted and Calculated.rda')

	# ### occupancy predictors
	# ########################
	
		# pika$latestOccStatus <- factor(pika$latestOccStatus, levels=c('0 never', '1 old', '2 occupied'))
		
		# preds <- predTable$var[predTable$useOccAbund]

		# for (occWindow in occWindows_y) {
			
			# figs <- list()
			# for (countPred in seq_along(preds)) {
			
				# pred <- preds[countPred]
				
				# predIndex <- which(predTable$var == pred)
			
				# predNice <- predTable$varNice[predIndex]
				# predDescriptorUnit <- paste0(predTable$unitDescriptor[predIndex], ' (', predTable$unit[predIndex], ')')
			
				# predWindow <- paste0('occVar_', pred, '_', occWindow, 'yrWindow')
			
				# mus <- data.frame(
					# latestOccStatus = c('0 never', '1 old', '2 occupied'),
					# mu = c(
						# mean(pika[pika$latestOccStatus == '0 never', predWindow]),
						# mean(pika[pika$latestOccStatus == '1 old', predWindow]),
						# mean(pika[pika$latestOccStatus == '2 occupied', predWindow])
					# )
				# )

				# say(predWindow)


				# thisData <- pika[ , c('latestOccStatus', predWindow)]
				# names(thisData)[2] <- 'value'
				
				# figs[[countPred]] <- ggplot(data=thisData, aes(x=value, col=latestOccStatus, fill=latestOccStatus)) +
					# geom_density(size=1) +
					# scale_color_manual(values=c('0 never'='firebrick3', '1 old'='darkgoldenrod3', '2 occupied'='darkgreen')) +
					# scale_fill_manual(values=alpha(c('0 never'='firebrick3', '1 old'='darkgoldenrod3', '2 occupied'='darkgreen'), 0.2)) +
					# labs(title=predNice, subtitle=pred, x=predDescriptorUnit, y=NULL) +
					# geom_vline(data=mus, aes(xintercept=mu, color=latestOccStatus), linetype='dotted', size=1) +
					# theme(legend.position='none', plot.title=element_text(face='bold'))
					
				# figs[[countPred]]
					
			# } # next predictor

			# main <- plot_grid(plotlist=figs, align='h', ncol=3, rel_widths=1)
		
			# main
		
			# ggsave(paste0('./Figures & Tables/Distributions of Occupancy Variables for ', occWindow, '-yr Window.png'), width=8.5, height=11, units='in')

		# } # next occupancy window

	# ### density predictors
	# ######################
	
		# preds <- names(pika)[grepl(names(pika), pattern='densVar_')]

		# figs <- list()
		# for (countPred in seq_along(preds)) {
		
			# predFull <- preds[countPred]
			# timePeriod <- if (grepl(predFull, pattern='1yrPrior')) { 'prior year' } else { 'same year'}
			
			# say(predFull)
			
			# pred <- gsub(predFull, pattern='densVar_', replacement='')
			# pred <- gsub(pred, pattern='_1yrPrior', replacement='')
			# pred <- gsub(pred, pattern='_0yrPrior', replacement='')
			
			# predIndex <- which(predTable$var == pred)
		
			# predNice <- predTable$varNice[predIndex]
			# predDescriptorUnit <- paste0(predTable$unitDescriptor[predIndex], ' (', predTable$unit[predIndex], ')')
		
			# thisData <- pika[ , c('latestDensity', predFull)]
			# names(thisData)[2] <- 'value'
			# thisData <- thisData[complete.cases(thisData), ]

			# subtitle <- paste0(timePeriod, ' (', pred, ')')
				
			# figs[[countPred]] <- ggplot(data=thisData, aes(x=value, y=latestDensity)) +
				# geom_point(size=2, shape=16, col='gray40') +
				# labs(title=predNice, subtitle=subtitle, x=predDescriptorUnit, y='Abundance') +
				# theme(
					# legend.position='none',
					# plot.title=element_text(face='bold', size=8),
					# plot.subtitle=element_text(size=6),
					# axis.text.x=element_text(size=7),
					# axis.text.y=element_text(size=7)
				# )
				
		# } # next predictor

		# main <- plot_grid(plotlist=figs, align='h', ncol=6, rel_widths=1)
	
		# ggsave(plot=main, paste0('./Figures & Tables/Distributions of Density Variables.png'), width=11, height=8.5, units='in')

	# ### correlations for occupancy variables: heatmap
	# #################################################
	
		# for (occWindow in occWindows_y) {
		
			# occVars <- names(pika)[grepl(names(pika), pattern='occVar_')]
			# occVars <- occVars[grepl(occVars, pattern=paste0(occWindow, 'yrWindow'))]

			# corr <- cor(pika[ , occVars])
			
			# # change variable names
			# rownames(corr) <- colnames(corr) <- makeNiceVars(colnames(corr), occOrDens='occupancy')

			# title <- paste0('Occupancy variables with a ', occWindow, '-yr window')
			# fig <- ggcorrplot(corr, hc.order=TRUE, outline.color='white', lab=TRUE, title=title, legend.title='Correlation', show.diag=FALSE)
			
			# ggsave(plot=fig, paste0('./Figures & Tables/Correlations between Occupancy Variables Using a ', occWindow, '-yr Window Heat Map.png'), width=10, height=10, units='in')
		
		# }
	
	# ### correlations for occupancy variables: spoke plot
	# ####################################################

	# for (occWindow in occWindows_y) {
	
		# occVars <- names(pika)[grepl(names(pika), pattern='occVar_')]
		# occVars <- occVars[grepl(occVars, pattern=paste0(occWindow, 'yrWindow'))]

		# corr <- cor(pika[ , occVars])
		
		# # change variable names
		# niceVars <- makeNiceVars(colnames(corr), occOrDens='occupancy', incTime=FALSE)

		# title <- paste0('Occupancy variables with a ', occWindow, '-yr window')
		# png(paste0('./Figures & Tables/Correlations between Occupancy Variables Using a ', occWindow, '-yr Window Spoke Plot.png'), width=1000, height=1000)

			# par(oma=c(3, 10, 3, 3))
			# spoke(pos=corr > thold, neg=corr < -thold, labels=niceVars, ltyNeg='solid', colNeg='red', main=title, cexLabel=1.6, cex.main=1.9)
			# legend('bottomright', inset=-0.01, legend=c(paste0('Positive >', thold), paste0('Negative <-', thold)), lwd=1, col=c('black', 'red'), bty='n', cex=1.8)
			
		# dev.off()

	# }


	# ### correlations for density variables: heatmap
	# ###############################################
	
		# densVars <- names(pika)[grepl(names(pika), pattern='densVar_')]

		# recordedDens <- which(!is.na(pika$latestDensSurveyYear))
		# corr <- cor(pika[recordedDens, densVars])
		
		# # change variable names
		# vars <- colnames(corr)
		# vars <- gsub(vars, pattern='densVar_', replacement='')
		
		# times <- rep(NA, length(vars))
		# times[grepl(vars, pattern=paste0('_1yrPrior'))] <- '(1 yr prior)'
		# times[grepl(vars, pattern=paste0('_0yrPrior'))] <- '(same year)'
		
		# vars <- gsub(vars, pattern=paste0('_1yrPrior'), replacement='')
		# vars <- gsub(vars, pattern=paste0('_0yrPrior'), replacement='')
		# preds <- predTable$varNice[match(vars, predTable$var)]
		
		# preds <- paste(preds, times)
		# colnames(corr) <- rownames(corr) <- preds

		# title <- paste0('Density variables')
		# fig <- ggcorrplot(corr, hc.order=TRUE, outline.color='white', lab=TRUE, title=title, legend.title='Correlation', show.diag=FALSE, lab_size=2.5)
			
		# ggsave(plot=fig, paste0('./Figures & Tables/Correlations between Density Variables Heat Map.png'), width=11, height=8.5, units='in')
		
	# ### correlations for density variables: spoke plot
	# ####################################################
	
		# densVars <- names(pika)[grepl(names(pika), pattern='densVar_')]
		# recordedDens <- which(!is.na(pika$latestDensSurveyYear))
		# corr <- cor(pika[recordedDens, densVars])
		
		# # change variable names
		# niceVars <- makeNiceVars(colnames(corr), occOrDens='density', incTime=TRUE, wrapTime=TRUE)

		# title <- paste0('Density variables')
		# png(paste0('./Figures & Tables/Correlations between Density Variables Spoke Plot.png'), width=1000, height=1000)

			# par(oma=c(3, 10, 3, 3))
			# spoke(pos=corr > thold, neg=corr < -thold, labels=niceVars, ltyNeg='solid', colNeg='red', main=title, cexLabel=1.6, cex.main=1.9)
			# legend('bottomright', inset=-0.01, legend=c(paste0('Positive >', thold), paste0('Negative <-', thold)), lwd=1, col=c('black', 'red'), bty='n', cex=1.8)
			
		# dev.off()

# say('################################')
# say('### define regions and folds ###')
# say('################################')

	# load('./Data/02 New Mexico Pika - Environmental Values Extracted and Calculated.rda')

	# # define regions
	# pika$region <- NA
	# pika$region[pika$longitude < -105.9448 & pika$latitude > 36.29301] <- 'northwest'
	# pika$region[pika$longitude < -105.9448 & pika$latitude < 36.29301] <- 'southwest'
	# pika$region[pika$longitude > -105.9448 & pika$latitude > 36.29301] <- 'northeast'
	# pika$region[pika$longitude > -105.9448 & pika$latitude < 36.29301] <- 'southeast'
		
	# # define folds
	# pika$fold <- NA
	# for (region in regions) {
	
		# lats <- pika$latitude[pika$region == region]
		# mid <- mean(lats)
		# pika$fold[pika$region == region & pika$latitude > mid] <- 1
		# pika$fold[pika$region == region & pika$latitude <= mid] <- 2
		# if (region == 'northwest') {

			# pika$fold[pika$region == region & ((pika$longitude > -106.1828 & pika$latitude > 36.58048) | pika$latitude > mid)] <- 1
		
		# }
	
	# }

	# # visual check
	# sp <- vect(as.matrix(pika[ , ll]), att=pika, crs=getCRS('wgs84'))
	# plot(sp, col=NA)
	# points(sp[sp$region == 'northwest' & sp$fold == 1, ], col='red', pch=1)
	# points(sp[sp$region == 'northwest' & sp$fold == 2, ], col='red', pch=2)

	# points(sp[sp$region == 'southwest' & sp$fold == 1, ], col='green', pch=1)
	# points(sp[sp$region == 'southwest' & sp$fold == 2, ], col='green', pch=2)

	# points(sp[sp$region == 'northeast' & sp$fold == 1, ], col='blue', pch=1)
	# points(sp[sp$region == 'northeast' & sp$fold == 2, ], col='blue', pch=2)

	# points(sp[sp$region == 'southeast' & sp$fold == 1, ], col='orange', pch=1)
	# points(sp[sp$region == 'southeast' & sp$fold == 2, ], col='orange', pch=2)

	# # sample size checks
	# for (region in regions) {

		# say(region, ' ============================')
		# for (fold in 1:2) {

			# for (occStatus in occStatuses) {

				# n <- sum(pika$latestOccStatus[pika$region == region & pika$fold == fold]==occStatus)
				# say('fold ', fold, ': ', occStatus, ' ', n)
			# }
			
		# }
		
	# }

	# # presence/absence
	# pika$presAbs <- 0
	# pika$presAbs[pika$latestOccStatus == '2 occupied'] <- 1

	# save(pika, file='./Data/03 New Mexico Pika - Assigned Folds.rda')
	
say('DONE!!!', level=1, deco='%')
