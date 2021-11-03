### NEW MEXICO PIKA -- Extract climate data for Bandelier Sites
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('C:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/XX New Mexico Pika Occupancy & Abundance Analysis - Extract Data for Bandelier.r')
### source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/XX New Mexico Pika Occupancy & Abundance Analysis - Extract Data for Bandelier.r')
###
### This script extracts climate data for select sites at Bandalier. Data is extracted for the 7 yr prior to sampling *each* survey at each site.
###
### CONTENTS ###
### process and clean pika data ###
### extract environmental data and calculate predictors ###

#############
### setup ###
#############

	# source('C:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r')

	source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r')

say('###################################')
say('### process and clean pika data ###')
say('###################################')

	# reformat Excel data into R
	file <- paste0('./Data/00 NM Database 2021 Master [Downloaded 2021-04-23 ORIGINAL].xlsx')
	pika <- read_excel(file, sheet='All Patches', skip=4, col_types='text', na=c('N/A', 'NA', ''))

	pika <- as.data.frame(pika)

	### nice names
	##############

		names(pika)[names(pika) == 'Polygon Name'] <- 'polygonName'
		names(pika)[names(pika) == 'Latitude'] <- 'latitude'
		names(pika)[names(pika) == 'Longitude'] <- 'longitude'
		names(pika)[names(pika) == 'Elevation (m)'] <- 'elevation_m'
		names(pika)[names(pika) == 'Elevation (ft)'] <- 'elevation_ft'
		names(pika)[names(pika) == 'Mean Aspect'] <- 'aspectMean_deg'
		names(pika)[names(pika) == 'Mean Slope'] <- 'slopeMean_deg'
		names(pika)[names(pika) == 'Size (#HR)'] <- 'numHomeRanges'
		names(pika)[names(pika) == 'Color Talus (most recent survey):'] <- 'origOccCodeMostRecent'
		names(pika)[names(pika) == '2016 Survey Results (Color)'] <- 'origOccCode2016'
		names(pika)[names(pika) == '2017 Survey Results (Color)'] <- 'origOccCode2017'
		names(pika)[names(pika) == '2018 Survey Results (Color)'] <- 'origOccCode2018'
		names(pika)[names(pika) == '2019 Survey Results (Color)'] <- 'origOccCode2019'
		names(pika)[names(pika) == '2020 Suvey Results (Color)'] <- 'origOccCode2020'
		names(pika)[names(pika) == 'Change in Occupancy'] <- 'changeInOcc'
		names(pika)[names(pika) == 'Patch Notes'] <- 'patchNotes'
		pika$`...17` <- NULL
		names(pika)[names(pika) == '2016 Survey Results (Density)'] <- 'origDensity2016'
		names(pika)[names(pika) == '2017 Survey Results (Density)'] <- 'origDensity2017'
		names(pika)[names(pika) == '2018 Survey Results (Density)'] <- 'origDensity2018'
		names(pika)[names(pika) == '2019 Survey Results (Density)'] <- 'origDensity2019'
		names(pika)[names(pika) == '2020 Survey Results (Density)'] <- 'origDensity2020'
		pika$`...23` <- NULL
		names(pika)[names(pika) == 'Grazing Status (Most Recent Survey; Yes=1, No=0)'] <- 'grazingMostRecent'
		names(pika)[names(pika) == '2016 Grazing Status (Yes=1, No=0)'] <- 'grazing2016'
		names(pika)[names(pika) == '2017 Grazing Status (Yes=1, No=0)'] <- 'grazing2017'
		names(pika)[names(pika) == '2018 Grazing Status (Yes=1, No=0)'] <- 'grazing2018'
		names(pika)[names(pika) == '2019 Grazing Status (Yes=1, No=0)'] <- 'grazing2019'
		names(pika)[names(pika) == '2020 Grazing Status (Yes=1, No=0)'] <- 'grazing2020'
		pika$`...30` <- NULL
		names(pika)[names(pika) == 'Percent Patch Perimeter Burned       (Most Recent Survey)'] <- 'perimBurnedMostRecent_perc'
		names(pika)[names(pika) == '2016 Percent Patch Perimeter Burned'] <- 'perimBurned2016_perc'
		names(pika)[names(pika) == '2017 Percent Patch Perimeter Burned'] <- 'perimBurned2017_perc'
		names(pika)[names(pika) == '2018 Percent Patch Perimeter Burned'] <- 'perimBurned2018_perc'
		names(pika)[names(pika) == '2019 Percent Patch Perimeter Burned'] <- 'perimBurned2019_perc'
		names(pika)[names(pika) == '2020 Percent Patch Perimeter Burned'] <- 'perimBurned2020_perc'
		names(pika)[names(pika) == 'Year of Most Recent Fire (Prior to 2017)'] <- 'mostRecentFireYear'
		names(pika)[names(pika) == 'Most Recent Fire Name (Prior to 2017)'] <- 'mostRecentFireName'
		names(pika)[names(pika) == '# Acres Buned in Most Recent Fire (Prior to 2017)'] <- 'mostRecentFireArea_ac'
		pika$`...40` <- NULL
		names(pika)[names(pika) == '2016 Patch-Wide Mean Grass Cover, all spp. (%)'] <- 'grass2016_perc'
		names(pika)[names(pika) == '2016 Patch-Wide Mean Forb Cover, all spp. (%)'] <- 'forb2016_perc'
		names(pika)[names(pika) == '2016 Patch-Wide Mean Total G+F, all spp. (%)'] <- 'grassForb2016_perc'
		names(pika)[names(pika) == '2017 Patch-Wide Mean Grass Cover, all spp. (%)'] <- 'grass2017_perc'
		names(pika)[names(pika) == '2017 Patch-Wide Mean Forb Cover, all spp. (%)'] <- 'forb2017_perc'
		names(pika)[names(pika) == '2017 Patch-Wide Mean Total G+F, all spp. (%)'] <- 'grassForb2017_perc'
		names(pika)[names(pika) == '2018 Patch-Wide Mean Grass Cover, all spp. (%)'] <- 'grass2018_perc'
		names(pika)[names(pika) == '2018 Patch-Wide Mean Forb Cover, all spp. (%)'] <- 'forb2018_perc'
		names(pika)[names(pika) == '2018 Patch-Wide Mean Total G+F, all spp. (%)'] <- 'grassForb2018_perc'
		names(pika)[names(pika) == '2019 Patch-Wide Mean Grass Cover, all spp. (%)'] <- 'grass2019_perc'
		names(pika)[names(pika) == '2019 Patch-Wide Mean Forb Cover, all spp. (%)'] <- 'forb2019_perc'
		names(pika)[names(pika) == '2019 Patch-Wide Mean Total G+F, all spp. (%)'] <- 'grassForb2019_perc'
		names(pika)[names(pika) == '2020 Patch-Wide Mean Grass Cover, all spp. (%)'] <- 'grass2020_perc'
		names(pika)[names(pika) == '2020 Patch-Wide Mean Forb Cover, all spp. (%)'] <- 'forb2020_perc'
		names(pika)[names(pika) == '2020 Patch-Wide Mean Total G+F, all spp. (%)'] <- 'grassForb2020_perc'
		pika$`...56` <- NULL
		names(pika)[names(pika) == '2016 Survey Date'] <- 'surveyDate2016'
		names(pika)[names(pika) == '2017 Survey Date'] <- 'surveyDate2017'
		names(pika)[names(pika) == '2018 Survey Date'] <- 'surveyDate2018'
		names(pika)[names(pika) == '2019 Survey Date'] <- 'surveyDate2019'
		names(pika)[names(pika) == '2020 Survey Date'] <- 'surveyDate2020'
		pika$`...62` <- NULL
		names(pika)[names(pika) == 'Area'] <- 'areaName'
		names(pika)[names(pika) == 'Best Access Notes'] <- 'bestAccessNotes'
		names(pika)[names(pika) == 'Added By'] <- 'addedBy'
		names(pika)[names(pika) == 'Historical Site'] <- 'historicalSite'
		pika$`...67` <- NULL

	### remove all "problem sites"
	##############################
	
		probsStart <- which(pika$polygonName == 'PROBLEM SITES')
		pika <- pika[1:(probsStart - 1), ]
		
		empty <- which(is.na(pika$polygonName))
		if (length(empty) > 0) pika <- pika[-empty, ]

		if (any(pika$polygonName == 'Misattributed Survey: Not Pecos5')) pika <- pika[-which(pika$polygonName == 'Misattributed Survey: Not Pecos5'), ]
		
	### column-specific changes and data corrections
	################################################
	
		pika$longitude <- as.numeric(pika$longitude)
		pika$latitude <- as.numeric(pika$latitude)
		pika$elevation_m <- as.numeric(pika$elevation_m)

	### select selected Bandelier sites
	###################################
	
		# NB excluding Band10 and BandUnsearched?
		bands <- c('Band 101', 'Band 18', 'Band 19', 'Band Fire', 'Band103', 'BAND11', 'BAND12', 'BAND13', 'BAND15', 'Band16', 'Band17', 'BAND3', 'BAND4', 'BAND5', 'BAND6', 'BAND8', 'BAND9', 'BANDnew', 'BandUnsearched2')

		pika <- pika[pika$polygonName %in% bands, ]

	### checks
	##########
	
		say('Number of unique polygon names is same as number of rows: ', length(unique(pika$polygonName)) == nrow(pika))
		say('Number of records with NA for longitude ........ ', sum(is.na(pika$longitude)))
		say('Number of records with NA for latitude ......... ', sum(is.na(pika$latitude)))
		say('Number of records with NA for elevation_m ...... ', sum(is.na(pika$origOccCodeMostRecent)))

say('###########################################################')
say('### extract environmental data and calculate predictors ###')
say('###########################################################')

	say('This step extracts environmental data from PRISM and calculates the derived climate variables as described in the Excel document created by Erik and Maria.', breaks=80)

	### generalization
	# fail <- FALSE # if TRUE then fail if predictor has any infinite or NA values!
	fail <- TRUE # if TRUE then fail if predictor has any infinite or NA values!
	
	say('Fail on NA/infinite value: ', fail)

	# PRISM daily base directory
	prDir <- 'I:/Ecology/Climate/PRISM/working/an81'
	
	# load airUpThere functions for extracting climate data
	fxs <- listFiles(paste0(drive, '/Ecology/Drive/R/airUpThere/R'), pattern='.r')
	for (fx in fxs) source(fx)
	
	fxs <- listFiles(paste0(drive, '/Ecology/Drive/R/airUpThere/data'), pattern='.rda')
	for (fx in fxs) load(fx)
	
	### extract PRISM for OCCCUPANCY variables
	##########################################

		for (surveyYear in surveyYears) {
		
			say(surveyYear, level=2, deco='*')
		
			pika$dummyDate <- as.Date(paste0(surveyYear, '-06-06'))

			for (mostPriorYear in occWindows_y) {

				### chronicCold_C: Average winter temperature (tmean), November-March
				#####################################################################
					
					# generalization
					thisPred <- 'chronicCold_C' # name of this variable
					prVar <- 'tmean' # variable to extract from PRISM/TerraClimate
					priorYears <- mostPriorYear:0 # years along window to extract
					
					# function to summarize values within a year
					summaryFx <- mean # may need to put arguments in the actual function call below!

					pika <- extractAndCalcVarForOcc(
						thisPred=thisPred,
						prVar=prVar,
						predTable=predTable,
						priorYears=priorYears,
						pika=pika,
						args=args,
						dateField='dummyDate',
						fail=fail,
						summaryFx=summaryFx,
						na.rm=TRUE
					)
					
					names(pika)[ncol(pika)] <- paste0(names(pika)[ncol(pika)], '_end', surveyYear)

				### subLethalCold_d: Length of longest run with winter minimum temperature (tmin) <5°C, September-March
				#######################################################################################################
					
					# generalization
					thisPred <- 'subLethalCold_d' # name of this variable
					prVar <- 'tmin' # variable to extract from PRISM/TerraClimate
					priorYears <- mostPriorYear:0 # years along window to extract
					
					# function to summarize values within a year
					summaryFx <- maxRuns # may need to put arguments in the actual function call below!
					
					fx <- function(x, th=th) { x < th }
					args <- list(th = 5)
					
					pika <- extractAndCalcVarForOcc(
						thisPred=thisPred,
						prVar=prVar,
						predTable=predTable,
						priorYears=priorYears,
						pika=pika,
						summaryFx=summaryFx,
						dateField='dummyDate',
						fail=fail,
						failIfAllNA=TRUE,
						fx=fx,
						args=args
					)	

					names(pika)[ncol(pika)] <- paste0(names(pika)[ncol(pika)], '_end', surveyYear)

				### acuteCold_d: Number of winter days with minimum temperature (tmin) <-10°C, September-March
				##############################################################################################
					
					# generalization
					thisPred <- 'acuteCold_d' # name of this variable
					prVar <- 'tmin' # variable to extract from PRISM/TerraClimate
					priorYears <- mostPriorYear:0 # years along window to extract
					
					# function to summarize values within a year
					summaryFx <- totalDays # may need to put arguments in the actual function call below!
					
					fx <- function(x, th=th) { x < th }
					args <- list(th = -10)
					
					pika <- extractAndCalcVarForOcc(
						thisPred=thisPred,
						prVar=prVar,
						predTable=predTable,
						priorYears=priorYears,
						pika=pika,
						summaryFx=summaryFx,
						dateField='dummyDate',
						fail=fail,
						failIfAllNA=TRUE,
						fx=fx,
						args=args
					)	

					names(pika)[ncol(pika)] <- paste0(names(pika)[ncol(pika)], '_end', surveyYear)

				### acuteHeat_d: Number of days with maximum temperature ≥26°C, January-June
				############################################################################
					
					# generalization
					thisPred <- 'acuteHeat_d' # name of this variable
					prVar <- 'tmax' # variable to extract from PRISM/TerraClimate
					priorYears <- mostPriorYear:1 # years along window to extract
					
					# function to summarize values within a year
					summaryFx <- totalDays # may need to put arguments in the actual function call below!
					
					fx <- function(x, th=th) { x >= th }
					args <- list(th = 26)
					
					pika <- extractAndCalcVarForOcc(
						thisPred=thisPred,
						prVar=prVar,
						predTable=predTable,
						priorYears=priorYears,
						pika=pika,
						summaryFx=summaryFx,
						dateField='dummyDate',
						fail=fail,
						failIfAllNA=TRUE,
						fx=fx,
						args=args
					)	

					names(pika)[ncol(pika)] <- paste0(names(pika)[ncol(pika)], '_end', surveyYear)

				### subLethalHeat22deg_d: Longest run of days with mean (tmean) temperature ≥22°C, January-June
				###############################################################################################
					
					# generalization
					thisPred <- 'subLethalHeat22deg_d' # name of this variable
					prVar <- 'tmean' # variable to extract from PRISM/TerraClimate
					priorYears <- mostPriorYear:1 # years along window to extract
					
					# function to summarize values within a year
					summaryFx <- maxRuns # may need to put arguments in the actual function call below!
					
					fx <- function(x, th=th) { x >= th }
					args <- list(th = 22)
					
					pika <- extractAndCalcVarForOcc(
						thisPred=thisPred,
						prVar=prVar,
						predTable=predTable,
						priorYears=priorYears,
						pika=pika,
						summaryFx=summaryFx,
						dateField='dummyDate',
						fail=fail,
						failIfAllNA=TRUE,
						fx=fx,
						args=args
					)	

					names(pika)[ncol(pika)] <- paste0(names(pika)[ncol(pika)], '_end', surveyYear)

				### subLethalHeat20deg_d: Longest run of days with mean (tmean) temperature ≥20°C, January-June
				###############################################################################################
					
					# generalization
					thisPred <- 'subLethalHeat20deg_d' # name of this variable
					prVar <- 'tmean' # variable to extract from PRISM/TerraClimate
					priorYears <- mostPriorYear:1 # years along window to extract
					
					# function to summarize values within a year
					summaryFx <- maxRuns # may need to put arguments in the actual function call below!
					
					fx <- function(x, th=th) { x >= th }
					args <- list(th = 20)
					
					pika <- extractAndCalcVarForOcc(
						thisPred=thisPred,
						prVar=prVar,
						predTable=predTable,
						priorYears=priorYears,
						pika=pika,
						summaryFx=summaryFx,
						dateField='dummyDate',
						fail=fail,
						failIfAllNA=TRUE,
						fx=fx,
						args=args
					)	

					names(pika)[ncol(pika)] <- paste0(names(pika)[ncol(pika)], '_end', surveyYear)

				### subLethalHeat18deg_d: Longest run of days with mean (tmean) temperature ≥18°C, January-June
				###############################################################################################
					
					# generalization
					thisPred <- 'subLethalHeat18deg_d' # name of this variable
					prVar <- 'tmean' # variable to extract from PRISM/TerraClimate
					priorYears <- mostPriorYear:1 # years along window to extract
					
					# function to summarize values within a year
					summaryFx <- maxRuns # may need to put arguments in the actual function call below!
					
					fx <- function(x, th=th) { x >= th }
					args <- list(th = 18)
					
					pika <- extractAndCalcVarForOcc(
						thisPred=thisPred,
						prVar=prVar,
						predTable=predTable,
						priorYears=priorYears,
						pika=pika,
						summaryFx=summaryFx,
						dateField='dummyDate',
						fail=fail,
						failIfAllNA=TRUE,
						fx=fx,
						args=args
					)	

					names(pika)[ncol(pika)] <- paste0(names(pika)[ncol(pika)], '_end', surveyYear)

				### chronicHeat_C: Mean of summer average temperature (tmean), June-September
				#############################################################################
					
					# generalization
					thisPred <- 'chronicHeat_C' # name of this variable
					prVar <- 'tmean' # variable to extract from PRISM/TerraClimate
					priorYears <- mostPriorYear:1 # years along window to extract
					
					# function to summarize values within a year
					summaryFx <- mean # may need to put arguments in the actual function call below!
					
					pika <- extractAndCalcVarForOcc(
						thisPred=thisPred,
						prVar=prVar,
						predTable=predTable,
						priorYears=priorYears,
						pika=pika,
						summaryFx=summaryFx,
						dateField='dummyDate',
						fail=fail,
						na.rm=TRUE
					)	

					names(pika)[ncol(pika)] <- paste0(names(pika)[ncol(pika)], '_end', surveyYear)

				### summerRespiteHeat_C: Mean of summer minimum temperature (tmin), June-September
				###################################################################################
					
					# generalization
					thisPred <- 'summerRespiteHeat_C' # name of this variable
					prVar <- 'tmin' # variable to extract from PRISM/TerraClimate
					priorYears <- mostPriorYear:1 # years along window to extract
					
					# function to summarize values within a year
					summaryFx <- mean # may need to put arguments in the actual function call below!
					
					pika <- extractAndCalcVarForOcc(
						thisPred=thisPred,
						prVar=prVar,
						predTable=predTable,
						priorYears=priorYears,
						pika=pika,
						summaryFx=summaryFx,
						dateField='dummyDate',
						fail=fail,
						na.rm=TRUE
					)	

					names(pika)[ncol(pika)] <- paste0(names(pika)[ncol(pika)], '_end', surveyYear)

				### chronicMoistStress_hPa: Mean of daily summer minimum vapor pressure deficit (vpdmin), June-September
				########################################################################################################
					
					# generalization
					thisPred <- 'chronicMoistStress_hPa' # name of this variable
					prVar <- 'vpdmin' # variable to extract from PRISM/TerraClimate
					priorYears <- mostPriorYear:1 # years along window to extract
					
					# function to summarize values within a year
					summaryFx <- mean # may need to put arguments in the actual function call below!
					
					pika <- extractAndCalcVarForOcc(
						thisPred=thisPred,
						prVar=prVar,
						predTable=predTable,
						priorYears=priorYears,
						pika=pika,
						summaryFx=summaryFx,
						dateField='dummyDate',
						fail=fail,
						na.rm=TRUE
					)	

					names(pika)[ncol(pika)] <- paste0(names(pika)[ncol(pika)], '_end', surveyYear)

				### gsVpdMin_hPa: Mean of daily summer minimum vapor pressure deficit (vpdmin), May-September
				#############################################################################################
					
					# generalization
					thisPred <- 'gsVpdMin_hPa' # name of this variable
					prVar <- 'vpdmin' # variable to extract from PRISM/TerraClimate
					priorYears <- mostPriorYear:1 # years along window to extract
					
					# function to summarize values within a year
					summaryFx <- mean # may need to put arguments in the actual function call below!
					
					pika <- extractAndCalcVarForOcc(
						thisPred=thisPred,
						prVar=prVar,
						predTable=predTable,
						priorYears=priorYears,
						pika=pika,
						summaryFx=summaryFx,
						dateField='dummyDate',
						fail=fail,
						na.rm=TRUE
					)	

					names(pika)[ncol(pika)] <- paste0(names(pika)[ncol(pika)], '_end', surveyYear)

				### peakMoistStress_hPa: Mean of summer maximum vapor pressure deficit (vpdmax), May-September
				##############################################################################################
					
					# generalization
					thisPred <- 'peakMoistStress_hPa' # name of this variable
					prVar <- 'vpdmax' # variable to extract from PRISM/TerraClimate
					priorYears <- mostPriorYear:1 # years along window to extract
					
					# function to summarize values within a year
					summaryFx <- mean # may need to put arguments in the actual function call below!
					
					pika <- extractAndCalcVarForOcc(
						thisPred=thisPred,
						prVar=prVar,
						predTable=predTable,
						priorYears=priorYears,
						pika=pika,
						summaryFx=summaryFx,
						dateField='dummyDate',
						fail=fail,
						na.rm=TRUE
					)	

					names(pika)[ncol(pika)] <- paste0(names(pika)[ncol(pika)], '_end', surveyYear)

				### gsPpt_mm: Total growing-season precipitation (ppt), May-September
				#####################################################################
					
					# generalization
					thisPred <- 'gsPpt_mm' # name of this variable
					prVar <- 'ppt' # variable to extract from PRISM/TerraClimate
					priorYears <- mostPriorYear:1 # years along window to extract
					
					# function to summarize values within a year
					summaryFx <- sum # may need to put arguments in the actual function call below!
					
					pika <- extractAndCalcVarForOcc(
						thisPred=thisPred,
						prVar=prVar,
						predTable=predTable,
						priorYears=priorYears,
						pika=pika,
						summaryFx=summaryFx,
						dateField='dummyDate',
						fail=fail,
						na.rm=TRUE
					)	

					names(pika)[ncol(pika)] <- paste0(names(pika)[ncol(pika)], '_end', surveyYear)

				### monsoonPpt_mm: Total monsoon precipitation (ppt), min-June-August
				#####################################################################
					
					# generalization
					thisPred <- 'monsoonPpt_mm' # name of this variable
					prVar <- 'ppt' # variable to extract from PRISM/TerraClimate
					priorYears <- mostPriorYear:1 # years along window to extract
					
					# function to summarize values within a year
					summaryFx <- sum # may need to put arguments in the actual function call below!
					
					pika <- extractAndCalcVarForOcc(
						thisPred=thisPred,
						prVar=prVar,
						predTable=predTable,
						priorYears=priorYears,
						pika=pika,
						summaryFx=summaryFx,
						dateField='dummyDate',
						fail=fail,
						na.rm=TRUE
					)	

					names(pika)[ncol(pika)] <- paste0(names(pika)[ncol(pika)], '_end', surveyYear)

				### winterSnow_mm: Total precipitation, December-March
				######################################################
					
					# generalization
					thisPred <- 'winterSnow_mm' # name of this variable
					prVar <- 'ppt' # variable to extract from PRISM/TerraClimate
					priorYears <- mostPriorYear:0 # years along window to extract
					
					# function to summarize values within a year
					summaryFx <- sum # may need to put arguments in the actual function call below!
					
					pika <- extractAndCalcVarForOcc(
						thisPred=thisPred,
						prVar=prVar,
						predTable=predTable,
						priorYears=priorYears,
						pika=pika,
						summaryFx=summaryFx,
						dateField='dummyDate',
						fail=fail,
						na.rm=TRUE
					)

					names(pika)[ncol(pika)] <- paste0(names(pika)[ncol(pika)], '_end', surveyYear)

				### annualPpt_mm: Total annual precipitation, June-May
				######################################################
					
					# generalization
					thisPred <- 'annualPpt_mm' # name of this variable
					prVar <- 'ppt' # variable to extract from PRISM/TerraClimate
					priorYears <- mostPriorYear:1 # years along window to extract
					
					# function to summarize values within a year
					summaryFx <- sum # may need to put arguments in the actual function call below!
					
					pika <- extractAndCalcVarForOcc(
						thisPred=thisPred,
						prVar=prVar,
						predTable=predTable,
						priorYears=priorYears,
						pika=pika,
						summaryFx=summaryFx,
						dateField='dummyDate',
						fail=fail,
						na.rm=TRUE
					)	

					names(pika)[ncol(pika)] <- paste0(names(pika)[ncol(pika)], '_end', surveyYear)

				### springPpt_mm: Total spring precipitation (ppt), mid-March-May
				#################################################################
					
					# generalization
					thisPred <- 'springPpt_mm' # name of this variable
					prVar <- 'ppt' # variable to extract from PRISM/TerraClimate
					priorYears <- mostPriorYear:0 # years along window to extract
					
					# function to summarize values within a year
					summaryFx <- sum # may need to put arguments in the actual function call below!
					
					pika <- extractAndCalcVarForOcc(
						thisPred=thisPred,
						prVar=prVar,
						predTable=predTable,
						priorYears=priorYears,
						pika=pika,
						summaryFx=summaryFx,
						dateField='dummyDate',
						fail=fail,
						na.rm=TRUE
					)
					
					names(pika)[ncol(pika)] <- paste0(names(pika)[ncol(pika)], '_end', surveyYear)

			} # next most prior year
			
		} # next survey year

	pika$dummyDate <- NULL

	write.csv(pika, './Data/XX New Mexico Pika - Environmental Values Extracted and Calculated for Bandelier Sites.csv', row.names=FALSE)
		
say('DONE!!!', level=1, deco='%')
