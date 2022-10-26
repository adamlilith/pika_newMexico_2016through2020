### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('C:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/04 Site Statistics.r')
### source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/04 Site Statistics.r')
###
### CONTENTS ###
### setup ###
### site environmental statistics ###
### mean elevations of each occupancy class overall and by region ###
### EOO ###

#############
### setup ###
#############

	# drive <- 'C:'
	# drive <- 'D:'
	drive <- 'E:'

	source(paste0(drive, '/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r'))

# say('#####################################')		
# say('### site environmental statistics ###')
# say('#####################################')

	# load('./Data/03 New Mexico Pika - Assigned Folds.rda')

	# say('MEAN ANNUAL PRECIPITATION (mm)')
	# x <- pika$occVar_annualPpt_mm_10yrWindow
	# say('Min: ', min(x))
	# say('Max: ', max(x))
	# say('Mean: ', mean(x))
	# say('Range: ', paste(range(x)))

	# say('SUMMER TEMPERATURE (AVERAGE June-Sept, deg C)')
	# x <- pika$occVar_chronicHeat_C_10yrWindow
	# say('Min (mm): ', min(x))
	# say('Max (mm): ', max(x))
	# say('Mean (mm): ', mean(x))
	# say('Range: ', paste(range(x)))

	# say('WINTER TEMPERATURE (AVERAGE Nov-March, deg C)')
	# x <- pika$occVar_chronicCold_C_10yrWindow
	# say('Min (mm): ', min(x))
	# say('Max (mm): ', max(x))
	# say('Mean (mm): ', mean(x))
	# say('Range: ', paste(range(x)))

# say('#####################################################################')
# say('### mean elevations of each occupancy class overall and by region ###')
# say('#####################################################################')

	# load('./Data/03 New Mexico Pika - Assigned Folds.rda')


	# never <- round(mean(pika$elevation_m[pika$latestOccStatus == '0 never']), 0)
	# old <- round(mean(pika$elevation_m[pika$latestOccStatus == '1 old']), 0)
	# occ <- round(mean(pika$elevation_m[pika$latestOccStatus == '2 occupied']), 0)

	# results <- data.frame(
		# region = 'all',
		# occ = occ,
		# old = old,
		# never = never
	# )

	# for (region in c('northwest', 'southwest', 'northeast', 'southeast')) {
	
		# occ <- round(mean(pika$elevation_m[pika$latestOccStatus == '2 occupied' & pika$region == region]), 0)
		# old <- round(mean(pika$elevation_m[pika$latestOccStatus == '1 old' & pika$region == region]), 0)
		# never <- round(mean(pika$elevation_m[pika$latestOccStatus == '0 never' & pika$region == region]), 0)

		# results <- rbind(
			# results,
				# data.frame(
				# region = region,
				# occ = occ,
				# old = old,
				# never = never
			# )
		# )
	
	# }

	# write.csv(results, './Figures & Tables/Mean Elevations by Site Status across All Sites and by Region.csv', row.names=FALSE)

say('###########')
say('### EOO ###')
say('###########')

	load('./Data/03 New Mexico Pika - Assigned Folds.rda')
	xy <- pika[ , c('longitude', 'latitude')]
	xy <- as.matrix(xy)
	
	pikaVect <- vect(xy, atts=pika, crs=getCRS('wgs84'))
	mcp <- convHull(pikaVect)
	area_km2 <- expanse(mcp) / 1000^2
	say('All-sites EOO (km2): ', area_km2)

	pikaVect <- pikaVect[pikaVect$latestOccStatus == '2 occupied', ]
	mcp <- convHull(pikaVect)
	area_km2 <- expanse(mcp) / 1000^2
	say('Occupied EOO (km2): ', area_km2)
	
say('DONE!!!', level=1, deco='%')
