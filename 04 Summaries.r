### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('C:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/04 Summaries.r')
### source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/04 Summaries.r')
###
### CONTENTS ###
### setup ###
### site environmental statistics ###
### mean elevations of each occupancy class overall and by region ###
### EOO ###
### plot of variable importance by model type ###

#############
### setup ###
#############

	rm(list=ls())

	drive <- 'C:'
	# drive <- 'D:'
	# drive <- 'E:'

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

# say('###########')
# say('### EOO ###')
# say('###########')

	# load('./Data/03 New Mexico Pika - Assigned Folds.rda')
	# xy <- pika[ , c('longitude', 'latitude')]
	# xy <- as.matrix(xy)
	
	# pikaVect <- vect(xy, atts=pika, crs=getCRS('wgs84'))
	# mcp <- convHull(pikaVect)
	# area_km2 <- expanse(mcp) / 1000^2
	# say('All-sites EOO (km2): ', area_km2)

	# pikaVect <- pikaVect[pikaVect$latestOccStatus == '2 occupied', ]
	# mcp <- convHull(pikaVect)
	# area_km2 <- expanse(mcp) / 1000^2
	# say('Occupied EOO (km2): ', area_km2)

say('#################################################')
say('### plot of variable importance by model type ###')
say('#################################################')

	# occWindow_y <- 7
	occWindow_y <- 10

	viOrdinal <- read.csv(paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Ordinal Models Using All Data - ', occWindow_y, '-yr Window - Var Import.csv'))
	viBinary <- read.csv(paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Binary Models Using All Data - ', occWindow_y, '-yr Window - Var Import.csv'))
	viDensity <- read.csv(paste0('./Figures & Tables/Density - Simple Models/Density - Simple GLMs - Var Import.csv'))
	
	viOrdinal <- viOrdinal[order(viOrdinal$meanWeight, decreasing=TRUE), ]
	viBinary <- viBinary[order(viBinary$meanWeight, decreasing=TRUE), ]
	viDensity <- viDensity[order(viDensity$meanWeight, decreasing=TRUE), ]

	### colors and fills
	####################
	
	cols <- data.frame(
		pred = c(
			'chronic heat',
			'GS precipitation',
			'Isolation',
			'sub-lethal heat',
			'sub-lethal cold',
			'monsoon precipitation',
			'chronic moisture deficit',
			'chronic cold',
			'peak moisture deficit',
			'acute cold',
			'GS moisture deficit',
			'acute heat',
			'winter snowfall',
			'summer respite',
			'spring precipitation',
			'annual precipitation'
		),
		fill = c(
			'#FF0000',
			'#00FF00',
			'#D3D3D3',
			'#FFFF00',
			'#008B8B',
			'#00008B',
			'#006400',
			'#0000FF',
			'gray30',
			'gray30',
			'gray30',
			'gray30',
			'#00FFFF',
			'gray30',
			'gray30',
			'gray30'
		)
	)
	

		fills = c(
			'chronic heat' = '#FF0000',
			'GS precipitation' = '#00FF00',
			'Isolation' = '#D3D3D3',
			'sub-lethal heat' = '#FFFF00',
			'sub-lethal cold' = '#008B8B',
			'monsoon precipitation' = '#00008B',
			'chronic moisture deficit' = '#006400',
			'chronic cold' = '#0000FF',
			'peak moisture deficit' = 'gray30',
			'acute cold' = 'gray30',
			'GS moisture deficit' = 'gray30',
			'acute heat' = 'gray30',
			'winter snowfall' = '#00FFFF',
			'summer respite' = 'gray30',
			'spring precipitation' = 'gray30',
			'annual precipitation' = 'gray30'
		)


	### ordinal occupancy
	#####################
	
	x <- viOrdinal
	
	x$sumWeightNice <- round(x$sumWeight, 2)
	x$sumWeightNice[x$sumWeightNice > 0] <- sprintf('%.2f', x$sumWeightNice[x$sumWeightNice > 0])
	
	predsNice <- gsub(x$niceVar, pattern=paste0(' \\(', occWindow_y, '-yr\\)'), replacement='')
	predsNice <- capIt(predsNice)
	x$predsNice <- predsNice
	
	x$predsNice <- factor(x$predsNice, levels = x$predsNice[order(x$meanWeight)])
	
	xSumWeight <- x
	xSumWeight <- xSumWeight[xSumWeight$sumWeight > 0.001, ]
	
	ylim <- c(0, 1.15 * max(x$meanWeight))
	
	x$fill <- NA
	x$fill <- cols$fill[match(tolower(x$predsNice), tolower(cols$pred))]
	
	ordinal <- ggplot(x, aes(x=predsNice, y=meanWeight, fill=predsNice)) +
		geom_bar(stat='identity') +
		scale_fill_manual(
			values = c(
				'Chronic Heat' = '#FF0000',
				'GS Precipitation' = '#00FF00',
				'Isolation' = '#D3D3D3',
				'Sub-lethal Heat' = '#FFFF00',
				'Sub-lethal Cold' = '#008B8B',
				'Monsoon Precipitation' = '#00008B',
				'Chronic Moisture Deficit' = '#006400',
				'Chronic Cold' = '#0000FF',
				'Peak Moisture Deficit' = 'gray30',
				'Acute Cold' = 'gray30',
				'GS Moisture Deficit' = 'gray30',
				'Acute Heat' = 'gray30',
				'Winter Snowfall' = '#00FFFF',
				'Summer Respite' = 'gray30',
				'Spring Precipitation' = 'gray30',
				'Annual Precipitation' = 'gray30'
			)
		) +
		geom_text(data=xSumWeight, aes(y=meanWeight, label=sumWeightNice), hjust=-0.2) +
		coord_flip() +
		ylim(ylim[1], ylim[2]) +
		xlab('') + ylab('Mean AICc Weight') +
		ggtitle('(a) Ordinal occupancy') +
		theme(
			legend.position = 'none'
		)
	
	### binary occupancy
	####################
	
	x <- viBinary
	
	x$sumWeightNice <- round(x$sumWeight, 2)
	x$sumWeightNice[x$sumWeightNice > 0] <- sprintf('%.2f', x$sumWeightNice[x$sumWeightNice > 0])
	
	predsNice <- gsub(x$niceVar, pattern=paste0(' \\(', occWindow_y, '-yr\\)'), replacement='')
	predsNice <- capIt(predsNice)
	x$predsNice <- predsNice
	
	x$predsNice <- factor(x$predsNice, levels = x$predsNice[order(x$meanWeight)])
	
	xSumWeight <- x
	xSumWeight <- xSumWeight[xSumWeight$sumWeight > 0.001, ]
	
	ylim <- c(0, 1.15 * max(x$meanWeight))
	
	binary <- ggplot(x, aes(x=predsNice, y=meanWeight, fill=predsNice)) +
		geom_bar(stat='identity') +
		scale_fill_manual(
			values = c(
				'Chronic Heat' = '#FF0000',
				'GS Precipitation' = '#00FF00',
				'Isolation' = '#D3D3D3',
				'Sub-lethal Heat' = '#FFFF00',
				'Sub-lethal Cold' = '#008B8B',
				'Monsoon Precipitation' = '#00008B',
				'Chronic Moisture Deficit' = '#006400',
				'Chronic Cold' = '#0000FF',
				'Peak Moisture Deficit' = 'gray30',
				'Acute Cold' = 'gray30',
				'GS Moisture Deficit' = 'gray30',
				'Acute Heat' = 'gray30',
				'Winter Snowfall' = '#00FFFF',
				'Summer Respite' = 'gray30',
				'Spring Precipitation' = 'gray30',
				'Annual Precipitation' = 'gray30'
			)
		) +
		geom_text(data=xSumWeight, aes(y=meanWeight, label=sumWeightNice), hjust=-0.2) +
		coord_flip() +
		ylim(ylim[1], ylim[2]) +
		xlab('') + ylab('Mean AICc Weight') +
		ggtitle('(b) Binary occupancy') +
		theme(
			legend.position = 'none'
		)
	
	occsVert <- plot_grid(ordinal, binary, align='h', ncol=1, rel_widths=1, labels=NULL)
	occsVarImpHoriz <- plot_grid(ordinal, binary, align='h', ncol=2, rel_widths=1, labels=NULL)

	### density
	###########
	
	x <- viDensity
	
	x$sumWeightNice <- round(x$sumWeight, 2)
	x$sumWeightNice[x$sumWeightNice > 0] <- sprintf('%.2f', x$sumWeightNice[x$sumWeightNice > 0])
	
	predsNice <- gsub(x$niceVar, pattern=paste0(' \\(', occWindow_y, '-yr\\)'), replacement='')
	predsNice <- capIt(predsNice)
	predsNice <- gsub(predsNice, pattern='Yr', replacement='yr')
	x$predsNice <- predsNice
	
	x$predsNice <- factor(x$predsNice, levels = x$predsNice[order(x$meanWeight)])
	
	xSumWeight <- x
	xSumWeight <- xSumWeight[xSumWeight$sumWeight > 0.001, ]
	
	ylim <- c(0, 1.225 * max(x$meanWeight))
	
	breaks <- pretty(x$meanWeight, 3)
	
	density <- ggplot(x, aes(x=predsNice, y=meanWeight, fill=predsNice)) +
		geom_bar(stat='identity') +
		scale_fill_manual(
			values = c(
				'Chronic Heat (1 yr)' = '#FF0000',
				'Chronic Heat (0 yr)' = '#FF0000',
				'GS Precipitation (1 yr)' = '#00FF00',
				'GS Precipitation (0 yr)' = '#00FF00',
				'Isolation' = '#D3D3D3',
				'Sub-lethal Heat (1 yr)' = '#FFFF00',
				'Sub-lethal Heat (0 yr)' = '#FFFF00',
				'Sub-lethal Cold (1 yr)' = '#008B8B',
				'Sub-lethal Cold (0 yr)' = '#008B8B',
				'Monsoon Precipitation (1 yr)' = '#00008B',
				'Monsoon Precipitation (0 yr)' = '#00008B',
				'Chronic Moisture Deficit (1 yr)' = '#006400',
				'Chronic Moisture Deficit (0 yr)' = '#006400',
				'Chronic Cold (1 yr)' = '#0000FF',
				'Chronic Cold (0 yr)' = '#0000FF',
				'Peak Moisture Deficit (1 yr)' = 'gray30',
				'Peak Moisture Deficit (0 yr)' = 'gray30',
				'Acute Cold (1 yr)' = 'gray30',
				'Acute Cold (0 yr)' = 'gray30',
				'GS Moisture Deficit (1 yr)' = 'gray30',
				'GS Moisture Deficit (0 yr)' = 'gray30',
				'Acute Heat (1 yr)' = '#8B0000',
				'Acute Heat (0 yr)' = '#8B0000',
				'Winter Snowfall (1 yr)' = '#00FFFF',
				'Winter Snowfall (0 yr)' = '#00FFFF',
				'Summer Respite (1 yr)' = '#800080',
				'Summer Respite (0 yr)' = '#800080',
				'Spring Precipitation (1 yr)' = 'gray30',
				'Spring Precipitation (0 yr)' = 'gray30',
				'Annual Precipitation (1 yr)' = 'gray30',
				'Annual Precipitation (0 yr)' = 'gray30'
			)
		) +
		geom_text(data=xSumWeight, aes(y=meanWeight, label=sumWeightNice), hjust=-0.2) +
		coord_flip() +
		xlab('') + ylab('Mean AICc Weight') +
		scale_y_continuous(limits = ylim, breaks = breaks) +
		ggtitle('(c) Density') +
		theme(
			legend.position = 'none'
		)
		

	occsVarImpHoriz <- occsVarImpHoriz +
		theme(
			plot.margin = unit(c(0, 0.4, 0, 0), 'cm')
		)

	occDensVarImp <- plot_grid(occsVert, density, ncol=2, rel_widths=1, labels=NULL)
	occDensVarImp <- occDensVarImp +
		theme(
			plot.margin = unit(c(0, 0.4, 0, 0), 'cm')
		)


	ggsave(occsVarImpHoriz, file=paste0('./Figures & Tables/Variable Importance for All Models (Using ', occWindow_y, '-yr Occupancy Window - Occupancy.pdf'), width=8, height=6)
	ggsave(occDensVarImp, file=paste0('./Figures & Tables/Variable Importance for All Models (Using ', occWindow_y, '-yr Occupancy Window.pdf'), width=8, height=8)
	
say('DONE!!!', level=1, deco='%')
