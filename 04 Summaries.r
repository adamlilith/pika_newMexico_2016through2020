### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('C:/Ecology/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/04 Summaries.r')
### source('E:/Adam/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/04 Summaries.r')
###
### CONTENTS ###
### setup ###
### site environmental statistics ###
### summarize occupancy status across regions ###
### mean elevations of each occupancy class overall and by region ###
### EOO ###
### plot of variable importance by model type ###

### distributions of predictors by region ###
### distributions of predictors by region and ordinal occupancy for main text ###
### distributions of predictors by region and binary occupancy for main text ###
### distributions of predictors by region and binary and ordinal occupancy for main text ###
### contingency table analysis of site status by region ###

### statistical comparison of environmental variables across regions ###


#############
### setup ###
#############

	rm(list=ls())

	drive <- 'C:/Ecology/'
	# drive <- 'E:/Adam/'

	source(paste0(drive, '/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r'))

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

# say('#################################################')
# say('### summarize occupancy status across regions ###')
# say('#################################################')

	# load('./Data/05 New Mexico Pika - Added PRISM Cell Number.rda')

	# sink('./Figures & Tables/Occupancy - Simple Models/Occupancy - Status by Region.txt')

		# say('Summary of occupancy status by region')
		# say(date(), post=2)

		# say('ALL REGIONS:')
		# n <- nrow(pika)
		# unis <- length(unique(pika$prismCellNum))
		# say('     Never occupied:                 ', sum(pika$latestOccStatus == '0 never'), ' (', round(sum(pika$latestOccStatus == '0 never') / n, 2), ')')
		# say('     Old evidence:                   ', sum(pika$latestOccStatus == '1 old'), ' (', round(sum(pika$latestOccStatus == '1 old') / n, 2), ')')
		# say('     Occupied:                       ', sum(pika$latestOccStatus == '2 occupied'), ' (', round(sum(pika$latestOccStatus == '2 occupied') / n, 2), ')')
		# say('     All sites:                      ', n)
		# y <- pika[!duplicated(pika$prismCellNum), ]
		# say('     Sites not in same PRISM cell:   ', n - nrow(y), pre = 1, post = 2)

		# say('NORTHEAST:')
		# region <- 'northeast'
		# n <- sum(pika$region == region)
		# say('     Never occupied:                 ', sum(pika$region == region & pika$latestOccStatus == '0 never'), ' (', round(sum(pika$region == region & pika$latestOccStatus == '0 never') / n, 2), ')')
		# say('     Old evidence:                   ', sum(pika$region == region & pika$latestOccStatus == '1 old'), ' (', round(sum(pika$region == region & pika$latestOccStatus == '1 old') / n, 2), ')')
		# say('     Occupied:                       ', sum(pika$region == region & pika$latestOccStatus == '2 occupied'), ' (', round(sum(pika$region == region & pika$latestOccStatus == '2 occupied') / n, 2), ')')
		# say('     All sites:                      ', n)
		# x <- pika[pika$region == region, ]
		# y <- x[!duplicated(x$prismCellNum), ]
		# say('     Sites not in same PRISM cell:   ', nrow(x) - nrow(y), pre = 1)

		# say('NORTHWEST:', pre=1)
		# region <- 'northwest'
		# n <- sum(pika$region == region)
		# say('     Never occupied:                 ', sum(pika$region == region & pika$latestOccStatus == '0 never'), ' (', round(sum(pika$region == region & pika$latestOccStatus == '0 never') / n, 2), ')')
		# say('     Old evidence:                   ', sum(pika$region == region & pika$latestOccStatus == '1 old'), ' (', round(sum(pika$region == region & pika$latestOccStatus == '1 old') / n, 2), ')')
		# say('     Occupied:                       ', sum(pika$region == region & pika$latestOccStatus == '2 occupied'), ' (', round(sum(pika$region == region & pika$latestOccStatus == '2 occupied') / n, 2), ')')
		# say('     All sites:                      ', n)
		# x <- pika[pika$region == region, ]
		# y <- x[!duplicated(x$prismCellNum), ]
		# say('     Sites not in same PRISM cell:   ', nrow(x) - nrow(y), pre = 1)

		# say('SOUTHEAST:', pre=1)
		# region <- 'southeast'
		# n <- sum(pika$region == region)
		# say('     Never occupied:                 ', sum(pika$region == region & pika$latestOccStatus == '0 never'), ' (', round(sum(pika$region == region & pika$latestOccStatus == '0 never') / n, 2), ')')
		# say('     Old evidence:                   ', sum(pika$region == region & pika$latestOccStatus == '1 old'), ' (', round(sum(pika$region == region & pika$latestOccStatus == '1 old') / n, 2), ')')
		# say('     Occupied:                 ', sum(pika$region == region & pika$latestOccStatus == '2 occupied'), ' (', round(sum(pika$region == region & pika$latestOccStatus == '2 occupied') / n, 2), ')')
		# say('     All sites:                      ', n)
		# x <- pika[pika$region == region, ]
		# y <- x[!duplicated(x$prismCellNum), ]
		# say('     Sites not in same PRISM cell:   ', nrow(x) - nrow(y), pre = 1)

		# say('SOUTHWEST:', pre=1)
		# region <- 'southwest'
		# n <- sum(pika$region == region)
		# say('     Never occupied:                 ', sum(pika$region == region & pika$latestOccStatus == '0 never'), ' (', round(sum(pika$region == region & pika$latestOccStatus == '0 never') / n, 2), ')')
		# say('     Old evidence:                   ', sum(pika$region == region & pika$latestOccStatus == '1 old'), ' (', round(sum(pika$region == region & pika$latestOccStatus == '1 old') / n, 2), ')')
		# say('     Occupied:                       ', sum(pika$region == region & pika$latestOccStatus == '2 occupied'), ' (', round(sum(pika$region == region & pika$latestOccStatus == '2 occupied') / n, 2), ')')
		# say('     All sites:                      ', n)
		# x <- pika[pika$region == region, ]
		# y <- x[!duplicated(x$prismCellNum), ]
		# say('     Sites not in same PRISM cell:   ', nrow(x) - nrow(y), pre = 1)

	# sink()

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

# say('#################################################')
# say('### plot of variable importance by model type ###')
# say('#################################################')

	# # occWindow_y <- 7
	# occWindow_y <- 10

	# viOrdinal <- read.csv(paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Ordinal Models Using All Data - ', occWindow_y, '-yr Window - Var Import.csv'))
	# viBinary <- read.csv(paste0('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Binary Models Using All Data - ', occWindow_y, '-yr Window - Var Import.csv'))
	# viDensity <- read.csv(paste0('./Figures & Tables/Density - Simple Models/Density - Simple GLMs - Var Import.csv'))
	
	# viOrdinal <- viOrdinal[order(viOrdinal$meanWeight, decreasing=TRUE), ]
	# viBinary <- viBinary[order(viBinary$meanWeight, decreasing=TRUE), ]
	# viDensity <- viDensity[order(viDensity$meanWeight, decreasing=TRUE), ]

	# ### colors and fills
	# ####################
	
	# cols <- data.frame(
		# pred = c(
			# 'chronic heat',
			# 'GS precipitation',
			# 'Isolation',
			# 'sub-lethal heat',
			# 'sub-lethal cold',
			# 'monsoon precipitation',
			# 'chronic moisture deficit',
			# 'chronic cold',
			# 'peak moisture deficit',
			# 'acute cold',
			# 'GS moisture deficit',
			# 'acute heat',
			# 'winter snowfall',
			# 'summer respite',
			# 'spring precipitation',
			# 'annual precipitation'
		# ),
		# fill = c(
			# '#FF070B', # chronic heat
			# '#A8D08D', # GS precip
			# '#E7E6E6', # isolation
			# '#FFC000',
			# '#DEEBF7',
			# '#538135',
			# '#BF9000',
			# '#9FC5D1',
			# '#FFD966',
			# '#4472C4',
			# '#C55A11',
			# '#C00000',
			# '#FFFFFF',
			# '#F8CBAD',
			# '#C5E0B4',
			# '#92D050'
		# )
	# )
	
	# fills <- cols$fill
	# names(fills) <- capIt(cols$pred)
	
	# ### ordinal occupancy
	# #####################
	
	# x <- viOrdinal
	
	# x$sumWeightNice <- round(x$sumWeight, 2)
	# x$sumWeightNice[x$sumWeightNice > 0] <- sprintf('%.2f', x$sumWeightNice[x$sumWeightNice > 0])
	
	# predsNice <- gsub(x$niceVar, pattern=paste0(' \\(', occWindow_y, '-yr\\)'), replacement='')
	# predsNice <- capIt(predsNice)
	# x$predsNice <- predsNice
	
	# x$predsNice <- factor(x$predsNice, levels = x$predsNice[order(x$meanWeight)])
	
	# xSumWeight <- x
	# xSumWeight <- xSumWeight[xSumWeight$sumWeight >= 0.01, ]
	
	# ylim <- c(0, 1.15 * max(x$meanWeight))
	
	# x$fill <- NA
	# x$fill <- cols$fill[match(tolower(x$predsNice), tolower(cols$pred))]
	
	# ordinal <- ggplot(x, aes(x=predsNice, y=meanWeight, fill=predsNice)) +
		# geom_bar(stat='identity', color = 'black', linewidth = 0.2) +
		# scale_fill_manual(
			# values = fills
		# ) +
		# geom_text(data = xSumWeight, aes( y = meanWeight, label=sumWeightNice), hjust=-0.2) +
		# coord_flip() +
		# ylim(ylim[1], ylim[2]) +
		# xlab('') + ylab('Mean AICc Weight') +
		# ggtitle('(a) Ordinal occupancy') +
		# theme(
			# legend.position = 'none'
		# )
	
	# ### binary occupancy
	# ####################
	
	# x <- viBinary
	
	# x$sumWeightNice <- round(x$sumWeight, 2)
	# x$sumWeightNice[x$sumWeightNice > 0] <- sprintf('%.2f', x$sumWeightNice[x$sumWeightNice > 0])
	
	# predsNice <- gsub(x$niceVar, pattern=paste0(' \\(', occWindow_y, '-yr\\)'), replacement='')
	# predsNice <- capIt(predsNice)
	# x$predsNice <- predsNice
	
	# x$predsNice <- factor(x$predsNice, levels = x$predsNice[order(x$meanWeight)])
	
	# xSumWeight <- x
	# xSumWeight <- xSumWeight[xSumWeight$sumWeight >= 0.01, ]
	
	# ylim <- c(0, 1.15 * max(x$meanWeight))
	
	# binary <- ggplot(x, aes(x=predsNice, y=meanWeight, fill=predsNice)) +
		# geom_bar(stat='identity', color = 'black', linewidth = 0.2) +
		# scale_fill_manual(
			# values = fills
		# ) +
		# geom_text(data=xSumWeight, aes(y=meanWeight, label=sumWeightNice), hjust=-0.2) +
		# coord_flip() +
		# ylim(ylim[1], ylim[2]) +
		# xlab('') + ylab('Mean AICc Weight') +
		# ggtitle('(b) Binary occupancy') +
		# theme(
			# legend.position = 'none'
		# )
	
	# occsVert <- plot_grid(ordinal, binary, align='h', ncol=1, rel_widths=1, labels=NULL)
	# occsVarImpHoriz <- plot_grid(ordinal, binary, align='h', ncol=2, rel_widths=1, labels=NULL)

	# ### density
	# ###########
	
	# x <- viDensity
	
	# x$sumWeightNice <- round(x$sumWeight, 2)
	# x$sumWeightNice[x$sumWeightNice > 0] <- sprintf('%.2f', x$sumWeightNice[x$sumWeightNice > 0])
	
	# predsNice <- gsub(x$niceVar, pattern=paste0(' \\(', occWindow_y, '-yr\\)'), replacement='')
	# predsNice <- capIt(predsNice)
	# predsNice <- gsub(predsNice, pattern='Yr', replacement='yr')
	# x$predsNice <- predsNice
	
	# x$predsNice <- factor(x$predsNice, levels = x$predsNice[order(x$meanWeight)])
	
	# xSumWeight <- x
	# xSumWeight <- xSumWeight[xSumWeight$sumWeight >= 0.01, ]
	
	# ylim <- c(0, 1.225 * max(x$meanWeight))
	
	# breaks <- pretty(x$meanWeight, 3)
	
	# density <- ggplot(x, aes(x=predsNice, y=meanWeight, fill=predsNice)) +
		# geom_bar(stat='identity', color = 'black', linewidth = 0.2) +
		# scale_fill_manual(
			# values = c(
				# 'Chronic Heat (1 yr)' = unname(fills['Chronic Heat']),
				# 'Chronic Heat (0 yr)' = unname(fills['Chronic Heat']),
				# 'GS Precipitation (1 yr)' = unname(fills['GS Precipitation']),
				# 'GS Precipitation (0 yr)' = unname(fills['GS Precipitation']),
				# 'Isolation' = unname(fills['Isolation']),
				# 'Sub-lethal Heat (1 yr)' = unname(fills['Sub-lethal Heat']),
				# 'Sub-lethal Heat (0 yr)' = unname(fills['Sub-lethal Heat']),
				# 'Sub-lethal Cold (1 yr)' = unname(fills['Sub-lethal Cold']),
				# 'Sub-lethal Cold (0 yr)' = unname(fills['Sub-lethal Cold']),
				# 'Monsoon Precipitation (1 yr)' = unname(fills['Monsoon Precipitation']),
				# 'Monsoon Precipitation (0 yr)' = unname(fills['Monsoon Precipitation']),
				# 'Chronic Moisture Deficit (1 yr)' = unname(fills['Chronic Moisture Deficit']),
				# 'Chronic Moisture Deficit (0 yr)' = unname(fills['Chronic Moisture Deficit']),
				# 'Chronic Cold (1 yr)' = unname(fills['Chronic Cold']),
				# 'Chronic Cold (0 yr)' = unname(fills['Chronic Cold']),
				# 'Peak Moisture Deficit (1 yr)' = unname(fills['Peak Moisture Deficit']),
				# 'Peak Moisture Deficit (0 yr)' = unname(fills['Peak Moisture Deficit']),
				# 'Acute Cold (1 yr)' = unname(fills['Acute Cold']),
				# 'Acute Cold (0 yr)' = unname(fills['Acute Cold']),
				# 'GS Moisture Deficit (1 yr)' = unname(fills['GS Moisture Deficit']),
				# 'GS Moisture Deficit (0 yr)' = unname(fills['GS Moisture Deficit']),
				# 'Acute Heat (1 yr)' = unname(fills['Acute Heat']),
				# 'Acute Heat (0 yr)' = unname(fills['Acute Heat']),
				# 'Winter Snowfall (1 yr)' = unname(fills['Winter Snowfall']),
				# 'Winter Snowfall (0 yr)' = unname(fills['Winter Snowfall']),
				# 'Summer Respite (1 yr)' = unname(fills['Summer Respite']),
				# 'Summer Respite (0 yr)' = unname(fills['Summer Respite']),
				# 'Spring Precipitation (1 yr)' = unname(fills['Spring Precipitation']),
				# 'Spring Precipitation (0 yr)' = unname(fills['Spring Precipitation']),
				# 'Annual Precipitation (1 yr)' = unname(fills['Annual Precipitation']),
				# 'Annual Precipitation (0 yr)' = unname(fills['Annual Precipitation'])
			# )
		# ) +
		# geom_text(data=xSumWeight, aes(y=meanWeight, label=sumWeightNice), hjust=-0.2) +
		# coord_flip() +
		# xlab('') + ylab('Mean AICc Weight') +
		# scale_y_continuous(limits = ylim, breaks = breaks) +
		# ggtitle('(c) Density') +
		# theme(
			# legend.position = 'none'
		# )
		

	# occsVarImpHoriz <- occsVarImpHoriz +
		# theme(
			# plot.margin = unit(c(0, 0.4, 0, 0), 'cm')
		# )

	# occDensVarImp <- plot_grid(occsVert, density, ncol=2, rel_widths=1, labels=NULL)
	# occDensVarImp <- occDensVarImp +
		# theme(
			# plot.margin = unit(c(0, 0.4, 0, 0), 'cm')
		# )


	# # ggsave(occsVarImpHoriz, file=paste0('./Figures & Tables/Variable Importance for All Models (Using ', occWindow_y, '-yr Occupancy Window) - Occupancy.pdf'), width=8, height=6)
	# ggsave(occDensVarImp, file=paste0('./Figures & Tables/Variable Importance for All Models (Using ', occWindow_y, '-yr Occupancy Window).pdf'), width=8, height=8)
	
# say('#############################################')
# say('### distributions of predictors by region ###')
# say('#############################################')

	# titleSize <- 8
	# subtitleSize <- 7
	# legendTitleSize <- 7
	# legendTextSize <- 7
	# axisLabelSize <- 7
	# axisTextSize <- 6
	
	# legendKeySize <- 0.4

	# ### data
	# load('./Data/03 New Mexico Pika - Assigned Folds.rda')

	# ### occupancy predictors
	# ########################
	
	# pika$latestOccStatus <- factor(pika$latestOccStatus, levels=c('0 never', '1 old', '2 occupied'))
	
	# preds <- predTable$var[predTable$useOccAbund]

	# for (occWindow in occWindows_y) {
		
		# figs <- list()
		# for (countPred in seq_along(preds)) {
		
			# pred <- preds[countPred]
			
			# predWindow <- paste0('occVar_', pred, '_', occWindow, 'yrWindow')
			# say(predWindow)
		
			# predIndex <- which(predTable$var == pred)
		
			# predNice <- predTable$varNice[predIndex]
			# predDescriptorUnit <- paste0(predTable$unitDescriptor[predIndex], ' (', predTable$unit[predIndex], ')')
		
			# # mus <- data.frame(
				# # latestOccStatus = c('0 never', '1 old', '2 occupied'),
				# # mu = c(
					# # mean(pika[pika$latestOccStatus == '0 never', predWindow]),
					# # mean(pika[pika$latestOccStatus == '1 old', predWindow]),
					# # mean(pika[pika$latestOccStatus == '2 occupied', predWindow])
				# # )
			# # )

			# thisData <- pika[ , c('latestOccStatus', 'region', predWindow)]
			# names(thisData)[3] <- 'value'

			# xlim <- range(thisData$value)
			
			# # all regions together
			# title <- paste0(predNice, ': regions')
			# figs[[length(figs) + 1]] <- ggplot(data=thisData, aes(x=value, col=latestOccStatus, fill=latestOccStatus)) +
				# geom_density(size=1) +
				# scale_color_manual(
					# labels = c('none', 'old', 'occ'),
					# values=c('0 never'='firebrick3', '1 old'='darkgoldenrod3', '2 occupied'='darkgreen')
				# ) +
				# scale_fill_manual(
					# labels = c('none', 'old', 'occ'),
					# values=alpha(c('0 never'='firebrick3', '1 old'='darkgoldenrod3', '2 occupied'='darkgreen'), 0.2)
				# ) +
				# labs(title=title, subtitle=pred, x=predDescriptorUnit, y=NULL) +
				# # geom_vline(data=mus, aes(xintercept=mu, color=latestOccStatus), linetype='dotted', size=1) +
				# guides(
					# color=guide_legend(title='Evid.'),
					# fill=guide_legend(title='Evid.')
				# ) +
				# xlim(xlim[1], xlim[2]) +
				# theme(
					# legend.key.size = unit(legendKeySize, 'cm'),
					# plot.title=element_text(size=titleSize, face='bold'),
					# plot.subtitle=element_text(size=subtitleSize),
					# legend.title=element_text(size=legendTitleSize),
					# legend.text=element_text(size=legendTextSize),
					# axis.title=element_text(size=axisLabelSize),
					# axis.text=element_text(size=axisTextSize)
				# )

			# # "never" by region
			# title <- paste0(predNice, ': none')
			# thisThisData <- thisData[thisData$latestOccStatus == '0 never', ]
			# figs[[length(figs) + 1]] <- ggplot(data=thisThisData, aes(x=value, col=region, fill=region)) +
				# geom_density(size=1) +
				# scale_color_manual(
					# labels=c('southwest'='SW', 'southeast'='SE', 'northwest'='NW', 'northeast'='NE'),
					# values=c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred')
				# ) +
				# scale_fill_manual(
					# labels=c('southwest'='SW', 'southeast'='SE', 'northwest'='NW', 'northeast'='NE'),
					# values=alpha(c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred'), 0.2),
				# ) +
				# labs(title=title, subtitle=pred, x=predDescriptorUnit, y=NULL) +
				# guides(
					# color=guide_legend(title='Region'),
					# fill=guide_legend(title='Region')
				# ) +
				# xlim(xlim[1], xlim[2]) +
				# theme(
					# legend.key.size = unit(legendKeySize, 'cm'),
					# plot.title=element_text(size=titleSize, face='bold'),
					# plot.subtitle=element_text(size=subtitleSize),
					# legend.title=element_text(size=legendTitleSize),
					# legend.text=element_text(size=legendTextSize),
					# axis.title=element_text(size=axisLabelSize),
					# axis.text=element_text(size=axisTextSize)
				# )

			# # "old" by region
			# title <- paste0(predNice, ': old')
			# thisThisData <- thisData[thisData$latestOccStatus == '1 old', ]
			# figs[[length(figs) + 1]] <- ggplot(data=thisThisData, aes(x=value, col=region, fill=region)) +
				# geom_density(size=1) +
				# scale_color_manual(
					# labels=c('southwest'='SW', 'southeast'='SE', 'northwest'='NW', 'northeast'='NE'),
					# values=c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred')
				# ) +
				# scale_fill_manual(
					# labels=c('southwest'='SW', 'southeast'='SE', 'northwest'='NW', 'northeast'='NE'),
					# values=alpha(c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred'), 0.2),
				# ) +
				# labs(title=title, subtitle=pred, x=predDescriptorUnit, y=NULL) +
				# guides(
					# color=guide_legend(title='Region'),
					# fill=guide_legend(title='Region')
				# ) +
				# xlim(xlim[1], xlim[2]) +
				# theme(
					# legend.key.size = unit(legendKeySize, 'cm'),
					# plot.title=element_text(size=titleSize, face='bold'),
					# plot.subtitle=element_text(size=subtitleSize),
					# legend.title=element_text(size=legendTitleSize),
					# legend.text=element_text(size=legendTextSize),
					# axis.title=element_text(size=axisLabelSize),
					# axis.text=element_text(size=axisTextSize)
				# )

			# # "occupied" by region
			# title <- paste0(predNice, ': occ')
			# thisThisData <- thisData[thisData$latestOccStatus == '2 occupied', ]
			# figs[[length(figs) + 1]] <- ggplot(data=thisThisData, aes(x=value, col=region, fill=region)) +
				# geom_density(size=1) +
				# scale_color_manual(
					# labels=c('southwest'='SW', 'southeast'='SE', 'northwest'='NW', 'northeast'='NE'),
					# values=c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred')
				# ) +
				# scale_fill_manual(
					# labels=c('southwest'='SW', 'southeast'='SE', 'northwest'='NW', 'northeast'='NE'),
					# values=alpha(c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred'), 0.2),
				# ) +
				# labs(title=title, subtitle=pred, x=predDescriptorUnit, y=NULL) +
				# guides(
					# color=guide_legend(title='Region'),
					# fill=guide_legend(title='Region')
				# ) +
				# xlim(xlim[1], xlim[2]) +
				# theme(
					# legend.key.size = unit(legendKeySize, 'cm'),
					# plot.title=element_text(size=titleSize, face='bold'),
					# plot.subtitle=element_text(size=subtitleSize),
					# legend.title=element_text(size=legendTitleSize),
					# legend.text=element_text(size=legendTextSize),
					# axis.title=element_text(size=axisLabelSize),
					# axis.text=element_text(size=axisTextSize)
				# )

				
		# } # next predictor

		# main <- plot_grid(plotlist=figs[1:(7 * 4)], align='h', ncol=4, rel_widths=1)
		# ggsave(paste0('./Figures & Tables/Distributions of Occupancy Variables for ', occWindow, '-yr Window by Region Set 1.png'), width=8, height=10, units='in')

		# main <- plot_grid(plotlist=figs[((7 * 4) + 1):length(figs)], align='h', ncol=4, rel_widths=1)
		# ggsave(paste0('./Figures & Tables/Distributions of Occupancy Variables for ', occWindow, '-yr Window by Region Set 2.png'), width=8, height=10, units='in')

	# } # next occupancy window

# say('#################################################################################')
# say('### distributions of predictors by region and ordinal occupancy for main text ###')
# say('#################################################################################')

	# titleSize <- 8
	# subtitleSize <- 7
	# legendTitleSize <- 7
	# legendTextSize <- 7
	# axisLabelSize <- 7
	# axisTextSize <- 6
	
	# legendKeySize <- 0.4
	
	# # occWindow <- 7
	# occWindow <- 10

	# ### data
	# load('./Data/03 New Mexico Pika - Assigned Folds.rda')

	# ### occupancy predictors
	# ########################
	
	# pika$latestOccStatus <- factor(pika$latestOccStatus, levels=c('0 never', '1 old', '2 occupied'))
	
	# preds <- c('chronicCold_C', 'chronicHeat_C', 'subLethalHeat18deg_d', 'gsPpt_mm')

	# figs <- list()
	# for (countPred in seq_along(preds)) {
	
		# pred <- preds[countPred]
		
		# predWindow <- paste0('occVar_', pred, '_', occWindow, 'yrWindow')
		# say(predWindow)
	
		# predIndex <- which(predTable$var == pred)
	
		# predNice <- predTable$varNice[predIndex]
		# predNice <- capIt(predNice)
		# predDescriptorUnit <- paste0(predTable$unitDescriptor[predIndex], ' (', predTable$unit[predIndex], ')')
	
		# # mus <- data.frame(
			# # latestOccStatus = c('0 never', '1 old', '2 occupied'),
			# # mu = c(
				# # mean(pika[pika$latestOccStatus == '0 never', predWindow]),
				# # mean(pika[pika$latestOccStatus == '1 old', predWindow]),
				# # mean(pika[pika$latestOccStatus == '2 occupied', predWindow])
			# # )
		# # )

		# thisData <- pika[ , c('latestOccStatus', 'region', predWindow)]
		# names(thisData)[3] <- 'value'

		# xlim <- range(thisData$value)
		
		# letter <- letters[countPred]
		# letter <- paste0('(', letter, ') ')
		
		# # all regions together
		# title <- paste0(letter, predNice)
		# figs[[length(figs) + 1]] <- ggplot(data=thisData, aes(x=value, col=latestOccStatus, fill=latestOccStatus)) +
			# geom_density(linewidth=1) +
			# scale_color_manual(
				# labels = c('none', 'old', 'occ'),
				# values=c('0 never'='firebrick3', '1 old'='darkgoldenrod3', '2 occupied'='darkgreen')
			# ) +
			# scale_fill_manual(
				# labels = c('none', 'old', 'occ'),
				# values=alpha(c('0 never'='firebrick3', '1 old'='darkgoldenrod3', '2 occupied'='darkgreen'), 0.2)
			# ) +
			# labs(title=title, subtitle='Regions Together', x=predDescriptorUnit, y='Density') +
			# # geom_vline(data=mus, aes(xintercept=mu, color=latestOccStatus), linetype='dotted', size=1) +
			# guides(
				# color=guide_legend(title='Evidence'),
				# fill=guide_legend(title='Evidence')
			# ) +
			# xlim(xlim[1], xlim[2]) +
			# theme(
				# legend.key.size = unit(legendKeySize, 'cm'),
				# plot.title=element_text(size=titleSize, face='bold'),
				# plot.subtitle=element_text(size=subtitleSize),
				# legend.title=element_text(size=legendTitleSize),
				# legend.text=element_text(size=legendTextSize),
				# axis.title=element_text(size=axisLabelSize),
				# axis.text=element_text(size=axisTextSize)
			# )

		# # "never" by region
		# title <- paste0('')
		# thisThisData <- thisData[thisData$latestOccStatus == '0 never', ]
		# figs[[length(figs) + 1]] <- ggplot(data=thisThisData, aes(x=value, col=region, fill=region)) +
			# geom_density(linewidth=1) +
			# scale_color_manual(
				# labels=c('southwest'='SW', 'southeast'='SE', 'northwest'='NW', 'northeast'='NE'),
				# values=c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred')
			# ) +
			# scale_fill_manual(
				# labels=c('southwest'='SW', 'southeast'='SE', 'northwest'='NW', 'northeast'='NE'),
				# values=alpha(c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred'), 0.2),
			# ) +
			# labs(title=title, subtitle='No Evidence', x=predDescriptorUnit, y='Density') +
			# guides(
				# color=guide_legend(title='Region'),
				# fill=guide_legend(title='Region')
			# ) +
			# xlim(xlim[1], xlim[2]) +
			# theme(
				# legend.key.size = unit(legendKeySize, 'cm'),
				# plot.title=element_text(size=titleSize, face='bold'),
				# plot.subtitle=element_text(size=subtitleSize),
				# legend.title=element_text(size=legendTitleSize),
				# legend.text=element_text(size=legendTextSize),
				# axis.title=element_text(size=axisLabelSize),
				# axis.text=element_text(size=axisTextSize)
			# )

		# # "old" by region
		# title <- paste0('')
		# thisThisData <- thisData[thisData$latestOccStatus == '1 old', ]
		# figs[[length(figs) + 1]] <- ggplot(data=thisThisData, aes(x=value, col=region, fill=region)) +
			# geom_density(linewidth=1) +
			# scale_color_manual(
				# labels=c('southwest'='SW', 'southeast'='SE', 'northwest'='NW', 'northeast'='NE'),
				# values=c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred')
			# ) +
			# scale_fill_manual(
				# labels=c('southwest'='SW', 'southeast'='SE', 'northwest'='NW', 'northeast'='NE'),
				# values=alpha(c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred'), 0.2),
			# ) +
			# labs(title=title, subtitle='Previously Occupied', x=predDescriptorUnit, y='Density') +
			# guides(
				# color=guide_legend(title='Region'),
				# fill=guide_legend(title='Region')
			# ) +
			# xlim(xlim[1], xlim[2]) +
			# theme(
				# legend.key.size = unit(legendKeySize, 'cm'),
				# plot.title=element_text(size=titleSize, face='bold'),
				# plot.subtitle=element_text(size=subtitleSize),
				# legend.title=element_text(size=legendTitleSize),
				# legend.text=element_text(size=legendTextSize),
				# axis.title=element_text(size=axisLabelSize),
				# axis.text=element_text(size=axisTextSize)
			# )

		# # "occupied" by region
		# title <- paste0('')
		# thisThisData <- thisData[thisData$latestOccStatus == '2 occupied', ]
		# figs[[length(figs) + 1]] <- ggplot(data=thisThisData, aes(x=value, col=region, fill=region)) +
			# geom_density(linewidth=1) +
			# scale_color_manual(
				# labels=c('southwest'='SW', 'southeast'='SE', 'northwest'='NW', 'northeast'='NE'),
				# values=c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred')
			# ) +
			# scale_fill_manual(
				# labels=c('southwest'='SW', 'southeast'='SE', 'northwest'='NW', 'northeast'='NE'),
				# values=alpha(c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred'), 0.2),
			# ) +
			# labs(title=title, subtitle='Currently Occupied', x=predDescriptorUnit, y='Density') +
			# guides(
				# color=guide_legend(title='Region'),
				# fill=guide_legend(title='Region')
			# ) +
			# xlim(xlim[1], xlim[2]) +
			# theme(
				# legend.key.size = unit(legendKeySize, 'cm'),
				# plot.title=element_text(size=titleSize, face='bold'),
				# plot.subtitle=element_text(size=subtitleSize),
				# legend.title=element_text(size=legendTitleSize),
				# legend.text=element_text(size=legendTextSize),
				# axis.title=element_text(size=axisLabelSize),
				# axis.text=element_text(size=axisTextSize)
			# )
			
	# } # next predictor

	# main <- plot_grid(plotlist=figs, align='h', ncol=4, rel_widths=1, labels=NULL, label_size=12)
	# ggsave(paste0('./Figures & Tables/Distributions of Occupancy Variables for ', occWindow, '-yr Window by Region MAIN TEXT Ordinal.pdf'), width=8, height=length(preds) * 1.8, units='in')

# say('################################################################################')
# say('### distributions of predictors by region and binary occupancy for main text ###')
# say('################################################################################')

	# titleSize <- 8
	# subtitleSize <- 7
	# legendTitleSize <- 6.5
	# legendTextSize <- 6.5
	# axisLabelSize <- 7
	# axisTextSize <- 6
	
	# legendKeySize <- 0.4
	
	# lw <- 0.5 # width of lines for density smoother
	# occWindow <- 10

	# ### data
	# load('./Data/03 New Mexico Pika - Assigned Folds.rda')

	# ### occupancy predictors
	# ########################
	
	# pika$latestOccStatus <- factor(pika$latestOccStatus, levels=c('0 never', '1 old', '2 occupied'))
	# pika$presAbs <- ifelse(pika$latestOccStatus %in% c('0 never', '1 old'), 'unoccupied', 'occupied')

	# preds <- c('chronicCold_C', 'chronicHeat_C', 'subLethalHeat18deg_d', 'gsPpt_mm')

	# figs <- list()
	# for (countPred in seq_along(preds)) {
	
		# pred <- preds[countPred]
		
		# predWindow <- paste0('occVar_', pred, '_', occWindow, 'yrWindow')
		# say(predWindow)
	
		# predIndex <- which(predTable$var == pred)
	
		# predNice <- predTable$varNice[predIndex]
		# predNice <- capIt(predNice)
		# predDescriptorUnit <- paste0(predTable$unitDescriptor[predIndex], ' (', predTable$unit[predIndex], ')')
	
		# # mus <- data.frame(
			# # latestOccStatus = c('0 never', '1 old', '2 occupied'),
			# # mu = c(
				# # mean(pika[pika$latestOccStatus == '0 never', predWindow]),
				# # mean(pika[pika$latestOccStatus == '1 old', predWindow]),
				# # mean(pika[pika$latestOccStatus == '2 occupied', predWindow])
			# # )
		# # )

		# thisData <- pika[ , c('presAbs', 'region', predWindow)]
		# names(thisData)[3] <- 'value'

		# xlim <- range(thisData$value)
		
		# letter <- letters[countPred]
		# letter <- paste0('(', letter, ') ')
		
		# # all regions together
		# title <- paste0(letter, predNice)
		# figs[[length(figs) + 1]] <- ggplot(data=thisData, aes(x=value, col=presAbs, fill=presAbs)) +
			# geom_density(linewidth=lw) +
			# scale_color_manual(
				# labels = c('occ.', 'unocc.'),
				# values=c('unoccupied'='firebrick3', 'occupied'='darkgreen')
			# ) +
			# scale_fill_manual(
				# labels = c('occ.', 'unocc.'),
				# values=alpha(c('unoccupied'='firebrick3', 'occupied'='darkgreen'), 0.2)
			# ) +
			# labs(title=title, subtitle='Regions Together', x=predDescriptorUnit, y='Density') +
			# # geom_vline(data=mus, aes(xintercept=mu, color=latestOccStatus), linetype='dotted', size=1) +
			# guides(
				# color=guide_legend(title='Status'),
				# fill=guide_legend(title='Status')
			# ) +
			# xlim(xlim[1], xlim[2]) +
			# theme(
				# legend.key.size = unit(legendKeySize, 'cm'),
				# plot.title=element_text(size=titleSize, face='bold'),
				# plot.subtitle=element_text(size=subtitleSize),
				# legend.title=element_text(size=legendTitleSize),
				# legend.text=element_text(size=legendTextSize),
				# axis.title=element_text(size=axisLabelSize),
				# axis.text=element_text(size=axisTextSize)
			# )

		# # "unoccupied" by region
		# title <- paste0('')
		# thisThisData <- thisData[thisData$presAbs == 'unoccupied', ]
		# figs[[length(figs) + 1]] <- ggplot(data=thisThisData, aes(x=value, col=region, fill=region)) +
			# geom_density(linewidth=lw) +
			# scale_color_manual(
				# labels=c('southwest'='SW', 'southeast'='SE', 'northwest'='NW', 'northeast'='NE'),
				# values=c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred')
			# ) +
			# scale_fill_manual(
				# labels=c('southwest'='SW', 'southeast'='SE', 'northwest'='NW', 'northeast'='NE'),
				# values=alpha(c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred'), 0.2),
			# ) +
			# labs(title=title, subtitle='Unoccupied', x=predDescriptorUnit, y='Density') +
			# guides(
				# color=guide_legend(title='Region'),
				# fill=guide_legend(title='Region')
			# ) +
			# xlim(xlim[1], xlim[2]) +
			# theme(
				# legend.key.size = unit(legendKeySize, 'cm'),
				# plot.title=element_text(size=titleSize, face='bold'),
				# plot.subtitle=element_text(size=subtitleSize),
				# legend.title=element_text(size=legendTitleSize),
				# legend.text=element_text(size=legendTextSize),
				# axis.title=element_text(size=axisLabelSize),
				# axis.text=element_text(size=axisTextSize)
			# )

		# # "occupied" by region
		# title <- paste0('')
		# thisThisData <- thisData[thisData$presAbs == 'occupied', ]
		# figs[[length(figs) + 1]] <- ggplot(data=thisThisData, aes(x=value, col=region, fill=region)) +
			# geom_density(linewidth=lw) +
			# scale_color_manual(
				# labels=c('southwest'='SW', 'southeast'='SE', 'northwest'='NW', 'northeast'='NE'),
				# values=c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred')
			# ) +
			# scale_fill_manual(
				# labels=c('southwest'='SW', 'southeast'='SE', 'northwest'='NW', 'northeast'='NE'),
				# values=alpha(c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred'), 0.2),
			# ) +
			# labs(title=title, subtitle='Currently Occupied', x=predDescriptorUnit, y='Density') +
			# guides(
				# color=guide_legend(title='Region'),
				# fill=guide_legend(title='Region')
			# ) +
			# xlim(xlim[1], xlim[2]) +
			# theme(
				# legend.key.size = unit(legendKeySize, 'cm'),
				# plot.title=element_text(size=titleSize, face='bold'),
				# plot.subtitle=element_text(size=subtitleSize),
				# legend.title=element_text(size=legendTitleSize),
				# legend.text=element_text(size=legendTextSize),
				# axis.title=element_text(size=axisLabelSize),
				# axis.text=element_text(size=axisTextSize)
			# )
			
	# } # next predictor

	# main <- plot_grid(plotlist=figs, align='h', ncol=3, rel_widths=1, labels=NULL, label_size=12)
	# ggsave(paste0('./Figures & Tables/Distributions of Occupancy Variables for ', occWindow, '-yr Window by Region MAIN TEXT Binary.pdf'), width=7, height=length(preds) * 1.8, units='in')

say('############################################################################################')
say('### distributions of predictors by region and binary and ordinal occupancy for main text ###')
say('############################################################################################')

	occWindow <- 10

	titleSize <- 8
	subtitleSize <- 7
	legendTitleSize <- 6.5
	legendTextSize <- 6.5
	axisLabelSize <- 7
	axisTextSize <- 6
	
	legendKeySize <- 0.4
	
	lw <- 0.5 # width of lines for density smoother

	### data
	load('./Data/03 New Mexico Pika - Assigned Folds.rda')

	### occupancy predictors
	########################
	
	pika$latestOccStatus <- factor(pika$latestOccStatus, levels=c('0 never', '1 old', '2 occupied'))
	pika$presAbs <- ifelse(pika$latestOccStatus %in% c('0 never', '1 old'), 'unoccupied', 'occupied')

	preds <- c('chronicCold_C', 'chronicHeat_C', 'subLethalHeat18deg_d', 'gsPpt_mm')

	### BINARY
	ordinals <- binaries <- list()
	for (countPred in seq_along(preds)) {
	
		pred <- preds[countPred]
		
		predWindow <- paste0('occVar_', pred, '_', occWindow, 'yrWindow')
		say(predWindow)
	
		predIndex <- which(predTable$var == pred)
	
		predNice <- predTable$varNice[predIndex]
		predNice <- capIt(predNice)
		predDescriptorUnit <- paste0(predTable$unitDescriptor[predIndex], ' (', predTable$unit[predIndex], ')')
	
		thisData <- pika[ , c('presAbs', 'latestOccStatus', 'region', predWindow)]
		names(thisData)[ncol(thisData)] <- 'value'

		xlim <- range(thisData$value)
		
		letter <- letters[countPred]
		letter <- paste0('(', letter, ') ')
		
		letterPlus <- letters[countPred + length(preds)]
		letterPlus <- paste0('(', letterPlus, ') ')
		
		# binary
		title <- paste0(letter, predNice, ' (Binary)')
		binaries[[length(binaries) + 1]] <- ggplot(data=thisData, aes(x=value, col=presAbs, fill=presAbs)) +
			geom_density(linewidth=lw) +
			scale_color_manual(
				labels = c('occupied', 'unoccupied'),
				values=c('unoccupied'='firebrick3', 'occupied'='darkgreen')
			) +
			scale_fill_manual(
				labels = c('occupied', 'unoccupied'),
				values=alpha(c('unoccupied'='firebrick3', 'occupied'='darkgreen'), 0.2)
			) +
			labs(title=title, x=predDescriptorUnit, y='Density') +
			guides(
				color=guide_legend(title='Status'),
				fill=guide_legend(title='Status')
			) +
			xlim(xlim[1], xlim[2]) +
			theme(
				legend.key.size = unit(legendKeySize, 'cm'),
				plot.title=element_text(size=titleSize, face='bold'),
				plot.subtitle=element_text(size=subtitleSize),
				legend.title=element_text(size=legendTitleSize),
				legend.text=element_text(size=legendTextSize),
				axis.title=element_text(size=axisLabelSize),
				axis.text=element_text(size=axisTextSize)
			)

		# ordinal
		title <- paste0(letterPlus, predNice, ' (Ordinal)')
		ordinals[[length(ordinals) + 1]] <- ggplot(data=thisData, aes(x=value, col=latestOccStatus, fill=latestOccStatus)) +
			geom_density(linewidth=lw) +
			scale_color_manual(
				labels = c('0 never' = 'no evidence', '1 old' = 'previously\n  occupied', '2 occupied' = 'currently\n  occupied'),
				values=c('0 never'='firebrick3', '1 old' = 'darkgoldenrod3', '2 occupied' = 'darkgreen')
			) +
			scale_fill_manual(
				labels = c('0 never' = 'no evidence', '1 old' = 'previously\n  occupied', '2 occupied' = 'currently\n  occupied'),
				values=alpha(c('0 never'='firebrick3', '1 old' = 'darkgoldenrod3', '2 occupied' = 'darkgreen'), 0.2)
			) +
			labs(title=title, x=predDescriptorUnit, y='Density') +
			guides(
				color=guide_legend(title='Status'),
				fill=guide_legend(title='Status')
			) +
			xlim(xlim[1], xlim[2]) +
			theme(
				legend.key.size = unit(legendKeySize, 'cm'),
				plot.title=element_text(size=titleSize, face='bold'),
				plot.subtitle=element_text(size=subtitleSize),
				legend.title=element_text(size=legendTitleSize),
				legend.text=element_text(size=legendTextSize),
				axis.title=element_text(size=axisLabelSize),
				axis.text=element_text(size=axisTextSize)
			)
			
	} # next predictor

	binaries <- plot_grid(plotlist = binaries, align='v', ncol=1)
	ordinals <- plot_grid(plotlist = ordinals, align='v', ncol=1)
	main <- binaries + ordinals
	
	ggsave(paste0('./Figures & Tables/Distributions of Occupancy Variables for ', occWindow, '-yr Window by Region MAIN TEXT Binary & Ordinal.pdf'), width=6, height = 8, units='in')

# say('###########################################################')
# say('### contingency table analysis of site status by region ###')
# say('###########################################################')

	# ### data
	# load('./Data/03 New Mexico Pika - Assigned Folds.rda')
	
	# cont <- table(pika$latestOccStatus, pika$region)
	# regions <- c('NE', 'NW', 'SE', 'SW')
	# statuses <- c('None', 'Previous', 'Occupied')
	# colnames(cont) <- regions
	# rownames(cont) <- statuses
	
	# contPrime <- t(cont)
	# contPrime <- contPrime[ , 3:1]
	
	# propPerRegion <- rowSums(contPrime) / sum(contPrime)
	# regionCenters <- propPerRegion[1] / 2
	# regionCenters <- c(
		# propPerRegion[1] / 2,
		# propPerRegion[1] + propPerRegion[2] / 2,
		# sum(propPerRegion[1:2]) + propPerRegion[3] / 2,
		# sum(propPerRegion[1:3]) + propPerRegion[4] / 2
	# )
	
	# yCont <- contPrime
	# sitesPerRegion <- rowSums(contPrime)
	# for (i in seq_along(regions)) yCont[i, ] <- yCont[i, ] / sitesPerRegion[i]
	# yContFull <- yCont
	# yCont[ , 'None'] <- yContFull[ , 'None'] / 2
	# yCont[ , 'Previous'] <- yContFull[ , 'None'] + yContFull[ , 'Previous'] / 2
	# yCont[ , 'Occupied'] <- yContFull[ , 'None'] + yContFull[ , 'Previous'] + yContFull[ , 'Occupied'] / 2
	
	# png('./Figures & Tables/Sites by Region.png', width=1200, height=900, res=600)
	
		# par(oma = c(0, 0, 0, 0), mar = c(2, 2, 1, 1), cex = 0.6, lwd = 0.6)
		# mosaicplot(
			# contPrime,
			# xlab='Region',
			# ylab='Evidence',
			# color=c('chartreuse3', 'darkgoldenrod3', 'firebrick'),
			# main=''
		# )
		
		# # frequencies
		# for (r in seq_along(regions)) {
		
			# x <- regionCenters[r]
			# for (s in seq_along(statuses)) {
				
				# y <- yCont[r, s]
				# n <- contPrime[r, s]
				# text(x, y, labels = n, adj = c(0.41, 1), cex = 0.8)
				
			# }
		
		# }
	
	# dev.off()
	
	# sink('./Figures & Tables/Sites by Region Risk Ratio Text.txt', split=TRUE)
		# say('Risk ratio test of site status by region:', post=2)
		# epitools::riskratio(cont, method = 'wald')
	# sink()

	# write.csv(cont, './Figures & Tables/Sites by Region.csv')
	
# say('########################################################################')
# say('### statistical comparison of environmental variables across regions ###')
# say('########################################################################')

	# load('./Data/05 New Mexico Pika - Added PRISM Cell Number & Cell-Based Weight.rda')
	
	# ### occupancy variables
	# vars <- getVars('occupancy')
	# results <- data.frame()
	# for (i in seq_along(vars)) {
	
		# var <- vars[i]
		
		# subPika <- pika[ , c('region', var)]
		# names(subPika)[2] <- 'y'
		# kw <- kruskal.test(y ~ region, data = subPika)
	
		# results <- rbind(
			# results,
			# data.frame(
				# type = 'occupancy',
				# variable = var,
				# varNice = makeNiceVars(var, occOrDens = 'occupancy'),
				# KWChiSquared = kw$statistic,
				# df = kw$parameter,
				# p = kw$p.value
			# )
		
		# )
	
	# }
	
	# ### density variables
	# vars <- getVars('density')
	# for (i in seq_along(vars)) {
	
		# var <- vars[i]
		
		# subPika <- pika[ , c('region', var)]
		# names(subPika)[2] <- 'y'
		# kw <- kruskal.test(y ~ region, data = subPika)
	
		# results <- rbind(
			# results,
			# data.frame(
				# type = 'density',
				# variable = var,
				# varNice = makeNiceVars(var, occOrDens = 'density'),
				# KWChiSquared = kw$statistic,
				# df = kw$parameter,
				# p = kw$p.value
			# )
		
		# )
	
	# }

	# # isolation: occupancy sites
	# subPika <- pika[ , c('region', 'meanDistToClosest4Patches')]
	# names(subPika)[2] <- 'y'
	# kw <- kruskal.test(y ~ region, data = subPika)
		# results <- rbind(
			# results,
			# data.frame(
				# type = 'occupancy (i.e., all sites)',
				# variable = 'meanDistToClosest4Patches',
				# varNice = 'mean distance to closest 4 patches',
				# KWChiSquared = kw$statistic,
				# df = kw$parameter,
				# p = kw$p.value
			# )
		
		# )

	# # isolation: density sites
	# subPika <- pika[!is.na(pika$latestDensity) , c('region', 'meanDistToClosest4Patches')]
	# names(subPika)[2] <- 'y'
	# kw <- kruskal.test(y ~ region, data = subPika)
		# results <- rbind(
			# results,
			# data.frame(
				# type = 'density (i.e., density > 0)',
				# variable = 'meanDistToClosest4Patches',
				# varNice = 'mean distance to closest 4 patches',
				# KWChiSquared = kw$statistic,
				# df = kw$parameter,
				# p = kw$p.value
			# )
		
		# )

	# # number of home ranges: occupancy sites
	# subPika <- pika[ , c('region', 'numHomeRanges')]
	# names(subPika)[2] <- 'y'
	# kw <- kruskal.test(y ~ region, data = subPika)
		# results <- rbind(
			# results,
			# data.frame(
				# type = 'occupancy (i.e., all sites)',
				# variable = 'numHomeRanges',
				# varNice = 'number of home ranges',
				# KWChiSquared = kw$statistic,
				# df = kw$parameter,
				# p = kw$p.value
			# )
		
		# )

	# # number of home ranges: density sites
	# subPika <- pika[!is.na(pika$latestDensity) , c('region', 'numHomeRanges')]
	# names(subPika)[2] <- 'y'
	# kw <- kruskal.test(y ~ region, data = subPika)
		# results <- rbind(
			# results,
			# data.frame(
				# type = 'density (i.e., density > 0)',
				# variable = 'numHomeRanges',
				# varNice = 'number of home ranges',
				# KWChiSquared = kw$statistic,
				# df = kw$parameter,
				# p = kw$p.value
			# )
		
		# )

	# rownames(results) <- NULL

	# write.csv(results, './Figures & Tables/Kruskal-Wallis Test for Each Environmental Variable by Region.csv', row.names = FALSE)
	
	
say('DONE!!!', level=1, deco='%')
