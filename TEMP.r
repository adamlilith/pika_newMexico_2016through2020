# source('C:/Ecology/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/TEMP.r')

		### map of change in future suitability
		#######################################

			# color ranges
			ranges <- c(Inf, -Inf)
			for (ssp in ssps) {
				for (period in periods) {
					mm <- minmax(deltaPredictions[[paste(ssp, period)]])
					ranges[1] <- min(ranges[1], mm[1, ])
					ranges[2] <- max(ranges[2], mm[2, ])
				}
			}

			colNeg <- colorRampPalette(c('red', 'gray90'))
			colPos <- colorRampPalette(c('gray90', 'green3'))
			negRanges <- abs(100 * round(ranges[1], 2))
			posRanges <- 100 * round(ranges[2], 2)
			
			maxRanges <- max(negRanges, posRanges)
			colNeg <- colNeg(maxRanges)
			colPos <- colPos(maxRanges)
			
			colNeg <- colNeg[1:negRanges]
			colPos <- colPos[1:posRanges]

			deltaCol <- c(colNeg, colPos)

			png(paste0('./Figures & Tables/Maps of Predictions/Binary Occurrence Futures.png'), res = 600, width = 2 * 1200, height = 2 * 1300)
			
				par(mfrow = c(2, 2), mar = c(0, 0, 0, 0))

				i <- 1
				for (ssp in ssps) {
				
					for (period in periods) {

						nicePeriod <- if (period == '2041_2070') { '2050s' } else { '2080s' }
						
						extent <- as.vector(ext(studyRegionFocus))
						xlim <- extent[1:2]
						ylim <- extent[3:4]

						if (i %in% 1:2) { mar <- c(0, 0.1, 0.1, 0.1) } else { mar <- c(0.1, 0.1, 0, 0.1) }

						main <- paste0('(', letters[i + 1], ') SSP ', ssp, ': ', nicePeriod)
						col <- paste0('gray', 30:70)
						plot(hs, col = col, legend = FALSE, axes = FALSE, ext = studyRegionFocus, box = FALSE, maxcell = maxcell, mar = c(0.4, 0.1, 0.1, 0.1))
						# plot(sqPredictions[['northwest']], legend = FALSE, axes = FALSE, ext = studyRegionFocus, box = FALSE, maxcell = maxcell, mar = mar)
						
						for (region in regions) {
						
							x <- deltaPredictions[[paste(ssp, period)]][[region]]
							
							plot(x, col = deltaCol, alpha = 0.65, legend = FALSE, axes = FALSE, ext = studyRegionFocus, box = FALSE, maxcell = maxcell, range = ranges, add = TRUE)
							
						}
						
						plot(focusCounties, lwd = 0.3, add = TRUE)
						plot(studyRegionFocus, add = TRUE, lwd = 0.7)

						if (i == 4) {

							legendGrad(
								x = 'bottom',
								inset = -0.02,
								vert = FALSE,
								width = 1,
								height = 0.06,
								adjX = c(-0.8, 1),
								adjY = c(-0.1, 1),
								col = deltaCol,
								title = 'Change',
								titleCex = 0.7,
								titleAdj = c(-0.92, 0.45),
								labels = '',
								boxBorder = NA,
								lwd = 0.8,
								xpd = NA
							)
							
							niceRanges <- roundTo(ranges, 0.01)
							
							y <- -1065937
							text(-825000, y, labels = niceRanges[1], cex = 0.63, xpd = NA)
							text(-575000, y, labels = 0, cex = 0.63, xpd = NA)
							text(-440000, y, labels = paste0('+', niceRanges[2]), cex = 0.63, xpd = NA)
							
						}
						
						labelFig(label = main, adj = c(-0.05, -0.015), cex = 0.75)
						
						i <- i + 1
						
					} # next period
				
				} # next ssp

			dev.off()
