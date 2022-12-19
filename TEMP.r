# source('C:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/TEMP.r')

	### colors
	
		# colors for elevation
		load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')
		
		minElev <- globalx(elev, 'min')
		lowestNeverOcc <- min(pika$elevation_m[pika$latestOccStatus == '0 never'])
		# medianNeverOcc <- median(pika$elevation_m[pika$latestOccStatus == '0 never'])
		medianPastOcc <- median(pika$elevation_m[pika$latestOccStatus == '1 old'])
		medianOcc <- median(pika$elevation_m[pika$latestOccStatus == '2 occupied'])
		maxElev <- globalx(elev, 'max')

		elevBreaks <- c(minElev, lowestNeverOcc, medianPastOcc, medianOcc, maxElev)
		names(elevBreaks) <- c('minElev', 'medianNeverOcc', 'medianPastOcc', 'medianOcc', 'maxElev')

		elevCols <- c('gray80', 'firebrick3', 'darkgoldenrod3', 'chartreuse')
		elevCols <- alpha(elevCols, 0.5)

		# hillshade colors
		hsCols <- colorRampPalette(c('gray0', 'gray100'))
		hsCols <- hsCols(20)
		
	### placement
	
		legendInset <- c(0.019, 0.02)
		mar <- c(0, 0.5, 1, 1)
		
	### plot!
	png('./Figures & Tables/Study Region with Sampling Sites.png', width=2300, height=2400, res=300)

		par(mfrow=c(2, 2), mar=rep(0, 4), mai=rep(0, 4), oma=rep(0, 4), mgp=c(0, 0, 0))

		### occurrences
		###############
			
			plot(hs, col=hsCols, legend=FALSE, axes=FALSE, mar=mar)
			plot(nam1LgCrop, lwd=3, add=TRUE)
			plot(usa2LgCrop, lwd=1, add=TRUE)
			plot(elev, legend = FALSE, axes = FALSE, add = TRUE, col = elevCols, breaks = elevBreaks)

			plot(occs, pch=21, bg=alpha('chartreuse', 0.5), cex=1.2, add=TRUE)

			legend('bottomleft', inset=legendInset, legend=c('Currently occupied'), pch=c(21), pt.bg=c( 'chartreuse'), bg='white', cex=1.05)
			
		### scale bar
		#############

			usr <- par('usr')
			length <- 50 # in km
			nudge <- 12000
			x <- c(usr[2] - length * 1000 - nudge, usr[2] - nudge)
			y <- rep(usr[3] + 0.03 * (usr[4] - usr[3]), 2)
			lines(x, y, lwd=9, lend=1)
			y <- usr[3] + 0.07 * (usr[4] - usr[3])
			text(mean(x), y, labels=paste(length, 'km'), cex=1.1)
			
		### old evidence
		################
			
			plot(hs, col=hsCols, legend=FALSE, axes=FALSE, mar=mar)
			plot(nam1LgCrop, lwd=3, add=TRUE)
			plot(usa2LgCrop, lwd=1, add=TRUE)
			plot(elev, legend = FALSE, axes = FALSE, add = TRUE, col = elevCols, breaks = elevBreaks)

			plot(oldEvid, pch=22, bg=alpha('darkgoldenrod3', 0.5), cex=1.2, add=TRUE)

			legend('bottomleft', inset=legendInset, legend=c('Previously occupied'), pch=c(22), pt.bg=c('darkgoldenrod3'), bg='white', cex=1.05)
			
		### no evidence
		###############
			
			plot(hs, col=hsCols, legend=FALSE, axes=FALSE, mar=mar)
			plot(nam1LgCrop, lwd=3, add=TRUE)
			plot(usa2LgCrop, lwd=1, add=TRUE)
			plot(elev, legend = FALSE, axes = FALSE, add = TRUE, col = elevCols, breaks = elevBreaks)

			plot(noEvid, pch=25, bg=alpha('firebrick2', 0.5), cex=1.2, add=TRUE)

			legend('bottomleft', inset=legendInset, legend=c('No evidence'), pch=c(25), pt.bg=c('firebrick2'), bg='white', cex=1.05)
			
		### all together
		################
			
			plot(hs, col=hsCols, legend=FALSE, axes=FALSE, mar=mar)
			plot(nam1LgCrop, lwd=3, add=TRUE)
			plot(usa2LgCrop, lwd=1, add=TRUE)
			plot(elev, legend = FALSE, axes = FALSE, add = TRUE, col = elevCols, breaks = elevBreaks)

			noEvid <- pikaVectProj[pika$latestOccStatus == '0 never']
			oldEvid <- pikaVectProj[pika$latestOccStatus == '1 old']
			occs <- pikaVectProj[pika$latestOccStatus == '2 occupied']
			
			plot(noEvid, pch=25, bg=alpha('firebrick2', 0.5), cex=1.2, add=TRUE)
			plot(oldEvid, pch=22, bg=alpha('darkgoldenrod3', 0.5), cex=1.2, add=TRUE)
			plot(occs, pch=21, bg=alpha('chartreuse', 0.5), cex=1.2, add=TRUE)

			legend('bottomleft', inset=legendInset, legend=c('Currently occupied', 'Previously occupied', 'No evidence'), pch=c(21, 22, 25), pt.bg=c('chartreuse', 'darkgoldenrod3', 'firebrick2'), bg='white', cex=1.05)
			
		### color ramp
		##############
		legendBreaks(
			x = 'bottomright',
			inset = legendInset,
			width = 0.2,
			height = 0.38,
			labels = roundTo(elevBreaks, 10),
			labAdjX = 0.7,
			labAdjY = c(0.05, 0.25, 0.5, 0.75, 0.98),
			col = elevCols,
			colBorder = 'black',
			title = 'Elevation\n(m)',
			titleAdj = c(0.5, 0.88),
			adjX = c(0.1, 0.4),
			adjY = c(0.05, 0.72)
		)
		
		### inset
		#########
			
			# par(fig = c(0.455, 0.715, 0.237, 0.537), bg='white', new=TRUE)
			par(fig = c(0.455, 0.715, 0.237 - 0.01, 0.537 - 0.01), bg='white', new=TRUE)
			
			ext <- ext(nam1XL)
			corners <- ext@ptr$vector
			ext <- raster::extent(corners)
			ext <- as(ext, 'SpatialPolygons')
			projection(ext) <- getCRS('naAlbers')
			ext <- vect(ext)

			nam1XXL <- crop(nam1XXL, ext)
			plot(ext, border='gray', col='gray', axes=FALSE)
			plot(nam1XXL, col='white', border='gray40', ann=FALSE, add=TRUE)
			
			# highlight study region in inset
			pikaBuffProjSmExt <- ext(pikaBuffProjSm)
			corners <- pikaBuffProjSmExt@ptr$vector
			focus <- raster::extent(corners)
			focus <- as(focus, 'SpatialPolygons')
			projection(focus) <- getCRS('naAlbers')
			focus <- vect(focus)
			plot(focus, lwd=1.9, border='black', add=TRUE)
			
			plot(ext, add=TRUE)
			
	dev.off()
