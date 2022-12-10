# source('C:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/TEMP.r')

	### data
	load('./Data/03 New Mexico Pika - Assigned Folds.rda')
	
	cont <- table(pika$latestOccStatus, pika$region)
	regions <- c('NE', 'NW', 'SE', 'SW')
	statuses <- c('None', 'Previous', 'Occupied')
	colnames(cont) <- regions
	rownames(cont) <- statuses
	
	contPrime <- t(cont)
	contPrime <- contPrime[ , 3:1]
	
	propPerRegion <- rowSums(contPrime) / sum(contPrime)
	regionCenters <- propPerRegion[1] / 2
	regionCenters <- c(
		propPerRegion[1] / 2,
		propPerRegion[1] + propPerRegion[2] / 2,
		sum(propPerRegion[1:2]) + propPerRegion[3] / 2,
		sum(propPerRegion[1:3]) + propPerRegion[4] / 2
	)
	
	yCont <- contPrime
	sitesPerRegion <- rowSums(contPrime)
	for (i in seq_along(regions)) yCont[i, ] <- yCont[i, ] / sitesPerRegion[i]
	yContFull <- yCont
	yCont[ , 'None'] <- yContFull[ , 'None'] / 2
	yCont[ , 'Previous'] <- yContFull[ , 'None'] + yContFull[ , 'Previous'] / 2
	yCont[ , 'Occupied'] <- yContFull[ , 'None'] + yContFull[ , 'Previous'] + yContFull[ , 'Occupied'] / 2
	
	png('./Figures & Tables/Sites by Region.png', width=1200, height=900, res=600)
	
		par(oma = c(0, 0, 0, 0), mar = c(2, 2, 1, 1), cex = 0.6, lwd = 0.6)
		mosaicplot(
			contPrime,
			xlab='Region',
			ylab='Evidence',
			color=c('chartreuse3', 'darkgoldenrod3', 'firebrick'),
			main=''
		)
		
		# frequencies
		for (r in seq_along(regions)) {
		
			x <- regionCenters[r]
			for (s in seq_along(statuses)) {
				
				y <- yCont[r, s]
				n <- contPrime[r, s]
				text(x, y, labels = n, adj = c(0.41, 1), cex = 0.8)
				
			}
		
		}
	
	dev.off()
