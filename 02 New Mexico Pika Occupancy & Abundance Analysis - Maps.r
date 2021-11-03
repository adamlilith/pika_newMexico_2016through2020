### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/Code/02 New Mexico Pika Occupancy & Abundance Analysis - Maps.r')
###
### CONTENTS ###
### setup ###
### map of sampling sites ###

#############
### setup ###
#############

	source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/Code/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r')

# say('#############################')
# say('### map of sampling sites ###')
# say('#############################')

	# # survey sites
	# load('./Data/02 New Mexico Pika - Environmental Values Extracted and Calculated.rda')
	# pikaVectUnproj <- vect(as.matrix(pika[ , ll]), 'points', crs=getCRS('wgs84'))
	
	# # plot extent
	# pikaBuffUnprojXXL <- terra::buffer(pikaVectUnproj, width=2000000)
	# pikaBuffUnprojXL <- terra::buffer(pikaVectUnproj, width=1000000)
	# pikaBuffUnprojLg <- terra::buffer(pikaVectUnproj, width=200000)
	# pikaBuffUnprojSm <- terra::buffer(pikaVectUnproj, width=40000)
	# pikaBuffProjSm <- project(pikaBuffUnprojSm, getCRS('albersNA'))

	# # GADM
	# mex1 <- getData('GADM', country='MEX', level=1, path=paste0(drive, '/ecology/!Scratch'))
	# usa1 <- getData('GADM', country='USA', level=1, path=paste0(drive, '/ecology/!Scratch'))
	# usa2 <- getData('GADM', country='USA', level=2, path=paste0(drive, '/ecology/!Scratch'))

	# nam1 <- rbind(usa1, mex1)
	# nam1 <- vect(nam1)
	# usa2 <- vect(usa2)

	# nam1XXL <- crop(nam1, pikaBuffUnprojXXL)
	# nam1XL <- crop(nam1, pikaBuffUnprojXL)
	# nam1Lg <- crop(nam1, pikaBuffUnprojLg)
	# usa2Lg <- crop(usa2, pikaBuffUnprojLg)
	
	# # elevation
	# srtm <- rast('F:/ecology/Topography/SRTM - CGIAR/Tiles/cut_n30w120.tif')
	# srtm <- crop(srtm, pikaBuffUnprojLg)
	
	# slope <- terrain(srtm, 'slope', unit='radians')
	# aspect <- terrain(srtm, 'aspect', unit='radians')
	# slopeR <- raster(slope)
	# aspectR <- raster(aspect)

	# hs <- hillShade(slopeR, aspectR, direction=45)
	# hs <- rast(hs)
	# hs <- project(hs, getCRS('albersNA'))
	
	# # project
	# pikaVectProj <- project(pikaVectUnproj, getCRS('albersNA'))
	# nam1XXL <- project(nam1XXL, getCRS('albersNA'))
	# nam1XL <- project(nam1XL, getCRS('albersNA'))
	# nam1Lg <- project(nam1Lg, getCRS('albersNA'))
	# usa2Lg <- project(usa2Lg, getCRS('albersNA'))
	
	# # split records
	# noEvid <- pikaVectProj[pika$latestOccStatus == '0 never']
	# oldEvid <- pikaVectProj[pika$latestOccStatus == '1 old']
	# occs <- pikaVectProj[pika$latestOccStatus == '2 occupied']

	# # plot extent
	# pikaBuffProjSmSp <- as(pikaBuffProjSm, 'Spatial')
	# ext <- extent(pikaBuffProjSmSp)
	# ext <- as(ext, 'SpatialPolygons')
	# projection(ext) <- getCRS('albersNA')

	# usa2LgCrop <- crop(usa2Lg, ext)
	# nam1LgCrop <- crop(nam1Lg, ext)
	# hs <- crop(hs, ext)
	# hs <- hs - global(hs, 'min', na.rm=TRUE)$min
	# hs <- hs / global(hs, 'max', na.rm=TRUE)$max
	# hsR <- raster(hs)

	# ### plot!
	# png('./Figures & Tables/Study Region with Sampling Sites.png', width=2400, height=2400, res=300)

		# par(mfrow=c(2, 2), mar=rep(0, 4), oma=rep(0, 4), mgp=c(0, 0, 0))
		
		# ### occurrences
		# ###############
			
			# plot(pikaBuffProjSmSp, border='white', ann=FALSE)
			# plot(hsR, col=paste0('gray', 0:100), legend=FALSE, ann=FALSE, add=TRUE)
			# plot(nam1LgCrop, lwd=3, add=TRUE)
			# plot(usa2LgCrop, lwd=1, add=TRUE)

			# plot(occs, pch=21, bg=alpha('chartreuse', 0.5), cex=1.2, add=TRUE)

			# legend('bottomleft', inset=c(0.045, 0.06), legend=c('Occupied'), pch=c(21), pt.bg=c( 'chartreuse'), bg='white', cex=1.05)
			
		# ### old evidence
		# ################
			
			# plot(pikaBuffProjSmSp, border='white', ann=FALSE)
			# plot(hsR, col=paste0('gray', 0:100), legend=FALSE, ann=FALSE, add=TRUE)
			# plot(nam1LgCrop, lwd=3, add=TRUE)
			# plot(usa2LgCrop, lwd=1, add=TRUE)

			# plot(oldEvid, pch=22, bg=alpha('darkgoldenrod3', 0.5), cex=1.2, add=TRUE)

			# legend('bottomleft', inset=c(0.045, 0.06), legend=c('Old evidence'), pch=c(22), pt.bg=c('darkgoldenrod3'), bg='white', cex=1.05)
			
		# ### no evidence
		# ###############
			
			# plot(pikaBuffProjSmSp, border='white', ann=FALSE)
			# plot(hsR, col=paste0('gray', 0:100), legend=FALSE, axes=FALSE, ann=FALSE, add=TRUE)
			# plot(nam1LgCrop, lwd=3, add=TRUE)
			# plot(usa2LgCrop, lwd=1, add=TRUE)

			# plot(noEvid, pch=25, bg=alpha('firebrick2', 0.5), cex=1.2, add=TRUE)

			# legend('bottomleft', inset=c(0.045, 0.06), legend=c('No evidence'), pch=c(25), pt.bg=c('firebrick2'), bg='white', cex=1.05)
			
		# ### all together
		# ################
			
			# plot(ext, border='black', ann=FALSE)
			# plot(hsR, col=paste0('gray', 0:100), ann=FALSE, axes=FALSE, legend=FALSE, add=TRUE)
			# plot(nam1LgCrop, lwd=3, add=TRUE)
			# plot(usa2LgCrop, lwd=1, add=TRUE)

			# noEvid <- pikaVectProj[pika$latestOccStatus == '0 never']
			# oldEvid <- pikaVectProj[pika$latestOccStatus == '1 old']
			# occs <- pikaVectProj[pika$latestOccStatus == '2 occupied']
			
			# plot(noEvid, pch=25, bg=alpha('firebrick2', 0.5), cex=1.2, add=TRUE)
			# plot(oldEvid, pch=22, bg=alpha('darkgoldenrod3', 0.5), cex=1.2, add=TRUE)
			# plot(occs, pch=21, bg=alpha('chartreuse', 0.5), cex=1.2, add=TRUE)

			# legend('bottomleft', inset=c(0.045, 0.06), legend=c('Occupied', 'Old evidence', 'No evidence'), pch=c(21, 22, 25), pt.bg=c('chartreuse', 'darkgoldenrod3', 'firebrick2'), bg='white', cex=1.05)
			
			# ### scale bar
			# usr <- par('usr')
			# length <- 50 # in km
			# nudge <- 12000
			# x <- c(usr[2] - length * 1000 - nudge, usr[2] - nudge)
			# y <- rep(usr[3] + 0.07 * (usr[4] - usr[3]), 2)
			# lines(x, y, lwd=9, lend=1)
			# y <- usr[3] + 0.105 * (usr[4] - usr[3])
			# text(mean(x), y, labels=paste(length, 'km'), cex=1.1)
			
			# ### inset
			# par(fig = c(0.455, 0.715, 0.237, 0.537), bg='white', new=TRUE)
			
			# ext <- ext(nam1XL)
			# corners <- ext@ptr$vector
			# ext <- raster::extent(corners)
			# ext <- as(ext, 'SpatialPolygons')
			# projection(ext) <- getCRS('albersNA')
			# ext <- vect(ext)

			# nam1XXL <- crop(nam1XXL, ext)
			# plot(ext, border='gray', col='gray', axes=FALSE)
			# plot(nam1XXL, col='white', border='gray40', ann=FALSE, add=TRUE)
			
			# # highlight study region in inset
			# pikaBuffProjSmExt <- ext(pikaBuffProjSm)
			# corners <- pikaBuffProjSmExt@ptr$vector
			# focus <- raster::extent(corners)
			# focus <- as(focus, 'SpatialPolygons')
			# projection(focus) <- getCRS('albersNA')
			# focus <- vect(focus)
			# plot(focus, lwd=1.9, border='black', add=TRUE)
			
			# plot(ext, add=TRUE)
			
	# dev.off()


say('DONE!!!', level=1, deco='%')
