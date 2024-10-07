### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('C:/Ecology/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/02 New Mexico Pika Occupancy & Abundance Analysis - Maps.r')
### source('E:/Adam/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/02 New Mexico Pika Occupancy & Abundance Analysis - Maps.r')
###
### CONTENTS ###
### setup ###
### map of sampling sites ###
### map of sampling sites & GBIF occurrences ###

#############
### setup ###
#############

	source('C:/Ecology/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r')
	# source('E:/Adam/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r')

# say('##############################')
# say('### fetch elevation raster ###')
# say('##############################')

	# library(elevatr)
	
	# # elevation from AWS Terrain Tiles
	# load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')
	# pika <- vect(pika, geom = ll, crs = getCRS('WGS84'))
	# focus <- buffer(pika, 50000)
	# focus <- ext(focus)
	# focus <- as.polygons(focus, crs = getCRS('WGS84'))
	# focus <- sf::st_as_sf(focus)
	
	# elev_fine_m <- get_elev_raster(focus, z = 10)
	# elev_fine_m <- round(elev_fine_m, 0L)
	# elev_fine_m <- setMinMax(elev_fine_m)
	# names(elev_fine_m) <- 'elev_fine_m'
	
	# writeRaster(elev_fine_m, './Data/elev_fine_m.tif', datatype = 'INT2S', overwrite = TRUE)
	
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
	# pikaBuffProjSm <- project(pikaBuffUnprojSm, getCRS('naAlbers'))

	# # GADM
	# mex1 <- gadm(country='MEX', level=1, path=paste0('C:/!Scratch'), version=4.1, resolution=1)
	# usa1 <- gadm(country='USA', level=1, path=paste0('C:/!Scratch'), version=4.1, resolution=1)
	# usa2 <- gadm(country='USA', level=2, path=paste0('C:/!Scratch'), version=4.1, resolution=1)

	# usa1 <- usa1[usa1$NAME_1 != 'Alaska', ]
	# usa1 <- usa1[usa1$NAME_1 != 'Hawaii', ]

	# usa2 <- usa2[usa2$NAME_1 != 'Alaska', ]
	# usa2 <- usa2[usa2$NAME_1 != 'Hawaii', ]

	# nam1 <- rbind(usa1, mex1)
	
	# nam1XXL <- crop(nam1, ext(pikaBuffUnprojXXL))
	# nam1XL <- crop(nam1, ext(pikaBuffUnprojXL))
	# nam1Lg <- crop(nam1, ext(pikaBuffUnprojLg))
	# usa2Lg <- crop(usa2, ext(pikaBuffUnprojLg))
	
	# # elevation
	# elev_fine_m <- rast('./Data/elev_fine_m.tif')
	# elev_fine_m <- crop(elev_fine_m, pikaBuffUnprojLg)
	
	# slope <- terrain(elev_fine_m, 'slope', unit='radians')
	# aspect <- terrain(elev_fine_m, 'aspect', unit='radians')
	# slopeR <- raster(slope)
	# aspectR <- raster(aspect)

	# hs <- hillShade(slopeR, aspectR, direction=45)
	# hs <- rast(hs)
	# hs <- project(hs, getCRS('naAlbers'))
	
	# # project
	# pikaVectProj <- project(pikaVectUnproj, getCRS('naAlbers'))
	# nam1XXL <- project(nam1XXL, getCRS('naAlbers'))
	# nam1XL <- project(nam1XL, getCRS('naAlbers'))
	# nam1Lg <- project(nam1Lg, getCRS('naAlbers'))
	# usa2Lg <- project(usa2Lg, getCRS('naAlbers'))
	
	# # split records
	# noEvid <- pikaVectProj[pika$latestOccStatus == '0 never']
	# oldEvid <- pikaVectProj[pika$latestOccStatus == '1 old']
	# occs <- pikaVectProj[pika$latestOccStatus == '2 occupied']

	# # plot extent
	# pikaBuffProjSmSp <- as(pikaBuffProjSm, 'Spatial')
	# ext <- extent(pikaBuffProjSmSp)
	# ext <- as(ext, 'SpatialPolygons')
	# projection(ext) <- getCRS('naAlbers')

	# usa2LgCrop <- crop(usa2Lg, ext)
	# nam1LgCrop <- crop(nam1Lg, ext)
	# hs <- crop(hs, ext)
	# hs <- hs - global(hs, 'min', na.rm=TRUE)$min
	# hs <- hs / global(hs, 'max', na.rm=TRUE)$max
	# hs <- crop(hs, pikaBuffProjSmSp)
	# hsR <- rast(hs)
	
	# elev <- project(elev_fine_m, getCRS('naAlbers'))
	# elev <- crop(elev, pikaBuffProjSmSp)

	# ### colors
	
		# # colors for elevation
		# load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')
		
		# minElev <- globalx(elev, 'min')
		# lowestNeverOcc <- min(pika$elevation_m[pika$latestOccStatus == '0 never'])
		# # medianNeverOcc <- median(pika$elevation_m[pika$latestOccStatus == '0 never'])
		# medianPastOcc <- median(pika$elevation_m[pika$latestOccStatus == '1 old'])
		# medianOcc <- median(pika$elevation_m[pika$latestOccStatus == '2 occupied'])
		# maxElev <- globalx(elev, 'max')

		# elevBreaks <- c(minElev, lowestNeverOcc, medianPastOcc, medianOcc, maxElev)
		# elevBreaks <- round(elevBreaks)
		# names(elevBreaks) <- c('minElev', 'medianNeverOcc', 'medianPastOcc', 'medianOcc', 'maxElev')

		# elevCols <- c('gray80', 'firebrick3', 'darkgoldenrod3', 'chartreuse')
		# elevCols <- alpha(elevCols, 0.5)

		# # hillshade colors
		# hsCols <- colorRampPalette(c('gray0', 'gray100'))
		# hsCols <- hsCols(20)
		
	# ### placement
	
		# legendInset <- c(0.019, 0.02)
		# mar <- c(0, 0.5, 1, 1)
		
	# ### plot!
	# png('./Figures & Tables/Study Region with Sampling Sites.png', width=2300, height=2400, res=300)

		# par(mfrow=c(2, 2), mar=rep(0, 4), mai=rep(0, 4), oma=rep(0, 4), mgp=c(0, 0, 0))

		# ### occurrences
		# ###############
			
			# plot(hs, col=hsCols, legend=FALSE, axes=FALSE, mar=mar)
			# plot(nam1LgCrop, lwd=3, add=TRUE)
			# plot(usa2LgCrop, lwd=1, add=TRUE)
			# plot(elev, legend = FALSE, axes = FALSE, add = TRUE, col = elevCols, breaks = elevBreaks)

			# plot(occs, pch=21, bg=alpha('chartreuse', 0.5), cex=1.2, add=TRUE)

			# legend('bottomleft', inset=legendInset, legend=c('Currently occupied'), pch=c(21), pt.bg=c( 'chartreuse'), bg='white', cex=1.05)
			
		# ### scale bar
		# #############

			# usr <- par('usr')
			# length <- 50 # in km
			# nudge <- 12000
			# x <- c(usr[2] - length * 1000 - nudge, usr[2] - nudge)
			# y <- rep(usr[3] + 0.03 * (usr[4] - usr[3]), 2)
			# lines(x, y, lwd=9, lend=1)
			# y <- usr[3] + 0.07 * (usr[4] - usr[3])
			# text(mean(x), y, labels=paste(length, 'km'), cex=1.1)
			
		# ### old evidence
		# ################
			
			# plot(hs, col=hsCols, legend=FALSE, axes=FALSE, mar=mar)
			# plot(nam1LgCrop, lwd=3, add=TRUE)
			# plot(usa2LgCrop, lwd=1, add=TRUE)
			# plot(elev, legend = FALSE, axes = FALSE, add = TRUE, col = elevCols, breaks = elevBreaks)

			# plot(oldEvid, pch=22, bg=alpha('darkgoldenrod3', 0.5), cex=1.2, add=TRUE)

			# legend('bottomleft', inset=legendInset, legend=c('Previously occupied'), pch=c(22), pt.bg=c('darkgoldenrod3'), bg='white', cex=1.05)
			
		# ### no evidence
		# ###############
			
			# plot(hs, col=hsCols, legend=FALSE, axes=FALSE, mar=mar)
			# plot(nam1LgCrop, lwd=3, add=TRUE)
			# plot(usa2LgCrop, lwd=1, add=TRUE)
			# plot(elev, legend = FALSE, axes = FALSE, add = TRUE, col = elevCols, breaks = elevBreaks)

			# plot(noEvid, pch=25, bg=alpha('firebrick2', 0.5), cex=1.2, add=TRUE)

			# legend('bottomleft', inset=legendInset, legend=c('No evidence'), pch=c(25), pt.bg=c('firebrick2'), bg='white', cex=1.05)
			
		# ### all together
		# ################
			
			# plot(hs, col=hsCols, legend=FALSE, axes=FALSE, mar=mar)
			# plot(nam1LgCrop, lwd=3, add=TRUE)
			# plot(usa2LgCrop, lwd=1, add=TRUE)
			# plot(elev, legend = FALSE, axes = FALSE, add = TRUE, col = elevCols, breaks = elevBreaks)

			# noEvid <- pikaVectProj[pika$latestOccStatus == '0 never']
			# oldEvid <- pikaVectProj[pika$latestOccStatus == '1 old']
			# occs <- pikaVectProj[pika$latestOccStatus == '2 occupied']
			
			# plot(noEvid, pch=25, bg=alpha('firebrick2', 0.5), cex=1.2, add=TRUE)
			# plot(oldEvid, pch=22, bg=alpha('darkgoldenrod3', 0.5), cex=1.2, add=TRUE)
			# plot(occs, pch=21, bg=alpha('chartreuse', 0.5), cex=1.2, add=TRUE)

			# legend('bottomleft', inset=legendInset, legend=c('Currently occupied', 'Previously occupied', 'No evidence'), pch=c(21, 22, 25), pt.bg=c('chartreuse', 'darkgoldenrod3', 'firebrick2'), bg='white', cex=1.05)
			
		# ### color ramp
		# ##############
		# legendBreaks(
			# x = 'bottomright',
			# inset = legendInset,
			# width = 0.2,
			# height = 0.38,
			# labels = roundTo(elevBreaks, 10),
			# labAdjX = 0.7,
			# labAdjY = c(0.05, 0.25, 0.5, 0.75, 0.98),
			# col = elevCols,
			# colBorder = 'black',
			# title = 'Elevation\n(m)',
			# titleAdj = c(0.5, 0.88),
			# adjX = c(0.1, 0.4),
			# adjY = c(0.05, 0.72)
		# )
		
		# ### inset
		# #########
			
			# # par(fig = c(0.455, 0.715, 0.237, 0.537), bg='white', new=TRUE)
			# par(fig = c(0.455, 0.715, 0.237 - 0.01, 0.537 - 0.01), bg='white', new=TRUE)
			
			# ext <- ext(nam1XL)
			# corners <- as.vector(ext)
			# ext <- raster::extent(corners)
			# ext <- as(ext, 'SpatialPolygons')
			# projection(ext) <- getCRS('naAlbers')
			# ext <- vect(ext)

			# nam1XXL <- crop(nam1XXL, ext)
			# plot(ext, border='gray', col='gray', axes=FALSE)
			# plot(nam1XXL, col='white', border='gray40', ann=FALSE, add=TRUE)
			
			# # highlight study region in inset
			# pikaBuffProjSmExt <- ext(pikaBuffProjSm)
			# corners <- as.vector(pikaBuffProjSmExt)
			# focus <- raster::extent(corners)
			# focus <- as(focus, 'SpatialPolygons')
			# projection(focus) <- getCRS('naAlbers')
			# focus <- vect(focus)
			# plot(focus, lwd=1.9, border='black', add=TRUE)
			
			# plot(ext, add=TRUE)
			
	# dev.off()

say('################################################')
say('### map of sampling sites & GBIF occurrences ###')
say('################################################')

	# survey sites
	load('./Data/02 New Mexico Pika - Environmental Values Extracted and Calculated.rda')
	pikaVectUnproj <- vect(as.matrix(pika[ , ll]), 'points', crs=getCRS('wgs84'))
	
	# plot extent
	pikaBuffUnprojXXL <- terra::buffer(pikaVectUnproj, width=2000000)
	pikaBuffUnprojXL <- terra::buffer(pikaVectUnproj, width=1000000)
	pikaBuffUnprojLg <- terra::buffer(pikaVectUnproj, width=200000)
	pikaBuffUnprojSm <- terra::buffer(pikaVectUnproj, width=40000)
	pikaBuffProjSm <- project(pikaBuffUnprojSm, getCRS('naAlbers'))

	# GADM
	mex1 <- gadm(country='MEX', level=1, path=paste0('C:!Scratch'), version=4.1, resolution=1)
	usa1 <- gadm(country='USA', level=1, path=paste0('C:/!Scratch'), version=4.1, resolution=1)
	usa2 <- gadm(country='USA', level=2, path=paste0('C:/!Scratch'), version=4.1, resolution=1)

	usa1 <- usa1[usa1$NAME_1 != 'Alaska', ]
	usa1 <- usa1[usa1$NAME_1 != 'Hawaii', ]

	usa2 <- usa2[usa2$NAME_1 != 'Alaska', ]
	usa2 <- usa2[usa2$NAME_1 != 'Hawaii', ]

	nam1 <- rbind(usa1, mex1)
	
	nam1XXL <- crop(nam1, ext(pikaBuffUnprojXXL))
	nam1XL <- crop(nam1, ext(pikaBuffUnprojXL))
	nam1Lg <- crop(nam1, ext(pikaBuffUnprojLg))
	usa2Lg <- crop(usa2, ext(pikaBuffUnprojLg))
	
	# elevation
	elev_fine_m <- rast('./Data/elev_fine_m.tif')
	elev_fine_m <- crop(elev_fine_m, pikaBuffUnprojLg)
	elev_fine_m <- crop(elev_fine_m, pikaBuffUnprojLg)
	
	slope <- terrain(elev_fine_m, 'slope', unit='radians')
	aspect <- terrain(elev_fine_m, 'aspect', unit='radians')
	slopeR <- raster(slope)
	aspectR <- raster(aspect)

	hs <- hillShade(slopeR, aspectR, direction=45)
	hs <- rast(hs)
	hs <- project(hs, getCRS('naAlbers'))
	
	# project
	pikaVectProj <- project(pikaVectUnproj, getCRS('naAlbers'))
	nam1XXL <- project(nam1XXL, getCRS('naAlbers'))
	nam1XL <- project(nam1XL, getCRS('naAlbers'))
	nam1Lg <- project(nam1Lg, getCRS('naAlbers'))
	usa2Lg <- project(usa2Lg, getCRS('naAlbers'))
	
	# split records
	noEvid <- pikaVectProj[pika$latestOccStatus == '0 never']
	oldEvid <- pikaVectProj[pika$latestOccStatus == '1 old']
	occs <- pikaVectProj[pika$latestOccStatus == '2 occupied']

	# plot extent
	pikaBuffProjSmSp <- as(pikaBuffProjSm, 'Spatial')
	ext <- extent(pikaBuffProjSmSp)
	ext <- as(ext, 'SpatialPolygons')
	projection(ext) <- getCRS('naAlbers')

	usa2LgCrop <- crop(usa2Lg, ext)
	nam1LgCrop <- crop(nam1Lg, ext)
	hs <- crop(hs, ext)
	hs <- hs - global(hs, 'min', na.rm=TRUE)$min
	hs <- hs / global(hs, 'max', na.rm=TRUE)$max
	hs <- crop(hs, pikaBuffProjSmSp)
	hsR <- rast(hs)
	
	elev <- project(elev_fine_m, getCRS('naAlbers'))
	elev <- crop(elev, pikaBuffProjSmSp)

	# GBIF
	gbif <- read.csv('./Data/GBIF 2022-12-19/ochotona_princeps.csv')
	gbif <- gbif[gbif$species == 'Ochotona princeps', ]
	gbif <- gbif[gbif$year <= 2015, ]
	
	ext <- as.vector(ext(project(usa2LgCrop, getCRS('wgs84'))))
	
	gbif <- vect(gbif, geom=c('decimalLongitude', 'decimalLatitude'), crs=getCRS('wgs84'), keepgeom = TRUE)
	crs(gbif) <- getCRS('wgs84')
	
	gbif <- project(gbif, getCRS('naAlbers'))
	
	ins <- extract(nam1LgCrop, gbif)
	ins$CC_1 <- NULL
	gbif <- gbif[ins$COUNTRY == 'United States', ]

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
	png('./Figures & Tables/Study Region with Sampling Sites with GBIF.png', width=2300, height=2400, res=300)

		par(mar=rep(0.1, 4), mai=rep(0.1, 4), oma=rep(0, 4), mgp=c(0, 0, 0))
			
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
			plot(gbif, pch=21, cex=1, bg='blue', add=TRUE)

			legend('bottomleft', inset=legendInset, legend=c('Currently occupied', 'Previously occupied', 'No evidence', 'GBIF'), pch=c(21, 22, 25, 21), pt.bg=c('chartreuse', 'darkgoldenrod3', 'firebrick2', 'blue'), bg='white', cex=1.05)
			
			
	dev.off()

say('DONE!!!', level=1, deco='%')
