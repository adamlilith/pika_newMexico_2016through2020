### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('C:/Kaji/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/02 New Mexico Pika Occupancy & Abundance Analysis - Maps.r')
###
### CONTENTS ###
### setup ###
### fetch elevation raster ###
### map of sampling sites ###
### map of sampling sites & GBIF occurrences ###
### map of sampling sites with color gradients for elevational bands ###
### map of sampling sites with solid colors for elevational bands ###

#############
### setup ###
#############

	rm(list = ls())

	source('C:/Kaji/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r')

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

# say('################################################')
# say('### map of sampling sites & GBIF occurrences ###')
# say('################################################')

# 	# survey sites
# 	load('./Data/02 New Mexico Pika - Environmental Values Extracted and Calculated.rda')
# 	pikaVectUnproj <- vect(as.matrix(pika[ , ll]), 'points', crs=getCRS('wgs84'))
	
# 	# plot extent
# 	pikaBuffUnprojXXL <- terra::buffer(pikaVectUnproj, width=2000000)
# 	pikaBuffUnprojXL <- terra::buffer(pikaVectUnproj, width=1000000)
# 	pikaBuffUnprojLg <- terra::buffer(pikaVectUnproj, width=200000)
# 	pikaBuffUnprojSm <- terra::buffer(pikaVectUnproj, width=40000)
# 	pikaBuffProjSm <- project(pikaBuffUnprojSm, getCRS('naAlbers'))

# 	# GADM
# 	mex1 <- gadm(country='MEX', level=1, path=paste0('C:!Scratch'), version=4.1, resolution=1)
# 	usa1 <- gadm(country='USA', level=1, path=paste0('C:/!Scratch'), version=4.1, resolution=1)
# 	usa2 <- gadm(country='USA', level=2, path=paste0('C:/!Scratch'), version=4.1, resolution=1)

# 	usa1 <- usa1[usa1$NAME_1 != 'Alaska', ]
# 	usa1 <- usa1[usa1$NAME_1 != 'Hawaii', ]

# 	usa2 <- usa2[usa2$NAME_1 != 'Alaska', ]
# 	usa2 <- usa2[usa2$NAME_1 != 'Hawaii', ]

# 	nam1 <- rbind(usa1, mex1)
	
# 	nam1XXL <- crop(nam1, ext(pikaBuffUnprojXXL))
# 	nam1XL <- crop(nam1, ext(pikaBuffUnprojXL))
# 	nam1Lg <- crop(nam1, ext(pikaBuffUnprojLg))
# 	usa2Lg <- crop(usa2, ext(pikaBuffUnprojLg))
	
# 	# elevation
# 	elev_fine_m <- rast('./Data/elev_fine_m.tif')
# 	elev_fine_m <- crop(elev_fine_m, pikaBuffUnprojLg)
# 	elev_fine_m <- crop(elev_fine_m, pikaBuffUnprojLg)
	
# 	slope <- terrain(elev_fine_m, 'slope', unit='radians')
# 	aspect <- terrain(elev_fine_m, 'aspect', unit='radians')
# 	slopeR <- raster(slope)
# 	aspectR <- raster(aspect)

# 	hs <- hillShade(slopeR, aspectR, direction=45)
# 	hs <- rast(hs)
# 	hs <- project(hs, getCRS('naAlbers'))
	
# 	# project
# 	pikaVectProj <- project(pikaVectUnproj, getCRS('naAlbers'))
# 	nam1XXL <- project(nam1XXL, getCRS('naAlbers'))
# 	nam1XL <- project(nam1XL, getCRS('naAlbers'))
# 	nam1Lg <- project(nam1Lg, getCRS('naAlbers'))
# 	usa2Lg <- project(usa2Lg, getCRS('naAlbers'))
	
# 	# split records
# 	noEvid <- pikaVectProj[pika$latestOccStatus == '0 never']
# 	oldEvid <- pikaVectProj[pika$latestOccStatus == '1 old']
# 	occs <- pikaVectProj[pika$latestOccStatus == '2 occupied']

# 	# plot extent
# 	pikaBuffProjSmSp <- as(pikaBuffProjSm, 'Spatial')
# 	ext <- extent(pikaBuffProjSmSp)
# 	ext <- as(ext, 'SpatialPolygons')
# 	projection(ext) <- getCRS('naAlbers')

# 	usa2LgCrop <- crop(usa2Lg, ext)
# 	nam1LgCrop <- crop(nam1Lg, ext)
# 	hs <- crop(hs, ext)
# 	hs <- hs - global(hs, 'min', na.rm=TRUE)$min
# 	hs <- hs / global(hs, 'max', na.rm=TRUE)$max
# 	hs <- crop(hs, pikaBuffProjSmSp)
# 	hsR <- rast(hs)
	
# 	elev <- project(elev_fine_m, getCRS('naAlbers'))
# 	elev <- crop(elev, pikaBuffProjSmSp)

# 	# GBIF
# 	gbif <- read.csv('./Data/GBIF 2022-12-19/ochotona_princeps.csv')
# 	gbif <- gbif[gbif$species == 'Ochotona princeps', ]
# 	gbif <- gbif[gbif$year <= 2015, ]
	
# 	ext <- as.vector(ext(project(usa2LgCrop, getCRS('wgs84'))))
	
# 	gbif <- vect(gbif, geom=c('decimalLongitude', 'decimalLatitude'), crs=getCRS('wgs84'), keepgeom = TRUE)
# 	crs(gbif) <- getCRS('wgs84')
	
# 	gbif <- project(gbif, getCRS('naAlbers'))
	
# 	ins <- extract(nam1LgCrop, gbif)
# 	ins$CC_1 <- NULL
# 	gbif <- gbif[ins$COUNTRY == 'United States', ]

# 	### colors
	
# 		# colors for elevation
# 		load('./Data/04 New Mexico Pika - Added Distance to Closest Patches.rda')
		
# 		minElev <- globalx(elev, 'min')
# 		lowestNeverOcc <- min(pika$elevation_m[pika$latestOccStatus == '0 never'])
# 		# medianNeverOcc <- median(pika$elevation_m[pika$latestOccStatus == '0 never'])
# 		medianPastOcc <- median(pika$elevation_m[pika$latestOccStatus == '1 old'])
# 		medianOcc <- median(pika$elevation_m[pika$latestOccStatus == '2 occupied'])
# 		maxElev <- globalx(elev, 'max')

# 		elevBreaks <- c(minElev, lowestNeverOcc, medianPastOcc, medianOcc, maxElev)
# 		names(elevBreaks) <- c('minElev', 'medianNeverOcc', 'medianPastOcc', 'medianOcc', 'maxElev')

# 		elevCols <- c('gray80', 'firebrick3', 'darkgoldenrod3', 'chartreuse')
# 		elevCols <- alpha(elevCols, 0.5)

# 		# hillshade colors
# 		hsCols <- colorRampPalette(c('gray0', 'gray100'))
# 		hsCols <- hsCols(20)
		
# 	### placement
	
# 		legendInset <- c(0.019, 0.02)
# 		mar <- c(0, 0.5, 1, 1)
		
# 	### plot!
# 	png('./Figures & Tables/Study Region with Sampling Sites with GBIF.png', width=2300, height=2400, res=300)

# 		par(mar=rep(0.1, 4), mai=rep(0.1, 4), oma=rep(0, 4), mgp=c(0, 0, 0))
			
# 		### all together
# 		################
			
# 			plot(hs, col=hsCols, legend=FALSE, axes=FALSE, mar=mar)
# 			plot(nam1LgCrop, lwd=3, add=TRUE)
# 			plot(usa2LgCrop, lwd=1, add=TRUE)
# 			plot(elev, legend = FALSE, axes = FALSE, add = TRUE, col = elevCols, breaks = elevBreaks)

# 			noEvid <- pikaVectProj[pika$latestOccStatus == '0 never']
# 			oldEvid <- pikaVectProj[pika$latestOccStatus == '1 old']
# 			occs <- pikaVectProj[pika$latestOccStatus == '2 occupied']
			
# 			plot(noEvid, pch=25, bg=alpha('firebrick2', 0.5), cex=1.2, add=TRUE)
# 			plot(oldEvid, pch=22, bg=alpha('darkgoldenrod3', 0.5), cex=1.2, add=TRUE)
# 			plot(occs, pch=21, bg=alpha('chartreuse', 0.5), cex=1.2, add=TRUE)
# 			plot(gbif, pch=21, cex=1, bg='blue', add=TRUE)

# 			legend('bottomleft', inset=legendInset, legend=c('Currently occupied', 'Previously occupied', 'No evidence', 'GBIF'), pch=c(21, 22, 25, 21), pt.bg=c('chartreuse', 'darkgoldenrod3', 'firebrick2', 'blue'), bg='white', cex=1.05)
			
			
# 	dev.off()

# say('########################################################################')
# say('### map of sampling sites with color gradients for elevational bands ###')
# say('########################################################################')

# 	# user-defined

# 		# 2935 m is the 25th quantile of currently occupied sites
# 		# 2565 m is the lowest elevation with current pika presence
# 		# 2303 m is the lowest elevation with past pika presence
# 		# 2280 m is the lowest elevation sampled

# 		min_high_elev_m <- 2935
# 		min_mid_elev_m <- 2565
# 		min_min_elev_m <- 2303

# 	# North America
# 	nam1 <- gadm(country = c('CAN', 'USA', 'MEX'), level = 1, path = paste0('C:/!Scratch/gadm'), version = 4.1, resolution = 1)
# 	nam1 <- nam1[nam1$NAME_1 %notin% c('Alaska', 'Hawaii')]
# 	nam1 <- nam1[nam1$NAME_1 %in% c('California', 'Oregon', 'Washington', 'British Columbia', 'Alberta', 'Saskatchewan', 'Manitoba', 'Yukon', 'Northwest Territories', 'Idaho', 'Montana', 'North Dakota', 'South Dakota', 'Colorado', 'Nebraska', 'Wyoming', 'Utah', 'Nevada', 'Arizona', 'New Mexico', 'Baja California', 'Sonora', 'Oklahoma', 'Texas', 'Kansas', 'Nebraska')]

# 	# IUCN range
# 	iucn <- vect('./Data/IUCN Range Map 2025-03-12/pika_range_map.gpkg')

# 	# survey sites
# 	load('./Data/02 New Mexico Pika - Environmental Values Extracted and Calculated.rda')
# 	pika <- vect(pika, geom = ll, crs = getCRS('NAD83'))

# 	# plot extent
# 	extent_nm <- buffer(pika, width = 20000)
# 	extent_nam <- buffer(iucn, width = 600000)

# 	extent_nm <- ext(extent_nm)
# 	extent_nam <- ext(extent_nam)

# 	extent_nm <- as.polygons(extent_nm, crs = getCRS('WGS84'))
# 	extent_nam <- as.polygons(extent_nam, crs = getCRS('WGS84'))

# 	nam1_nm <- crop(nam1, extent_nm)

# 	# cities
# 	cities <- data.frame(
# 		name = c('Santa Fe', 'Los Alamos'),
# 		x = c(-105.964444, -106.263889),
# 		y = c(35.667222, 35.891111)
# 	)
# 	cities <- vect(cities, geom = c('x', 'y'), crs = getCRS('NAD83'))
	
# 	# major rivers and lakes
# 	rivers <- vect('C:/Kaji/Research Data/Rivers and Lakes - North America USGS/hydrography_l_rivers_v2.shp') # Adjust path to your rivers data
# 	rivers <- makeValid(rivers)
# 	extent_nm_proj <- project(extent_nm, rivers)
# 	rivers <- crop(rivers, extent_nm_proj)
	
# 	# elevation
# 	elev_fine_m <- rast('./Data/elev_fine_m.tif')
# 	elev_fine_m <- crop(elev_fine_m, extent_nm)
	
# # elev_fine_m <- aggregate(elev_fine_m, 16, mean)

# 	# hillshade
# 	slope <- terrain(8 * elev_fine_m, 'slope', unit = 'radians')
# 	aspect <- terrain(8 * elev_fine_m, 'aspect', unit = 'radians')
# 	hs <- shade(slope, aspect, angle = 45, direction = 315)

# 	# elevation above given threshold
# 	elev_high_m <- elev_mid_m <- elev_low_m <- elev_fine_m
# 	elev_high_m[elev_high_m < min_high_elev_m] <- NA
# 	elev_mid_m[elev_mid_m < min_mid_elev_m | elev_mid_m >= min_high_elev_m] <- NA
# 	elev_low_m[elev_low_m < min_min_elev_m | elev_low_m >= min_mid_elev_m] <- NA

# 	max_high_elev_m <- globalx(elev_high_m, 'max')

# 	# project
# 	hs <- project(hs, getCRS('North America Lambert'))
# 	pika <- project(pika, getCRS('North America Lambert'))
# 	nam1 <- project(nam1, getCRS('North America Lambert'))
# 	cities <- project(cities, getCRS('North America Lambert'))
# 	rivers <- project(rivers, getCRS('North America Lambert'))
# 	iucn <- project(iucn, getCRS('North America Lambert'))
# 	extent_nm_proj <- project(extent_nm_proj, getCRS('North America Lambert'))
# 	extent_nam_proj <- project(extent_nam, getCRS('North America Lambert'))

# 	nam1_iucn <- iucn * nam1 

# 	extent_nam_proj_vect <- as.vector(ext(extent_nam_proj))

# 	iucn_nm <- crop(iucn, extent_nm_proj)

# 	extent_nam_proj <- ext(extent_nam_proj)
# 	extent_nam_proj <- as.vector(extent_nam_proj)

# 	extent_nm_proj_vect <- as.vector(ext(extent_nm_proj))
# 	say('Extent along x-axis is ', (extent_nm_proj_vect[2] - extent_nm_proj_vect[1]) / 1000, ' km')

# 	# split records
# 	present <- pika[pika$latestOccStatus %in% c('2 occupied')]
# 	absent <- pika[pika$latestOccStatus %in% c('0 never', '1 old')]

# 	# hillshade colors
# 	hs_cols <- colorRampPalette(c('gray30', 'gray100'))(20)
# 	high_cols <- colorRampPalette(c('forestgreen', 'green2'))(5)
# 	mid_cols <- colorRampPalette(c('goldenrod3', 'gold'))(5)
# 	low_cols <- colorRampPalette(c('indianred4', 'indianred1'))(5)

# 	high_cols <- alpha(high_cols, 0.7)
# 	mid_cols <- alpha(mid_cols, 0.7)
# 	low_cols <- alpha(low_cols, 0.7)

# 	### range map	
# 	extent_nam_proj[1] <- extent_nam_proj[1] + 1300000
# 	extent_nam_proj[2] <- extent_nam_proj[2] - 600000
# 	extent_nam_proj[3] <- extent_nam_proj[3] + 600000
# 	extent_nam_proj[4] <- extent_nam_proj[4] - 900000
	
# 	range_map <- ggplot() +
# 		layer_spatial(nam1, fill = 'gray85') +
# 		layer_spatial(iucn, fill = 'gray40') +
# 		layer_spatial(extent_nm_proj, color = 'black', fill = NA, linewidth = 1) +
# 		layer_spatial(nam1_iucn, color = 'gray60', fill = NA) +
# 		xlim(extent_nam_proj[1], extent_nam_proj[2]) + ylim(extent_nam_proj[3], extent_nam_proj[4]) +
# 		theme(
# 			axis.text = element_text(size = 16)
# 		)

# 	### study region map
# 	sr_map <- ggplot() +
# 		layer_spatial(hs, aes(fill = stat(band1))) +
# 		scale_fill_gradientn(colors = hs_cols, guide = 'none', na.value = 'transparent') +
# 		new_scale_fill() +
# 		layer_spatial(elev_low_m, aes(fill = stat(band1))) +
# 		scale_fill_gradientn(colors = low_cols, name = 'Low\nElevation (m)', na.value = 'transparent', breaks = c(min_min_elev_m, min_mid_elev_m), guide = guide_colorbar(order = 3, label = FALSE)) +
# 		new_scale_fill() +
# 		layer_spatial(elev_mid_m, aes(fill = stat(band1))) +
# 		scale_fill_gradientn(colors = mid_cols, name = 'Middle\nElevation (m)', na.value = 'transparent', breaks = c(min_mid_elev_m, min_high_elev_m), guide = guide_colorbar(order = 2, label = FALSE)) +
# 		new_scale_fill() +
# 		layer_spatial(elev_high_m, aes(fill = stat(band1))) +
# 		scale_fill_gradientn(colors = high_cols, name = 'High\nElevation (m)', na.value = 'transparent', breaks = c(min_high_elev_m, max_high_elev_m), guide = guide_colorbar(order = 1, label = FALSE)) +
# 		layer_spatial(rivers, color = 'blue', size = 0.5) +
# 		layer_spatial(iucn_nm, fill = NA, color = 'black', linewidth = 1.2, linetype = 'dashed') +
# 		layer_spatial(absent, pch = 2, size = 4.1, alpha = 1, color = 'red') +
# 		layer_spatial(present, pch = 1, size = 4.2, alpha = 1, color = 'black') +
# 		layer_spatial(cities, pch = 19, size = 5) +
# 		layer_spatial(nam1_nm) +
# 		geom_sf_text(data = st_as_sf(cities), aes(label = name), 
# 			nudge_x = c(-10000, 10000), nudge_y = c(-5000, -5000),
# 			size = 4.5, fontface = 'bold'
# 		) +
# 		xlim(extent_nm_proj_vect[1], extent_nm_proj_vect[2]) +
# 		ylim(extent_nm_proj_vect[3], extent_nm_proj_vect[4]) +
# 		coord_sf(expand = FALSE) +
# 		theme_void() +
# 		theme(
# 			legend.position = c(-0.08, 0.3),
# 			legend.key.height = unit(0.5, 'cm'),
# 			legend.title = element_text(size = 16),
# 			legend.text = element_text(size = 16),
# 			plot.margin = margin(t = 1, r = 1, b = 1, l = 100, unit = 'pt')
# 		)
	
# 	ggsave(sr_map, filename = './Figures & Tables/Study Region with Sampling Sites V2 Study Region.png', width = 10, height = 9, dpi = 600, bg = 'white')
# 	ggsave(range_map, filename = './Figures & Tables/Study Region with Sampling Sites V2 Range Map.png', width = 6, height = 8, dpi = 600, bg = 'white')

say('#####################################################################')
say('### map of sampling sites with solid colors for elevational bands ###')
say('#####################################################################')

	# user-defined

		# 2935 m is the 25th quantile of currently occupied sites
		# 2565 m is the lowest elevation with current pika presence
		# 2303 m is the lowest elevation with past pika presence
		# 2280 m is the lowest elevation sampled

		min_high_elev_m <- 2935
		min_mid_elev_m <- 2565

	# North America
	nam1 <- gadm(country = c('CAN', 'USA', 'MEX'), level = 1, path = paste0('C:/!Scratch/gadm'), version = 4.1, resolution = 1)
	nam1 <- nam1[nam1$NAME_1 %notin% c('Alaska', 'Hawaii')]
	nam1 <- nam1[nam1$NAME_1 %in% c('California', 'Oregon', 'Washington', 'British Columbia', 'Alberta', 'Saskatchewan', 'Manitoba', 'Yukon', 'Northwest Territories', 'Idaho', 'Montana', 'North Dakota', 'South Dakota', 'Colorado', 'Nebraska', 'Wyoming', 'Utah', 'Nevada', 'Arizona', 'New Mexico', 'Baja California', 'Sonora', 'Oklahoma', 'Texas', 'Kansas', 'Nebraska')]

	# IUCN range
	iucn <- vect('./Data/IUCN Range Map 2025-03-12/pika_range_map.gpkg')

	# survey sites
	load('./Data/02 New Mexico Pika - Environmental Values Extracted and Calculated.rda')
	pika <- vect(pika, geom = ll, crs = getCRS('NAD83'))

	# plot extent
	extent_nm <- buffer(pika, width = 20000)
	extent_nam <- buffer(iucn, width = 600000)

	extent_nm <- ext(extent_nm)
	extent_nam <- ext(extent_nam)

	extent_nm <- as.polygons(extent_nm, crs = getCRS('WGS84'))
	extent_nam <- as.polygons(extent_nam, crs = getCRS('WGS84'))

	nam1_nm <- crop(nam1, extent_nm)

	# cities
	cities <- data.frame(
		name = c('Santa Fe', 'Los Alamos'),
		x = c(-105.964444, -106.263889),
		y = c(35.667222, 35.891111)
	)
	cities <- vect(cities, geom = c('x', 'y'), crs = getCRS('NAD83'))
	
	# major rivers and lakes
	rivers <- vect('C:/Kaji/Research Data/Rivers and Lakes - North America USGS/hydrography_l_rivers_v2.shp') # Adjust path to your rivers data
	rivers <- makeValid(rivers)
	extent_nm_proj <- project(extent_nm, rivers)
	rivers <- crop(rivers, extent_nm_proj)
	
	# elevation
	elev_fine_m <- rast('./Data/elev_fine_m.tif')
	elev_fine_m <- crop(elev_fine_m, extent_nm)
	
# elev_fine_m <- aggregate(elev_fine_m, 16, mean)

	# hillshade
	slope <- terrain(20 * elev_fine_m, 'slope', unit = 'radians')
	aspect <- terrain(20 * elev_fine_m, 'aspect', unit = 'radians')
	hs <- shade(slope, aspect, angle = 45, direction = 315)
	hs <- hs - globalx(hs, 'mean')
	
	# elevation above given threshold
	elev_high_m <- elev_mid_m <- elev_fine_m
	elev_high_m[elev_high_m < min_high_elev_m] <- NA
	elev_mid_m[elev_mid_m < min_mid_elev_m | elev_mid_m >= min_high_elev_m] <- NA

	elev_high <- elev_high_m * 0 + 1
	elev_mid <- elev_mid_m * 0 + 1

	# project
	hs <- project(hs, getCRS('North America Lambert'))
	pika <- project(pika, getCRS('North America Lambert'))
	nam1 <- project(nam1, getCRS('North America Lambert'))
	cities <- project(cities, getCRS('North America Lambert'))
	rivers <- project(rivers, getCRS('North America Lambert'))
	iucn <- project(iucn, getCRS('North America Lambert'))
	extent_nm_proj <- project(extent_nm_proj, getCRS('North America Lambert'))
	extent_nam_proj <- project(extent_nam, getCRS('North America Lambert'))
	elev_high <- project(elev_high, getCRS('North America Lambert'), method = 'near')
	elev_mid <- project(elev_mid, getCRS('North America Lambert'), method = 'near')

	nam1_iucn <- iucn * nam1 

	extent_nam_proj_vect <- as.vector(ext(extent_nam_proj))

	iucn_nm <- crop(iucn, extent_nm_proj)

	extent_nam_proj <- ext(extent_nam_proj)
	extent_nam_proj <- as.vector(extent_nam_proj)

	extent_nm_proj_vect <- as.vector(ext(extent_nm_proj))
	say('Extent along x-axis is ', (extent_nm_proj_vect[2] - extent_nm_proj_vect[1]) / 1000, ' km')

	# split records
	present <- pika[pika$latestOccStatus %in% c('2 occupied')]
	absent <- pika[pika$latestOccStatus %in% c('0 never', '1 old')]

	# hillshade colors
	hs_cols <- colorRampPalette(c('gray10', 'gray100'))(20)
	high_cols <- 'forestgreen'
	mid_cols <- 'springgreen3'

	### range map	
	extent_nam_proj[1] <- extent_nam_proj[1] + 1300000
	extent_nam_proj[2] <- extent_nam_proj[2] - 600000
	extent_nam_proj[3] <- extent_nam_proj[3] + 600000
	extent_nam_proj[4] <- extent_nam_proj[4] - 900000
	
	range_map <- ggplot() +
		layer_spatial(nam1, fill = 'gray85') +
		layer_spatial(iucn, color = 'black', fill = '#6F4E37', linetype = 'solid') +
		layer_spatial(extent_nm_proj, color = 'black', fill = NA, linewidth = 1) +
		layer_spatial(nam1_iucn, color = 'gray60', fill = NA) +
		xlim(extent_nam_proj[1], extent_nam_proj[2]) + ylim(extent_nam_proj[3], extent_nam_proj[4]) +
		theme(
			axis.text = element_text(size = 16)
		)

	### study region map
	sr_map <- ggplot() +
		layer_spatial(hs, aes(fill = stat(band1))) +
		scale_fill_gradientn(colors = hs_cols, guide = 'none', na.value = 'transparent') +
		new_scale_fill() +
		layer_spatial(elev_mid, aes(fill = stat(band1)), alpha = 0.5) +
		scale_fill_gradientn(
			colours = c(mid_cols[1], mid_cols[1]),
			breaks = 1,
			labels = 'Mid-elevation\n(2303-2829 m)',
			name = NULL,
			na.value = 'transparent'
		) +
		new_scale_fill() +
		layer_spatial(elev_high, aes(fill = stat(band1)), alpha = 0.5) +
		scale_fill_gradientn(
			colours = c(high_cols[1], high_cols[1]),
			breaks = 1,
			labels = 'High-elevation\n(>2829 m)',
			name = NULL,
			na.value = 'transparent'
		) +
		layer_spatial(rivers, color = 'blue', size = 0.8) +
		layer_spatial(iucn_nm, fill = NA, color = 'black', linewidth = 1.3, linetype = 'dashed') +
		layer_spatial(absent, pch = 2, size = 3.8, alpha = 1, color = 'red') +
		layer_spatial(present, pch = 1, size = 4.1, alpha = 1, color = 'black') +
		layer_spatial(cities, pch = 19, size = 5) +
		layer_spatial(nam1_nm) +
		geom_sf_text(data = st_as_sf(cities), aes(label = name), 
			nudge_x = c(-10000, 10000), nudge_y = c(-5000, -5000),
			size = 4.5, fontface = 'bold'
		) +
		coord_sf(xlim = c(extent_nm_proj_vect[1], extent_nm_proj_vect[2]), ylim = c(extent_nm_proj_vect[3], extent_nm_proj_vect[4]), expand = FALSE) +
		theme_void() +
		theme(
			legend.position = 'none'
			# legend.key.height = unit(0.7, 'cm'),
			# legend.title = element_text(size = 16),
			# legend.text = element_text(size = 16),
			# plot.margin = margin(t = 1, r = 1, b = 5, l = 5, unit = 'pt')
		)
	
	ggsave(sr_map, filename = './Figures & Tables/Study Region with Sampling Sites V3 Study Region.png', width = 10, height = 9, dpi = 600, bg = 'white')
	ggsave(range_map, filename = './Figures & Tables/Study Region with Sampling Sites V3 Range Map.png', width = 6, height = 8, dpi = 600, bg = 'white')


say('DONE!!!', level=1, deco='%')
