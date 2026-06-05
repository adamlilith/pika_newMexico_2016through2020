# source('C:/Kaji/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/TEMP.r')

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
