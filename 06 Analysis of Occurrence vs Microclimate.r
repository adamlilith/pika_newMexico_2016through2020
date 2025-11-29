### NEW MEXICO PIKA ANALYSIS
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2021-04
###
### source('C:/Kaji/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/06 Analysis of Occurrence vs Microclimate.r')
###
### CONTENTS ###
### setup ###
### custom functions ###
### calculate correlations between microclimate variables for single-year data (Marie Westover) ###
### calculate correlations between microclimate variables for multi-year data (Erik Beever) ###

### exploratory data analysis for single-year data (Marie Westover) ###
### analyze patterns of occupancy for single-year data (Maria Westover) ###

### exploratory data analysis for multi-year data (Erik Beever) ###
### analyze patterns of occupancy for multi-year data (Erik Beever) ###

### joint exploratory data analysis for single- and multi year data ###

#############
### setup ###
#############

	rm(list=ls())

	drive <- 'C:/Kaji/'

	source(paste0(drive, '/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/00 New Mexico Pika Occupancy & Abundance Analysis - Shared Functions & Constants.r'))

	dirCreate('./Figures & Tables/Microclimate Analysis')

########################
### custom functions ###
########################

	### Make nice versions of single-year variables
	# vars		Vector of variable names
	# units		If TRUE, include units
	make_nice_microclimate_vars_single_year <- function(vars, units = TRUE) {
	
		vars[vars == 'chroniccoldavgNov-Mar'] <- paste0('Chronic cold', ifelse(units, ' (°C)', ''))
		vars[vars == 'chroniccoldavgNovMar'] <- paste0('Chronic cold', ifelse(units, ' (°C)', ''))
		vars[vars == 'ChronicHeatavgJunSept'] <- paste0('Chronic heat', ifelse(units, ' (°C)', ''))
		vars[vars == 'Hot22Days'] <- paste0('Hot days', ifelse(units, ' (log₁₀ + 1)', ''))
		vars[vars == 'InferredSnowDays'] <- paste0('Snow days', ifelse(units, ' (d)', ''))
		vars[vars == 'MeanJulyTemp'] <- paste0('July temperature', ifelse(units, ' (°C)', ''))
		vars[vars == 'Min.temp.x'] <- paste0('Minimum temperature', ifelse(units, ' (°C)', ''))
		vars[vars == 'meanDistToClosest4Patches_m'] <- paste0('Isolation', ifelse(units, ' (log₁₀, m)', ''))
		vars[vars == 'numHomeRanges'] <- paste0('Number of home ranges', ifelse(units, ' (log₁₀ + 1)', ''))
		vars[vars == 'Min.temp.x:Hot22Days'] <- paste0('Min. temp. × Hot days', ifelse(units, ' (°C d)', ''))
		vars[vars == 'InferredSnowDays:Min.temp.x'] <- paste0('Min. temp × Snow days', ifelse(units, ' (°C d)', ''))
		vars
	
	}

	make_nice_microclimate_vars_multi_year <- function(vars, units = TRUE) {

		vars[vars == 'subLethalCold'] <- paste0('Sub-lethal cold', ifelse(units, ' (°C)', ''))
		vars[vars == 'meltRefreeze'] <- paste0('Melt/re-freeze', ifelse(units, ' (log₁₀ + 1, d)', ''))
		vars[vars == 'chronicCold'] <- paste0('Chronic cold', ifelse(units, ' (°C)', ''))
		vars[vars == 'coldNoSnow'] <- paste0('Cold, no snow', ifelse(units, ' (log₁₀ + 1, d)', ''))
		vars[vars == 'acuteCold'] <- paste0('Acute cold', ifelse(units, ' (log₁₀ + 1, d)', ''))
		vars[vars == 'acuteHeat'] <- paste0('Acute heat', ifelse(units, ' (log₁₀ + 1, d)', ''))
		vars[vars == 'subLethalHeat'] <- paste0('Sub-lethal heat', ifelse(units, ' (log₁₀ + 1, d)', ''))
		vars[vars == 'chronicHeat'] <- paste0('Chronic heat', ifelse(units, ' (°C)', ''))
		vars[vars == 'summerRespite'] <- paste0('Summer respite', ifelse(units, ' (°C)', ''))
		vars[vars == 'meanDistToClosest4Patches'] <- paste0('Isolation', ifelse(units, ' (log₁₀, m)', ''))
		vars[vars == 'numHomeRanges'] <- paste0('Number of home ranges', ifelse(units, ' (log₁₀ + 1)', ''))
		vars[vars == 'subLethalHeat:subLethalCold'] <- paste0('Sub-lethal cold × Sub-lethal heat', ifelse(units, ' (°C²)', ''))
		vars[vars == 'acuteCold:acuteHeat'] <- paste0('Acute cold × Acute heat', ifelse(units, ' (d²)', ''))
		vars[vars == 'occy_prev'] <- paste0('Previous-year occurrence', ifelse(units, '', ''))
		vars	
	}

# say('###################################################################################################')
# say('### calculate correlations between microclimate variables for single-year data (Marie Westover) ###')
# say('###################################################################################################')

# 	# This chunk creates heatmaps and dendrograms of microclimate variables for the single-year microclimate data set (i.e., not the time series of microclimate data). The correlations will be used to construct models with variables of minimal correlation.

# 	x <- fread('./Data/Microclimate/Single Year 2025-08-07 Marie Westover/iButton Data MLW matching site names 8.7.25.csv')

# 	vars <- x[ , c('chroniccoldavgNov-Mar', 'ChronicHeatavgJunSept', 'Hot22Days', 'MeanJulyTemp', 'InferredSnowDays', 'Min.temp.x', 'numHomeRanges', 'meanDistToClosest4Patches_m')]

# 	vars[ , log10_numHomeRanges := log10(numHomeRanges)]
# 	vars[ , log10_meanDistToClosest4Patches_m := log10(meanDistToClosest4Patches_m)]

# 	vars[ , numHomeRanges := NULL]
# 	vars[ , meanDistToClosest4Patches_m := NULL]

# 	library(gplots)

# 	cors <- cor(vars, method = 'spearman', use = 'complete.obs')
	
# 	# labels for each cell with rounded correlation values
# 	cell_labels <- matrix(sprintf('%.2f', cors), nrow = nrow(cors), ncol = ncol(cors))
	
# 	png('./Figures & Tables/Microclimate Analysis - Single Year/Marie - Variable Correlation Heatmap.png', width = 1000, height = 1000, res = 120)
# 	heatmap.2(
# 		cors,
# 		main = '',
# 		col = colorRampPalette(c('#4575b4', '#ffffbf', '#d73027'))(100),
# 		margins = c(16, 16),
# 		symm = TRUE,
# 		trace = 'none',
# 		density.info = 'none',
# 		cellnote = cell_labels,
# 		notecol = 'black',
# 		notecex = 1.2,
# 		key = FALSE
# 	)
# 	title('Marie\'s Data - Spearman Rank Correlations', cex.main = 1.2) # Reduce title size
# 	dev.off()

# 	dists <- 1 - abs(cors)
# 	dists <- as.dist(dists)

# 	dendro <- hclust(dists)

# 	png('./Figures & Tables/Microclimate Analysis - Single Year/Marie - Variable Correlation Dendrogram.png', width = 1000, height = 1000, res = 120)
# 	plot(dendro, ylab = '1 - abs(correlation)', xlab = '', main = 'Marie\'s Data - Spearman Rank Correlations')
# 	abline(h = 0.3, col = 'red')
# 	dev.off()

# say('###############################################################################################')
# say('### calculate correlations between microclimate variables for multi-year data (Erik Beever) ###')
# say('###############################################################################################')

# 	# This chunk creates heatmaps and dendrograms of microclimate variables for the multi-year microclimate data set (i.e., the time series of microclimate data). The correlations will be used to construct models with variables of minimal correlation.

# 	# Bandolier sites with microclimate over > 1 year
# 	band <- fread('./Data/Microclimate/Multi-Year 2025-10-23 Erik Beever & Dylan Ryals/band_clim_14.csv')
# 	band_sites <- tolower(band$patch)
# 	band_sites[band_sites %in% c('band101a', 'band101b')] <- 'band101'

# 	# add isolation and home range size of patches
# 	load('./Data/05 New Mexico Pika - Added PRISM Cell Number & Cell-Based Weight.rda')
# 	all_sites <- as.data.table(pika)
# 	all_sites_sites <- tolower(all_sites$polygonName)
# 	all_sites_sites <- gsub(all_sites_sites, pattern = ' ', replacement = '')

# 	matches <- match(band_sites, all_sites_sites)
# 	band[ , meanDistToClosest4Patches := all_sites$meanDistToClosest4Patches[matches]]
# 	band[ , numHomeRanges := all_sites$numHomeRanges[matches]]

# 	vars <- band[ , c('patch', 'subLethalCold', 'meltRefreeze', 'chronicCold', 'coldNoSnow', 'acuteCold', 'acuteHeat', 'subLethalHeat', 'chronicHeat', 'summerRespite', 'meanDistToClosest4Patches', 'numHomeRanges', 'occy', 'occy_prev')]

# 	vars[ , log10_numHomeRanges := log10(numHomeRanges)]
# 	vars[ , log10_meanDistToClosest4Patches := log10(meanDistToClosest4Patches)]

# 	vars[ , numHomeRanges := NULL]
# 	vars[ , meanDistToClosest4Patches := NULL]

# 	# calculate site means
# 	vars <- vars[ , lapply(.SD, mean, na.rm = TRUE), by = 'patch', .SDcols = setdiff(names(vars), 'patch')]
# 	vars[ , patch := NULL]


# 	cors <- cor(vars, method = 'spearman', use = 'complete.obs')
	
# 	# labels for each cell with rounded correlation values
# 	cell_labels <- matrix(sprintf('%.2f', cors), nrow = nrow(cors), ncol = ncol(cors))
	
# 	png('./Figures & Tables/Microclimate Analysis - Multiple Years/Microclimate Analysis - Multiple Years - Variable Correlation Heatmap.png', width = 1000, height = 1000, res = 120)
# 	heatmap.2(
# 		cors,
# 		main = '',
# 		col = colorRampPalette(c('#4575b4', '#ffffbf', '#d73027'))(100),
# 		margins = c(16, 16),
# 		symm = TRUE,
# 		trace = 'none',
# 		density.info = 'none',
# 		cellnote = cell_labels,
# 		notecol = 'black',
# 		notecex = 0.8,
# 		key = FALSE
# 	)
# 	title(paste0('Multi-Year Data - Spearman Rank Correlations ', date()), cex.main = 1.2) # Reduce title size
# 	dev.off()

# 	dists <- 1 - abs(cors)
# 	dists <- as.dist(dists)

# 	dendro <- hclust(dists)

# 	png('./Figures & Tables/Microclimate Analysis - Multiple Years/Microclimate Analysis - Multiple Years - Variable Correlation Dendrogram.png', width = 1000, height = 1000, res = 120)
# 	plot(dendro, ylab = '1 - abs(correlation)', xlab = '', main = paste0('Multi-year Data - Spearman Rank Correlations ', date()))
# 	abline(h = 0.3, col = 'red')
# 	dev.off()

# say('#######################################################################')
# say('### exploratory data analysis for single-year data (Marie Westover) ###')
# say('#######################################################################')

# 	data <- fread('./Data/Microclimate/Single Year 2025-08-07 Marie Westover/iButton Data MLW matching site names 8.7.25.csv')
# 	model_forms <- read_xlsx('./Data/Microclimate/Single Year 2025-08-07 Marie Westover/PikaiButtonModelsMLW.xlsx', sheet = 'Sheet1')

# 	vars <- c(model_forms$term1, model_forms$term2, model_forms$term3)
# 	vars <- vars[!is.na(vars)]
# 	vars <- sort(unique(vars))
# 	vars <- vars[!grepl(vars, pattern = '\\*')]
# 	vars <- c(vars, 'numHomeRanges', 'meanDistToClosest4Patches_m')

# 	names(data)[names(data) == 'chroniccoldavgNov-Mar'] <- 'chroniccoldavgNovMar'
# 	vars[vars == 'chroniccoldavgNov-Mar'] <- 'chroniccoldavgNovMar'

# 	vars_nice <- make_nice_microclimate_vars_single_year(vars, units = FALSE)
# 	vars_nice_units <- make_nice_microclimate_vars_single_year(vars, units = TRUE)

# 	data$binaryStatus <- factor(data$binaryStatus)

# 	hists <- list()
# 	for (i in seq_along(vars)) {
	
# 		var <- vars[i]
# 		x <- data[[var]]
# 		if (var == 'numHomeRanges') x <- log10(x + 1)
# 		if (var == 'meanDistToClosest4Patches_m') x <- log10(x)
# 		data$var <- x

# 		var_nice <- vars_nice[i]
# 		var_nice_units <- vars_nice_units[i]

# 		# means for each group
# 		means <- data[ , .(mean = mean(var, na.rm = TRUE)), by = binaryStatus]

# 		hists[[i]] <- ggplot(data, aes(x = var, fill = binaryStatus)) +
# 			geom_histogram(
# 				bins = 30, 
# 				color = 'black', 
# 				position = 'stack', 
# 				alpha = 0.5
# 			) +
# 			geom_vline(
# 				data = means,
# 				aes(xintercept = mean, color = binaryStatus),
# 				lwd = 1.1,
# 				show.legend = FALSE
# 			) +
# 			scale_color_manual(
# 				values = c('0' = 'red', '1' = 'forestgreen')
# 			) +
# 			scale_fill_manual(
# 				name = 'Occupancy',
# 				values = c('0' = 'red', '1' = 'forestgreen'),
# 				labels = c('0' = 'Absent', '1' = 'Present')
# 			) +
# 			xlab(var_nice_units) +
# 			ggtitle(var_nice) +
# 			theme_minimal() +
# 			theme(
# 				axis.title.y = element_blank()
# 			)

# 	} # next variable

# 	hists <- plot_grid(plotlist = hists, nrow = 2)
# 	ggsave(hists, filename = './Figures & Tables/Microclimate Analysis - Single Year/Histograms of Predictors.png', width = 14, height = 8, bg = 'white')

# say('###########################################################################')
# say('### analyze patterns of occupancy for single-year data (Maria Westover) ###')
# say('###########################################################################')

# 	data <- fread('./Data/Microclimate/Single Year 2025-08-07 Marie Westover/iButton Data MLW matching site names 8.7.25.csv')
# 	names(data)[names(data) == 'chroniccoldavgNov-Mar'] <- 'chroniccoldavgNovMar'

# 	data$Hot22Days <- log10(1 + data$Hot22Days)
# 	data$numHomeRanges <- log10(1 + data$numHomeRanges)
# 	data$meanDistToClosest4Patches_m <- log10(data$meanDistToClosest4Patches_m)

# 	model_forms_raw <- read_xlsx('./Data/Microclimate/Single Year 2025-08-07 Marie Westover/PikaiButtonModelsMLW.xlsx', sheet = 'Sheet1')
# 	model_forms <- list()
# 	for (i in 1:nrow(model_forms_raw)) {
	
# 		terms <- model_forms_raw[i, ]
# 		terms <- unlist(terms)
# 		terms <- terms[!is.na(terms)]
# 		terms <- gsub(terms, pattern = '\\*', replacement = ':')

# 		model_forms[[i]] <- terms

# 	}
	
# 	model_forms_iso <- model_forms_hr <- model_forms_iso_hr <- model_forms

# 	for (i in seq_along(model_forms)) {
		
# 		model_forms_iso[[i]] <- c(model_forms_iso[[i]], 'meanDistToClosest4Patches_m')
# 		model_forms_hr[[i]] <- c(model_forms_hr[[i]], 'numHomeRanges')
# 		model_forms_iso_hr[[i]] <- c(model_forms_iso_hr[[i]], 'meanDistToClosest4Patches_m + numHomeRanges')

# 	}

# 	model_forms <- c(model_forms, model_forms_iso, model_forms_hr, model_forms_iso_hr)
# 	model_forms <- lapply(model_forms, paste, collapse = ' + ')
# 	model_forms <- c(model_forms, list('numHomeRanges', 'meanDistToClosest4Patches_m', 'numHomeRanges + meanDistToClosest4Patches_m'))
# 	for (i in seq_along(model_forms)) model_forms[[i]] <- c('binaryStatus ~ 1', model_forms[[i]])
# 	model_forms <- c(model_forms, list('binaryStatus ~ 1'))
# 	model_forms <- lapply(model_forms, paste, collapse = ' + ')
# 	model_forms <- lapply(model_forms, as.formula)

# 	# scale terms
# 	terms <- character()
# 	for (i in seq_along(model_forms)) {
# 		this_terms <- terms(model_forms[[i]])
# 		this_terms <- attr(this_terms, 'term.labels')
# 		terms <- c(terms, this_terms)
# 	}
# 	terms_with_ias <- unique(terms)
# 	terms_sans_ias <- terms_with_ias[!grepl(terms_with_ias, pattern = '\\:')]

# 	# table for terms
# 	coeff_table_empty <- as.data.table(t(as.data.table(terms_with_ias)))
# 	colnames(coeff_table_empty) <- terms_with_ias
# 	coeff_table_empty[] <- NA_real_
# 	coeff_table_empty[ , '(Intercept)' := NA_real_]

# 	data[ , (terms_sans_ias) := lapply(.SD, function(x) as.numeric(scale(x))), .SDcols = terms_sans_ias]
# 	data <- data[complete.cases(data[ , ..terms_sans_ias])]
# 	n <- nrow(data)

# 	# log likelihood for intercept-only model for calculation of Nagelkerke's R2
# 	intercept_only_model <- glm(binaryStatus ~ 1, data = data, family = binomial())
# 	intercept_only_model_lik <- exp(as.numeric(logLik(intercept_only_model)))

# 	results <- data.table()
# 	for (i in seq_along(model_forms)) {
	
# 		form <- model_forms[[i]]
# 		terms <- terms(form)
# 		terms <- attr(terms, 'term.labels')
# 		terms <- terms[!grepl(terms, pattern = '\\:')]

# 		model <- glm(form, data = data, family = binomial())
	
# 		coeffs <- coefficients(model)
# 		aicc <- AICc(model)
# 		ll <- as.numeric(logLik(model))
# 		lik <- exp(ll)
# 		nr2 <- nagelR2(likeNull = intercept_only_model_lik, likeFull = lik, n = n)

# 		form_char <- as.character(form)
# 		form_char <- form_char[form_char %notin% c('binaryStatus', '~', '1')]
# 		if (length(form_char) == 0) form_char <- '1'

# 		this_results <- data.table(
# 			model = form_char,
# 			aicc = aicc,
# 			delta_aicc = NA_real_,
# 			weight = NA_real_,
# 			pseudo_r2 = nr2
# 		)

# 		this_coeff_table <- coeff_table_empty
# 		for (count_term in seq_along(coeffs)) {

# 			term <- names(coeffs[count_term])
# 			this_coeff_table[1, (term)] <- coeffs[count_term]
		
# 		}

# 		this_results <- cbind(this_results, this_coeff_table)

# 		results <- rbind(results, this_results)

# 	}

# 	min_aicc <- min(results$aicc)
# 	results[ , delta_aicc := aicc - min_aicc]
# 	results[ , weight := exp(-0.5 * delta_aicc) / sum(exp(-0.5 * delta_aicc))]

# 	results <- results[order(weight, decreasing = TRUE)]

# 	### variable importance
# 	coeffs <- names(coeff_table_empty)
# 	coeffs <- coeffs[coeffs != '(Intercept)']

# 	var_imp <- data.table(
# 		coefficient = coeffs,
# 		aicc_weight_sum = 0,
# 		coeff_sum = 0,
# 		aicc_weighted_coeff = NA_real_,
# 		n_models = 0,
# 		n_models_positive = 0,
# 		n_models_negative = 0
# 	)

# 	for (i in 1:nrow(results)) {

# 		for (j in 1:nrow(var_imp)) {
		
# 			coeff <- var_imp$coefficient[j]
# 			val <- results[i, ..coeff]
# 			val <- val[[1]]

# 			if (!is.na(val)) {
			
# 				aicc_weight <- results$weight[i]

# 				var_imp[coefficient == coeff, 'aicc_weight_sum'] <- as.numeric(var_imp[coefficient == coeff, 'aicc_weight_sum']) + aicc_weight
# 				var_imp[coefficient == coeff, 'coeff_sum'] <- as.numeric(var_imp[coefficient == coeff, 'coeff_sum']) + val

# 				var_imp[coefficient == coeff, 'n_models'] <- as.numeric(var_imp[coefficient == coeff, n_models]) + 1
# 				if (val < 0) var_imp[coefficient == coeff, 'n_models_negative'] <- as.numeric(var_imp[coefficient == coeff, 'n_models_negative']) + 1
# 				if (val > 0) var_imp[coefficient == coeff, 'n_models_positive'] <- as.numeric(var_imp[coefficient == coeff, 'n_models_positive']) + 1
			
# 			}
		
# 		} # next variable

# 	} # next model

# 	var_imp[ , mean_aicc_weight := aicc_weight_sum / n_models]
# 	var_imp[ , aicc_weighted_coeff := coeff_sum / n_models]

# 	var_imp[ , 'coeff_sum' := NULL]
# 	var_imp <- var_imp[order(mean_aicc_weight, decreasing = TRUE)]

# 	fwrite(results, './Figures & Tables/Microclimate Analysis - Single Year/Microclimate Analysis - Single Year - Models.csv')
# 	fwrite(var_imp, './Figures & Tables/Microclimate Analysis - Single Year/Microclimate Analysis - Single Year - Variable Importance.csv')

# say('###################################################################')
# say('### exploratory data analysis for multi-year data (Erik Beever) ###')
# say('###################################################################')

# 	# Bandolier sites with microclimate over > 1 year
# 	band <- fread('./Data/Microclimate/Multi-Year 2025-10-23 Erik Beever & Dylan Ryals/band_clim_14.csv')
# 	band_sites <- tolower(band$patch)
# 	band_sites[band_sites %in% c('band101a', 'band101b')] <- 'band101'

# 	# add isolation and home range size of patches
# 	load('./Data/05 New Mexico Pika - Added PRISM Cell Number & Cell-Based Weight.rda')
# 	all_sites <- as.data.table(pika)
# 	all_sites_sites <- tolower(all_sites$polygonName)
# 	all_sites_sites <- gsub(all_sites_sites, pattern = ' ', replacement = '')

# 	matches <- match(band_sites, all_sites_sites)
# 	band[ , meanDistToClosest4Patches := all_sites$meanDistToClosest4Patches[matches]]
# 	band[ , numHomeRanges := all_sites$numHomeRanges[matches]]

# 	data <- band[ , c('patch', 'occy', 'subLethalCold', 'meltRefreeze', 'chronicCold', 'coldNoSnow', 'acuteCold', 'acuteHeat', 'subLethalHeat', 'chronicHeat', 'summerRespite', 'meanDistToClosest4Patches', 'numHomeRanges')]

# 	# # calculate site means
# 	# data <- data[ , lapply(.SD, mean, na.rm = TRUE), by = 'patch', .SDcols = setdiff(names(data), 'patch')]
# 	# data[ , patch := NULL]

# 	vars <- names(data)
# 	vars <- vars[vars %notin% c('occy', 'patch')]
# 	vars_nice <- make_nice_microclimate_vars_multi_year(vars, units = FALSE)
# 	vars_nice_units <- make_nice_microclimate_vars_multi_year(vars, units = TRUE)

# 	data[ , meltRefreeze := log10(meltRefreeze + 1)]
# 	data[ , coldNoSnow := log10(coldNoSnow + 1)]
# 	data[ , acuteCold := log10(acuteCold + 1)]
# 	data[ , acuteHeat := log10(acuteHeat + 1)]
# 	data[ , subLethalHeat := log10(subLethalHeat + 1)]
# 	data[ , subLethalHeat := log10(subLethalHeat + 1)]
# 	data[ , meanDistToClosest4Patches := log10(meanDistToClosest4Patches)]
# 	data[ , numHomeRanges := log10(numHomeRanges)]

# 	data$occy <- factor(data$occy)

# 	means <- data[ , lapply(.SD, mean, na.rm = TRUE), by = occy, .SDcols = vars]

# 	hists <- list()
# 	for (i in seq_along(vars)) {
	
# 		var <- vars[i]
# 		x <- data[[var]]
# 		data$var <- x

# 		var_nice <- vars_nice[i]
# 		var_nice_units <- vars_nice_units[i]

# 		this_means <- data.frame(occy = factor(c(0, 1)), xintercept = means[[var]])

# 		hists[[i]] <- ggplot(data, aes(x = var, fill = occy)) +
# 			geom_histogram(
# 				bins = 30, 
# 				color = NA, 
# 				position = 'stack', 
# 				alpha = 0.9
# 			) +
# 			geom_vline(
# 				data = this_means,
# 				aes(xintercept = xintercept, color = occy),
# 				lwd = 1.1,
# 				show.legend = FALSE
# 			) +
# 			scale_color_manual(
# 				values = c('0' = 'red', '1' = 'forestgreen')
# 			) +
# 			scale_fill_manual(
# 				name = 'Occupancy',
# 				values = c('0' = 'coral2', '1' = 'forestgreen'),
# 				labels = c('0' = 'Absent', '1' = 'Present')
# 			) +
# 			xlab(var_nice_units) +
# 			ggtitle(var_nice) +
# 			theme_minimal() +
# 			theme(
#             legend.position = 'none',
# 				axis.title.y = element_blank()
# 			)

# 	} # next variable

# 	hists <- plot_grid(plotlist = hists, nrow = 2)
# 	ggsave(hists, filename = './Figures & Tables/Microclimate Analysis - Multiple Years/Microclimate Analysis - Multiple Years - Histograms of Predictors.png', width = 14, height = 6, bg = 'white')

say('#######################################################################')
say('### analyze patterns of occupancy for multi-year data (Erik Beever) ###')
say('#######################################################################')

	# model terms
	terms <- sort(c('subLethalCold', 'subLethalHeat', 'acuteHeat', 'acuteCold', 'summerRespite', 'chronicCold', 'chronicHeat', 'coldNoSnow', 'meltRefreeze'))
	terms_univariate <- sort(c('subLethalCold', 'subLethalHeat', 'acuteHeat', 'acuteCold', 'summerRespite', 'chronicCold', 'chronicHeat', 'coldNoSnow', 'meltRefreeze'))

	terms <- c(terms, 'meanDistToClosest4Patches', 'numHomeRanges', 'occy_prev')
	terms_univariate <- c(terms_univariate, 'meanDistToClosest4Patches', 'numHomeRanges', 'occy_prev')

	# Bandolier sites with microclimate over > 1 year
	band <- fread('./Data/Microclimate/Multi-Year 2025-10-23 Erik Beever & Dylan Ryals/band_clim_14.csv')
	band_sites <- tolower(band$patch)
	band_sites[band_sites %in% c('band101a', 'band101b')] <- 'band101'

	# add isolation and home range size of patches
	load('./Data/05 New Mexico Pika - Added PRISM Cell Number & Cell-Based Weight.rda')
	all_sites <- as.data.table(pika)
	all_sites_sites <- tolower(all_sites$polygonName)
	all_sites_sites <- gsub(all_sites_sites, pattern = ' ', replacement = '')

	matches <- match(band_sites, all_sites_sites)
	band[ , meanDistToClosest4Patches := all_sites$meanDistToClosest4Patches[matches]]
	band[ , numHomeRanges := all_sites$numHomeRanges[matches]]

	# amalgamate BAND101A and BAND101B
	years <- sort(unique(band$year))
	for (year in years) {
	
		idx <- which(band$year == year & band$patch %in% c('BAND101A', 'BAND101B'))

		if (length(idx) > 1) {
		
			this <- band[idx]
			this_new <- this[1]

			this_new$patch <- 'BAND101'
			this_new$subLethalCold <- mean(this$subLethalCold, na.rm = TRUE)
			this_new$meltRefreeze <- mean(this$meltRefreeze, na.rm = TRUE)
			this_new$chronicCold <- mean(this$chronicCold, na.rm = TRUE)
			this_new$coldNoSnow <- mean(this$coldNoSnow, na.rm = TRUE)
			this_new$acuteCold <- mean(this$acuteCold, na.rm = TRUE)
			this_new$acuteHeat <- mean(this$acuteHeat, na.rm = TRUE)
			this_new$subLethalHeat <- mean(this$subLethalHeat, na.rm = TRUE)
			this_new$chronicHeat <- mean(this$chronicHeat, na.rm = TRUE)
			this_new$summerRespite <- mean(this$summerRespite, na.rm = TRUE)

			this_new$occy <- max(this$occy, na.rm = TRUE)

			band <- band[-idx]
			band <- rbind(band, this_new)

		}

	
	}

	band$patch[band$patch == 'BAND101A'] <- 'BAND101'

	# center and scale selected columns in 'band'
	cols_to_scale <- c('subLethalCold', 'meltRefreeze', 'chronicCold', 'coldNoSnow', 'acuteCold', 'acuteHeat', 'subLethalHeat', 'chronicHeat', 'summerRespite', 'meanDistToClosest4Patches', 'numHomeRanges')

	band[ , occy := as.numeric(occy)]

	# # remove rows with NA in any of the selected columns
	# band <- band[complete.cases(band[, ..cols_to_scale])]

	# dummy variable for AR1
	band[ , year_recode := year - min(year) + 1]
	band[ , year_recode := as.integer(year_recode)]

	# ### tally patch by year
	# #######################

	# # table with 0 = recorded patch status in that year, NA = did not
	# patches <- sort(unique(band$patch))
	# years <- sort(unique(band$year))

	# tallies <- data.table(patch = patches)
	# for (year in years) {
	
	# 	tallies[ , this_year := NA_integer_]
	# 	for (patch in patches) {
		
	# 		n <- sum(band$year == year & band$patch == patch)
	# 		if (n > 0) tallies$this_year[tallies$patch == patch] <- 0
		
	# 	}
	
	# 	setnames(tallies, old = 'this_year', new = paste0('year_', year))
	# 	# colnames(tallies)[ncol(tallies)] <- paste0('year_', year)
	# 	# n <- ncol(tallies)

	# }

	# # re-code non-NA values to 0/1 if detection or not
	# for (year in years) {
	# 	idx <- which(!is.na(tallies[[paste0('year_', year)]]))
	# 	for (i in idx) {

	# 		patch_name <- tallies$patch[i]
	# 		occy_val <- band$occy[band$patch == patch_name & band$year == year]
	# 		if (length(occy_val) > 0 && any(occy_val == 1, na.rm = TRUE)) {
	# 			tallies[[paste0('year_', year)]][i] <- 1
	# 		}
		
	# 	}
	# }

	# year_cols <- names(tallies)[names(tallies) != 'patch']
	# tallies[ , n_years_obs := rowSums(!is.na(.SD)), .SDcols = year_cols]

	# # # remove any patches that have just one year of observation
	# # if (any(tallies$n_years_obs <= 1)) {
	
	# # 	bads <- tallies$patch[tallies$n_years_obs <= 1]
	# # 	band <- band[patch %notin% bads]
	# # 	tallies <- tallies[patch %notin% bads]
	
	# # }

	# # add predictors
	# patches <- tallies$patch
	# for (i in seq_along(terms_univariate)) {
	
	# 	term <- terms_univariate[i]
	# 	tallies[ , DUMMY := NA_real_]
		
	# 	for (this_patch in patches) {
			
	# 		y <- band[band$patch == this_patch, ..term]
	# 		y <- y[[1]]
	# 		y <- mean(y)
	# 		tallies$DUMMY[tallies$patch == this_patch] <- y

	# 	}
		
	# 	setnames(tallies, 'DUMMY', term)

	# }

	# tallies <- tallies[order(meanDistToClosest4Patches)]
	# fwrite(tallies, './Figures & Tables/Microclimate Analysis - Multiple Years/Microclimate Analysis - Multiple Years - Tallies of Detections by Site and Year.csv')

	# model formulae
	model_forms <- c(
		'1',
		# '1 + subLethalCold + subLethalHeat',
		'subLethalHeat',
		# 'acuteHeat + acuteCold',
		'acuteHeat',
		# 'acuteCold + acuteHeat + acuteCold:acuteHeat',
		'subLethalHeat + summerRespite',
		'acuteHeat + summerRespite',
		# 'acuteHeat + subLethalCold',
		'summerRespite',
		'chronicCold',
		'acuteCold',
		'subLethalCold',
		'chronicHeat + subLethalHeat',
		'chronicHeat',
		'coldNoSnow',
		'coldNoSnow + subLethalCold',
		'coldNoSnow + acuteCold',
		'coldNoSnow + chronicCold',
		# 'coldNoSnow + summerRespite',
		# 'coldNoSnow + chronicHeat',
		# 'coldNoSnow + subLethalCold + chronicHeat',
		# 'coldNoSnow + subLethalCold + summerRespite',
		# 'subLethalHeat + subLethalCold + subLethalHeat:subLethalCold',
		# 'meltRefreeze + acuteHeat',
		# 'chronicHeat + subLethalCold'
		'meltRefreeze',
		# 'meltRefreeze + summerRespite',
		# 'meltRefreeze + chronicHeat',
		'meltRefreeze + acuteCold',
		'meltRefreeze + subLethalCold'
	)

	model_forms_season <- c(
		'1' = NA,
		# '1 + subLethalCold + subLethalHeat',
		'subLethalHeat' = 'summer',
		# 'acuteHeat + acuteCold',
		'acuteHeat' = 'summer',
		# 'acuteCold + acuteHeat + acuteCold:acuteHeat',
		'subLethalHeat + summerRespite' = 'summer',
		'acuteHeat + summerRespite' = 'summer',
		# 'acuteHeat + subLethalCold',
		'summerRespite' = 'summer',
		'chronicCold' = 'winter',
		'acuteCold' = 'winter',
		'subLethalCold' = 'winter',
		'chronicHeat + subLethalHeat' = 'summer',
		'chronicHeat' = 'summer',
		'coldNoSnow' = 'winter',
		'coldNoSnow + subLethalCold' = 'winter',
		'coldNoSnow + acuteCold' = 'winter',
		'coldNoSnow + chronicCold' = 'winter',
		# 'coldNoSnow + summerRespite',
		# 'coldNoSnow + chronicHeat',
		# 'coldNoSnow + subLethalCold + chronicHeat',
		# 'coldNoSnow + subLethalCold + summerRespite',
		# 'subLethalHeat + subLethalCold + subLethalHeat:subLethalCold',
		# 'meltRefreeze + acuteHeat',
		# 'chronicHeat + subLethalCold'
		'meltRefreeze' = 'winter',
		# 'meltRefreeze + summerRespite',
		# 'meltRefreeze + chronicHeat',
		'meltRefreeze + acuteCold' = 'winter',
		'meltRefreeze + subLethalCold' = 'winter'
	)

	control <- glmmTMBControl(optCtrl = list(iter.max = 2e4, maxfun = 1e5))
	# control <- glmmTMBControl(optCtrl = list(iter.max = 2e4), optimizer = optim, optArgs = list(method = 'BFGS'))
	# control <- glmmTMBControl(optCtrl = list(iter.max = 2e4), optimizer = optim, optArgs = list(method = 'L-BFGS-B'))
	# control <- glmmTMBControl(optCtrl = list(iter.max = 2e4), optimizer = 'bobyqa')

	# temporal autocorrelation formula... site-specific
	ar_form <- ~ ar1(year_recode + 0 | patch)

	# center/scale
	band[ , (cols_to_scale) := lapply(.SD, scale), .SDcols = cols_to_scale]

	# null model for Nagelkerke's R2
	null_model <- glmmTMB(occy ~ 1, family = binomial, data = band, control = control)
	intercept_only_model_ll <- as.numeric(logLik(null_model))
	intercept_only_model_lik <- exp(intercept_only_model_ll)

		# function to remember multi-year model
		remember_multi_year_model <- function(model, form, results, intercept_only_model_lik, site_intercept, tac, season) {

			# if (model$fit$convergence != 0) stop('Convergence failure on model ', form, '. Value ', model$fit$convergence, '.')
			if (model$fit$convergence != 0) return(results)

			# remember model
			n <- nrow(band)
			ll <- as.numeric(logLik(model))
			lik_mean <- ll / n
			lik <- exp(ll)

			nr2 <- nagelR2(likeNull = intercept_only_model_lik, likeFull = lik, n = n)

			# aicc <- AICc(model)

			# get deviance from fitted glmmTMB model
			deviance <- -2 * as.numeric(logLik(model))

			this_results <- data.table(
				model = form,
				season = season,
				site_intercept = site_intercept,
				# tac = tac,
				n_obs = n,
				log_like = ll,
				lik_mean = lik_mean,
				pseudo_r2 = nr2,
				deviance = deviance
			)

			# remember coefficients			
			for (count_term in seq_along(terms)) {
				this_results[ , DUMMY := NA_real_]
				colnames(this_results)[ncol(this_results)] <- terms[count_term]
			}

			coeffs <- fixef(model)$cond
			names(coeffs) <- gsub(names(coeffs), pattern = 'cond\\(\\(', replacement = '')
			names(coeffs) <- gsub(names(coeffs), pattern = 'cond\\(', replacement = '')
			names(coeffs) <- gsub(names(coeffs), pattern = '\\)', replacement = '')

			coeffs <- coeffs[names(coeffs) != 'Int']

			for (count_coeff in seq_along(coeffs)) {
				
				col <- which(names(this_results) == names(coeffs)[count_coeff])
				this_results[1, col] <- coeffs[count_coeff]

			}

			results <- rbind(results, this_results)
			results

		}

	results <- data.table()
	n_forms <- length(model_forms)
	for (previous_year_effect in c(FALSE, TRUE)) {
	
		for (i in 1:n_forms) {
		
			this_form <- model_forms[i]
			if (previous_year_effect) this_form <- paste0(this_form, ' + occy_prev')

			this_season <- model_forms_season[i]
			say(this_form)

			this_form <- paste('occy ~', this_form)
			this_form_iso <- paste(this_form, ' + meanDistToClosest4Patches')
			this_form_hr <- paste(this_form, ' + numHomeRanges')
			this_form_iso_hr <- paste(this_form, ' + meanDistToClosest4Patches + numHomeRanges')

			### models: no TAC and no site effects
			######################################

			site_intercept <- FALSE
			tac <- FALSE

			# model with no TAC, no site effect, no HR, no isolation
			form <- this_form
			model <- glmmTMB(as.formula(form),
				family = binomial,
				data = band,
				control = control
			)

			results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# model with no TAC, no site effect, no HR, YES isolation
			form <- this_form_iso
			model <- glmmTMB(as.formula(form),
				family = binomial,
				data = band,
				control = control
			)

			results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# model with no TAC, no site effect, YES HR, NO isolation
			form <- this_form_hr
			model <- glmmTMB(as.formula(form),
				family = binomial,
				data = band,
				control = control
			)

			results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# model with no TAC, no site effect, YES HR, YES isolation
			form <- this_form_iso_hr
			model <- glmmTMB(as.formula(form),
				family = binomial,
				data = band,
				control = control
			)

			results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			### models: no TAC and YES site effects
			#######################################

			site_intercept <- TRUE
			tac <- FALSE

			# model with no TAC, YES site effect, no HR, no isolation
			form <- this_form
			form <- paste(form, '+ (1 | patch)')
			model <- glmmTMB(as.formula(form),
				family = binomial,
				data = band,
				control = control
			)

			results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# model with no TAC, YES site effect, no HR, YES isolation
			form <- this_form_iso
			form <- paste(form, '+ (1 | patch)')
			model <- glmmTMB(as.formula(form),
				family = binomial,
				data = band,
				control = control
			)

			results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# model with no TAC, YES site effect, YES HR, NO isolation
			form <- this_form_hr
			form <- paste(form, '+ (1 | patch)')
			model <- glmmTMB(as.formula(form),
				family = binomial,
				data = band,
				control = control
			)

			results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# model with no TAC, YES site effect, YES HR, YES isolation
			form <- this_form_iso_hr
			form <- paste(form, '+ (1 | patch)')
			model <- glmmTMB(as.formula(form),
				family = binomial,
				data = band,
				control = control
			)

			results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# ### models: YES TAC and NO site effects
			# #######################################

			# site_intercept <- FALSE
			# tac <- TRUE

			# # model with YES TAC, no site effect, no HR, no isolation
			# form <- this_form
			# model <- glmmTMB(as.formula(form),
			# 	family = binomial,
			# 	data = band,
			# 	control = control,
			# 	dispformula = ar_form
			# )

			# results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# # model with YES TAC, no site effect, no HR, YES isolation
			# form <- this_form_iso
			# model <- glmmTMB(as.formula(form),
			# 	family = binomial,
			# 	data = band,
			# 	control = control,
			# 	dispformula = ar_form
			# )

			# results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# # model with YES TAC, no site effect, YES HR, NO isolation
			# form <- this_form_hr
			# model <- glmmTMB(as.formula(form),
			# 	family = binomial,
			# 	data = band,
			# 	control = control,
			# 	dispformula = ar_form
			# )

			# results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# # model with YES TAC, no site effect, YES HR, YES isolation
			# form <- this_form_iso_hr
			# model <- glmmTMB(as.formula(form),
			# 	family = binomial,
			# 	data = band,
			# 	control = control,
			# 	dispformula = ar_form
			# )

			# results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# ### models: YES TAC and YES site effects
			# ########################################

			# site_intercept <- TRUE
			# tac <- TRUE

			# # model with YES TAC, YES site effect, no HR, no isolation
			# form <- this_form
			# form <- paste(form, '+ (1 | patch)')
			# model <- glmmTMB(as.formula(form),
			# 	family = binomial,
			# 	data = band,
			# 	control = control,
			# 	dispformula = ar_form
			# )

			# results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# # model with YES TAC, YES site effect, no HR, YES isolation
			# form <- this_form_iso
			# form <- paste(form, '+ (1 | patch)')
			# model <- glmmTMB(as.formula(form),
			# 	family = binomial,
			# 	data = band,
			# 	control = control,
			# 	dispformula = ar_form
			# )

			# results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# # model with YES TAC, YES site effect, YES HR, NO isolation
			# form <- this_form_hr
			# form <- paste(form, '+ (1 | patch)')
			# model <- glmmTMB(as.formula(form),
			# 	family = binomial,
			# 	data = band,
			# 	control = control,
			# 	dispformula = ar_form
			# )

			# results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

			# # model with YES TAC, YES site effect, YES HR, YES isolation
			# form <- this_form_iso_hr
			# form <- paste(form, '+ (1 | patch)')
			# model <- glmmTMB(as.formula(form),
			# 	family = binomial,
			# 	data = band,
			# 	control = control,
			# 	dispformula = ar_form
			# )

			# results <- remember_multi_year_model(model = model, form = form, results, intercept_only_model_lik, site_intercept = site_intercept, tac = tac, season = this_season)

		} # next model form

	} # include previous year site status as predictor?

	results <- results[order(lik_mean, decreasing = TRUE)]
	# results[ , delta_aicc := aicc - min(aicc)]
	# results[ , weight := exp(-0.5 * delta_aicc) / sum(exp(-0.5 * delta_aicc))]

	fwrite(results, './Figures & Tables/Microclimate Analysis - Multiple Years/Microclimate Analysis - Multiple Years - Models.csv')

	### variable importance
	#######################

	var_imp <- data.table(
		coefficient = terms,
		n_models = NA_integer_,
		n_models_positive = NA_integer_,
		n_models_negative = NA_integer_,
		mean_rank = NA_real_
	)

	results$rank <- 1:nrow(results)
	for (i in seq_along(terms)) {
	
		term <- terms[i]
		col <- which(names(results) == term)
		idx <- which(complete.cases(results[ , ..col]))

		if (length(idx) > 0) {

			this_results <- results[idx]
			coeff <- this_results[[term]]

			var_imp$n_models[i] <- nrow(this_results)
			var_imp$n_models_positive[i] <- sum(coeff > 0)
			var_imp$n_models_negative[i] <- sum(coeff < 0)
			var_imp$mean_rank[i] <- mean(this_results$rank)

		}
		
	}

	var_imp <- var_imp[order(mean_rank, decreasing = FALSE)]

	fwrite(var_imp, './Figures & Tables/Microclimate Analysis - Multiple Years/Microclimate Analysis - Multiple Years - Variable Importance.csv')
 
# 	# fwrite(band, './Figures & Tables/Microclimate Analysis - Multiple Years/Microclimate Analysis - Multiple Years - Data Centered & Scaled.csv')

say('###################################################################')
say('### plots of variable importance for single and multi-year data ###')
say('###################################################################')

	vi_single <- fread('./Figures & Tables/Microclimate Analysis - Single Year/Microclimate Analysis - Single Year - Variable Importance.csv')
	vi_multi <- fread('./Figures & Tables/Microclimate Analysis - Multiple Years/Microclimate Analysis - Multiple Years - Variable Importance.csv')

	vi_single$var_nice <- make_nice_microclimate_vars_single_year(vi_single$coefficient, units = FALSE)
	vi_multi$var_nice <- make_nice_microclimate_vars_multi_year(vi_multi$coefficient, units = FALSE)

	### single year
	###############

	# ensure coefficients are ordered by mean_aicc_weight
	vi_single$coefficient <- factor(vi_single$coefficient, levels = vi_single$coefficient[order(vi_single$mean_aicc_weight, decreasing = FALSE)])

	vi_single$var_nice <- factor(vi_single$var_nice, levels = vi_single$var_nice[order(vi_single$mean_aicc_weight, decreasing = FALSE)])

	vi_single$bar_fill <- 'steelblue'
	vi_single$bar_fill[vi_single$var_nice %in% c('Number of home ranges', 'Isolation')] <- 'gray'
	vi_single$bar_fill[vi_single$var_nice %in% c('Chronic cold', 'Snow days', 'Minimum temperature', 'Min. temp × Snow days')] <- '#9fc5d1'
	vi_single$bar_fill[vi_single$var_nice %in% c('Hot days', 'July temperature', 'Chronic heat')] <- '#ff070b'

	xmax <- 1.1 * max(vi_single$mean_aicc_weight)

	single <- ggplot(vi_single, aes(x = mean_aicc_weight, y = var_nice)) +
		geom_bar(stat = 'identity', aes(fill = bar_fill)) +
		scale_fill_identity() +
		geom_text(
			aes(
				label = sprintf('%.2f', round(aicc_weight_sum, 2)), 
				x = mean_aicc_weight + 0.0005
			),
			hjust = 0,
			size = 5
		) +
		coord_cartesian(xlim = c(0, xmax)) +
		labs(
			x = 'Mean AICc Weight\n(higher = more evidence)',
			title = 'a) Microclimate variable importance:\n    Single year'
		) +
		theme(
			plot.title = element_text(size = 14),
			axis.title.x = element_text(size = 12),
			axis.title.y = element_blank(),
			axis.text.x = element_text(size = 10),
			axis.text.y = element_text(size = 12)
		)

	### multi-year
	##############
	
	# ensure coefficients are ordered by mean rank
	vi_multi$var_nice <- factor(vi_multi$var_nice, levels = vi_multi$var_nice[order(vi_multi$mean_rank, decreasing = TRUE)])

	vi_multi$bar_fill <- 'steelblue'
	vi_multi$bar_fill[vi_multi$var_nice %in% c('Number of home ranges', 'Isolation', 'Previous-year occupancy')] <- 'gray'
	vi_multi$bar_fill[vi_multi$var_nice %in% c('Acute cold', 'Melt/re-freeze', 'Cold, no snow', 'Chronic cold', 'Sub-lethal cold')] <- '#9fc5d1'
	vi_multi$bar_fill[vi_multi$var_nice %in% c('Acute heat', 'Summer respite', 'Chronic heat', 'Sub-lethal heat')] <- '#ff070b'

	xmax <- 1.17 * max(vi_multi$mean_rank)

	multi <- ggplot(vi_multi, aes(x = mean_rank, y = var_nice)) +
		geom_bar(stat = 'identity', aes(fill = bar_fill)) +
		scale_fill_identity() +
		geom_text(
			aes(
				label = sprintf('%.2f', round(mean_rank, 2)), 
				x = mean_rank + 3
			),
			hjust = 0,
			size = 5
		) +
		coord_cartesian(xlim = c(0, xmax)) +
		labs(
			x = 'Mean Rank\n(lower = more evidence)',
			title = 'b) Microclimate variable importance:\n    Multiple years'
		) +
		theme(
			plot.title = element_text(size = 14),
			axis.title.x = element_text(size = 12),
			axis.title.y = element_blank(),
			axis.text.x = element_text(size = 10),
			axis.text.y = element_text(size = 12)
		)

		# combo <- single + multi + plot_layout(ncol = 2) & theme(plot.margin = margin(5, 30, 5, 5))
		combo <- single + multi

	ggsave(combo, filename = './Figures & Tables/Microclimate Analysis/Microclimate Variable Importance.pdf', width = 12, height = 8, bg = 'white')

# say('#######################################################################')
# say('### joint exploratory data analysis for single- and multi year data ###')
# say('#######################################################################')

# 	# compare distributions of environmental predictors shared by single- and multi-year data sets

# 	### multi-year (Bandolier)

# 		multi <- fread('./Data/Microclimate/Multi-Year 2025-10-23 Erik Beever & Dylan Ryals/band_clim_14.csv')
# 		multi_sites <- tolower(multi$patch)
# 		multi_sites[multi_sites %in% c('multi101a', 'multi101b')] <- 'multi101'

# 		# add isolation and home range size of patches
# 		load('./Data/05 New Mexico Pika - Added PRISM Cell Number & Cell-Based Weight.rda')
# 		all_sites <- as.data.table(pika)
# 		all_sites_sites <- tolower(all_sites$polygonName)
# 		all_sites_sites <- gsub(all_sites_sites, pattern = ' ', replacement = '')

# 		matches <- match(multi_sites, all_sites_sites)
# 		multi[ , meanDistToClosest4Patches := all_sites$meanDistToClosest4Patches[matches]]
# 		multi[ , numHomeRanges := all_sites$numHomeRanges[matches]]

# 		multi <- multi[ , c('patch', 'occy', 'subLethalCold', 'meltRefreeze', 'chronicCold', 'coldNoSnow', 'acuteCold', 'acuteHeat', 'subLethalHeat', 'chronicHeat', 'summerRespite', 'meanDistToClosest4Patches', 'numHomeRanges')]

# 		vars <- names(multi)
# 		vars <- vars[vars %notin% c('occy', 'patch')]
# 		vars_nice <- make_nice_microclimate_vars_multi_year(vars, units = FALSE)
# 		vars_nice_units <- make_nice_microclimate_vars_multi_year(vars, units = TRUE)

# 		multi[ , meltRefreeze := log10(meltRefreeze + 1)]
# 		multi[ , coldNoSnow := log10(coldNoSnow + 1)]
# 		multi[ , acuteCold := log10(acuteCold + 1)]
# 		multi[ , acuteHeat := log10(acuteHeat + 1)]
# 		multi[ , subLethalHeat := log10(subLethalHeat + 1)]
# 		multi[ , subLethalHeat := log10(subLethalHeat + 1)]
# 		multi[ , meanDistToClosest4Patches := log10(meanDistToClosest4Patches)]
# 		multi[ , numHomeRanges := log10(numHomeRanges)]

# 		multi$occy <- factor(multi$occy)

# 	### single-year
		
# 		single <- fread('./Data/Microclimate/Single Year 2025-08-07 Marie Westover/iButton Data MLW matching site names 8.7.25.csv')
# 		model_forms <- read_xlsx('./Data/Microclimate/Single Year 2025-08-07 Marie Westover/PikaiButtonModelsMLW.xlsx', sheet = 'Sheet1')

# 		vars <- c(model_forms$term1, model_forms$term2, model_forms$term3)
# 		vars <- vars[!is.na(vars)]
# 		vars <- sort(unique(vars))
# 		vars <- vars[!grepl(vars, pattern = '\\*')]
# 		vars <- c(vars, 'numHomeRanges', 'meanDistToClosest4Patches_m')

# 		names(single)[names(single) == 'chroniccoldavgNov-Mar'] <- 'chroniccoldavgNovMar'
# 		vars[vars == 'chroniccoldavgNov-Mar'] <- 'chroniccoldavgNovMar'

# 		vars_nice <- make_nice_microclimate_vars_single_year(vars, units = FALSE)
# 		vars_nice_units <- make_nice_microclimate_vars_single_year(vars, units = TRUE)

# 		single$binaryStatus <- factor(single$binaryStatus)

# 	### plotting function

# 		joint_hist <- function(collated, means, title, xlab) {

# 			collated$status_time_span <- apply(collated[ , c('time_span', 'status')], 1, paste, collapse = ' - ')
# 			collated$status_time_span <- factor(collated$status_time_span )

# 			means <- collated[ , lapply(.SD, mean, na.rm = TRUE), by = status_time_span, .SDcols = 'x']
# 			means$status_time_span <- factor(means$status_time_span )

# 			out <- ggplot() +
# 				geom_histogram(
# 					data = collated,
# 					mapping = aes(x = x, fill = status_time_span),
# 					bins = 50, 
# 					color = NA, 
# 					position = 'stack', 
# 					alpha = 0.5
# 				) +
# 				geom_vline(
# 					data = means,
# 					aes(xintercept = x, color = status_time_span),
# 					lwd = 1.1,
# 					show.legend = FALSE
# 				) +
# 				scale_fill_manual(
# 					values = c('Multi-year - 0' = 'firebrick4', 'Multi-year - 1' = 'chartreuse4', 'Single year - 0' = 'firebrick1', 'Single year - 1' = 'chartreuse')
# 				) +
# 				scale_color_manual(
# 					values = c('Multi-year - 0' = 'firebrick4', 'Multi-year - 1' = 'chartreuse4', 'Single year - 0' = 'firebrick1', 'Single year - 1' = 'chartreuse')
# 				) +
# 				labs(
# 					title = title,
# 					fill = 'Dataset/\nOccurrence'
# 				) +
# 				xlab(xlab) +
# 				ylab('Number of patches') +
# 				theme_minimal()

# 			out

# 		}


# 	hists <- list()
# 	### chronic cold

# 		title <- 'Chronic cold'
# 		xlab <- 'Chronic cold (°C)'

# 		collated <- data.table(
# 			x = c(multi$chronicCold, single$chroniccoldavgNovMar),
# 			status = c(multi$occy, single$binaryStatus),
# 			time_span = c(rep('Multi-year', nrow(multi)), rep('Single year', nrow(single)))
# 		)

# 		hists[[length(hists) + 1]] <- joint_hist(collated = collated, means = means, title = title, xlab = xlab)

# 	### chronic heat

# 		title <- 'Chronic heat'
# 		xlab <- 'Chronic heat (°C)'

# 		collated <- data.table(
# 			x = c(multi$chronicHeat, single$ChronicHeatavgJunSept),
# 			status = c(multi$occy, single$binaryStatus),
# 			time_span = c(rep('Multi-year', nrow(multi)), rep('Single year', nrow(single)))
# 		)

# 		hists[[length(hists) + 1]] <- joint_hist(collated = collated, means = means, title = title, xlab = xlab)

# 	### patch size (home ranges)

# 		title <- 'Patch size'
# 		xlab <- 'Patch size (log10, home ranges)'

# 		collated <- data.table(
# 			x = c(multi$numHomeRanges, single$numHomeRanges),
# 			status = c(multi$occy, single$binaryStatus),
# 			time_span = c(rep('Multi-year', nrow(multi)), rep('Single year', nrow(single)))
# 		)
# 		collated[ , x := log10(x)]

# 		hists[[length(hists) + 1]] <- joint_hist(collated = collated, means = means, title = title, xlab = xlab)

# 	### isolation

# 		title <- 'Isolation'
# 		xlab <- 'Isolation (log10, km)'

# 		collated <- data.table(
# 			x = c(multi$ meanDistToClosest4Patches, single$ meanDistToClosest4Patches),
# 			status = c(multi$occy, single$binaryStatus),
# 			time_span = c(rep('Multi-year', nrow(multi)), rep('Single year', nrow(single)))
# 		)
# 		collated[ , x := log10(x)]

# 		hists[[length(hists) + 1]] <- joint_hist(collated = collated, means = means, title = title, xlab = xlab)


# 	hists <- plot_grid(plotlist = hists, nrow = 2)
# 	ggsave(hists, filename = './Figures & Tables/Microclimate Analysis/Microclimate Analysis - Single & Multiple Years - Joint Histograms of Predictors.png', width = 14, height = 10, bg = 'white')

say('DONE', level = 1)

