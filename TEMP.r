			
			figs <- list()
			
			xlim <- range(thisData$value)
			
			# all regions together
			title <- paste0(predNice, ': all regions')
			figs[[length(figs) + 1]] <- ggplot(data=thisData, aes(x=value, col=latestOccStatus, fill=latestOccStatus)) +
				geom_density(size=1) +
				scale_color_manual(
					labels = c('never', 'old', 'occupied'),
					values=c('0 never'='firebrick3', '1 old'='darkgoldenrod3', '2 occupied'='darkgreen')
				) +
				scale_fill_manual(
					labels = c('never', 'old', 'occupied'),
					values=alpha(c('0 never'='firebrick3', '1 old'='darkgoldenrod3', '2 occupied'='darkgreen'), 0.2)
				) +
				xlim(xlim[1], xlim[2]) +
				labs(title=title, subtitle=pred, x=predDescriptorUnit, y=NULL) +
				# geom_vline(data=mus, aes(xintercept=mu, color=latestOccStatus), linetype='dotted', size=1) +
				theme(
					# legend.position='none',
					plot.title=element_text(face='bold')
				) +
				guides(
					color=guide_legend(title='Status'),
					fill=guide_legend(title='Status')
				)

			# "never" by region
			title <- paste0(predNice, ': no evidence')
			thisThisData <- thisData[thisData$latestOccStatus == '0 never', ]
			figs[[length(figs) + 1]] <- ggplot(data=thisThisData, aes(x=value, col=region, fill=region)) +
				geom_density(size=1) +
				scale_color_manual(
					values=c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred')
				) +
				scale_fill_manual(
					values=alpha(c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred'), 0.2),
				) +
				xlim(xlim[1], xlim[2]) +
				labs(title=title, subtitle=pred, x=predDescriptorUnit, y=NULL) +
				theme(
					plot.title=element_text(face='bold')
				) +
				guides(
					color=guide_legend(title='Region'),
					fill=guide_legend(title='Region')
				)

			# "old" by region
			title <- paste0(predNice, ': old evidence')
			thisThisData <- thisData[thisData$latestOccStatus == '1 old', ]
			figs[[length(figs) + 1]] <- ggplot(data=thisThisData, aes(x=value, col=region, fill=region)) +
				geom_density(size=1) +
				scale_color_manual(
					values=c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred')
				) +
				scale_fill_manual(
					values=alpha(c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred'), 0.2),
				) +
				xlim(xlim[1], xlim[2]) +
				labs(title=title, subtitle=pred, x=predDescriptorUnit, y=NULL) +
				theme(
					plot.title=element_text(face='bold')
				) +
				guides(
					color=guide_legend(title='Region'),
					fill=guide_legend(title='Region')
				)

			# "occupied" by region
			title <- paste0(predNice, ': occupied')
			thisThisData <- thisData[thisData$latestOccStatus == '2 occupied', ]
			figs[[length(figs) + 1]] <- ggplot(data=thisThisData, aes(x=value, col=region, fill=region)) +
				geom_density(size=1) +
				scale_color_manual(
					values=c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred')
				) +
				scale_fill_manual(
					values=alpha(c('southwest'='darkgoldenrod3', 'southeast'='darkgreen', 'northwest'='navyblue', 'northeast'='darkred'), 0.2),
				) +
				xlim(xlim[1], xlim[2]) +
				labs(title=title, subtitle=pred, x=predDescriptorUnit, y=NULL) +
				theme(
					plot.title=element_text(face='bold')
				) +
				guides(
					color=guide_legend(title='Region'),
					fill=guide_legend(title='Region')
				)
