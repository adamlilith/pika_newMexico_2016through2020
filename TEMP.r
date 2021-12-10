# source('E:/Ecology/Drive/Research Active/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/TEMP.r')

	p <- ggplot(data=results, aes(x=modelNice, y=auc, col=region, pch=numHomeRanges)) +
		geom_point(size=2) +
		labs(shape='Number of Home\nRanges as Covariate', color='Region Covariate') +
		xlab(NULL) + ylab('AUC') +
		# ylim(ylim[1], ylim[2]) +
		theme(
			axis.text.y=element_text(size=8),
			legend.title=element_text(size=6),
			legend.text=element_text(size=6),
			legend.position=c(-1.1, 0.9),
			plot.margin = unit(c(0.5, 0.5, 0.5, 2), 'cm')
		) +
		coord_flip()
	
	pdf('./Figures & Tables/Occupancy - Simple Models/Occupancy - Simple Binary Models Cross-validation Results.pdf', width=8, height=10.5)
		print(p)
	dev.off()
