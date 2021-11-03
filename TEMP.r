# source('E:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/TEMP.r')


			png(paste0('./Figures & Tables/OCCUPANCY Correlations between Occupancy Variables Using a ', occWindow, '-yr Window Spoke Plot.png'), width=1000, height=1000)

				par(oma=c(3, 10, 3, 3))
				spoke(pos=corr > thold, neg=corr < -thold, labels=niceVars, ltyNeg='solid', colNeg='red', main=title, cexLabel=1.6, cex.main=1.9)
				legend('bottomright', inset=-0.01, legend=c(paste0('Positive >', thold), paste0('Negative <-', thold)), lwd=1, col=c('black', 'red'), bty='n', cex=1.8)
				
			dev.off()
			