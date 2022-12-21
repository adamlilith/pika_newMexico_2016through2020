# source('C:/Ecology/Drive/Research/Pikas - New Mexico 2016-2020 (Erik Beever et al)/pika_newMexico_2016through2020/TEMP.r')

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
