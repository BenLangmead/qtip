args <- commandArgs(trailingOnly = TRUE)

# Take table of alignments and parallel vector of mapping qualities and
# return a data frame containing columns for MAPQs (sorted in
# descending order), cumulative number of correct alignments, and
# cumulative number of incorrect alignments.
roc_table <- function(correct, mapq) {
	xord <- order(mapq)
	correct <- factor(correct[xord], levels=c(F, T))
	tab <- table(mapq[xord], correct)
	tab <- data.frame(mapq=as.numeric(rownames(tab)), cor=tab[,2], incor=tab[,1])
	tab <- tab[order(tab$mapq, decreasing=T),]
	tab$corCum <- cumsum(tab$cor)
	tab$incorCum <- cumsum(tab$incor)
	return(tab)
}

#
# Make a plot that bins alignments by nearest MAPQ and plots predicted versus
# actual MAPQ.
#
bucketError <- function(
	correct,
	mapq,
	col="red",
	binby=1,
	max_mapq=100,
	new_plot=T,
	main="")
{
	mapq <- round(pmin(mapq, max_mapq))
	if(binby > 1) {
		mapq <- trunc((mapq+1)/binby)*binby
	}
	tab <- roc_table(correct, mapq)
	mapqEx <- -10 * log10(tab$incor / (tab$cor + tab$incor))
	mapqEx <- round(pmin(mapqEx, max_mapq))
	mx <- max(max(tab$mapq), max(mapqEx))
	frac = ((tab$cor + tab$incor) / sum(tab$cor + tab$incor))
	lfrac = log(frac)
	lfrac <- lfrac - min(lfrac)
	lfrac <- lfrac / max(lfrac)
	if(new_plot) {
		plot(tab$mapq, mapqEx, xlim=c(0, mx), ylim=c(0, mx),
			col=rgb(1.0, 0.0, 0.0, lfrac), pch=16,
			xlab="Predicted MAPQ", ylab="Actual MAPQ",
			main=paste("Predicted versus actual MAPQ", main))
	} else {
		points(tab$mapq, mapqEx, col=rgb(1.0, 0.0, 0.0, frac), pch=21)
	}
	points(tab$mapq, mapqEx, col="red", pch=1)
	abline(0, 1)
}

mapqize <- function(pcor) { return(-10.0 * log10(1.0 - pcor)) }
tab <- read.csv(args[1])

pdf(file=args[2])
bucketError(tab$correct == 1, round(mapqize(tab$pcor)))
dev.off()

pdf(file=args[3])
bucketError(tab$correct == 1, round(mapqize(tab$pcor)), binby=5)
dev.off()

pdf(file=args[4])
bucketError(tab$correct == 1, round(mapqize(tab$pcor)), binby=10)
dev.off()
