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

# Plot a roc-like curve (except using absolute numbers of true
# positives and false negatives) for the given alignments (x) and
# mapping-quality ranking (mapq).
roc_plot <- function(
		correct,
		mapq,
		color="blue",
		first=T,
		xlim=NULL,
		ylim=NULL,
		main="MAPQ ROCs",
		xlab="Incorrect alignments",
		ylab="Correct alignments",
		type="l")
{
	tab <- roc_table(correct, mapq)
	if(first) {
		plot(cumsum(tab$incor), cumsum(tab$cor), col=color, type=type, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)
	} else {
		lines(cumsum(tab$incor), cumsum(tab$cor), col=color, type=type)
	}
}

# Plot two roc-like curves (except using absolute numbers of true
# positives and false negatives) for the given alignments (x) and the
# given two mapping-quality rankings (mapq1 and mapq2).
roc_2plot <- function(
	correct,
	mapq1,
	mapq2,
	colors=c("blue", "red"),
	sca=1,
	round_digits=2,
	mapq1_lab="MAPQ 1",
	mapq2_lab="MAPQ 2",
	main="MAPQ ROCs",
	xlim=NULL,
	ylim=NULL)
{
	if(round_digits > 0) {
		mapq1 <- round(mapq1, digits=round_digits)
		mapq2 <- round(mapq2, digits=round_digits)
	}
	xroc <- roc_table(correct, mapq1)
	xroc2 <- roc_table(correct, mapq2)
	mx_x <- max(max(cumsum(xroc$incor)), max(cumsum(xroc2$incor)))
	mn_x <- min(min(cumsum(xroc$incor)), min(cumsum(xroc2$incor)))
	mx_y <- max(max(cumsum(xroc$cor)), max(cumsum(xroc2$cor)))
	mn_y <- min(min(cumsum(xroc$cor)), min(cumsum(xroc2$cor)))
	if(is.null(xlim)) {
		xlim <- c(mn_x, mn_x + (mx_x - mn_x) * sca)
	}
	if(is.null(ylim)) {
		ylim <- c(mn_y + (mx_y - mn_y) * (1-sca), mx_y)
	}
	roc_plot(correct, mapq1, color=colors[1], first=T, xlim=xlim, ylim=ylim, main=main)
	roc_plot(correct, mapq2, color=colors[2], first=F)
	legend("bottomright", c(mapq1_lab, mapq2_lab), col=colors, pch=c(1, 1), lty=c(1, 1))
}

mapqize <- function(pcor) { return(-10.0 * log10(1.0 - pcor)) }

tab <- read.csv(args[1])
pdf(file=args[2])
roc_2plot(tab$correct == 1, round(mapqize(tab$pcor)), tab$mapq, mapq1_lab="Precicted", mapq2_lab="Aligner-reported")
dev.off()
