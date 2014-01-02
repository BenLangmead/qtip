args <- commandArgs(trailingOnly = TRUE)

plotDropRate <- function(pcor, cor, minus.log=T, pcor2=NULL, pdf.fn=NULL) {
	if(!is.null(pdf.fn)) {
		pdf(file=pdf.fn)
	}
	ordr <- order(jitterPcor(pcor), decreasing=T)
	if(minus.log) {
    	x <- seq(1, length(cor)) / length(cor)
    	maxx <- 11
    	labs <- sprintf("%0.3f", seq(0.0, 1, 0.1))
    } else {
    	x <- -log2(1.0 - seq(1, length(cor)) / length(cor))
    	maxx <- ceiling(log2(length(cor)))
    	labs <- sprintf("%0.5f", (2.0 ^ seq(-1, -maxx)))
    }
	y <- cumsum(!cor[ordr])
	plot(x, y, xlab="Drop rate", xaxt="n", ylab="Cumulative # incorrect alignments", type="l")
	if(!is.null(pcor2)) {
		ordr2 <- order(jitterPcor(pcor2), decreasing=T)
		y2 <- cumsum(!cor[ordr2])
		lines(x, y2, type="l", col="red")
	}
	axis(1, at=seq(1, maxx), labels=labs, las=2, cex.axis=0.7)
	legend("bottomright", c('Predicted MAPQ', 'Original MAPQ'), col=c('black', 'red'), lty=c(1, 1))
	if(!is.null(pdf.fn)) {
		dev.off()
	}
}

# I add jitter so that ties are broken randomly
jitterPcor <- function(pcor) {
	# Add jitter that is small compared to the smallest difference
	# between any two unequal pcors
	ordr <- order(pcor)
	df <- diff(pcor[ordr])
	mn <- min(df[df > 0])
	return(jitter(pcor, amount=mn/10.0))
}

plotSmoothedComparison <- function(pcor, cor, minus.log=T, pdf.fn=NULL, spar=0.5, nbootstrap=1) {
	if(!is.null(pdf.fn)) {
		pdf(file=pdf.fn)
	}
	ordr <- order(jitterPcor(pcor), decreasing=T)
	xfracs <- seq(1, length(cor)+1) / (length(cor)+1)
	xfracs <- xfracs[1:length(xfracs)-1]
	if(minus.log) {
    	x <- xfracs
    	maxx <- 11
    	labs <- sprintf("%0.3f", seq(0.0, 1, 0.1))
    } else {
    	x <- -log2(1.0 - xfracs)
    	maxx <- ceiling(log2(length(cor)))
    	labs <- sprintf("%0.5f", (2.0 ^ seq(-1, -maxx)))
    }
	y.pred <- pcor[ordr]
	plot(x, y.pred, xaxt="n", xlab="Drop rate", ylab="Correct", ylim=c(0, 1), col="white", type="l")
	first.below.50 <- which(pcor[ordr] < 0.99999)[1]
	first.below.40 <- which(pcor[ordr] < 0.99990)[1]
	first.below.30 <- which(pcor[ordr] < 0.99900)[1]
	first.below.20 <- which(pcor[ordr] < 0.99000)[1]
	first.below.10 <- which(pcor[ordr] < 0.90000)[1]
	# Draw colored segments
	lines(x[1:first.below.50-1], y.pred[1:first.below.50-1], col="dodgerblue4", type="l")
	lines(x[first.below.50:first.below.40], y.pred[first.below.50:first.below.40], col="orange", type="l")
	lines(x[first.below.40:first.below.30], y.pred[first.below.40:first.below.30], col="purple", type="l")
	lines(x[first.below.30:first.below.20], y.pred[first.below.30:first.below.20], col="blue", type="l")
	lines(x[first.below.20:first.below.10], y.pred[first.below.20:first.below.10], col="green", type="l")
	lines(x[first.below.10:length(pcor)], y.pred[first.below.10:length(pcor)], col="dodgerblue", type="l")
	# Draw X axis
	axis(1, at=seq(1, maxx), labels=labs, las=2, cex.axis=0.7)
	for(i in 1:nbootstrap) {
		sm <- sample(length(pcor), replace=T)
		ordr.bs <- order(jitterPcor(pcor[sm]), decreasing=T)
		y.act <- smooth.spline(cor[sm][ordr.bs], spar=spar)$y
		lines(x, y.act, col="red")
	}
	if(!is.null(pdf.fn)) {
		dev.off()
	}
}

plotSmoothedDifference <- function(pcor, cor, minus.log=T, pdf.fn=NULL, spar=0.5, nbootstrap=1) {
	if(!is.null(pdf.fn)) {
		pdf(file=pdf.fn)
	}
	ordr <- order(jitterPcor(pcor), decreasing=T)
	xfracs <- seq(1, length(cor)) / length(cor)
	if(minus.log) {
    	x <- xfracs
    	maxx <- 11
    	labs <- sprintf("%0.3f", seq(0.0, 1, 0.1))
    } else {
        x <- -log2(1.0 - xfracs)
        maxx <- ceiling(log2(length(cor)))
    	labs <- sprintf("%0.5f", (2.0 ^ seq(-1, -maxx)))
    }
	y <- pcor[ordr]
	ydiff <- y - smooth.spline(cor[ordr], spar=spar)$y
	plot(x, ydiff, xaxt="n", xlab="Drop rate", ylab="Correct", ylim=c(-1, 1), col="white")
	first.below.50 <- which(pcor[ordr] < 0.99999)[1]
	first.below.40 <- which(pcor[ordr] < 0.99990)[1]
	first.below.30 <- which(pcor[ordr] < 0.99900)[1]
	first.below.20 <- which(pcor[ordr] < 0.99000)[1]
	first.below.10 <- which(pcor[ordr] < 0.90000)[1]
	lines(x[1:first.below.50-1], ydiff[1:first.below.50-1], col="dodgerblue4", type="l")
	lines(x[first.below.50:first.below.40], ydiff[first.below.50:first.below.40], col="orange", type="l")
	lines(x[first.below.40:first.below.30], ydiff[first.below.40:first.below.30], col="purple", type="l")
	lines(x[first.below.30:first.below.20], ydiff[first.below.30:first.below.20], col="blue", type="l")
	lines(x[first.below.20:first.below.10], ydiff[first.below.20:first.below.10], col="green", type="l")
	lines(x[first.below.10:length(pcor)], ydiff[first.below.10:length(pcor)], col="dodgerblue", type="l")
	axis(1, at=seq(1, maxx), labels=labs, las=2, cex.axis=0.7)
	for(i in 1:nbootstrap) {
		sm <- sample(length(pcor), replace=T)
		ordr.bs <- order(jitterPcor(pcor[sm]), decreasing=T)
		ydiff <- pcor[sm][ordr.bs] - smooth.spline(cor[sm][ordr.bs], spar=spar)$y
		lines(x, ydiff, col="red")
	}
	if(!is.null(pdf.fn)) {
		dev.off()
	}
}

sdiff <- function(pcor, cor, spar=0.5, n=10, normalize=F) {
	errs <- c()
	for(i in 1:n) {
		ordr <- order(jitterPcor(pcor), decreasing=T)
		diffs <- smooth.spline(cor[ordr], spar=spar)$y - pcor[ordr]
		err = sum(diffs ^ 2)
		if(normalize) {
			ordr.null <- sample(length(pcor))
			diffs.null <- smooth.spline(cor[ordr.null], spar=spar)$y - pcor[ordr.null]
			err.null = sum(diffs.null ^ 2)
			errs <- append(errs, err / err.null)
		} else {
			errs <- append(errs, err)
		}
	}
	return(mean(errs))
}

rankerr <- function(pcor, cor, n=50, normalize=F) {
	errs <- c()
	for(i in 1:n) {
		ordr <- order(jitterPcor(pcor), decreasing=T)
		err <- sum(as.numeric(which(!cor[ordr])))
		if(normalize) {
			ordr.null <- sample(length(pcor))
			err.null = sum(as.numeric(which(!cor[ordr.null])))
			errs <- append(errs, err / err.null)
		} else {
			errs <- append(errs, err)
		}
		errs <- append(errs, err)
	}
	return(mean(errs))
}

tab <- read.csv(args[1], colClasses=c('numeric', 'numeric', 'integer', 'logical'))
sdiff.tab <- sdiff(tab$pcor, tab$correct)
rankerr.tab <- rankerr(tab$pcor, tab$correct)

if(all(!is.na(tab$correct))) {

    print(paste('sdiff=', format(sdiff.tab, scientific=T)))
    print(paste('rankerr=', format(rankerr.tab, scientific=T)))

    plotSmoothedDifference(tab$pcor, tab$correct, minus.log=F, paste0(args[1], '_diffcmp.pdf'))
    plotSmoothedDifference(tab$pcor, tab$correct, minus.log=T, paste0(args[1], '_diffcmplog.pdf'))

    plotSmoothedComparison(tab$pcor, tab$correct, minus.log=F, paste0(args[1], '_smoothcmp.pdf'))
    plotSmoothedComparison(tab$pcor, tab$correct, minus.log=T, paste0(args[1], '_smoothcmplog.pdf'))

    plotDropRate(tab$pcor, tab$correct, minus.log=F, pcor2=tab$orig, paste0(args[1], '_drop.pdf'))
    plotDropRate(tab$pcor, tab$correct, minus.log=T, pcor2=tab$orig, paste0(args[1], '_droplog.pdf'))

} else {
    print("Warning: no info about which alignments are correct")
}
