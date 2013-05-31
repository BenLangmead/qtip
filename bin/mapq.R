##
# mapq.R
#
# Load alignment and training data output by ts.py, and use it to predict
# better mapping qualities for the alignments.  We're trying to estimate pcor,
# the probability that the alignment is correct.  Mapping quality (mapq)
# equals -10 * log10(1 - pcor). 
#
# A mixed dataset of unpaired and paired-end alignments consist of alignments
# from these categories:
#
# 1. Unpaired alignments
# 2. Paired-end mates that aligned in an unpaired fashion
# 3. Paired-end mates that aligned in a discordant fashion
# 4. Paired-end mates that aligned in a concordant fashion
#
# Cases 1, 2 and 3
# ================
#
# One approach for categories 1, 2 and 3 is to lump them all together,
# discarding any paired-end information from the mates and treating them all
# as unpaired alignments.
#
# Other approaches are to treat 1 separately from 2 and 3, or to treat all
# three separately.
#
# Case 4
# ======
#
# One approach for category 4 is to obtain estimates for pcor values for the
# two mates separately using the model fit for case 1 (perhaps also including
# cases 2 and/or 3).  Then do a query in a third model where the features are
# (1) pcor for mate 1, (2) pcor for mate 2, and (3) fragment length.
#
#
# TODO:
# - Weights for logistic regression
#

library(bitops)
library(zoo)
library(plyr)

myblues <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6")

# Given data frame of alignments, return vector of booleans indicating which aligned
selectAligned <- function(x) { return(bitAnd(x$flag, 4) == 0) }

# Given data frame of alignments, return vector of booleans indicating which are unpaired
selectUnpaired <- function(x) { return(x$YT.Z == "UU") }

# Given data frame of alignments, return vector of booleans indicating which are concordant
selectConcordant <- function(x) { return(x$YT.Z == "CP") }

# Take table of alignments and parallel vector of mapping qualities and
# return a data frame containing columns for MAPQs (sorted in
# descending order), cumulative number of correct alignments, and
# cumulative number of incorrect alignments.
roc_table <- function(correct, mapq) {
	xord <- order(mapq)
	correct <- factor(correct[xord], levels=c(0, 1))
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
	xroc <- roc_table(correct, mapq1, correct=correct)
	xroc2 <- roc_table(correct, mapq2, correct=correct)
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
	err1 <- rankingError(correct, mapq1, correct=correct)
	err2 <- rankingError(correct, mapq2, correct=correct)
	mapq1_lab <- paste0(mapq1_lab, " (err=", err1 / err1, ")")
	mapq2_lab <- paste0(mapq2_lab, " (err=", err2 / err1, ")")
	legend("bottomright", c(mapq1_lab, mapq2_lab), col=colors, pch=c(1, 1), lty=c(1, 1))
}

# Penalty is sum of ranks of all incorrectly aligned reads.  Rank is lower for
# reads with lower MAPQ.
rankingError <- function(correct, mapq, correct=NULL) {
	ordr <- order(mapq)
	return(sum(as.numeric(which(correct[ordr] == 0))))
}

quantileError <- function(correct, mapq, nquantile=10) {
	mxmapq <- max(mapq)
	xeq = x[x$predMapq == mxmapq,]
	xlt = x[x$predMapq < mxmapq,]
	ord <- order(xlt$predPcor)
	xlt <- xlt[ord,]
	ct <- cut(xlt$predPcor, nquantile, labels=F)
	pcorOld <- 1 - (10 ^ (-0.1 * xlt$mapq[ord]))
	pred <- c()
	actual <- c()
	predOld <- c()
	actualOld <- c()
	mnPcor <- c()
	mxPcor <- c()
	mnMapq <- c()
	mxMapq <- c()
	for(i in 1:nquantile) {
		n <- sum(ct == i)
		predOld <- append(predOld, sum(pcorOld[ct == i] / n))
		pred <- append(pred, sum(xlt$predPcor[ct == i]) / n)
		actual <- append(actual, sum(xlt$ZC.i[ct == i]) / n)
		mnMapq <- append(mnMapq, -10.0 * log10(1.0 - min(xlt$predPcor[ct == i])))
		mxMapq <- append(mxMapq, -10.0 * log10(1.0 - max(xlt$predPcor[ct == i])))
		mnPcor <- append(mnPcor, min(xlt$predPcor[ct == i]))
		mxPcor <- append(mxPcor, max(xlt$predPcor[ct == i]))
	}
	n <- sum(x$predMapq == mxmapq)
	predOld <- append(predOld, sum( 1 - (10 ^ (-0.1 * x$mapq[x$predMapq == mxmapq])) ) / n)
	pred <- append(pred, sum(xeq$predPcor) / n)
	actual <- append(actual, sum(xeq$ZC.i) / n)
	mnMapq <- append(mnMapq, mxmapq)
	mxMapq <- append(mxMapq, mxmapq)
	mnPcor <- append(mnPcor, max(xlt$predPcor))
	mxPcor <- append(mxPcor, max(xlt$predPcor))
	predPhred <- pmin(-10 * log10(1.0 - pred), 99.9, na.rm=T)
	predOldPhred <- pmin(-10 * log10(1.0 - predOld), 99.9, na.rm=T)
	actualPhred <- pmin(-10 * log10(1.0 - actual), 99.9, na.rm=T)
	mx <- max(max(predPhred[is.finite(predPhred)]), max(actualPhred[is.finite(actualPhred)]))
	if(T) {
		{
			quartz(width=6, height=6)
			ymn <- min(min(pred), min(predOld), min(actual))
			ymx <- max(max(pred), max(predOld), max(actual))
			plot(seq(1, nquantile+1), pred, type="o", col="red", ylim=c(ymn, ymx), xlab="", ylab="Fraction correct", main="Fraction predicted/acutally correct by decile", xaxt='n')
			axis(1, at=seq(1, nquantile+1), labels=paste0("[", sprintf("%0.2f", mnPcor), ", ", sprintf("%0.2f", mxPcor), "]"), las=2, cex.axis=0.75)
			points(seq(1, nquantile+1), predOld, type="o", col="blue")
			points(seq(1, nquantile+1), actual, type="o", col="green")
			legend("bottomright", lty=c(1, 1, 1), col=c("red", "blue", "green"), pch=c(1, 1, 1), legend=c("SimQ predicted", "Bowtie 2 predicted", "Actual"))
		}
		{
			quartz(width=6, height=6)
			ymn <- min(min(predPhred), min(predOldPhred), min(actualPhred))
			ymx <- max(max(predPhred), max(predOldPhred), max(actualPhred))
			plot(seq(1, nquantile+1), predPhred, type="o", col="red", ylim=c(ymn, ymx), xlab="", ylab="MAPQ", main="MAPQ predicted/acutal by decile", xaxt='n')
			axis(1, at=seq(1, nquantile+1), labels=paste0("[", sprintf("%0.2f", mnMapq), ", ", sprintf("%0.2f", mxMapq), "]"), las=2, cex.axis=0.75)
			points(seq(1, nquantile+1), predOldPhred, type="o", col="blue")
			points(seq(1, nquantile+1), actualPhred, type="o", col="green")
			legend("topleft", lty=c(1, 1, 1), col=c("red", "blue", "green"), pch=c(1, 1, 1), legend=c("SimQ predicted", "Bowtie 2 predicted", "Actual"))
		}
	} else {
		{
			quartz(width=6, height=6)
			ymn <- min(min(pred), min(actual))
			ymx <- max(max(pred), max(actual))
			plot(seq(1, nquantile+1), pred, type="o", col="red", ylim=c(ymn, ymx), xlab="", ylab="Fraction correct", main="Fraction predicted/acutally correct by decile", xaxt='n')
			axis(1, at=seq(1, nquantile+1), labels=paste0("[", sprintf("%0.2f", mnPcor), ", ", sprintf("%0.2f", mxPcor), "]"), las=2, cex.axis=0.75)
			points(seq(1, nquantile+1), actual, type="o", col="blue")
			legend("bottomright", lty=c(1, 1), col=c("red", "blue"), pch=c(1, 1), legend=c("SimQ predicted", "Actual"))
		}
		{
			quartz(width=6, height=6)
			ymn <- min(min(predPhred), min(actualPhred))
			ymx <- max(max(predPhred), max(actualPhred))
			plot(seq(1, nquantile+1), predPhred, type="o", col="red", ylim=c(ymn, ymx), xlab="", ylab="MAPQ", main="MAPQ predicted/acutal by decile", xaxt='n')
			axis(1, at=seq(1, nquantile+1), labels=paste0("[", sprintf("%0.2f", mnMapq), ", ", sprintf("%0.2f", mxMapq), "]"), las=2, cex.axis=0.75)
			points(seq(1, nquantile+1), actualPhred, type="o", col="blue")
			legend("topleft", lty=c(1, 1), col=c("red", "blue"), pch=c(1, 1), legend=c("SimQ predicted", "Actual"))
		}
	}
}

bucketError <- function(correct, mapq) {
	tab <- roc_table(correct, round(mapq))
	mapqEx <- -10 * log10(tab$incor / (tab$cor + tab$incor))
	mx <- max(max(tab$mapq), mapqEx[is.finite(mapqEx)])
	plot(tab$mapq, mapqEx, xlim=c(0, mx), ylim=c(0, mx))
	points(seq(1, mx), seq(1, mx), col="red")
}

openTrainingAndSam <- function(x, merge123=T) {
	
	options(stringsAsFactors = FALSE)
	
	# Read SAM output
	print(paste("  Scanning sat file", format(Sys.time(), "%X"), "..."))
	samTsvFn <- paste0(x, ".sat")
	samTab <- read.table(samTsvFn, comment.char="", quote="", header=T)
	if(! "len" %in% colnames(samTab)) {
		if(! "seq" %in% colnames(samTab)) {
			stop("no way to get read length")
		}
		samTab$len <- nchar(samTab$seq)
	}
	ret <- list()
	
	# Compile alignment data frames
	print(paste("  Splitting sat table by category", format(Sys.time(), "%X"), "..."))
	al <- samTab[selectAligned(samTab),]
	ret$alU <- al[selectUnpaired(al),]
	ret$alM <- al[!selectUnpaired(al) & !selectConcordant(al),]
	ret$alC <- al[selectConcordant(al),]
	
	# Read training data
	trainUFn <- paste0(x, ".unp.training.tsv")
	if(file.exists(trainUFn)) {
		print(paste("  Reading unpaired training data", format(Sys.time(), "%X"), "..."))
		ret$trainU <- read.table(trainUFn, comment.char="", quote="", header=T)
	}
	trainMFn <- paste0(x, ".m.training.tsv")
	if(file.exists(trainMFn)) {
		print(paste("  Reading paired-end M training data", format(Sys.time(), "%X"), "..."))
		ret$trainM <- read.table(trainMFn, comment.char="", quote="", header=T)
	}
	trainCFn <- paste0(x, ".conc.training.tsv")
	if(file.exists(trainCFn )) {
		print(paste("  Reading paired-end concordant training data", format(Sys.time(), "%X"), "..."))
		ret$trainC <- read.table(trainCFn , comment.char="", quote="", header=T)
	}
	return(ret)
}

##
# Stratify given data by given column.  Merge strata so that no two encompass
# identical intervals of data.
#
stratify <- function(data, column, strata=NULL) {
	if(is.null(strata)) {
		qs <- quantile(data[,column], probs=seq(0, 1, 0.1))
		strata <- qs[diff(qs) > 0]
		if(length(strata) == 0) {
			strata <- c(-Inf)
		}
	}
	data$stratum <- findInterval(data[,column], strata)
	data$stratum <- pmax(data$stratum, 1)
	return(list(data=data, strata=strata))
}

fit2DPlot <- function(
	fit,
	model,
	Xorig,
	Yorig,
	correct,
	Xex=NULL,
	Yex=NULL,
	alpha=NULL,
	main="",
	xlab="",
	ylab="",
	plot.fn=NULL)
{
	mnx <- min(Xorig); mxx <- max(Xorig)
	mny <- min(Yorig); mxy <- max(Yorig)
	if(!is.null(Xex)) {
		mnx <- min(mnx, min(Xex)); mxx <- max(mxx, max(Xex))
		mny <- min(mny, min(Yex)); mxy <- max(mxy, max(Yex))
	}
	Xs <- seq(mnx, mxx, length.out=200)
	Ys <- seq(mny, mxy, length.out=200)
	grid <- expand.grid(list(X=Xs, Y=Ys))
	if(model == "loess") {
		p <- pmax(pmin(predict(fit, newdata=grid), 1.0), 0.0)
	} else {
		p <- pmax(pmin(predict(fit, newdata=grid, type="response"), 1.0), 0.0)
	}
	p <- -10 * log10(1.0 - p)
	p <- pmin(p, 100.0)
	p <- matrix(p, nrow=200)
	if(!is.null(plot.fn)) {
		pdf(file=plot.fn, width=6, height=6)
	} else {
		quartz(width=6, height=6)
	}
	image(
		Xs, Ys, p, col=myblues, zlim=c(0, 100.0),
		main=paste0(main, ": ", model, " model"),
		xlab=xlab, ylab=ylab)
	jit <- function(x, factor=2) { jitter(x, factor=factor) }
	points(jit(Xorig[correct]), jit(Yorig[correct]), col=rgb(0.0, 0.0, 0.0, 0.1))
	points(jit(Xorig[!correct]), jit(Yorig[!correct]), col="red")
	if(!is.null(Xex)) {
		if(!is.null(alpha)) {
			points(jit(Xex), jit(Yex), col=rgb(0.0, 0.0, 1.0, alpha))
		} else {
			points(jit(Xex), jit(Yex), col=rgb(0.0, 0.0, 1.0))
		}
	}
	if(!is.null(plot.fn)) { dev.off() }
}

fit1DPlot <- function(
	fit,
	model,
	Xorig,
	correct,
	Xex=NULL,
	Yex=NULL,
	alpha=NULL,
	main="",
	xlab="",
	plot.fn=NULL)
{
	mnx <- min(Xorig); mxx <- max(Xorig)
	if(!is.null(Xex)) {
		mnx <- min(mnx, min(Xex)); mxx <- max(mxx, max(Xex))
	}
	Xs <- seq(mnx, mxx, length.out=200)
	Ys <- seq(0, 2, length.out=200)
	grid <- expand.grid(list(X=Xs, Y=Ys))
	if(model == "loess") {
		p <- pmax(pmin(predict(fit, newdata=grid), 1.0), 0.0)
	} else {
		p <- pmax(pmin(predict(fit, newdata=grid, type="response"), 1.0), 0.0)
	}
	p <- -10 * log10(1.0 - p)
	lim <- 50.0
	p <- pmin(p, lim)
	p <- matrix(p, nrow=200)
	if(!is.null(plot.fn)) {
		pdf(file=plot.fn, width=6, height=6)
	} else {
		quartz(width=6, height=6)
	}
	image(
		Xs, Ys, p, col=myblues, zlim=c(0, lim),
		main=paste0(main, ": ", model, " model"),
		yaxt='n', ylab="", xlab=xlab)
	jit <- function(x, factor=2) { jitter(x, factor=factor) }
	points(jit(Xorig[correct]), jit(rep(1, length(Xorig[correct])), factor=15), col=rgb(0.0, 0.0, 0.0, 0.1))
	points(jit(Xorig[!correct]), jit(rep(1, length(Xorig[!correct])), factor=15), col="red")
	if(!is.null(Xex)) {
		if(!is.null(alpha)) {
			points(jit(Xex), jit(rep(1, length(Xex)), factor=15), col=rgb(0.0, 0.0, 1.0, alpha))
		} else {
			points(jit(Xex), jit(rep(1, length(Xex)), factor=15), col=rgb(0.0, 0.0, 1.0))
		}
	}
	if(!is.null(plot.fn)) { dev.off() }
}

fitUnpaired <- function(
	train,               # unpaired training matrix
	stratum=1,
	model="loess",
	span=1.0,            # loess span
	degree=1.0,          # loess degree
	scaleXrep=NULL,
	scaleYrep=NULL,
	scaleXuni=NULL,
	use.sqrt=F,
	do.plot=F,
	plot.fn=NULL)
{
	train$diffValid <- train$maxValid - train$minValid
	train$bestNorm <- (train$best - train$maxValid) / train$diffValid

	trRep <- train[train$diff < 999999,]
	trUni <- train[train$diff == 999999,]
	
	trRep$diffNorm <- trRep$diff / trRep$diffValid
	
	#
	# Fit model for Case 1: repetitive
	#
	summ <- ddply(trRep, .(bestNorm, diffNorm), summarise, mean=mean(correct), n=length(correct))
	# Gather Xs (best scores) and Ys (differences)
	X <- summ$bestNorm
	mxxRep <- max(X); mnxRep <- min(X)
	if(!is.null(scaleXrep)) { X <- scaleXrep(X, mnxRep, mxxRep) }
	Y <- summ$diffNorm
	mxyRep <- max(Y); mnyRep <- min(Y)
	if(!is.null(scaleYrep)) { Y <- scaleYrep(Y, mnyRep, mxyRep) }
	weights <- ifelse(rep(use.sqrt, nrow(summ)), sqrt(summ$n), summ$n)
	if(model == "loess") {
		fitRep <- loess(summ$mean ~ X * Y, span=span, degree=degree, weight=weights)
		pcorRep <- predict(fitRep, newdata=data.frame(X))
	} else if(model == "logit") {
		fitRep <- glm(summ$mean ~ X * Y, family=binomial("logit"), weight=weights)
		pcorRep <- predict(fitRep, newdata=data.frame(X), type="response")
	}
	pcorRep <- pmin(pmax(pcorRep, 0.0), 1.0)
	mapqRep <- -10.0 * log10(1.0 - pcorRep)
	XorigRep <- trRep$bestNorm
	YorigRep <- trRep$diffNorm
	correctRep <- trRep$correct == 1
	if(!is.null(scaleXrep)) { XorigRep <- scaleXrep(XorigRep, mnxRep, mxxRep) }
	if(!is.null(scaleYrep)) { YorigRep <- scaleYrep(YorigRep, mnyRep, mxyRep) }
	
	plotRep <- function(Xex=NULL, Yex=NULL, alpha=NULL) {
		my.plot.fn <- plot.fn
		if(!is.null(plot.fn)) {
			my.plot.fn <- paste0(plot.fn, "unp.", stratum, ".fit.rep.pdf")
		}
		fit2DPlot(
			fitRep,
			model,
			XorigRep,
			YorigRep,
			correctRep,
			Xex=Xex,
			Yex=Yex,
			alpha=alpha,
			xlab="Best aln score",
			ylab="Absolute diff b/t best, 2nd best aln score",
			main="Unpaired case 1",
			plot.fn=my.plot.fn)
		if(!is.null(plot.fn)) {
			pdf(file=paste0(plot.fn, "unp.", x, ".trainroc.rep.pdf"), width=6, height=6)
		} else {
			quartz(width=6, height=6)
		}
		roc_plot(trRep$correct, mapqRep, main="Unpaired case 1 (rep) training-data ROC")
		if(!is.null(plot.fn)) { dev.off() }
	}
	if(do.plot) { plotRep() }
	
	#
	# Fit model for Case 2: unique
	#
	summ <- ddply(trUni, .(bestNorm), summarise, mean=mean(correct), n=length(correct))
	X <- summ$bestNorm
	mxxUni <- max(X); mnxUni <- min(X)
	if(!is.null(scaleXuni)) { X <- scaleXuni(X, mnxUni, mxxUni) }
	weights <- ifelse(rep(use.sqrt, nrow(summ)), sqrt(summ$n), summ$n)
	if(model == "loess") {
		fitUni <- loess(summ$mean ~ X, span=span, degree=degree, weight=weights)
		pcorUni <- predict(fitUni, newdata=data.frame(X))
	} else if(model == "logit") {
		fitUni <- glm(summ$mean ~ X, family=binomial("logit"), weight=weights)
		pcorUni <- predict(fitUni, newdata=data.frame(X), type="response")
	}
	pcorUni <- pmin(pmax(pcorUni, 0.0), 1.0)
	mapqUni <- -10.0 * log10(1.0 - pcorUni)
	XorigUni <- trUni$bestNorm
	if(!is.null(scaleXuni)) { XorigUni <- scaleXuni(XorigUni, mnxUni, mxxUni) }
	correctUni <- trUni$correct == 1
	plotUni <- function(Xex=NULL, alpha=NULL) {
		my.plot.fn <- plot.fn
		if(!is.null(plot.fn)) {
			my.plot.fn <- paste0(plot.fn, "unp.", stratum, ".fit.uni.pdf")
		}
		fit1DPlot(
			fitUni,
			model,
			XorigUni,
			correctUni,
			Xex=Xex,
			alpha=alpha,
			xlab="Best aln score",
			main="Unpaired case 2",
			plot.fn=my.plot.fn)
		if(!is.null(plot.fn)) {
			pdf(file=paste0(plot.fn, "unp.", x, ".trainroc.uni.pdf"), width=6, height=6)
		} else {
			quartz(width=6, height=6)
		}
		roc_plot(trUni$correct, mapqUni, main="Unpaired case 2 (uni) training-data ROC")
		if(!is.null(plot.fn)) { dev.off() }
	}
	if(do.plot) { plotUni() }
	
	return(list(
		rep=fitRep,
		repBounds=c(mnxRep, mxxRep, mnyRep, mxyRep),
		repX=XorigRep,
		repY=YorigRep,
		repCorrect=correctRep,
		uni=fitUni,
		uniBounds=c(mnxUni, mxxUni),
		uniX=trUni$bestNorm,
		uniCorrect=correctUni,
		model=model,
		plotRep=plotRep,
		plotUni=plotUni))
}

fitUnpairedStrata <- function(
	train,
	model="loess",
	span=2.0,
	scaleXrep=NULL,
	scaleYrep=NULL,
	scaleXuni=NULL,
	use.sqrt=F,
	do.plot=F,
	plot.fn=NULL)
{
	return(lapply(unique(train$stratum), function(x) {
		fitUnpaired(
			train[train$stratum == x,],
			stratum=x,
			model=model,
			span=span,
			scaleXrep=scaleXrep,
			scaleYrep=scaleYrep,
			scaleXuni=scaleXuni,
			use.sqrt=use.sqrt,
			do.plot=do.plot,
			plot.fn=plot.fn)
	}))
}

##
# Fit concordantly-aligned training data using three models covering
# three cases:
#
# 1. There was more than one concordant alignment and more than one
#    alignment for the mate whose MAPQ is being estimated
# 2. There was one concordant alignment but more than one alignment for
#    the mate whose MAPQ is being estimated
# 3. There was just one concordant alignment and just one alignment for
#    the mate whose MAPQ is being estimated
#
fitPaired <- function(
	train,               # training
	stratum=1,
	model="loess",
	span=1.0,            # loess span
	degree=1.0,          # loess degree
	scaleXcrep=NULL,
	scaleYcrep=NULL,
	scaleXcuni=NULL,
	scaleYcuni=NULL,
	scaleXuni=NULL,
	use.sqrt=F,
	do.plot=F,
	plot.fn=NULL)
{
	train$diffValidA  <- train$maxValidA - train$minValidA
	train$diffValidCA <- train$diffValidA + (train$maxValidB - train$minValidB)
	train$bestANorm   <- (train$bestA - train$maxValidA) / train$diffValidA
	train$bestCANorm  <- ((train$bestA + train$bestB) - (train$maxValidA + train$maxValidB)) / train$diffValidCA
	
	trCrep <- train[train$diffCA  < 999999,]
	trCuni <- train[train$diffCA == 999999 & train$diffA  < 999999,]
	trUni  <- train[train$diffCA == 999999 & train$diffA == 999999,]
	
	trCrep$diffANorm  <- ifelse(trCrep$diffA == 999999, 1.1 * trCrep$diffValidA, pmin(trCrep$diffA, trCrep$diffValidA))  / trCrep$diffValidA
	trCrep$diffCANorm <- trCrep$diffCA / trCrep$diffValidCA
	
	trCuni$diffANorm  <- trCuni$diffA  / trCuni$diffValidA
	
	#
	# Fit model for Case 1: repetitive
	#
	summ <- ddply(trCrep, .(diffANorm, diffCANorm), summarise, mean=mean(correctA), n=length(correctA))
	# Gather Xs (best scores) and Ys (differences)
	X <- summ$diffANorm
	mxxCrep <- max(X); mnxCrep <- min(X)
	if(!is.null(scaleXcrep)) { X <- scaleXcrep(X, mnxCrep, mxxCrep) }
	Y <- summ$diffCANorm
	mxyCrep <- max(Y); mnyCrep <- min(Y)
	if(!is.null(scaleYcrep)) { Y <- scaleYcrep(Y, mnyCrep, mxyCrep) }
	weights <- ifelse(rep(use.sqrt, nrow(summ)), sqrt(summ$n), summ$n)
	if(model == "loess") {
		fitCrep <- loess(mean ~ X * Y, span=span, degree=degree, weight=weights, data=summ)
		pcorCrep <- predict(fitCrep, newdata=data.frame(X, Y))
	} else if(model == "logit") {
		fitCrep <- glm(mean ~ X * Y, family=binomial("logit"), weight=weights, data=summ)
		pcorCrep <- predict(fitCrep, newdata=data.frame(X, Y), type="response")
	}
	pcorCrep <- pmin(pmax(pcorCrep, 0.0), 1.0)
	mapqCrep <- -10.0 * log10(1.0 - pcorCrep)
	XorigCrep <- trCrep$diffANorm
	YorigCrep <- trCrep$diffCANorm
	correctCrep <- trCrep$correctA == 1
	if(!is.null(scaleXcrep)) { XorigCrep <- scaleXcrep(XorigCrep, mnxCrep, mxxCrep) }
	if(!is.null(scaleYcrep)) { YorigCrep <- scaleYcrep(YorigCrep, mnyCrep, mxyCrep) }
	plotCrep <- function(Xex=NULL, Yex=NULL, alpha=NULL) {
		my.plot.fn <- plot.fn
		if(!is.null(plot.fn)) {
			my.plot.fn <- paste0(plot.fn, "pair.", stratum, ".crep.pdf")
		}
		fit2DPlot(
			fitCrep,
			model,
			XorigCrep,
			YorigCrep,
			correctCrep,
			Xex=Xex,
			Yex=Yex,
			alpha=alpha,
			xlab="Score of mate aln - best unchosen mate aln",
			ylab="Best concordant score - 2nd best",
			main="Paired-end case 1",
			plot.fn=my.plot.fn)
		if(!is.null(plot.fn)) {
			pdf(file=paste0(plot.fn, "pair.", x, ".trainroc.crep.pdf"), width=6, height=6)
		} else {
			quartz(width=6, height=6)
		}
		roc_plot(trCrep$correctA, mapqCrep, main="Paired case 1 (crep) training-data ROC")
		if(!is.null(plot.fn)) { dev.off() }
	}
	if(do.plot) { plotCrep() }
	
	#
	# Fit model for Case 2: unique condordant, repetitive mate
	#
	summ <- ddply(trCuni, .(diffANorm), summarise, mean=mean(correctA), n=length(correctA))
	# Gather Xs (best concordant scores)
	X <- summ$diffANorm
	mxxCuni <- max(X); mnxCuni <- min(X)
	if(!is.null(scaleXcuni)) { X <- scaleXcuni(X, mnxCuni, mxxCuni) }
	mxyCuni <- NULL; mnyCuni <- NULL
	weights <- ifelse(rep(use.sqrt, nrow(summ)), sqrt(summ$n), summ$n)
	if(model == "loess") {
		fitCuni <- loess(mean ~ X, span=span, degree=degree, weight=weights, data=summ)
		pcorCuni <- predict(fitCuni, newdata=data.frame(X))
	} else if(model == "logit") {
		fitCuni <- glm(mean ~ X, family=binomial("logit"), weight=weights, data=summ)
		pcorCuni <- predict(fitCuni, newdata=data.frame(X), type="response")
	}
	pcorCuni <- pmin(pmax(pcorCuni, 0.0), 1.0)
	mapqCuni <- -10.0 * log10(1.0 - pcorCuni)
	XorigCuni <- trCuni$diffANorm
	YorigCuni <- NULL
	correctCuni <- trCuni$correctA == 1
	if(!is.null(scaleXcuni)) { XorigCuni <- scaleXcuni(XorigCuni, mnxCuni, mxxCuni) }
	plotCuni <- function(Xex=NULL, alpha=NULL) {
		my.plot.fn <- plot.fn
		if(!is.null(plot.fn)) {
			my.plot.fn <- paste0(plot.fn, "pair.", stratum, ".fit.cuni.pdf")
		}
		fit1DPlot(
			fitCuni,
			model,
			XorigCuni,
			correctCuni,
			Xex=Xex,
			alpha=alpha,
			xlab="Score of mate aln - best unchosen mate aln",
			main="Paired-end case 2",
			plot.fn=my.plot.fn)
		if(!is.null(plot.fn)) {
			pdf(file=paste0(plot.fn, "pair.", x, ".trainroc.cuni.pdf"), width=6, height=6)
		} else {
			quartz(width=6, height=6)
		}
		roc_plot(trCuni$correctA, mapqCuni, main="Paired case 2 (cuni) training-data ROC")
		if(!is.null(plot.fn)) { dev.off() }
	}
	if(do.plot) { plotCuni() }
	
	#
	# Fit model for Case 3: unique
	#
	summ <- ddply(trUni, .(bestANorm), summarise, mean=mean(correctA), n=length(correctA))
	# Gather Xs (best concordant scores)
	X <- summ$bestANorm
	mxxUni <- max(X); mnxUni <- min(X)
	if(!is.null(scaleXuni)) { X <- scaleXuni(X, mnxUni, mxxUni) }
	weights <- ifelse(rep(use.sqrt, nrow(summ)), sqrt(summ$n), summ$n)
	if(model == "loess") {
		fitUni <- loess(mean ~ X, span=span, degree=degree, weight=weights, data=summ)
		pcorUni <- predict(fitUni, newdata=data.frame(X))
	} else if(model == "logit") {
		fitUni <- glm(mean ~ X, family=binomial("logit"), weight=weights, data=summ)
		pcorUni <- predict(fitUni, newdata=data.frame(X), type="response")
	}
	pcorUni <- pmin(pmax(pcorUni, 0.0), 1.0)
	mapqUni <- -10.0 * log10(1.0 - pcorUni)
	XorigUni <- trUni$bestANorm
	if(!is.null(scaleXuni)) { XorigUni <- scaleXuni(XorigUni, mnxUni, mxxUni) }
	correctUni <- trUni$correctA == 1
	plotUni <- function(Xex=NULL, alpha=NULL) {
		my.plot.fn <- plot.fn
		if(!is.null(plot.fn)) {
			my.plot.fn <- paste0(plot.fn, "pair.", stratum, ".fit.uni.pdf")
		}
		fit1DPlot(
			fitUni,
			model,
			XorigUni,
			correctUni,
			Xex=Xex,
			alpha=alpha,
			xlab="Best mate aln score",
			main="Paired-end case 3",
			plot.fn=my.plot.fn)
		if(!is.null(plot.fn)) {
			pdf(file=paste0(plot.fn, "pair.", x, ".trainroc.uni.pdf"), width=6, height=6)
		} else {
			quartz(width=6, height=6)
		}
		roc_plot(trUni$correctA, mapqUni, main="Paired case 3 (uni) training-data ROC")
		if(!is.null(plot.fn)) { dev.off() }
	}
	if(do.plot) { plotUni() }
	
	return(list(
		crep=fitCrep,
		crepBounds=c(mnxCrep, mxxCrep, mnyCrep, mxyCrep),
		crepX=XorigCrep,
		crepY=YorigCrep,
		crepCorrect=correctCrep,
		cuni=fitCuni,
		cuniBounds=c(mnxCuni, mxxCuni, mnyCuni, mxyCuni),
		cuniX=XorigCuni,
		cuniY=YorigCuni,
		cuniCorrect=correctCuni,
		uni=fitUni,
		uniBounds=c(mnxUni, mxxUni),
		uniX=XorigUni,
		uniCorrect=correctUni,
		model=model,
		plotCrep=plotCrep,
		plotCuni=plotCuni,
		plotUni=plotUni))
}

fitPairedStrata <- function(
	train,
	model="loess",
	span=2.0,
	scaleXcrep=NULL,
	scaleYcrep=NULL,
	scaleXcuni=NULL,
	scaleYcuni=NULL,
	scaleXuni=NULL,
	use.sqrt=F,
	do.plot=F,
	plot.fn=NULL)
{
	return(lapply(unique(train$stratum), function(x) {
		fitPaired(
			train[train$stratum == x,],
			stratum=x,
			model=model,
			span=span,
			scaleXcrep=scaleXcrep,
			scaleYcrep=scaleYcrep,
			scaleXcuni=scaleXcuni,
			scaleYcuni=scaleYcuni,
			scaleXuni=scaleXuni,
			use.sqrt=use.sqrt,
			do.plot=do.plot,
			plot.fn=plot.fn)
	}))
}

predictUnpaired <- function(
	al,
	fitU,
	scaleXrep=NULL,
	scaleYrep=NULL,
	scaleXuni=NULL)
{
	fitRep <- fitU$rep; fitUni <- fitU$uni
	repBounds <- fitU$repBounds
	uniBounds <- fitU$uniBounds
	
	#
	# Case 1: Repeat
	#
	alNextBest <- al$XS.i
	if("Xs.i" %in% colnames(al)) {
		alNextBest  <- pmax(alNextBest , al$Xs.i, na.rm=T)
	}
	al_rep_idx <- !is.na(alNextBest)
	al_rep <- al[al_rep_idx,]
	nextBest <- alNextBest[al_rep_idx]
	mx <- al_rep$Yn.i
	mn <- al_rep$YN.i
	
	best <- al_rep$AS.i
	diff <- best - nextBest
	X <- (best - mx) / (mx - mn)
	X <- pmin(pmax(X, repBounds[1]), repBounds[2])
	if(!is.null(scaleXrep)) { X <- scaleXrep(X, repBounds[1], repBounds[2]) }
	Y <- diff / (mx - mn)
	Y <- pmin(pmax(Y, repBounds[3]), repBounds[4])
	if(!is.null(scaleYrep)) { Y <- scaleYrep(Y, repBounds[3], repBounds[4]) }
	
	if(fitU$model == "loess") {
		pcorRep <- predict(fitRep, newdata=data.frame(X, Y))
	} else if(fitU$model == "logit") {
		pcorRep <- predict(fitRep, newdata=data.frame(X, Y), type="response")
	} else {
		print(paste0("Error: bad model: ", fitU$model))
	}
	if(any(is.na(pcorRep))) {
		print(pcorRep)
		stop("ERROR: NAs in unpaired pcorRep")
	}
	pcorRep <- pmin(1.0, pmax(0.0, pcorRep))
	mapqRep <- -10 * log10(1.0 - pcorRep)
	if("ZC.i" %in% colnames(al_rep)) {
		#alpha <- pmin(pmax(mapqRep - 5.0, 0.0) / 30.0, 1.0)
		badidx <- al_rep$ZC.i == 0 & mapqRep >= 5.0
		if(any(badidx)) { fitU$plotRep(X[badidx], Y[badidx]) }
	}
	
	#
	# Case 2: Unique
	#
	al_uni_idx <- is.na(alNextBest)
	al_uni <- al[al_uni_idx,]
	mx <- al_uni$Yn.i
	mn <- al_uni$YN.i
	
	best <- al_uni$AS.i
	X <- (best - mx) / (mx - mn)
	X <- pmin(pmax(X, uniBounds[1]), uniBounds[2])
	if(!is.null(scaleXuni)) { X <- scaleXuni(X, uniBounds[1], uniBounds[2]) }
	if(fitU$model == "loess") {
		pcorUni <- predict(fitUni, newdata=data.frame(X))
	} else if(fitU$model == "logit") {
		pcorUni <- predict(fitUni, newdata=data.frame(X), type="response")
	}
	if(any(is.na(pcorUni))) {
		print(pcorUni)
		stop("ERROR: NAs in unpaired pcorUni")
	}
	pcorUni <- pmin(1.0, pmax(0.0, pcorUni))
	mapqUni <- -10 * log10(1.0 - pcorUni)
	if("ZC.i" %in% colnames(al_uni)) {
		#alpha <- pmin(pmax(mapqUni - 5.0, 0.0) / 30.0, 1.0)
		#alpha <- pmin(pmax(mapqUni - 5.0, 0.0) / 30.0, 1.0)
		badidx <- al_uni$ZC.i == 0 & mapqUni >= 5.0
		if(any(badidx)) { fitU$plotUni(X[badidx]) }
	}
	
	# Compose result
	pcor <- rep(NA, nrow(al))
	pcor[al_uni_idx] <- pcorUni
	pcor[al_rep_idx] <- pcorRep
	if(any(is.na(pcor))) {
		print(pcor)
		stop("ERROR: Still some NAs")
	}
	mapq <- -10 * log10(1.0 - pcor)
	
	return(list(pcor=pcor, pcorRep=pcorRep, pcorUni=pcorUni,
	            mapq=mapq, mapqRep=mapqRep, mapqUni=mapqUni,
	            repIdx=al_rep_idx, uniIdx=al_uni_idx))
}

predictPaired <- function(
	al,
	fitP,
	scaleXcrep=NULL,
	scaleYcrep=NULL,
	scaleXcuni=NULL,
	scaleYcuni=NULL,
	scaleXuni=NULL,
	do.plot=F,
	plot.fn=NULL)
{
	fitCrep <- fitP$crep; fitCuni <- fitP$cuni; fitUni <- fitP$uni
	crepBounds <- fitP$crepBounds
	cuniBounds <- fitP$cuniBounds
	uniBounds <- fitP$uniBounds
	
	#
	# Case 1: Concordant repeat, mate repeat
	#
	al_crep_idx <- !is.na(al$Zp.i)
	al_crep <- al[al_crep_idx,]
	pcorCrep <- c()
	mapqCrep <- c()
	if(nrow(al_crep) > 0) {
		mxp <- al_crep$Yn.i + al_crep$Zn.i
		mnp <- al_crep$YN.i + al_crep$ZN.i
		mx <- al_crep$Yn.i
		mn <- al_crep$YN.i
		
		bestConc <- al_crep$ZP.i
		nextConc <- al_crep$Zp.i
		diffConc <- bestConc - nextConc
		
		best <- al_crep$AS.i
		bestUnchosen <- ifelse(is.na(al_crep$XS.i), mn, al_crep$XS.i)
		if("Xs.i" %in% colnames(al_crep)) {
			bestUnchosen <- pmax(bestUnchosen, al_crep$Xs.i, na.rm=T)
		}
		diff <- best - bestUnchosen
		
		X <- diff / (mx - mn)
		X <- pmin(pmax(X, crepBounds[1]), crepBounds[2])
		if(!is.null(scaleXcrep)) { X <- scaleXcrep(X, crepBounds[1], crepBounds[2]) }
		Y <- diffConc / (mxp - mnp)
		Y <- pmin(pmax(Y, crepBounds[3]), crepBounds[4])
		if(!is.null(scaleYcrep)) { Y <- scaleYcrep(Y, crepBounds[3], crepBounds[4]) }
		
		if(fitP$model == "loess") {
			pcorCrep <- predict(fitCrep, newdata=data.frame(X, Y))
		} else if(fitP$model == "logit") {
			pcorCrep <- predict(fitCrep, newdata=data.frame(X, Y), type="response")
		} else {
			print(paste0("Error: bad model: ", fitP$model))
		}
		pcorCrep <- pmin(1.0, pmax(0.0, pcorCrep))
		mapqCrep <- -10 * log10(1.0 - pcorCrep)
		if(any(is.na(pcorCrep))) {
			stop("ERROR: NAs in pcorCrep")
		}
		if("ZC.i" %in% colnames(al_crep)) {
			#alpha <- pmin(pmax(mapqCrep - 5.0, 0.0) / 30.0, 1.0)
			badidx <- al_crep$ZC.i == 0 & mapqCrep >= 5.0
			if(any(badidx)) { fitP$plotCrep(X[badidx], Y[badidx]) }
		}
	}

	#
	# Case 2: Mate repeat, no concordant repeat
	#
	bestUnchosen <- al$XS.i
	if("Xs.i" %in% colnames(al)) {
		bestUnchosen <- pmax(bestUnchosen, al$Xs.i, na.rm=T)
	}
	al_cuni_idx <- is.na(al$Zp.i) & !is.na(bestUnchosen)
	al_cuni <- al[al_cuni_idx,]
	pcorCuni <- c()
	mapqCuni <- c()
	if(nrow(al_cuni) > 0) {
		mxp <- al_cuni$Yn.i + al_cuni$Zn.i
		mnp <- al_cuni$YN.i + al_cuni$ZN.i
		mx <- al_cuni$Yn.i
		mn <- al_cuni$YN.i
		
		best <- al_cuni$AS.i
		bestUnchosenCuni <- al_cuni$XS.i
		if("Xs.i" %in% colnames(al_cuni)) {
			bestUnchosenCuni <- pmax(bestUnchosenCuni, al_cuni$Xs.i, na.rm=T)
		}
		diff <- best - bestUnchosenCuni

		bestConc <- al_cuni$ZP.i
		
		X <- diff / (mx - mn)
		X <- pmin(pmax(X, cuniBounds[1]), cuniBounds[2])
		if(!is.null(scaleXcuni)) { X <- scaleXcuni(X, cuniBounds[1], cuniBounds[2]) }
		if(fitP$model == "loess") {
			pcorCuni <- predict(fitCuni, newdata=data.frame(X))
		} else if(fitP$model == "logit") {
			pcorCuni <- predict(fitCuni, newdata=data.frame(X), type="response")
		}
		pcorCuni <- pmin(1.0, pmax(0.0, pcorCuni))
		mapqCuni <- -10 * log10(1.0 - pcorCuni)
		if(any(is.na(pcorCuni))) { stop("ERROR: NAs in pcorCuni") }
		if("ZC.i" %in% colnames(al_cuni)) {
			#alpha <- pmin(pmax(mapqCuni - 5.0, 0.0) / 30.0, 1.0)
			badidx <- al_cuni$ZC.i == 0 & mapqCuni >= 5.0
			if(any(badidx)) { fitP$plotCuni(X[badidx]) }
		}
	}
	
	#
	# Case 3: No mate repeat, no concordant repeat
	#
	al_uni_idx <- is.na(al$Zp.i) & is.na(bestUnchosen)
	al_uni <- al[al_uni_idx,]
	pcorUni <- c()
	mapqUni <- c()
	if(nrow(al_uni) > 0) {
		mx <- al_uni$Yn.i
		mn <- al_uni$YN.i
		best <- al_uni$AS.i
		X <- (best - mx) / (mx - mn)
		X <- pmin(pmax(X, uniBounds[1]), uniBounds[2])
		if(!is.null(scaleXuni)) { X <- scaleXuni(X, uniBounds[1], uniBounds[2]) }
		if(fitP$model == "loess") {
			pcorUni <- predict(fitUni, newdata=data.frame(X))
		} else if(fitP$model == "logit") {
			pcorUni <- predict(fitUni, newdata=data.frame(X), type="response")
		}
		pcorUni <- pmin(1.0, pmax(0.0, pcorUni))
		mapqUni <- -10 * log10(1.0 - pcorUni)
		if(any(is.na(pcorUni))) {
			stop("ERROR: NAs in pcorUni")
		}
		if("ZC.i" %in% colnames(al_uni)) {
			#alpha <- pmin(pmax(mapqUni - 5.0, 0.0) / 30.0, 1.0)
			badidx <- al_uni$ZC.i == 0 & mapqUni >= 5.0
			if(any(badidx)) { fitP$plotUni(X[badidx]) }
		}
	}
	
	# Compose result
	pcor <- rep(NA, nrow(al))
	pcor[al_crep_idx] <- pcorCrep
	pcor[al_cuni_idx] <- pcorCuni
	pcor[al_uni_idx] <- pcorUni
	if(any(is.na(pcor))) {
		print(pcor)
		stop("ERROR: Still some NAs")
	}
	mapq <- -10 * log10(1.0 - pcor)
	return(list(
		pcor=pcor,
		pcorCrep=pcorCrep,
		pcorCuni=pcorCuni,
		pcorUni=pcorUni,
		mapq=mapq,
		mapqCrep=mapqCrep,
		mapqCuni=mapqCuni,
		mapqUni=mapqUni,
		crepIdx=al_crep_idx,
		cuniIdx=al_cuni_idx,
		uniIdx=al_uni_idx))
}

predictUnpairedStrata <- function(
	al,
	fit,
	scaleXrep=NULL,
	scaleYrep=NULL,
	scaleXuni=NULL,
	do.plot=F,
	plot.fn=NULL)
{
	al$predStratum <- rep(Inf, nrow(al))
	al$predPcor <- rep(Inf, nrow(al))
	al$predMapq <- rep(Inf, nrow(al))
	for(x in unique(al$stratum)) {
		wh <- which(al$stratum == x)
		alst <- al[wh,]
		pred <- predictUnpaired(
			alst,
			fit[[x]],
			scaleXrep=scaleXrep,
			scaleYrep=scaleYrep,
			scaleXuni=scaleXuni)
		al$predStratum[wh] <- x
		al$predPcor[wh] <- pred$pcor
		al$predMapq[wh] <- pred$mapq
		if(do.plot) {
			if(!is.null(plot.fn)) {
				pdf(file=paste0(plot.fn, "unp.st", x, ".roc.comb.pdf"), width=6, height=6)
			} else {
				quartz(width=6, height=6)
			}
			roc_2plot(
				alst$ZC.i,
				alst$mapq,
				pred$mapq,
				mapq1_lab="Bowtie 2", mapq2_lab="SimQ",
				main=paste("Unpaired alignments (all), stratum", x))
			if(!is.null(plot.fn)) {
				dev.off()
				pdf(file=paste0(plot.fn, "unp.st", x, ".roc.rep.pdf"), width=6, height=6)
			} else {
				quartz(width=6, height=6)
			}
			roc_2plot(
				alst[pred$repIdx,]$ZC.i,
				alst$mapq[pred$repIdx],
				pred$mapq[pred$repIdx],
				mapq1_lab="Bowtie 2", mapq2_lab="SimQ",
				main=paste("Unpaired alignments case 1 (rep), stratum", x))
			if(!is.null(plot.fn)) {
				dev.off()
				pdf(file=paste0(plot.fn, "unp.st", x, ".roc.uni.pdf"), width=6, height=6)
			} else {
				quartz(width=6, height=6)
			}
			roc_2plot(
				alst[pred$uniIdx,]$ZC.i,
				alst$mapq[pred$uniIdx],
				pred$mapq[pred$uniIdx],
				mapq1_lab="Bowtie 2", mapq2_lab="SimQ",
				main=paste("Unpaired alignments case 2 (uni), stratum", x))
			if(!is.null(plot.fn)) { dev.off() }
		}
	}
	if(any(!is.finite(al$predStratum)) ||
	   any(!is.finite(al$predPcor)))
	{
		print("At least one non-finite in loess predictions!!!")
	}
	return(al)
}

predictPairedStrata <- function(
	al,
	fit,
	scaleXcrep=NULL,
	scaleYcrep=NULL,
	scaleXcuni=NULL,
	scaleYcuni=NULL,
	scaleXuni=NULL,
	do.plot=F,
	plot.fn=NULL)
{
	al$predStratum <- rep(Inf, nrow(al))
	al$predPcor <- rep(Inf, nrow(al))
	al$predMapq <- rep(Inf, nrow(al))
	for(x in unique(al$stratum)) {
		wh <- which(al$stratum == x)
		alst <- al[wh,]
		pred <- predictPaired(
			alst,
			fit[[x]],
			scaleXcrep=scaleXcrep,
			scaleYcrep=scaleYcrep,
			scaleXcuni=scaleXcuni,
			scaleYcuni=scaleYcuni,
			scaleXuni=scaleXuni)
		al$predStratum[wh] <- x
		al$predPcor[wh] <- pred$pcor
		al$predMapq[wh] <- pred$mapq
		if(do.plot) {
			if(!is.null(plot.fn)) {
				pdf(file=paste0(plot.fn, "conc.st", x, ".roc.comb.pdf"), width=6, height=6)
			} else {
				quartz(width=6, height=6)
			}
			roc_2plot(
				alst$ZC.i,
				alst$mapq,
				pred$mapq,
				mapq1_lab="Bowtie 2", mapq2_lab="SimQ",
				main="Concordant alignments (all)")
			if(!is.null(plot.fn)) {
				dev.off()
				pdf(file=paste0(plot.fn, "conc.st", x, ".roc.crep.pdf"), width=6, height=6)
			} else {
				quartz(width=6, height=6)
			}
			roc_2plot(
				alst[pred$crepIdx,]$ZC.i,
				alst$mapq[pred$crepIdx],
				pred$mapq[pred$crepIdx],
				mapq1_lab="Bowtie 2", mapq2_lab="SimQ",
				main="Concordant alignments case 1 (rep, rep)")
			if(!is.null(plot.fn)) {
				dev.off()
				pdf(file=paste0(plot.fn, "conc.st", x, ".roc.cuni.pdf"), width=6, height=6)
			} else {
				quartz(width=6, height=6)
			}
			roc_2plot(
				alst[pred$cuniIdx,]$ZC.i,
				alst$mapq[pred$cuniIdx],
				pred$mapq[pred$cuniIdx],
				mapq1_lab="Bowtie 2", mapq2_lab="SimQ",
				main="Concordant alignments case 2 (uni, rep)")
			if(!is.null(plot.fn)) {
				dev.off()
				pdf(file=paste0(plot.fn, "conc.st", x, ".roc.uni.pdf"), width=6, height=6)
			} else {
				quartz(width=6, height=6)
			}
			roc_2plot(
				alst[pred$uniIdx,]$ZC.i,
				alst$mapq[pred$uniIdx],
				pred$mapq[pred$uniIdx],
				mapq1_lab="Bowtie 2", mapq2_lab="SimQ",
				main="Concordant alignments case 3 (uni, uni)")
			if(!is.null(plot.fn)) {
				dev.off()
			}
		}
	}
	if(any(!is.finite(al$predStratum)) ||
	   any(!is.finite(al$predPcor)))
	{
		print("At least one non-finite in loess predictions!!!")
	}
	return(al)
}

loadResults <- function(
	len,
	sensitivity,
	paired=F,
	model="logit",
	span=1.0,
	do.plot=F,
	plot.fn=NULL,
	use.sqrt=F)
{
	pairedstr = ifelse(paired, "12", "0")
	fn <- paste0("r", pairedstr, "_ill_mason_", len, "_500k.bt2_", sensitivity)
	
	# Read in alignments and training data
	print(paste("  Opening training and SAM", format(Sys.time(), "%X"), "..."))
	data <- openTrainingAndSam(fn)
	
	# Let all of the dimensions representing differences be log-scaled
	scaleXrep_u  <- function(x, mnx, mxx) { log2(-x - mnx + 0.1) } # best
	scaleYrep_u  <- function(y, mny, mxy) { log2( y - mny + 0.1) } # diff
	scaleXuni_u  <- NULL                                     # best
	scaleXcrep_p <- function(x, mnx, mxx) { pmax(x, -0.25) } # diff concordant
	scaleYcrep_p <- NULL # diff mate
	scaleXcuni_p <- function(x, mnx, mxx) { log2( x - mnx + 0.1) } # diff mate
	#scaleXcuni_p <- NULL # diff mate
	scaleYcuni_p <- function(y, mny, mxy) { 1 } # best concordant
	scaleXuni_p  <- NULL                                     # best mate
	
	ret <- list()
	
	# Convert training data to vectors suitable for handing to loess
	if(!is.null(data$trainU)) {
		print(paste("  Stratifying training data (unpaired)", format(Sys.time(), "%X"), "..."))
		strataU <- stratify(data$trainU, "len")
		print(paste("  Fitting (unpaired)", format(Sys.time(), "%X"), "..."))
		fitU <- fitUnpairedStrata(
			strataU$data,
			model=model,
			span=span,
			scaleXrep=scaleXrep_u,
			scaleYrep=scaleYrep_u,
			scaleXuni=scaleXuni_u,
			do.plot=do.plot,
			plot.fn=plot.fn,
			use.sqrt=use.sqrt)
		print(paste("  Stratifying real data (unpaired)", format(Sys.time(), "%X"), "..."))
		resU <- stratify(data$alU, "len", strata=strataU$strata)
		alU <- resU$data; alU$stratum <- pmax(alU$stratum, 1)
		print(paste("  Predicting (unpaired)", format(Sys.time(), "%X"), "..."))
		predU <- predictUnpairedStrata(
			alU,
			fitU,
			scaleXrep_u,
			scaleYrep_u,
			scaleXuni_u,
			do.plot=do.plot,
			plot.fn=plot.fn)
		predU$predMapq <- pmin(predU$predMapq, 99.0)
		ret$predU <- predU
		ret$mapqU <- predU$predMapq
		ret$errU.old <- rankingError(predU, predU$mapq)
		ret$errU <- rankingError(predU, predU$predMapq)
	}
	if(!is.null(data$trainM)) {
		print(paste("  Stratifying training data (paired)", format(Sys.time(), "%X"), "..."))
		strataM <- stratify(data$trainM, "len")
		strataC <- stratify(data$trainC, "lenA")
		print(paste("  Fitting (paired)", format(Sys.time(), "%X"), "..."))
		fitM <- fitUnpairedStrata(
			strataM$data,
			model=model,
			span=span,
			scaleXrep=scaleXrep_u,
			scaleYrep=scaleYrep_u,
			scaleXuni=scaleXuni_u,
			do.plot=do.plot,
			plot.fn=plot.fn,
			use.sqrt=use.sqrt)
		fitC <- fitPairedStrata(
			strataC$data,
			model=model,
			span=span,
			scaleXcrep=scaleXcrep_p,
			scaleYcrep=scaleYcrep_p,
			scaleXcuni=scaleXcuni_p,
			scaleYcuni=scaleYcuni_p,
			scaleXuni=scaleXuni_p,
			do.plot=do.plot,
			plot.fn=plot.fn,
			use.sqrt=use.sqrt)
		print(paste("  Stratifying real data (paired)", format(Sys.time(), "%X"), "..."))
		resM <- stratify(data$alM, "len", strata=strataM$strata)
		resC <- stratify(data$alC, "len", strata=strataC$strata)
		alM <- resM$data; alM$stratum <- pmax(alM$stratum, 1)
		alC <- resC$data; alC$stratum <- pmax(alC$stratum, 1)
		print(paste("  Predicting (paired)", format(Sys.time(), "%X"), "..."))
		predM <- predictUnpairedStrata(
			alM,
			fitM,
			scaleXrep_u,
			scaleYrep_u,
			scaleXuni_u,
			do.plot=do.plot,
			plot.fn=plot.fn)
		predM$predMapq <- pmin(predM$predMapq, 99.0)
		predC <- predictPairedStrata(
			alC,
			fitC,
			scaleXcrep=scaleXcrep_p,
			scaleYcrep=scaleYcrep_p,
			scaleXcuni=scaleXcuni_p,
			scaleYcuni=scaleYcuni_p,
			scaleXuni=scaleXuni_p,
			do.plot=do.plot,
			plot.fn=plot.fn)
		predC$predMapq <- pmin(predC$predMapq, 99.0)
		ret$predM <- predM
		ret$mapqM <- predM$predMapq
		ret$errM.old <- rankingError(predM, predM$mapq)
		ret$errM <- rankingError(predM, predM$predMapq)
		
		ret$predC <- predC
		ret$mapqC <- predC$predMapq
		ret$errC.old <- rankingError(predC, predC$mapq)
		ret$errC <- rankingError(predC, predC$predMapq)
	}
	return(ret)
}

if(FALSE) {
	setwd("~/Documents/workspace/mapq/examples/mason/")
	for(paired in c(F, T)) {
		for(len in seq(50, 300, 25)) {
			for(sens in c("vf", "f", "s", "vs")) {
				id <- paste(ifelse(paired, "pair", "unp"), len, sens, sep=".")
				print(id)
				dir.create(id)
				loadResults(len, sens, paired, model="loess", span=1.0, do.plot=T, plot.fn=paste0(id, "/"), use.sqrt=F)
			}
		}
	}
}
