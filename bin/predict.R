#!/usr/bin/env Rscript

library(plyr)

# Linear regression with best * diff
fit.unpaired.linear.times_log <- function(tab) {
	summ <- ddply(tab, .(bestNorm, diffNorm), summarise, mean=mean(correct), n=length(correct))
	fit <- glm(mean ~ bestNorm * I(log(diffNorm - min(diffNorm) + 0.1)), data=summ, family=gaussian(), weight=n)
	return (function(tab) { predict(fit, tab, type="response") })
}

fit.paired.linear.times_log <- function(tab) {
	summ <- ddply(tab, .(bestNorm1, diffNorm1), summarise, mean=mean(correct), n=length(correct))
	fit <- glm(mean ~ bestNorm1 * I(log(diffNorm1 - min(diffNorm1) + 0.1)), data=summ, family=gaussian(), weight=n)
	return (function(tab) { predict(fit, tab, type="response") })
}

# Logistic regression with best * diff
fit.unpaired.logit.times <- function(tab) {
	summ <- ddply(tab, .(bestNorm, diffNorm), summarise, mean=mean(correct), n=length(correct))
	fit <- glm(mean ~ bestNorm * diffNorm, data=summ, family=binomial("logit"), weight=n)
	return (function(tab) { predict(fit, tab, type="response") })
}

fit.paired.logit.times <- function(tab) {
	summ <- ddply(tab, .(bestNorm1, diffNorm1), summarise, mean=mean(correct), n=length(correct))
	fit <- glm(mean ~ bestNorm1 * diffNorm1, data=summ, family=binomial("logit"), weight=n)
	return (function(tab) { predict(fit, tab, type="response") })
}

# Logistic regression after log-transforming diff
fit.unpaired.logit.times_log <- function(tab) {
	summ <- ddply(tab, .(bestNorm, diffNorm), summarise, mean=mean(correct), n=length(correct))
	fit <- glm(mean ~ bestNorm * I(log(diffNorm - min(diffNorm) + 0.1)), data=summ, family=binomial("logit"), weight=n)
	return (function(tab) { predict(fit, tab, type="response") })
}

fit.paired.logit.times_log <- function(tab) {
	summ <- ddply(tab, .(bestNorm1, diffNorm1), summarise, mean=mean(correct), n=length(correct))
	fit <- glm(mean ~ bestNorm1 * I(log(diffNorm1 - min(diffNorm1) + 0.1)), data=summ, family=binomial("logit"), weight=n)
	return (function(tab) { predict(fit, tab, type="response") })
}

# Logistic regression with best + diff
fit.unpaired.logit.plus <- function(tab) {
	summ <- ddply(tab, .(bestNorm, diffNorm), summarise, mean=mean(correct), n=length(correct))
	fit <- glm(mean ~ bestNorm + diffNorm, data=summ, family=binomial("logit"), weight=n)
	return (function(tab) { predict(fit, tab, type="response") })
}

fit.paired.logit.plus <- function(tab) {
	summ <- ddply(tab, .(bestNorm1, diffNorm1), summarise, mean=mean(correct), n=length(correct))
	fit <- glm(mean ~ bestNorm1 + diffNorm1, data=summ, family=binomial("logit"), weight=n)
	return (function(tab) { predict(fit, tab, type="response") })
}

# Random forest with best + diff
fit.unpaired.random_forest <- function(tab) {
	library(randomForest)
	summ <- ddply(tab, .(bestNorm, diffNorm), summarise, mean=mean(correct), n=length(correct))
	fit <- randomForest(mean ~ bestNorm + diffNorm, data=summ, weight=n)
	return (function(tab) { predict(fit, tab, type="response") })
}

fit.paired.random_forest <- function(tab) {
	library(randomForest)
	#summ <- ddply(tab, .(bestNorm1, diffNorm1), summarise, mean=mean(correct), n=length(correct))
	#fit <- randomForest(mean ~ bestNorm1 + diffNorm1, data=summ, weight=n)
	#fit <- randomForest(correct ~ bestNorm1 + diffNorm1 + diffConc, data=tab)
	summ <- ddply(tab, .(bestNorm1, diffNorm1, diffConc), summarise, mean=mean(correct), n=length(correct))
	fit <- randomForest(mean ~ bestNorm1 + diffNorm1 + diffConc, data=summ, weight=n)
	return (function(tab) { predict(fit, tab, type="response") })
}

# Read a training or test CSV file.  The arguments are so that we can
# ensure they're both scaled the same way.
read.ts.unpaired <- function(fn) {
	tab <- read.csv(fn)
	minbest.default <- min(tab$best)
	maxbest.default <- max(tab$best)
	minbest <- ifelse(is.na(tab$minv), minbest.default, tab$minv)
	maxbest <- ifelse(is.na(tab$maxv), maxbest.default, tab$maxv)
	tab$bestNorm <- (tab$best - minbest) / (maxbest - minbest)
	tab$secbestNorm <- (tab$secbest  - minbest) / (maxbest - minbest)
	minsecbest.norm <- min(tab$secbestNorm, na.rm=T)
	tab$secbestNorm <- ifelse(is.na(tab$secbestNorm), minsecbest.norm, tab$secbestNorm)
	tab$diffNorm <- tab$bestNorm - tab$secbestNorm
	return(tab)
}

# Read a training or test CSV file.  The arguments are so that we can
# ensure they're both scaled the same way.
read.ts.paired <- function(fn) {
	tab <- read.csv(fn)
	minbest.default1 <- min(tab$best1)
	maxbest.default1 <- max(tab$best1)
	minbest.default2 <- min(tab$best2)
	maxbest.default2 <- max(tab$best2)
	minbest1 <- ifelse(is.na(tab$minv1), minbest.default1, tab$minv1)
	maxbest1 <- ifelse(is.na(tab$maxv1), maxbest.default1, tab$maxv1)
	minbest2 <- ifelse(is.na(tab$minv2), minbest.default2, tab$minv2)
	maxbest2 <- ifelse(is.na(tab$maxv2), maxbest.default2, tab$maxv2)
	minconc <- minbest1 + minbest2
	maxconc <- maxbest1 + maxbest2
	
	# Mate 1s
	tab$bestNorm1 <- (tab$best1 - minbest1) / (maxbest1 - minbest1)
	tab$secbestNorm1 <- (tab$secbest1 - minbest1) / (maxbest1 - minbest1)
	minsecbest.norm1 <- min(tab$secbestNorm1, na.rm=T)
	tab$secbestNorm1 <- ifelse(is.na(tab$secbestNorm1), minsecbest.norm1, tab$secbestNorm1)
	
	# Mate 2s
	tab$bestNorm2 <- (tab$best2 - minbest2) / (maxbest2 - minbest2)
	tab$secbestNorm2 <- (tab$secbest2 - minbest2) / (maxbest2 - minbest2)
	minsecbest.norm2 <- min(tab$secbestNorm2, na.rm=T)
	tab$secbestNorm2 <- ifelse(is.na(tab$secbestNorm2), minsecbest.norm2, tab$secbestNorm2)
	
	# Concordant pairs
	tab$bestConc <- (tab$bestconc - minconc) / (maxconc - minconc)
	tab$secbestConc <- (tab$secbestconc - minconc) / (maxconc - minconc)
	minsecbest.conc <- min(tab$secbestConc, na.rm=T)
	tab$secbestConc <- ifelse(is.na(tab$secbestConc), minsecbest.conc, tab$secbestConc)
	
	tab$diffNorm1 <- tab$bestNorm1 - tab$secbestNorm1
	tab$diffNorm2 <- tab$bestNorm2 - tab$secbestNorm2
	tab$diffConc <- tab$bestConc - tab$secbestConc
	
	return(tab)
}

plotDropRateMLog <- function(pcor, cor, pcor2=NULL, pdf.fn=NULL) {
	if(!is.null(pdf.fn)) {
		pdf(file=pdf.fn)
	}
	ordr <- order(jitterPcor(pcor), decreasing=T)
	x <- -log2(1.0 - seq(1, length(cor)) / length(cor))
	maxx <- ceiling(log2(length(cor)))
	y <- cumsum(!cor[ordr])
	plot(x, y, xlab="Drop rate", xaxt="n", ylab="Cumulative # incorrect alignments", type="l")
	if(!is.null(pcor2)) {
		ordr2 <- order(jitterPcor(pcor2), decreasing=T)
		y2 <- cumsum(!cor[ordr2])
		lines(x, y2, type="l", col="red")
	}
	axis(1, at=seq(1, maxx), labels=sprintf("%0.5f", (2.0 ^ seq(-1, -maxx))), las=2, cex.axis=0.7)
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

plotSmoothedComparison <- function(pcor, cor, pdf.fn=NULL, spar=0.5, nbootstrap=1) {
	if(!is.null(pdf.fn)) {
		pdf(file=pdf.fn)
	}
	ordr <- order(jitterPcor(pcor), decreasing=T)
	xfracs <- seq(1, length(cor)+1) / (length(cor)+1)
	xfracs <- xfracs[1:length(xfracs)-1]
	maxx <- ceiling(log2(length(cor)))
	x <- -log2(1.0 - xfracs)
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
	xlab <- append("1.00000", sprintf("%0.5f", (2.0 ^ seq(-1, -maxx))))
	axis(1, at=seq(0, maxx), labels=xlab, las=2, cex.axis=0.7)
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

plotSmoothedDifference <- function(pcor, cor, pdf.fn=NULL, spar=0.5, nbootstrap=1) {
	if(!is.null(pdf.fn)) {
		pdf(file=pdf.fn)
	}
	ordr <- order(jitterPcor(pcor), decreasing=T)
	xfracs <- seq(1, length(cor)) / length(cor)
	maxx <- ceiling(log2(length(cor)))
	x <- -log2(1.0 - xfracs)
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
	xlab <- append("1.00000", sprintf("%0.5f", (2.0 ^ seq(-1, -maxx))))
	axis(1, at=seq(0, maxx), labels=xlab, las=2, cex.axis=0.7)
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

predictors <- list(
	#'svm'               = list(fit.unpaired=fit.unpaired.svm, fit.paired=fit.paired.svm)
	'random_forest'    = list(fit.unpaired=fit.unpaired.random_forest, fit.paired=fit.paired.random_forest)
	#,'logit_times_log'  = list(fit.unpaired=fit.unpaired.logit.times_log, fit.paired=fit.paired.logit.times_log)
)

inputs <- list(
	#'unp' = list(train='training_unp.csv.gz', test='test_unp.csv.gz', paired=F, out='test_unp_mapq')
	#,'mates' = list(train='training_mates.csv.gz', test='test_mates.csv.gz', read=read.ts.unpaired, fit=fit.unpaired, out='test_mates_mapq.csv')
	'conc' = list(train='training_conc.csv.gz', test='test_conc.csv.gz', paired=T, out='test_conc_mapq.csv')
)

#fractions <- c(0.5, 0.4, 0.3, 0.25, 0.225, 0.2, 0.175, 0.15, 0.125, 0.1, 0.075, 0.05, 0.025, 0.02, 0.015, 0.01)
fractions <- c(0.5, 0.25, 0.2, 0.15, 0.1)

predictMapq <- function(inp, predict.training=T) {
	for(inpname in names(inputs)) {
		inp <- inputs[[inpname]]
		if(inp$paired) {
			train <- read.ts.paired(inp$train)
			test <- read.ts.paired(inp$test)
		} else {
			train <- read.ts.unpaired(inp$train)
			test <- read.ts.unpaired(inp$test)
		}
		sdiff.summ <- list()
		rankerr.summ <- list()
		for(predname in names(predictors)) {
			print(paste("Predictor:", predname, ", input:", inpname))
			predictor <- predictors[[predname]]
			sdiff.summ[[paste(predname, inpname)]] <- c()
			for(fraction in fractions) {
				print(paste("  Fraction:", fraction))
				sm <- sample(nrow(train), size=round(fraction*nrow(train)))
				if(inp$paired) {
					pcor <- predictor$fit.paired(train[sm,])(test)
					orig.mapq <- test$mapq1
				} else {
					pcor <- predictor$fit.unpaired(train[sm,])(test)
					orig.mapq <- test$mapq
				}
				mapq = -10.0 * log10(1.0 - pcor)
				cor = test$correct
				sdiff.test <- sdiff(pcor, cor)
				rankerr.test <- rankerr(pcor, cor)
				sdiff.summ[[paste(predname, inpname)]] <- append(sdiff.summ[[paste(predname, inpname)]], sdiff.test)
				rankerr.summ[[paste(predname, inpname)]] <- append(rankerr.summ[[paste(predname, inpname)]], rankerr.test)
				if(fraction == fractions[1]) {
					write.csv(data.frame(pcor=pcor, mapq=mapq, orig=orig.mapq, correct=cor), file=paste0(inp$out, "_test_", predname, ".csv"), quote=F, row.names=F)
					print(paste('test sdiff=', format(sdiff.test, scientific=T)))
					print(paste('test rankerr=', format(rankerr.test, scientific=T)))
					if(all(!is.na(test$correct))) {
						plotSmoothedDifference(pcor, test$correct, paste0(inp$out, "_test_", predname, "_smoothdiff.pdf"))
						plotSmoothedComparison(pcor, test$correct, paste0(inp$out, "_test_", predname, "_smoothcmp.pdf"))
						plotDropRateMLog(pcor, test$correct, orig.mapq, paste0(inp$out, "_test_", predname, "_drop.pdf"))
					}
				}
				if(predict.training) {
					if(inp$paired) {
						pcor <- predictor$fit.paired(train[sm,])(train)
						orig.mapq <- train$mapq1
					} else {
						pcor <- predictor$fit.unpaired(train[sm,])(train)
						orig.mapq <- train$mapq
					}
					mapq = -10.0 * log10(1.0 - pcor)
					cor <- train$correct
					if(fraction == fractions[1]) {
						print(paste('training sdiff=', format(sdiff(pcor, cor), scientific=T)))
						print(paste('training rankerr=', format(rankerr(pcor, cor), scientific=T)))
						write.csv(data.frame(pcor=pcor, mapq=mapq, orig=orig.mapq, correct=cor), file=paste0(inp$out, "_train_", predname, ".csv"), quote=F, row.names=F)
						plotSmoothedDifference(pcor, train$correct, paste0(inp$out, "_train_", predname, "_smoothdiff.pdf"))
						plotSmoothedComparison(pcor, train$correct, paste0(inp$out, "_test_", predname, "_smoothcmp.pdf"))
						plotDropRateMLog(pcor, train$correct, orig.mapq, paste0(inp$out, "_train_", predname, "_drop.pdf"))
					}
				}
			}
		}
		write.csv(data.frame(fraction=fractions, sdiff.summ), file=paste0(inp$out, "_sdiff.csv"), row.names=F)
		write.csv(data.frame(fraction=fractions, rankerr.summ), file=paste0(inp$out, "_rankerr.csv"), row.names=F)
	}
	return(NULL)
}

predictMapq()
