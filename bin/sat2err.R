#!/usr/bin/env Rscript

library(bitops)

# Open a tabulated sam file
openTabulatedSam <- function(fn) {
	return(read.table(fn, comment.char="", quote="", header=T, stringsAsFactors=F))
}

# Given data frame of alignments, return data frame with just aligned reads
selectAligned <- function(x) {
	return(bitAnd(x$flag, 4) == 0)
}

# Penalty is sum of ranks of all incorrectly aligned reads.  Rank is lower for
# reads with lower MAPQ.
rankingError <- function(x, mapq) {
	ordr <- order(mapq)
	return(sum(as.numeric(which(x$ZC.i[ordr] == 0))))
}

# Take table of alignments and parallel vector of mapping qualities and
# return a data frame containing columns for MAPQs (sorted in
# descending order), cumulative number of correct alignments, and
# cumulative number of incorrect alignments.
roc_table <- function(x, mapq) {
	xord <- order(mapq)
	tab <- table(mapq[xord], x$ZC.i[xord])
	tab <- data.frame(mapq=as.numeric(rownames(tab)), cor=tab[,2], incor=tab[,1])
	tab <- tab[order(tab$mapq, decreasing=T),]
	return(tab)
}

# Plot a roc-like curve (except using absolute numbers of true
# positives and false negatives) for the given alignments (x) and
# mapping-quality ranking (mapq).
roc_plot <- function(x, mapq, color="blue", first=T, xlim=NULL, ylim=NULL, main="MAPQ ROCs", xlab="Incorrect alignments", ylab="Correct alignments") {
	tab <- roc_table(x, mapq)
	if(first) {
		plot(cumsum(tab$incor), cumsum(tab$cor), col="blue", type="o", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)
	} else {
		lines(cumsum(tab$incor), cumsum(tab$cor), col="red", type="o")
	}
}

# Plot two roc-like curves (except using absolute numbers of true
# positives and false negatives) for the given alignments (x) and the
# given two mapping-quality rankings (mapq1 and mapq2).
roc_2plot <- function(x, mapq1, mapq2, colors=c("blue", "red"), sca=1, round_digits=2, mapq1_lab="MAPQ 1", mapq2_lab="MAPQ 2", main="MAPQ ROCs") {
	if(round_digits > 0) {
		mapq1 <- round(mapq1, digits=round_digits)
		mapq2 <- round(mapq2, digits=round_digits)
	}
	xroc <- roc_table(x, mapq1)
	xroc2 <- roc_table(x, mapq2)
	mx_x <- max(max(cumsum(xroc$incor)), max(cumsum(xroc2$incor)))
	mn_x <- min(min(cumsum(xroc$incor)), min(cumsum(xroc2$incor)))
	mx_y <- max(max(cumsum(xroc$cor)), max(cumsum(xroc2$cor)))
	mn_y <- min(min(cumsum(xroc$cor)), min(cumsum(xroc2$cor)))
	roc_plot(x, mapq1,  color=colors[0], first=T, xlim=c(mn_x, mn_x + (mx_x - mn_x) * sca), ylim=c(mn_y + (mx_y - mn_y) * (1-sca), mx_y), main=main)
	roc_plot(x, mapq2, color=colors[1], first=F, xlim=c(mn_x, mn_x + (mx_x - mn_x) * sca), ylim=c(mn_y + (mx_y - mn_y) * (1-sca), mx_y))
	legend("bottomright", c(mapq1_lab, mapq2_lab), col=colors, pch=c(1, 1), lty=c(1, 1))
}

compareToEmpiricalMap <- function(x) {
	
}

args <- commandArgs()
fn <- args[length(args)]
tab <- openTabulatedSam(fn)
tab <- tab[selectAligned(tab),]

rocFn <- sub(".sat$", ".old.roc", fn)
rocNewFn <- sub(".sat$", ".new.roc", fn)
rocPdfFn <- sub(".sat$", ".roc.pdf", fn)

# Calculate ranking errors
errOrig    <- rankingError(tab, tab$mapq)
errTs      <- rankingError(tab, tab$XQ.f)
errTsRound <- rankingError(tab, round(tab$XQ.f))

errDiff <- errOrig - errTsRound

# Make a ROC table and corresponding ROC plot comparing Bowtie 2's native MAPQ
# and the MAPQ predicted in our tandem simulation scheme.
pdf(file=rocPdfFn)
roc.old <- roc_table(tab, tab$mapq)
roc.new <- roc_table(tab, round(tab$XQ.f))
write.table(roc.old, rocFn, quote=F, sep="\t", row.names=F)
write.table(roc.new, rocNewFn, quote=F, sep="\t", row.names=F)
roc_2plot(tab, round(tab$XQ.f), tab$mapq, mapq1_lab="Tandem simulation (rounded)", mapq2_lab="Bowtie 2", main=paste(fn, ", ranking error (old-new) =", errDiff, round_digits=0))
dev.off()

# Print ranking error info
if(errTsRound > errOrig) {
  print("***********************")
}
print(paste("mapq", errOrig))
print(paste("ts", errTs))
print(paste("ts rounded", errTsRound))
print(paste("# distinct mapqs", length( unique( round(tab$XQ.f) ) )   ))
if(errTsRound > errOrig) {
	print("***********************")
}
