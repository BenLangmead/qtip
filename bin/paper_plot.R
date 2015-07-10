plot_cums <- function(ncol, sse_cum, sse_incor1, sse_incor2, incor_cum, incor_incor1, incor_incor2, mapq, mapq_orig)
{
    options(scipen=10)

    plot_mapq <- function(cum, mapq, y=0) {
        m40 <- append(cum[mapq > 40.0], max(cum))
        m30 <- append(cum[mapq <= 40.0 & mapq > 30.0], max(cum))
        m20 <- append(cum[mapq <= 30.0 & mapq > 20.0], max(cum))
        m10 <- append(cum[mapq <= 20.0 & mapq > 10.0], max(cum))
        m0 <- append(cum[mapq <= 10.0], max(cum))
        yheight <- 10.0
        xl <- 0
        xr <- max(cum[mapq > 40.0], 0)
        rect(xl, y, xr, y + yheight, col='blue', border=NA)
        xl <- xr
        xr <- max(cum[mapq > 30.0], 0)
        rect(xl, y, xr, y + yheight, col='dodgerblue1', border=NA)
        xl <- xr
        xr <- max(cum[mapq > 20.0], 0)
        rect(xl, y, xr, y + yheight, col='green', border=NA)
        xl <- xr
        xr <- max(cum[mapq > 10.0], 0)
        rect(xl, y, xr, y + yheight, col='red', border=NA)
        xl <- xr
        xr <- max(cum[mapq <= 10.0], 0)
        rect(xl, y, xr, y + yheight, col='brown', border=NA)
    }

    logmod <- function(x) {
        return(sign(x) * log10(abs(x) + 1))
    }

    inv_logmod <- function(x) {
        sx <- sign(x)
        return(sx * (10.0 ** (sx * x) - 1))
    }

    line_width <- 3

    # Plot the top plot
    par(mgp=c(3, 0.6, 0), cex.axis=0.7)
    sse_diff <- sse_incor1 - sse_incor2
    sse_diff <- sign(sse_diff) * log10(abs(sse_diff) + 1)
    cidp_diff <- incor_incor1 - incor_incor2
    cidp_diff <- sign(cidp_diff) * log10(abs(cidp_diff) + 1)
    par(mar = c(0.1, 2.1, 2.1, 2.1))
    plot(incor_cum, cidp_diff, type='l', col='red', xaxt='n', yaxt='n', bty='n', xlab='', ylab='')
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(0.95, 0.95, 0.95))
    lines(incor_cum, cidp_diff, type='l', col='red', lwd=line_width)
    axis(2, at=logmod(c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)),
            labels=c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000), las=2)
    abline(a=0, b=0)
    legend("bottomleft", legend=c(">40", "30s", "20s", "10s", "<10"),
           col=c("blue", "dodgerblue1", "green", "red", "brown"), pch=15, horiz=T, cex=1.2,
           bty="n", adj = c(0.0, NA))

    #
    # Plot the MAPQ bars
    #
    par(xpd=NA)

    # Aligner MAPQ bar
    par(mar = c(0, 2.1, 0, 2.1))
    plot(incor_cum, rep(0, length(incor_cum)), type="n", axes=F, xlab="", ylab="", ylim=c(0, 10))
    plot_mapq(incor_cum, mapq_orig, 0)
    text(0.03*(max(incor_cum) - min(incor_cum)), 5, "Aligner", col='white', adj = c(0.0, NA))

    # Qsim MAPQ bar
    par(mar = c(0, 2.1, 0, 2.1))
    plot(incor_cum, rep(0, length(incor_cum)), type="n", axes=F, xlab="", ylab="", ylim=c(0, 10))
    plot_mapq(incor_cum, mapq, 0)
    text(0.03*(max(incor_cum) - min(incor_cum)), 5, "Qsim", col='white', adj = c(0.0, NA))

    # Difference MAPQ bar - a bit harder since it's a bunch of rectangles
    #par(mar = c(0, 2.1, 0, 2.1))
    #plot(incor_cum, rep(0, length(incor_cum)), type="n", axes=F, xlab="", ylab="", ylim=c(0, 10))

    #plot(seq(1, length(mapq)), mapq - mapq_orig)

    #text(0.03*(max(incor_cum) - min(incor_cum)), 5, "Qsim - aligner", col='white', adj = c(0.0, NA))
    par(xpd=F)

    # Plot the bottom plot
    par(mar = c(2.1, 2.1, 0.1, 2.1))
    mn <- min(sse_diff)
    mx <- max(sse_diff)
    rn = mx - mn
    plot(sse_cum, sse_diff, type='l', col='red', bty='n', yaxt='n',
         ylim=c(mn, mn + rn), cex.axis=0.7, xlab='', ylab='')
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(0.95, 0.95, 0.95))
    lines(sse_cum, sse_diff, type='l', col='red', lwd=line_width)
    axis(2, at=logmod(c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)),
            labels=c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000), las=2)
    abline(a=0, b=0)

    auc <- sum(incor_incor1) - sum(incor_incor2)
    auc_pct_diff <- 100.0 * auc / sum(incor_incor1)
    sse <- sse_incor1[length(sse_incor1)] - sse_incor2[length(sse_incor2)]
    sse_pct_diff <- 100.0 * sse / sse_incor1[length(sse_incor1)]
    legend("bottomleft", sprintf("Relative change in AUC: %0.2f%%, SSE: %0.2f%%", auc_pct_diff, sse_pct_diff),
           cex=1.2, bty="n", adj = c(0.0, NA))

}
