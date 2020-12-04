#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-08-12
#
# HT functions for spatial coherence on big datasets.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-08-12
#
#~~~~~~~~~~~~~~~~~~~~~~~


#### FUNCTIONS ####-----------------------------------------
plottexmex <- function (x, plots = "gpd", main = "", plotnum=NULL, ...) {
  if (casefold(plots) == "gpd") {
    which <- 1
  }
  else {
    which <- 2
  }
  d2 <- dim(x$boot[[1]][[which]])
  pointEst <- x$simpleDep
  condVar <- names(x$simpleMar$data)[x$which]
  margins <- x$margins
  x <- x$boot
  co <- unlist(lapply(x, function(z, wh) z[[wh]], wh = which))
  co <- array(co, dim = c(d2[1], d2[2], length(co)/prod(d2)))
  lco <- list(length = prod(d2))
  for (i in 1:d2[2]) {
    for (j in 1:d2[1]) {
      lco[[j + d2[1] * (i - 1)]] <- co[j, i, ]
    }
  }
  cn <- colnames(x[[1]][[which]])
  rn <- rownames(x[[1]][[which]])
  if (which == 2) {
    cn <- paste(cn, "|", condVar)
  }
  labs <- paste(rep(rn, length(cn)),
                rep(cn, each = switch(which, 2, 6)), sep = "  ")
  fun <- function(X, z, label, ...) {
    hist(z[[X]], prob = TRUE, xlab = label[X], main = main, ...)
    lines(density(z[[X]], n = 100))
    invisible()
  }
  if (which == 1) {
    if(is.null(plotnum)){
      plotnum <- 1:d2[2]
    }
    lapply(plotnum, fun, z = lco, label = labs, ...)
  }
  if (which == 2) {
    fun <- function(X, z, label, ...) {
      offset <- (X - 1) * 6
      plot(lco[[offset + 1]], lco[[offset + 2]],
           xlab = labs[offset + 1], ylab = labs[offset + 2],
           main = main, ...)
      points(pointEst[1, X], pointEst[2, X], pch = "@", col = 2)
      if (margins[[1]] == "gumbel") {
        plot(lco[[offset + 3]], lco[[offset + 4]],
             xlab = labs[offset + 3], ylab = labs[offset + 4],
             main = main, ...)
        points(pointEst[3, X], pointEst[4, X], pch = "@", col = 2)
      }
    }
    if(is.null(plotnum)){
      plotnum <- 1:d2[2]
    }
    lapply(plotnum, fun, z = lco, label = labs, ...)
  }
  invisible()
}


plotpredmex <- function(x, pch = c(1, 3, 20), col = c(2, 8, 3),
                        cex = c(1, 1, 1), ask = TRUE, plotnum=NULL, ...){
  d <- dim(x$data$simulated)[[2]] - 1
  if (prod(par("mfrow")) < d) {
    if (ask) {
      op <- par(ask = TRUE)
      on.exit(par(op))
    }
  }
  xdat <- x$data$real[, 1]
  upts <- seq(from = 0.001, to = 1 - 1e-04, len = 100)
  xpts <- texmex:::revTransform(upts, data = x$data$real[, 1],
                                qu = mean(x$data$real[, 1] < x$mth[1]),
                                th = x$mth[1],
                                sigma = x$gpd.coef[3, 1],
                                xi = x$gpd.coef[4, 1])
  if(is.null(plotnum)){
    plotnum <- 2:(dim(x$data$real)[[2]])
  }
  for (i in plotnum) {
    ydat <- x$data$real[, i]
    xlimits <- range(xdat, x$data$simulated[, 1])
    ylimits <- range(ydat, x$data$simulated[, i])
    plot(xdat, ydat, xlim = xlimits, ylim = ylimits,
         xlab = names(x$data$simulated)[1], 
         ylab = names(x$data$simulated)[i], type = "n", ...)
    points(x$data$simulated[x$data$CondLargest, 1],
           x$data$simulated[x$data$CondLargest, i],
           col = col[3], pch = pch[3], cex = cex[3])
    points(x$data$simulated[!x$data$CondLargest, 1],
           x$data$simulated[!x$data$CondLargest, i],
           col = col[2], pch = pch[2], cex = cex[2])
    points(xdat, ydat, pch = pch[1], col = col[1], cex = cex[1])
    abline(v = x$data$pth, lty = 2, col = 3)
    ypts <- texmex:::revTransform(upts, data = x$data$real[, i],
                                  qu = mean(x$data$real[, i] < x$mth[i]),
                                  th = x$mth[i],
                                  sigma = x$gpd.coef[3, i],
                                  xi = x$gpd.coef[4, i])
    lines(xpts, ypts, col = 3)
  }
  invisible()
}


mexMonteCarlo2 <- function (nSample, mexList, mult = 10)
  # edited from texmex::mexMonteCarlo to get the return periods out.
{
  d <- length(mexList)
  data <- mexList[[1]]$margins$data
  margins <- mexList[[1]]$dependence$margins
  nData <- dim(data)[1]
  which <- sample(1:nData, size = nSample, replace = TRUE)
  MCsampleOriginal <- data[which, ]
  dataLaplace <- texmex:::mexTransform(mexList[[1]]$margins, margins = margins, 
                                       method = "mixture")$transformed
  MCsampleLaplace <- dataLaplace[which, ]
  whichMax <- apply(MCsampleLaplace, 1, which.max)
  dth <- sapply(mexList, function(l) l$dependence$dth)
  dqu <- sapply(mexList, function(l) l$dependence$dqu)
  whichMaxAboveThresh <- sapply(1:nSample,
                                function(i) MCsampleLaplace[i, whichMax[i]] >= dth[whichMax[i]])
  
  
  mexKeep <- vector("list", d)
  mexLap <- vector("list", d)
  
  for(i in 1:d){
    mc <- texmex:::predict.mex(mexList[[i]], pqu = dqu[i], 
                               nsim = nSample * d * mult)
    
    mexKeep[[i]] <- mc$data$simulated[mc$data$CondLargest,
                                      order(c(i, c(1:d)[-i]))]
    
    mexLap[[i]] <- mc$data$transformed[mc$data$CondLargest,
                                       order(c(i, c(1:d)[-i]))]
    
  }
  nR <- rep(0, d)
  names(nR) <- names(data)
  for (i in 1:d) {
    replace <- whichMax == i & whichMaxAboveThresh
    nReplace <- sum(replace)
    if (nReplace > 0) {
      nR[i] <- nReplace
      MCsampleOriginal[replace, ] <- as.matrix(mexKeep[[i]])[1:nReplace, ]
      MCsampleLaplace[replace, ]  <- as.matrix(mexLap[[i]])[1:nReplace,  ]
    }
  }
  res <- list(nR = nR,
              MCsample = MCsampleOriginal,
              MCsample_L = MCsampleLaplace,
              whichMax = whichMax, 
              whichMaxAboveThresh = whichMaxAboveThresh,
              q2pp = mexList[[1]]$dependence$margins$q2p,
              MCsample_P = mexList[[1]]$dependence$margins$q2p(MCsampleLaplace)
  )
  oldClass(res) <- "mexMC"
  res
}


mexMonteCarloBIG <- function (nSample, mexList, di = 2, dj = 3, mult = 10) 
{
  d <- length(mexList)
  data <- mexList[[1]]$margins$data
  margins <- mexList[[1]]$dependence$margins
  nData <- dim(data)[1]
  which <- sample(1:nData, size = nSample, replace = TRUE)
  MCsampleOriginal <- data[which, ]
  dataLaplace <- texmex:::mexTransform(mexList[[1]]$margins, margins = margins,
                                       method = "mixture")$transformed
  MCsampleLaplace <- dataLaplace[which, ]
  whichMax <- apply(MCsampleLaplace, 1, which.max)
  dth <- sapply(mexList, function(l) l$dependence$dth)
  dqu <- sapply(mexList, function(l) l$dependence$dqu)
  whichMaxAboveThresh <- sapply(1:nSample, function(i) MCsampleLaplace[i,whichMax[i]] >= dth[whichMax[i]])
  mexKeep <- lapply(di:dj, function(i) {
    mc <- predict.mex(mexList[[i]], pqu = dqu[i], nsim = nSample * d * mult)
    mc$data$simulated[mc$data$CondLargest, order(c(i, c(1:d)[-i]))]
  })
  nR <- rep(0, (dj-di)+1)
  names(nR) <- names(data)
  for (i in di:dj) {
    replace <- whichMax == i & whichMaxAboveThresh
    nReplace <- sum(replace)
    if (nReplace > 0) {
      nR[i+(di-1)] <- nReplace
      MCsampleOriginal[replace, ] <- as.matrix(mexKeep[[i]])[1:nReplace, ]
    }
  }
  res <- list(nR = nR, MCsample = MCsampleOriginal, whichMax = whichMax, 
              whichMaxAboveThresh = whichMaxAboveThresh)
  oldClass(res) <- "mexMC"
  res
}


