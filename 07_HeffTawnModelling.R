#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-06-17
#
# Using HT model for spatial coherence.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-06-17
#
#~~~~~~~~~~~~~~~~~~~~~~~

##### SETUP #####-----------------------------------------------------------

library(texmex)
library(readr)
library(dplyr)
library(ncdf4)
library(extRemes)

if (substr(osVersion,1,3) == "Win") {
  ncname <- "S:/CodeABG/InterimData/dmflow_copy.nc"  # rechunked for spaceslices
  ncname2 <- "S:/run_hmfg2g/outputs/dmflow_RCM01_198012_201011_out.nc" 
  # pre-rechunk
  wd <- "S:/"
  wd_id <- "S:/CodeABG/InterimData/"
  wd_cd <- "S:/CodeABG/"
  
} else {
  ncname <- "/prj/aquacat/CodeABG/InterimData/dmflow_copy.nc"
  ncname2 <- "S:/run_hmfg2g/outputs/dmflow_RCM01_198012_201011_out.nc"
  wd <- "/prj/aquacat/"
  wd_id <- "/prj/aquacat/CodeABG/InterimData/"
  wd_cd <- "/prj/aquacat/CodeABG/"
}

##### DATA #####-----------------------------------------------------------

ND <- 10800 # Number of days

# river network
rn <- read_csv(paste0(wd_id, "hasData2.csv"))
NH <- nrow(rn)


# threshold for inundation
threshVal <- c(5/365, 2/365, 1/365, 0.2/365, 0.1/365)
threshName <- c("POT5", "POT2", "POT1", "Q5", "Q10")
NT <- length(threshVal)



# cut-off bound for widespread event
wsBound <- c(0.05, 0.02, 0.01, 0.005, 0.001)
wsName <- c("pc5", "pc2", "pc1", "pc05", "pc01")
NW <- length(wsBound)



# lists of which days different thresholds were exceeded at different points
# NT lists of NW lists
threshDayExcList <- readRDS(paste0(wd_id, "threshDayExcList2.rds"))


# matrix of threshold value (col) at a given cell (row)
threshMat <- read_csv(paste0(wd_id, "threshMat2.csv"))
#dim(threshMat)  =  19914 x 5


# eventLList length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
load(paste0(wd_id, "eventLists03.RDa")) 
NE <- length(eventDayList[[2]][[4]]) # POT2, 2% inun.


# timewise maxima at each cell for each event ((NE + 2) x NH)
eventDF <- readr::read_csv(paste0(wd,"Data/eventdf_POT2_pc05.csv"))


# PoE under different computations with extra data. Tidy format.
present <- readr::read_csv(paste0(wd,"Data/present2_returnlevels.csv"))
dummy <- present %>% filter(loc < 4)

#dummy <- readr::read_csv(paste0(wd, "Data/dummy2_returnlevels.csv"))





##### HEFFTAWN CODING #####------------------------------------------------

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


#### TRIAL FIT ####-------------------------------------------

dummy_melt <- dcast(dummy, eventNo~loc, value.var="gpp")
rownames(dummy_melt) <- paste0("E",dummy_melt$eventNo)
dummy_melt <- dummy_melt[,-1]
colnames(dummy_melt) <- paste0("L", 1:3)

ST <- Sys.time()
  # fits everything else conditional on location 1.
  mmex <- mex(dummy_melt, mqu=.7, penalty="none", dqu=.7, which=1)
print(Sys.time() - ST)

ST <- Sys.time()
  myboot <- bootmex(mmex, R=100)  # R=10 is quite low.
print(Sys.time() - ST)
plottexmex(myboot, plots="dependence")

ST <- Sys.time()
  mypred <- predict(myboot,  pqu=.99)
print(Sys.time() - ST)
plotpredmex(mypred)


#### GET THE COEFFICIENTS OUT ####---------------------------
ab_mmex <- mmex$dependence$coefficients
marg_mmex <- sapply(mmex$margins$models, function(L){L$coefficients})





