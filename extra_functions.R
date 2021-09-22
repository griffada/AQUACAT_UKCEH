cdfPrimer <- function(RCM, period, method, NE, NH, thresh1, ws1, rn, savepath, chunks=T){
  
  if(file.exists(savepath)){ file.remove(savepath) }
  
  if(chunks){
    cs <- c(4,250)
  }else{
    cs <- NA
  }
  
  loc_dim <- ncdim_def(name="loc", units="Num", longname="Location",
                       vals=(1:NH))
  event_dim <- ncdim_def("event", units="Num", vals=1:NE, 
                         unlim=T, longname="Event")
  dpe_var <- ncvar_def("dpe","ProbOfExc", list(loc_dim, event_dim), -9999,
                       "Daily Probability of Exceedance", prec="double",
                       compression=2, chunksizes=cs)
  ape_var <- ncvar_def("ape","ProbOfExc",list(loc_dim, event_dim), -9999,
                       "Annual Probability of Exceedance",prec="double",
                       compression=2, chunksizes=cs)
  #ape_var20 <- ncvar_def("ape20","ProbOfExc",list(loc_dim, event_dim), -9999,
  #                   "Annual Probability of Exceedance (20/yr)",prec="double",
  #                   compression=2, chunksizes=cs)
  flow_var <- ncvar_def("flow", "cumecs",list(loc_dim, event_dim), -9999,
                        "Peak Flow",prec="double",
                        compression=2, chunksizes=cs)
  event_var <- ncvar_def("eventNo", "Matching OBS event", list(event_dim), -1,
                         "", prec="integer")
  row_var <- ncvar_def("row", "Num", list(loc_dim), -9999,
                       "Row", prec="integer")
  col_var <- ncvar_def("col", "Num", list(loc_dim), -9999,
                       "Column", prec="integer")
  north_var <- ncvar_def("northing", "metres", list(loc_dim), -9999,
                         "Northing", prec="double")
  east_var <- ncvar_def("easting", "metres", list(loc_dim), -9999,
                        "Easting", prec="double")
  
  nc.file <- nc_create(savepath, list(row_var, col_var, north_var, east_var,
                                      flow_var, dpe_var, ape_var,# ape_var20, 
                                      event_var),
                       force_v4 = T)
  
  ncvar_put(nc.file, row_var, unlist(rn[,1]))
  ncvar_put(nc.file, col_var, unlist(rn[,2]))
  ncvar_put(nc.file, north_var, unlist(rn[,3]))
  ncvar_put(nc.file, east_var, unlist(rn[,4]))
  
  ncatt_put(nc.file,0, "RCM", paste("RCM", RCM))
  ncatt_put(nc.file,0, "period", period)
  ncatt_put(nc.file,0, "event threshold", thresh1)
  ncatt_put(nc.file,0, "area lower limit", ws1)
  ncatt_put(nc.file,0, "Method", method)
  
  print(nc.file)
  nc_close(nc.file)
}

dpeApeComputer2 <- function(h, vals, obs_events, pars, thresh_val, ncin, 
                           events_per_year=2, rare_limit=1e-3){
  # Calculates daily and annual probability of exceedence based on observed flow
  # and GLO parameters using a hybrid empirical/modelled distribution.
  # Calculates 360-day conversion and 20-events-per-year conversion.
  # Saves directly to the netcdf file.
  # Returns number of extreme events, and annual RP of those events.
  
  
  ecd <- ecdf(c(-1,vals,1e8))

  gpa_poe <- 1 - lmomco::cdfglo(obs_events, pars)
  ecd_poe <- 1 - ecd(obs_events)
  gpa_poe[is.na(gpa_poe)] <- 1
  wh_ext <- obs_events > thresh_val
  wh_ext[is.na(wh_ext)] <- FALSE

  gpa_poe[wh_ext] <- gpa_poe[wh_ext]*(1 - ecd(thresh_val))
  gpa_poe[!wh_ext] <- ecd_poe[!wh_ext]

  valsape <- 1 - exp(-360*gpa_poe)
  valsape[is.na(valsape)] <- 1
  

  valsape[valsape < 1/5000] <- 1/5000
  gpa_poe[gpa_poe < 5e-7] <- 5e-7
  
  if(any(valsape < rare_limit)){
    print(paste("****", h))
    print(paste("obs beyond rare limit",
                paste(1/valsape[valsape < rare_limit], collapse=" ")))
  }
  
  ncvar_put(nc=ncin, varid="dpe", vals=gpa_poe, start=c(h,1),
            count=c(1,length(gpa_poe)))
  ncvar_put(nc=ncin, varid="ape", vals=valsape, start=c(h,1),
            count=c(1,length(gpa_poe)))
  
  # if("ape20" %in% names(ncin$var)){
  #   valsape20 <- 1 - exp(-20*gpa_poe)
  #   valsape20[is.na(valsape)] <- 1
  #   
  #   ncvar_put(nc=ncin, varid="ape20", vals=valsape20, start=c(h,1),
  #         count=c(1, length(gpa_poe)))
  #   
  #   if(any(valsape20 < rare_limit)){
  #     print(paste("****", h))
  #     print(paste("obs beyond rare limit (20/yr)",
  #               paste(1/valsape20[valsape20 < rare_limit], collapse=" ")))
  #   }
  #}

  
  return(c(sum(valsape < rare_limit), 1/min(valsape[valsape > (1/5000)])))
}

dpeApeComputer <- function(h, vals, obs_events, pars, thresh_val, ncin, 
                           events_per_year=2, rare_limit=1e-3){
  # Calculates daily and annual probability of exceedence based on observed flow
  # and GPA parameters using a hybrid empirical/modelled distribution.
  # Calculates 360-day conversion and 20-events-per-year conversion.
  # Saves directly to the netcdf file.
  # Returns number of extreme events, and annual RP of those events.
  
  parsGL <- vec2par(unlist(pars[3:5]), type='glo')
  parsGP <- vec2par(unlist(pars[6:8]), type='gpa')
  ecd <- ecdf(c(-1,vals,1e8))

  glo_poe <- 1 - lmomco::cdfgpa(obs_events, parsGL)
  gpa_poe <- 1 - lmomco::cdfgpa(obs_events, parsGP)
  
  gpa_poe <- pmin(gpa_poe, glo_poe)  #take the rarer of glo_poe and gpa_poe
  
  ecd_poe <- 1 - ecd(obs_events)
  gpa_poe[is.na(gpa_poe)] <- 1
  wh_ext <- (obs_events > thresh_val) & (gpa_poe < ecd_poe)
  wh_ext[is.na(wh_ext)] <- FALSE

  gpa_poe[wh_ext] <- gpa_poe[wh_ext]*(1 - ecd(thresh_val))
  gpa_poe[!wh_ext] <- ecd_poe[!wh_ext]

  valsape <- 1 - exp(-360*gpa_poe)
  valsape[is.na(valsape)] <- 1
  

  valsape[valsape < 1/5000] <- 1/5000
  gpa_poe[gpa_poe < -log(1-(1/5000))/360] <- -log(1-(1/5000))/360
  
  if(any(valsape < rare_limit)){
    print(paste("****", h))
    print(paste("obs beyond rare limit",
                paste(1/valsape[valsape < rare_limit], collapse=" ")))
  }
  
  ncvar_put(nc=ncin, varid="dpe", vals=gpa_poe, start=c(h,1),
            count=c(1,length(gpa_poe)))
  ncvar_put(nc=ncin, varid="ape", vals=valsape, start=c(h,1),
            count=c(1,length(gpa_poe)))
  
  return(c(sum(valsape < rare_limit), 1/min(valsape[valsape > (1/5000)])))
}


eventSummaryLine <- function(obs_events, i, threshold, r1=1:19914){
  ni <- eventNo[i]
  vvec <- sum(
    (ncvar_get(obs_events, "flow", start=c(1,i), count=c(-1,1)) > threshold),
    na.rm=T)
  avec <- ncvar_get(obs_events, "ape", start=c(1,i), count=c(-1,1))
  min_dvec <- min(
    ncvar_get(obs_events, "dpe", start=c(1,i), count=c(-1,1)),
    na.rm=T)
  D <- eventDayList[[jT]][[jW]][ni]
  L <- eventLList[[jT]][[jW]][ni]
  #ncl <- nclustevent(avec < (1-exp(-2*20/360)))
  #pk <- peakyfun(avec, (1-exp(-2*20/360)))

  return(list(i, D, L, vvec, min(avec), min_dvec, season(D), 0,0))
}

POT_extract_UKFE2 <- function (x, div = NULL, thresh = 0.975, Plot = TRUE){
    Low.Func <- function(TS) {
      L <- length(TS) - 2
      L1 <- length(TS) - 1
      L2 <- length(TS)
      Vec1 <- TS[1:L]
      Vec2 <- TS[2:L1]
      Vec3 <- TS[3:L2]
      P1 <- ifelse(Vec2 <= Vec1 & Vec2 <= Vec3 & Vec1 != Vec2,
                   Vec2, NA)
      return(P1)
    }
    P.Func <- function(TS) {
      L <- length(TS) - 2
      L1 <- length(TS) - 1
      L2 <- length(TS)
      Vec1 <- TS[1:L]
      Vec2 <- TS[2:L1]
      Vec3 <- TS[3:L2]
      P1 <- ifelse(Vec2 >= Vec1 & Vec2 >= Vec3 & Vec1 != Vec2,
                   Vec2, NA)
      return(P1)
    }
    VP <- function(j, mu) {
      maxll <- suppressWarnings(max(which(lows[1:j] <= mu),
                                    na.rm = T))
      if (maxll == -Inf) {
        maxll <- j
      } else {
        maxll <- maxll
      }
      minlr <- suppressWarnings(min(which(lows[j:length(lows)] <=
                                            mu), na.rm = T))
      if (minlr == Inf) {
        minlr <- j
      } else {
        minlr <- j + (minlr - 1)
      }
      if (peaks[j] == max(peaks[maxll:minlr], na.rm = T)) {
        vp <- peaks[j]
      } else {
        vp <- NA
      }
      return(vp)
    }
    NAs <- FALSE
    if (class(x) != "data.frame") {
      if (is.null(div)) {
        mu <- mean(x, na.rm = TRUE)
      }
      else {
        mu <- div
      }
      if (mu >= quantile(x, thresh, na.rm = TRUE))
        stop("The event division must be significantly lower than the event threshold")
      QThresh <- as.numeric(quantile(x, thresh, na.rm = TRUE))
      MinMuP <- min(which(x <= mu))
      MaxMuP <- max(which(x <= mu))
      PkBegin <- max(x[1:MinMuP])
      PkEnd <- max(x[MaxMuP:length(x)])
      x <- x[MinMuP:MaxMuP]
      lows <- Low.Func(x)
      peaks <- P.Func(x)
      MinMuL <- min(which(lows <= mu))
      MaxMuL <- max(which(lows <= mu))
      pt.ind <- which(peaks > QThresh)
      pt <- peaks[pt.ind]
      l <- length(pt.ind)
      POT <- NULL
      {
        for (i in 1:l) {
          POT[i] <- VP(pt.ind[i], mu)
        }
      }
      if (PkBegin > QThresh) {
        POT <- append(PkBegin, POT)
      }
      if (PkEnd > QThresh) {
        POT <- append(POT, PkEnd)
      }
      if (NAs == TRUE) {
        POT <- POT
      }
      else {
        POT <- POT[which(is.na(POT) == FALSE)]
      }
      return(POT)
    }
    else {
      if (is.null(div)) {
        mu <- mean(x[, 2], na.rm = TRUE)
      }
      else {
        mu <- div
      }
      if (mu >= quantile(x[, 2], thresh, na.rm = TRUE))
        stop("The event division must be significantly lower than the event threshold")
      QThresh <- as.numeric(quantile(x[, 2], thresh, na.rm = TRUE))
      MinMuP <- min(which(x[, 2] <= mu), na.rm = TRUE)
      MaxMuP <- max(which(x[, 2] <= mu), na.rm = TRUE)
      PkBegin <- which(x[1:MinMuP, 2] == max(x[1:MinMuP, 2],
                                             na.rm = TRUE))
      PkEnd <-
        which(x[MaxMuP:length(x[, 2]), 2] == max(x[MaxMuP:length(x[,
                                                                   2]), 2], na.rm = TRUE))
      DBegin <- x[PkBegin,]
      DEnd <- x[((MaxMuP - 1) + PkEnd),]
      colnames(DBegin) <- c("Date", "peak")
      colnames(DEnd) <- c("Date", "peak")
      xUse <- x[MinMuP:MaxMuP,]
      lows <- Low.Func(xUse[, 2])
      peaks <- P.Func(xUse[, 2])
      pt.ind <- which(peaks > QThresh)
      pt <- peaks[pt.ind]
      L <- length(pt.ind)
      POT <- NULL
      {
        for (i in 1:L) {
          POT[i] <- VP(pt.ind[i], mu)
        }
      }
      POT.Dates <- (xUse[, 1][pt.ind]) + 1
      res <- data.frame(POT.Dates, POT)
      colnames(res) <- c("Date", "peak")
      if (DBegin$peak > QThresh) {
        res <- rbind(DBegin, res)
      }
      if (DEnd$peak > QThresh) {
        res <- rbind(res, DEnd)
      }
      rownames(res) <- seq(1, length(res[, 2]))
      if (NAs == TRUE) {
        res <- res
      }
      else {
        res <- res[which(is.na(res$peak) == FALSE),]
      }
      if (Plot == TRUE) {
        plot(
          x,
          type = "l",
          main = "Peaks over threshold",
          ylab = "Quantile",
          xlab = "Date"
        )
        abline(h = quantile(x[, 2], thresh, na.rm = TRUE),
               col = "blue")
        points(res, col = "red")
        abline(h = mu, col = rgb(0, 0.7, 0.3))
      }
      Years <- as.numeric((x[length(x[, 1]), 1] - x[1, 1]) / 365.25)
      PPY <- length(res[, 1]) / Years
      # print(paste("Peaks per year:", format(PPY, trim = TRUE),
      #     sep = " "))
      return(res)
    }
}
