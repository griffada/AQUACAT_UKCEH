#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-02-21
#
# Rewritten functions to perform HT modelling on big dataset.
# Takes texmex and pulled it apart, making use of a small number of internal
# texmex functions.
# Generally, this is an aim to reduce overheads by removing copying of large
# objects and large list structures.
#
# This is tested in slimline_testing2.R
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-08-14, edited 2020-09-03
#
# OUTPUTS: None
#
#~~~~~~~~~~~~~~~~~~~~~~~

library(texmex)


migpd_slim <- function(mth, mqu, penalty = "gaussian", maxit = 10000, 
          trace = 0, verbose = FALSE, priorParameters = NULL, cov = "observed", 
          family = gpd){
  # Extension of **migpd** from texmex, removing data from the arguments.
  # Reduces the number of outputs to lower overheads.
  # fetches DATA from the global environment.
  
  if(!exists("DATA")){
    stop("Needs DATA object; one event per row, one location per column.")
  }
  
  theCall <- match.call()
  if (is.null(colnames(DATA))) {
    colnames(DATA) <- paste(rep("Column", ncol(DATA)), 
                            1:ncol(DATA), sep = "")
  }
  d <- dim(DATA)[2]
  if (missing(mth) & missing(mqu)) 
    stop("you must provide one of mth or mqu")
  if (!missing(mth) & !missing(mqu)) 
    stop("you must provide precisely one of mth or mqu")
  if (!(family$name %in% c("GPD", "CGPD"))) {
    stop("family should be either gpd or cgpd")
  }
  if (!missing(mth)) 
    mth <- rep(mth, length = d)
  if (!missing(mqu)) 
    mqu <- rep(mqu, length = d)
  if (missing(mqu)) 
    mqu <- sapply(1:d, function(i, mth) 1 - mean(DATA[, i] > mth[i]), mth = mth)
  if (missing(mth)) 
    mth <- sapply(1:d, function(i, prob) quantile(DATA[,  i], prob = prob[i]), prob = mqu)
  if (penalty %in% c("quadratic", "gaussian") & 
      is.null(priorParameters)) {
    gp = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2))
    priorParameters <- vector("list", length = length(mth))
    for (i in 1:length(mth)) priorParameters[[i]] <- gp
    names(priorParameters) <- dimnames(DATA)[[2]]
  }
  else if (penalty %in% c("quadratic", "gaussian")) {
    nm <- names(priorParameters)
    if (is.null(nm)) {
      stop("priorParameters must be a named list")
    }
    else if (any(!is.element(nm, dimnames(DATA)[[2]]))) {
      stop("the names of priorParameters must match the column names of the DATA")
    }
  }
  #####
  wrapgpd_slim <- function(i, mth, penalty, maxit, verbose, trace, 
                      priorParameters) {
    if (verbose) 
      cat("Fitting model", i, "\n")
    if (!is.null(priorParameters)) 
      priorParameters <- 
        priorParameters[[(1:length(priorParameters))[names(priorParameters) == dimnames(DATA)[[2]][i]]]]
    x <- c(DATA[, i])
    mth <- mth[i]
    evm(x, th = mth, penalty = penalty, priorParameters = priorParameters, 
        maxit = maxit, trace = trace, cov = cov, family = family)
  }
  #####
  modlist <- lapply(1:d, wrapgpd_slim, penalty = penalty, 
                    mth = mth, verbose = verbose, priorParameters = priorParameters, 
                    maxit = maxit, trace = trace)
  if (length(dimnames(DATA)[[2]]) == dim(DATA)[[2]]) {
    names(modlist) <- dimnames(DATA)[[2]]
  }
  names(mth) <- names(mqu) <- dimnames(DATA)[[2]]
  res <- list(models = modlist, mth = mth, mqu = mqu)
  oldClass(res) <- "migpd"
  invisible(res)
}


mexTransform_slim <- function(marginfns, mth, r=NULL, method = "mixture",
                              divisor = "n+1", na.rm = TRUE){ 
  
  # if required, r output from makeReferenceMarginalDistribution
  # Extension of **mexTransform** from texmex, removing data and models from the arguments.
  # Reduces the number of outputs to lower overheads.
  # fetches MODELS and DATA from the global environment
  
  if(!exists("DATA")){
    stop("Needs DATA object; one event per row, one location per column.")
  }
  if(!exists("MODELS")){
    stop("Needs MODELS object; a list of one marginal model per location, output from migpd, $models item.")
  }
  
    
  if (!is.element(method, c("mixture", "empirical")))
    stop("method should be either 'mixture' or 'empirical'")
  if (!is.element(divisor, c("n", "n+1")))
    stop("divisor can be 'n' or 'n+1'")
  
  if (is.null(r)){
    r <- list(mth=mth)
    r$transData <- lapply(1:dim(DATA)[2], function(i)DATA[,i])
  }
  
  #####
  transFun <- function(i, th, divisor, method){
    x <- DATA[, i]
    r <- DATA[, i]
    mod <- MODELS[[i]]
    th <- th[i]
    
    if (divisor == "n") divisor <- length(r)
    else if (divisor == "n+1") divisor <- length(r) + 1
    
    ox <- order(x)
    r <- sort(r)
    run <- rle(r)
    p <- cumsum(run$lengths) / divisor
    p <- rep(p, run$lengths) # this calculated from r
    
    Femp <- p[sapply(x,function(y) which.min(abs(r-y)))]
    if (method == "mixture"){
      sigma <- exp(mod$coefficients[1])
      xi <- mod$coefficients[2]
      
      Para <- (1 + xi * (x - th) / sigma) ^ (-1 / xi) 
      # this calculated from model fitted to r but data values are x
      Para <- 1 - mean(r > th) * Para
      Para[Para==1] <- 1 - 1e-8  # log doesn't like zeroes.
      res <- ifelse((x <= th), Femp, Para)
    }
    else res <- Femp
    
    res[ox] <- sort(res)
    res
    
  } # Close transfun
  #####
  res <- sapply(1:ncol(DATA), transFun, th = r$mth,
                divisor = divisor, method=method)
  dimnames(res) <- list(NULL, names(MODELS))
  
  x <- list(transformed=marginfns$p2q(res))
  invisible(x)
}

PosGumb.Laplace.negloglik <- function (yex, ydep, a, b, m, s, constrain, v, aLow){
  BigNumber <- 10^40
  WeeNumber <- 10^(-10)
  if (a < aLow[1] | s < WeeNumber | a > 1 - WeeNumber | b > 1 - WeeNumber) {
    res <- BigNumber
  }
  else {
    mu <- a * yex + m * yex^b
    sig <- s * yex^b
    res <- sum(0.5 * log(2 * pi) + log(sig) + 0.5 * ((ydep -  mu)/sig)^2)
    if (is.infinite(res)) {
      if (res < 0) {
        res <- -BigNumber
      }
      else {
        res <- BigNumber
      }
      warning("Infinite value of Q in mexDependence")
    }
    else if (constrain) {
      zpos <- range(ydep - yex)
      z <- range((ydep - yex * a)/(yex^b))
      zneg <- range(ydep + yex)
      if (!texmex:::ConstraintsAreSatisfied(a, b, z, zpos, zneg, v)) {
        res <- BigNumber
      }
    }
  }
  res
}

PosGumb.Laplace.negProfileLogLik <- function (yex, ydep, a, b, constrain, v, aLow) {
  Z <- (ydep - yex * a)/(yex^b)
  m <- mean(Z)
  s <- sd(Z)
  res <- PosGumb.Laplace.negloglik(yex, ydep, a, b, m = m, 
                                   s = s, constrain, v, aLow = aLow)
  res <- list(profLik = res, m = m, s = s)
  res
}

mexDependence_slim <- function (which, dqu, mth, mqu=0.7, margins = "laplace",
                                constrain = TRUE, v = 10, maxit = 1e+06,
                                start = c(0.01, 0.01), marTransform = "mixture", 
                                referenceMargin = NULL, marginsTransformed = NULL,
                                nOptim = 1, zspot=FALSE){
  
  # Extension of **mexDependence** from texmex, removing data from the arguments.
  # Reduces the number of outputs to lower overheads.
  # fetches DATA from the global environment
  
  if(!exists("DATA")){
    stop("Needs DATA object; one event per row, one location per column.")
  }
  
  theCall <- match.call()
  errcode <- 0
  marginfns <- list(casefold(margins),
      p2q = switch(casefold(margins),
                   gumbel = function(p) -log(-log(p)),
                   laplace = function(p) ifelse(p <  0.5, log(2 * p), -log(2 * (1 - p)))),
      q2p = switch(casefold(margins),
                   gumbel = function(q) exp(-exp(-q)),
                   laplace = function(q) ifelse(q <  0, exp(q)/2, 1 - 0.5 * exp(-q))))
  
  if(!is.null(marginsTransformed)){
    x <- list(transformed = marginsTransformed)
  }else{
    x <- mexTransform_slim(marginfns = marginfns, mth = mth, method = marTransform, 
                      r = referenceMargin)
  }
  x$referenceMargin <- referenceMargin
  if (marginfns[[1]] == "gumbel" & constrain) {
    warning("With Gumbel margins, you can't constrain, setting constrain=FALSE")
    constrain <- FALSE
  }
  if (missing(which)) {
    message("Missing 'which'. Conditioning on", dimnames(x$transformed)[[2]][1], 
            "\n")
    which <- 1
  }
  else if (length(which) > 1) 
    stop("which must be of length 1")
  else if (is.character(which)) 
    which <- match(which, dimnames(x$transformed)[[2]])
  if (missing(dqu)) {
    message(paste("Assuming same quantile for dependence thesholding as was used\n",
                  "to fit corresponding marginal model...\n"))
    dqu <- mqu[which]
  }
  dth <- quantile(x$transformed[, which], dqu)
  dependent <- (1:(dim(DATA)[[2]]))[-which]
  if (length(dqu) < length(dependent)) 
    dqu <- rep(dqu, length = length(dependent))
  aLow <- ifelse(margins[[1]] == "gumbel", 10^(-10), 
                 -1 + 10^(-10))
  if (missing(start)) {
    start <- c(0.01, 0.01)
  }
  else if (inherits(start, "mex")) {
    start <- start$dependence$coefficients[1:2, ]
  }
  if (length(start) == 2) {
    start <- matrix(rep(start, length(dependent)), nrow = 2)
  }
  if (length(start) != 2 * length(dependent)) {
    stop(paste("start should be of type 'mex' or be a vector of length 2, or be a matrix",
    "with 2 rows and ncol equal to the number of dependence models to be estimated"))
  }
  #####
  qfun <- function(X, yex, wh, aLow, margins, constrain, v, maxit, start) {
    #print(X)
    Qpos <- function(param, yex, ydep, constrain, v, aLow) {
      a <- param[1]
      b <- param[2]
      res <- PosGumb.Laplace.negProfileLogLik(yex, ydep, a, b, constrain, v, aLow)
      res$profLik
    }
    o <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit), 
                   yex = yex[wh], ydep = X[wh], constrain = constrain, 
                   v = v, aLow = aLow))
    if (inherits(o, "try-error")) {
      if(interactive()) browser()
      warning("Error in optim call from mexDependence")
      errcode <- 101
      o <- as.list(o)
      o$par <- rep(NA, 6)
      o$value <- NA
    }
    else if (o$convergence != 0) {
      warning("Non-convergence in mexDependence")
      errcode <- 102
      o <- as.list(o)
      o$par <- rep(NA, 6)
    }
    else if (nOptim > 1) {
      for (i in 2:nOptim) {
        o <- try(optim(par = o$par, fn = Qpos, control = list(maxit = maxit), 
                       yex = yex[wh], ydep = X[wh], constrain = constrain, 
                       v = v, aLow = aLow), silent = TRUE)
        if (inherits(o, "try-error")) {
          warning("Error in optim call from mexDependence")
          errcode <- 103
          o <- as.list(o)
          o$par <- rep(NA, 6)
          o$value <- NA
          (break)()
        }
        else if (o$convergence != 0) {
          warning("Non-convergence in mexDependence")
          errcode <- 104
          o <- as.list(o)
          o$par <- rep(NA, 6)
          (break)()
        }
      }
    }
    if (!is.na(o$par[1])) {
      if (margins == "gumbel" & o$par[1] <= 10^(-5) & 
          o$par[2] < 0) {
        lo <- c(10^(-10), -Inf, -Inf, 10^(-10), -Inf, 
                10^(-10))
        Qneg <- function(yex, ydep, param) {
          param <- param[-1]
          b <- param[1]
          cee <- param[2]
          d <- param[3]
          m <- param[4]
          s <- param[5]
          obj <- function(yex, ydep, b, cee, d, m, s) {
            mu <- cee - d * log(yex) + m * yex^b
            sig <- s * yex^b
            log(sig) + 0.5 * ((ydep - mu)/sig)^2
          }
          res <- sum(obj(yex, ydep, b, cee, d, m, s))
          res
        }
        o <- try(optim(c(0, 0, 0, 0, 0, 1), Qneg, method = "L-BFGS-B", 
                       lower = lo, 
                       upper = c(1, 1 - 10^(-10), Inf, 1 - 10^(-10), Inf, Inf),
                       yex = yex[wh], ydep = X[wh]), 
                 silent = TRUE)
        if (inherits(o, "try-error") || o$convergence != 
            0) {
          warning("Non-convergence in mexDependence")
          errcode <- 105
          o <- as.list(o)
          o$par <- rep(NA, 6)
        }
      }
      else {
        Z <- (X[wh] - yex[wh] * o$par[1])/(yex[wh]^o$par[2])
        o$par <- c(o$par[1:2], 0, 0, mean(Z), sd(Z))
      }
    }
    c(o$par[1:6], o$value)
  }
  #####
  yex <- c(x$transformed[, which])
  wh <- yex > unique(dth)
  
  res <- sapply(1:length(dependent), 
                function(X, dat, yex, wh, aLow, margins, constrain, v, maxit, start){
                  qfun(dat[, X], yex, wh, aLow, margins, constrain, v, maxit, start[, X])},
                dat = as.matrix(x$transformed[, dependent]), yex = yex, 
                wh = wh, aLow = aLow, margins = marginfns[[1]], constrain = constrain, 
                v = v, maxit = maxit, start = start)
  
  loglik <- -res[7, ]
  res <- matrix(res[1:6, ], nrow = 6)
  dimnames(res)[[1]] <- c(letters[1:4], "m", "s")
  dimnames(res)[[2]] <- dimnames(x$transformed)[[2]][dependent]
  gdata <- as.matrix(x$transformed[wh, -which])
  ####
  tfun <- function(i, data_temp, yex, a, b, cee, d) {
    data_temp <- data_temp[, i]
    a <- a[i]
    b <- b[i]
    cee <- cee[i]
    d <- d[i]
    if (is.na(a)) 
      rep(NA, length(data_temp))
    else {
      if (a < 10^(-5) & b < 0) 
        a <- cee - d * log(yex)
      else a <- a * yex
      (data_temp - a)/(yex^b)
    }
  }
  ####
  z <- try(sapply(1:(dim(gdata)[[2]]), tfun, data_temp = gdata, 
                  yex = yex[wh], a = res[1, ], b = res[2, ],
                  cee = res[3,], d = res[4, ]))
  if (inherits(z, c("Error", "try-error"))) {
    errcode <- 106
    z <- matrix(nrow = 0, ncol = dim(DATA)[[2]] - 1)
  }
  else if (!is.array(z)) {
    z <- matrix(nrow = 0, ncol = dim(DATA)[[2]] - 1)
    errcode <- 107
  }
  dimnames(z) <- list(NULL, dimnames(x$transformed)[[2]][dependent])
  if(zspot){
    print(dim(z))
  }
  COEFFS[,,which] <<- res
  if(is.array(z) && !inherits(z, c("Error", "try-error"))){
    Z[1:(dim(z)[1]),1:(dim(z)[2]),which] <<- z
  }
  
  res2 <- list(dth = unique(dth), 
               dqu = unique(dqu), which = which,
               conditioningVariable = colnames(DATA)[which], 
               loglik = loglik, marginfns = marginfns, constrain = constrain, 
               v = v)
  
  #oldClass(res2) <- "mexDependence"  # A bit of a lie, but helps things work.
  
  output <- list(margins = list(#transformed=x$transformed,
                                referenceMargin=x$referenceMargin),
                 dependence = res2,
                 errcode = errcode,
                 zspot = dim(z))
  #oldClass(output) <- "mex" # A bit of a lie but helps things work.
  output
}




mexMonteCarlo_slim <- function(marginfns, referenceMargin=NULL,
                               mth, mqu, nSample, mexList, mult = 10){
  
  # Extension of **mexMonteCarlo** from texmex, removing data from the arguments.
  # Reduces the number of outputs to lower overheads.
  # fetches DATA from the global environment
  # fetches TRANSFORMED from the global environment
  
  d <- length(mexList)

  nData <- dim(DATA)[1]
  which <- sample(1:nData, size = nSample, replace = TRUE)
  MCsampleOriginal <- DATA[which, ]
  MCsampleLaplace <- TRANSFORMED[which, ]
  whichMax <- apply(MCsampleLaplace, 1, which.max)
  dth <- sapply(mexList, function(l) l$dependence$dth)
  dqu <- sapply(mexList, function(l) l$dependence$dqu)
  whichMaxAboveThresh <- sapply(1:nSample, 
                          function(i) MCsampleLaplace[i, whichMax[i]] >= dth[whichMax[i]])
  
  nReplace <- sapply(1:d, function(i){sum(whichMax==i & whichMaxAboveThresh)})
  nR <- rep(0, d)
  names(nR) <- names(DATA)
  for (i in 1:d) {
    if(d < 50 | i %% 10 == 0){
      print(paste0("Currently processing at ", (i/d)*100, "%"))
    }
    replace <- whichMax == i & whichMaxAboveThresh
    if (nReplace[i] > 0) {
      MCsampleOriginal[replace, ] <- predict.mex_slim(which=i,
                                      referenceMargin=referenceMargin,
                                      marginfns=marginfns,
                                      constrain=mexList[[i]]$dependence$constrain,
                                      coeffs_in = COEFFS[,,i],
                                      z_in = Z[,,i],
                                      pqu = dqu[i],
                                      mth=mth,
                                      mqu=mqu,
                                      nsim = nSample * d * mult,
                                      d=d)[1:nReplace[i],]
    }
  }
  res <- list(nR = nReplace, MCsample = MCsampleOriginal, whichMax = whichMax, 
              whichMaxAboveThresh = whichMaxAboveThresh)
  # oldClass(res) <- "mexMC" # A bit of a lie, but keeps things smooth
  res
}


coef.migpd_slim <- function(mth, mqu, ...){
  
  # Extension of **coef.migpd** from texmex, removing data from the arguments.
  # Reduces the number of outputs to lower overheads.
  # fetches MODELS from the global environment
  
  if(!exists("MODELS")){
    stop("Needs MODELS object; a list of one marginal model per location, output from migpd, $models item.")
  }
  
  
  
    co <- sapply(MODELS, coef)
    up <- sapply(MODELS, endPoint, verbose = FALSE)
    co <- rbind(mth, mqu, co, up)
    dimnames(co) <- list(c("Threshold", "P(X < threshold)", 
                           "sigma", "xi", "Upper end point"), 
                         names(MODELS))
    co[3, ] <- exp(co[3, ])
    co
}

revTransform_slim <- function (x, data_temp, qu, th = 0, sigma = 1, xi = 0, 
                               method = c("mixture", "empirical")){
  
  # Extension of **revTransform** from texmex, removing data from the arguments.
  # Reduces the number of outputs to lower overheads.
  # fetches DATA from the global environment
  
  if(!exists("DATA")){
    stop("Needs DATA object; one event per row, one location per column.")
  }
  
  method <- match.arg(method)
  n <- length(data_temp)
  probs <- (1:n)/(n + 1)
  # px <- vapply(x, p[[which.min(abs(x-p)))]]}, 0, p=probs) 
  # 
  # px <- as.integer(round(px * (1 + n)))
  
  px <- as.integer(pmax( pmin( round( x*(n+1) ), n), 1))
  
  res <- sort(data_temp)[px]
  if (method == "mixture") {
    i.x <- x >= qu
    i.r <- res > th
    i.rx <- i.x & i.r # apply(cbind(i.x, i.r), 1, all)
    if(is.na(i.rx[1])){
      print(str(i.rx))
      if(interactive()) browser()
    }
    if (sum(i.rx > 0)) {
      wh <- texmex:::u2gpd(x[i.rx], p = 1 - qu, th = th, sigma = sigma, 
                  xi = xi)
      rth <- res[i.rx]
      o <- order(rth)
      rth <- rth[o]

      
      rth[length(rth):(length(rth) - length(wh) + 1)] <- rev(sort(wh))
      rth <- rth[order(o)]
      res[i.rx] <- rth
    }
  }
  res[order(x)] <- sort(res)
  res
}


predict.mex_slim <- function(which, referenceMargin=NULL, marginfns,
                             constrain, coeffs_in, z_in, 
                             mth, mqu, pqu = .99, nsim = 1000, trace=10,
                             smoothZdistribution=FALSE, d, ...){
  
  # Extension of **predict.mex** from texmex, removing data from the arguments.
  # Reduces the number of outputs to lower overheads.
  # fetches DATA and TRANSFORMED from the global environment
  
  if(!exists("DATA")){
    stop("Needs DATA object; one event per row, one location per column.")
  }
  if(!exists("TRANSFORMED")){
    stop("Needs TRANSFORMED object; laplace-transformed output from mexTransform, $transformed item")
  }
    
    if(is.null(referenceMargin)){
      migpd <- list(transformed=TRANSFORMED)
    } else {
      migpd <- referenceMargin
    }
      marginfns <- marginfns
      constrain <- constrain
    
    ################################################################
    MakeThrowData <- function(dco,z,coxi,coxmi){
      ui <- runif(nsim , min = max(c(mqu[which], pqu)))
      y <-  marginfns$p2q(ui)
      distFun <- marginfns$q2p
      z <- z[!is.na(z[,1]), !is.na(z[1,])]
      z <- as.matrix(z[ sample( 1:( dim( z )[ 1 ] ), size=nsim, replace=TRUE ) ,])
      if(smoothZdistribution){
        z <- apply(z,2,function(x)x + rnorm(length(x),0,bw.nrd(x)))
      }
      ymi <- sapply( 1:( dim( z )[[ 2 ]] ) , makeYsubMinusI, z=z, v=dco , y=y )
      #print(range(ymi[,11]))
      xmi <- apply( ymi, 2, distFun )
      xmi[xmi > (1 - 1e-8)] <- (1 - 1e-8) # ! # FUDGE TO AVOID EXACTLY 1
      #print(range(xmi[,11]))
      xi <- texmex:::u2gpd( ui, p = 1 - mqu[which], th = mth[which],
                   sigma = coxi[1], xi = coxi[2] )
      
      for( i in 1:( dim(xmi)[[2]] ) ){
        #print(paste("I=", i))
        xmi[, i] <- revTransform_slim(
                                   xmi[, i], as.matrix(DATA[, -which])[, i],
                                   th = mth[-which][i],
                                   qu = mqu[-which][i],
                                   sigma=coxmi[1, i], xi=coxmi[2, i])
      }
      sim <- data.frame( xi , xmi)
      names( sim ) <- c( colnames( DATA )[ which ],
                         colnames( DATA )[ -which ])
      sim[,dim(sim)[2]+1] <- y > apply(ymi,1,max) # condlargest extra column
      sim
    }
    
    ################################################################
    makeYsubMinusI <- function( i, z, v , y ){
      v <- v[ , i ]
      z <- z[ , i ]
      if ( !is.na( v[ 1 ] ) ){
        if( v[ 1 ] < 10^(-5) & v[ 2 ] < 0 ){
          if( v[ 4 ] < 10^(-5 ) ) d <- 0
          else d <- v[ 4 ]
          a <- v[ 3 ] - d * log( y )
        }
        else a <- v[ 1 ] * y
      } # close if( !is.na...
      else a <- NA
      a + ( y^v[ 2 ] ) * z
    }
    bootRes <- NULL
    
    cox <- coef.migpd_slim(mth, mqu)[3:4, which]
    coxmi <- as.matrix(coef.migpd_slim(mth, mqu)[3:4, -which])
    
    sim <- MakeThrowData(dco=coeffs_in,
                         z=z_in,
                         coxi=cox,
                         coxmi=coxmi)
    CondLargest <- sim[,dim(sim)[2]]
    sim <- sim[,-(ncol(sim))]
    
    m <- 1 / ( 1 - pqu ) # Need to estimate pqu quantile
    zeta <- 1 - mqu[ which ] # Coles, page 81
    pth <- mth[ which ] + cox[ 1 ] / cox[ 2 ] * ( ( m*zeta )^cox[ 2 ] - 1 )
    # 
    # datafit <- list( #real = data.frame( data[, which], data[, -which] ),
    #               simulated = sim,
    #               CondLargest=CondLargest)
    # 
    #res <- list(datafit = datafit)
    as.matrix(sim[CondLargest, order(c(which, c(1:d)[-which]))])
    #oldClass( res ) <- "predict.mex" # A bit of a lie, but keeps things smooth.
    
    #res
}

