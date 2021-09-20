library(fevd)
library(texmex)
library(extRemes)
library(parallel)
library(foreach)
library(readr)

OBS <- t(as.data.frame(readr::read_csv("./Data/RCM01_198012_201011/eventflow_OBS_POT2_pc01_RCM01_198012_201011.csv"))[,-(1:4)])
APE <- t(as.data.frame(readr::read_csv("./Data/RCM01_198012_201011/NW/eventape_OBS_region_NW_RCM01_198012_201011.csv"))[,-(1:4)])
DPE <- t(as.data.frame(readr::read_csv("./Data/RCM01_198012_201011/eventdpe_OBS_POT2_pc01_RCM01_198012_201011.csv"))[,-(1:4)])

partable <-  as.data.frame(readr::read_csv(paste0("./Data/RCM01_198012_201011/",
                           "paramtableG_POT2_RCM01_198012_201011.csv"),
                            col_types=cols(.default= col_double())))

d <- 10
n <- 606
RankOBS <- apply(OBS[1:n, 1:d], 2, rank, ties.method = "random")
M <- 1000
dm <- 5

system.time({
  cl <- parallel::makeCluster(3, outfile = "out1.out")
  doParallel::registerDoParallel(cl)
  Vnew <- foreach(m = 1:M , .combine = rbind) %dopar% {
    print(m)
    U <- matrix(runif(d * n), nrow = n, ncol = d)
    V <- apply(U, 2, sort)  # This is the slow bit.
    w <- sample.int(n, dm)
    Vnew0 <- matrix(NA, nrow=dm, ncol=d)
    for (i in 1:d) {
      Vnew0[,i] <- V[RankOBS[w, i], i]
    }
    Vnew0
  }
  stopCluster(cl)
})
pp <- 2/360

OBS0 <- OBS[1:n,1:d]

qnewB <- sapply(1:d, function(i){quantile(OBS0[,i], probs=Vnew[,i])})

vnewC <- Vnew[, 1][Vnew[,1]>(1-pp)]
f <- function(v, pp=0.95, params, OBS=OBS0){
 qout <- rep(0, length(v))
 up <- v > (1-pp)
 vsub <- v[up]
 qout[up] <- qglo(1 - (1-vsub)/pp, loc=params[1], scale=params[2], sh=params[3])
 vsub <- v[!up]
 qout[!up] <- quantile(OBS, probs=vsub)
 qout
}

qnewC <- sapply(1:d, function(i){f(Vnew[,i], pp=2/360, params=partable[i,3:5], OBS=OBS0[,i])})

qnewC <- qglo(1 - ((1-vnewC)/pp), loc=31.2, scale=2.75, sh=-0.523)
qnewD <- qnewB
qnewD[Vnew[,1]>(1-pp)] <- qnewC


#flow_glo <- function(v, params, OBS = OBS0) {
v=Vnew[,1]
params=partable[1,]
OBS=OBS0[,1]

  thresh0 <- params$threshold
  logsp <- logspline(OBS, lbound=0)
  threshquan <- plogspline(thresh0, logsp)
  #threshquan <- ecdf(OBS)(thresh0)
  qout <- rep(0, length(v))
  up <- v > threshquan
  vsub <- v[!up]
  qout[!up] <- quantile(OBS, probs = vsub)
  range(qlogspline(vsub, logspline(OBS0, lbound=0)))
  vsub <- v[up]
  dp <- pglo(thresh0, loc = params$loc,
    scale = params$sca,
    sh = params$shape
  )
  
  range(1/((1-vsub)*(1-dp)/(1-threshquan)))
  
  qout[up] <- qglo(
    1 - (((1-vsub)*(1-dp))/(1-threshquan)),
    loc = params$loc,
    scale = params$sca,
    sh = params$shape
  )
  range(qout)
#}
debug(flow_glo)
flow_glo(Vnew[,1], partable[1,], OBS0[,1])

dpeA <- 1-ecdf(tSlice)(qout)






