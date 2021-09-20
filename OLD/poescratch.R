h <- 81

thr <- threshMat[h,jT]
meanInt <- partable$meanint[h]
scaleH <- partable$scale[h]
shapeH <- partable$shape[h]

vals <- ncvar_get(ncin_pres, "dmflow",
                       start=c(rn$row[h], rn$col[h], 1),
                       count=c(1, 1, -1))

soo <- sort(obs_events[h,], decreasing=T)[1:10]

ep1 <- (ilaprosUtils::extractPeaks(vecObs=vals, mintimeDiff=7) == 1) &
            (vals > thr)
pgl <- pargpa(lmoms(vals[ep1]), xi=thr)
ub <- thr + (pgl$para[2]/pgl$para[3])

ecd <- ecdf(c(-1,vals,1e8))
valsdpe <- 1 - ecd(unlist(obs_events[h,]))
valsape <- 1 - exp(-valsdpe*360)


XX <- 1/(1 - cdfgpa(as.numeric(soo), pargpa(lmoms(vals[ep1]), xi=thr)))
XY <- 1/(1 - pevd(as.numeric(soo), threshold=thr, scale=scaleH, shape=shapeH, type='GP'))
XZ <- meanInt*(1:10 - 0.44)/(sum(ep1) + 0.12)

XW <- do.call(cbind, list(soo=unlist(soo), lmomms=XX,mli=XY, gri=1/XZ, vals=1/(sort(valsape)[1:10])))
              