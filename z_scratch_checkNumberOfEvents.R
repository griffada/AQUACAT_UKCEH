library(ncdf4)
data_wd <- "S:/Data/"
REG <- "NW"
ws1 <- "pc01"
thresh1 <- "POT2"
FF_flag <- FALSE
RCMS <- c("01","04","05","06","07","08","09","10","11","12","13","15")
for(RCM in RCMS){
for(j in 1:2){
period <- c("present","future")[j]
suffix <- c("_198012_201011", "_205012_208011")[j]
if(FF_flag){
subfold <- paste0("RCM", RCM, suffix, "_FF/")
}else{
subfold <- paste0("RCM", RCM, suffix, "/")
}
savepath <- paste0(data_wd, subfold, REG, "/eventHT_",REG,"_",
thresh1,"_", ws1, "_RCM", RCM, suffix, ".nc")
try({
ht_events <- nc_open(savepath, write=T)
evNo <- ncvar_get(ht_events, "eventNo")
print(savepath)
print(range(evNo, na.rm=T))
})
}
}
