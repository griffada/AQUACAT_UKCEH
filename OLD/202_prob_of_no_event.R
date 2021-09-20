rcm <- c("01","04","05","06","07","08","09","10","11","12","13","15")
period <- c("present","future")

df <- expand.grid(rcm, period)
colnames(df) <- c("rcm","period")
df$noevent <- 0
df$numevents <- 0


if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

for(i in seq_len(nrow(df))){
print(df[i,])
RCM <- df$rcm[i]
period <- df$period[i]
if(period == "present"){
  suffix <- "_198012_201011"
}else if (period == "future"){
  suffix <- "_205012_208011"
}
subfold <- paste0("RCM", RCM, suffix, "/")

thresh1 <- "POT2"
ws1 <- "pc05"
# Only using POT2 and 2% inundation minimums.
jT <- which(threshName==thresh1)
jW <- which(wsName == ws1)

load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))

eventz <- eventLList[[jT]][[jW]]

df$noevent[i] <- 1-sum(eventz)/ND

df$numevents[i] <- length(eventz)/30

df$q100[i] <- 1 - exp(-0.01*length(eventz)/30)

df$q100b[i] <- 1- (1-0.01)^(length(eventz)/30)

}
df2 <- df %>% group_by(period) %>% summarise(noeventav = mean(noevent))
df3 <- df %>% group_by(rcm) %>% summarise(diff = max(noevent)-min(noevent))
write.csv(df, "./prob_of_no_event.csv")



##### THE WRONG CALCULATION ###
1 - pbinom(ceiling(0.005*NH), NH, prob=2/360, lower.tail=FALSE)
# This assumes all the sites are independent


