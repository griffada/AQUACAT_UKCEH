library(dplyr)
st <- 19801200
en <- 20101200


# L <- unzip(zipfile = "/prj/aquacat/Other Data/Misc/BODC_data.zip", list = TRUE)
# LS <- sapply(L[,1], function(x){strsplit(x, split="[\\._]+")[[1]][c(2,3)]})
# LS[,1] <- as.numeric(LS[,1])
# LS[,2] <- as.numeric(LS[,2])
# LS <- as.data.frame(t(LS), stringsAsFactors=F)
# st <- 19801200
# en <- 20101200
# LS$within <- (LS[,1] < st & LS[,2] > st)|(LS[,1] > st & LS[,2] < en)|(LS[,1] < en & LS[,2] > en)
# which.max(LS$within)
# LS$within[is.na(LS$within)] <- FALSE
# 
# 
# for(l in L[which(LS$within),1]){
#   print(l)
#   unzip(zipfile = "/prj/aquacat/Other Data/Misc/BODC_data.zip",
#       files= l,
#       exdir="/prj/aquacat/Other Data/SeaLevels")
# }
#   unzip(zipfile = "/prj/aquacat/Other Data/Misc/BODC_data.zip",
#       files= L[which(LS$within)[2],1],
#       exdir="/prj/aquacat/Other Data/SeaLevels")


L <- dir("S:/Other Data/SeaLevels")
LS <- sapply(L, function(x){strsplit(x, split="_")[[1]][1]})
U <- unique(LS)
for(l in L){

 df <- as.data.frame(data.table::fread(paste0("/prj/aquacat/Other Data/SeaLevels/",l)))
 df <- df[,c(9,11)]
}


for(u in U){
 print(u)
 Lsub <- L[LS==u]
 dflist <- vector("list",length(Lsub))
 for(i in 1:length(Lsub)){
   dflist[[i]] <- as.data.frame(
     data.table::fread(paste0("S:/Other Data/SeaLevels/",Lsub[i])))[,c(9,11)]

 }
  df <- do.call(rbind,dflist)
  df$day <- as.Date(df[,1])
  df <- df[((df$day > as.Date("1980-11-30")) & (df$day < as.Date("2010-12-01"))),]
  df <- df %>% group_by(day) %>% summarise(max(V11))
  print(dim(df))
  write.csv(df, paste0("/prj/aquacat/Other Data/SeaLevels/",u,".csv"))
}

for(u in U){
 print(u)
 Lsub <- L[LS==u]
 dftemp <- as.data.frame(
     data.table::fread(paste0("S:/Other Data/SeaLevels/",Lsub[1])))

info <- dftemp[1,c(1,4,5)]
  df <- do.call(rbind,dflist)
  df$day <- as.Date(df[,1])
  df <- df[((df$day > as.Date("1980-11-30")) & (df$day < as.Date("2010-12-01"))),]
  df <- df %>% group_by(day) %>% summarise(max(V11))
  print(dim(df))
  write.csv(df, paste0("/prj/aquacat/Other Data/SeaLevels/",u,".csv"))
}