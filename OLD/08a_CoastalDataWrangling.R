#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-06-17
#
# Data wrangling for datasets from BODC
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-06-19
#
# DONT TRY THIS ON LINUX, THE UNZIP IS NOT RELIABLE AND TAKES TOO LONG.
#
#~~~~~~~~~~~~~~~~~~~~~~~

if (substr(osVersion,1,3) == "Win") {
  wd <- "S:/"
  wd_id <- "S:/Data/"
  wd_cd <- "S:/CodeABG/"
  
} else {
  wd <- "/prj/aquacat/"
  wd_id <- "/prj/aquacat/Data/"
  wd_cd <- "/prj/aquacat/CodeABG/"
}

##### SETUP #####--------------------------------------------------

library(dplyr)

library(lubridate)

# Get the .csvs from the zip folder
UList <- unzip(zipfile=paste0(wd_id,"BODC_data.zip"), list=T)$Name
UList <- UList[grep(".csv", UList)]



locations <- sapply(strsplit(UList,"_"),
                   function(x){paste(x[1:(length(x)-2)], collapse="_")})
Ulocations <- unique(locations)

AMAXlist <- vector("list", length(Ulocations))

AMAXtemp <- data.frame()


allAmaxDays <- c()

#PT <- proc.time()
dateFull <- c()
for(i in 1:length(Ulocations)){  # for each location
  
  print(Ulocations[i])
  loclist <- UList[locations == Ulocations[i]]
  
  for(k in 1:length(loclist)){  # for each file for that location
    
    
    # get the file
    U <- unzip(zipfile=paste0(wd_id,"BODC_data.zip"), files=loclist[k])
    U2 <- readr::read_csv(loclist[k],
                          col_types="cccddcccccdc")
    colnames(U2) <- sub(" ", "_",colnames(U2))
    
    
    #get AMAX, based on DMAX
    UV <- U2 %>%
      mutate(Date2=ymd_hms(Date)) %>%
      mutate(month = month(Date2),
             year=year(Date2), date3=as_date(Date2)) %>%
      group_by(date3, year, Site_Name) %>%
      summarise(DMAX = max(Data_value)) %>%
      group_by(year, Site_Name) %>%
      summarise(maxm = max(DMAX),
                maxd = date3[which(DMAX == max(DMAX))[1]])
    
    allAmaxDays <- c(allAmaxDays, UV$maxd)
    
    if(k==1){
      AMAXtemp <- UV
    }else{
      AMAXtemp <- rbind(AMAXtemp, UV)
    }
    
    dateFull <- c(dateFull, UV$maxd)
    #AMAXtemp <- ifelse(k == 1, UV, rbind(AMAXtemp, UV))
    #print(proc.time() - PT)
    
    file.remove(loclist[k])  # don't keep hold of all the GBs of files.
  }
  readr::write_csv(x=AMAXtemp,
                  path=paste0(wd_id,"CoastalSumm/",Ulocations[i],".csv"))
}

Udate <- unique(dateFull)
saveRDS(Udate, paste0(wd_id,"Udate.RDS"))


##### STEP 2 #####
# Go back and get the other days for the extreme events.
