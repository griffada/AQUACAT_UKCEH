#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-06-25
#
# Trial SQL setup for event sets.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-06-15
#
#~~~~~~~~~~~~~~~~~~~~~~~

library("readr")
library("dplyr")
library("dbplyr")
library("RSQLite")
library("reshape2")


setwd("S:")
rl_db_file <- "./Data/present_return_levels2.sqlite"
rl_db <- src_sqlite(rl_db_file, create=TRUE)

PRES <- readr::read_csv("./Data/present2_returnlevels.csv")

load(paste0(wd,"/CodeABG/InterimData/eventLists03.RDa"))

generated_events <- readr::read_csv("./Data/NewEventPresentEC_POT2_pc05.csv")
colnames(generated_events) <- paste0("L", 1:19914)
rownames(generated_events) <- paste0("D", 1:nrow(generated_events))

colnames(PRES)[6] <- "startDay"

locationThresh <- PRES %>%
  dplyr::filter(eventNo == 1) %>%
  dplyr::select(loc, Easting, Northing, thresh)

eventDuration <- data.frame(eventNo     = 1:285,
                            startDay    = eventDayList[[2]][[2]],
                            eventLength = eventLList[[2]][[2]])

presmin <- PRES %>%
  dplyr::select(eventNo, loc, val, gpp, ecdf, gev)

copy_to(rl_db, locationThresh, temporary=FALSE, overwrite=TRUE)
copy_to(rl_db, eventDuration,  temporary=FALSE, overwrite=TRUE)
copy_to(rl_db, presmin,        temporary=FALSE, overwrite=TRUE)

rl_db

eventVloc <- acast(presmin, eventNo ~ loc, value.var="gpp")

colnames(eventVloc) <- paste0("L",1:19914)
rownames(eventVloc) <- paste0("E",1:285)

eventMags <- data.frame(t(eventVloc))

copy_to(rl_db,
        eventMags,
        temporary=FALSE,
        overwrite=TRUE)

copy_to(rl_db,
        generated_events,
        temporary=FALSE,
        overwrite=TRUE)



# my_db <- DBI::dbConnect(RSQLite::SQLite(), "./Data/present_return_levels.sqlite")
# src_dbi(my_db)
# dbListTables(my_db)
# dbListFields(my_db, "eventDuration")
# dbListFields(my_db, "locationThresh")
# dbListFields(my_db, "presmin")

##### TUTORIAL STUFF #####----------------------------------------------------
if(FALSE){
  setwd("~/AQUACAT_C")
  
  my_db_file <- "./bodc_sea_level.sqlite"
  my_db <- src_sqlite(my_db_file, create=TRUE)
  my_db
  
  SH <- readr::read_csv("./Sheerness_19861024_20180430.csv",
                        col_types="cccddcccccdc")
  
  copy_to(my_db, SH, temporary=FALSE, overwrite=TRUE)
  
  my_db
  
  my_sea_level <- DBI::dbConnect(RSQLite::SQLite(), "./bodc_sea_level.sqlite")
  src_dbi(my_sea_level)
  DBI::dbDisconnect(my_sea_level)
  
  DBI::dbDisconnect(my_db)
}