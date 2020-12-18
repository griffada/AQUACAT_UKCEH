AQUACAT - NEC07441
Flood risk estimates combining the bbest of techniques from catastrophe modelling and future flood risk assessment tools.

Project in association with Sayers and Partners

Files:

102_Threshold_Extract.R: Determine per-cell threshold for event extraction based on one output of 30-years.
# OUTPUTS: threshDayExcList.rda,
#          thresMat.rda,
#          threshGrid2.asc

103_Regional_Splitting: Script to split GB into hydrological regions based on sets of NRFA Hydrological Areas.
# OUTPUTS: hasData_Regions.csv

03_Event_Extract.R/.Rmd: Different sizes of widespread events and amounts of grid-cell threshold limits. 

104_Event_Extract.R: Using threshiold matrices from 102, event extraction from daily flow, no test for independence.
# OUTPUTS: eventLists***.Rda = eventLList, eventDayList

105_Event_Summary.R: Script for summarising events (pointwise event maxima), and calculating size and duration of extreme events.
# OUTPUTS: eventdf_***.csv

106_PoE_Estimation.R: Estimating probability of exceedence along time series using ecdf and plotting positions. Uses inputs from 104_Event_Extract and 105_Event_Summary.
# OUTPUTS: present_returnlevels***.csv

107_Empirical_Copula.R: The Empirical Copula functions for simulating new events from given data.
# OUTPUTS: NewEventPresentEC_***.csv

108_Event_Splitting.R: Splitting events by region defined in 103.
# OUTPUTS: regionalEvents_***.csv, eventdf_region_***.csv

109_HeffTawn_Modelling.R: Using Heffernan and Tawn model for spatial coherence to generate new events.
# OUTPUTS: NewEventHT_***.csv, coefficients.rds, zscores.rds, depStruct.rds

110_HT_PoEEstimation.R: Estimating probability of exceedence for Heffernan-Tawn events.
# OUTPUTS: HTraritydf_***.csv

111_SQL_Compilation.R: SQL database construction for HT/EC modelled events.
# OUTPUTS: eventdb_EC_***.sqlite, eventdb_OBS_***.sqlite, eventdb_HT_***.sqlite

199_RLplotting1.R: Function for plotting extracted values in terms of return period.
# OUTPUTS: .png plots


### DEV FILES

04_timeDF_compile.R: Summarises sizes of events and extract events for EC/HT.

04b_timeDF_compile.R: as in 04, but just for 2% spread at POT2 level.

05_PoEestimation.R: Estimating PoE along time series using ecdf and plotting positions.

06_EmpiricalCopula.R: Functions for EC and initial set of events.

06_EmpiricalTailProbabilities.Rmd: Trying to decide on a tail probability function.

07_HeffTawnModelling.R: functions for HT and initial set of events.

08_CoastalDataWrangling.R: take BODC and extract AMAX and full day data.

08a_CoastalDataWrangling.R: (NOT for use on Linux) take BODC and extract AMAX and full day data.
