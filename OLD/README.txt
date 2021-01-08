AQUACAT - NEC07441
Flood risk estimates combining the bbest of techniques from catastrophe modelling and future flood risk assessment tools.

Project in association with Sayers and Partners

Files:
00_prelimPOE.R - Initial investigation of 2-year slice of time from VB: methods of identifying spatial events.

01_texmex_scratch.Rmd - R markdown outlining texmex functions, and basic
outline of data from G2G outputs.

02_Threshold_Extract.R: Determine per-cell threshold for event extraction based on one output of 30-years.

02b_Threshold_Plotting.R: plotting the maps of thresholds from 02

02c_Threshold_Plotting2.R: different version of plotting as in 02b

03_Event_Extract.R/.Rmd: Different sizes of widespread events and amounts of gric-cell threshold limits. 

04_timeDF_compile.R: Summarises sizes of events and extract events for EC/HT.

04b_timeDF_compile.R: as in 04, but just for 2% spread at POT2 level.

05_PoEestimation.R: Estimating PoE along time series using ecdf and plotting positions.

06_EmpiricalCopula.R: Functions for EC and initial set of events.

06_EmpiricalTailProbabilities.Rmd: Trying to decide on a tail probability function.

07_HeffTawnModelling.R: functions for HT and initial set of events.

08_CoastalDataWrangling.R: take BODC and extract AMAX and full day data.

08a_CoastalDataWrangling.R: (NOT for use on Linux) take BODC and extract AMAX and full day data.
