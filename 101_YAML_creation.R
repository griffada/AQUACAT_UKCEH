#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2021-08-12
#
# Construction of setting YAML file, which tracks steps of the pipeline
# successfully completed, along with number of EC and HT simulations.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-08-12
#
# OUTPUTS: settings.yaml
#
#~~~~~~~~~~~~~~~~~~~~~~~

if(interactive()){commandArgs <- function(...){c("04","future", "NW")}}

#### SETUP ####----------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) %in% c("Cen","Fed")){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

L <- list(Msims=25000,
  nSampleHT=400,
  region='NW',
  RCM=RCM,
  period=period,
  thresh="POT2",
  wsname="pc01",
  paramtable=FALSE,
  eventlist=FALSE,
  OBSflow=FALSE,
  OBSdpe=FALSE,
  OBSape=FALSE,
  OBSsumm=FALSE,
  EC2flow=FALSE,
  EC2dpe=FALSE,
  EC2ape=FALSE,
  EC2summ=FALSE,
  HTsplit=FALSE,
  HTstructure=FALSE,
  HTflow=FALSE,
  HTdpe=FALSE,
  HTape=FALSE,
  HTsumm=FALSE
)
write_yaml(L, settingspath)

print("YAML created. 101 done.")