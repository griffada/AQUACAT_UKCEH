print("looking for packages")

print(R.version)
print(Sys.info())

PA <- c("dplyr", "texmex", "ncdf4", "raster", "fields", "rgdal", "readr", "extRemes", "reshape2")

print(PA[PA %in% installed.packages()$Package])

install.packages(pkgs = PA[!(PA %in% installed.packages()$Package)],
                 lib = "/home/users/adagri/R/x86_64-conda_cos6-linux-gnu-library/3.6",
                 repos = "https://www.stats.bris.ac.uk/R/")
require(dplyr)
require(texmex)
require(ncdf4)
require(raster)
require(fields)
require(rgdal)
require(readr)

print(sessionInfo())