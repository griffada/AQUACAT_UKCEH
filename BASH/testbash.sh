#!/bin/bash

# This bash script runs the bare essentials to get the Regional split and the
# point-by-point thresholds for RCM01, present and future.

echo "102 01 present"
Rscript --vanilla 102_Threshold_Extract.R 01 present > pipe102_pr.out 2> pipe102_pr.out
tail pipe102_pr.out
echo "103 01 present"
Rscript --vanilla 103_Regional_Splitting.R > pipe103.out 2> pipe103.out 
tail pipe103.out
echo "102 01 future"
Rscript --vanilla 102_Threshold_Extract.R 01 future > pipe102_fu.out 2> pipe102_fu.out 
tail pipe102_fu.out

echo "Complete"