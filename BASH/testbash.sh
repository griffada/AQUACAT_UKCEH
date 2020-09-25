#!/bin/bash

echo "01 future run of 102-103"

Rscript --vanilla 102_Threshold_Extract.R 01 present > pipe102_pr.err 2> pipe102_pr.out
Rscript --vanilla 103_Regional_Splitting.R 01 present > pipe103_pr.err 2> pipe103_pr.out

Rscript --vanilla 102_Threshold_Extract.R 01 future > pipe102_fu.err 2> pipe102_fu.out
Rscript --vanilla 103_Regional_Splitting.R 01 future > pipe103_fu.err 2> pipe103_fu.out