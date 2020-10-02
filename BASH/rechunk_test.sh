#!/bin/bash

echo "rechunking RCM01"

nccopy -c Northing/19,Easting/13,Time/4 ../Data/RCM01_205012_208011/dmflow_copy_RCM01_205012_208011.nc ../Data/RCM01_205012_208011/dmflow_rechunked_RCM01_205012_208011.nc 
