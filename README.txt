#AQUACAT - NEC07441
README

### Flood risk estimates combining the bbest of techniques from catastrophe modelling and future flood risk assessment tools.

### Project in association with Sayers and Partners

Last updated: 2021-02-01 Adam Griffin

## Folders:
>> RCMXX_YYYY12_ZZZZ11: 
        RCM ensemble member XX (between 1 and 15, doesn't include 2,3 or 14), for the period 1980-2010 or 2050-2080, with the year starting in December.

## Files:

>> hasData_primary.csv
    Table for converting netCDF row and column to Northing and Easting (GB National Grid). ALL tables with N rows correspond to this ordering (which is North to South, East to West)
    > contains:
        N rows: one location per row.
        4 columns: row, col, easting, northing.

>> hasData_Regions.csv
    Table for obtaining regions from row and column.
    > contains:
        N rows: one location per row.
        5 columns: row, col, HA_NUM (Hydrological Area Number), HA_NAME (Hydrological Area Name), REGION (region code)

    Regions:
    ANG - Anglia
    ESC - East Scotland
    NE - North East England
    NSC - North Scotland
    NW - North West England
    SE - South East England
    SEV - Severn Basin
    SSC - South Scotland
    SW - South West England
    THA - Thames Basin
    TRE - Trent Basin
    WAL - Wales

>> eventOBS_POT2_pc01_RCM**_****12_****11.nc
>> eventEC_POT2_pc01_RCM**_****12_****11.nc
    Legend:
    8 variables (excluding dimension variables):
        int row[loc]   (Contiguous storage)  
            units: Num
            _FillValue: -9999
            long_name: Row
        int col[loc]   (Contiguous storage)  
            units: Num
            _FillValue: -9999
            long_name: Column
        double northing[loc]   (Contiguous storage)  
            units: metres
            _FillValue: -9999
            long_name: Northing
        double easting[loc]   (Contiguous storage)  
            units: metres
            _FillValue: -9999
            long_name: Easting
        double flow[loc,event]   (Chunking: [4,250])  (Compression: level 2)
            units: cumecs
            _FillValue: -9999
            long_name: Peak Flow
        double dpe[loc,event]   (Chunking: [4,250])  (Compression: level 2)
            units: ProbOfExc
            _FillValue: -9999
            long_name: Daily Probability of Exceedance
        double ape[loc,event]   (Chunking: [4,250])  (Compression: level 2)
            units: ProbOfExc
            _FillValue: -9999
            long_name: Annual Probability of Exceedance
        int eventNo[event]   (Chunking: [599])  
            units: Matching OBS event
            _FillValue: -1

     2 dimensions:
        loc  Size:19914
            units: Num
            long_name: Location
        event  Size:599   *** is unlimited ***
            units: Num
            long_name: Event

    5 global attributes:
        RCM: RCM 07
        period: present
        event threshold: POT2
        area lower limit: pc01
        Method: OBS

>> eventOBS_region_ZZZ_RCM**_****12_****11.nc
>> eventEC_region_ZZZ_RCM**_****12_****11.nc
>> eventHT_region_ZZZ_RCM**_****12_****11.nc
    Legend as above except:
    _region_WAL: just gridpoints in region ZZZ, just events for which at least one point in ZZZ was affected.
    Regionalised data is just for POT2, pc01 observations.

    > both contain:
        N rows: one location per row - rows correspond to hasData_primary.csv (for national) or hasData_regions.csv (for regional).
        E columns: one column per event.


>> threshMat_RCM01_198012_2010.csv
    Threshold matrix, showing threshold per gridcell for each threshold POT5, POT2, POT1, Q2, Q5
    > contains:
        N rows: one location per row.
        5 columns: one threshold per column.

>> paramtableG_POT2_RCM01_198012_201011
    Generalised Pareto distribution parameters fitted to threshold exceedences (POT2) for each gridcell.
    > contains:
        N rows: one location per row.
        6 columns:  meanint - mean interval between events (measured in years). meanint of 1 equates to AMAX.
                    threshold - flow value of threshold
                    location - GPA location parameter, similar to threshold.
                    scale - GPA scale parameter > 0.
                    shape - GPA shape parameter.
                    threshquan - quantile of daily flow corresponding to threshold.

>> coefficients_NW_POT2_pc01_RCM**_****12_****11.csv
    Parameters from Heffernan-Tawn model - a,b,c,d,mu,sigma, for each pair of points in the region.







