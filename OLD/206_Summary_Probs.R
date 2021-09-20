library(dplyr)

files_in <- c(#"S:/Data/eventSumm_OBS_POT2_pc01_ALL.csv",
            #"S:/Data/eventSumm_EC_POT2_pc01_ALL.csv",
            #"S:/Data/eventSumm_OBS_region_NW_POT2_pc01_ALL.csv",
            #"S:/Data/eventSumm_EC_region_NW_POT2_pc01_ALL.csv",
            "S:/Data/eventSumm_HT_NW_POT2_pc01_ALL.csv")
files_out <- stringr::str_replace(files_in, "ALL", "SEASON_PROBS")

files_yr <- stringr::str_replace(files_in, "ALL", "EVENTS_PER_YEAR")

for(i in 1:length(files_in)){
df1 <- readr::read_csv(files_in[i])
df2 <- df1 %>% dplyr::count(rcm, period, season) %>%
  dplyr::group_by(rcm, period) %>% 
  dplyr::mutate(seasonprob=n/sum(n)) %>%
  tidyr::pivot_wider(id_cols=c("rcm","period"), names_from="season", 
                                  values_from="seasonprob")
readr::write_csv(df2, files_out[i])



df4 <- df1 %>% dplyr::mutate(startyr = ifelse(period=="present",1980,2050),
                             year = floor(eventDay/360) + startyr) %>%
  dplyr::count(rcm, period, year) %>%
  dplyr::group_by(rcm, period) %>%
  dplyr::summarise(events_per_year=round(mean(n)), .groups="keep")
readr::write_csv(df4, files_yr[i])

}