# devtools::install_github("https://github.com/Cdevenish/hardRain")

# load analysis tools
library(warbleR)
library(suncalc)
library(weathercan)
library(naniar) # explore missing data
library(sf)
library(mapview)
library(hardRain)
library(tuneR)
library(parallel)
library(progress)
library(furrr)

# load wrangling tools
library(tidyverse)
library(lubridate)
library(tidylog)
library(readxl)
library(patchwork) # combining plots
library(scales)
library(ggrepel) # text labels
library(GGally) # pair plots
library(skimr)
library(tictoc)

# set options
options("scipen" = 100) # prevent display with scientific notation
set.seed(123) # set random number generator seed
theme_set(theme_bw())

# set working directory
setwd("C:/Users/Jonathan/Dropbox/SS/EWPW")



#0. Wrangle recording metadata ------------------------------------------------

# note: deleted three 0 KB files: 
# SM12011_20170718_051200.wav, SM13520_20190512_073000.wav, & SM13520_20190516_101639.wav

# recording length (duration) (takes a couple minutes to run)
df_durations <- 
  # list directories
  c(list.dirs("D:/GRASSLAND2017/")[-1], # -1 removes root directory
    list.dirs("D:/PEPNWA_2019/")[-1]) %>% 
  # iterate over directories to read durations
  map_dfr(~{warbleR::wav_dur(path = .x)}) %>% 
  # tidy up
  rename(Recording = sound.files, Duration_sec = duration) %>%
  mutate(Duration_min = round(Duration_sec/60, 0))

table(df_durations$Duration_min)

save(df_durations, file = "Recording list with durations.Rdata")
# load("Recording list with durations.Rdata")

# make list of all recordings
df_recordings <- df_durations %>% select(Recording)


# site names and coordinates
df_locations <- 
  read_excel("GrassCommAllLocations2016-19.xlsx", sheet = "AllSites") %>% 
  mutate(SongMeter = paste0("SM", SongMeter))


# combine recording metadata for joining
df_recording_meta <-
  df_recordings %>% 
  # extract metadata from recording name
  separate(Recording, into = c("SongMeter", "Date", "Time_rec", "Filetype"), remove = FALSE) %>%
  select(-Filetype) %>% 
  # wrangle dates & times
  mutate(Date = ymd(Date),
         Year = year(Date),
         Yday = yday(Date),
         Time_rec = ymd_hms(paste(Date, Time_rec)),
         .after = Date) %>% 
  mutate(Time_rec = force_tz(Time_rec, "America/Toronto"))



#1. Read in all validated files and bind together -----------------------------

# list validated results files
detect_file_list <- list.files(path = "Processing_finished/", full.names = TRUE, pattern = "_validated.txt")

# read and combine
df_detect_raw <- map_dfr(detect_file_list,
                         read_tsv,
                         col_names = c("File_name", "Time_offset", "Duration", "Level", "Quality", "Score", "Recognizer", "Comments"),
                         col_types = c("ccddddcc"))

# quick peek at data
skimr::skim(df_detect_raw)
sort(unique(df_detect_raw$Comments))



#2. Tidy detection data -------------------------------------------------------

df_detect_clean <- 
  df_detect_raw %>% 
  # clean up file names
  mutate(File_name = str_replace(File_name, "D:\\\\", "")) %>% 
  mutate(File_name = str_replace(File_name, "E:\\\\", "")) %>% 
  mutate(File_name = str_replace(File_name, "/Volumes/Backup Plus/", "")) %>% 
  mutate(File_name = str_replace_all(File_name, "\\\\", "/")) %>% 
  # extract metadata from file path
  separate(File_name, into = c("Project", "SongMeter", "Recording"), sep = "/") %>% 
  # extract metadata from recording name
  separate(Recording, into = c("SongMeter2", "Date", "Time_rec", "Filetype"), remove = FALSE) %>% 
  # wrangle dates & times
  mutate(Date = ymd(Date),
         Year = year(Date),
         Yday = yday(Date),
         Time_rec = ymd_hms(paste(Date, Time_rec)),
         .after = Date) %>% 
  mutate(Time_rec = force_tz(Time_rec, "America/Toronto")) %>%
  mutate(Time_voc = Time_rec + hms(Time_offset), .after = Time_rec) %>%
  # calculate minute interval of each detection
  mutate(Rec_minute = (minute(hms(Time_offset)) + 1), .after = Time_offset) %>% 
  # tidy up columns
  select(-SongMeter2, -Filetype) %>% 
  # tidy up SongMeter IDs
  mutate(SongMeter = if_else(Year == 2019, paste0("SM", SongMeter), SongMeter)) %>% 
  # tidy up comments
  mutate(Comments_clean = 
           case_when(Comments == "?" ~ NA_character_,
                     Comments == "11" ~ "1",
                     Comments == "2 here" ~ "2",
                     Comments == "3?" ~ "3",
                     Comments == "Minimum 2" ~ "2",
                     Comments == "Minimum 3" ~ "3",
                     Comments == "na" ~ NA_character_,
                     Comments == "Na" ~ NA_character_,
                     TRUE ~ Comments)) %>% 
  mutate(Abundance = if_else(is.na(Comments_clean), 0, as.numeric(Comments_clean)),
         Detection = if_else(Abundance > 0, 1, 0))


table(df_detect_clean$Rec_minute)



#3. occupancy per every 1-minute recording interval ----

# complete list of 1-minute recording intervals
df_1mins <- 
  df_durations %>% 
  uncount(Duration_min, .id = "Rec_minute") %>% 
  select(-Duration_sec)

# detections in minute intervals
df_1min_detects <- 
  df_detect_clean %>% 
  group_by(Recording, Rec_minute) %>% 
  summarise(Occupied = max(Detection)) %>% 
  filter(Occupied == 1)

table(df_1min_detects$Rec_minute)


#3b. number of detections per 1 min
df_1min_n_detects <- 
  df_detect_clean %>% 
  filter(Abundance >= 1) %>% 
  group_by(Recording, Rec_minute) %>% 
  summarise(Num_detects = n())


#3c. abundance per recording
df_abundance_recording <- 
  df_detect_clean %>% 
  filter(Abundance >= 1) %>% 
  group_by(Recording) %>% 
  summarise(Abundance_recording = max(Abundance))

hist(df_abundance_recording$Abundance_recording)


#3d. abundance per site
df_abundance_site <- 
  df_detect_clean %>% 
  filter(Abundance >= 1) %>% 
  group_by(Project, SongMeter, Year) %>% 
  summarise(Abundance_site = max(Abundance))

hist(df_abundance_recording$Abundance_recording)


# join summaries of detections with complete interval list
df_1min_occ_abd <- 
  df_1mins %>% 
  left_join(df_1min_detects) %>% 
  left_join(df_1min_n_detects) %>% 
  left_join(df_abundance_recording) %>% 
  left_join(df_recording_meta, by = "Recording") %>% 
  left_join(df_locations, by = c("SongMeter", "Year")) %>% 
  left_join(df_abundance_site) %>% 
  select(Recording, Rec_minute, Occupied, Num_detects, Abundance_recording, Abundance_site) %>% 
  replace_na(list(Occupied           = 0, 
                  Num_detects          = 0, 
                  Abundance_recording = 0, 
                  Abundance_site     = 0)) 
  

table(df_1min_occupancy$Rec_minute, df_1min_occupancy$Occupied)

df_1min_occupancy %>% 
  filter(Rec_minute <= 15) %>% 
ggplot() + 
  geom_histogram(aes(x = Rec_minute), bins = 29) +
  scale_x_continuous(breaks = c(5, 10, 15))


#5. Calculate lunar & solar covariates for each day at each location ----------

# list of every recording minute and location ---
df_covar_prep <- 
  df_1mins %>% 
  left_join(df_recording_meta, by = "Recording") %>% 
  left_join(df_locations, by = c("SongMeter", "Year")) %>% 
  mutate(Time_rec_1min = Time_rec + minutes(Rec_minute - 1),
         Hour_rec_1min = hour(Time_rec_1min)) %>% 
  mutate(Date_EWPW = if_else(Hour_rec_1min < 12, Date - days(1), Date), .after = Time_rec)

# sun
df_sun_times <- 
  df_covar_prep %>% 
  distinct(SongMeter, Date, Latitude, Longitude) %>% 
  rename(date = Date, lat = Latitude, lon = Longitude) %>% 
  getSunlightTimes(data = ., tz = "America/Toronto") %>%
  rename(Date = date, Latitude = lat, Longitude = lon)


# moon

# list of every date x every location
location_x_dates <- 
  df_covar_prep %>% 
  group_by(SiteName) %>% 
  summarise(first = min(Date),
            last  = max(Date)) %>% 
  mutate(start = first - days(1),
         end   = last  + days(1)) %>% 
  group_by(SiteName) %>% 
  complete(start = seq.Date(start, end, by = "day")) %>% 
  select(SiteName, Date = start) %>% 
  left_join(df_locations)

df_moon_times <- 
  location_x_dates %>% 
  distinct(SongMeter, Date, Latitude, Longitude) %>% 
  rename(date = Date, lat = Latitude, lon = Longitude) %>%
  getMoonTimes(data = ., tz = "America/Toronto", keep = c("rise", "set")) %>% 
  mutate(date = date - days(1)) %>% # correction for date offset error in suncalc
  rename(Date = date, Latitude = lat, Longitude = lon,
         moonrise = rise, moonset = set)


df_moon_position <- 
  df_covar_prep %>% 
  distinct(SongMeter, Time_rec_1min, Latitude, Longitude) %>% 
  rename(date = Time_rec_1min, lat = Latitude, lon = Longitude) %>%
  getMoonPosition(data = ., keep = c("altitude", "azimuth")) %>% 
  mutate(moon_altitude_deg = (altitude * 180) / pi,
         moon_azimuth_deg  = (pi + azimuth) * 180 / pi %% 360,
         moonUp = if_else(moon_altitude_deg >= 0, "Y", "N")) %>% 
  rename(Time_rec_1min = date, Latitude = lat, Longitude = lon,
         moon_altitude = altitude, moon_azimuth = azimuth)


df_moon_phase <- 
  df_covar_prep %>% 
  distinct(Date) %>% 
  rename(date = Date) %>%
  pull(date) %>% 
  getMoonIllumination(keep = c("fraction")) %>% 
  rename(Date = date, moon_fraction = fraction)


# all astro together
df_sunmoon <- 
  df_covar_prep %>% 
  left_join(df_sun_times, by = c("Date", "Latitude", "Longitude")) %>% 
  left_join(df_moon_times, by = c("Date", "Latitude", "Longitude")) %>% 
  left_join(df_moon_phase, by = "Date") %>% 
  left_join(df_moon_position, by = c("Latitude", "Longitude", "Time_rec_1min"))



#6. Pull nearest (best) weather covariates for each location ------------------

# find nearby stations
stn200_all <- stations_search(coords = c(44.4, -76.75), dist = 200) %>% 
  filter(start <= 2017, end >= 2019, !(interval == "month"))

# only stations with *both* daily and hourly data
stn_dayhour <- 
  table(stn200_all$station_name, stn200_all$interval) %>%
  as.data.frame() %>% 
  mutate(station_name = as.character(Var1)) %>% 
  pivot_wider(names_from = Var2, values_from = Freq) %>% 
  filter(day == 1 & hour == 1)

stn200 <- 
  stn200_all %>% 
  filter(station_name %in% c(stn_dayhour$station_name))
  

# map locations
stn200.sf <- st_as_sf(stn200, coords = c("lon", "lat"), crs = 4326)
SM.sf <- st_as_sf(df_locations, coords = c("Longitude", "Latitude"), crs = 4326)

wx_map <- ggplot() +
  geom_sf(data = SM.sf, aes(colour = factor(Year)), size = 3) +
  geom_sf(data = stn200.sf) +
  geom_text_repel(data = stn200, aes(x = lon, y = lat, label = station_name), size = 2) +
  facet_grid(~interval) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

wx_map

ggsave("Weather station proximity.pdf", 
       plot = wx_map, width = 11, height = 8.5, units = "in", dpi = 220, device = cairo_pdf)

mapview(stn200.sf, co.regions = "red") +
mapview(SM.sf, zcol = "Year")

# find nearest weather station to each SongMeter
df_wx_nearest <- st_join(SM.sf, stn200.sf, join = st_nearest_feature)

table(df_wx_nearest$station_name)

wx_nearest_IDs <- as.character(sort(unique(df_wx_nearest$station_id)))

# find date ranges
df_recording_meta %>% 
  group_by(Year) %>% 
  summarise(start = min(Date), end = max(Date))

# daily data ---
wx_daily_2017 <- weather_dl(station_ids = wx_nearest_IDs,
                            start = as_date("2017-04-15"),
                            end = as_date("2017-08-03"),
                            interval = "day")

wx_daily_2019 <- weather_dl(station_ids = wx_nearest_IDs,
                            start = as_date("2019-05-06"),
                            end = as_date("2019-07-17"),
                            interval = "day")

wx_daily <- 
  bind_rows(wx_daily_2017, wx_daily_2019)


# hourly data ---
wx_hourly_2017 <- weather_dl(station_ids = wx_nearest_IDs,
                             start = as_date("2017-04-15"),
                             end = as_date("2017-08-03"),
                             interval = "hour")

wx_hourly_2019 <- weather_dl(station_ids = wx_nearest_IDs,
                             start = as_date("2019-05-06"),
                             end = as_date("2019-07-17"),
                             interval = "hour")

wx_hourly <- 
  bind_rows(wx_hourly_2017, wx_hourly_2019) %>% 
  mutate(time = force_tz(time, "Etc/GMT+5"))


# check for missing data ---
p_wx_na_h <- 
  wx_hourly %>% 
  select(station_name, weather:wind_spd) %>% 
  select(! contains("flag")) %>% 
  naniar::gg_miss_var(show_pct = TRUE, facet = station_name) +
  labs(title = "Hourly data", y = NULL)

p_wx_na_d <- 
  wx_daily %>% 
  select(station_name, cool_deg_days:total_snow) %>% 
  select(! contains("flag")) %>% 
  naniar::gg_miss_var(show_pct = TRUE, facet = station_name)+
  labs(title = "Daily data")
  
p_wx_na <- p_wx_na_h + p_wx_na_d
p_wx_na

ggsave("Weather station data availability.pdf", 
       plot = p_wx_na, width = 11, height = 8.5, units = "in", dpi = 220, device = cairo_pdf)


# Point Petre is closest, but missing a lot of precip data
# let's see how similar Trenton and Pt Pet are
wx_compare_daily <- function(VAR) {
  wx_daily %>% 
    select(date, station_name, all_of(VAR)) %>% 
    pivot_wider(names_from = station_name,
                values_from = VAR) %>% 
    ggpairs(columns = 2:7)+
  labs(title = paste(VAR, "(daily)"))
}

wx_compare_daily("total_precip")
wx_compare_daily("max_temp")
wx_compare_daily("min_temp")


# make table of which wx station to use for which SM site
NE_sites <- c("Limerick Rd", "Weedmark Rd", "Burnt Lands", "Panmvre Rd")

wx_stn_index <- 
  df_locations %>%
  mutate(station_name = if_else(SiteName %in% NE_sites, "OTTAWA CDA RCS", "TRENTON A"),
         station_id   = if_else(SiteName %in% NE_sites, "30578", "5126"))


# join hourly weather data together
df_weather_hourly_interp_trenton <- 
  df_covar_prep %>% 
  left_join(wx_stn_index %>% select(SiteName, station_name, station_id), by = "SiteName") %>% 
  #left_join(wx_daily, by = c("Date" = "date", "station_name", "station_id")) %>%
  mutate(time = Time_rec_1min) %>%
  filter(station_name == "TRENTON A") %>% 
  weather_interp(weather = wx_hourly %>% filter(station_name == "TRENTON A"))

df_weather_hourly_interp_ottawa <-
  df_covar_prep %>% 
  left_join(wx_stn_index %>% select(SiteName, station_name, station_id), by = "SiteName") %>% 
  #left_join(wx_daily, by = c("Date" = "date", "station_name", "station_id")) %>%
  mutate(time = Time_rec_1min) %>%
  filter(station_name == "OTTAWA CDA RCS") %>% 
  weather_interp(weather = wx_hourly %>% filter(station_name == "OTTAWA CDA RCS"))
  
df_weather_hourly_interp <- 
  bind_rows(df_weather_hourly_interp_trenton, df_weather_hourly_interp_ottawa) 


# join hourly and daily data
df_weather <-  
  df_weather_hourly_interp %>% 
  left_join(wx_daily, by = c("Date" = "date", "station_name", "station_id"))
  


#7. Run Hard Rain for each minute of each recording ----

####
# takes about a day to scan everything
# just load the results instead:
####
load("hardRain results.Rdata")
####

# create list of files
file_list <- list.files(path = "D:", recursive = TRUE, full.names = TRUE, pattern = ".wav")
file_list_trunc <- list.files(path = "D:", recursive = TRUE, full.names = TRUE, pattern = "_trunc60.wav")

file_list <- setdiff(file_list, file_list_trunc)

# load duration data
load("Recording list with durations.Rdata")

file_df <- 
  data.frame(paths = file_list) %>% 
  separate(paths, into = c("folder", "Songmeter", "Recording"), sep = "/", remove = FALSE) %>% 
  left_join(df_durations) %>% 
  arrange(paths)

# check how many recording minutes
file_df %>% 
  group_by(Duration_min) %>% 
  summarise(n = n(),
            mins = sum(Duration_min)) %>% 
  mutate(proc_hrs_est = (mins/5.5)/60/60)

# truncate 5 hr recordings to 60 mins
wav_trunc <- function(WAV) {
   wav <-  tuneR::readWave(WAV, from = 0, to = 3600, units = "seconds")
   wav %>% tuneR::writeWave(str_replace(WAV, ".wav", "_trunc60.wav"))
   rm(wav)
}

wav_long <- 
  file_df %>% 
  filter(Duration_min > 60) %>% 
  pull(paths)

wav_done <- 
  list.files(path = "D:", recursive = TRUE, full.names = TRUE, pattern = "_trunc60.wav") %>% 
  str_replace(., "_trunc60.wav", ".wav")

wav_todo <- base::setdiff(wav_long, wav_done)

tic()
pb <- progress_bar$new(
  format = "Crunching... [:bar] :percent in :elapsed, ETA: :eta",
  total = length(wav_todo), clear = FALSE, width = 80)

map(wav_todo, ~{pb$tick(); wav_trunc(.)})
beepr::beep(2)
toc()


# make a "to-do" list subset for hardRain
todo <- 
  file_df %>% 
  filter(Duration_min == 31) %>% 
  #slice(1:13890) %>% 
  pull(paths)

todo_trunc <- file_list_trunc

todo_misc  <- 
  file_df %>% 
  filter(Recording %in% c("SM7687_20170724_051800.wav", "SM13520_20190515_163000.wav")) %>% 
  pull(paths)
  

# split to-do list into chunks of 6 recordings 
# (number of processing cores available)
create_chunks <- function(x, n) {
  split(x, cut(seq_along(x), length(x)/n, labels = FALSE))
}

chunk_list <- create_chunks(todo_trunc, 6)

# create custom hardRain-running function
do_hardRain <- function(TODO) {
  pb$tick()
  hardRain::getMetrics(wav = TODO, 
                       t.step = 60, 
                       parallel = TRUE) %>% 
  as_tibble(rownames = "Recording") %>%
  group_by(Recording) %>% 
  mutate(Rec_minute = seq(1, length(Recording)), .after = Recording) %>% 
  write_excel_csv("hardRain results.csv", append = TRUE)
}  

# run function on list of chunks
detach("package:tidylog", unload = TRUE)

tic() 

pb <- progress_bar$new(
  format = "Crunching... [:bar] :percent in :elapsed, ETA: :eta",
  total = length(chunk_list), clear = FALSE, width = 80)

hr <- map_dfr(todo_misc, do_hardRain)

beepr::beep(2)

toc()

# read in csv created by function
library(tidylog)

csv <- read_csv("hardRain results.csv") %>% 
  distinct()

write_excel_csv(csv, path = "hardRain results.csv")

sort(unique(csv$Recording)) %>% length()

df_hardRain <- csv %>% 
  mutate(Recording = str_replace(Recording, "_trunc60.wav", ".wav"))

save(df_hardRain, file = "hardRain results.Rdata")



#8. Combine everything with list of 1 minute recording intervals 
#   (validation, lunar/solar covs, weather covs, hard rain)
#   (and then calculate time-since and time-until vars)

df_analysis_occupancy <-
  # combine
  df_1min_occ_abd %>% 
  left_join(df_covar_prep, by = c("Recording", "Rec_minute")) %>% 
  left_join(df_sunmoon) %>% 
  left_join(df_weather_hourly_interp) %>% 
  left_join(wx_daily, by = c("Date_EWPW" = "date", "station_name", "station_id")) %>% 
  left_join(df_hardRain) %>% 
  # tidy
  distinct() %>% 
  select(-time) %>% 
  # calculate derivative covariates
  mutate(sunsetTimeSince  = round(difftime(Time_rec_1min, sunset, units = "hours"), 2),
         sunsetTimeSince  = if_else(Hour_rec_1min < 12, sunsetTimeSince + 24, sunsetTimeSince),
         sunriseTimeSince = round(difftime(Time_rec_1min, sunrise, units = "hours"), 2),
         sunriseTimeSince = if_else(Hour_rec_1min > 12, sunriseTimeSince - 24, sunriseTimeSince),
         .after = goldenHour)

# save final data and full workspace
save(df_analysis_occupancy, file = "EWPW occupancy analysis data.Rdata")

save.image("EWPW wrangling v1.RData")

####
# shortcut to the end:
# load("EWPW wrangling v1.RData")
####