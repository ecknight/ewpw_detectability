#title: "Refining surveys for nocturnal birds using acoustic data: Eastern Whip-poor-will"
#author: "Elly C. Knight"
#date: "December 17, 2020"

#Preamble#####

options(scipen = 9999, max.print=2000)

#Load required R packages
library(tidyverse)
library(lubridate)
library(usdm)
library(GGally)
library(unmarked)
library(MuMIn)
library(AICcmodavg)
library(data.table)
library(readxl)
library(gridExtra)
library(formatR)
library(sf)
library(ggmap)
library(inauguration)
library(cowplot)
library(patchwork)
library(lme4)
library(merTools)
library(mgcv)

my.theme <- theme_classic() +
  theme(text=element_text(size=12, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        plot.title=element_text(size=12, hjust = 0.5))

#Load data
load("Analysisv7_workspace.Rdata")

#SECTION 1. Recognizer evaluation####

##1a. Prepare data####

#Read in and summarize benchmark data (human listener processed)
path <- "/Users/ellyknight/Documents/UoA/Projects/Projects/EWPW/"
bench <- read_excel(paste0(path, "Data/TestDatasetCallCount.xlsx"),
                    col_types = c("text", "text", "numeric", "numeric", "text", "numeric", "numeric", "numeric", "text", "text")) %>% 
  mutate(File=paste0(File, ".wav")) %>% 
  rename(file=File) %>% 
  mutate(LevelNo = case_when(Level=="barely perceptible" ~ 1,
                             Level %in% c("barely perceptible to medium", "barely perceptible to quiet", "quiet to barely perceptible") ~ 2,
                             Level=="quiet" ~ 3,
                             Level=="quiet to medium" ~ 4,
                             Level=="medium" ~ 5,
                             Level=="medium to loud" ~ 6,
                             Level== "loud" ~7,
                             Level=="very loud" ~ 8)) %>% 
  mutate(starttime = 60*StartMin + StartSec,
         endtime = 60*EndMin + EndSec,
         totalmins = endtime - starttime + 1) %>% 
  dplyr::filter(StartMin < 5) %>% 
  as.data.frame() %>% 
  mutate(id = row_number())

bench.sum <- bench %>% 
  group_by(LevelNo) %>% 
  summarize(bouts = n(),
            duration = sum(totalmins)/60,
            calls = sum(Count)) %>% 
  ungroup()

write.csv(bench.sum, "LoudnessRankingSummary.csv", row.names=FALSE)

bench.rec <- bench %>% 
  dplyr::select(file) %>% 
  unique()

call <- sum(bench$Count)

bench.level <- bench %>% 
  group_by(LevelNo) %>% 
  summarize(calls = sum(Count))

bench.slice <- as.data.frame(lapply(bench, rep, bench$totalmins + 2)) %>% 
  group_by(id) %>% 
  mutate(row = row_number(),
         time = starttime + row - 2) %>% 
  ungroup() %>% 
  dplyr::select(id, file, Count, LevelNo, time)

#Read in raw validated data
dat.test <- read.table("/Users/ellyknight/Documents/UoA/Projects/Projects/EWPW/Processing/Test_EWPWV1_0_20_results_validated.txt",
                       sep="\t",
                       header=FALSE,
                       col.names = c("filepath", "start", "duration", "rsl", "quality", "score", "recognizer", "validation")) %>% 
  separate(filepath, into=c("f1", "f2", "f3", "f4", "presence", "file"), sep="/", remove=FALSE) %>% 
  separate(file, into=c("aru", "datetext", "timetext", "filetype"), remove=FALSE) %>% 
  mutate(ewpw=ifelse(validation=="y", 1, 0)) %>% 
  separate(start, into=c("starthour", "startminute", "startsecond"), sep=":", remove=FALSE) %>% 
  mutate(time = round(60*as.numeric(startminute) + as.numeric(startsecond))) %>% 
  dplyr::filter(as.numeric(startminute) < 5) %>% 
  left_join(bench.slice)


##1b. Evaluate for detection####

#score <- seq(round(min(dat.test$score)), round(max(dat.test$score)))
score <- seq(20, 80)
results.call <- data.frame()
for(i in 1:length(score)){
  
  score.i <- as.numeric(score[i])
  
  dat.i <- dat.test %>% 
    dplyr::filter(score >= score.i)
  
  level.i <- dat.i %>% 
    group_by(ewpw) %>% 
    summarize(rslsd = sd(rsl),
              rsl = mean(rsl),
              level = mean(LevelNo, na.rm=TRUE),
              levelsd = sd(LevelNo, na.rm=TRUE)) %>% 
    dplyr::filter(ewpw==1)
  
  results.i <- dat.i %>% 
    dplyr::summarize(det=sum(ewpw), hit=n()) %>% 
    mutate(r=det/call, p=det/hit, fp=(hit-det)/hit, tp=det/hit, fn=(call-det)/call, acc=(det-(hit-det))/hit, leveltp = level.i$level, rsltp = level.i$rsl, levelmintp = level.i$levelmin, levelsdtp = level.i$levelsd, rslsdtp = level.i$rslsd) %>% 
    mutate(f=(2*p*r)/(p+r),
           score=score[i]) %>% 
    ungroup()

  
  results.call <- rbind(results.call, results.i)
  
}

##1c. Evaluate for presence/absence per recording####

results.rec <- data.frame()
for(i in 1:length(score)){
  score.i <- as.numeric(score[i])
  
  dat.i <- dat.test %>% 
    dplyr::filter(score >= score.i) %>% 
    group_by(file, presence) %>% 
    summarize(det=sum(ewpw), hit=n()) %>% 
    ungroup() %>% 
    full_join(bench.rec) %>% 
    mutate(det = ifelse(is.na(det), 0, det)) %>% 
    mutate(hit = ifelse(is.na(hit), 0, hit)) %>% 
    mutate(pres.rec = ifelse(det > 0, 1, 0),
           score=score.i)
  
  results.rec <- rbind(results.rec, dat.i)
  
}

r.rec <- results.rec %>% 
  group_by(score) %>% 
  summarize(rpres=sum(pres.rec)/nrow(bench.rec)) %>% 
  ungroup()

##1d. Put together----
r.eval <- left_join(results.call, r.rec)

write.csv(r.eval, "RecognizerEvaluationResults.csv", row.names = FALSE)

#SECTION 2. Summary of automated processing####

#Read in clean dataset summarized by 1 minute of recording
load("/Users/ellyknight/Dropbox/SS/EWPW/EWPW occupancy analysis data.Rdata") 

#ID sites with detections
sites <- df_analysis_occupancy %>% 
  dplyr::filter(as.numeric(sunsetTimeSince) > 0,
                as.numeric(sunriseTimeSince) < 0, 
                Rec_minute <= 5) %>% 
  group_by(SiteName, Year) %>% 
  summarize(detections=sum(Occupied)) %>% 
  dplyr::filter(detections > 0) %>% 
  ungroup() %>% 
  left_join(df_analysis_occupancy) %>% 
  dplyr::select(SiteName, Year, SongMeter, Latitude, Longitude, detections) %>% 
  unique()

#filter to just the data used in analysis
df_analysis <- df_analysis_occupancy %>% 
  mutate(ID=paste0(Year, "-", SiteName)) %>% 
  dplyr::filter(as.numeric(sunsetTimeSince) > -1,
                as.numeric(sunriseTimeSince) < 1, 
                Rec_minute <= 5) %>% 
  group_by(Recording) %>% 
  mutate(Rec_length = max(Rec_minute)) %>% 
  ungroup() %>% 
  dplyr::filter(Rec_length==5) %>% 
  dplyr::filter(SiteName %in% sites$SiteName)

## 2a. Precision of recognizer####

#Wrangle raw validated text files from SongScope
files <- data.frame(file=list.files("/Users/ellyknight/Dropbox/SS/EWPW/Processing_finished", pattern="*validated.txt", full.names=TRUE))
val.list <- list()
for(i in 1:nrow(files)){
  val.list[[i]] <- read.table(as.character(files$file[i]),
                              sep="\t",
                              col.names = c("filepath", "start", "duration", "rsl", "quality", "score", "recognizer", "validation"))
}
val.raw <- rbindlist(val.list) %>% 
  mutate(validation=case_when(is.na(validation) ~ 0,
                              validation %in% c("1", "11", "2", "2 here", "3?", "Minimum 2", "Minimum 3", "3", "4") ~ 1,
                              validation %in% c("0", "na", "Na", "?", "") ~ 0)) %>% 
  mutate(filepath = str_replace(filepath, "D:\\\\", "")) %>% 
  mutate(filepath = str_replace(filepath, "E:\\\\", "")) %>% 
  mutate(filepath = str_replace(filepath, "/Volumes/Backup Plus/", "")) %>% 
  mutate(filepath = str_replace_all(filepath, "\\\\", "/")) %>% 
  separate(filepath, into = c("Project", "SongMeter", "Recording"), sep = "/") %>% 
  dplyr::filter(Recording %in% df_analysis$Recording)

#Summarize raw validated data
#True positives
tp <- nrow(val.raw %>% 
             dplyr::filter(validation==1))
tp

# False positives
fp <- nrow(val.raw %>% 
             dplyr::filter(validation==0))
fp

# Precision
tp/(fp+tp)

#Interrogate false positives
val.time <- val.raw %>% 
  separate(Recording, into=c("ID", "datename", "timename"), sep="_", remove=FALSE) %>% 
  mutate(DateTime = ymd_hms(paste0(datename, str_sub(timename, 1, 6))),
         doy = yday(DateTime),
         hour = hour(DateTime))

ggplot(val.time) +
  geom_histogram(aes(x=doy, colour=factor(validation)))

ggplot(val.time) +
  geom_histogram(aes(x=hour, colour=factor(validation)))

#Number of false positives at 4 & 5 am (i.e., dawn chorus)
fp.dawn <- val.time %>% 
  dplyr::filter(validation==0, hour %in% c(4,5)) %>% 
  nrow()

fp.dawn/fp

#Precision without dawn chorus
tp.dawn <- val.time %>% 
  dplyr::filter(validation==1, hour %in% c(4,5)) %>% 
  nrow()

(tp - tp.dawn)/((fp-fp.dawn)+(tp-tp.dawn))

##2b. Recording availability####

#Length and number of recordings per location
recordings <- df_analysis %>% 
  mutate(ID=paste0(Year, "-", SiteName)) %>% 
  group_by(Year, ID, SiteName, Recording) %>% 
  summarize(Duration = max(Rec_minute)) %>% 
  ungroup()
table(recordings$ID, recordings$Duration)

#Time and date of recordings per location
recordings.time <- df_analysis %>% 
  mutate(ID=paste0(Year, "-", SiteName)) %>% 
  #  dplyr::filter(Rec_minute <= 5) %>% 
  group_by(Year, ID, SiteName, Recording, Hour_rec_1min, Date) %>% 
  summarize(Duration = max(Rec_minute),
            Occupied = ifelse(sum(Occupied) > 0, 1, 0)) %>% 
  ungroup()
table(recordings.time$ID, recordings.time$Hour_rec_1min)

#Total number of recording minutes per location
recordings.sum <- recordings %>% 
  group_by(Year, ID, SiteName) %>% 
  summarize(Duration.hr=sum(Duration)/60) %>% 
  ungroup() %>% 
  arrange(Duration.hr)
recordings.sum

#Make appendix 1
appendix1 <- table(recordings.time$SiteName, recordings.time$Hour_rec_1min) %>% 
  data.frame() %>% 
  rename(SiteName=Var1, Hour=Var2, Recordings=Freq) %>% 
  mutate(Year = ifelse(str_sub(SiteName, 1, 3)=="PEP", 2019, 2017)) %>% 
  pivot_wider(id_cols=c(SiteName, Year), names_from=Hour, names_prefix="Hr", values_from=Recordings) %>% 
  mutate(TotalRecordings=Hr0 + Hr1 + Hr2 + Hr3 + Hr4 + Hr5 + Hr6 + Hr20 + Hr21 + Hr22 + Hr23,
         TotalMinutes=TotalRecordings*5) %>% 
  arrange(Year, TotalMinutes) %>% 
  dplyr::select(-TotalRecordings)

write.csv(appendix1, "Appendix1.csv", row.names = FALSE)

#Total # of recording minutes
sum(appendix1$TotalMinutes)
sum(appendix1$TotalMinutes)/60

##2c. Summary of detections####

summary(df_analysis)

#Detections per recording
detections.rec <- df_analysis %>% 
  group_by(Year, SiteName, Recording) %>% 
  summarize(DetectionMinutes = sum(Occupied),
            Num_detects = sum(Num_detects)) %>% 
  ungroup()
summary(detections.rec)
sd(detections.rec$Num_detects)

#Summarize by site
detections.site <- df_analysis %>% 
  group_by(Year, SiteName, Recording) %>% 
  summarize(DetectionMinutes = sum(Occupied),
            Num_detects = sum(Num_detects),
            DetectionRecordings = ifelse(DetectionMinutes > 0, 1, 0)) %>% 
  group_by(Year, SiteName) %>% 
  summarize(Num_detects = sum(Num_detects),
            DetectionMinutes = sum(DetectionMinutes),
            DetectionRecordings = sum(DetectionRecordings),
            Recordings=n()) %>% 
  ungroup() %>% 
  left_join(recordings.sum) %>% 
  mutate(Occupied = ifelse(DetectionMinutes > 0, 1, 0),
         DetectsPerOccuMinute = Num_detects/Duration.hr,
         DetectsPerRecMinute = Num_detects/Recordings*5)

hist(detections.site$DetectionRecordings)

#Detections per minute of occupied recording at occupied sites
detections.site %>% 
  dplyr::filter(Occupied==1) %>% 
  dplyr::select(DetectsPerOccuMinute) %>% 
  summary()

#Detections per minute of recording at occupied sites
detections.site %>% 
  dplyr::filter(Occupied==1) %>% 
  dplyr::select(DetectsPerRecMinute) %>% 
  summary()
sd(detections.site$DetectsPerRecMinute)

#Specific detections for sites with few detections
detections.summary <- detections.rec %>% 
  dplyr::filter(Num_detects > 0) %>% 
  dplyr::select(Year, SiteName, Recording) %>% 
  left_join(detections.site)

#SECTION 3: Preliminary detectability analysis####

##3a. Wrangle data####

#Summarize by recording
dat <- df_analysis %>% 
  inner_join(duration) %>% 
  dplyr::filter(Rec_minute <= 5) %>% 
  mutate(moonUp = ifelse(moonUp=="Y", 1, 0)) %>% 
  group_by(SiteName, Latitude, Longitude, Year, Yday, Date, Recording, Duration, Time_rec) %>% 
  summarize(detectionsRecording = sum(Occupied),
            sunsetTimeSince = mean(as.numeric(sunsetTimeSince)),
            moon_fraction = mean(moon_fraction),
            moon_altitude = mean(moon_altitude),
            temp = mean(temp),
            wind_spd = mean(wind_spd),
            total_rain = mean(total_rain),
            moonUp = round(mean(moonUp)),
            band1_psd = mean(band1_psd),
            band1_s2n = mean(band1_s2n),
            band2_psd = mean(band2_psd),
            band2_s2n = mean(band2_s2n),
            Abundance_recording=mean(Abundance_recording),
            Abundance_site=mean(Abundance_site)) %>% 
  mutate(Occupied = ifelse(detectionsRecording > 0, 1, 0),
         moonUp_altitude = ifelse(moonUp==1, moon_altitude, 0)) %>% 
  dplyr::filter(!is.na(total_rain),
                !is.na(band1_psd)) %>% 
  group_by(SiteName, Year) %>% 
  arrange(Recording) %>% 
  mutate(n=row_number()) %>% 
  ungroup()

#Visualize
ggplot(dat) +
  geom_jitter(aes(x=sunsetTimeSince, y=Occupied)) + 
  geom_smooth(aes(x=sunsetTimeSince, y=Occupied))

ggplot(dat) +
  geom_jitter(aes(x=Yday, y=Occupied)) + 
  geom_smooth(aes(x=Yday, y=Occupied))

ggplot(dat, aes(x=moon_fraction, y=Occupied)) +
  geom_jitter() + 
  geom_smooth()

ggplot(dat, aes(x=moon_altitude, y=Occupied)) +
  geom_jitter() + 
  geom_smooth()

ggplot(dat, aes(x=temp, y=Occupied)) +
  geom_jitter() + 
  geom_smooth()

ggplot(dat, aes(x=total_rain, y=Occupied)) +
  geom_jitter() + 
  geom_smooth()

ggplot(dat, aes(x=wind_spd, y=Occupied)) +
  geom_jitter() + 
  geom_smooth()

ggplot(dat, aes(x=band2_psd, y=Occupied)) +
  geom_jitter() + 
  geom_smooth()

ggplot(dat, aes(x=band2_s2n, y=Occupied)) +
  geom_hex() + 
  geom_smooth()

##3b. Check for VIF & correlation####

#Check VIF & correlation
covs <- dat %>% 
  dplyr::select(Yday, sunsetTimeSince, moon_fraction, moon_altitude, temp, total_rain, wind_spd, band1_psd, band2_psd, band1_s2n, band2_s2n) %>% 
  data.frame()

usdm::vif(covs)
cor(covs)
#Remove band1_psd

covs <- dat %>% 
  dplyr::select(Yday, sunsetTimeSince, moon_fraction, moon_altitude, temp, total_rain, wind_spd, band2_psd, band1_s2n, band2_s2n) %>% 
  data.frame()

usdm::vif(covs)
cor(covs)
#all good

##3c. Covariate selection####

#Create lists for saving results
temp.dredge.list <- list()
weather.dredge.list <- list()

#Set number of bootstraps
boot <- 100

#Start for loop
for(i in 1:boot){
  
  #Subsample data
  dat.sub <- dat %>% 
    group_by(SiteName) %>% 
    sample_n(50, replace=TRUE) %>%
    ungroup() %>% 
    unique()
  
  #Format for occupancy
  detections <- dat.sub %>% 
    dplyr::select(SiteName, n, Occupied) %>% 
    spread(key = n, value = Occupied) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  doy <- dat.sub %>% 
    dplyr::select(SiteName, n, Yday) %>% 
    spread(key = n, value = Yday) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  suntime <- dat.sub %>% 
    dplyr::select(SiteName, n, sunsetTimeSince) %>% 
    spread(key = n, value = sunsetTimeSince) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  moonfrac <- dat.sub %>% 
    dplyr::select(SiteName, n, moon_fraction) %>% 
    spread(key = n, value = moon_fraction) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  moonalt <- dat.sub %>% 
    dplyr::select(SiteName, n, moon_altitude) %>% 
    spread(key = n, value = moon_altitude) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  moonupalt <- dat.sub %>% 
    dplyr::select(SiteName, n, moonUp_altitude) %>% 
    spread(key = n, value = moonUp_altitude) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  temp <- dat.sub %>% 
    dplyr::select(SiteName, n, temp) %>% 
    spread(key = n, value = temp) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  rain <- dat.sub %>% 
    dplyr::select(SiteName, n, total_rain) %>% 
    spread(key = n, value = total_rain) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  wind <- dat.sub %>% 
    dplyr::select(SiteName, n, wind_spd) %>% 
    spread(key = n, value = wind_spd) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  psd2 <- dat.sub %>% 
    dplyr::select(SiteName, n, band2_psd) %>% 
    spread(key = n, value = band2_psd) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  s2n1 <- dat.sub %>% 
    dplyr::select(SiteName, n, band1_s2n) %>% 
    spread(key = n, value = band1_s2n) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  s2n2 <- dat.sub %>% 
    dplyr::select(SiteName, n, band2_s2n) %>% 
    spread(key = n, value = band2_s2n) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  obs.cov <- list(doy=doy, suntime=suntime, moonupalt=moonupalt, moonfrac=moonfrac, moonalt=moonalt, temp=temp, rain=rain, wind=wind, psd2=psd2, s2n1=s2n1, s2n2=s2n2)
  
  site.cov <- dat.sub %>% 
    dplyr::select(SiteName, Latitude, Longitude) %>% 
    unique() %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  dat.occ <- unmarkedFrameOccu(detections, siteCovs=site.cov, obsCovs=obs.cov)
  
  #Fit models
  
  #Temporal first
  mod.temp <- occu(~ moonfrac*moonalt + doy + suntime + I(suntime^2)
                   ~1,
                   data=dat.occ)
  temp.dredge <-  try(dredge(mod.temp))
  if(class(temp.dredge)=="model.selection"){
    temp.dredge.list[[i]] <- temp.dredge %>% 
      data.frame() %>% 
      mutate(boot=i,
             mod=row.names(temp.dredge))
  }
  
  #Weather next
  mod.weather <- occu(~ temp + I(temp^2) + rain + wind + psd2 + s2n2
                      ~1,
                      data=dat.occ)
  weather.dredge <-  try(dredge(mod.weather))
  if(class(weather.dredge)=="model.selection"){
    weather.dredge.list[[i]] <- weather.dredge %>% 
      data.frame() %>% 
      mutate(boot=i,
             mod=row.names(weather.dredge))
  }
  
  
  print(paste0("Finished bootstrap ", i, " of ", boot))
  
}

#Summarize bootstrapped models
temp.dredge.raw <- rbindlist(temp.dredge.list)

write.csv(temp.dredge.raw, "TemporalCovariatesModelSelection_raw.csv", row.names=FALSE)

temp.dredge.mn <- rbindlist(temp.dredge.list) %>% 
  group_by(mod) %>% 
  summarize_all(~mean(.x)) %>% 
  ungroup() %>% 
  arrange(-weight) %>% 
  mutate(metric="mean")
temp.dredge.mn

temp.dredge.sd <- rbindlist(temp.dredge.list) %>% 
  group_by(mod) %>% 
  summarize_all(~sd(.x)) %>% 
  ungroup() %>% 
  arrange(-weight) %>% 
  mutate(metric="sd")
temp.dredge.sd

temp.dredge.all <- rbind(temp.dredge.mn, temp.dredge.sd)

write.csv(temp.dredge.all, "TemporalCovariatesModelSelection.csv", row.names=FALSE)

temp.dredge.top <- temp.dredge.all %>% 
  head(1) %>% 
  pivot_longer(p.doy.:p.moonalt.moonfrac., names_to="variable", values_to="coefficient") %>% 
  dplyr::select(variable, coefficient, delta, weight) %>% 
  dplyr::filter(!is.na(coefficient))
temp.dredge.top
#Select Doy, suntime, suntime^2, moonalt (everything but moon frac)

weather.dredge.raw <- rbindlist(weather.dredge.list)

write.csv(weather.dredge.raw, "WeatherCovariatesModelSelection_raw.csv", row.names=FALSE)

weather.dredge.mn <- rbindlist(weather.dredge.list) %>% 
  group_by(mod) %>% 
  summarize_all(~mean(.x)) %>% 
  ungroup() %>% 
  arrange(-weight) %>% 
  mutate(metric="mean")
weather.dredge.mn

weather.dredge.sd <- rbindlist(weather.dredge.list) %>% 
  group_by(mod) %>% 
  summarize_all(~sd(.x)) %>% 
  ungroup() %>% 
  arrange(-weight) %>% 
  mutate(metric="sd")
weather.dredge.sd

weather.dredge.all <- rbind(weather.dredge.mn, weather.dredge.sd)

write.csv(weather.dredge.all, "WeatherCovariatesModelSelection.csv", row.names=FALSE)

weather.dredge.top <- weather.dredge.all %>% 
  head(1) %>% 
  pivot_longer(p.psd2.:p.wind., names_to="variable", values_to="coefficient") %>% 
  dplyr::select(variable, coefficient, delta, weight) %>% 
  dplyr::filter(!is.na(coefficient))
weather.dredge.top
#Select psd2, ps2n2, temp, temp^2, wind (everything but rain)

##3d. Fit final model####

#Create new data for model prediction in loop
moondat <- expand.grid(moonalt=seq(round(min(dat$moon_altitude), -1), max(dat$moon_altitude), 0.1),
                       suntime=mean(dat$sunsetTimeSince),
                       doy=min(dat$Yday),
                       wind=min(dat$wind_spd),
                       rain=min(dat$total_rain),
                       temp=mean(dat$temp),
                       psd2=min(dat$band2_psd),
                       s2n2=min(dat$band2_s2n))

suntimedat <- expand.grid(moonalt=max(dat$moon_altitude),
                          suntime=seq(round(min(dat$sunsetTimeSince), -1), max(dat$sunsetTimeSince), 1),
                          doy=min(dat$Yday),
                          wind=min(dat$wind_spd),
                          rain=min(dat$total_rain),
                          temp=mean(dat$temp),
                          psd2=min(dat$band2_psd),
                          s2n2=min(dat$band2_s2n))

doydat <- expand.grid(moonalt=max(dat$moon_altitude),
                      suntime=mean(dat$sunsetTimeSince),
                      doy=seq(min(dat$Yday), max(dat$Yday), 1),
                      wind=min(dat$wind_spd),
                      rain=min(dat$total_rain),
                      temp=mean(dat$temp),
                      psd2=min(dat$band2_psd),
                      s2n2=min(dat$band2_s2n))

winddat <- expand.grid(moonalt=max(dat$moon_altitude),
                       suntime=mean(dat$sunsetTimeSince),
                       doy=min(dat$Yday),
                       wind=seq(round(min(dat$wind_spd)), max(dat$wind_spd), 1),
                       rain=min(dat$total_rain),
                       temp=mean(dat$temp),
                       psd2=min(dat$band2_psd),
                       s2n2=min(dat$band2_s2n))

tempdat <- expand.grid(moonalt=max(dat$moon_altitude),
                       suntime=mean(dat$sunsetTimeSince),
                       doy=min(dat$Yday),
                       wind=min(dat$wind_spd),
                       rain=min(dat$total_rain),
                       temp=seq(min(round(dat$temp)), max(dat$temp), 1),
                       psd2=min(dat$band2_psd),
                       s2n2=min(dat$band2_s2n))


psddat <- expand.grid(moonalt=max(dat$moon_altitude),
                      suntime=mean(dat$sunsetTimeSince),
                      doy=min(dat$Yday),
                      wind=min(dat$wind_spd),
                      rain=min(dat$total_rain),
                      temp=mean(dat$temp),
                      psd2=seq(0, round(max(dat$band2_psd), 3), 0.001),
                      s2n2=min(dat$band2_s2n))

s2ndat <- expand.grid(moonalt=max(dat$moon_altitude),
                      suntime=mean(dat$sunsetTimeSince),
                      doy=min(dat$Yday),
                      wind=min(dat$wind_spd),
                      rain=min(dat$total_rain),
                      temp=mean(dat$temp),
                      psd2=min(dat$band2_psd),
                      s2n2=seq(round(min(dat$band2_s2n), 2), round(max(dat$band2_s2n), 2), 0.01))

#Create lists to save out results
mod.final.coeff.list <- list()
mod.final.pred.list <- list()

#Set number of bootstraps  
boot <- 100

#Start for loop
for(i in 1:boot){
  
  #Sample data
  dat.sub <- dat %>% 
    group_by(SiteName) %>% 
    sample_n(50, replace=TRUE) %>%
    ungroup() %>% 
    unique()
  
  #Format for occupancy  
  detections <- dat.sub %>% 
    dplyr::select(SiteName, n, Occupied) %>% 
    spread(key = n, value = Occupied) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  abundance.site <- dat.sub %>% 
    dplyr::select(SiteName, n, Abundance_site) %>% 
    spread(key = n, value = Abundance_site) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  abundance.rec <- dat.sub %>% 
    dplyr::select(SiteName, n, Abundance_recording) %>% 
    spread(key = n, value = Abundance_recording) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  year <- dat.sub %>% 
    dplyr::select(SiteName, n, Year) %>% 
    spread(key = n, value = Year) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  doy <- dat.sub %>% 
    dplyr::select(SiteName, n, Yday) %>% 
    spread(key = n, value = Yday) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  suntime <- dat.sub %>% 
    dplyr::select(SiteName, n, sunsetTimeSince) %>% 
    spread(key = n, value = sunsetTimeSince) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  moonalt <- dat.sub %>% 
    dplyr::select(SiteName, n, moon_altitude) %>% 
    spread(key = n, value = moon_altitude) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  temp <- dat.sub %>% 
    dplyr::select(SiteName, n, temp) %>% 
    spread(key = n, value = temp) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  rain <- dat.sub %>% 
    dplyr::select(SiteName, n, total_rain) %>% 
    spread(key = n, value = total_rain) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  wind <- dat.sub %>% 
    dplyr::select(SiteName, n, wind_spd) %>% 
    spread(key = n, value = wind_spd) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  psd2 <- dat.sub %>% 
    dplyr::select(SiteName, n, band2_psd) %>% 
    spread(key = n, value = band2_psd) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  s2n2 <- dat.sub %>% 
    dplyr::select(SiteName, n, band2_s2n) %>% 
    spread(key = n, value = band2_s2n) %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  obs.cov <- list(doy=doy, suntime=suntime, moonalt=moonalt, temp=temp, rain=rain, wind=wind, psd2=psd2, s2n2=s2n2, abundance.site=abundance.site, abundance.rec=abundance.rec, year=year)
  
  site.cov <- dat.sub %>% 
    dplyr::select(SiteName, Latitude, Longitude) %>% 
    unique() %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  dat.occ <- unmarkedFrameOccu(detections, siteCovs=site.cov, obsCovs=obs.cov)
  
  #Fit model 
  mod.final <- occu(~ moonalt + suntime + I(suntime^2) + doy + temp + I(temp^2) + wind + psd2 + s2n2
                    ~1,
                    data=dat.occ)
  
  #Save coefficient estimates
  invisible({capture.output({
    mod.final.coeff.list[[i]] <- summary(mod.final)[['det']] %>% 
      mutate(boot=i,
             var=row.names(data.frame(mod.final@estimates@estimates[["det"]]@estimates)))
  })})
  
  #Predict on new data
  mod.final.pred.list[[i]] <- moondat %>% 
    cbind(predict(mod.final, type="det", newdata=moondat)) %>% 
    dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
    mutate(plot="Moon altitude") %>% 
    rbind(suntimedat %>% 
            cbind(predict(mod.final, type="det", newdata=suntimedat)) %>% 
            dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
            mutate(plot="Time since sunset (hrs)")) %>% 
    rbind(winddat %>% 
            cbind(predict(mod.final, type="det", newdata=winddat)) %>% 
            dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
            mutate(plot="Wind speed (km/h)")) %>% 
    rbind(tempdat %>% 
            cbind(predict(mod.final, type="det", newdata=tempdat)) %>% 
            dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
            mutate(plot="Temperature (C)")) %>% 
    rbind(doydat %>% 
            cbind(predict(mod.final, type="det", newdata=doydat)) %>% 
            dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
            mutate(plot="Day of year")) %>% 
    rbind(psddat %>% 
            cbind(predict(mod.final, type="det", newdata=psddat)) %>% 
            dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
            mutate(plot="Power spectrum density")) %>% 
    rbind(s2ndat %>% 
            cbind(predict(mod.final, type="det", newdata=s2ndat)) %>% 
            dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
            mutate(plot="Signal to noise ratio")) %>% 
    mutate(boot=i)
  
  print(paste0("Finished bootstrap ", i, " of ", boot))
  
}

#Summarize bootstrapped coefficients
mod.final.coeff <- rbindlist(mod.final.coeff.list)
mod.final.coeff.sum <- mod.final.coeff %>% 
  group_by(var) %>% 
  summarize(est.mn=mean(Estimate),
            est.sd=sd(Estimate),
            se.mn=mean(SE),
            se.sd=sd(SE)) %>% 
  ungroup()
mod.final.coeff.sum

write.csv(mod.final.coeff, "FinalModelCoefficients.csv", row.names = FALSE)


#Summarize model predictions
mod.final.pred <- rbindlist(mod.final.pred.list)
mod.final.pred.sum <- mod.final.pred %>% 
  group_by(plot, moonalt, suntime, wind, temp) %>% 
  summarize(det.mn = mean(Det),
            se.mn = mean(DetSE),
            up.mn = mean(DetUpper), 
            lw.mn = mean(DetLower)) %>% 
  ungroup()

write.csv(mod.final.pred, "FinalModelPredictions.csv", row.names = FALSE)

##3e. Plot final model####
plot.moon <- ggplot(subset(mod.final.pred, plot=="Moon altitude")) +
  geom_line(aes(x=moonalt, y=det.mn), colour="blue") +
  geom_ribbon(aes(x=moonalt, ymin=lw.mn, ymax=up.mn), alpha=0.5) +
  ylim(c(0,1)) +
  xlab("Moon altitude") +
  ylab("EWPW detectability")

plot.suntime <- ggplot(subset(mod.final.pred, plot=="Time since sunset (hrs)")) +
  geom_line(aes(x=suntime, y=det.mn), colour="blue") +
  geom_ribbon(aes(x=suntime, ymin=lw.mn, ymax=up.mn), alpha=0.5) +
  ylim(c(0,1)) +
  xlab("Time since sunset (hrs)") +
  ylab("EWPW detectability")

plot.wind <- ggplot(subset(mod.final.pred, plot=="Wind speed (km/h)")) +
  geom_line(aes(x=wind, y=det.mn), colour="blue") +
  geom_ribbon(aes(x=wind, ymin=lw.mn, ymax=up.mn), alpha=0.5) +
  ylim(c(0,1)) +
  xlab("Wind speed (km/h)") +
  ylab("EWPW detectability")

plot.temp <- ggplot(subset(mod.final.pred, plot=="Temperature (C)")) +
  geom_line(aes(x=temp, y=det.mn), colour="blue") +
  geom_ribbon(aes(x=temp, ymin=lw.mn, ymax=up.mn), alpha=0.5) +
  ylim(c(0,1)) +
  xlab("Temperature (C)") +
  ylab("EWPW detectability")

#Put all four plots together
plot.pred <- gridExtra::grid.arrange(plot.moon, plot.suntime, plot.wind, plot.temp,
                                     ncol=2, nrow=2)

#SECTION 4: Survey effort analysis####

#4a. No constraints on covariates----

#Set number of bootstraps  
boot <- 100

#Set protocol options
length <- c(1:5)
#length <- 5
visit <- c(1:30, 40, 50, 60, 70, 80, 90, 100)
#visit <- c(1:30)
#visit <- 20

protocol <- expand.grid(length=length, visit=visit)

#Create lists to save out results
mod.any.null.pred.list1 <- list()
mod.any.null.pred.list2 <- list()
mod.any.cov.pred.list1 <- list()
mod.any.cov.pred.list2 <- list()
summary.any <- data.frame()

for(h in 1:nrow(protocol)){
  
  #Set protocol for this set
  length.h <- protocol$length[h]
  visit.h <- protocol$visit[h]
  
  #Wrangle data
  dat.h <- df %>% 
    inner_join(duration) %>% 
    dplyr::filter(Rec_minute <= length.h) %>% 
    mutate(moonUp = ifelse(moonUp=="Y", 1, 0)) %>% 
    group_by(SiteName, Latitude, Longitude, Year, Yday, Date, Recording, Duration, Time_rec) %>% 
    summarize(detectionsRecording = sum(Occupied),
              sunsetTimeSince = mean(as.numeric(sunsetTimeSince)),
              moon_fraction = mean(moon_fraction),
              moon_altitude = mean(moon_altitude),
              temp = mean(temp),
              wind_spd = mean(wind_spd),
              total_rain = mean(total_rain),
              moonUp = round(mean(moonUp)),
              band1_psd = mean(band1_psd),
              band1_s2n = mean(band1_s2n),
              band2_psd = mean(band2_psd),
              band2_s2n = mean(band2_s2n),
              Abundance_recording=mean(Abundance_recording),
              Abundance_site=mean(Abundance_site)) %>% 
    mutate(Occupied = ifelse(detectionsRecording > 0, 1, 0),
           moonUp_altitude = ifelse(moonUp==1, moon_altitude, 0)) %>% 
    dplyr::filter(!is.na(total_rain),
                  !is.na(band1_psd)) %>% 
    group_by(SiteName, Year) %>% 
    arrange(Recording) %>% 
    mutate(n=row_number()) %>% 
    ungroup()
  
  #Start for loop
  for(i in 1:boot){
    
    #Sample data
    dat.sub <- dat.h %>% 
      group_by(SiteName) %>% 
      sample_n(visit.h, replace=TRUE) %>%
      ungroup() %>% 
      unique() %>% 
      data.frame()
    
    #Format for occupancy  
    Occupied <- dat.sub %>% 
      dplyr::select(SiteName, n, Occupied) %>% 
      spread(key = n, value = Occupied) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    Abundance_site <- dat.sub %>% 
      dplyr::select(SiteName, n, Abundance_site) %>% 
      spread(key = n, value = Abundance_site) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    Abundance_recording <- dat.sub %>% 
      dplyr::select(SiteName, n, Abundance_recording) %>% 
      spread(key = n, value = Abundance_recording) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    Yday <- dat.sub %>% 
      dplyr::select(SiteName, n, Yday) %>% 
      spread(key = n, value = Yday) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    sunsetTimeSince <- dat.sub %>% 
      dplyr::select(SiteName, n, sunsetTimeSince) %>% 
      spread(key = n, value = sunsetTimeSince) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    moon_altitude <- dat.sub %>% 
      dplyr::select(SiteName, n, moon_altitude) %>% 
      spread(key = n, value = moon_altitude) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    temp <- dat.sub %>% 
      dplyr::select(SiteName, n, temp) %>% 
      spread(key = n, value = temp) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    total_rain <- dat.sub %>% 
      dplyr::select(SiteName, n, total_rain) %>% 
      spread(key = n, value = total_rain) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    wind_spd <- dat.sub %>% 
      dplyr::select(SiteName, n, wind_spd) %>% 
      spread(key = n, value = wind_spd) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    band2_psd <- dat.sub %>% 
      dplyr::select(SiteName, n, band2_psd) %>% 
      spread(key = n, value = band2_psd) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    band2_s2n <- dat.sub %>% 
      dplyr::select(SiteName, n, band2_s2n) %>% 
      spread(key = n, value = band2_s2n) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    obs.cov <- list(Yday=Yday, sunsetTimeSince=sunsetTimeSince, moon_altitude=moon_altitude, temp=temp, total_rain=total_rain, wind_spd=wind_spd, band2_psd=band2_psd, band2_s2n=band2_s2n, Abundance_site=Abundance_site, Abundance_recording=Abundance_recording)
    
    site.cov <- dat.sub %>% 
      dplyr::select(SiteName, Latitude, Longitude) %>% 
      unique() %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    dat.occ <- unmarkedFrameOccu(Occupied, siteCovs=site.cov, obsCovs=obs.cov)
    
    #Fit null model
    mod.any.null <- try(occu(~ 1
                             ~1,
                             data=dat.occ))
    
    #Predict on new data
    if(class(mod.any.null)[1]=="unmarkedFitOccu"){
      
      mod.any.null.pred.list1[[i]] <- dat.sub %>% 
        cbind(predict(mod.any.null, type="det", newdata=dat.sub)) %>% 
        dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
        cbind(predict(mod.any.null, type="state", newdata=dat.sub)) %>% 
        dplyr::rename(Occu=Predicted, OccuSE=SE, OccuLower=lower, OccuUpper=upper) %>% 
        unique() %>% 
        summarize(occu.mn.mn=mean(Occu, na.rm=TRUE),
                  occu.mn.sd=sd(Occu, na.rm=TRUE),
                  occu.lwr.mn=mean(OccuLower, na.rm=TRUE),
                  occu.upr.mn=mean(OccuUpper, na.rm=TRUE),
                  det.mn.mn=mean(Det, na.rm=TRUE),
                  det.mn.sd=sd(Det, na.rm=TRUE),
                  det.lwr.mn=mean(DetLower, na.rm=TRUE),
                  det.upr.mn=mean(DetUpper, na.rm=TRUE)) %>% 
        mutate(boot=i)
    }
    
    #Fit covariate model
    mod.any.cov <- try(occu(~ moon_altitude + sunsetTimeSince + I(sunsetTimeSince^2) + Yday + temp + I(temp^2) + band2_psd + band2_s2n
                            ~1,
                            data=dat.occ))
    
    #Predict on new data
    if(class(mod.any.cov)[1]=="unmarkedFitOccu"){
      
      mod.any.cov.pred.list1[[i]] <- dat.sub %>% 
        cbind(predict(mod.any.cov, type="det", newdata=dat.sub)) %>% 
        dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
        cbind(predict(mod.any.cov, type="state", newdata=dat.sub)) %>% 
        dplyr::rename(Occu=Predicted, OccuSE=SE, OccuLower=lower, OccuUpper=upper) %>% 
        unique() %>% 
        summarize(occu.mn.mn=mean(Occu, na.rm=TRUE),
                  occu.mn.sd=sd(Occu, na.rm=TRUE),
                  occu.lwr.mn=mean(OccuLower, na.rm=TRUE),
                  occu.upr.mn=mean(OccuUpper, na.rm=TRUE),
                  det.mn.mn=mean(Det, na.rm=TRUE),
                  det.mn.sd=sd(Det, na.rm=TRUE),
                  det.lwr.mn=mean(DetLower, na.rm=TRUE),
                  det.upr.mn=mean(DetUpper, na.rm=TRUE)) %>% 
        mutate(boot=i)
      
    }
    
    #Summarize detections by site
    summary.any <- dat.sub %>% 
      group_by(SiteName) %>% 
      summarize(Occupied_minutes = sum(Occupied),
                Occupied = ifelse(Occupied_minutes > 0, 1, 0)) %>% 
      ungroup() %>% 
      mutate(boot=i,
             length=length.h,
             visit=visit.h) %>% 
      rbind(summary.any)
    
    #Report status
    print(paste0("Finished bootstrap ", i, " of ", boot))
  }
  
  #Summarize occupancy & detectability estimates
  mod.any.null.pred.list2[[h]] <- rbindlist(mod.any.null.pred.list1) %>% 
    mutate(length=length.h,
           visit=visit.h)
  
  mod.any.cov.pred.list2[[h]] <- rbindlist(mod.any.cov.pred.list1) %>% 
    mutate(length=length.h,
           visit=visit.h)
  
  print(paste0("Finished protocol combination ", h, " of ", nrow(protocol)))
  
}

mod.any.null.pred <- rbindlist(mod.any.null.pred.list2) %>% 
  dplyr::filter(!is.na(occu.mn.mn)) %>% 
  group_by(length, visit) %>% 
  summarize(occu.mn.mn=mean(occu.mn.mn, na.rm=TRUE),
            occu.mn.sd=sd(occu.mn.mn, na.rm=TRUE),
            occu.lwr.mn=mean(occu.lwr.mn, na.rm=TRUE),
            occu.upr.mn=mean(occu.upr.mn, na.rm=TRUE),
            det.mn.mn=mean(det.mn.mn, na.rm=TRUE),
            det.mn.sd=sd(det.mn.mn, na.rm=TRUE),
            det.lwr.mn=mean(det.lwr.mn, na.rm=TRUE),
            det.upr.mn=mean(det.upr.mn, na.rm=TRUE),
            bootstraps=n()) %>% 
  ungroup() %>% 
  mutate(model="null") 

mod.any.cov.pred <- rbindlist(mod.any.cov.pred.list2) %>% 
  dplyr::filter(!is.na(occu.mn.mn)) %>% 
  group_by(length, visit) %>% 
  summarize(occu.mn.mn=mean(occu.mn.mn, na.rm=TRUE),
            occu.mn.sd=sd(occu.mn.mn, na.rm=TRUE),
            occu.lwr.mn=mean(occu.lwr.mn, na.rm=TRUE),
            occu.upr.mn=mean(occu.upr.mn, na.rm=TRUE),
            det.mn.mn=mean(det.mn.mn, na.rm=TRUE),
            det.mn.sd=sd(det.mn.mn, na.rm=TRUE),
            det.lwr.mn=mean(det.lwr.mn, na.rm=TRUE),
            det.upr.mn=mean(det.upr.mn, na.rm=TRUE),
            bootstraps=n()) %>% 
  ungroup() %>% 
  mutate(model="cov") 

mod.any.pred <- rbind(mod.any.cov.pred, mod.any.null.pred)

write.csv(mod.any.pred, "UnrestrictedModelResults.csv", row.names = FALSE)

#4b. Constrain covariates####

#Set number of bootstraps  
boot <- 100

#Set protocol options
length <- c(1:5)
#length <- 5
visit <- c(1:10)
#visit <- 20

protocol <- expand.grid(length=length, visit=visit)

#Create lists to save out results
mod.prot.null.pred.list1 <- list()
mod.prot.null.pred.list2 <- list()
mod.prot.cov.pred.list1 <- list()
mod.prot.cov.pred.list2 <- list()
summary.prot <- data.frame()

for(h in 1:nrow(protocol)){
  
  #Set protocol for this set
  length.h <- protocol$length[h]
  visit.h <- protocol$visit[h]
  
  #Wrangle data
  dat.h <- df %>% 
    inner_join(duration) %>% 
    dplyr::filter(Rec_minute <= length.h) %>% 
    dplyr::filter(moon_altitude >= 0.6,
                  temp <= 20,
                  temp >= 7,
                  sunsetTimeSince >= 0,
                  sunsetTimeSince <= 7,
                  wind_spd <= 19) %>% 
    mutate(moonUp = ifelse(moonUp=="Y", 1, 0)) %>% 
    group_by(SiteName, Latitude, Longitude, Year, Yday, Date, Recording, Duration, Time_rec) %>% 
    summarize(detectionsRecording = sum(Occupied),
              sunsetTimeSince = mean(as.numeric(sunsetTimeSince)),
              moon_fraction = mean(moon_fraction),
              moon_altitude = mean(moon_altitude),
              temp = mean(temp),
              wind_spd = mean(wind_spd),
              total_rain = mean(total_rain),
              moonUp = round(mean(moonUp)),
              band1_psd = mean(band1_psd),
              band1_s2n = mean(band1_s2n),
              band2_psd = mean(band2_psd),
              band2_s2n = mean(band2_s2n),
              Abundance_recording=mean(Abundance_recording),
              Abundance_site=mean(Abundance_site)) %>% 
    mutate(Occupied = ifelse(detectionsRecording > 0, 1, 0),
           moonUp_altitude = ifelse(moonUp==1, moon_altitude, 0)) %>% 
    dplyr::filter(!is.na(total_rain),
                  !is.na(band1_psd)) %>% 
    group_by(SiteName, Year) %>% 
    arrange(Recording) %>% 
    mutate(n=row_number()) %>% 
    ungroup()
  
  #Start for loop
  for(i in 1:boot){
    
    #Sample data
    dat.sub <- dat.h %>% 
      group_by(SiteName) %>% 
      sample_n(visit.h, replace=TRUE) %>%
      ungroup() %>% 
      unique() %>% 
      data.frame()
    
    #Format for occupancy  
    Occupied <- dat.sub %>% 
      dplyr::select(SiteName, n, Occupied) %>% 
      spread(key = n, value = Occupied) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    Abundance_site <- dat.sub %>% 
      dplyr::select(SiteName, n, Abundance_site) %>% 
      spread(key = n, value = Abundance_site) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    Abundance_recording <- dat.sub %>% 
      dplyr::select(SiteName, n, Abundance_recording) %>% 
      spread(key = n, value = Abundance_recording) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    Yday <- dat.sub %>% 
      dplyr::select(SiteName, n, Yday) %>% 
      spread(key = n, value = Yday) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    sunsetTimeSince <- dat.sub %>% 
      dplyr::select(SiteName, n, sunsetTimeSince) %>% 
      spread(key = n, value = sunsetTimeSince) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    moon_altitude <- dat.sub %>% 
      dplyr::select(SiteName, n, moon_altitude) %>% 
      spread(key = n, value = moon_altitude) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    temp <- dat.sub %>% 
      dplyr::select(SiteName, n, temp) %>% 
      spread(key = n, value = temp) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    total_rain <- dat.sub %>% 
      dplyr::select(SiteName, n, total_rain) %>% 
      spread(key = n, value = total_rain) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    wind_spd <- dat.sub %>% 
      dplyr::select(SiteName, n, wind_spd) %>% 
      spread(key = n, value = wind_spd) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    band2_psd <- dat.sub %>% 
      dplyr::select(SiteName, n, band2_psd) %>% 
      spread(key = n, value = band2_psd) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    band2_s2n <- dat.sub %>% 
      dplyr::select(SiteName, n, band2_s2n) %>% 
      spread(key = n, value = band2_s2n) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    obs.cov <- list(Yday=Yday, sunsetTimeSince=sunsetTimeSince, moon_altitude=moon_altitude, temp=temp, total_rain=total_rain, wind_spd=wind_spd, band2_psd=band2_psd, band2_s2n=band2_s2n, Abundance_site=Abundance_site, Abundance_recording=Abundance_recording)
    
    site.cov <- dat.sub %>% 
      dplyr::select(SiteName, Latitude, Longitude) %>% 
      unique() %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    dat.occ <- unmarkedFrameOccu(Occupied, siteCovs=site.cov, obsCovs=obs.cov)
    
    #Fit null model
    mod.prot.null <- try(occu(~ 1
                              ~1,
                              data=dat.occ))
    
    #Predict on new data
    if(class(mod.prot.null)[1]=="unmarkedFitOccu"){
      
      mod.prot.null.pred.list1[[i]] <- dat.sub %>% 
        cbind(predict(mod.prot.null, type="det", newdata=dat.sub)) %>% 
        dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
        cbind(predict(mod.prot.null, type="state", newdata=dat.sub)) %>% 
        dplyr::rename(Occu=Predicted, OccuSE=SE, OccuLower=lower, OccuUpper=upper) %>% 
        unique() %>% 
        summarize(occu.mn.mn=mean(Occu, na.rm=TRUE),
                  occu.mn.sd=sd(Occu, na.rm=TRUE),
                  occu.lwr.mn=mean(OccuLower, na.rm=TRUE),
                  occu.upr.mn=mean(OccuUpper, na.rm=TRUE),
                  det.mn.mn=mean(Det, na.rm=TRUE),
                  det.mn.sd=sd(Det, na.rm=TRUE),
                  det.lwr.mn=mean(DetLower, na.rm=TRUE),
                  det.upr.mn=mean(DetUpper, na.rm=TRUE)) %>% 
        mutate(boot=i)
    }
    
    #Fit covariate model
    mod.prot.cov <- try(occu(~ moon_altitude + sunsetTimeSince + I(sunsetTimeSince^2) + Yday + temp + I(temp^2) + band2_psd + band2_s2n
                             ~1,
                             data=dat.occ))
    
    #Predict on new data
    if(class(mod.prot.cov)[1]=="unmarkedFitOccu"){
      
      mod.prot.cov.pred.list1[[i]] <- dat.sub %>% 
        cbind(predict(mod.prot.cov, type="det", newdata=dat.sub)) %>% 
        dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
        cbind(predict(mod.prot.cov, type="state", newdata=dat.sub)) %>% 
        dplyr::rename(Occu=Predicted, OccuSE=SE, OccuLower=lower, OccuUpper=upper) %>% 
        unique() %>% 
        summarize(occu.mn.mn=mean(Occu, na.rm=TRUE),
                  occu.mn.sd=sd(Occu, na.rm=TRUE),
                  occu.lwr.mn=mean(OccuLower, na.rm=TRUE),
                  occu.upr.mn=mean(OccuUpper, na.rm=TRUE),
                  det.mn.mn=mean(Det, na.rm=TRUE),
                  det.mn.sd=sd(Det, na.rm=TRUE),
                  det.lwr.mn=mean(DetLower, na.rm=TRUE),
                  det.upr.mn=mean(DetUpper, na.rm=TRUE)) %>% 
        mutate(boot=i)
      
    }
    
    #Summarize detections by site
    summary.prot <- dat.sub %>% 
      group_by(SiteName) %>% 
      summarize(Occupied_minutes = sum(Occupied),
                Occupied = ifelse(Occupied_minutes > 0, 1, 0)) %>% 
      ungroup() %>% 
      mutate(boot=i,
             length=length.h,
             visit=visit.h) %>% 
      rbind(summary.prot)
    
    #Report status
    print(paste0("Finished bootstrap ", i, " of ", boot))
  }
  
  #Summarize occupancy & detectability estimates
  mod.prot.null.pred.list2[[h]] <- rbindlist(mod.prot.null.pred.list1) %>% 
    mutate(length=length.h,
           visit=visit.h)
  
  mod.prot.cov.pred.list2[[h]] <- rbindlist(mod.prot.cov.pred.list1) %>% 
    mutate(length=length.h,
           visit=visit.h)
  
  print(paste0("Finished protocol combination ", h, " of ", nrow(protocol)))
  
}

mod.prot.null.pred <- rbindlist(mod.prot.null.pred.list2) %>% 
  dplyr::filter(!is.na(occu.mn.mn)) %>% 
  group_by(length, visit) %>% 
  summarize(occu.mn.mn=mean(occu.mn.mn, na.rm=TRUE),
            occu.mn.sd=sd(occu.mn.mn, na.rm=TRUE),
            occu.lwr.mn=mean(occu.lwr.mn, na.rm=TRUE),
            occu.upr.mn=mean(occu.upr.mn, na.rm=TRUE),
            det.mn.mn=mean(det.mn.mn, na.rm=TRUE),
            det.mn.sd=sd(det.mn.mn, na.rm=TRUE),
            det.lwr.mn=mean(det.lwr.mn, na.rm=TRUE),
            det.upr.mn=mean(det.upr.mn, na.rm=TRUE),
            bootstraps=n()) %>% 
  ungroup() %>% 
  mutate(model="null")

mod.prot.cov.pred <- rbindlist(mod.prot.cov.pred.list2) %>% 
  dplyr::filter(!is.na(occu.mn.mn)) %>% 
  group_by(length, visit) %>% 
  summarize(occu.mn.mn=mean(occu.mn.mn, na.rm=TRUE),
            occu.mn.sd=sd(occu.mn.mn, na.rm=TRUE),
            occu.lwr.mn=mean(occu.lwr.mn, na.rm=TRUE),
            occu.upr.mn=mean(occu.upr.mn, na.rm=TRUE),
            det.mn.mn=mean(det.mn.mn, na.rm=TRUE),
            det.mn.sd=sd(det.mn.mn, na.rm=TRUE),
            det.lwr.mn=mean(det.lwr.mn, na.rm=TRUE),
            det.upr.mn=mean(det.upr.mn, na.rm=TRUE),
            bootstraps=n()) %>% 
  ungroup() %>% 
  mutate(model="cov")

mod.prot.pred <- rbind(mod.prot.cov.pred, mod.prot.null.pred) %>% 
  mutate(cov="constrained") %>% 
  dplyr::select(cov, model, length, visit, occu.mn.mn, occu.lwr.mn, occu.upr.mn, det.mn.mn, det.lwr.mn, det.upr.mn)
mod.any.pred <- rbind(mod.any.cov.pred, mod.any.null.pred) %>% 
  mutate(cov="all") %>% 
  dplyr::select(cov, model, length, visit, occu.mn.mn, occu.lwr.mn, occu.upr.mn, det.mn.mn, det.lwr.mn, det.upr.mn)

mod.pred <- rbind(mod.prot.pred, mod.any.pred)

write.csv(mod.pred, "SamplingEffortResults.csv", row.names = FALSE)

#Summary results
summary <- rbind(summary.any %>% 
                   mutate(cov="all"),
                 summary.prot %>% 
                   mutate(cov="constrained")) %>% 
  mutate(DetsPerVisit = Occupied_minutes/visit,
         DetsPerMinute=DetsPerVisit/length)

#ggplot(summary, aes(x=visit, y=DetsPerMinute, colour=factor(length), linetype=cov)) +
#  geom_point() +
#  geom_smooth()

summary.sum <- summary %>% 
  group_by(cov, length, visit, boot) %>% 
  summarize(sites = sum(Occupied)/32) %>% 
  ungroup() %>% 
  mutate(minutes = length*visit)

summary.sum.any <- summary.sum %>% 
  dplyr::filter(cov=="all")

ggplot(summary.sum, aes(x=visit, y=sites, colour=factor(length))) +
  geom_jitter() +
  geom_smooth(method="loess") +
  facet_wrap(~cov, scales="free_x") +
  xlim(c(0,10))

ggplot(summary.sum, aes(x=visit, y=sites, colour=factor(length), linetype=cov)) +
  #  geom_jitter() +
  geom_smooth(method="loess") +
  xlim(c(0,10))

write.csv(summary, "SamplingEffortResults_Summary.csv", row.names = FALSE)
write.csv(summary.sum.any, "NLSData.csv", row.names = FALSE)

#4c. Model cumulative probability of detection----

#4ci. Logistic glmm----

#Compare constrained and unconstrained
summary.10 <- summary %>% 
  dplyr::filter(visit<=10)

table(summary.10$cov) #Good, equal samples

mod.glmm <- glmer(Occupied ~ visit*length + visit*cov + length*cov + (1|SiteName), data=summary.10, family="binomial")
summary(mod.glmm)

newdat <- expand.grid(length = c(1:5), visit = c(1:10), cov=c("all", "constrained"),
                      SiteName=unique(summary.10$SiteName))

pred <- predictInterval(mod.glmm, newdata = newdat, which="fixed", n.sims = 1000, type="probability", level=0.95) %>% 
  cbind(newdat)

pred.glmm <- data.frame(pred=predict(mod.glmm, newdata = newdat, type="response", re.form=NA)) %>% 
  cbind(newdat) %>% 
  cbind(predictInterval(mod.glmm, newdata = newdat, which="fixed", n.sims = 1000, type="probability", level=0.95))

ggplot(pred.glmm) +
  geom_ribbon(aes(x=visit, ymin=lwr, ymax=upr, group=factor(length)), alpha=0.2) +
  geom_line(aes(x=visit, y=pred, colour=factor(length), linetype=cov)) +
  facet_wrap(~cov)

write.csv(pred.glmm, "GLMMPredictions.csv", row.names = FALSE)

#4cii. Logistic growth curve for unconstrained by # visits----
length <- c(1:5)
visit <- c(1:30, 40, 50, 60, 70, 80, 90, 100)

mod.cumu.any.visit.list <- list()
pred.any.visit.list <- list()
pred.asym.list <- list()
for(i in 1:length(length)){
  
  length.i <- length[i]
  
  summary.sum.any.i <- summary.sum.any %>% 
    dplyr::filter(length==length.i)
  
  mod.cumu.any.i <- nls(sites ~ SSasymp(visit, Asym, R0, lrc), data = summary.sum.any.i)
  #mod.cumu.any.i <- drm(sites ~ visit, fct = DRC.asymReg(), data = summary.sum.any.i)
  
  mod.cumu.any.visit.list[[i]] <- data.frame(summary(mod.cumu.any.i)$coefficients) %>% 
    mutate(length=length[i],
           var=row.names(summary(mod.cumu.any.i)$coefficients))
  
  #Fit growth curve to new data
  newdat <- data.frame(visit=seq(min(summary.sum.any.i$visit), max(summary.sum.any.i$visit),
                                 length.out = 100))
  
  pred.any.visit.list[[i]] <- data.frame(
    r=predict(newdata = newdat, object = mod.cumu.any.i),
    visit=newdat$visit) %>% 
    mutate(length=length[i])
  
  #Find asymptote
  asym <- round(environment(mod.cumu.any.i[["m"]][["fitted"]])[["env"]][["Asym"]], 3)
  ls.sum.99 <- pred.any.visit.list[[i]] %>%
    mutate(r = round(r, digits = 3)) %>% filter(r >= 0.99 * asym)
  ls.sum.95 <- pred.any.visit.list[[i]] %>%
    mutate(r = round(r, digits = 3)) %>% filter(r >= 0.95 * asym)
  pred.asym.list[[i]] <- data.frame(asym=asym,
                                    n = round(min(ls.sum.99$visit)),
                                    length=length[i])
  
}

mod.cumu.any.visit <- rbindlist(mod.cumu.any.visit.list)
pred.any.visit <- rbindlist(pred.any.visit.list)
pred.asym <- rbindlist(pred.asym.list)

ggplot() +
  geom_jitter(data=summary.sum.any, aes(x=visit, y=sites, colour=factor(length))) +
  geom_line(data=pred.any.visit, aes(x=visit, y=r, colour=factor(length))) +
  geom_vline(data=pred.asym, aes(xintercept=n99, colour=factor(length)), linetype="dashed") +
  geom_vline(data=pred.asym, aes(xintercept=n95, colour=factor(length)), linetype="dotted") +
  ylim(c(0,1)) + 
#  facet_wrap(~length) +
  scale_fill_viridis_c()

pred.asym.visit <- pred.asym %>% 
  mutate(minutes = n*length)

write.csv(pred.any.visit, "NLSPredictions_Visits.csv", row.names = FALSE)
write.csv(pred.asym.visit, "NLSAsymptotes_Visits.csv", row.names = FALSE)

#4ciii. Logistic growth curve for unconstrained by recording length----
length <- c(1:5)
visit <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

#For each recording length
mod.cumu.any.length.list <- list()
pred.any.length.list <- list()
pred.asym.list <- list()
for(i in 1:length(visit)){
  
  visit.i <- visit[i]
  
  summary.sum.any.i <- summary.sum.any %>% 
    dplyr::filter(visit==visit.i)
  
  mod.cumu.any.i <- nls(sites ~ SSasymp(length, Asym, R0, lrc), data = summary.sum.any.i)
  mod.cumu.any.length.list[[i]] <- data.frame(summary(mod.cumu.any.i)$coefficients) %>% 
    mutate(visit=visit[i],
           var=row.names(summary(mod.cumu.any.i)$coefficients))
  
  #Fit growth curve to new data
  newdat <- data.frame(length=seq(min(summary.sum.any.i$length), max(summary.sum.any.i$length),
                                   length.out = 100))
  
  pred.any.length.list[[i]] <- data.frame(
    r=predict(newdata = newdat, object = mod.cumu.any.i),
    length=newdat$length) %>% 
    mutate(visit=visit[i])
  
  #Find asymptote
  asym <- round(environment(mod.cumu.any.i[["m"]][["fitted"]])[["env"]][["Asym"]], 3)
  ls.sum.99 <- pred.any.length.list[[i]] %>%
    mutate(r = round(r, digits = 3)) %>% filter(r >= 0.99 * asym)
  ls.sum.95 <- pred.any.length.list[[i]] %>%
    mutate(r = round(r, digits = 3)) %>% filter(r >= 0.95 * asym)
  pred.asym.list[[i]] <- data.frame(asym=asym,
                                    n = round(min(ls.sum.99$length)),
                                    visit=visit[i])
  
}

mod.cumu.any.length <- rbindlist(mod.cumu.any.length.list)
pred.any.length <- rbindlist(pred.any.length.list)
pred.asym <- rbindlist(pred.asym.list)

ggplot() +
  geom_jitter(data=summary.sum.any, aes(x=length, y=sites, group=factor(visit))) +
  geom_line(data=pred.any.length, aes(x=length, y=r, colour=factor(visit))) +
  geom_vline(data=pred.asym, aes(xintercept=n99, colour=factor(visit)), linetype="dashed") +
  ylim(c(0,1))

pred.asym.length <- pred.asym %>% 
  mutate(visit95 = n/length)

write.csv(pred.any.length, "NLSPredictions_Length.csv", row.names = FALSE)
write.csv(pred.asym.length, "NLSAsymptotes_Length.csv", row.names = FALSE)

#SECTION 5. Recording length----
#5a. Model----

table(duration$SiteName, duration$Duration)

#Set number of bootstraps  
boot <- 100

#Set protocol options
length <- c(1:10)
visit <- seq(1, 100, 10)

protocol <- expand.grid(length=length, visit=visit)

#Create lists to save out results
mod.10.null.pred.list1 <- list()
mod.10.null.pred.list2 <- list()
mod.10.cov.pred.list1 <- list()
mod.10.cov.pred.list2 <- list()
summary.10 <- data.frame()

for(h in 1:nrow(protocol)){
  
  #Set protocol for this set
  length.h <- protocol$length[h]
  visit.h <- protocol$visit[h]
  
  #Wrangle data
  dat.h <- df %>% 
    inner_join(duration) %>% 
    dplyr::filter(str_sub(SiteName, 1, 3)!="PEP", Rec_minute <= length.h) %>% 
    mutate(moonUp = ifelse(moonUp=="Y", 1, 0)) %>% 
    group_by(SiteName, Latitude, Longitude, Year, Yday, Date, Recording, Duration, Time_rec) %>% 
    summarize(detectionsRecording = sum(Occupied),
              sunsetTimeSince = mean(as.numeric(sunsetTimeSince)),
              moon_fraction = mean(moon_fraction),
              moon_altitude = mean(moon_altitude),
              temp = mean(temp),
              wind_spd = mean(wind_spd),
              total_rain = mean(total_rain),
              moonUp = round(mean(moonUp)),
              band1_psd = mean(band1_psd),
              band1_s2n = mean(band1_s2n),
              band2_psd = mean(band2_psd),
              band2_s2n = mean(band2_s2n),
              Abundance_recording=mean(Abundance_recording),
              Abundance_site=mean(Abundance_site)) %>% 
    mutate(Occupied = ifelse(detectionsRecording > 0, 1, 0),
           moonUp_altitude = ifelse(moonUp==1, moon_altitude, 0)) %>% 
    dplyr::filter(!is.na(total_rain),
                  !is.na(band1_psd)) %>% 
    group_by(SiteName, Year) %>% 
    arrange(Recording) %>% 
    mutate(n=row_number()) %>% 
    ungroup()
  
  #Start for loop
  for(i in 1:boot){
    
    #Sample data
    dat.sub <- dat.h %>% 
      group_by(SiteName) %>% 
      sample_n(visit.h, replace=TRUE) %>%
      ungroup() %>% 
      unique() %>% 
      data.frame()
    
    #Format for occupancy  
    Occupied <- dat.sub %>% 
      dplyr::select(SiteName, n, Occupied) %>% 
      spread(key = n, value = Occupied) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    Abundance_site <- dat.sub %>% 
      dplyr::select(SiteName, n, Abundance_site) %>% 
      spread(key = n, value = Abundance_site) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    Abundance_recording <- dat.sub %>% 
      dplyr::select(SiteName, n, Abundance_recording) %>% 
      spread(key = n, value = Abundance_recording) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    Yday <- dat.sub %>% 
      dplyr::select(SiteName, n, Yday) %>% 
      spread(key = n, value = Yday) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    sunsetTimeSince <- dat.sub %>% 
      dplyr::select(SiteName, n, sunsetTimeSince) %>% 
      spread(key = n, value = sunsetTimeSince) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    moon_altitude <- dat.sub %>% 
      dplyr::select(SiteName, n, moon_altitude) %>% 
      spread(key = n, value = moon_altitude) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    temp <- dat.sub %>% 
      dplyr::select(SiteName, n, temp) %>% 
      spread(key = n, value = temp) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    total_rain <- dat.sub %>% 
      dplyr::select(SiteName, n, total_rain) %>% 
      spread(key = n, value = total_rain) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    wind_spd <- dat.sub %>% 
      dplyr::select(SiteName, n, wind_spd) %>% 
      spread(key = n, value = wind_spd) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    band2_psd <- dat.sub %>% 
      dplyr::select(SiteName, n, band2_psd) %>% 
      spread(key = n, value = band2_psd) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    band2_s2n <- dat.sub %>% 
      dplyr::select(SiteName, n, band2_s2n) %>% 
      spread(key = n, value = band2_s2n) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    obs.cov <- list(Yday=Yday, sunsetTimeSince=sunsetTimeSince, moon_altitude=moon_altitude, temp=temp, total_rain=total_rain, wind_spd=wind_spd, band2_psd=band2_psd, band2_s2n=band2_s2n, Abundance_site=Abundance_site, Abundance_recording=Abundance_recording)
    
    site.cov <- dat.sub %>% 
      dplyr::select(SiteName, Latitude, Longitude) %>% 
      unique() %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    dat.occ <- unmarkedFrameOccu(Occupied, siteCovs=site.cov, obsCovs=obs.cov)
    
    #Fit null model
    mod.10.null <- try(occu(~ 1
                             ~1,
                             data=dat.occ))
    
    #Predict on new data
    if(class(mod.10.null)[1]=="unmarkedFitOccu"){
      
      mod.10.null.pred.list1[[i]] <- dat.sub %>% 
        cbind(predict(mod.10.null, type="det", newdata=dat.sub)) %>% 
        dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
        cbind(predict(mod.10.null, type="state", newdata=dat.sub)) %>% 
        dplyr::rename(Occu=Predicted, OccuSE=SE, OccuLower=lower, OccuUpper=upper) %>% 
        unique() %>% 
        summarize(occu.mn.mn=mean(Occu, na.rm=TRUE),
                  occu.mn.sd=sd(Occu, na.rm=TRUE),
                  occu.lwr.mn=mean(OccuLower, na.rm=TRUE),
                  occu.upr.mn=mean(OccuUpper, na.rm=TRUE),
                  det.mn.mn=mean(Det, na.rm=TRUE),
                  det.mn.sd=sd(Det, na.rm=TRUE),
                  det.lwr.mn=mean(DetLower, na.rm=TRUE),
                  det.upr.mn=mean(DetUpper, na.rm=TRUE)) %>% 
        mutate(boot=i)
    }
    
    #Fit covariate model
    mod.10.cov <- try(occu(~ moon_altitude + sunsetTimeSince + I(sunsetTimeSince^2) + Yday + temp + I(temp^2) + band2_psd + band2_s2n
                            ~1,
                            data=dat.occ))
    
    #Predict on new data
    if(class(mod.10.cov)[1]=="unmarkedFitOccu"){
      
      mod.10.cov.pred.list1[[i]] <- dat.sub %>% 
        cbind(predict(mod.10.cov, type="det", newdata=dat.sub)) %>% 
        dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
        cbind(predict(mod.10.cov, type="state", newdata=dat.sub)) %>% 
        dplyr::rename(Occu=Predicted, OccuSE=SE, OccuLower=lower, OccuUpper=upper) %>% 
        unique() %>% 
        summarize(occu.mn.mn=mean(Occu, na.rm=TRUE),
                  occu.mn.sd=sd(Occu, na.rm=TRUE),
                  occu.lwr.mn=mean(OccuLower, na.rm=TRUE),
                  occu.upr.mn=mean(OccuUpper, na.rm=TRUE),
                  det.mn.mn=mean(Det, na.rm=TRUE),
                  det.mn.sd=sd(Det, na.rm=TRUE),
                  det.lwr.mn=mean(DetLower, na.rm=TRUE),
                  det.upr.mn=mean(DetUpper, na.rm=TRUE)) %>% 
        mutate(boot=i)
      
    }
    
    #Summarize detections by site
    summary.10 <- dat.sub %>% 
      group_by(SiteName) %>% 
      summarize(Occupied_minutes = sum(Occupied),
                Occupied = ifelse(Occupied_minutes > 0, 1, 0)) %>% 
      ungroup() %>% 
      mutate(boot=i,
             length=length.h,
             visit=visit.h) %>% 
      rbind(summary.10)
    
    #Report status
    print(paste0("Finished bootstrap ", i, " of ", boot))
  }
  
  #Summarize occupancy & detectability estimates
  mod.10.null.pred.list2[[h]] <- rbindlist(mod.10.null.pred.list1) %>% 
    mutate(length=length.h,
           visit=visit.h)
  
  mod.10.cov.pred.list2[[h]] <- rbindlist(mod.10.cov.pred.list1) %>% 
    mutate(length=length.h,
           visit=visit.h)
  
  print(paste0("Finished protocol combination ", h, " of ", nrow(protocol)))
  
}

mod.10.null.pred <- rbindlist(mod.10.null.pred.list2) %>% 
  dplyr::filter(!is.na(occu.mn.mn)) %>% 
  group_by(length, visit) %>% 
  summarize(occu.mn.mn=mean(occu.mn.mn, na.rm=TRUE),
            occu.mn.sd=sd(occu.mn.mn, na.rm=TRUE),
            occu.lwr.mn=mean(occu.lwr.mn, na.rm=TRUE),
            occu.upr.mn=mean(occu.upr.mn, na.rm=TRUE),
            det.mn.mn=mean(det.mn.mn, na.rm=TRUE),
            det.mn.sd=sd(det.mn.mn, na.rm=TRUE),
            det.lwr.mn=mean(det.lwr.mn, na.rm=TRUE),
            det.upr.mn=mean(det.upr.mn, na.rm=TRUE),
            bootstraps=n()) %>% 
  ungroup() %>% 
  mutate(model="null") 

mod.10.cov.pred <- rbindlist(mod.10.cov.pred.list2) %>% 
  dplyr::filter(!is.na(occu.mn.mn)) %>% 
  group_by(length, visit) %>% 
  summarize(occu.mn.mn=mean(occu.mn.mn, na.rm=TRUE),
            occu.mn.sd=sd(occu.mn.mn, na.rm=TRUE),
            occu.lwr.mn=mean(occu.lwr.mn, na.rm=TRUE),
            occu.upr.mn=mean(occu.upr.mn, na.rm=TRUE),
            det.mn.mn=mean(det.mn.mn, na.rm=TRUE),
            det.mn.sd=sd(det.mn.mn, na.rm=TRUE),
            det.lwr.mn=mean(det.lwr.mn, na.rm=TRUE),
            det.upr.mn=mean(det.upr.mn, na.rm=TRUE),
            bootstraps=n()) %>% 
  ungroup() %>% 
  mutate(model="cov") 

mod.10.pred <- rbind(mod.10.cov.pred, mod.10.null.pred)

write.csv(mod.10.pred, "RecordingLengthModelResults.csv", row.names = FALSE)

mod.10.pred.boot <- rbindlist(mod.10.null.pred.list2) %>% 
  dplyr::filter(!is.na(occu.mn.mn)) %>% 
  mutate(minutes = length*visit)

#Visualize
ggplot(mod.10.pred) +
  geom_line(aes(x=length, y=occu.mn.mn, colour=factor(visit))) +
  facet_wrap(~model)

ggplot(mod.10.pred) +
  geom_line(aes(x=length, y=det.mn.mn, colour=factor(visit))) +
  facet_wrap(~model)

ggplot(mod.10.pred.boot) +
#  geom_point(aes(x=length, y=occu.mn.mn, colour=factor(visit))) +
  geom_smooth(aes(x=length, y=occu.mn.mn, colour=factor(visit)))

ggplot(mod.10.pred.boot) +
  #  geom_point(aes(x=length, y=occu.mn.mn, colour=factor(visit))) +
  geom_smooth(aes(x=length, y=det.mn.mn, colour=factor(visit)))

ggplot(mod.10.pred.boot) +
  geom_point(aes(x=minutes, y=occu.mn.mn)) +
  geom_smooth(aes(x=minutes, y=occu.mn.mn, colour=factor(length)))

ggplot(mod.10.pred.boot) +
  geom_point(aes(x=minutes, y=det.mn.mn)) +
  geom_smooth(aes(x=minutes, y=det.mn.mn, colour=factor(length)))

#5b. NLS for recording length----
summary.sum.10 <- summary.10 %>% 
  group_by(length, visit, boot) %>% 
  summarize(sites = sum(Occupied)/12) %>% 
  ungroup()

length <- c(1:10)
visit <- seq(1, 100, 10)

mod.10.visit.list <- list()
pred.10.visit.list <- list()
pred.10.asym.list <- list()
for(i in 1:length(visit)){
  
  visit.i <- visit[i]
  
  dat.i <- summary.sum.10 %>% 
    dplyr::filter(visit==visit.i)
  
  mod.10.i <- nls(sites ~ SSlogis(length, Asym, xmid, scal), data = dat.i)
  mod.10.visit.list[[i]] <- data.frame(summary(mod.10.i)$coefficients) %>% 
    mutate(length=length[i],
           var=row.names(summary(mod.10.i)$coefficients))
  
  #Fit growth curve to new data
  newdat <- data.frame(length=seq(min(dat.i$length), max(dat.i$length),
                                 length.out = 100))
  
  pred.10.visit.list[[i]] <- data.frame(
    r=predict(newdata = newdat, object = mod.10.i),
    length=newdat$length) %>% 
    mutate(visit=visit[i])
  
  #Find asymptote
  asym <- round(environment(mod.10.i[["m"]][["fitted"]])[["env"]][["Asym"]], 3)
  ls.sum <- pred.10.visit.list[[i]] %>%
    mutate(r = round(r, digits = 3)) %>% filter(r >= 0.99 * asym)
  pred.10.asym.list[[i]] <- data.frame(asym=asym,
                                    n = round(min(ls.sum$length)),
                                    visit=visit[i])
  
}

mod.10.visit <- rbindlist(mod.10.visit.list)
pred.10.visit <- rbindlist(pred.10.visit.list)
pred.10.asym <- rbindlist(pred.10.asym.list)

ggplot() +
  geom_jitter(data=summary.sum.10, aes(x=length, y=sites, colour=factor(visit))) +
  geom_line(data=pred.10.visit, aes(x=length, y=r, colour=factor(visit))) +
  geom_vline(data=pred.10.asym, aes(xintercept=n, colour=factor(visit)), linetype="dashed") +
  ylim(c(0,1))

write.csv(pred.10.visit, "NLSPredictions_Lengths.csv", row.names = FALSE)
write.csv(summary.sum.10, "NLSData_Lengths.csv", row.names = FALSE)
write.csv(pred.10.asym, "NLSAsymptotes_Lengths.csv", row.names = FALSE)

#Conclusion: detectability does not asymptote within 10 minutes

#SECTION 6. Canadian Nightjar Survey analysis----

#6a. Occupancy models-----
#Set number of bootstraps  
boot <- 100

#Set protocol options
length <- c(1:5)
visit <- c(1)

protocol <- expand.grid(length=length, visit=visit)

#Create lists to save out results
mod.cns.pred.list1 <- list()
mod.cns.pred.list2 <- list()
summary.cns <- data.frame()

for(h in 1:nrow(protocol)){
  
  #Set protocol for this set
  length.h <- protocol$length[h]
  visit.h <- protocol$visit[h]
  
  #Wrangle data
  dat.h <- df %>% 
    inner_join(duration) %>% 
    dplyr::filter(Rec_minute <= length.h) %>% 
    dplyr::filter(moon_altitude >= 0.6,
                  temp <= 20,
                  temp >= 7,
                  sunsetTimeSince >= 0,
                  sunsetTimeSince <= 7,
                  wind_spd <= 19) %>% 
    mutate(moonUp = ifelse(moonUp=="Y", 1, 0)) %>% 
    group_by(SiteName, Latitude, Longitude, Year, Yday, Date, Recording, Duration, Time_rec) %>% 
    summarize(detectionsRecording = sum(Occupied),
              sunsetTimeSince = mean(as.numeric(sunsetTimeSince)),
              moon_fraction = mean(moon_fraction),
              moon_altitude = mean(moon_altitude),
              temp = mean(temp),
              wind_spd = mean(wind_spd),
              total_rain = mean(total_rain),
              moonUp = round(mean(moonUp)),
              band1_psd = mean(band1_psd),
              band1_s2n = mean(band1_s2n),
              band2_psd = mean(band2_psd),
              band2_s2n = mean(band2_s2n),
              Abundance_recording=mean(Abundance_recording),
              Abundance_site=mean(Abundance_site)) %>% 
    mutate(Occupied = ifelse(detectionsRecording > 0, 1, 0),
           moonUp_altitude = ifelse(moonUp==1, moon_altitude, 0)) %>% 
    dplyr::filter(!is.na(total_rain),
                  !is.na(band1_psd)) %>% 
    group_by(SiteName, Year) %>% 
    arrange(Recording) %>% 
    mutate(n=row_number()) %>% 
    ungroup()
  
  #Start for loop
  for(i in 1:boot){
    
    #Sample data
    dat.sub <- dat.h %>% 
      group_by(SiteName) %>% 
      sample_n(visit.h, replace=TRUE) %>%
      ungroup() %>% 
      unique() %>% 
      data.frame()
    
    #Format for occupancy  
    Occupied <- dat.sub %>% 
      dplyr::select(SiteName, n, Occupied) %>% 
      spread(key = n, value = Occupied) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    Abundance_site <- dat.sub %>% 
      dplyr::select(SiteName, n, Abundance_site) %>% 
      spread(key = n, value = Abundance_site) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    Abundance_recording <- dat.sub %>% 
      dplyr::select(SiteName, n, Abundance_recording) %>% 
      spread(key = n, value = Abundance_recording) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    Yday <- dat.sub %>% 
      dplyr::select(SiteName, n, Yday) %>% 
      spread(key = n, value = Yday) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    sunsetTimeSince <- dat.sub %>% 
      dplyr::select(SiteName, n, sunsetTimeSince) %>% 
      spread(key = n, value = sunsetTimeSince) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    moon_altitude <- dat.sub %>% 
      dplyr::select(SiteName, n, moon_altitude) %>% 
      spread(key = n, value = moon_altitude) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    temp <- dat.sub %>% 
      dplyr::select(SiteName, n, temp) %>% 
      spread(key = n, value = temp) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    total_rain <- dat.sub %>% 
      dplyr::select(SiteName, n, total_rain) %>% 
      spread(key = n, value = total_rain) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    wind_spd <- dat.sub %>% 
      dplyr::select(SiteName, n, wind_spd) %>% 
      spread(key = n, value = wind_spd) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    band2_psd <- dat.sub %>% 
      dplyr::select(SiteName, n, band2_psd) %>% 
      spread(key = n, value = band2_psd) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    band2_s2n <- dat.sub %>% 
      dplyr::select(SiteName, n, band2_s2n) %>% 
      spread(key = n, value = band2_s2n) %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    obs.cov <- list(Yday=Yday, sunsetTimeSince=sunsetTimeSince, moon_altitude=moon_altitude, temp=temp, total_rain=total_rain, wind_spd=wind_spd, band2_psd=band2_psd, band2_s2n=band2_s2n, Abundance_site=Abundance_site, Abundance_recording=Abundance_recording)
    
    site.cov <- dat.sub %>% 
      dplyr::select(SiteName, Latitude, Longitude) %>% 
      unique() %>% 
      arrange(SiteName) %>% 
      dplyr::select(-SiteName) %>% 
      data.frame()
    
    dat.occ <- unmarkedFrameOccu(Occupied, siteCovs=site.cov, obsCovs=obs.cov)
    
    #Fit covariate model
    mod.prot.cov <- try(occu(~ moon_altitude + sunsetTimeSince + I(sunsetTimeSince^2) + Yday + temp + I(temp^2)
                             ~1,
                             data=dat.occ))
    
    #Predict on new data
    if(class(mod.prot.cov)[1]=="unmarkedFitOccu"){
      
      mod.cns.pred.list1[[i]] <- dat.sub %>% 
        cbind(predict(mod.prot.cov, type="det", newdata=dat.sub)) %>% 
        dplyr::rename(Det=Predicted, DetSE=SE, DetLower=lower, DetUpper=upper) %>% 
        cbind(predict(mod.prot.cov, type="state", newdata=dat.sub)) %>% 
        dplyr::rename(Occu=Predicted, OccuSE=SE, OccuLower=lower, OccuUpper=upper) %>% 
        unique() %>% 
        summarize(occu.mn.mn=mean(Occu, na.rm=TRUE),
                  occu.mn.sd=sd(Occu, na.rm=TRUE),
                  occu.lwr.mn=mean(OccuLower, na.rm=TRUE),
                  occu.upr.mn=mean(OccuUpper, na.rm=TRUE),
                  det.mn.mn=mean(Det, na.rm=TRUE),
                  det.mn.sd=sd(Det, na.rm=TRUE),
                  det.lwr.mn=mean(DetLower, na.rm=TRUE),
                  det.upr.mn=mean(DetUpper, na.rm=TRUE)) %>% 
        mutate(boot=i)
      
    }
    
    #Summarize detections by site
    summary.cns <- dat.sub %>% 
      group_by(SiteName) %>% 
      summarize(Occupied_minutes = sum(Occupied),
                Occupied = ifelse(Occupied_minutes > 0, 1, 0)) %>% 
      ungroup() %>% 
      mutate(boot=i,
             length=length.h,
             visit=visit.h) %>% 
      rbind(summary.cns)
    
    #Report status
    print(paste0("Finished bootstrap ", i, " of ", boot))
  }
  
  #Summarize occupancy & detectability estimates
  mod.cns.pred.list2[[h]] <- rbindlist(mod.cns.pred.list1) %>% 
    mutate(length=length.h,
           visit=visit.h)
  
  print(paste0("Finished protocol combination ", h, " of ", nrow(protocol)))
  
}

mod.cns.pred <- rbindlist(mod.cns.pred.list2) %>% 
  dplyr::filter(!is.na(occu.mn.mn)) %>% 
  group_by(length, visit) %>% 
  summarize(occu.mn.mn=mean(occu.mn.mn, na.rm=TRUE),
            occu.mn.sd=sd(occu.mn.mn, na.rm=TRUE),
            occu.lwr.mn=mean(occu.lwr.mn, na.rm=TRUE),
            occu.upr.mn=mean(occu.upr.mn, na.rm=TRUE),
            det.mn.mn=mean(det.mn.mn, na.rm=TRUE),
            det.mn.sd=sd(det.mn.mn, na.rm=TRUE),
            det.lwr.mn=mean(det.lwr.mn, na.rm=TRUE),
            det.upr.mn=mean(det.upr.mn, na.rm=TRUE),
            bootstraps=n()) %>% 
  ungroup()
View(mod.cns.pred)

ggplot(mod.cns.pred) +
  geom_ribbon(aes(x=length, ymin=det.lwr.mn, ymax=det.upr.mn), alpha=0.3) +
  geom_line(aes(x=length, y=det.mn.mn), size=1.5) +
  xlab("") +
  ylab("Probability of detection") +
  #  xlim(c(1,10)) +
  ylim(c(0,1)) +
  my.theme

ggplot(mod.cns.pred) +
  geom_ribbon(aes(x=length, ymin=occu.lwr.mn, ymax=occu.upr.mn), alpha=0.3) +
  geom_line(aes(x=length, y=occu.mn.mn), size=1.5) +
  xlab("") +
  ylab("Probability of occupancy") +
  #  xlim(c(1,10)) +
  ylim(c(0,1)) +
  my.theme

summary.cns.sum <- summary.cns %>% 
  group_by(length, visit, boot) %>% 
  summarize(sites = sum(Occupied)/32) %>% 
  ungroup()

ggplot(summary.cns.sum, aes(x=length, y=sites)) +
  geom_jitter() +
  geom_smooth(method="loess")

#6b: Correction factors----
summary.cns.sum <- summary.sum %>% 
  dplyr::filter(cov=="constrained",
                visit==1)

ggplot(summary.cns.sum, aes(x=length, y=sites)) +
  geom_jitter() +
  geom_smooth(method="loess")


#TO DO: FIGURE OUT OFFSETS

mod.cns <- lm(sites ~ length, data=summary.cns.sum)
pred.cns <- data.frame(predict(mod.cns, newdata=data.frame(length=6), se.fit=TRUE)) %>% 
  mutate(offset.rec = 20/17,
         offset.aru = (1/0.75)^2,
         p.adj = fit*offset.rec*offset.aru,
         lwr.adj = p.adj - se.fit*1.96*offset.rec*offset.aru,
         upr.adj = p.adj + se.fit*1.96*offset.rec*offset.aru)
pred.cns

save.image("Analysisv7_workspace.Rdata")


#APPENDIX 3: Loudness analysis----
dat.test.1 <- dat.test %>% 
  dplyr::filter(ewpw==1) %>% 
  rename(level=LevelNo)

lm.rsl.1 <- lmer(score ~ rsl + (1|file), data=dat.test.1)
lm.rsl.0 <- lmer(score ~ 1 + (1|file), data=dat.test.1)
aictab(list(lm.rsl.1, lm.rsl.0))

lm.level.1 <- lmer(score ~ level + (1|file), data=dat.test.1)
lm.level.0 <- lmer(score ~ 1 + (1|file), data=dat.test.1)
aictab(list(lm.level.1, lm.level.0))
