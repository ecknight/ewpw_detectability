#' ---
#' title: "Refining surveys for nocturnal birds using acoustic data: Eastern Whip-poor-will interim report"
#' author: "Elly C. Knight"
#' date: "December 17, 2020"
#' output: pdf_document
#' ---
#' 
#+ echo=FALSE
options(scipen = 9999)
knitr::opts_chunk$set(tidy.opts=list(width.cutoff=60), tidy=TRUE)
knitr::opts_chunk$set(message=FALSE, warning=FALSE)
#' #Preamble####
#+ eval=TRUE
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

#' #Section 1. Recognizer evaluation####
#' 
#'  We followed Knight et al. 2017 to evaluate the precision and recall of our recognizer and select a score threshold for processing. 
#' 
#' ##1a. Prepare data####
#' 
#' First, we randomly selected 20 night time recordings with EWPW present and 20 night time recordings with EWPW absent to evaluate our recognizer. We visually and/or aurally processed those recordings and counted the number of EWPW calls in each calling bout in each recording. We defined the start and end time of bouts by gaps of at least 2 seconds between other calls.  We subjectively evaluated the amplitude of the calls in each bout on a scale of 1-5. We also noted the maximum abundance of EWPW calling simultaneously in each recording.
#' 
#+ eval=TRUE
#Read in and summarize benchmark data (human listener processed)
path <- "/Users/ellyknight/Documents/UoA/Data/AutomatedProcessing/EWPW/"
bench <- read_excel(paste0(path, "TestDatasetCallCount.xlsx"),
                    col_types = c("text", "text", "text", "numeric", "text", "text")) %>% 
  mutate(File=paste0(File, ".wav")) %>% 
  rename(file=File)

bench.rec <- bench %>% 
  dplyr::select(file) %>% 
  unique()

call <- sum(bench$Count)

#' We then used our trained SongScope recognizer to process those 40 selected recordings with a score threshold of 20 and a quality threshold of 20. We visually validated the recognizer results to separate true and false positives.
#' 
#+ eval=TRUE
#Read in raw validated data
dat.test <- read.table("/Users/ellyknight/Documents/UoA/Data/AutomatedProcessing/EWPW/Test_EWPWV1_0_20_results_validated.txt",
                       sep="\t",
                       header=FALSE,
                       col.names = c("filepath", "start", "duration", "rsl", "quality", "score", "recognizer", "validation")) %>% 
  separate(filepath, into=c("f1", "f2", "f3", "f4", "presence", "file"), sep="/", remove=FALSE) %>% 
  separate(file, into=c("aru", "datetext", "timetext", "filetype"), remove=FALSE) %>% 
  mutate(ewpw=ifelse(validation=="y", 1, 0))

#' ##1b. Evaluate for detection####
#' 
#' Then we compared number of detections in the recognizer results to the number of detections in the benchmark data across all possible score thresholds to determine precision and recall.
#' 
#+ eval=TRUE
score <- seq(round(min(dat.test$score)), round(max(dat.test$score)))
results.call <- data.frame()
for(i in 1:length(score)){
  
  score.i <- as.numeric(score[i])
  
  dat.i <- dat.test %>% 
    dplyr::filter(score >= score.i)
  
  results.i <- dat.i %>% 
    dplyr::summarize(det=sum(ewpw), hit=n()) %>% 
    mutate(r=det/call, p=det/hit, fp=(hit-det)/hit, tp=det/hit, fn=(call-det)/call, acc=(det-(hit-det))/hit) %>% 
    mutate(f=(2*p*r)/(p+r),
           score=score[i]) %>% 
    ungroup()
  
  results.call <- rbind(results.call, results.i)
  
}

#' ##1c. Evaluate for presence/absence per recording####
#' 
#' We also determined the recall of EWPW presence per recording across all possible score thresholds.
#' 
#+ eval=TRUE
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
  summarize(r=sum(pres.rec)/nrow(bench.rec)) %>% 
  ungroup()

#' ##1d. Plot evaluation results####
#' 
#' Based on our results, we selected a score threshold of 60 to minnimize false positives, while still detecting EWPW in 90% of recordings where they were detected by human listener.
#' 
#+ eval=TRUE
#Plot recall
r <- ggplot(results.call) +
  geom_smooth(aes(x=score, y=r), colour="blue") +
  geom_vline(aes(xintercept=60)) +
  xlab("") +
  ylab("Recall") +
  ylim(0, 1)

#+ eval=TRUE
#Plot precision
p <- ggplot(results.call) +
  geom_smooth(aes(x=score, y=p), colour="blue") +
  geom_vline(aes(xintercept=60)) + 
  xlab("") +
  ylab("Precision") +
  ylim(0, 1)

#+ eval=TRUE
#Plot presence-absence recall
pres <- ggplot(r.rec) +
  geom_smooth(aes(x=score, y=r)) +
  geom_vline(aes(xintercept=60)) +
  xlab("Score threshold") +
  ylab("Presence-absence recall") +
  ylim(0, 1)

#+ eval=TRUE
#+ fig.width=8, fig.height=3
#Put all three plots together
plot.eval <- gridExtra::grid.arrange(p, r, pres, ncol=3)

#' # Section 2. Summary of automated processing
#' 
#' ## 2a. Precision of recognizer####
#+ eval=TRUE
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
                              validation %in% c("0", "na", "Na", "?") ~ 0))
#'
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

#' ##2b. Recording availability
#' 
#Read in clean dataset summarized by 1 minute of recording
load("/Users/ellyknight/Dropbox/SS/EWPW/EWPW occupancy analysis data.Rdata") 

#Length and number of recordings per location
recordings <- df_analysis_occupancy %>% 
  mutate(ID=paste0(Year, "-", SiteName)) %>% 
  group_by(Year, ID, SiteName, Recording) %>% 
  summarize(Duration = max(Rec_minute)) %>% 
  ungroup()
table(recordings$ID, recordings$Duration)

#Total number of recording minutes per location
recordings.sum <- recordings %>% 
  group_by(Year, ID, SiteName) %>% 
  summarize(Duration.hr=sum(Duration)/60) %>% 
  ungroup() %>% 
  arrange(Duration.hr)

#+ fig.width=4, fig.height=4
ggplot(recordings.sum) +
  geom_histogram(aes(x=Duration.hr, fill=factor(Year))) +
  scale_fill_viridis_d("Year")

#' ##2c. Summary of detections
#' 
#Summarize by site
detections <- df_analysis_occupancy %>% 
  group_by(Year, SiteName) %>% 
  summarize(DetectionMinutes = sum(Occupied),
            Num_detects = sum(Num_detects)) %>% 
  ungroup() %>% 
  left_join(recordings.sum) %>% 
  mutate(Occupied = ifelse(DetectionMinutes > 0, 1, 0),
         DetectsPerMinute = Num_detects/Duration.hr)

#Number of sites with detections
table(detections$Occupied)

#Detections per minute of recording at occupied sites
detections %>% 
  dplyr::filter(Occupied==1) %>% 
  dplyr::select(DetectsPerMinute) %>% 
  summary()

#' ##2d. Summary of abundance
#'
#Summarize by recording
abundance <- df_analysis_occupancy %>% 
  dplyr::select(Year, SiteName, Recording, Abundance_site, Abundance_recording) %>% 
  unique()

#Abundance per occupied recording
abundance %>% 
  dplyr::filter(Abundance_recording > 0) %>% 
  dplyr::select(Abundance_recording) %>% 
  summary()

#Distribution of abundance across occupied sites
abundance.site <- abundance %>% 
  dplyr::select(Year, SiteName, Abundance_site) %>% 
  unique()

#+ fig.width=4, fig.height=4
ggplot(abundance.site %>% dplyr::filter(Abundance_site > 0)) +
  geom_histogram(aes(x=Abundance_site))

#' #Section 3: Preliminary detectability analysis
#' 
#' First, we used occupancy modelling to determine which temporal, solar, lunar, and weather covariates are important predictors of EWPW detectability so that we could constrain further analyses to recordings with optimal conditions. 
#' 
#' Because we were interested in predicting EWPW detectability, we used only those sites for which EWPW presence was known. In other words, we only used data from the 34 sites where EWPW were detected in our dataset, and removed the 7 sites were no EWPW were detected.
#' 
#' ##3a. Wrangle data
#' 
#' To standardize sampling effort per visit (i.e., recording), we used only the first five minutes of every recording. One site from 2019 had 244 three-minute recordings, which we removed from the dataset.
#+ eval=TRUE
#Identify sites with no detections
df <- df_analysis_occupancy %>% 
  group_by(SiteName, Year) %>% 
  summarize(detections=sum(Occupied)) %>% 
  dplyr::filter(detections > 0) %>% 
  ungroup() %>% 
  left_join(df_analysis_occupancy)

#Remove three minute recordings
duration <- df %>% 
  group_by(SiteName, Recording) %>% 
  summarize(Duration = max(Rec_minute)) %>% 
  dplyr::filter(Duration!=3) %>% 
  ungroup()
table(duration$SiteName, duration$Duration)

#Summarize by recording
dat <- df %>% 
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

#' ##3b. Check for VIF & correlation
#' We checked all potential covariates for variance inflation and correlation and removed any variables with VIF > 5 or correlation > 0.7.
#' 
#Check VIF & correlation
covs <- dat %>% 
  dplyr::select(Yday, sunsetTimeSince, moon_fraction, moon_altitude, temp, total_rain, wind_spd, band1_psd, band2_psd, band1_s2n, band2_s2n) %>% 
  data.frame()

vif(covs)
cor(covs)
#Remove band1_psd

covs <- dat %>% 
  dplyr::select(Yday, sunsetTimeSince, moon_fraction, moon_altitude, temp, total_rain, wind_spd, band2_psd, band1_s2n, band2_s2n) %>% 
  data.frame()

vif(covs)
cor(covs)
#all good

#' ##3c. Covariate selection
#' 
#'   We then randomly sampled 200 5-minute visits from each of the 34 study sites to even out sampling across the study sites. We fit two sets of occupancy models: one with all potential solar, lunar, and temporal covariates, and one with all weather covariates. For each set, we fit one global model and all potential combinations of the covariates in that model. We bootstrapped this visit selection and model fitting process 100 times. 
#'   
#+ eval=TRUE
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
    sample_n(200, replace=TRUE) %>%
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
  temp.dredge <-  dredge(mod.temp)
  temp.dredge.list[[i]] <- temp.dredge %>% 
    data.frame() %>% 
    mutate(boot=i,
           mod=row.names(temp.dredge))
  
  #Weather next
  mod.weather <- occu(~ temp + I(temp^2) + rain + wind + psd2 + s2n2
                      ~1,
                      data=dat.occ)
  weather.dredge <-  dredge(mod.weather)
  weather.dredge.list[[i]] <- weather.dredge %>% 
    data.frame() %>% 
    mutate(boot=i,
           mod=row.names(weather.dredge))
  
}

#'We selected the best fitting model from each set as the model with the highest mean model weight across the 100 bootstraps.
#' 
#+ eval=TRUE  
#Summarize bootstrapped models
temp.dredge.all <- rbindlist(temp.dredge.list) %>% 
  group_by(mod) %>% 
  summarize_all(~mean(.x)) %>% 
  ungroup() %>% 
  arrange(-weight)  %>% 
  head(1) %>% 
  pivot_longer(p.doy.:p.moonalt.moonfrac., names_to="variable", values_to="coefficient") %>% 
  dplyr::select(variable, coefficient, delta, weight) %>% 
  dplyr::filter(!is.na(coefficient))
temp.dredge.all
#Select Doy, suntime, suntime^2, moonalt

weather.dredge.all <- rbindlist(weather.dredge.list) %>% 
  group_by(mod) %>% 
  summarize_all(~mean(.x)) %>% 
  ungroup() %>% 
  arrange(-weight) %>% 
  head(1) %>% 
  pivot_longer(p.psd2.:p.wind., names_to="variable", values_to="coefficient") %>% 
  dplyr::select(variable, coefficient, delta, weight) %>% 
  dplyr::filter(!is.na(coefficient))
weather.dredge.all 
#Select psd2, rain, ps2n2, temp, temp^2, wind (everything)

#' ##3d. Fit final model
#' 
#'   We then combined all the covariates from the two best fitting models into a final model. We again randomly sampled 200 5-minute visits and fit them to this final model. We bootstrapped visit sampling and model fitting 100 times.
#' 

#Create new data for model prediction in loop
moondat <- expand.grid(moonalt=seq(round(min(dat$moon_altitude), -1), max(dat$moon_altitude), 0.1),
                      suntime=1,
                      doy=mean(dat$Yday),
                      wind=min(dat$wind_spd),
                      rain=min(dat$total_rain),
                      temp=mean(dat$temp),
                      psd2=max(dat$band2_psd),
                      s2n2=min(dat$band2_s2n))

suntimedat <- expand.grid(moonalt=mean(dat$moon_altitude),
                          suntime=seq(round(min(dat$sunsetTimeSince), -1), max(dat$sunsetTimeSince), 1),
                          doy=mean(dat$Yday),
                          wind=min(dat$wind_spd),
                          rain=min(dat$total_rain),
                          temp=mean(dat$temp),
                          psd2=max(dat$band2_psd),
                          s2n2=min(dat$band2_s2n))

winddat <- expand.grid(moonalt=mean(dat$moon_altitude),
                       suntime=1,
                       doy=mean(dat$Yday),
                       wind=seq(round(min(dat$wind_spd)), max(dat$wind_spd), 1),
                       rain=min(dat$total_rain),
                       temp=mean(dat$temp),
                       psd2=max(dat$band2_psd),
                       s2n2=min(dat$band2_s2n))

tempdat <- expand.grid(moonalt=mean(dat$moon_altitude),
                       suntime=1,
                       doy=mean(dat$Yday),
                       wind=min(dat$wind_spd),
                       rain=min(dat$total_rain),
                       temp=seq(min(round(dat$temp)), max(dat$temp), 1),
                       psd2=max(dat$band2_psd),
                       s2n2=min(dat$band2_s2n))

#+ eval=TRUE
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
    sample_n(200, replace=TRUE) %>%
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
  
  obs.cov <- list(doy=doy, suntime=suntime, moonalt=moonalt, temp=temp, rain=rain, wind=wind, psd2=psd2, s2n2=s2n2, abundance.site=abundance.site, abundance.rec=abundance.rec)
  
  site.cov <- dat.sub %>% 
    dplyr::select(SiteName, Latitude, Longitude) %>% 
    unique() %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  dat.occ <- unmarkedFrameOccu(detections, siteCovs=site.cov, obsCovs=obs.cov)
  
#Fit model 
  mod.final <- occu(~ moonalt + suntime + I(suntime^2) + doy + temp + I(temp^2) + wind + rain + psd2 + s2n2
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
    mutate(boot=i)
  
}

#'We summarized the 100 bootstraps to obtain mean and standard error estimates for the coefficients in our final model.
#'
#+ eval=TRUE
#Summarize bootstrapped coefficients
mod.final.coeff <- rbindlist(mod.final.coeff.list) %>% 
  group_by(var) %>% 
  summarize(est.mn=mean(Estimate),
            est.sd=sd(Estimate),
            se.mn=mean(SE),
            se.sd=sd(SE)) %>% 
  ungroup()
mod.final.coeff

#Summarize model predictions
mod.final.pred <- rbindlist(mod.final.pred.list) %>% 
  group_by(plot, moonalt, suntime, wind, temp) %>% 
  summarize(det.mn = mean(Det),
            se.mn = mean(DetSE),
            up.mn = mean(DetUpper), 
            lw.mn = mean(DetLower)) %>% 
  ungroup()

#' ##3e. Plot final model
#' 
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

#+ eval=TRUE
#+ fig.width=8, fig.height=8
#Put all four plots together
plot.pred <- gridExtra::grid.arrange(plot.moon, plot.suntime, plot.wind, plot.temp,
                                     ncol=2, nrow=2)
