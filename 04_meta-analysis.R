#title: "Refining surveys for nocturnal birds using acoustic data: Eastern Whip-poor-will - supplementary results for meta-analysis by Matt Broadway"
#author: "Elly C. Knight"
#date: "May 13, 2025

#Preamble#####

library(tidyverse)

options(scipen = 9999, max.print=2000)

#Load data
load("Analysisv7_workspace.Rdata")

#Create lists to save out results
mod.moonalt.list <- list()
mod.moonfrac.list <- list()

#Set number of bootstraps  
boot <- 100

#Start for loop
for(i in 1:boot){
  
  #Sample data
  set.seed(i)
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
  
  moonfrac <- dat.sub %>% 
    dplyr::select(SiteName, n, moon_fraction) %>% 
    spread(key = n, value = moon_fraction) %>% 
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
  
  obs.cov <- list(doy=doy, suntime=suntime, moonalt=moonalt, moonfrac=moonfrac, temp=temp, rain=rain, wind=wind, psd2=psd2, s2n2=s2n2, abundance.site=abundance.site, abundance.rec=abundance.rec, year=year)
  
  site.cov <- dat.sub %>% 
    dplyr::select(SiteName, Latitude, Longitude) %>% 
    unique() %>% 
    arrange(SiteName) %>% 
    dplyr::select(-SiteName) %>% 
    data.frame()
  
  dat.occ <- unmarkedFrameOccu(detections, siteCovs=site.cov, obsCovs=obs.cov)
  
  #Fit model 
  #model structure selected via results described in Knight et al. 2022 (https://ace-eco.org/vol17/iss1/art21/)
  mod.moonalt.list[[i]]  <- occu(~ moonalt + suntime + I(suntime^2) + doy + temp + I(temp^2) + wind + psd2 + s2n2
                    ~1,
                    data=dat.occ)
  
  mod.moonfrac.list[[i]] <- occu(~ moonfrac + suntime + I(suntime^2) + temp + I(temp^2) + wind + psd2 + s2n2
                      ~1,
                      data=dat.occ)
  
  
  
  print(paste0("Finished bootstrap ", i, " of ", boot))
  
}

#Save
save(mod.moonalt.list, mod.moonfrac.list, file="MoonFrac&MoonAltModels.Rdata")
