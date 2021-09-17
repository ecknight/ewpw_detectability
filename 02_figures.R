#title: "Figures for refining surveys for nocturnal birds using acoustic data: Eastern Whip-poor-will"
#author: "Elly C. Knight"
#date: "February 15, 2021"

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
library(nord)
library(cowplot)
library(patchwork)
library(grid)
library(ggspatial)
library(ggsn)
library(scales)

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

load("Analysisv7_workspace.Rdata")


#Figure 1. Study area----

#Remove sites with no detections
load("/Users/ellyknight/Dropbox/SS/EWPW/EWPW occupancy analysis data.Rdata")

sites <- df_analysis_occupancy %>% 
  dplyr::filter(as.numeric(sunsetTimeSince) > 0,
                as.numeric(sunriseTimeSince) < 0, 
                Rec_minute <= 5) %>% 
  group_by(SiteName, Year) %>% 
  summarize(detections=sum(Occupied)) %>%  
  dplyr::filter(detections > 0) %>% 
  ungroup() %>% 
  left_join(df_analysis_occupancy) %>% 
  dplyr::select(SiteName, Year, SongMeter, Latitude, Longitude) %>% 
  unique()

#1a. Map of North America----
nam <- map_data("world", region=c("Canada", 
                                  "USA", 
                                  "Mexico",
                                  "Guatemala", 
                                  "Belize", 
                                  "El Salvador",
                                  "Honduras", 
                                  "Nicaragua", 
                                  "Costa Rica",
                                  "Panama", 
                                  "Jamaica", 
                                  "Cuba", 
                                  "The Bahamas",
                                  "Haiti", 
                                  "Dominican Republic", 
                                  "Antigua and Barbuda",
                                  "Dominica", 
                                  "Saint Lucia", 
                                  "Saint Vincent and the Grenadines", 
                                  "Barbados",
                                  "Grenada",
                                  "Trinidad and Tobago")) %>% 
  dplyr::filter(!group%in%c(258:264))

nam.eq <- nam %>% 
  st_as_sf(coords=c("long", "lat"), crs=4326) %>% 
  st_transform(crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs") %>%
  st_coordinates() %>% 
  data.frame() %>% 
  cbind(nam)

area.eq.center <- sites %>% 
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326) %>% 
  st_transform(crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs") %>%
  st_coordinates() %>% 
  as.data.frame() %>% 
  dplyr::summarize(X=mean(X),
                   Y=mean(Y))
#Read in EWPW range
range <- read_sf("/Volumes/ECK004/GIS/Misc/BirdLifeRanges/EWPW.shp") %>% 
  st_transform(crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs") %>% 
  dplyr::filter(SEASONAL==2)

map.nam <- ggplot() +
  geom_polygon(data=nam.eq, aes(x=X, y=Y, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_sf(data=range, aes(fill=SCINAME), fill="grey50", colour="gray85") +
  geom_text(label="★", aes(x=X, y=Y), size=6, family = "HiraKakuPro-W3", data=area.eq.center, colour="black") +
  xlim(c(-4000000, 3000000)) +
  coord_sf(datum = NA) +
  xlab("") +
  ylab("") +
  my.theme +
  theme(plot.margin = unit(c(0,0,-1,-1), "cm"),
        legend.position="bottom")
map.nam

#1b. Study sites----
sites.sf <- sites %>% 
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)

#area.center <- sites.sf %>% 
#  st_coordinates() %>% 
#  as.data.frame() %>% 
#  dplyr::summarize(X=mean(X),
#                   Y=mean(Y))
area.center <- data.frame(X=-77, Y=44.6)

#Get background data
register_google(key="AIzaSyCta9P4x7jGNELznpwlx07VZkkLVk3FP4M")

map <- get_map(maptype="satellite", location=area.center, zoom=8, force=TRUE, color="color")

map_attributes <- attributes(map)

map_transparent <- matrix(adjustcolor(map, 
                                      alpha.f = 0.8), 
                          nrow = nrow(map))
attributes(map_transparent) <- map_attributes

#Weather stations
weather <- data.frame(name=c("“OTTAWA CDA RCS”, “TRENTON A"),
                      Latitude = c(45.3825, 44.1189),
                      Longitude = c(-75.7141, -77.5281),
                      id=rep("Weather station",2))

#Plot site map
map.sites <- ggmap(map_transparent) +
  geom_spatial_point(aes(x = Longitude, y = Latitude,
                         fill=factor(Year)),
                     data = sites, 
                     crs=4326,
                     alpha = 1,
                     colour="grey85",
                     size=4,
                     shape=21) +
  geom_spatial_point(aes(x = Longitude, y = Latitude, colour=id),
                     data = weather, 
                     crs=4326,
                     size=4,
                     shape=23,
                     fill="black") +
  geom_text(label="Lake Ontario", aes(x=-77.5, y=43.6), size=5, colour="black") +
  ggspatial::annotation_north_arrow(location = "bl",
                                    style = ggspatial::north_arrow_orienteering(fill = c("grey80", "grey20"), line_col = "grey20")) +
  ggsn::scalebar(x.min = -76.5, x.max = -75.6, 
                 y.min = 43.55, y.max = 43.8, 
                 transform=TRUE, model="WGS84",
                 dist=25, dist_unit="km",
                 box.fill=c("grey80", "grey20"),
                 box.color="grey20",
                 height=0.1,
                 st.bottom=TRUE, st.dist=0.15) +
  xlim(c(-78.1, -75.5)) +
  ylim(c(43.5, 45.5)) +
  scale_fill_manual(values=c("darkgoldenrod1", "tomato3"), name="Year surveyed") +
  scale_colour_manual(values=c("grey85"), name="") +
  my.theme +
  xlab("") +
  ylab("") +
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        legend.position = "bottom")
map.sites

#1c. Put together----
plot.sa <- map.sites +
  inset_element(map.nam,
                right=0.48,
                bottom=0.6,
                left=0.02,
                top=0.98)
#plot.sa

ggsave(file="Figures/Figure1StudyArea.jpeg", height=8, width=6, units="in", device="jpeg",
       plot=plot.sa)

#Figure 3. Evaluation results####

results.eval <- read.csv("RecognizerEvaluationResults.csv") %>% 
  dplyr::select(score, p, r, rpres) %>% 
  pivot_longer(cols=p:rpres, names_to="metric", values_to="value")

plot.eval <- ggplot(results.eval) +
  geom_line(aes(x=score, y=value, colour=metric), size=1) +
  geom_vline(aes(xintercept=60), linetype="dashed") +
  scale_colour_manual(values=c("darkgoldenrod3", "#A8A8A8", "tomato4"), name="Evaluation metric",
                      labels=c("Precision", "Recall", "Presence-absence recall\n(per 5-minute recording)")) +
  xlab("Score threshold") +
  ylab("Evaluation metric value") +
  ylim(0, 1) +
  my.theme
plot.eval

ggsave(plot.eval, file="Figures/Fig3Evaluation.jpg", width=7, height=4, units="in", device="jpeg")

#1e. Summary statistics----

#Number of EWPW calls in test dataset
sum(bench$Count)

#Number of EWPW calls detected by recognizer
sum(dat.test$ewpw)

#Number of EWPW calls detected by recognizer with score threshold of 60
dat.test %>% 
  dplyr::filter(score >= 60) %>% 
  summarize(ewpw=sum(ewpw))

#Proportion of calls detected at score threshold of 60
dat.test %>% 
  dplyr::filter(score >= 60) %>% 
  summarize(ewpw=sum(ewpw)/sum(bench$Count))

#Number of recordings that EWPW is detected in
dat.test %>% 
  group_by(file) %>% 
  summarize(ewpw=sum(ewpw),
            ewpw=ifelse(ewpw>0, 1, 0)) %>% 
  ungroup() %>% 
  summarize(rec=sum(ewpw))

#Number of recordings that EWPW is detected in at score threshold of 60
dat.test %>% 
  dplyr::filter(score >= 60) %>% 
  group_by(file) %>% 
  summarize(ewpw=sum(ewpw),
            ewpw=ifelse(ewpw>0, 1, 0)) %>% 
  ungroup() %>% 
  summarize(rec=sum(ewpw))

20/17 #Offset

#Table 1. Model selection results----
temp.dredge.all <- read.csv("TemporalCovariatesModelSelection.csv")

temp.dredge.raw <- read.csv("TemporalCovariatesModelSelection_raw.csv") %>% 
  arrange(mod, boot) %>% 
  dplyr::filter(mod %in% c(27:32))

temp.dredge.top <- temp.dredge.raw %>% 
  group_by(boot) %>% 
  top_n(1, wt=weight) %>% 
  ungroup() %>% 
  group_by(mod) %>% 
  summarize(percent = n()/100) %>% 
  arrange(-percent)

temp.dredge.summary <- temp.dredge.raw %>% 
  dplyr::select(mod, delta, weight) %>% 
  group_by(mod) %>% 
  summarize(wt.mn = mean(weight),
            wt.sd = sd(weight),
            d.mn = mean(delta),
            d.sd = sd(delta)) %>% 
  ungroup() %>% 
  left_join(temp.dredge.top) %>% 
  arrange(-wt.mn)
temp.dredge.summary

weather.dredge.all <- read.csv("WeatherCovariatesModelSelection.csv")

weather.dredge.raw <- read.csv("WeatherCovariatesModelSelection_raw.csv") %>% 
  arrange(mod, boot) %>% 
  dplyr::filter(mod %in% c(61, 62, 64, 34, 58, 57))

weather.dredge.top <- weather.dredge.raw %>% 
  group_by(boot) %>% 
  top_n(1, wt=weight) %>% 
  ungroup() %>% 
  group_by(mod) %>% 
  summarize(percent = n()/100) %>% 
  arrange(-percent)

weather.dredge.summary <- weather.dredge.raw %>% 
  dplyr::select(mod, delta, weight) %>% 
  group_by(mod) %>% 
  summarize(wt.mn = mean(weight),
            wt.sd = sd(weight),
            d.mn = mean(delta),
            d.sd = sd(delta)) %>% 
  ungroup() %>% 
  left_join(weather.dredge.top) %>% 
  arrange(-wt.mn)
View(weather.dredge.summary)


#Figure 4. Covariate effects----

mod.final.pred <- read.csv("FinalModelPredictions.csv")
mod.final.pred.sum <- mod.final.pred %>% 
  group_by(plot, moonalt, suntime, doy, wind, temp, psd2, s2n2) %>% 
  summarize(det.mn = mean(Det),
            se.mn = mean(DetSE),
            up.mn = mean(DetUpper), 
            lw.mn = mean(DetLower)) %>% 
  ungroup()

plot.moon <- ggplot(subset(mod.final.pred.sum, plot=="Moon altitude")) +
  geom_line(aes(x=moonalt, y=det.mn), colour=nord_palettes$baie_mouton[1]) +
  geom_ribbon(aes(x=moonalt, ymin=lw.mn, ymax=up.mn), alpha=0.5, fill=nord_palettes$baie_mouton[1]) +
  ylim(c(0,1)) +
  xlab("Moon altitude") +
  my.theme +
  theme(axis.title.y=element_blank())

plot.suntime <- ggplot(subset(mod.final.pred.sum, plot=="Time since sunset (hrs)")) +
  geom_line(aes(x=suntime, y=det.mn), colour=nord_palettes$baie_mouton[2]) +
  geom_ribbon(aes(x=suntime, ymin=lw.mn, ymax=up.mn), alpha=0.5, fill=nord_palettes$baie_mouton[2]) +
  ylim(c(0,1)) +
  xlab("Time since sunset (hrs)") +
  my.theme +
  theme(axis.title.y=element_blank())

plot.day <- ggplot(subset(mod.final.pred.sum, plot=="Day of year")) +
  geom_line(aes(x=doy, y=det.mn), colour=nord_palettes$baie_mouton[3]) +
  geom_ribbon(aes(x=doy, ymin=lw.mn, ymax=up.mn), alpha=0.5, fill=nord_palettes$baie_mouton[3]) +
  ylim(c(0,1)) +
  xlab("Day of year") +
  my.theme +
  theme(axis.title.y=element_blank())

plot.wind <- ggplot(subset(mod.final.pred.sum, plot=="Wind speed (km/h)")) +
  geom_line(aes(x=wind, y=det.mn), colour=nord_palettes$baie_mouton[4]) +
  geom_ribbon(aes(x=wind, ymin=lw.mn, ymax=up.mn), alpha=0.5, fill=nord_palettes$baie_mouton[4]) +
  ylim(c(0,1)) +
  xlab("Wind speed (km/h)") +
  my.theme +
  theme(axis.title.y=element_blank())

plot.temp <- ggplot(subset(mod.final.pred.sum, plot=="Temperature (C)")) +
  geom_line(aes(x=temp, y=det.mn), colour=nord_palettes$baie_mouton[5]) +
  geom_ribbon(aes(x=temp, ymin=lw.mn, ymax=up.mn), alpha=0.5, fill=nord_palettes$baie_mouton[5]) +
  ylim(c(0,1)) +
  xlab("Temperature (C)") +
  my.theme +
  theme(axis.title.y=element_blank())

plot.psd <- ggplot(subset(mod.final.pred.sum, plot=="Power spectrum density")) +
  geom_line(aes(x=psd2, y=det.mn), colour=nord_palettes$baie_mouton[6]) +
  geom_ribbon(aes(x=psd2, ymin=lw.mn, ymax=up.mn), alpha=0.5, fill=nord_palettes$baie_mouton[6]) +
  ylim(c(0,1)) +
  xlab("Power spectrum density") +
  my.theme +
  theme(axis.title.y=element_blank())

plot.s2n <- ggplot(subset(mod.final.pred.sum, plot=="Signal to noise ratio")) +
  geom_line(aes(x=s2n2, y=det.mn), colour=nord_palettes$baie_mouton[7]) +
  geom_ribbon(aes(x=s2n2, ymin=lw.mn, ymax=up.mn), alpha=0.5, fill=nord_palettes$baie_mouton[7]) +
  ylim(c(0,1)) +
  xlab("Signal to noise ratio") +
  my.theme +
  theme(axis.title.y=element_blank())

#Put plots together
plot.pred <- gridExtra::grid.arrange(plot.moon, plot.suntime, plot.day, plot.wind, plot.temp, plot.psd, plot.s2n,
                                     ncol=4, nrow=2,
                                     left = textGrob("EWPW detectability", rot = 90, vjust = 1))

ggsave(plot.pred, file="Figures/Fig4CovariatePredictions.jpg", width=12, height=6, units="in", device="jpeg")

#Figure 5. Sampling effects----

mod.pred <- read.csv("SamplingEffortResults.csv") %>% 
  dplyr::filter(model=="cov")

summary <- read.csv("SamplingEffortResults_summary.csv")

summary.sum <- summary %>% 
  group_by(cov, length, visit, boot) %>% 
  summarize(sites = sum(Occupied)/32) %>% 
  ungroup()

vline <- data.frame(xint=c(10, NA), cov=c("all", "constrained"))

mod.pred$cov <- factor(mod.pred$cov, levels=c("all", "constrained"), labels=c("Unconstrained", "Constrained"))
summary.sum$cov <- factor(summary.sum$cov, levels=c("all", "constrained"), labels=c("Unconstrained", "Constrained"))
vline$cov <- factor(vline$cov, levels=c("all", "constrained"), labels=c("Unconstrained", "Constrained"))

plot.occu <- ggplot(mod.pred) +
  geom_ribbon(aes(x=visit, ymin=occu.lwr.mn, ymax=occu.upr.mn, fill=factor(length)), alpha=0.3) +
  geom_line(aes(x=visit, y=occu.mn.mn, colour=factor(length)), size=1) +
  geom_vline(aes(xintercept=xint), linetype="dashed", data=vline) +
  scale_colour_nord("aurora", name="Recording\nlength\n(minutes)") + 
  scale_fill_nord("aurora", name="Recording\nlength\n(minutes)") +
  xlab("Number of recordings") +
  ylab("Probability of occupancy") +
#  xlim(c(1,10)) +
  ylim(c(0,1)) +
  my.theme +
  facet_wrap(~cov, scales="free_x") +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_x_continuous(breaks= pretty_breaks())
plot.occu

plot.det <- ggplot(mod.pred) +
  geom_ribbon(aes(x=visit, ymin=det.lwr.mn, ymax=det.upr.mn, fill=factor(length)), alpha=0.3) +
  geom_line(aes(x=visit, y=det.mn.mn, colour=factor(length)), size=1) +
  geom_vline(aes(xintercept=xint), linetype="dashed", data=vline) +
  scale_colour_nord("aurora", name="Recording\nlength\n(minutes)") + 
  scale_fill_nord("aurora", name="Recording\nlength\n(minutes)") +
  xlab("") +
  ylab("Probability of detection") +
#  xlim(c(1,10)) +
  ylim(c(0,1)) +
  my.theme +
  facet_wrap(~cov, scales="free_x") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks= pretty_breaks())
plot.det

plot.legend <- ggplot(summary.sum) +
  geom_point(aes(x=visit, y=sites, colour=factor(length), fill=factor(length))) +
  geom_smooth(aes(x=visit, y=sites, colour=factor(length), fill=factor(length))) +
  scale_colour_nord("aurora", name="Recording length (minutes)") + 
  scale_fill_nord("aurora", name="Recording length (minutes)") +
  xlab("Number of recordings") +
  ylab("Cumulative probability of detection") +
  #  xlim(c(1,10)) +
  ylim(c(0,1)) +
  my.theme +
  facet_wrap(~cov, scales="free_x") +
  theme(legend.position = "bottom")

legend <- get_legend(plot.legend)

ggsave(file="Figures/Fig5OccupancySamplingEffort.jpeg", height=8, width=8, units="in", device="jpeg",
       grid.arrange(plot.det, plot.occu, legend,
                    widths = c(8),
                    heights = c(4, 4, 1),
                    layout_matrix = rbind(c(1),
                                          c(2),
                                          c(3))))

#Figure 6. Cumulative detection probability comparison----

pred.glmm <- read.csv("GLMMPredictions.csv")

plot.glmm <- ggplot(pred.glmm) +
  geom_ribbon(aes(x=visit, ymin=lwr, ymax=upr, fill=factor(length)), alpha=0.2, data=subset(pred.glmm, cov=="all")) +
  geom_line(aes(x=visit, y=pred, colour=factor(length), linetype=cov)) +
  scale_colour_nord("aurora", name="Recording length\n(minutes)") + 
  scale_fill_nord("aurora", name="Recording length\n(minutes)") +
  scale_linetype_discrete(name="Recording selection", labels=c("Unconstrained", "Constrained")) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10)) +
  xlab("Number of recordings") +
  ylab("Cumulative probability of detection") +
  my.theme

ggsave(plot.glmm, file="Figures/Fig6CumulativeProbabilityComparison.jpeg", height=4, width=6, units="in", device="jpeg")

#Figure 7. Logistic Growth curve----
summary.sum.any <- read.csv("NLSData.csv")

pred.any.visit <- read.csv("NLSPredictions_Visits.csv")
pred.asym.visit <- read.csv("NLSAsymptotes_Visits.csv") %>% 
  rbind(data.frame(asym=rep(0,5), n=read.csv("NLSAsymptotes_Visits.csv")$n, length=c(1:5), minutes=n*length)) %>% 
  mutate(linetype="Sample size for asymptote",
         asym.99 = asym*0.99) %>% 
  dplyr::filter(n!="Inf")

plot.nls.visit <- ggplot() +
  geom_jitter(data=summary.sum.any, aes(x=visit, y=sites, colour=factor(length))) +
  geom_line(data=pred.any.visit, aes(x=visit, y=r, colour=factor(length))) +
  geom_line(data=pred.asym.visit, aes(x=n, y=asym.99, colour=factor(length), linetype=linetype)) +
  scale_colour_nord("aurora", name="Recording length\n(minutes)") + 
  scale_linetype_manual(name="", labels=c("Sample size\nfor asymptote"), values=c("dashed")) +
  ylim(c(0,1)) +
  xlab("Number of recordings") +
  ylab("Cumulative proportion of sites with detections") +
  my.theme
plot.nls.visit

pred.any.length <- read.csv("NLSPredictions_Length.csv")
summary.sum.any.length <- summary.sum.any %>% 
  dplyr::filter(visit %in% unique(pred.any.length$visit))

plot.nls.length <- ggplot() +
  geom_jitter(data=summary.sum.any.length, aes(x=length, y=sites, colour=factor(visit))) +
  geom_line(data=pred.any.length, aes(x=length, y=r, colour=factor(visit))) +
  scale_colour_nord("aurora", name="Number of\nrecordings") + 
  scale_linetype_manual(name="", labels=c("Sample size\nfor asymptote"), values=c("dashed")) +
  ylim(c(0,1)) +
  xlab("Recording length (minutes)") +
  ylab("") +
  my.theme
#plot.nls.length

ggsave(grid.arrange(plot.nls.visit, plot.nls.length, ncol=2, widths=c(7,6)), file="Figures/Fig7CumulativeProbabilityNLS.jpeg", height=4, width=10, units="in", device="jpeg")

View(pred.asym.visit)


#Summary stats
mod.pred %>% 
  group_by(cov) %>% 
  summarize(mean=mean(det.mn.mn))

mod.pred %>% 
  group_by(cov) %>% 
  summarize(mean=mean(occu.mn.mn))


#Appendix 3. RSL----

r.eval <- read.csv("RecognizerEvaluationResults.csv")

plot.rsl1 <- ggplot(r.eval) +
  geom_ribbon(aes(x=r, ymax=rsltp + rslsdtp, ymin=rsltp - rslsdtp), fill="grey85", colour="grey85") +
  geom_line(aes(x=r, y=rsltp), colour="grey30") +
  geom_point(aes(x=r, y=rsltp, fill=score), shape=21, colour="grey30", size=3) +
  scale_fill_viridis_c(name="Score threshold", breaks=c(20, 50, 80)) +
  xlab("Recall") +
  ylab("Relative sound level (RSL)") +
  my.theme +
  theme(legend.position = "bottom",
        axis.title.y = element_text(vjust = 10),
        plot.margin = margin(l=30))
plot.rsl1

legend <- get_legend(plot.rsl1)

plot.rsl <- plot.rsl1 + theme(legend.position="none")

plot.level <- ggplot(r.eval) +
  geom_ribbon(aes(x=r, ymax=leveltp + levelsdtp, ymin=leveltp - levelsdtp), fill="grey85", colour="grey85") +
  geom_line(aes(x=r, y=leveltp), colour="grey30") +
  geom_point(aes(x=r, y=leveltp, fill=score), shape=21, colour="grey30", size=3) +
  scale_fill_viridis_c(name="Score threshold") +
  xlab("Recall") +
  scale_y_continuous(name="Loudness ranking", breaks=c(5, 6, 7, 8),
                     labels = c("Medium", "Medium-loud", "Loud", "Very loud")) +
  my.theme + 
  theme(axis.text.y=element_text(angle=45),
        legend.position="none",
        axis.title.y = element_text(vjust = -5),
        plot.margin = margin(l=-10))
#plot.level


ggsave(grid.arrange(plot.level, plot.rsl, legend, nrow=3, heights=c(4,4,1)), file="Figures/FigA3.1Loudness.jpeg", device="jpeg", height = 9, width = 6, unit="in")
