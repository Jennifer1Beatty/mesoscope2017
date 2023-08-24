### 6/5/23 MAPS
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringi)  # for manipulating strings
library(maps)
library(cowplot)

##### fix the PIT data ####
pit1 <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/pit1.csv", header=F)
names(pit1)=c("Time","Lat","Long","Station")
pit2 <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/pit2.csv", header=F)

names(pit2)=c("Time","Lat","Long")
pit2$Station="PIT2"
pit3 <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/pit3.csv", header=F)
names(pit3)=c("Time","Lat","Long")
pit3$Station="PIT3"
pit4 <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/pit4.csv", header=F)
names(pit4)=c("Time","Lat","Long")
pit4$Station="PIT4"
pit5 <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/pit5.csv", header=F)
names(pit5)=c("Time","Lat","Long")
pit5$Station="PIT5"
pit6 <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/pit6.csv", header=F)
names(pit6)=c("Time","Lat","Long")
pit6$Station="PIT6"
pit7 <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/pit7.csv", header=F)
names(pit7)=c("Time","Lat","Long")
pit7$Station="PIT7"
pit8 <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/pit8.csv", header=F)
names(pit8)=c("Time","Lat","Long")
pit8$Station="PIT8"
pit9 <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/pit9.csv", header=F)
names(pit9)=c("Time","Lat","Long")
pit9$Station="PIT9"
pit10 <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/pit10.csv", header=F)
names(pit10)=c("Time","Lat","Long")
pit10$Station="PIT10"
pit11 <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/pit11.csv", header=F)
names(pit11)=c("Time","Lat","Long")
pit11$Station="PIT11"
pit12 <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/pit12.csv", header=F)
names(pit12)=c("Time","Lat","Long")
pit12$Station="PIT12"

pit_tracks<-rbind(pit1,pit2,pit3,pit4,pit5,pit6,pit7,pit8,pit9,pit10,pit11,pit12)
write.csv(pit_tracks,"~/Desktop/Chapter1/MESOSCOPE2017/Data/pit_tracks.csv")

####### input the pit tracks, station points #####
pit_tracks <-read.csv("~/github/mesoscope2017/raw-data/pit_tracks.csv", header=T)
pit_tracks$Station<-factor(as.character(pit_tracks$Station),levels=c("PIT12","PIT11","PIT10","PIT9","PIT8","PIT7","PIT6","PIT5","PIT4","PIT3","PIT2","PIT1"))
pit_points<-read.csv("~/github/mesoscope2017/raw-data/pit_points.csv", header=T)
pit_points$Station<-factor(as.character(pit_points$Station),levels=c("PIT12","PIT11","PIT10","PIT9","PIT8","PIT7","PIT6","PIT5","PIT4","PIT3","PIT2","PIT1"))
stn_points <-read.csv("~/github/mesoscope2017/raw-data/station_points.csv", header=T)
stn_points$Station<-factor(as.character(stn_points$Station),levels=c("S4","S6","S8","S10","S12","S14","Stn. ALOHA"))

# input sla information
sla_df_lat <- read.csv("~/github/mesoscope2017/raw-data/SLAcorr_20170701.csv", header=F) 
sla_long <- read.csv("~/github/mesoscope2017/raw-data/SLA_long.csv", header=F)
# edit it
names(sla_df_lat) = unlist(sla_long)
sla_tidy <- as_tibble(sla_df_lat, .name_repair = "universal") %>% 
  pivot_longer(!Lat, names_to = "Long", values_to = "SLA")
sla_tidy$Long <- as.numeric(stri_replace_all_fixed(sla_tidy$Long, "..", ""))
sla_tidy_subset <- sla_tidy %>%
  subset(., ((Long < 161 & Long > 156.5) & (Lat < 28 & Lat > 23)))

####### to just plot Hawaii #######
mp1 <- fortify(map(fill=TRUE, plot=FALSE))
mp2 <- mp1
mp2$long <- mp2$long
mp2$group <- mp2$group + max(mp2$group) + 1
mp <- rbind(mp1, mp2)
stn_ALOHA=subset(stn_points,Station=="Stn. ALOHA")

hawaii_map<-ggplot()+
  geom_path(data=mp, aes(x=long, y= lat, group = group),color="black")+
  scale_x_continuous(limits=c(-161,-150),breaks=c(-160,-158,-156,-154,-152,-150),label=c("160°W","158°W","156°W","154°W","152°W","150°W"))+
  scale_y_continuous(limits=c(18,28), breaks=c(18,20,22,24,26,28),label=c("18°N","20°N","22°N","24°N","26°N","28°N"))+
  geom_point(data=stn_ALOHA,aes(x=-Long,y=Lat), shape=8, color="black", size=5)+
  geom_text(data=stn_ALOHA, aes(x=-Long+.5,y=Lat-.25,label=Station))+
  geom_rect(aes(ymin = 23.25, ymax = 28, xmin = -156.5, xmax = -160.5,),alpha = .2,color="black")+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(linetype=2,colour = "grey"))
hawaii_map

###### to plot map of stations with SLA contour #####

subset_map<-ggplot() +
  geom_raster(data=sla_tidy_subset,aes(x=Long, y=Lat, z= SLA,fill = SLA), interpolate=TRUE) +
  scale_x_continuous(trans="reverse",limits=c(160.5,156.5),breaks=c(160,159,158,157),label=c("160°W","159°W","158°W","157°W"),)+
  scale_y_continuous(limits=c(23.25,28), breaks=c(24,25,26,27,28),label=c("24°N","25°N","26°N","27°N","28°N"))+
  scale_fill_gradient2(low = "#0818A8",
                       mid = "light gray",
                       high = "#880808",
                       midpoint = 0,
                       limits=c(-28,28))+
  geom_path(data=pit_tracks, aes(x=Long, y=Lat, group=Station),color="black",alpha=0.7)+
  geom_point(data=pit_points, aes(x=Long,y=Lat),size=2,shape=20,color="black",alpha=0.7)+
  geom_point(data=stn_points, aes(x=Long,y=Lat), color="black",shape=18,size=3)+
  geom_text(data=stn_points,aes(x=Long-.15, y=Lat+0.12,label=Station),color="black", size=2)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.margin=unit(c(0,0,0,0), "pt"))+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour="white"))
subset_map

#### put map together #####
plot.with.inset <-
  ggdraw() +
  draw_plot(hawaii_map) +
  draw_plot(subset_map, x = .5, y = .35, width = 0.5, height = 0.5)
plot.with.inset

ggsave(filename = "~/github/mesoscope2017/figures/map_combined.pdf", 
       plot = plot.with.inset,
       width = 5, 
       height = 6,
       units = "in",
       dpi = 300)
