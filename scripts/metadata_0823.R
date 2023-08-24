# J. Beatty 1/31/23

### Metadata analysis
library(plyr)
library(corrplot)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(patchwork)

######depth profiles ########

# Using the original document from scope website
metadata_fr <- read.csv("~/github/mesoscope2017/raw-data/mesoscope_complete_metadata.csv", header=T)
head(metadata_fr)
# has all the stations, I only want a few, and I want to rename the columns and take the mean if necessary per depth
metadata_stn <- filter(metadata_fr, (Station==4 | Station==6 | Station==8 | Station==10 | Station==12 | Station==14)) %>%
  select(., -c("Time","RosPos","Alkalinity","pH")) %>%
  group_by(Station, Cast, Depth) %>% 
  summarise(lat=mean(Latitude,na.rm=T),long=mean(Longitude,na.rm=T),
           CTD_Temp=mean(CTD_Temperature, na.rm=T),CTD_Salinity=mean(CTD_Salinity, na.rm=T),CTD_Oxygen=mean(CTD_Oxygen, na.rm=T),CTD_Chloropigment=mean(CTD_Chloropigment, na.rm=T),
          Potential_Temp=mean(Potential_Temperature, na.rm=T),Potential_Density=mean(Potential_Density, na.rm=T),Bottle_O2=mean(Bottle_Oxygen, na.rm=T),
          Si =mean(SiO4, na.rm=T), PSi = mean(PSi, na.rm=T),
          DIC=mean(DIC, na.rm=T),PC=mean(PC, na.rm=T),PN=mean(PN, na.rm=T),NO3_NO2=mean(NO3.NO2, na.rm=T),LLN=mean(LLN, na.rm=T),PP = mean(PP, na.rm=T),PO4=mean(PO4,na.rm=T),LLP=mean(LLP, na.rm=T),
          Chlorophyll=mean(Chlorophyll, na.rm=T), Pheopigment=mean(Pheopigment, na.rm=T),
          Het_Bac=mean(Heterotrophic_Bacteria,na.rm=T),Pro=mean(Prochlorococcus,na.rm = T),Syn=mean(Synechococcus,na.rm=T),Euk=mean(Eukaryotes,na.rm = T)) %>% 
  mutate(Stn = sub("^","S",Station)) # Need to add an S before Station so that we can combine the key_all info
head(metadata_stn)

write.csv(metadata_stn, "~/github/mesoscope2017/processed-data/metadata_bystation.csv")



metadata.m <- melt(metadata_stn, id.vars=c("Stn","Depth") , na.rm=TRUE)
metadata.m<-metadata.m[order(metadata.m[1],metadata.m[2]),]
tail(metadata.m)

metadata.m$Stn<-factor(as.character(metadata.m$Stn),levels=c("S4","S6","S8","S10","S12","S14"))
 
color_sla<-c("S6"="#982c1d", "S8"="#bd6a59","S14"="#e8c1b8","S10"="#d2cdca","S4"="#e4dcf4", "S12"="#8e76cf")                     
# put in a melted dataframe for plotting, it will make a plot for each of the different variables
meta_plot<- function(df){
  topics <- unique(df$variable)
  plots = list()
  i = 1
  for (topic in topics){
    df_2<-filter(df, variable==topic)
    plot<-ggplot(df_2,aes(x=value, y=Depth,color=Stn)) + 
      geom_path(linewidth=0.5) + # geom_point()+ 
      scale_y_reverse()+scale_x_continuous(position = "top") + 
      labs(x = topic, y="Depth (m)")+scale_color_manual(values=color_sla) +
      theme_bw() + theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"), text = element_text(size=8))+
      theme(
        legend.position = c(.05, .99),
        legend.justification = c("left", "top"),
        legend.box.just = "right")
    plots[[i]]=plot
    i= i+1
  }
  return(plots)
}

plots <- meta_plot(metadata.m)
unique(metadata.m$variable)

(plots[[5]] + plots[[8]] + plots[[10]] + plots[[11]])+ plot_annotation(tag_levels = "a")
ggplot2::ggsave("~/github/mesoscope2017/figures/metadata_0823.pdf", width=5, height = 6, units="in", dpi=300)

