library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library("RColorBrewer")


#Meso-scope metavariables
all_data<-read.csv("2017 Microbial Ecology of the Surface Ocean - SCOPE, Water Column Data.csv", header=T)
names(all_data)
para<-c("lat","lon","CTD_Temperature","CTD_Salinity","CTD_Oxygen","CTD_Chloropigment","Potential_Temperature", "Potential_Density", "Bottle_Oxygen","DIC","PO4","NO3_NO2","LLN","LLP","PC","PN","Chlorophyll", "Pheopigment",
        "Heterotrophic_Bacteria","Prochlorococcus", "Synechococcus","Eukaryotes")
S4_all = all_data[all_data$Station == "4",] #dim(S4_all) 69 rows and 29 columns
write.csv(S4_all,"S4_all.csv")
S6_all = all_data[all_data$Station == "6",] #dim(S6_all) 72 rows and 29 columns
write.csv(S6_all,"S6_all.csv")
S8_all = all_data[all_data$Station == "8",] #dim(S8_all) 48 rows and 29 columns
write.csv(S8_all,"S8_all.csv")
S10_all = all_data[all_data$Station == "10",] #dim(S10_all) 48 rows and 29 columns
write.csv(S10_all,"S10_all.csv")
S12_all = all_data[all_data$Station == "12",] #dim(S4_all) 44 rows and 29 columns
write.csv(S12_all,"S12_all.csv")
S14_all = all_data[all_data$Station == "14",] #dim(S4_all) 42 rows and 29 columns
write.csv(S14_all,"S14_all.csv")


#I'm going to summarize the data based on depth, so it's one row per depth with all the variables in that row, I have to add the station number separately 
#I manually edited the depths to combine the different niskin bottles depths, I took the average DCM depth from the casts at the station
S4.depth<-read.csv("S4_meta.csv",header=T)
S4.all<-S4.depth %>% group_by(depth) %>% summarise(lat=mean(lat,na.rm=T),lon=mean(lon,na.rm=T),CTD_Temperature=mean(CTD_Temperature, na.rm=T),CTD_Salinity=mean(CTD_Salinity, na.rm=T),CTD_Oxygen=mean(CTD_Oxygen, na.rm=T),CTD_Chloropigment=mean(CTD_Chloropigment, na.rm=T),
                                                   Potential_Temperature=mean(Potential_Temperature, na.rm=T),Potential_Density=mean(Potential_Density, na.rm=T),Bottle_Oxygen=mean(Bottle_Oxygen, na.rm=T),DIC=mean(DIC, na.rm=T),PO4=mean(PO4,na.rm=T),
                                                   NO3_NO2=mean(NO3_NO2, na.rm=T),LLN=mean(LLN, na.rm=T),LLP=mean(LLP, na.rm=T),PC=mean(PC, na.rm=T),Pheopigment=mean(Pheopigment, na.rm=T),PN=mean(PN, na.rm=T),Chlorophyll=mean(Chlorophyll, na.rm=T),
                                                   Heterotrophic_Bacteria=mean(Heterotrophic_Bacteria,na.rm=T),Prochlorococcus=mean(Prochlorococcus,na.rm = T),Synechococcus=mean(Synechococcus,na.rm=T),Eukaryotes=mean(Eukaryotes,na.rm = T))
S4<-rep(4,12)
S4.all$Station<-S4

S6.depth<-read.csv("S6_meta.csv",header=T)
S6.all<-S6.depth %>% group_by(depth) %>% summarise(lat=mean(lat,na.rm=T),lon=mean(lon,na.rm=T),CTD_Temperature=mean(CTD_Temperature, na.rm=T),CTD_Salinity=mean(CTD_Salinity, na.rm=T),CTD_Oxygen=mean(CTD_Oxygen, na.rm=T),CTD_Chloropigment=mean(CTD_Chloropigment, na.rm=T),
                                                   Potential_Temperature=mean(Potential_Temperature, na.rm=T),Potential_Density=mean(Potential_Density, na.rm=T),Bottle_Oxygen=mean(Bottle_Oxygen, na.rm=T),DIC=mean(DIC, na.rm=T),PO4=mean(PO4,na.rm=T),
                                                   NO3_NO2=mean(NO3_NO2, na.rm=T),LLN=mean(LLN, na.rm=T),LLP=mean(LLP, na.rm=T),PC=mean(PC, na.rm=T),Pheopigment=mean(Pheopigment, na.rm=T),PN=mean(PN, na.rm=T),Chlorophyll=mean(Chlorophyll, na.rm=T),
                                                   Heterotrophic_Bacteria=mean(Heterotrophic_Bacteria,na.rm=T),Prochlorococcus=mean(Prochlorococcus,na.rm = T),Synechococcus=mean(Synechococcus,na.rm=T),Eukaryotes=mean(Eukaryotes,na.rm = T))
S6<-rep(6,12)
S6.all$Station<-S6

S8.depth<-read.csv("S8_meta.csv",header=T)
S8.all<-S8.depth %>% group_by(depth) %>% summarise(lat=mean(lat,na.rm=T),lon=mean(lon,na.rm=T),CTD_Temperature=mean(CTD_Temperature, na.rm=T),CTD_Salinity=mean(CTD_Salinity, na.rm=T),CTD_Oxygen=mean(CTD_Oxygen, na.rm=T),CTD_Chloropigment=mean(CTD_Chloropigment, na.rm=T),
                                                   Potential_Temperature=mean(Potential_Temperature, na.rm=T),Potential_Density=mean(Potential_Density, na.rm=T),Bottle_Oxygen=mean(Bottle_Oxygen, na.rm=T),DIC=mean(DIC, na.rm=T),PO4=mean(PO4,na.rm=T),
                                                   NO3_NO2=mean(NO3_NO2, na.rm=T),LLN=mean(LLN, na.rm=T),LLP=mean(LLP, na.rm=T),PC=mean(PC, na.rm=T),Pheopigment=mean(Pheopigment, na.rm=T),PN=mean(PN, na.rm=T),Chlorophyll=mean(Chlorophyll, na.rm=T),
                                                   Heterotrophic_Bacteria=mean(Heterotrophic_Bacteria,na.rm=T),Prochlorococcus=mean(Prochlorococcus,na.rm = T),Synechococcus=mean(Synechococcus,na.rm=T),Eukaryotes=mean(Eukaryotes,na.rm = T))
S8<-rep(8,12)
S8.all$Station<-S8

S10.depth<-read.csv("S10_meta.csv",header=T)
S10.all<-S10.depth %>% group_by(depth) %>% summarise(lat=mean(lat,na.rm=T),lon=mean(lon,na.rm=T),CTD_Temperature=mean(CTD_Temperature, na.rm=T),CTD_Salinity=mean(CTD_Salinity, na.rm=T),CTD_Oxygen=mean(CTD_Oxygen, na.rm=T),CTD_Chloropigment=mean(CTD_Chloropigment, na.rm=T),
                                                   Potential_Temperature=mean(Potential_Temperature, na.rm=T),Potential_Density=mean(Potential_Density, na.rm=T),Bottle_Oxygen=mean(Bottle_Oxygen, na.rm=T),DIC=mean(DIC, na.rm=T),PO4=mean(PO4,na.rm=T),
                                                   NO3_NO2=mean(NO3_NO2, na.rm=T),LLN=mean(LLN, na.rm=T),LLP=mean(LLP, na.rm=T),PC=mean(PC, na.rm=T),Pheopigment=mean(Pheopigment, na.rm=T),PN=mean(PN, na.rm=T),Chlorophyll=mean(Chlorophyll, na.rm=T),
                                                   Heterotrophic_Bacteria=mean(Heterotrophic_Bacteria,na.rm=T),Prochlorococcus=mean(Prochlorococcus,na.rm = T),Synechococcus=mean(Synechococcus,na.rm=T),Eukaryotes=mean(Eukaryotes,na.rm = T))
S10<-rep(10,12)
S10.all$Station<-S10

S12.depth<-read.csv("S12_meta.csv",header=T)
S12.all<-S12.depth %>% group_by(depth) %>% summarise(lat=mean(lat,na.rm=T),lon=mean(lon,na.rm=T),CTD_Temperature=mean(CTD_Temperature, na.rm=T),CTD_Salinity=mean(CTD_Salinity, na.rm=T),CTD_Oxygen=mean(CTD_Oxygen, na.rm=T),CTD_Chloropigment=mean(CTD_Chloropigment, na.rm=T),
                                                   Potential_Temperature=mean(Potential_Temperature, na.rm=T),Potential_Density=mean(Potential_Density, na.rm=T),Bottle_Oxygen=mean(Bottle_Oxygen, na.rm=T),DIC=mean(DIC, na.rm=T),PO4=mean(PO4,na.rm=T),
                                                   NO3_NO2=mean(NO3_NO2, na.rm=T),LLN=mean(LLN, na.rm=T),LLP=mean(LLP, na.rm=T),PC=mean(PC, na.rm=T),Pheopigment=mean(Pheopigment, na.rm=T),PN=mean(PN, na.rm=T),Chlorophyll=mean(Chlorophyll, na.rm=T),
                                                   Heterotrophic_Bacteria=mean(Heterotrophic_Bacteria,na.rm=T),Prochlorococcus=mean(Prochlorococcus,na.rm = T),Synechococcus=mean(Synechococcus,na.rm=T),Eukaryotes=mean(Eukaryotes,na.rm = T))
S12<-rep(12,12)
S12.all$Station<-S12

S14.depth<-read.csv("S14_meta.csv",header=T)
S14.all<-S14.depth %>% group_by(depth) %>% summarise(lat=mean(lat,na.rm=T),lon=mean(lon,na.rm=T),CTD_Temperature=mean(CTD_Temperature, na.rm=T),CTD_Salinity=mean(CTD_Salinity, na.rm=T),CTD_Oxygen=mean(CTD_Oxygen, na.rm=T),CTD_Chloropigment=mean(CTD_Chloropigment, na.rm=T),
                                                   Potential_Temperature=mean(Potential_Temperature, na.rm=T),Potential_Density=mean(Potential_Density, na.rm=T),Bottle_Oxygen=mean(Bottle_Oxygen, na.rm=T),DIC=mean(DIC, na.rm=T),PO4=mean(PO4,na.rm=T),
                                                   NO3_NO2=mean(NO3_NO2, na.rm=T),LLN=mean(LLN, na.rm=T),LLP=mean(LLP, na.rm=T),PC=mean(PC, na.rm=T),Pheopigment=mean(Pheopigment, na.rm=T),PN=mean(PN, na.rm=T),Chlorophyll=mean(Chlorophyll, na.rm=T),
                                                   Heterotrophic_Bacteria=mean(Heterotrophic_Bacteria,na.rm=T),Prochlorococcus=mean(Prochlorococcus,na.rm = T),Synechococcus=mean(Synechococcus,na.rm=T),Eukaryotes=mean(Eukaryotes,na.rm = T))
S14<-rep(14,12)
S14.all$Station<-S14

#To combine all the dataframes back together
Sall.all<-bind_rows(S4.all,S6.all,S8.all,S10.all,S12.all,S14.all)
write.csv(Sall.all,"AllS_Meta_curated.csv")

###TO start with saved datafame
Sall.all<-read.csv("AllS_Meta_curated.csv",header=T)
Sall.all[1]=NULL

N<-c("depth","S4","Anti-Cyclone","S8","S10","Cyclone","S14")
Stats=c("S4","S6","S8","S10","S12","S14")
Col= c("#bdbdbd","#de2d26","#969696", "#737373","#525252","#2c7fb8","#252525")
label.all<-as.character(Stats)
colScale<-scale_color_manual(values=Col)
Sall.all$order<-factor(Sall.all$Station, levels = label.all)

###To plot- 
#from Sam G
#Chlorophyll plot
chl= ggplot(Sall.all, aes(x=CTD_Chloropigment, y=depth, color=Station)) +
  geom_point(size=2)+geom_path(size=1)+
  scale_y_reverse()+scale_x_continuous(position = "top")+theme_classic()+
  labs(x = expression("Chlorophyll ("*mu*"g L"^{"-1"}*")", y="Depth (m)"))
print(chl)
ggsave("chl_eddies.pdf", height=10, width=5)

Temp<-dcast(Sall.all[c(1,4,24)], depth~Station, value.var="CTD_Temperature") #restructure dataframe to only include depth oxygen
names(Temp)=N
Temp.mm<-reshape2::melt(Temp, id.var='depth')
Temperature<-ggplot(Temp.mm,aes(x=value, y=depth, col=variable)) + 
  geom_line() + geom_point(size=2)+ 
  scale_y_reverse()+scale_x_continuous(position = "top") + 
  labs(x = expression("Temperature ("*degree*"C)", y="Depth (m)"))+scale_color_manual(labels=c("S4","Anti-Cyclone","S8","S10","Cyclone","S14"), values=c("#bdbdbd","#de2d26","#969696", "#737373","#2c7fb8","#252525")) + 
  theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
Temperature
ggsave("ctd_temp.pdf", height=10, width=5)

Chla<-dcast(Sall.all[c(1,7,24)], depth~Station, value.var="CTD_Chloropigment") #restructure dataframe to only include depth oxygen
names(Chla)=N
Chla.mm<-reshape2::melt(Chla, id.var='depth')
Chla.plot<-ggplot(Chla.mm,aes(x=value, y=depth,color=variable)) + 
  geom_path() + geom_point()+ 
  scale_y_reverse()+scale_x_continuous(position = "top") + 
  labs(x = expression("Chlorophyll ("*mu*"g L"^{"-1"}*")", y="Depth (m)"))+scale_color_manual(labels=c("S4","Anti-Cyclone","S8","S10","Cyclone","S14"), values=c("#bdbdbd","#de2d26","#969696", "#737373","#2c7fb8","#252525"))
  theme_bw() + theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
Chla.plot
ggsave("ctd_chla.pdf", height=10, width=5)

O2<-dcast(Sall.all[c(1,6,24)], depth~Station, value.var="CTD_Oxygen") #restructure dataframe to only include depth oxygen
names(O2)=N
O2.mm<-reshape2::melt(O2, id.var='depth')
O2.plot<-ggplot(O2.mm,aes(x=value, y=depth,color=variable)) + 
  geom_path() + geom_point()+ 
  scale_y_reverse()+scale_x_continuous(position = "top") + 
  labs(x = expression("Oxygen ("*mu*"mol kg"^{"-1"}*")", y="Depth (m)"))+scale_color_manual(labels=c("S4","Anti-Cyclone","S8","S10","Cyclone","S14"), values=c("#bdbdbd","#de2d26","#969696", "#737373","#2c7fb8","#252525"))+
  theme_bw() + theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
O2.plot
ggsave("ctd_O2.pdf", height=10, width=5)

DIC<-dcast(Sall.all[c(1,11,24)],depth~Station, value.var="DIC") 
DIC[2]=NULL #No data for S4
names(DIC)=c("depth","S6","S8","S10","S12","S14")
DIC.mm<-reshape2::melt(DIC, id.var='depth')
DIC.plot<-ggplot(DIC.mm,aes(x=value, y=depth, col=variable)) + 
  geom_line() + geom_point()+ 
  scale_y_reverse()+scale_x_continuous(position = "top") + 
  labs(x = expression("Dissolved Inorganic Carbon ("*mu*"mol kg"^{"-1"}*")", y="Depth (m)"))+scale_color_manual(labels=c("Anti-Cyclone","S8","S10","Cyclone","S14"), values=c("#de2d26","#969696", "#737373","#2c7fb8","#252525")) + 
  theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
DIC.plot
ggsave("DIC.pdf", height=10, width=5)

NO<-dcast(Sall.all[c(1,13,24)],depth~Station, value.var="NO3_NO2")
names(NO)=N
NO.mm<-reshape2::melt(NO, id.var='depth')
NO.plot<-ggplot(NO.mm,aes(x=value, y=depth, col=variable)) + 
  geom_line() + geom_point()+ 
  scale_y_reverse()+scale_x_continuous(position = "top") + 
  labs(x = expression("Nitrate plus Nitrite ("*mu*"mol kg"^{"-1"}*")", y="Depth (m)"))+scale_color_manual(labels=c("S4","Anti-Cyclone","S8","S10","Cyclone","S14"), values=c("#bdbdbd","#de2d26","#969696", "#737373","#2c7fb8","#252525")) + 
  theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
NO.plot
ggsave("NO3_NO2.pdf", height=10, width=5)

PO4<-dcast(Sall.all[c(1,12,24)],depth~Station, value.var="PO4")
names(PO4)=N
PO4.mm<-reshape2::melt(PO4, id.var='depth')
PO4.plot<-ggplot(PO4.mm,aes(x=value, y=depth, col=variable)) + 
  geom_line() + geom_point()+ 
  scale_y_reverse()+scale_x_continuous(position = "top") + 
  labs(x = expression("Phosphate ("*mu*"mol kg"^{"-1"}*")", y="Depth (m)"))+scale_color_manual(labels=c("S4","Anti-Cyclone","S8","S10","Cyclone","S14"), values=c("#bdbdbd","#de2d26","#969696", "#737373","#2c7fb8","#252525")) + 
  theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
PO4.plot
ggsave("PO4.pdf", height=10, width=5)

LLN<-dcast(Sall.all[c(1,14,24)],depth~Station, value.var="LLN")
names(LLN)=N
LLN.mm<-reshape2::melt(LLN, id.var='depth')
LLN.plot<-ggplot(LLN.mm,aes(x=value, y=depth, col=variable)) + 
   geom_point()+ 
  scale_y_reverse()+scale_x_continuous(position = "top") + 
  labs(x = expression("Low-level Nitrogen ("*mu*"mol kg"^{"-1"}*")", y="Depth (m)"))+scale_color_manual(labels=c("S4","Anti-Cyclone","S8","S10","Cyclone","S14"), values=c("#bdbdbd","#de2d26","#969696", "#737373","#2c7fb8","#252525")) + 
  theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
LLN.plot
ggsave("LLN.pdf", height=10, width=5)

LLP<-dcast(Sall.all[c(1,15,24)],depth~Station, value.var="LLP")
names(LLP)=N
LLP.mm<-reshape2::melt(LLP, id.var='depth')
LLP.plot<-ggplot(LLP.mm,aes(x=value, y=depth, col=variable)) + 
  geom_path() +geom_point()+ 
  scale_y_reverse()+scale_x_continuous(position = "top") + 
  labs(x = expression("Low-level Phosphate ("*mu*"mol kg"^{"-1"}*")", y="Depth (m)"))+scale_color_manual(labels=c("S4","Anti-Cyclone","S8","S10","Cyclone","S14"), values=c("#bdbdbd","#de2d26","#969696", "#737373","#2c7fb8","#252525")) + 
  theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
LLP.plot
ggsave("LLP.pdf", height=10, width=5)

PC<-dcast(Sall.all[c(1,16,24)],depth~Station, value.var="PC")
names(PC)=N
PC.mm<-reshape2::melt(PC, id.var='depth')
PC.plot<-ggplot(PC.mm,aes(x=value, y=depth, col=variable)) + 
  geom_point()+ 
  scale_y_reverse()+scale_x_continuous(position = "top") + 
  labs(x = expression("Particulate Carbon ("*mu*"mol kg"^{"-1"}*")", y="Depth (m)"))+scale_color_manual(labels=c("S4","Anti-Cyclone","S8","S10","Cyclone","S14"), values=c("#bdbdbd","#de2d26","#969696", "#737373","#2c7fb8","#252525")) + 
  theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
PC.plot
ggsave("PC.pdf", height=10, width=5)

PN<-dcast(Sall.all[c(1,18,24)],depth~Station, value.var="PN")
names(PN)=N
PN.mm<-reshape2::melt(PN, id.var='depth')
PN.plot<-ggplot(PN.mm,aes(x=value, y=depth, col=variable)) + 
  geom_point()+ 
  scale_y_reverse()+scale_x_continuous(position = "top") + 
  labs(x = expression("Particulate Nitrogen ("*mu*"mol kg"^{"-1"}*")", y="Depth (m)"))+scale_color_manual(labels=c("S4","Anti-Cyclone","S8","S10","Cyclone","S14"), values=c("#bdbdbd","#de2d26","#969696", "#737373","#2c7fb8","#252525")) + 
  theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
PN.plot
ggsave("PN.pdf", height=10, width=5)

Pheo<-dcast(Sall.all[c(1,17,24)],depth~Station, value.var="Pheopigment")
names(Pheo)=N
Pheo.mm<-reshape2::melt(Pheo, id.var='depth')
Pheo.plot<-ggplot(Pheo.mm,aes(x=value, y=depth, col=variable)) + 
  geom_path()+geom_point()+ 
  scale_y_reverse()+scale_x_continuous(position = "top") + 
  labs(x = expression("Pheopigment ("*mu*"mol kg"^{"-1"}*")", y="Depth (m)"))+scale_color_manual(labels=c("S4","Anti-Cyclone","S8","S10","Cyclone","S14"), values=c("#bdbdbd","#de2d26","#969696", "#737373","#2c7fb8","#252525")) + 
  theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
Pheo.plot
ggsave("Pheo.pdf", height=10, width=5)

Chl<-dcast(Sall.all[c(1,19,24)],depth~Station, value.var="Chlorophyll")
names(Chl)=N
Chl.mm<-reshape2::melt(Chl, id.var='depth')
Chl.plot<-ggplot(Chl.mm,aes(x=value, y=depth, col=variable)) + 
  geom_path()+geom_point()+ 
  scale_y_reverse()+scale_x_continuous(position = "top") + 
  labs(x = expression("Chlorophyll ("*mu*"mol kg"^{"-1"}*")", y="Depth (m)"))+scale_color_manual(labels=c("S4","Anti-Cyclone","S8","S10","Cyclone","S14"), values=c("#bdbdbd","#de2d26","#969696", "#737373","#2c7fb8","#252525")) + 
  theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
Chl.plot
ggsave("Chl.pdf", height=10, width=5)

Pro<-dcast(Sall.all[c(1,21,24)],depth~Station, value.var="Prochlorococcus")
names(Pro)=N
Pro.mm<-reshape2::melt(Pro, id.var='depth')
Pro.plot<-ggplot(Pro.mm,aes(x=value, y=depth, col=variable)) + 
  geom_path()+geom_point()+ 
  scale_y_reverse()+scale_x_continuous(position = "top") + 
  labs(x = expression("Prochlorococcus (# ml"^{"-1"}*")", y="Depth (m)"))+scale_color_manual(labels=c("S4","Anti-Cyclone","S8","S10","Cyclone","S14"), values=c("#bdbdbd","#de2d26","#969696", "#737373","#2c7fb8","#252525")) + 
  theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
Pro.plot
ggsave("Pro.pdf", height=10, width=5)

Syn<-dcast(Sall.all[c(1,22,24)],depth~Station, value.var="Synechococcus")
names(Syn)=N
Syn.mm<-reshape2::melt(Syn, id.var='depth')
Syn.plot<-ggplot(Syn.mm,aes(x=value, y=depth, col=variable)) + 
  geom_path()+geom_point()+ 
  scale_y_reverse()+scale_x_continuous(position = "top") + 
  labs(x = expression("Synechococcus (# ml"^{"-1"}*")", y="Depth (m)"))+scale_color_manual(labels=c("S4","Anti-Cyclone","S8","S10","Cyclone","S14"), values=c("#bdbdbd","#de2d26","#969696", "#737373","#2c7fb8","#252525")) + 
  theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
Syn.plot
ggsave("Syn.pdf", height=10, width=5)

Hbac<-dcast(Sall.all[c(1,20,24)],depth~Station, value.var="Heterotrophic_Bacteria")
names(Hbac)=N
Hbac.mm<-reshape2::melt(Hbac, id.var='depth')
Hbac.plot<-ggplot(Hbac.mm,aes(x=value, y=depth, col=variable)) + 
  geom_path(directi)+geom_point()+ 
  scale_y_reverse()+scale_x_continuous(position = "top") + 
  labs(x = expression("Heterotrophic_Bacteria (# ml"^{"-1"}*")", y="Depth (m)"))+scale_color_manual(labels=c("S4","Anti-Cyclone","S8","S10","Cyclone","S14"), values=c("#bdbdbd","#de2d26","#969696", "#737373","#2c7fb8","#252525")) + 
  theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
Hbac.plot
ggsave("Hbac.pdf", height=10, width=5)

Euk<-dcast(Sall.all[c(1,23,24)],depth~Station, value.var="Eukaryotes")
names(Euk)=N
Euk.mm<-reshape2::melt(Euk, id.var='depth')
Euk.plot<-ggplot(Euk.mm,aes(x=value, y=depth, col=variable)) + 
    geom_path()+geom_point()+ 
   scale_y_reverse()+scale_x_continuous(position = "top") + 
   labs(x = expression("Eukaryotes (# ml"^{"-1"}*")", y="Depth (m)"))+scale_color_manual(labels=c("S4","Anti-Cyclone","S8","S10","Cyclone","S14"), values=c("#bdbdbd","#de2d26","#969696", "#737373","#2c7fb8","#252525")) + 
    theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
Euk.plot
ggsave("Euk.pdf", height=10, width=5)


###Trying to just flip the x/y axis to see if that will make plotting better
Euk.plot<-ggplot(Euk.mm,aes(x=depth, y=value, col=variable)) + 
  geom_path()+geom_point(na.rm=TRUE)+ coord_flip()+
  scale_x_reverse()+scale_y_continuous(position = "top") + 
  labs(x = expression("Eukaryotes (# ml"^{"-1"}*")", y="Depth (m)"))+scale_color_manual(labels=c("S4","Anti-Cyclone","S8","S10","Cyclone","S14"), values=c("#bdbdbd","#de2d26","#969696", "#737373","#2c7fb8","#252525")) + 
  theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
Euk.plot


Prot.plot<-ggplot(Pro.mm,aes(x=depth, y=value, col=variable)) + 
  geom_line(direction=hv)+geom_point()+ coord_flip()+
  scale_x_reverse()+ sc
  labs(y = expression("Prochlorococcus (# ml"^{"-1"}*")", x="Depth (m)"))+scale_color_manual(labels=c("S4","Anti-Cyclone","S8","S10","Cyclone","S14"), values=c("#bdbdbd","#de2d26","#969696", "#737373","#2c7fb8","#252525")) + 
  theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
Prot.plot
ggsave("Prot.pdf", height=10, width=5)

