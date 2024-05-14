## Figure 10. Flux from microscopy counts of sediment trap material
## J. L. Beatty
## Last updated: 05/2024

# load packages
library(plyr)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(lmodel2)

# for plotting
tax_full.pit<-c("Alveolata-Ciliates","Alveolata-Dinoflagellates","Other","Rhizaria-Foraminiferans","Rhizaria-Other","Rhizaria-Radiolarians","Stramenopila-Diatoms","Stramenopila-Silicoflagellates")
pit_colors <- c("Alveolata-Ciliates"='firebrick4',"Alveolata-Dinoflagellates"='indianred1',"Other"= 'grey',
                "Rhizaria-Foraminiferans"='mediumvioletred',"Rhizaria-Radiolarians"='lightpink',"Rhizaria-Other"='#C291A4',
                "Stramenopila-Diatoms"='#DDAD4B',"Stramenopila-Silicoflagellates"="#5C4033")
key_all<-read.delim("~/github/mesoscope2017/raw-data/variables_all.txt")
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}
# load data
flux<-read.csv("~/github/mesoscope2017/raw-data/CellsPit_Final.csv",header=T)
# rename groups
rownames(flux)=flux$Type
flux$Plotting = flux$Type
flux$Plotting[flux$Type == "Tintinnid"]="Alveolata-Ciliates"
flux$Plotting[flux$Type == "Ciliate"]="Alveolata-Ciliates"
flux$Plotting[flux$Type == "Dinoflagellate"]="Alveolata-Dinoflagellates"
flux$Plotting[flux$Type == "Rhizaria-Formanifera"]="Rhizaria-Formaniferans"
flux$Plotting[flux$Type == "Rhizaria-Radiolaria"]="Rhizaria-Radiolarians"
# combine acanths and rads to one category
flux$Plotting[flux$Type == "Rhizaria-Acantharians"]="Rhizaria-Radiolarians"
flux$Plotting[flux$Type == "Rhizaria-other"]="Rhizaria-Other"
flux$Plotting[flux$Type == "Diatom"]="Stramenopila-Diatoms"
flux$Plotting[flux$Type == "Dictyocha"]="Stramenopila-Silicoflagellates"
# reshape data
flux.m<-melt(select(flux, -c("Type")), id.vars = "Plotting") %>%
  separate(variable, c("Station", "Rep"), sep = "_", remove = TRUE)
head(flux.m)

#To calculate mean values of taxa by PIT
flux.mean<- flux.m %>% 
  group_by(Plotting, Station) %>%
  summarise(mean=mean(as.numeric(value)),stdev=sd(as.numeric(value))) %>%
  as.data.frame
flux.mean

flux.mean$tax.order<-factor(flux.mean$Plotting,levels=(tax_full.pit), labels=(tax_full.pit))
flux_sla <- join(flux.mean, select(key_all, c(Sample, Station, SLA)), by="Station",type="left")

# create a type 2 regression for all the factors
model_type2 <- flux_sla %>%
  group_by(Plotting) %>%
  do(model = lmodel2(mean ~ SLA,data=.))

# rename the model
model_cil <-model_type2[1,]$model[[1]][[3]]
names(model_cil) <- c("method", "intercept", "slope", "angle", "p-value")
model_cil$Plotting<-"Alveolata-Ciliates"
model_dino<-model_type2[2,]$model[[1]][[3]]
names(model_dino) <- c("method", "intercept", "slope", "angle", "p-value")
model_dino$Plotting<-"Alveolata-Dinoflagellates"
model_foram <-model_type2[3,]$model[[1]][[3]]
names(model_foram) <- c("method", "intercept", "slope", "angle", "p-value")
model_foram$Plotting<-"Rhizaria-Foraminiferans"
model_rhizother <-model_type2[4,]$model[[1]][[3]]
names(model_rhizother) <- c("method", "intercept", "slope", "angle", "p-value")
model_rhizother$Plotting<-"Rhizaria-Other"
model_rad <-model_type2[5,]$model[[1]][[3]]
names(model_rad) <- c("method", "intercept", "slope", "angle", "p-value")
model_rad$Plotting<-"Rhizaria-Radiolarians"
model_diatom <-model_type2[6,]$model[[1]][[3]]
names(model_diatom) <- c("method", "intercept", "slope", "angle", "p-value")
model_diatom$Plotting<-"Stramenopila-Diatoms"
model_sili <-model_type2[7,]$model[[1]][[3]]
names(model_sili) <- c("method", "intercept", "slope", "angle", "p-value")
model_sili$Plotting<-"Stramenopila-Silicoflagellates"

model_plot<-rbind(model_cil, model_dino) %>%
  rbind(.,model_foram)%>%
  rbind(.,model_rhizother) %>%
  rbind(., model_rad)%>%
  rbind(., model_diatom) %>%
  rbind(., model_sili)

# R values: 
scatter_sla_model2_rvalue<-ggplot(flux_sla, aes(y=mean,ymin= mean-stdev, ymax= mean+stdev, x=SLA,fill=tax.order))+
  geom_errorbar(width = 0.2) +
  geom_point(size=2, shape=21,color="black")+
  scale_fill_manual(values=pit_colors, name="Taxa")+
  theme_bw()+
  scale_y_continuous(labels = scientific_10)+
  theme(panel.grid.minor = element_line(colour="white"),panel.grid.major = element_line(colour="white")) +
  labs(y=expression("Mean flux by microscopy of trap material (Cells m"^{"-2"}*" d"^{"-1"}*")"), x="SLAcorr(cm)")+
  theme(
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0),
    axis.text.y = element_text(size = 8, angle = 0, hjust = 1, vjust = 0),
    axis.title.y = element_text(size = 10, angle = 90, hjust = .5, vjust = .5),
    text = element_text(size=8),
    strip.background = element_blank(),
    legend.position = c(.8, .05))+
  geom_abline(data = subset(model_plot, method=="SMA"), aes(intercept = intercept, slope = slope)) +
  stat_cor(label.x = .1, label.y = 1500)
scatter_sla_model2_rvalue %+% facet_wrap(~Plotting, ncol = 2, scales="free_y")

#ggplot2::ggsave("~/github/mesoscope2017/figures/fig10.pdf", width=5, height = 6, units="in", dpi=300)
