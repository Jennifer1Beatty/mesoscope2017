# J. Beatty 1/31/23

### Metadata analysis
library(plyr)
library(corrplot)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(patchwork)

key_all<-read.delim("~/Desktop/Chapter1/MESOSCOPE2017/Data/variables_all.txt")
key_all$Station<-factor(as.character(key_all$Station),levels=c("S4","S6","S8","S10","S12","S14",
                                                               "PIT12","PIT11","PIT10","PIT9","PIT8","PIT7","PIT6","PIT5","PIT4","PIT3","PIT2","PIT1"))
key_all$Station.sla<-factor(as.character(key_all$Station), levels=rev(c("S6", "S8", "S14","S10","S4","S12",
                                                                        "PIT12","PIT11","PIT10","PIT9","PIT8","PIT7","PIT6","PIT5","PIT4","PIT3","PIT2","PIT1")))
key_all$Loc <- factor(as.character(key_all$Loc), levels=c("AC-Outside","AC-Center","AC-Edge","Front","C-Center","C-Outside"))
key_all$Color.Loc <- factor(as.character(key_all$Color.Loc), levels=c("#bdbdbd","#de2d26","#969696","#737373","#2c7fb8"))
key_all$Color.sla <- factor(as.character(key_all$Color.sla), levels = c("#7ea6c3","#de2d26","#d77773","#cfcdcd","#2c7fb8", "#d2acab"))
key_all$Depth<-factor(as.character(key_all$Depth),levels=c("PIT","15m","DCM","175m","500m"))
key_all$Material<-factor(key_all$Material,levels=c("DNA","RNA","PIT"))
key_all$Color.depth<-factor(as.character(key_all$Color.depth),levels=c("#d55e00","#f0e442","#009e73","#cc79a7","#0072b2"))
key_all$Shape.type<-factor(as.numeric(key_all$Shape.type),levels=c("1","5"))
key_all$Shape.mat<-factor(as.numeric(key_all$Shape.mat),levels=c(25,21,22))
key_all$Line.mat<-factor(as.numeric(key_all$Line.mat),levels=c("6", "1", "3"))
summary(key_all)

setup_df <- select(key_all, c("Station","Depth","Material","Station.sla"))

metadata <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/metadata_edited.csv", header=T)
metadata_df <- select(metadata, -c("Station", "Depth_var"))
####### 1. Correlation plot ########
# http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram#computing-the-p-value-of-correlations 

M= cor(metadata_df, use="na.or.complete")
corrplot(M, type="upper")

# adding sig figs
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(metadata_df)
head(p.mat[, 1:5])

# plot with sig figs
corrplot(M, type="upper", p.mat = p.mat, sig.level = 0.01, insig = "blank")

## subsetting for the significant ones
metadata_sig <- select(metadata_df, -c("SLA", "NO3_NO2","LLN"))
M_2 = cor(metadata_sig, use="na.or.complete")
p.mat_2 = cor.mtest(metadata_sig)
corrplot(M_2, type="upper", p.mat = p.mat_2, sig.level = 0.01, insig = "blank")


###### 2. depth profiles ########

# Using the original document from scope website
metadata_fr <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/mesoscope_complete_metadata.csv", header=T)
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

write.csv(metadata_stn, "metadata_bystation.csv")

color_sla<-c("S6"="#982c1d", "S8"="#bd6a59","S14"="#e8c1b8","S10"="#d2cdca","S4"="#e4dcf4", "S12"="#8e76cf")
metadata_stn$Stn <- factor(metadata_stn$Stn, levels=c("S4","S6", "S8","S10","S12", "S14")) # in 05/23 I changed the orcerd back to S4-S14

metadata.m <- melt(metadata_stn, id.vars=c("Stn","Depth") , na.rm=TRUE)
tail(metadata.m)
key_all$Stn = key_all$Station
metadata_plot <- join(metadata.m, select(key_all, c("Stn","Col.sla")), by="Stn")
                      
# put in a melted dataframe for plotting, it will make a plot for each of the different variables
meta_plot<- function(df){
  topics <- unique(df$variable)
  plots = list()
  i = 1
  for (topic in topics){
    df_2<-filter(df, variable==topic)
    plot<-ggplot(df_2,aes(x=value, y=Depth,color=Stn)) + 
      geom_path() + geom_point()+ 
      scale_y_reverse()+scale_x_continuous(position = "top") + 
      labs(x = topic, y="Depth (m)")+scale_color_manual(values=color_sla) +
      theme_bw() + theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
    plots[[i]]=plot
    i= i+1
  }
  return(plots)
}

plots <- meta_plot(metadata_plot)
unique(metadata.m$variable)
ggarrange(plots[[4]],plots[[5]],plots[[6]],plots[[7]])
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/metadata/ctd_variables_0323.pdf", width=10, height = 10)
ggarrange(plots[[8]],plots[[9]])
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/metadata/potential_variables_0323.pdf", width=10, height = 10)
ggarrange(plots[[11]],plots[[12]])
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/metadata/carbon_0323.pdf", width=10, height = 10)
ggarrange(plots[[13]],plots[[14]],plots[[15]])
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/metadata/nitrogen_0323.pdf", width=10, height = 10)
ggarrange(plots[[16]],plots[[17]])
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/metadata/phosphate_0323.pdf", width=10, height = 10)
ggarrange(plots[[18]],plots[[19]])
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/metadata/pigments_0323.pdf", width=10, height = 10)
ggarrange(plots[[20]],plots[[21]],plots[[22]],plots[[23]])
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/metadata/flowcytometry_0323.pdf", width=10, height = 10)

ggarrange(plots[[4]],plots[[7]],plots[[9]],plots[[14]], common.legend = TRUE, legend="right")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/metadata/topmeta_0323.pdf", width=6, height = 10)

(plots[[4]] + plots[[7]] + plots[[9]] + plots[[14]]) + plot_layout(guides="collect") + plot_annotation(tag_levels = "a")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/metadata/topmeta_0323.pdf", width=8, height = 9)

ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Final_Figs/meta_0523.pdf", width=8, height = 8, units="in", dpi=300)

###### 3. make a table with variables ######
