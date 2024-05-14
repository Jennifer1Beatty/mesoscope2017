# Figure 4- Taxa barplot of trap material
# J. L. Beatty
# Last updated: 05/2024

# load required packages
library(plyr)
library(tidyverse)
library(reshape2)

# load metadata
key_all<-read.delim("~/github/mesoscope2017/raw-data/variables_all.txt")

# info for plotting
tax_color_phyla2_nometa<- c("Alveolata-Ciliates"='firebrick4',"Alveolata-Dinoflagellates"='indianred1',"Alveolata-Syndiniales"='tomato3', "Archaeplastida-Chlorophytes"='forestgreen',
                            "Excavata-Discobids"='yellowgreen', "Hacrobia-Cryptophytes"='darkblue',"Hacrobia-Haptophytes"='lightblue',
                            "Opisthokonta-Choanoflagellates"='moccasin',"Other/unknown"='grey',
                            "Rhizaria-Acantharians"='mediumvioletred',"Rhizaria-Cercozoans"="#C291A4","Rhizaria-Other"='#db6ea0',"Rhizaria-Polycystines"='lightpink', 
                            "Stramenopila-Chrysophytes"="#DECDBE","Stramenopila-Diatoms"='#DDAD4B',"Stramenopila-MAST"='tan2',"Stramenopila-Ochrophytes"='tan3',"Stramenopila-Other"='tan4',"Stramenopila-Pelagophytes"="#5C4033")
tax_full.pit<-c("Alveolata-Ciliates","Alveolata-Dinoflagellates","Other","Rhizaria-Foraminiferans","Rhizaria-Other","Rhizaria-Radiolarians","Stramenopila-Diatoms","Stramenopila-Silicoflagellates")
pit_colors <- c("Alveolata-Ciliates"='firebrick4',"Alveolata-Dinoflagellates"='indianred1',"Other"= 'grey',
                "Rhizaria-Foraminiferans"='mediumvioletred',"Rhizaria-Radiolarians"='lightpink',"Rhizaria-Other"='#C291A4',
                "Stramenopila-Diatoms"='#DDAD4B',"Stramenopila-Silicoflagellates"="#5C4033")

### load asv data
asv_df <- read.csv("~/github/mesoscope2017/processed-data/MS_Clean_No1_Means_ALL_Norm_NewTax_0424.csv", header=T)
asv_df[1]=NULL
row.names(asv_df) = asv_df$ASV1.ID
nometa <- asv_df[asv_df$Level3 != "Metazoa",]

data.m<-melt(nometa) #melt
head(data.m)
names(data.m)[14]="Sample" # rename the Sample column
names(data.m)
data.m.key <- join(data.m, key_all, by="Sample", type="left") # join metadata
head(data.m.key)
data.m.pit <- subset(data.m.key, Depth=="PIT") # isolate just the PIT samples
# Rename so it matches the microscopy taxonomy
data.m.pit$Plotting<-"Other"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Alveolata-Ciliates"]="Alveolata-Ciliates"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Alveolata-Dinoflagellates"]="Alveolata-Dinoflagellates"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Stramenopila-Diatoms"]="Stramenopila-Diatoms"
data.m.pit$Plotting[data.m.pit$Level4 == "Dictyochophyceae"]="Stramenopila-Silicoflagellates"
# combine acanths and rads to one category
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Rhizaria-Polycystines"]="Rhizaria-Radiolarians"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Rhizaria-Acantharians"]="Rhizaria-Radiolarians"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Rhizaria-Cercozoans"] = "Rhizaria-Other"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Rhizaria-Other"] = "Rhizaria-Other"

# aggregate data for plotting
data.agg.pit<-aggregate(data.m.pit$value, by=list(Taxa=data.m.pit$Plotting,Sample=data.m.pit$Station),sum) #sum sequences by taxonomic group
data.agg.pit$tax.order<-factor(data.agg.pit$Taxa,levels=(tax_full.pit), labels=(tax_full.pit))

### load microscopy data and edit to match asvs
counts<-read.csv("~/github/mesoscope2017/raw-data/CellsPit_Final.csv",header=T)
# Renaming the organisms from my microscopy dataframe to match the molecular names
rownames(counts)=counts$Type
# combine ciliates and tintinnid to the same category
counts$Plotting = counts$Type
counts$Plotting[counts$Type == "Tintinnid"]="Alveolata-Ciliates"
counts$Plotting[counts$Type == "Ciliate"]="Alveolata-Ciliates"
counts$Plotting[counts$Type == "Dinoflagellate"]="Alveolata-Dinoflagellates"
# combine acanths and rads to one category
counts$Plotting[counts$Type == "Rhizaria-Acantharians"]="Rhizaria-Radiolarians"
counts$Plotting[counts$Type == "Rhizaria-other"]="Rhizaria-Other"
counts$Plotting[counts$Type == "Diatom"]="Stramenopila-Diatoms"
counts$Plotting[counts$Type == "Dictyocha"]="Stramenopila-Silicoflagellates"

# melt the dataframe and identify reps
counts.m <- counts %>%
  melt(.,id.vars = c("Type","Plotting")) %>%
  separate(variable, c("Station", "Rep"), sep = "_", remove = TRUE)

#To calculate average values of taxa for plotting by PIT
counts.mean<- counts.m %>% 
  group_by(Plotting, Station) %>%
  summarise(mean=mean(value),stdev=sd(value)) %>%
  as.data.frame
counts.mean

# factoring for plotting
counts.mean$tax.order<-factor(counts.mean$Plotting,levels=(tax_full.pit), labels=(tax_full.pit))

#### combine data for plotting
pit_asvs_plotting <- data.agg.pit 
names(pit_asvs_plotting)= c("Taxa", "Sample", "value","tax.order") # rename the columns
pit_asvs_plotting$sampletype = "ASV" # add a sampletype column and assign asv
head(pit_asvs_plotting)
pit_micro_plotting <- select(counts.mean, -stdev) # remove the standard deviation 
names(pit_micro_plotting)= c("Taxa", "Sample", "value","tax.order") # rename columns to match asvs
pit_micro_plotting$sampletype = "Micro_count" # add a sampletype column and assign microscope
head(pit_micro_plotting)
pit_combo <- bind_rows(pit_asvs_plotting,pit_micro_plotting) # combine data
pit_combo$tax.order<-factor(pit_combo$tax.order,levels=(tax_full.pit), labels=(tax_full.pit))
pit_combo$Sample<-factor(pit_combo$Sample,levels=c("PIT12","PIT11","PIT10","PIT9","PIT8","PIT7","PIT6","PIT5","PIT4","PIT3","PIT2","PIT1"))

## plot
pit_plot<-ggplot(pit_combo, aes(y=value, fill=tax.order, x=Sample))+
  geom_bar(position = position_fill(),stat="identity",color="black",aes(fill=tax.order), width=0.8)+
  scale_fill_manual(values=pit_colors, name="Taxa",labels=tax_full.pit)+
  theme_bw()+
  labs(title="", x="",y="Relative abundance")+
  guides(fill = guide_legend(ncol = 1))+
  theme(axis.ticks = element_blank())+
  scale_x_discrete(labels = c("PIT12" = "P12","PIT11"="P11","PIT10"="P10","PIT9"="P9","PIT8"="P8",
                              "PIT7"="P7","PIT6"="P6","PIT5"="P5","PIT4"="P4","PIT3"="P3","PIT2"="P2","PIT1"="P1"))+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle=90),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
        legend.text = element_text(size = 10))+
  facet_wrap(sampletype~.,ncol=1) +
  theme(
    strip.background = element_blank())
pit_plot
# ggplot2::ggsave("~/github/mesoscope2017/figures/Figure4_pit_barplot.pdf", width=2, height= 3.5, units="in", dpi=300, scale=2)

#### statistics for individual groups ####
# Step 1. need to sum the number of reads per material type and depth and station
percent_micro <- counts.mean %>%
  group_by(Station) %>%
  summarise(total_counts = sum(mean))

# Step 2. sum the number of reads in the group I'm interested
sum_cil_micro <- counts.mean %>%
  group_by(Station) %>%
  subset(., (Plotting == "Alveolata-Ciliates")) %>%
  summarise(total_ciliates = sum(mean))

sum_alv_micro <- counts.mean %>%
  group_by(Station) %>%
  subset(., (Plotting == "Alveolata-Ciliates"| Plotting == "Alveolata-Dinoflagellates" )) %>%
  summarise(total_Alveolata = sum(mean))

sum_rhiz_micro <- counts.mean %>%
  group_by(Station) %>%
  subset(., (Plotting == "Rhizaria-Foraminifera"| Plotting == "Rhizaria-Other" | Plotting =="Rhizaria-Radiolaria")) %>%
  summarise(total_rhiz = sum(mean))

sum_stram_micro <- counts.mean %>%
  group_by(Station) %>%
  subset(., (Plotting == "Stramenopila-Diatoms"| Plotting == "Stramenopila-Dictyocha")) %>%
  summarise(total_stram = sum(mean))

sum_diatom_micro <- counts.mean %>%
  group_by(Station) %>%
  subset(., (Plotting == "Stramenopila-Diatoms")) %>%
  summarise(total_diatom = sum(mean))

# Step 3. divide the sum of the group reads by the total numbers of reads
percent_all_micro <- data.frame(percent_micro, sum_alv_micro$total_Alveolata, sum_rhiz_micro$total_rhiz,
                                sum_stram_micro$total_stram, sum_diatom_micro$total_diatom, sum_cil_micro$total_ciliates) %>%
  mutate(percent_alv=sum_alv_micro.total_Alveolata/total_counts,
         percent_rhiz= sum_rhiz_micro.total_rhiz/total_counts,
         percent_stram= sum_stram_micro.total_stram/total_counts,
         percent_diatom= sum_diatom_micro.total_diatom/total_counts,
         percent_ciliate = sum_cil_micro.total_ciliates/total_counts)

percent_pit_micro <- percent_all_micro %>%
  summarise(mean_alv = mean(percent_alv),
            stdev_alv = sd(percent_alv),
            mean_rhiz = mean(percent_rhiz),
            stdev_rhiz = sd(percent_rhiz),
            mean_stram = mean(percent_stram),
            stdev_stram = sd(percent_stram),
            mean_diatom = mean(percent_diatom),
            stdev_diatom = sd(percent_diatom),
            mean_cil = mean(percent_ciliate),
            stdev_cil = sd(percent_ciliate))


percent_pit_cylclonic <- filter(percent_all_micro, (Station=="PIT1" | Station =="PIT2")) %>%
  summarise(mean_alv = mean(percent_alv),
            stdev_alv = sd(percent_alv),
            mean_rhiz = mean(percent_rhiz),
            stdev_rhiz = sd(percent_rhiz),
            mean_stram = mean(percent_stram),
            stdev_stram = sd(percent_stram),
            mean_diatom = mean(percent_diatom),
            stdev_diatom = sd(percent_diatom),
            mean_cil = mean(percent_ciliate),
            stdev_cil = sd(percent_ciliate))

percent_pit_anticylclonic <- filter(percent_all_micro, (Station=="PIT12" | Station =="PIT11" | Station=="PIT10" | Station =="PIT9" | Station=="PIT8" | Station =="PIT7" | Station=="PIT6")) %>%
  summarise(mean_alv = mean(percent_alv),
            stdev_alv = sd(percent_alv),
            mean_rhiz = mean(percent_rhiz),
            stdev_rhiz = sd(percent_rhiz),
            mean_stram = mean(percent_stram),
            stdev_stram = sd(percent_stram),
            mean_diatom = mean(percent_diatom),
            stdev_diatom = sd(percent_diatom),
            mean_cil = mean(percent_ciliate),
            stdev_cil = sd(percent_ciliate))
