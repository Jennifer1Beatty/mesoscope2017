## Analysis of the MESO-SCOPE 2017 Cruise 18S data, final figures
## Author: Jennifer Beatty, jlbeatty@usc.edu


# Load necessary packages

library(plyr)
library(tidyverse)
library(reshape2)
# library(decontam)  # for removing possible contaminated ASVs, doesn't work with updated R
library(vegan)  # for statistical analysis
library(RColorBrewer)
#library(edgeR)  # for statistical normalization, doesn't work with udpated R
library(ggpubr)
library(scico) # for scientifically approved colors
library(ggupset)  # for upset plots
library(stringi)  # for maniuplating strings
library(VennDiagram)
library(ggvenn)
library(patchwork)
library(rstatix)
# for full size images: 5" w and 6" height, 300 dpi
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Final_Figs/Taxa.pdf", width=5, height = 6, units="in", dpi=300)

# for column size images: 3.5" w and, 300 dpi
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Final_Figs/Taxa.pdf", width=3.5, height = 6, units="in", dpi=300)

#### Input all the extra variables and factor them #####
key_all<-read.delim("~/Desktop/Chapter1/MESOSCOPE2017/Data/variables_all.txt")
key_all$Station<-factor(as.character(key_all$Station),levels=c("S4","S6","S8","S10","S12","S14",
                                                               "PIT12","PIT11","PIT10","PIT9","PIT8","PIT7","PIT6","PIT5","PIT4","PIT3","PIT2","PIT1"))
key_all$Station.sla<-factor(as.character(key_all$Station), levels=rev(c("S6", "S8", "S14","S10","S4","S12",
                                                                        "PIT12","PIT11","PIT10","PIT9","PIT8","PIT7","PIT6","PIT5","PIT4","PIT3","PIT2","PIT1")))
key_all$Color.sla <- factor(as.character(key_all$Color.sla), levels = c("#7ea6c3","#de2d26","#d77773","#cfcdcd","#2c7fb8", "#d2acab"))
key_all$Depth<-factor(as.character(key_all$Depth),levels=c("PIT","15m","DCM","175m","500m"))
key_all$Material<-factor(key_all$Material,levels=c("RNA","DNA","PIT"))
key_all$Color.depth<-factor(as.character(key_all$Color.depth),levels=c("#d55e00","#f0e442","#009e73","#cc79a7","#0072b2"))
key_all$Shape.type<-factor(as.numeric(key_all$Shape.type),levels=c(1,5))
key_all$Shape.mat<-as.numeric(key_all$Shape.mat, levels=c(22,25,21))
key_all$Line.mat<-factor(as.numeric(key_all$Line.mat),levels=c(6, 1, 3))
key_all$Station.density <- paste(key_all$Station, key_all$Potential_Density, sep="_")
key_all$Stat.depth.density <- paste(key_all$Station.density, key_all$Depth, sep="_")
key_all$Sample.type <- factor(as.character(key_all$Sample.type), levels=c( "Water column RNA","Water column DNA", "Trap DNA"))
summary(key_all)
# plotting variables
tax_order_nometa<-c("Alveolates-Ciliates","Alveolates-Dinophyceae","Alveolates-Syndiniales", "Archaeplastids-Chlorophytes",
             "Excavates-Discoba", "Hacrobia-Cryptophytes","Hacrobia-Haptophytes",
             "Opisthokont-Choanoflagellida","Other/unknown",
             "Rhizaria-Acantharia","Rhizaria-Cercozoa", "Rhizaria-Other","Rhizaria-Polycystines",
             "Stramenopiles-Chrysophytes","Stramenopiles-Diatoms","Stramenopiles-MAST","Stramenopiles-Ochrophyta","Stramenopiles-Other","Stramenopiles-Pelagophytes")
tax_color_phyla2_nometa<-c('firebrick4','indianred1','tomato3','forestgreen','yellowgreen','darkblue',
                           'lightblue','moccasin','grey','mediumvioletred',"#C291A4",'#db6ea0','lightpink',"#DECDBE",'#DDAD4B','tan2','tan3','tan4',"#5C4033")
tax_color_2 <- c("Alveolates-Ciliates"='firebrick4',"Alveolates-Dinophyceae"='indianred1',"Alveolates-Syndiniales"='tomato3', "Archaeplastids-Chlorophytes"='forestgreen',
                 "Excavates-Discoba"='yellowgreen', "Hacrobia-Cryptophytes"='darkblue',"Hacrobia-Haptophytes"='lightblue',
                 "Opisthokont-Choanoflagellida"='moccasin',"Other/unknown"='grey',
                 "Rhizaria-Acantharia"='mediumvioletred',"Rhizaria-Cercozoa"="#C291A4","Rhizaria-Other"='#db6ea0',"Rhizaria-Polycystines"='lightpink', 
                 "Stramenopiles-Chrysophytes"="#DECDBE","Stramenopiles-Diatoms"='#DDAD4B',"Stramenopiles-MAST"='tan2',"Stramenopiles-Ochrophyta"='tan3',"Stramenopiles-Other"='tan4',"Stramenopiles-Pelagophytes"="#5C4033")
WC <- c("RNA","DNA")
col_depth <-c("15m" = "#fee08b", "DCM" = "#7fbc41", "175m" = "#74add1", "500m" = "#542788")
col_material<- c("#7fcdbb","#2c7fb8","#fc8d59")
col_material_2<- c("Water column RNA"="#7fcdbb","Water column DNA"="#2c7fb8","Trap DNA"="#fc8d59")

## START HERE FOR ASV TABLE WITH UPDATED TAXONOMY
asv_df <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Output/MS_Clean_No1_Means_Complete_Norm_NewTax_1122.csv", header=T)
asv_df[1]=NULL
row.names(asv_df) = asv_df$ASV1.ID
nometa <-filter(asv_df, Taxa!="Opisthokont-Metazoa")
####### Figure 1. Map ########

####### Figure 2. Environmental features #######
####### Figure 3. WC Taxa Barplots #####
#### Taxa barplots for water column
data.nometa.m<-melt(nometa) #melt
head(data.nometa.m)
data.nometa.agg<-aggregate(data.nometa.m$value, by=list(Taxa=data.nometa.m$TaxaPlot,Sample=data.nometa.m$variable),sum) #sum sequences by taxonomic group
data.nometa.agg$Sample<-factor(data.nometa.agg$Sample,levels=names(nometa[2:61]))
nometa_bar_plot_df <- join(data.nometa.agg, key_all, by="Sample", type="left", match="first") %>%
  arrange(.,Station,Material,Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))
#Bar plot of community composition ordered by Material, Station and Depth
bars_nometa<-ggplot(nometa_bar_plot_df, aes(y=x, fill=tax, x=Station))+
  geom_bar(position = position_fill(reverse=TRUE),stat="identity",color="black",aes(fill=Taxa))+
  scale_fill_manual(values=tax_color_phyla2_nometa, name="Taxa")+
  theme_bw()+
  scale_x_discrete(limits=rev(factor(nometa_bar_plot_df$Station[1:912]))) +
  labs(title="", x="",y="")+
  theme(legend.position="right",
        legend.text = element_text(size = 8),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 10, angle = 45, hjust = .9, vjust = .9),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5, vjust = .5),
        legend.key.size = unit(0.5,'cm'))
bars_2 <-bars_nometa %+% subset(nometa_bar_plot_df, Material %in% WC) + coord_flip() + facet_grid(Depth~Material) + ylab("Relative abundance")
# tag_facet_outside(bars_2, open = c("",""), close = c("","."))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Final_Figs/Taxa_WaterColumn.pdf", width=8, height = 6, units="in", dpi=300)
# ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/Taxa_All_NOMeta_0323WC_DepthSLA.pdf", width=8, height = 6)

####### Fig 3, Stats/percentages of specific taxa #######
# Step 1. need to sum the number of reads per material type and depth and station
percent <- nometa_bar_plot_df %>%
  group_by(Material, Depth, Station) %>%
  summarise(total_reads = sum(x))

# Step 2. sum the number of reads in the group I'm interested
sum_alv <- nometa_bar_plot_df %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Alveolates-Ciliates"| Taxa == "Alveolates-Dinophyceae" | Taxa == 
               "Alveolates-Syndiniales")) %>%
  summarise(total_alveolates = sum(x))

sum_stram <- nometa_bar_plot_df %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Stramenopiles-Chrysophytes" | Taxa == "Stramenopiles-Diatoms"|Taxa == "Stramenopiles-MAST"|
               Taxa == "Stramenopiles-Ochrophyta" | Taxa == "Stramenopiles-Other" | Taxa == "Stramenopiles-Pelagophytes")) %>%
  summarise(total_stramenopiles = sum(x))

sum_chloro <- nometa_bar_plot_df %>%
  group_by(Material, Depth, Station) %>%
  subset(., Taxa == "Archaeplastids-Chlorophytes") %>%
  summarise(total_chlorophytes = sum(x))

sum_pelago<- nometa_bar_plot_df %>%
  group_by(Material, Depth, Station) %>%
  subset(., Taxa == "Stramenopiles-Pelagophytes") %>%
  summarise(total_pelago = sum(x))

sum_rhizaria<- nometa_bar_plot_df %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Rhizaria-Acantharia"| Taxa == "Rhizaria-Cercozoa" | Taxa == "Rhizaria-Other" | Taxa == "Rhizaria-Polycystines")) %>%
  summarise(total_rhizaria = sum(x))

sum_hapto <- nometa_bar_plot_df %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Hacrobia-Haptophytes")) %>%
  summarise(total_hapto = sum(x))

sum_acanth<- nometa_bar_plot_df %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Rhizaria-Acantharia")) %>%
  summarise(total_acanth = sum(x))

sum_poly<- nometa_bar_plot_df %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Rhizaria-Polycystines")) %>%
  summarise(total_poly = sum(x))

sum_rhiz_other<- nometa_bar_plot_df %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Rhizaria-Other")) %>%
  summarise(total_rhiz_other = sum(x))

# alveolates, chlorophytes, stramenopiles, rhizaria, pelagophytes
# want to talk about how the composition of alveolates was difference between the DNA and RNA
total_ciliates <- nometa_bar_plot_df %>%
  group_by(Material,Depth,Station) %>%
  subset(., Taxa=="Alveolates-Ciliates") %>%
  summarise(total_ciliates = sum(x))
total_dino <- nometa_bar_plot_df %>%
  group_by(Material,Depth,Station) %>%
  subset(., Taxa=="Alveolates-Dinophyceae") %>%
  summarise(total_dinos = sum(x))
total_synd <- nometa_bar_plot_df %>%
  group_by(Material,Depth,Station) %>%
  subset(., Taxa=="Alveolates-Syndiniales") %>%
  summarise(total_synd = sum(x))


# Step 3. divide the sum of the group reads by the total numbers of reads
percent_all <- data.frame(percent, sum_alv$total_alveolates, sum_chloro$total_chlorophytes, sum_pelago$total_pelago, sum_rhizaria$total_rhizaria, sum_stram$total_stramenopiles,
                          total_ciliates$total_ciliates, total_dino$total_dinos, total_synd$total_synd,
                          sum_hapto$total_hapto, sum_acanth$total_acanth, sum_poly$total_poly, sum_rhiz_other$total_rhiz_other) %>%
  mutate(percent_alv=sum_alv.total_alveolates/total_reads,
         percent_chloro= sum_chloro.total_chlorophytes/total_reads,
         percent_pelago= sum_pelago.total_pelago/total_reads,
         percent_rhizaria= sum_rhizaria.total_rhizaria/total_reads,
         percent_stram= sum_stram.total_stramenopiles/total_reads,
         percent_ciliates = total_ciliates.total_ciliates/total_reads,
         percent_dino = total_dino.total_dinos/total_reads,
         percent_synd = total_synd.total_synd/total_reads,
         percent_hapto = sum_hapto.total_hapto/total_reads,
         percent_acanth = sum_acanth.total_acanth/total_reads,
         percent_poly = sum_poly.total_poly/total_reads,
         percent_rhiz_other = sum_rhiz_other.total_rhiz_other/total_reads)

percent_wc_stats <- filter(percent_all, (Material=="DNA" | Material=="RNA")) %>%
  summarise(mean_alv = mean(percent_alv),
            stdev_alv = sd(percent_alv),
            mean_rhiz = mean(percent_rhizaria),
            stdev_rhiz = sd(percent_rhizaria))          

percent_dnanodepth_stats <- filter(percent_all, (Material=="DNA")) %>%
  summarise(mean_ciliates = mean(percent_ciliates),
            stdev_ciliates = sd(percent_ciliates),
            mean_dino = mean(percent_dino),
            stdev_dino = sd(percent_dino),
            mean_synd = mean(percent_synd),
            stdev_synd = sd(percent_synd),
            mean_hapto = mean(percent_hapto),
            stdev_hapto = sd(percent_hapto)) 

percent_rnanodepth_stats <- filter(percent_all, (Material=="RNA")) %>%
  summarise(mean_ciliates = mean(percent_ciliates),
            stdev_ciliates = sd(percent_ciliates),
            mean_dino = mean(percent_dino),
            stdev_dino = sd(percent_dino),
            mean_synd = mean(percent_synd),
            stdev_synd = sd(percent_synd),
            mean_stram = mean(percent_stram),
            stdev_stram = sd(percent_stram),
            mean_hapto = mean(percent_hapto),
            stdev_hapto = sd(percent_hapto))

percent_dna_stats <- filter(percent_all, (Material=="DNA")) %>%
  group_by(Depth) %>%
  summarise(mean_chloro = mean(percent_chloro),
            stdev_chloro = sd(percent_chloro),
            mean_stram = mean(percent_stram),
            stdev_stram = sd(percent_stram),
            mean_rhiz = mean(percent_rhizaria),
            stdev_rhiz = sd(percent_rhizaria),
            mean_ciliates = mean(percent_ciliates),
            stdev_ciliates = sd(percent_ciliates),
            mean_dino = mean(percent_dino),
            stdev_dino = sd(percent_dino),
            mean_synd = mean(percent_synd),
            stdev_synd = sd(percent_synd),
            mean_hapto = mean(percent_hapto),
            stdev_hapto = sd(percent_hapto),
            mean_acanth = mean(percent_acanth),
            stdev_acanth = sd(percent_acanth),
            mean_rhizother = mean(percent_rhiz_other),
            stdev_rhizother = sd(percent_rhiz_other))

percent_rna_stats <- filter(percent_all, (Material=="RNA")) %>%
  group_by(Depth) %>%
  summarise(mean_pelago = mean(percent_pelago),
            stdev_pelago = sd(percent_pelago),
            mean_rhiz = mean(percent_rhizaria),
            stdev_rhiz = sd(percent_rhizaria),
            mean_ciliates = mean(percent_ciliates),
            stdev_ciliates = sd(percent_ciliates),
            mean_dino = mean(percent_dino),
            stdev_dino = sd(percent_dino),
            mean_synd = mean(percent_synd),
            stdev_synd = sd(percent_synd),
            mean_hapto = mean(percent_hapto),
            stdev_hapto = sd(percent_hapto),
            mean_stram = mean(percent_stram),
            stdev_stram = sd(percent_stram),
            mean_acanth = mean(percent_acanth),
            stdev_acanth = sd(percent_acanth),
            mean_rhizother = mean(percent_rhiz_other),
            stdev_rhizother = sd(percent_rhiz_other))


percent_pit_stats <- filter(percent_all, (Material=="PIT")) %>%
  summarise(mean_stram = mean(percent_stram),
            stdev_stram = sd(percent_stram),
            mean_rhiz = mean(percent_rhizaria),
            stdev_rhiz = sd(percent_rhizaria),
            mean_alv = mean(percent_alv),
            stdev_alv = sd(percent_alv))

####### Figure 4. Trap barplots ######
# Microscopy counts
counts<-read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/CellsPit_Final.csv",header=T)
# Renaming the organisms from my microscopy dataframe to match the molecular names
rownames(counts)=counts$Type
# combine ciliates and tintinnid to the same category
counts$Plotting = counts$Type
counts$Plotting[counts$Type == "Tintinnid"]="Alveolates-Ciliates"
counts$Plotting[counts$Type == "Ciliate"]="Alveolates-Ciliates"
counts$Plotting[counts$Type == "Dinoflagellate"]="Alveolates-Dinoflagellates"
# combine acanths and rads to one category
counts$Plotting[counts$Type == "Rhizaria-Acantharia"]="Rhizaria-Radiolaria"
counts$Plotting[counts$Type == "Rhizaria-other"]="Rhizaria-Other"
counts$Plotting[counts$Type == "Diatom"]="Stramenopiles-Diatoms"
counts$Plotting[counts$Type == "Dictyocha"]="Stramenopiles-Dictyocha"
# make the color match the other type

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
tax_order.counts<-c("Alveolates-Ciliates","Alveolates-Dinoflagellates","Rhizaria-Foraminifera","Rhizaria-Other","Rhizaria-Radiolaria", "Stramenopiles-Diatoms","Stramenopiles-Dictyocha")
col.pit.counts <-c('firebrick4','indianred1','#fde0dd','#db6ea0','lightpink','#DDAD4B','gold1')
counts.mean$tax.order<-factor(counts.mean$Plotting,levels=(tax_order.counts), labels=(tax_order.counts))
counts.mean$Station <- factor(counts.mean$Station, levels=c("PIT12","PIT11","PIT10","PIT9","PIT8","PIT7","PIT6","PIT5","PIT4","PIT3","PIT2","PIT1"))

microscopy_barplot<-ggplot(counts.mean, aes(y=mean, fill=tax.order, x=Station))+
  geom_bar(position = position_fill(reverse=TRUE),stat="identity",color="black",aes(fill=tax.order))+
  scale_fill_manual(values=col.pit.counts,name="Taxa")+
  theme_bw()+
  labs(title="", x="",y="Relative abundance of trap microscopy counts")+
  theme(legend.position="right",
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 10, angle = 45, hjust = .9, vjust = .9),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5))
pit_micro_barplot <-microscopy_barplot + coord_flip() + scale_x_discrete(limits=rev)
# ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/MicroscopyCellMean_PIT_relav.pdf", width=8, height = 6)

#Now to restructure the ASV PIT data to look like the microscopy data
#no meta sequences
data.m<-melt(nometa) #melt
head(data.m)
names(data.m)[14]="Sample"
names(data.m)
data.m.key <- join(data.m, key_all, by="Sample", type="left")
head(data.m.key)
data.m.pit <- subset(data.m.key, Depth=="PIT")
# Rename so it matches the microscopy taxonomy
data.m.pit$Plotting<-"Other"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Alveolates-Ciliates"]="Alveolates-Ciliates"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Alveolates-Dinophyceae"]="Alveolates-Dinoflagellates"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Stramenopiles-Diatoms"]="Stramenopiles-Diatoms"
data.m.pit$Plotting[data.m.pit$Level4 == "Dictyochophyceae"]="Stramenopiles-Dictyocha"
# combine acanths and rads to one category
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Rhizaria-Polycystines"]="Rhizaria-Radiolaria"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Rhizaria-Acantharia"]="Rhizaria-Radiolaria"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Rhizaria-Cercozoa"] = "Rhizaria-Other"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Rhizaria-Other"] = "Rhizaria-Other"

tax_order.pit<-c("Alveolates-Ciliates","Alveolates-Dinoflagellates","Other","Rhizaria-Other","Rhizaria-Radiolaria", "Stramenopiles-Diatoms","Stramenopiles-Dictyocha")
col.pit <-c('firebrick4','indianred1', 'grey','#db6ea0','lightpink','#DDAD4B','gold1')

data.agg.pit<-aggregate(data.m.pit$value, by=list(Taxa=data.m.pit$Plotting,Sample=data.m.pit$Station),sum) #sum sequences by taxonomic group
data.agg.pit$tax.order<-factor(data.agg.pit$Taxa,levels=(tax_order.pit), labels=(tax_order.pit))

#Bar plot of community composition ordered by Material, Station and Depth
data.agg.pit$Sample<-factor(data.agg.pit$Sample,levels=c("PIT12","PIT11","PIT10","PIT9","PIT8","PIT7","PIT6","PIT5","PIT4","PIT3","PIT2","PIT1"))

asv_barplot<-ggplot(data.agg.pit, aes(y=x, fill=tax.order, x=Sample))+
  geom_bar(position = position_fill(reverse=TRUE),stat="identity",color="black",aes(fill=tax.order))+
  scale_fill_manual(values=col.pit, name="Taxa")+
  theme_bw()+
  labs(title="", x="",y="Relative abundance of trap DNA reads")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 10, angle = 45, hjust = .9, vjust = .9),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5))
asv_plot <-asv_barplot + coord_flip() + scale_x_discrete(limits=rev) # +  guides(fill = guide_legend(override.aes=list(color=col_full.pit, labels=tax_full.pit)))
# ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/ASVPIT_relav.pdf", width=8, height = 6)
# elements for the combined plot
tax_full.pit<-c("Alveolates-Ciliates","Alveolates-Dinoflagellates","Other","Rhizaria-Foraminifera","Rhizaria-Other","Rhizaria-Radiolaria","Stramenopiles-Diatoms","Stramenopiles-Dictyocha")
col_full.pit <-c('firebrick4','indianred1', 'grey','#fde0dd','lightpink','#db6ea0','#DDAD4B',"gold1")

trial_colors <- c("Alveolates-Ciliates"='firebrick4',"Alveolates-Dinoflagellates"='indianred1',"Other"= 'grey',
                  "Rhizaria-Foraminifera"='mediumvioletred',"Rhizaria-Radiolaria"='lightpink',"Rhizaria-Other"='#C291A4',
                  "Stramenopiles-Diatom"='#DDAD4B',"Stramenopiles-Dictyocha"="#5C4033")
asv_plot + pit_micro_barplot + plot_layout(guides="collect") 
# ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/PIT_relav_combo_newcolor.pdf", width=10, height = 6)

# try to combine the dataframes to use facet wrap instead.
pit_asvs_plotting <- data.agg.pit 
names(pit_asvs_plotting)= c("Taxa", "Sample", "value","tax.order")
pit_asvs_plotting$sampltype = "ASV"
head(pit_asvs_plotting)
pit_micro_plotting <- select(counts.mean, -stdev)
names(pit_micro_plotting)= c("Taxa", "Sample", "value","tax.order")
pit_micro_plotting$sampltype = "Micro_count"
head(pit_micro_plotting)
pit_combo <- bind_rows(pit_asvs_plotting,pit_micro_plotting)
tax_full.pit<-c("Alveolates-Ciliates","Alveolates-Dinoflagellates","Other","Rhizaria-Foraminifera","Rhizaria-Other","Rhizaria-Radiolaria","Stramenopiles-Diatoms","Stramenopiles-Dictyocha")
col_full.pit <-c('firebrick4','indianred1', 'grey','#fde0dd','lightpink','#db6ea0','#DDAD4B',"gold1")
pit_combo$tax.order<-factor(pit_combo$tax.order,levels=(tax_full.pit), labels=(tax_full.pit))

pit_plot<-ggplot(pit_combo, aes(y=value, fill=tax.order, x=Sample))+
  geom_bar(position = position_fill(),stat="identity",color="black",aes(fill=tax.order))+
  scale_fill_manual(values=col_full.pit, name="Taxa", labels=tax_full.pit)+
  theme_bw()+
  labs(title="", x="",y="Relative abundance")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.text = element_text(size = 8),
        axis.title.x = element_text(size = 10, angle = 0, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 8, angle = 90, hjust = .9, vjust = .9),
        axis.text.y = element_text(size = 8, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 10, angle = 90, hjust = .5, vjust = .5),
        legend.key.size = unit(0.25,'cm'))+
  facet_wrap(sampltype~.,ncol=1) 
pit_plot
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Final_Figs/Taxa_pits.pdf", width=3.5, height = 5, units="in", dpi=300)

####### Figure 4, stats/percentages #####
# Step 1. need to sum the number of reads per material type and depth and station
percent_micro <- counts.mean %>%
  group_by(Station) %>%
  summarise(total_counts = sum(mean))

# Step 2. sum the number of reads in the group I'm interested
sum_alv_micro <- counts.mean %>%
  group_by(Station) %>%
  subset(., (Plotting == "Alveolates-Ciliates"| Plotting == "Alveolates-Dinophyceae" )) %>%
  summarise(total_alveolates = sum(mean))

sum_rhiz_micro <- counts.mean %>%
  group_by(Station) %>%
  subset(., (Plotting == "Rhizaria-Foraminifera"| Plotting == "Rhizaria-Other" | Plotting =="Rhizaria-Radiolaria")) %>%
  summarise(total_rhiz = sum(mean))

sum_stram_micro <- counts.mean %>%
  group_by(Station) %>%
  subset(., (Plotting == "Stramenopiles-Diatoms"| Plotting == "Stramenopiles-Dictyocha")) %>%
  summarise(total_stram = sum(mean))

sum_diatom_micro <- counts.mean %>%
  group_by(Station) %>%
  subset(., (Plotting == "Stramenopiles-Diatoms")) %>%
  summarise(total_diatom = sum(mean))

# Step 3. divide the sum of the group reads by the total numbers of reads
percent_all_micro <- data.frame(percent_micro, sum_alv_micro$total_alveolates, sum_rhiz_micro$total_rhiz,
                                sum_stram_micro$total_stram, sum_diatom_micro$total_diatom) %>%
  mutate(percent_alv=sum_alv_micro.total_alveolates/total_counts,
         percent_rhiz= sum_rhiz_micro.total_rhiz/total_counts,
         percent_stram= sum_stram_micro.total_stram/total_counts,
         percent_diatom= sum_diatom_micro.total_diatom/total_counts)

percent_pit_micro <- percent_all_micro %>%
  summarise(mean_alv = mean(percent_alv),
            stdev_alv = sd(percent_alv),
            mean_rhiz = mean(percent_rhiz),
            stdev_rhiz = sd(percent_rhiz),
            mean_stram = mean(percent_stram),
            stdev_stram = sd(percent_stram),
            mean_diatom = mean(percent_diatom),
            stdev_diatom = sd(percent_diatom))


percent_pit_cylclonic <- filter(percent_all_micro, (Station=="PIT1" | Station =="PIT2")) %>%
  summarise(mean_alv = mean(percent_alv),
            stdev_alv = sd(percent_alv),
            mean_rhiz = mean(percent_rhiz),
            stdev_rhiz = sd(percent_rhiz),
            mean_stram = mean(percent_stram),
            stdev_stram = sd(percent_stram),
            mean_diatom = mean(percent_diatom),
            stdev_diatom = sd(percent_diatom))

percent_pit_anticylclonic <- filter(percent_all_micro, (Station=="PIT12" | Station =="PIT11" | Station=="PIT10" | Station =="PIT9" | Station=="PIT8" | Station =="PIT7" | Station=="PIT6")) %>%
  summarise(mean_alv = mean(percent_alv),
            stdev_alv = sd(percent_alv),
            mean_rhiz = mean(percent_rhiz),
            stdev_rhiz = sd(percent_rhiz),
            mean_stram = mean(percent_stram),
            stdev_stram = sd(percent_stram),
            mean_diatom = mean(percent_diatom),
            stdev_diatom = sd(percent_diatom))

####### Figure 5.  NMDS ASV  ######

nmdsPoints <- function(df,key_df){
  new_df<-dcast(df, variable~ASV.ID, fill=0) #restructure dataframe to only include ASV1.ID and samples
  row.names(new_df)<-new_df$variable; new_df$variable<-NULL
  NMDS_bray=metaMDS(new_df,distance="bray",k=2,trymax=100,engine=c("monoMDS"),autotransform=FALSE)
  stressplot(NMDS_bray)
  bray_pts <- NMDS_bray$points[1:nrow(NMDS_bray$points),]
  bray_pts<-as.data.frame(bray_pts) #makes nmds points into a dataframe
  plot(bray_pts$MDS1, bray_pts$MDS2) #quick view
  bray_pts$Sample<-row.names(bray_pts) #makes samples row names
  bray_df<-join(bray_pts, key_df, by="Sample", type="left", match="first")
  return(bray_df)
}
Norm.rel_nometa <-decostand(nometa[,2:61], MARGIN=2, method = "total") #this does it by rows
colSums(Norm.rel_nometa) # check! should all equal 1.
Norm.rel_nometa$ASV.ID<-rownames(Norm.rel_nometa)
melt_norm_nometa<- melt(Norm.rel_nometa)
vars_nometa<-colsplit(melt_norm_nometa$variable, "_", c("Station","Depth","Material"))
new_df_nometa<-data.frame(melt_norm_nometa, vars_nometa)

all_nometa<- nmdsPoints(new_df_nometa, key_all)
nmds_all_nometa<-ggplot(data=all_nometa, aes(x = MDS1, y = MDS2, fill=SLA,shape=Sample.type)) + 
  geom_point(color="black",size=3,position="jitter") + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = c(22,25,21))+ 
  # scale_fill_manual(values=all_nometa$Col.sla)+
  labs(title="", x="NMDS1",y="NMDS2",size=10)+
  scale_y_continuous(breaks = seq(-3, 3, by = 1))+
  scale_x_continuous(breaks = seq(-4, 2, by = 1))+
  scale_fill_gradient2(low = "#0818A8",
                       mid = "light gray",
                       high = "#880808",
                       midpoint = 0,
                       limits=c(-28,28))
# guides(fill = guide_legend(override.aes=list(values=all_nometa$Col.sla)))
nmds_all <- nmds_all_nometa +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Sample.type),linewidth=0.5) + scale_color_manual(values=col_material_2) #this adds in 95% confidence interval ellipses from ggplot
# ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_all_nometa_sla_depth.pdf", width=6, height = 6)

dna_nometa <- filter(new_df_nometa, Material=="DNA")
dna_df_nometa<-nmdsPoints(dna_nometa, key_all)
dna_nmds_nometa<-ggplot(data=dna_df_nometa, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
  geom_point(color="black",size=3,position="jitter") + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = unique(dna_df_nometa$Shape.depth))+ 
  # scale_fill_manual(values=all_nometa$Col.sla)+
  labs(title="", x="NMDS1",y="NMDS2",size=10)+
  scale_y_continuous(breaks = seq(-3, 3, by = 1))+
  scale_x_continuous(breaks = seq(-4, 2, by = 1))+
  scale_fill_gradient2(low = "#0818A8",
                       mid = "light gray",
                       high = "#880808",
                       midpoint = 0,
                       limits=c(-28,28))
# guides(fill = guide_legend(override.aes=list(values=all_nometa$Col.sla)))
dna_nmds <- dna_nmds_nometa +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth),linewidth=0.5, linetype="dashed") + scale_color_manual(values=col_depth) #this adds in 95% confidence interval ellipses from ggplot
dna_nmds

####### See if the DNA NMDS changes without alveolates ########
names(dna_nometa)[1]="ASV1.ID"
dna_nometa_int<-join(dna_nometa, select(nometa, c("ASV1.ID","TaxaPlot","Level2")), by="ASV1.ID", type="left")
dna_nometa_noalv<-filter(dna_nometa_int, Level2!="Alveolata")
names(dna_nometa_noalv)[1]="ASV.ID"
dna_df_nometa_noalv<-nmdsPoints(dna_nometa_noalv, key_all)
dna_nmds_nometa_noalv<-ggplot(data=dna_df_nometa_noalv, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
  geom_point(color="black",size=3,position="jitter") + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = unique(dna_df_nometa$Shape.depth))+ 
  # scale_fill_manual(values=all_nometa$Col.sla)+
  labs(title="DNA-no alveolates", x="NMDS1",y="NMDS2",size=10)+
  scale_y_continuous(breaks = seq(-3, 3, by = 1))+
  scale_x_continuous(breaks = seq(-4, 2, by = 1))+
  scale_fill_gradient2(low = "#0818A8",
                       mid = "light gray",
                       high = "#880808",
                       midpoint = 0,
                       limits=c(-28,28))
# guides(fill = guide_legend(override.aes=list(values=all_nometa$Col.sla)))
dna_nmds_noalv <- dna_nmds_nometa_noalv +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth),linewidth=0.5, linetype="dashed") + scale_color_manual(values=col_depth) #this adds in 95% confidence interval ellipses from ggplot
dna_nmds_noalv
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Final_Figs/nmds_dna_noalv.pdf", width=3.5, height = 3.5, units="in", dpi=300)

dna_nmds + dna_nmds_noalv + plot_layout(guides="collect")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Final_Figs/nmds_dna_and_noalv.pdf", width=6, height = 5, units="in", dpi=300)

########
rna_nometa <- filter(new_df_nometa, Material=="cDNA")
rna_df_nometa<-nmdsPoints(rna_nometa, key_all)
rna_nmds_nometa<-ggplot(data=rna_df_nometa, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
  geom_point(color="black",size=3,position="jitter") + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = unique(rna_df_nometa$Shape.depth))+ 
  # scale_fill_manual(values=all_nometa$Col.sla)+
  labs(title="", x="NMDS1",y="NMDS2",size=10)+
  scale_y_continuous(breaks = seq(-3, 3, by = 1))+
  scale_x_continuous(breaks = seq(-4, 2, by = 1))+
  scale_fill_gradient2(low = "#0818A8",
                       mid = "light gray",
                       high = "#880808",
                       midpoint = 0,
                       limits=c(-28,28))
# guides(fill = guide_legend(override.aes=list(values=all_nometa$Col.sla)))
rna_nmds <- rna_nmds_nometa +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth),linewidth=0.5, linetype="dashed") + scale_color_manual(values=col_depth) #this adds in 95% confidence interval ellipses from ggplot

(nmds_all / rna_nmds / dna_nmds) + plot_layout(guides="collect") +  plot_annotation(tag_levels = "a")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Final_Figs/nmds_all.pdf", width=5, height = 8, units="in", dpi=300)


####### Figure 6. Venn diagram & UpsetR plots ######
#### Venn diagram of ASVs found in different material type ####
venn_data <- list(
  "Water column DNA" = row.names(binary_wide_material_nometa[binary_wide_material_nometa$`Water column DNA`>0,]),
  "Water column RNA" = row.names(binary_wide_material_nometa[binary_wide_material_nometa$`Water column RNA`>0,]),
  "Trap DNA" = row.names(binary_wide_material_nometa[binary_wide_material_nometa$`Trap DNA`>0,])
)
venn<-ggvenn(data=venn_data, 
       columns=c("Water column RNA","Water column DNA","Trap DNA"),
       fill_color = col_material,
       #stroke_color = col_material,
       stroke_alpha=0.5,
       stroke_size = 0.5,
       set_name_size = 4,
       show_percentage = TRUE
)

#### UpsetR ####
melt_ASV_nometa<-melt(nometa)
names(melt_ASV_nometa)[14]="Sample"
melt_all_nometa<-join(melt_ASV_nometa,key_all,by="Sample", type="left", match="first")
head(melt_all_nometa)
toUpset_nometa<-melt_all_nometa[c(1:2,14:18,33)] # We want the ASV ID, Sample, TaxaPlot, value, Station, Depth, Material
head(toUpset_nometa) # Use below to add annotation
toUpset_nometa$Uniq<-paste(toUpset_nometa$ASV1.ID, toUpset_nometa$TaxaPlot, sep="_")

# Change to binary
summed_material_nometa<-aggregate(toUpset_nometa$value, by=list(Sample_type=toUpset_nometa$Sample.type, Uniq=toUpset_nometa$Uniq),sum)
summed_material_nometa$bin<-ifelse(summed_material_nometa$x > 0, 1, 0)
binary_wide_material_nometa<-dcast(summed_material_nometa[c(1,2,4)], Uniq~Sample_type, fill=0)
row.names(binary_wide_material_nometa)<-binary_wide_material_nometa$Uniq; binary_wide_material_nometa$Uniq<-NULL
head(binary_wide_material_nometa)

# Let's make "Intersect" column to tell us which months a given ASV shows up in 
binary_wide_material_nometa$Intersect <- apply(binary_wide_material_nometa > 0, 1, function(x){toString(names(binary_wide_material_nometa)[x])})
head(binary_wide_material_nometa)
binary_wide_nometa_forotherplots <- binary_wide_material_nometa
# Get rid of spaces after common in "Intersect" column
binary_wide_material_nometa$Intersect <- stri_replace_all_fixed(binary_wide_material_nometa$Intersect, " ", "")
head(binary_wide_material_nometa)
# Convert "Intersect" column to list format
binary_wide_material_nometa$Intersect <- strsplit(binary_wide_material_nometa$Intersect, ",")
head(binary_wide_material_nometa)

binary_tax_material_nometa <- binary_wide_material_nometa %>%
  rownames_to_column(var = "uniq") %>% 
  separate(uniq, c("ASV.ID", "Taxonomy"), sep = "_", remove = FALSE) %>% 
  # inner_join(tmp, by = "Taxonomy") %>% 
  column_to_rownames(var = "uniq") %>% 
  data.frame
head(binary_tax_material_nometa[1:2,])
# unique(binary_tax_material$tax_compiled)

totals_material_nometa<- binary_tax_material_nometa %>%
  group_by(Intersect) %>%
  summarise(total=sum(length(ASV.ID))) %>%
  data.frame
write.csv(as.data.frame(totals_material_nometa), "~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/totals_material.csv")
capture.output(summary(totals_material_nometa), file = "~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/totals_material.txt")
lapply(totals_material_nometa, function(x) write.table( data.frame(x), "~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/totals_material.csv"  , append= T, sep=',' ))

totals.num_material_nometa<-order(totals_material_nometa$total,decreasing=TRUE)
total.data_material_nometa<-data.frame(as.list(totals_material_nometa$Intersect),totals_material_nometa$total)

upset_plot<-binary_tax_material_nometa%>% ggplot(aes(x=Intersect, fill = Taxonomy)) +
  geom_bar(position="stack",stat = "count", color="black") + 
  # geom_text(stat="count", aes(label=after_stat(count), vjust=-1))+
  geom_text(data=totals_material_nometa ,aes(x=Intersect,y=total,label=total,fill=NULL),
         nudge_y = 10)+
  scale_x_upset()+scale_fill_manual(name="Taxa",values=c(tax_color_2))+
  theme_classic(base_size = 12)+xlab("")+ylab("Number of ASVs")+
  theme(legend.position="right",
        legend.text = element_text(size = 8),
        plot.margin = margin(10, 10, 10, 100))+
  guides(fill = guide_legend(ncol = 2))
# ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/UpSetR_nometa_sample_0323.pdf", width=12, height = 10)

#### ASV shared/unique, etc.####
# for shared between all- filter out the ASVs that are only found in all three samples
shared_asvs <- filter(binary_wide_nometa_forotherplots, Intersect=="Water column RNA, Water column DNA, Trap DNA")
shared_asvs$ASV1.ID <- row.names(shared_asvs) # add back the ASVs names
shared_asvs_names <- colsplit(shared_asvs$ASV1.ID, "_", names=c("ASV1.ID_real", "Taxon")) # split the ASV names and Taxon group
# filter the whole dataset by the shared ASVs
shared_asvs_only <- filter(nometa, ASV1.ID %in% shared_asvs_names$ASV1.ID_real)
sum_shared_each <- colSums(shared_asvs_only[,2:61]) # this gets the number of shared reads in all the samples
sum_total_each <-colSums(nometa[,2:61]) # this gets the total number of reads in all the samples
percent_shared <- sum_shared_each/sum_total_each
comparison_all <- data.frame(names(percent_shared),percent_shared)
names(comparison_all)[1]<-"Sample"
head(comparison_all)
comparison_plot <- join(comparison_all, key_all, by="Sample")
head(comparison_plot)

kruskal.test(comparison_plot$percent_shared,comparison_plot$Sample.type)
pairwise.wilcox.test(comparison_plot$percent_shared, comparison_plot$Sample.type, p.adjust="none")

percent_shared<-ggplot(comparison_plot, aes(x=Sample.type, y=percent_shared, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material)+
  theme_bw()+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  ylab("")+
  xlab("")+
  coord_cartesian(ylim=c(0,1))+
  labs(title="")+
  guides(fill=guide_legend(title="Sample Type"))+
  theme(#axis.text.x = element_text(hjust=1,vjust=0.5,size=8),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=12),legend.position = "right")
percent_shared + stat_compare_means(method = "kruskal.test")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/percent_shared_byall_nometa.pdf", width=4, height = 6)

pairwise.wilcox.test(comparison_plot$easy_comparison, comparison_plot$Depth, p.adjust="none")
# sig dif between 500m and all other depths

# plot of the number of ASVs per group
shared_barplot<-ggplot(shared_asvs_names, aes(x="", fill = Taxon)) +
  geom_bar(position="stack",stat = "count", color="black") + 
  scale_fill_manual(name="Taxa",values=c(tax_color_2))+
  scale_y_continuous(limits=c(0,1000))+
  theme_classic(base_size = 12)+xlab("")+ylab("Number of ASVs")+
  theme(legend.position = "none")

ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/shared_reads.pdf",  width=5, height = 7, units="in", dpi=300)

# plot the bar plot within the upsetr plot
plot.with.inset <-
  ggdraw() +
  draw_plot(upset_plot) +
  draw_plot(shared_barplot, x = .45, y = .3, width = .15, height = .5)
plot.with.inset
ggsave(filename = "plot.with.inset.pdf", 
       plot = plot.with.inset,
       width = 7, 
       height = 6,
       units = "in",
       dpi = 300)
shared_asvs_only$TAX_ORDER<-factor(shared_asvs_only$Taxonomy, levels = (tax_order), labels = (tax_order))
#### find the percentage of each shared group
shared.m<-melt(shared_asvs_names) #melt
head(unique.m)
unique.m$Sample<-factor(unique.m$variable,levels=names(unique_ASVs[2:61]))
unique_df_trial <- join(unique.m, key_all, by="Sample", type="left", match="first") %>%
  arrange(.,Station,Material,Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))
# Grouping by material, depth, station and taxa
asv_shared_total <- shared_asvs_names %>%
  group_by(Taxon) %>%
  count()# %>%
shared_asv_num=sum(asv_shared_total$n)

alv_shared  <- asv_shared_total %>%
  subset(., (Taxon == "Alveolates-Ciliates"| Taxon == "Alveolates-Dinophyceae" | Taxon=="Alveolates-Syndiniales"  | Taxon=="Alveolates-Dinoflagellates" ))
sum_alv_shared=sum(alv_shared$n)
percent_alv_shared = sum_alv_shared/shared_asv_num

stram_shared <- asv_shared_total %>%
  subset(., (Taxon == "Stramenopiles-Chrysophytes" | Taxon == "Stramenopiles-Diatoms"|Taxon == "Stramenopiles-MAST"|
               Taxon == "Stramenopiles-Ochrophyta" | Taxon == "Stramenopiles-Other" | Taxon == "Stramenopiles-Pelagophytes"))
sum_stram_shared=sum(stram_shared$n)
percent_stram_shared = sum_stram_shared/shared_asv_num

rhiz_shared <- asv_shared_total %>%
  subset(., (Taxon == "Rhizaria-Acantharia"| Taxon == "Rhizaria-Cercozoa" | Taxon == "Rhizaria-Other" | Taxon == "Rhizaria-Polycystines"))
sum_rhiz_shared=sum(rhiz_shared$n)
percent_rhiz_shared = sum_rhiz_shared/shared_asv_num



totals_material_nometa<- binary_tax_material_nometa %>%
  group_by(Intersect) %>%
  summarise(total=sum(length(ASV.ID))) %>%
  data.frame
write.csv(as.data.frame(totals_material_nometa), "~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/totals_material.csv")
capture.output(summary(totals_material_nometa), file = "~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/totals_material.txt")
lapply(totals_material_nometa, function(x) write.table( data.frame(x), "~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/totals_material.csv"  , append= T, sep=',' ))

totals.num_material_nometa<-order(totals_material_nometa$total,decreasing=TRUE)
total.data_material_nometa<-data.frame(as.list(totals_material_nometa$Intersect),totals_material_nometa$total)

ggplot(shared_asvs_only, aes(x=, fill = Taxonomy)) +
  geom_bar(position="stack",stat = "count", color="black") + 
  #geom_text(aes(y=sum(stat(count)), group=Intersect), vjust=0)+
  #geom_text(data=totals_material_nometa ,aes(x=Intersect,y=total,label=total,fill=NULL),
  #  nudge_y = 10)+
  scale_x_upset()+scale_fill_manual(name="Taxa",values=c(tax_color_phyla2_nometa))+
  theme_classic(base_size = 12)+xlab("")+ylab("Number of Shared ASVs")+
  theme(legend.position="right",
        legend.text = element_text(size = 8),
        plot.margin = margin(10, 10, 10, 100))

### now want to do the opposite, just the unique ones
DNA_only <- filter(binary_wide_nometa_forotherplots, Intersect=="Water column DNA")
RNA_only <- filter(binary_wide_nometa_forotherplots, Intersect=="Water column RNA")
PIT_only <- filter(binary_wide_nometa_forotherplots, Intersect=="Trap DNA")
all_unique<- rbind(DNA_only, RNA_only, PIT_only)
all_unique$ASV1.ID <- row.names(all_unique) # add back the ASVs names
unique_names <- colsplit(all_unique$ASV1.ID, "_", names=c("ASV1.ID_real", "Taxon")) # split the ASV names and Taxon group
# filter the whole dataset by the unique ASVs
unique_ASVs <- filter(nometa, ASV1.ID %in% unique_names$ASV1.ID_real)
sum_unique_each <- colSums(unique_ASVs[,2:61]) # this gets the number of shared reads in all the samples
percent_unique_each <- sum_unique_each/sum_total_each
comparison_unique <- data.frame(names(percent_unique_each),percent_unique_each)
names(comparison_unique)[1]<-"Sample"
head(comparison_unique)
comparison_unique_plot <- join(comparison_unique, key_all, by="Sample")
head(comparison_unique_plot)

percent_unique<-ggplot(comparison_unique_plot, aes(x=Sample.type, y=percent_unique_each, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material)+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  ylab("")+
  xlab("")+
  labs(title="")+
  guides(fill=guide_legend(title="Sample Type"))+
  theme(#axis.text.x = element_text(hjust=1,vjust=0.5,size=8),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=12),legend.position = "right")
percent_unique + stat_compare_means(method = "kruskal.test")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/percent_unique_byall_nometa.pdf", width=4, height = 6)
kruskal.test(comparison_unique_plot$percent_unique_each,comparison_unique_plot$Sample.type)
pairwise.wilcox.test(comparison_unique_plot$percent_unique_each, comparison_unique_plot$Sample.type, p.adjust="none")
##### percentage of ASVS within unique ######
unique.m<-melt(unique_ASVs) #melt
head(unique.m)
unique.m$Sample<-factor(unique.m$variable,levels=names(unique_ASVs[2:61]))
unique_df_trial <- join(unique.m, key_all, by="Sample", type="left", match="first") %>%
  arrange(.,Station,Material,Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))
# Grouping by material, depth, station and taxa
asv_unique_total <- unique_df_trial %>%
  group_by(Material, Depth, Station,Taxa) %>%
  summarise(count = sum(value > 0)) %>% # get the number of rows (ASVs) where there is > 0
  group_by(Material,Depth,Station) %>% 
  summarise(total_asv = sum(count)) # sum the total ASVs > 0
  
asv_trial <- unique_df_trial %>%
  group_by(Material, Depth, Station,Taxa) %>%
  summarise(count = sum(value > 0)) # get the number of ASVs where there is > 0

sum_alv_unique  <- asv_trial %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Alveolates-Ciliates"| Taxa == "Alveolates-Dinophyceae" | Taxa=="Alveolates-Syndiniales"  | Taxa=="Alveolates-Dinoflagellates" )) %>%
  summarise(total_alv = sum(count))

sum_stram_unique <- asv_trial %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Stramenopiles-Chrysophytes" | Taxa == "Stramenopiles-Diatoms"|Taxa == "Stramenopiles-MAST"|
               Taxa == "Stramenopiles-Ochrophyta" | Taxa == "Stramenopiles-Other" | Taxa == "Stramenopiles-Pelagophytes")) %>%
  summarise(total_stram = sum(count))

sum_rhiz_unique <- asv_trial %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Rhizaria-Acantharia"| Taxa == "Rhizaria-Cercozoa" | Taxa == "Rhizaria-Other" | Taxa == "Rhizaria-Polycystines")) %>%
  summarise(total_rhiz = sum(count))

percent_all_unique_2 <- data.frame(asv_unique_total, sum_alv_unique$total_alv, sum_stram_unique$total_stram,
                                   sum_rhiz_unique$total_rhiz) %>%
  mutate(percent_alv=sum_alv_unique.total_alv/total_asv,
         percent_stram = sum_stram_unique.total_stram/total_asv,
         percent_rhiz=sum_rhiz_unique.total_rhiz/total_asv)

percent_unique_means_2 <- percent_all_unique_2 %>%
  group_by(Material)%>%
  summarise(mean_alv = mean(percent_alv),
            stdev_alv = sd(percent_alv),
            mean_stram = mean(percent_stram),
            stdev_stram = sd(percent_stram),
            mean_rhiz = mean(percent_rhiz),
            stdev_rhiz = sd(percent_rhiz))

##### percentage of reads within unique, not used #######
unique.m<-melt(unique_ASVs) #melt
head(unique.m)
unique.agg<-aggregate(unique.m$value, by=list(Taxa=unique.m$TaxaPlot,Sample=unique.m$variable),sum) #sum sequences by taxonomic group
unique.agg$Sample<-factor(unique.agg$Sample,levels=names(unique_ASVs[2:61]))
unique_df <- join(unique.agg, key_all, by="Sample", type="left", match="first") %>%
  arrange(.,Station,Material,Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))

# Step 1. need to sum the number of reads per material type and depth and station
percent_unique_sample <- unique_df %>%
  group_by(Material, Depth, Station) %>%
  summarise(total_unique = sum(x))

# Step 2. sum the number of reads in the group I'm interested
sum_alv_unique  <- unique_df %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Alveolates-Ciliates"| Taxa == "Alveolates-Dinophyceae" | Taxa=="Alveolates-Syndiniales" )) %>%
  summarise(total_alveolates = sum(x))

sum_rhiz_unique  <- unique_df %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Rhizaria-Acantharia"| Taxa == "Rhizaria-Other" | Taxa =="Rhizaria-Polycystines" | Taxa=="Rhizaria-Cercozoa")) %>%
  summarise(total_rhiz = sum(x))

sum_other_unique  <- unique_df %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Other/unknown")) %>%
  summarise(total_other = sum(x))

sum_cil_unique  <- unique_df %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Alveolates-Ciliates" )) %>%
  summarise(total_cil = sum(x))

sum_dino_unique  <- unique_df %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Alveolates-Dinophyceae")) %>%
  summarise(total_dino = sum(x))

sum_synd_unique  <- unique_df %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa =="Alveolates-Syndiniales")) %>%
  summarise(total_synd = sum(x))

sum_stram_unique <- unique_df %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Stramenopiles-Chrysophytes" | Taxa == "Stramenopiles-Diatoms"|Taxa == "Stramenopiles-MAST"|
               Taxa == "Stramenopiles-Ochrophyta" | Taxa == "Stramenopiles-Other" | Taxa == "Stramenopiles-Pelagophytes")) %>%
  summarise(total_stram = sum(x))

# Step 3. divide the sum of the group reads by the total numbers of reads
percent_all_unique <- data.frame(percent_unique_sample, sum_alv_unique$total_alveolates, sum_rhiz_unique$total_rhiz,
                                sum_other_unique$total_other, sum_stram_unique$total_stram,
                                sum_cil_unique$total_cil, sum_synd_unique$total_synd, sum_dino_unique$total_dino) %>%
  mutate(percent_alv=sum_alv_unique.total_alveolates/total_unique,
         percent_rhiz= sum_rhiz_unique.total_rhiz/total_unique,
         percent_other= sum_other_unique.total_other/total_unique,
         percent_cil= sum_cil_unique.total_cil/total_unique,
         percent_dino= sum_dino_unique.total_dino/total_unique,
         percent_synd= sum_synd_unique.total_synd/total_unique,
         percent_stram = sum_stram_unique.total_stram/total_unique
         )

percent_unique_means <- percent_all_unique %>%
  group_by(Material)%>%
  summarise(mean_alv = mean(percent_alv),
            stdev_alv = sd(percent_alv),
            mean_rhiz = mean(percent_rhiz),
            stdev_rhiz = sd(percent_rhiz),
            mean_other = mean(percent_other),
            stdev_other = sd(percent_other),
            mean_cil = mean(percent_cil),
            stdev_cil = sd(percent_cil),
            mean_synd = mean(percent_synd),
            stdev_synd = sd(percent_synd),
            mean_dino = mean(percent_dino),
            stdev_dino = sd(percent_dino),
            mean_stram = mean(percent_stram),
            stdev_stram = sd(percent_stram))


percent_pit_cylclonic <- filter(percent_all_micro, (Station=="PIT1" | Station =="PIT2")) %>%
  summarise(mean_alv = mean(percent_alv),
            stdev_alv = sd(percent_alv),
            mean_rhiz = mean(percent_rhiz),
            stdev_rhiz = sd(percent_rhiz),
            mean_stram = mean(percent_stram),
            stdev_stram = sd(percent_stram),
            mean_diatom = mean(percent_diatom),
            stdev_diatom = sd(percent_diatom))

percent_pit_anticylclonic <- filter(percent_all_micro, (Station=="PIT12" | Station =="PIT11" | Station=="PIT10" | Station =="PIT9" | Station=="PIT8" | Station =="PIT7" | Station=="PIT6")) %>%
  summarise(mean_alv = mean(percent_alv),
            stdev_alv = sd(percent_alv),
            mean_rhiz = mean(percent_rhiz),
            stdev_rhiz = sd(percent_rhiz),
            mean_stram = mean(percent_stram),
            stdev_stram = sd(percent_stram),
            mean_diatom = mean(percent_diatom),
            stdev_diatom = sd(percent_diatom))


###### finish plotting the shared/unique reads
patch <- percent_shared / percent_unique + plot_layout(guides="collect")
patch | shared_barplot + plot_annotation(tag_levels="a")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Final_Figs/shared_unique_2.pdf",  width=7, height = 5, units="in", dpi=300)

shared_by_two_ASVs <- filter(nometa, !(ASV1.ID %in% c(unique_names$ASV1.ID_real,shared_asvs_names$ASV1.ID_real))) 
sum_two_each <- colSums(shared_by_two_ASVs[,2:61]) # this gets the number of shared reads in all the samples
percent_two_each <- sum_two_each/sum_total_each
comparison_two <- data.frame(names(percent_two_each),percent_two_each)
names(comparison_two)[1]<-"Sample"
head(comparison_two)
comparison_two_plot <- join(comparison_two, key_all, by="Sample")
head(comparison_two_plot)

percent_two<-ggplot(comparison_two_plot, aes(x=Sample.type, y=percent_two_each, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material)+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  ylab("")+
  xlab("")+
  labs(title="")+
  guides(fill=guide_legend(title="Sample Type"))+
  theme(#axis.text.x = element_text(hjust=1,vjust=0.5,size=8),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=12),legend.position = "right")
percent_two + stat_compare_means(method = "kruskal.test")
# trying to get the percentage of shared just between the water columns
DNA_RNAonly <- filter(binary_wide_nometa_forotherplots, Intersect=="Water column RNA, Water column DNA")
DNA_RNAonly$ASV1.ID <- row.names(DNA_RNAonly) # add back the ASVs names
wc_shared_names <- colsplit(DNA_RNAonly$ASV1.ID, "_", names=c("ASV1.ID_real", "Taxon")) # split the ASV names and Taxon group
# filter the whole dataset by the unique ASVs
WC_only_ASVs <- filter(nometa, ASV1.ID %in% wc_shared_names$ASV1.ID_real)
wc.m<-melt(WC_only_ASVs ) #melt
head(wc.m)
wc.m$Sample<-factor(wc.m$variable,levels=names(WC_only_ASVs[2:61]))
wc_df <- join(wc.m, key_all, by="Sample", type="left", match="first") %>%
  arrange(.,Station,Material,Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))
# Grouping by material, depth, station and taxa
asv_wc_total <- wc_df %>%
  group_by(Material, Depth, Station,Taxa) %>%
  summarise(count = sum(value > 0)) %>% # get the number of rows (ASVs) where there is > 0
  group_by(Material,Depth,Station) %>% 
  summarise(total_asv = sum(count)) # sum the total ASVs > 0

asv_wc <- wc_df%>%
  group_by(Material, Depth, Station,Taxa) %>%
  summarise(count = sum(value > 0)) # get the number of ASVs where there is > 0

sum_alv_wc  <- asv_wc %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Alveolates-Ciliates"| Taxa == "Alveolates-Dinophyceae" | Taxa=="Alveolates-Syndiniales"  | Taxa=="Alveolates-Dinoflagellates" )) %>%
  summarise(total_alv = sum(count))

sum_stram_wc <- asv_wc %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Stramenopiles-Chrysophytes" | Taxa == "Stramenopiles-Diatoms"|Taxa == "Stramenopiles-MAST"|
               Taxa == "Stramenopiles-Ochrophyta" | Taxa == "Stramenopiles-Other" | Taxa == "Stramenopiles-Pelagophytes")) %>%
  summarise(total_stram = sum(count))

sum_rhiz_wc <- asv_wc %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Rhizaria-Acantharia"| Taxa == "Rhizaria-Cercozoa" | Taxa == "Rhizaria-Other" | Taxa == "Rhizaria-Polycystines")) %>%
  summarise(total_rhiz = sum(count))

percent_wc <- data.frame(asv_wc_total, sum_alv_wc$total_alv, sum_stram_wc$total_stram,
                                   sum_rhiz_wc$total_rhiz) %>%
  mutate(percent_alv=sum_alv_wc.total_alv/total_asv,
         percent_stram = sum_stram_wc.total_stram/total_asv,
         percent_rhiz=sum_rhiz_wc.total_rhiz/total_asv)

percent_wc_means<- percent_wc %>%
  subset(., Material %in% WC)%>%
  summarise(mean_alv = mean(percent_alv),
            stdev_alv = sd(percent_alv),
            mean_stram = mean(percent_stram),
            stdev_stram = sd(percent_stram),
            mean_rhiz = mean(percent_rhiz),
            stdev_rhiz = sd(percent_rhiz))

# trying to combine them
unique_two <- join(comparison_unique,comparison_two, by="Sample", type="left")
all_groups <- join(unique_two, comparison_all, by="Sample", type="left")
head(all_groups)
combo_melt <- melt(all_groups)
all_plotting <- join(combo_melt, key_all, by="Sample", type="left")
all_plotting$Sample.type <- factor(all_plotting$Sample.type, levels=c("Water column RNA","Water column DNA","Trap DNA") )
box_shared<-ggplot(all_plotting, aes(x=Sample.type, y=value, fill=variable))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=c("#e0ecf4","#9ebcda","#8856a7"), labels=c("Unique", "Shared by two", "Shared by all"))+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  labs(title="", x="", y="Percent")+
  guides(fill=guide_legend(title=""))+
  theme(axis.text.x = element_text(angle=45,hjust=.9,vjust=.9,size=8),
    #axis.text.x=element_blank(),
    axis.text.y=element_text(size=12),legend.position = "right")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Final_Figs/percent_shared.pdf",  width=7, height = 5, units="in", dpi=300)

#dunn test
dunn_test_asvs<-all_plotting %>%
  group_by(variable)%>%
  dunn_test(value~Material)
write.csv(dunn_test_asvs, "~/Desktop/Chapter1/MESOSCOPE2017/Output/Final_Figs/dunn_shared_boxplot.csv")
 
plot.with.inset / (venn | box_shared) + plot_annotation(tag_levels="a")
 ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/fig6_2.pdf", width=10, height = 11)
 
####### Figure 7. Microscopy flux ######
flux<-read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/CellsPit_Final.csv",header=T)

rownames(flux)=flux$Type
# combine ciliates and tintinnid to the same category
flux$Plotting = flux$Type
flux$Plotting[flux$Type == "Tintinnid"]="Alveolates-Ciliates"
flux$Plotting[flux$Type == "Ciliate"]="Alveolates-Ciliates"
flux$Plotting[flux$Type == "Dinoflagellate"]="Alveolates-Dinoflagellates"
# combine acanths and rads to one category
flux$Plotting[flux$Type == "Rhizaria-Acantharia"]="Rhizaria-Radiolaria"
flux$Plotting[flux$Type == "Rhizaria-other"]="Rhizaria-Other"
flux$Plotting[flux$Type == "Diatom"]="Stramenopiles-Diatoms"
flux$Plotting[flux$Type == "Dictyocha"]="Stramenopiles-Dictyocha"

flux.m<-melt(select(flux, -c("Type")), id.vars = "Plotting") %>%
  separate(variable, c("Station", "Rep"), sep = "_", remove = TRUE)
head(flux.m)

#To calculate mean values of taxa by PIT
flux.mean<- flux.m %>% 
  group_by(Plotting, Station) %>%
  summarise(mean=mean(as.numeric(value)),stdev=sd(as.numeric(value))) %>%
  as.data.frame
flux.mean

# make the color match the other type
tax_order.counts<-c("Alveolates-Ciliates","Alveolates-Dinoflagellates","Rhizaria-Foraminifera","Rhizaria-Other","Rhizaria-Radiolaria", "Stramenopiles-Diatoms","Stramenopiles-Dictyocha")
col.pit.counts <-c('firebrick4','indianred1','#fde0dd','#db6ea0','lightpink','#DDAD4B','gold1')

flux.mean$tax.order<-factor(flux.mean$Plotting,levels=(tax_order.counts), labels=(tax_order.counts))
flux.mean$Station<-factor(flux.mean$Station, levels=c("PIT12","PIT11","PIT10","PIT9","PIT8","PIT7","PIT6","PIT5","PIT4","PIT3","PIT2","PIT1"))

pit_colors <- c("Alveolates-Ciliates"='firebrick4',"Alveolates-Dinoflagellates"='indianred1',"Other"= 'grey',
                "Rhizaria-Foraminifera"='mediumvioletred',"Rhizaria-Radiolaria"='lightpink',"Rhizaria-Other"='#C291A4',
                "Stramenopiles-Diatom"='#DDAD4B',"Stramenopiles-Dictyocha"="#5C4033")

flux_mean_ggplot <- flux.m %>% 
  group_by(Plotting, Station) %>%
  mean_se(value)

#To plot flux values
scatter<-ggplot(flux.mean, aes(y=mean,ymin= mean-stdev, ymax= mean + stdev, x=Station,fill=tax.order))+
  geom_errorbar(width = 0.2) +
  geom_point(size=3.5, shape=21,color="black")+
  scale_fill_manual(values=col.pit.counts, name="Taxa")+
  # scale_y_continuous(trans = 'log10',labels = function(x) format(x, scientific = TRUE))+
  theme_bw()+
  theme(panel.grid.minor = element_line(colour="white")) +
  labs(y=expression("Mean flux by microscopy of trap material (Cells m"^{"-2"}*" d"^{"-1"}*")"), x="")+
  # theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        # legend.title = element_blank(),
        legend.text = element_text(size = 10),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5))+
  facet_wrap(~tax.order, ncol=2, scales="free_y")
scatter %+% facet_wrap(~tax.order, ncol = 2) 
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/mean_flux_facet_stdev_samescale_0323.pdf", width=6.5, height = 7.5)
scatter %+% facet_wrap(~tax.order, ncol = 2, scales="free_y") 
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/mean_flux_facet_stdev_difscale_0323.pdf", width=6.5, height = 7.5)

scatter
ggplot2::ggsave("mean_flux_stdev.pdf", width=5, height = 5)
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/mean_flux_stdev.pdf", width=10, height = 4)

flux_sla <- join(flux.mean, select(key_all, c(Sample, Station, SLA)), by="Station",type="left")

scatter_sla<-ggplot(flux_sla, aes(y=mean,ymin= mean-stdev, ymax= mean + stdev, x=SLA,fill=tax.order))+
  geom_errorbar(width = 0.2) +
  geom_point(size=3.5, shape=21,color="black")+
  geom_smooth(method="lm",color="black",linewidth=0.5)+
  scale_fill_manual(values=col.pit.counts, name="Taxa")+
  scale_y_continuous(trans = 'log10',labels = function(x) format(x, scientific = TRUE))+
  theme_bw()+
  theme(panel.grid.minor = element_line(colour="white")) +
  labs(y=expression("Mean flux by microscopy of trap material (Cells m"^{"-2"}*" d"^{"-1"}*")"), x="SLAcorr(cm)")+
  # theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        # legend.title = element_blank(),
        legend.text = element_text(size = 10),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5))
scatter_sla %+% facet_wrap(~tax.order, ncol = 2, scales="free_y") 
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Final_Figs/flux_sla.pdf", width=6.5, height = 7.5)
# to get the model information for alll the types
model_type1 <- flux_sla %>%
  group_by(Plotting) %>%
  do(model = lm(SLA ~ mean,data=.))

# https://broom.tidymodels.org/reference/tidy.lmodel2.html

# create a type 2 regression for all the factors
model_type2 <- flux_sla %>%
  group_by(Plotting) %>%
  do(model = lmodel2(mean ~ SLA,data=.))
model_cil <-model_type2[1,]$model[[1]][[3]]

# trying to plot just the ciliates for proof of concept 
names(model_cil) <- c("method", "intercept", "slope", "angle", "p-value")
model_cil$Plotting<-"Alveolates-Ciliates"
flux_ciliates<- filter(flux_sla, Plotting=="Alveolates-Ciliates")
ggplot(flux_ciliates, aes(x=SLA, y=mean)) + 
  xlab("SLA") + ylab("flux (cells/L/day)") + 
  geom_point(aes(col=Plotting, size=5))+
  #geom_abline(data = model_cil, aes(intercept = intercept, slope = slope), show.legend = TRUE)
  geom_abline(data = subset(model_3, method=="SMA"), aes(intercept = intercept, slope = slope), show_guide = TRUE)+
  stat_regline_equation(label.x = .1, label.y = 1000) + stat_cor(label.x = .1, label.y = 1500)

# now do the rest
model_dino<-model_type2[2,]$model[[1]][[3]]
names(model_dino) <- c("method", "intercept", "slope", "angle", "p-value")
model_dino$Plotting<-"Alveolates-Dinoflagellates"
model_foram <-model_type2[3,]$model[[1]][[3]]
names(model_foram) <- c("method", "intercept", "slope", "angle", "p-value")
model_foram$Plotting<-"Rhizaria-Foraminifera"
model_rhizother <-model_type2[4,]$model[[1]][[3]]
names(model_rhizother) <- c("method", "intercept", "slope", "angle", "p-value")
model_rhizother$Plotting<-"Rhizaria-Other"
model_rad <-model_type2[5,]$model[[1]][[3]]
names(model_rad) <- c("method", "intercept", "slope", "angle", "p-value")
model_rad$Plotting<-"Rhizaria-Radiolaria"
model_diatom <-model_type2[6,]$model[[1]][[3]]
names(model_diatom) <- c("method", "intercept", "slope", "angle", "p-value")
model_diatom$Plotting<-"Stramenopiles-Diatoms"
model_dyct <-model_type2[7,]$model[[1]][[3]]
names(model_dyct) <- c("method", "intercept", "slope", "angle", "p-value")
model_dyct$Plotting<-"Stramenopiles-Dictyocha"

model_plot<-rbind(model_cil, model_dino) %>%
  rbind(.,model_foram)%>%
  rbind(.,model_rhizother) %>%
  rbind(., model_rad)%>%
  rbind(., model_diatom) %>%
  rbind(., model_dyct)

scatter_sla_model2<-ggplot(flux_sla, aes(y=mean,ymin= mean-stdev, ymax= mean+stdev, x=SLA,fill=tax.order))+
  geom_errorbar(width = 0.2) +
  geom_point(size=3.5, shape=21,color="black")+
  scale_fill_manual(values=col.pit.counts, name="Taxa")+
  theme_bw()+
  theme(panel.grid.minor = element_line(colour="white")) +
  labs(y=expression("Mean flux by microscopy of trap material (Cells m"^{"-2"}*" d"^{"-1"}*")"), x="SLAcorr(cm)")+
  theme(legend.position="right",
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5))+
  geom_abline(data = subset(model_plot, method=="SMA"), aes(intercept = intercept, slope = slope), show_guide = TRUE) # +
  stat_regline_equation(label.x = .1, label.y = 1000) + stat_cor(label.x = .1, label.y = 1500)
scatter_sla_model2 %+% facet_wrap(~Plotting, ncol = 2, scales="free_y")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Final_Figs/flux_sla_2.pdf", width=6.5, height = 7.5)

# R values: 
scatter_sla_model2_rvalue<-ggplot(flux_sla, aes(y=mean,ymin= mean-stdev, ymax= mean+stdev, x=SLA,fill=tax.order))+
  geom_errorbar(width = 0.2) +
  geom_point(size=3.5, shape=21,color="black")+
  scale_fill_manual(values=col.pit.counts, name="Taxa")+
  theme_bw()+
  theme(panel.grid.minor = element_line(colour="white")) +
  labs(y=expression("Mean flux by microscopy of trap material (Cells m"^{"-2"}*" d"^{"-1"}*")"), x="SLAcorr(cm)")+
  theme(legend.position="right",
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5))+
  geom_abline(data = subset(model_plot, method=="SMA"), aes(intercept = intercept, slope = slope), show_guide = TRUE) +
  #stat_regline_equation(label.x = .1, label.y = 1000) + 
  stat_cor(label.x = .1, label.y = 1500)
scatter_sla_model2_rvalue %+% facet_wrap(~Plotting, ncol = 2, scales="free_y")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Final_Figs/flux_sla_2_r.pdf", width=6.5, height = 7.5)


####### Figure S2. Diversity metrics & associated stats#######
# calculate metrics
shannon_nometa<-diversity(nometa[,2:61], index = "shannon",MARGIN = 2)
invsimp_nometa<-diversity(nometa[,2:61], index = "invsimpson",MARGIN = 2)
ASV_count_nometa<-colSums(nometa[,2:61]>0)
ASV_number_nometa<-colSums(nometa[,2:61])
Num_Seqs_stats <- summary(ASV_number_nometa)
num_seqs_sd <- sd(ASV_number_nometa)

#combine metrics
alpha_nometa<-data.frame(ASV_number_nometa, ASV_count_nometa, shannon_nometa, invsimp_nometa)
head(alpha_nometa)
# add in other variables
alpha_nometa$Sample<-row.names(alpha_nometa)
alpha.df_nometa<- alpha_nometa %>%
  melt(.) %>%
  join(.,key_all, by="Sample", type="left", match="first")
head(alpha.df_nometa)
alpha_nometa_means <- alpha.df_nometa %>%
  group_by(Material, variable)%>%
  summarise(mean = mean(value), stdev=sd(value))
write.csv(alpha_nometa_means, "~/Desktop/Chapter1/MESOSCOPE2017/Output/Final_Figs/Supplemental2_means_sd.csv")
dunn_test<-alpha.df_nometa %>%
  group_by(variable)%>%
  dunn_test(value~Material)
write.csv(dunn_test, "~/Desktop/Chapter1/MESOSCOPE2017/Output/Final_Figs/Supplemental2_stats_dunntest.csv")

station_dunn<-alpha.df_nometa %>%
  group_by(variable, Material)%>%
  dunn_test(value~Station)

#plot the number of reads
seq_num <-ggplot(alpha.df_nometa[alpha.df_nometa$variable=="ASV_number_nometa",], aes(x=Sample.type, y=value, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material, name="Sample Type")+
  theme_bw()+
  labs(title="", x="", y="Number of reads")+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  theme(#axis.text.x = element_text(angle=45,hjust = .9, vjust = .9,size=8),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=12),legend.position = "right")
seq_num+ stat_compare_means(method = "kruskal.test")

# plot the number of ASVs
alpha<-ggplot(alpha.df_nometa[alpha.df_nometa$variable=="ASV_count_nometa",], aes(x=Sample.type, y=value, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material, name="Sample Type")+
  theme_bw()+
  labs(title="", x="", y="Number of ASVs")+
  scale_y_continuous(breaks=seq(0,3000,1000), limits=c(0,3000)) +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  theme(#axis.text.x = element_text(angle=45,hjust = .9, vjust = .9,size=8),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=12),legend.position = "right")
alpha + stat_compare_means(method = "kruskal.test")

# plot shannon's index
shannons<-ggplot(alpha.df_nometa[alpha.df_nometa$variable=="shannon_nometa",], aes(x=Sample.type, y=value, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material, name="Sample Type")+
  theme_bw()+
  labs(title="", x="", y="Shannon's Index")+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  theme(#axis.text.x = element_text(angle=45,hjust = .9, vjust = .9,size=8),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=12),legend.position = "right")
shannons+ stat_compare_means(method = "kruskal.test")

# plot inverse simpson's
inv_simp<-ggplot(alpha.df_nometa[alpha.df_nometa$variable=="invsimp_nometa",], aes(x=Sample.type, y=value, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material, name="Sample Type")+
  theme_bw()+
  labs(title="", x="", y="Inverse Simpon's")+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  theme(#axis.text.x = element_text(angle=45,hjust = .9, vjust = .9,size=8),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=12),
        legend.position = "right")
inv_simp + stat_compare_means(method = "kruskal.test")


# combine plots
((seq_num + alpha) / (shannons + inv_simp)) + plot_layout(guides="collect") + plot_annotation(tag_levels="a")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Final_Figs/Supplemental2.pdf",  width=7, height = 7, units="in", dpi=300)
