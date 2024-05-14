## Figure 9- shared and unique ASVs, UpsetR, venn and boxplot
## Use eulerr for venn diagram, ggupset for upset plot
## J. L. Beatty, inspiration from Samantha Gleich, Gerid Ollison, Sarah Hu
## Last updated: 05/2024

# load packages
library(plyr)
library(tidyverse)
library(reshape2)
library(eulerr) # for venn diagram
library(stringi) 
library(ggupset) 
library(rstatix)
library(cowplot)

# load data and metadata
asv_df <- read.csv("~/github/mesoscope2017/processed-data/MS_Clean_No1_Means_ALL_Norm_NewTax_0424.csv", header=T)
asv_df[1]=NULL
row.names(asv_df) = asv_df$ASV1.ID
nometa <- asv_df[asv_df$Level3 != "Metazoa",] # remove metazoan asvs

#metadata
key_all<-read.delim("~/github/mesoscope2017/raw-data/variables_all.txt")
key_all$Sample.type <- factor(as.character(key_all$Sample.type), levels=c( "Water RNA", "Water DNA", "Trap DNA"))

# for plotting
WC <- c("RNA","DNA")
col_depth <-c("15m" = "#fee08b", "DCM" = "#7fbc41", "175m" = "#74add1", "500m" = "#542788")
col_material<- c("Water RNA"="#7fcdbb","Water DNA"="#2c7fb8","Trap DNA"="#fc8d59")
tax_color_phyla2_nometa<- c("Alveolata-Ciliates"='firebrick4',"Alveolata-Dinoflagellates"='indianred1',"Alveolata-Syndiniales"='tomato3', "Archaeplastida-Chlorophytes"='forestgreen',
                            "Excavata-Discobids"='yellowgreen', "Hacrobia-Cryptophytes"='darkblue',"Hacrobia-Haptophytes"='lightblue',
                            "Opisthokonta-Choanoflagellates"='moccasin',"Other/unknown"='grey',
                            "Rhizaria-Acantharians"='mediumvioletred',"Rhizaria-Cercozoans"="#C291A4","Rhizaria-Other"='#db6ea0',"Rhizaria-Polycystines"='lightpink', 
                            "Stramenopila-Chrysophytes"="#DECDBE","Stramenopila-Diatoms"='#DDAD4B',"Stramenopila-MAST"='tan2',"Stramenopila-Ochrophytes"='tan3',"Stramenopila-Other"='tan4',"Stramenopila-Pelagophytes"="#5C4033")


# first reshape data and add new labels
melt_ASV_nometa<-melt(nometa)
names(melt_ASV_nometa)[14]="Sample"
melt_all_nometa<-join(melt_ASV_nometa,key_all,by="Sample", type="left", match="first")
names(melt_all_nometa)
toUpset_nometa<-melt_all_nometa[c(1:2,14:18,22)] # We want the ASV ID, Sample, TaxaPlot, value, Station, Depth, Material, Sample type
head(toUpset_nometa) # Use below to add annotation
toUpset_nometa$Uniq<-paste(toUpset_nometa$ASV1.ID, toUpset_nometa$TaxaPlot, sep="_") # creating a unique combination of ASV and taxa ids

# Change to binary
summed_material_nometa<-aggregate(toUpset_nometa$value, by=list(Sample_type=toUpset_nometa$Sample.type, Uniq=toUpset_nometa$Uniq),sum)
summed_material_nometa$bin<-ifelse(summed_material_nometa$x > 0, 1, 0) # assign a 1 value if any reads are present, if not zero
binary_wide_material_nometa<-dcast(summed_material_nometa[c(1,2,4)], Uniq~Sample_type, fill=0) # change the data shape
row.names(binary_wide_material_nometa)<-binary_wide_material_nometa$Uniq; binary_wide_material_nometa$Uniq<-NULL # assign rownames and remove null values
head(binary_wide_material_nometa)

#### Figure 8b- Venn diagram of ASVs found in different material type ####
# create lists of lists that are the rownames (aka individual asvs) bserved in each separate group 
venn_data <- list(
  "Water RNA" = row.names(binary_wide_material_nometa[binary_wide_material_nometa$`Water RNA`>0,]),
  "Water DNA" = row.names(binary_wide_material_nometa[binary_wide_material_nometa$`Water DNA`>0,]),
  "Trap DNA" = row.names(binary_wide_material_nometa[binary_wide_material_nometa$`Trap DNA`>0,])
)
venn1<-euler(venn_data)
venn2<- plot(venn1, 
              fill = c("#7fcdbb","#2c7fb8","#fc8d59"), 
              alpha=0.5,
              quantities=list(type=c("counts","percent")))
#####
# Make "Intersect" column to tell us which material set a given ASV shows up in
binary_wide_material_nometa$Intersect <- apply(binary_wide_material_nometa > 0, 1, function(x){toString(names(binary_wide_material_nometa)[x])})
head(binary_wide_material_nometa)
binary_wide_nometa_forotherplots <- binary_wide_material_nometa

### Plot Figure 8a- UpsetR plot ####
# Get rid of spaces after common in "Intersect" column
binary_wide_material_nometa$Intersect <- stri_replace_all_fixed(binary_wide_material_nometa$Intersect, " ", "")
head(binary_wide_material_nometa)
# Convert "Intersect" column to list format
binary_wide_material_nometa$Intersect <- strsplit(binary_wide_material_nometa$Intersect, ",")
head(binary_wide_material_nometa)

binary_tax_material_nometa <- binary_wide_material_nometa %>%
  rownames_to_column(var = "uniq") %>% 
  separate(uniq, c("ASV.ID", "Taxonomy"), sep = "_", remove = FALSE) %>% 
  column_to_rownames(var = "uniq") %>% 
  data.frame
head(binary_tax_material_nometa[1:2,])

totals_material_nometa<- binary_tax_material_nometa %>%
  group_by(Intersect) %>%
  summarise(total=sum(length(ASV.ID))) %>%
  data.frame

upset_plot<-binary_tax_material_nometa%>% ggplot(aes(x=Intersect, fill = Taxonomy)) +
  geom_bar(position="stack",stat = "count", color="black") + 
  geom_text(data=totals_material_nometa ,aes(x=Intersect,y=total,label=total,fill=NULL),
            nudge_y = 10)+
  scale_x_upset()+scale_fill_manual(name="Taxa",values=c(tax_color_phyla2_nometa))+
  theme_classic(base_size = 12)+xlab("")+ylab("Number of ASVs") +
  theme(legend.position="right",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.25,'cm'),
        text = element_text(size=8),
        plot.margin = margin(0,0,0,0),
        legend.margin=margin(0,0,0,0))
# add the inset bar later

#### Figure 8c- shared and unique ASV boxplot ####
# filter out the ASVs that are only found in all three samples
shared_asvs <- filter(binary_wide_nometa_forotherplots, Intersect=="Water RNA, Water DNA, Trap DNA")
shared_asvs$ASV1.ID <- row.names(shared_asvs) # add back the ASVs names
shared_asvs_names <- colsplit(shared_asvs$ASV1.ID, "_", names=c("ASV1.ID_real", "Taxon")) # split the ASV names and Taxon group
#plot of the relative abundance of each group within the reads, to be inset later 
shared_barplot<-ggplot(shared_asvs_names, aes(x="", fill = Taxon)) +
  geom_bar(position="stack",stat = "count", color="black") + 
  scale_fill_manual(name="Taxa",values=c(tax_color_phyla2_nometa))+
  scale_y_continuous(limits=c(0,1000))+
  theme_classic(base_size = 12)+xlab("")+ylab("Number of ASVs")+
  theme(legend.position = "none")

# filter the whole dataset by the shared ASVs
shared_asvs_only <- filter(nometa, ASV1.ID %in% shared_asvs_names$ASV1.ID_real)
sum_shared_each <- colSums(shared_asvs_only[,1:60]) # this gets the number of shared reads in all the samples
sum_total_each <-colSums(nometa[,1:60]) # this gets the total number of reads in all the samples
percent_shared <- sum_shared_each/sum_total_each
comparison_shared <- data.frame(names(percent_shared),percent_shared)
names(comparison_shared)[1]<-"Sample"
head(comparison_shared)

### now want to do the opposite, just the unique ones
DNA_only <- filter(binary_wide_material_nometa, Intersect=="WaterDNA")
RNA_only <- filter(binary_wide_material_nometa, Intersect=="WaterRNA")
PIT_only <- filter(binary_wide_material_nometa, Intersect=="TrapDNA")
all_unique<- rbind(DNA_only, RNA_only, PIT_only)
all_unique$ASV1.ID <- row.names(all_unique) # add back the ASVs names
unique_names <- colsplit(all_unique$ASV1.ID, "_", names=c("ASV1.ID_real", "Taxon")) # split the ASV names and Taxon group
# filter the whole dataset by the unique ASVs
unique_ASVs <- filter(nometa, ASV1.ID %in% unique_names$ASV1.ID_real)
sum_unique_each <- colSums(unique_ASVs[,1:60]) # this gets the number of shared reads in all the samples
percent_unique_each <- sum_unique_each/sum_total_each
comparison_unique <- data.frame(names(percent_unique_each),percent_unique_each)
names(comparison_unique)[1]<-"Sample"
head(comparison_unique)

shared_by_two_ASVs <- filter(nometa, !(ASV1.ID %in% c(unique_names$ASV1.ID,shared_asvs_names$ASV1.ID_real))) 
sum_two_each <- colSums(shared_by_two_ASVs[,1:60]) # this gets the number of shared reads in all the samples
percent_two_each <- sum_two_each/sum_total_each
comparison_two <- data.frame(names(percent_two_each),percent_two_each)
names(comparison_two)[1]<-"Sample"
head(comparison_two)
### combine to for boxplot
unique_two <- join(comparison_unique,comparison_two, by="Sample", type="left")
all_groups <- join(unique_two, comparison_shared, by="Sample", type="left")
head(all_groups)
# check that they equal 1 
rowSums(all_groups[,2:4])
combo_melt <- melt(all_groups)
all_plotting <- join(combo_melt, key_all, by="Sample", type="left")
all_plotting$Sample.type <- factor(all_plotting$Sample.type, levels=c("Water RNA","Water DNA","Trap DNA") )

box_shared<-ggplot(all_plotting, aes(x=Sample.type, y=value, fill=variable))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.75, outlier.size = 1)+
  scale_fill_manual(values=c("#e0ecf4","#9ebcda","#8856a7"), labels=c("Unique", "Shared by two", "Shared by all"))+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  labs(title="", x="", y="Percent of reads")+
  guides(fill=guide_legend(title=""))+
  theme(axis.text.x = element_text(angle=45,hjust=.9,vjust=.9,size=8),
        axis.text.y=element_text(size=8),legend.position = "bottom",
        legend.text = element_text(size = 8),
        plot.margin = margin(0,0,0,0),
        text = element_text(size=8),
        legend.key.size = unit(0.25,'cm'))+
  theme(legend.title=element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))

#dunn test
dunn_test_asvs<-all_plotting %>%
  group_by(variable)%>%
  dunn_test(value~Material,p.adjust.method = "bonferroni")
#write.csv(dunn_test_asvs, "~/github/mesoscope2017/results/dunn_shared_boxplot.csv")

#### combining all plots ####
# plot the bar plot within the upsetr plot
upset.inset <-
  ggdraw() +
  draw_plot(upset_plot) +
  draw_plot(shared_barplot, x = .45, y = .3, width = .3, height = .5)
upset.inset

layout<-"
  AA
  AA
  BC"
wrap_plots(A=upset.inset, B=venn2, C=box_shared, design=layout) + plot_layout(widths = c())
#ggplot2::ggsave("~/github/mesoscope2017/figures/fig9.pdf", width=5, height = 6, units="in", dpi=300)

#### Plot Figure S3- nMDS of shared ASVs only ####
DNA_RNAonly <- filter(binary_wide_nometa_forotherplots, (Intersect=="Water RNA, Water DNA"|Intersect=="Water DNA, Water RNA"|Intersect=="Water RNA, Water DNA, Trap DNA")) # need to also include those found between all three
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

#First need to generate relative abundance
Norm.rel_nometa_shared <-decostand(WC_only_ASVs[,2:60], MARGIN=2, method = "total") #this does it by rows
colSums(Norm.rel_nometa_shared)# check! should all equal 1.
Norm.rel_nometa_shared$ASV.ID<-rownames(Norm.rel_nometa_shared)
melt_norm_nometa_shared<- melt(Norm.rel_nometa_shared)
vars_nometa<-colsplit(melt_norm_nometa_shared$variable, "_", c("Station","Depth","Material"))
new_df_nometa_shared<-data.frame(melt_norm_nometa_shared, vars_nometa)

# function for creating nmds
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

shared_all_nometa<- nmdsPoints(new_df_nometa_shared, key_all)
shared_nmds_all_nometa<-ggplot(data=shared_all_nometa, aes(x = MDS1, y = MDS2, fill=SLA,shape=Sample.type)) + 
  geom_point(color="black",size=2,position="jitter") + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = c(22,25,21))+ 
  labs(title="", x="NMDS1",y="NMDS2",size=10)+
  scale_y_continuous(breaks = seq(-3, 3, by = 1))+
  scale_x_continuous(breaks = seq(-4, 2, by = 1))+
  scale_fill_gradient2(low = "#0818A8",
                       mid = "light gray",
                       high = "#880808",
                       midpoint = 0,
                       limits=c(-28,28))+
  theme(text = element_text(size=8),plot.margin = unit(c(0, 0, 0, 0), "cm"))
shared_nmds_all <- shared_nmds_all_nometa +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Sample.type),linewidth=0.5) + scale_color_manual(values=col_material) #this adds in 95% confidence interval ellipses from ggplot
shared_nmds_all

dna_nometa_shared <- filter(new_df_nometa_shared, Material=="DNA")
dna_df_nometa_shared<-nmdsPoints(dna_nometa_shared, key_all)
dna_nmds_nometa_shared<-ggplot(data=dna_df_nometa_shared, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
  geom_point(color="black",size=2,position="jitter") + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = c(21,22,23,24))+ 
  labs(title="", x="NMDS1",y="NMDS2",size=10)+
  scale_y_continuous(breaks = seq(-3, 3, by = 1))+
  scale_x_continuous(breaks = seq(-4, 2, by = 1))+
  scale_fill_gradient2(low = "#0818A8",
                       mid = "light gray",
                       high = "#880808",
                       midpoint = 0,
                       limits=c(-28,28))+
  theme(text = element_text(size=8),plot.margin = unit(c(0, 0, 0, 0), "cm"))
dna_nmds_shared <- dna_nmds_nometa_shared +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth),linewidth=0.5, linetype="dashed") + scale_color_manual(values=col_depth) #this adds in 95% confidence interval ellipses from ggplot
dna_nmds_shared

rna_nometa_shared <- filter(new_df_nometa_shared, Material=="cDNA")
rna_df_nometa_shared<-nmdsPoints(rna_nometa_shared, key_all)
rna_nmds_nometa_shared<-ggplot(data=rna_df_nometa_shared, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
  geom_point(color="black",size=2,position="jitter") + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = c(21,22,23,24))+ 
  labs(title="", x="NMDS1",y="NMDS2",size=10)+
  scale_y_continuous(breaks = seq(-3, 3, by = 1))+
  scale_x_continuous(breaks = seq(-4, 2, by = 1))+
  scale_fill_gradient2(low = "#0818A8",
                       mid = "light gray",
                       high = "#880808",
                       midpoint = 0,
                       limits=c(-28,28))+
  theme(text = element_text(size=8),plot.margin = unit(c(0, 0, 0, 0), "cm"))
rna_nmds_shared <- rna_nmds_nometa_shared +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth),linewidth=0.5, linetype="dashed") + scale_color_manual(values=col_depth) #this adds in 95% confidence interval ellipses from ggplot
rna_nmds_shared
(shared_nmds_all / dna_nmds_shared / rna_nmds_shared ) + plot_layout(guides="collect") 
# ggplot2::ggsave("~/github/mesoscope2017/figures/S2_nmds.pdf", width=3.5, height = 6, units="in", dpi=300)

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
  subset(., (Taxa == "Alveolata-Ciliates"| Taxa == "Alveolata-Dinoflagellates" | Taxa=="Alveolata-Syndiniales"  | Taxa=="Alveolata-Dinoflagellates" )) %>%
  summarise(total_alv = sum(count))

sum_stram_wc <- asv_wc %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Stramenopila-Chrysophytes" | Taxa == "Stramenopila-Diatoms"|Taxa == "Stramenopila-MAST"|
               Taxa == "Stramenopila-Ochrophytes" | Taxa == "Stramenopila-Other" | Taxa == "Stramenopila-Pelagophytes")) %>%
  summarise(total_stram = sum(count))

sum_rhiz_wc <- asv_wc %>%
  group_by(Material, Depth, Station) %>%
  subset(., (Taxa == "Rhizaria-Acantharians"| Taxa == "Rhizaria-Cercozoans" | Taxa == "Rhizaria-Other" | Taxa == "Rhizaria-Polycystines")) %>%
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
