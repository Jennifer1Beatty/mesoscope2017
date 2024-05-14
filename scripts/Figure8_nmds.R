## Figure 8- nMDS from molecular communities
## J. L. Beatty
## Last update: 05/2024

# load packages
library(plyr)
library(tidyverse)
library(reshape2)
library(vegan) 
library(patchwork)
# load dataframe and metadata
asv_df <- read.csv("~/github/mesoscope2017/processed-data/MS_Clean_No1_Means_ALL_Norm_NewTax_0424.csv", header=T)
asv_df[1]=NULL
row.names(asv_df) = asv_df$ASV1.ID
nometa <- asv_df[asv_df$Level3 != "Metazoa",] # remove metazoan asvs

key_all<-read.delim("~/github/mesoscope2017/raw-data/variables_all.txt")
key_all$Station<-factor(as.character(key_all$Station),levels=c("S4","S6","S8","S10","S12","S14",
                                                               "PIT12","PIT11","PIT10","PIT9","PIT8","PIT7","PIT6","PIT5","PIT4","PIT3","PIT2","PIT1"))
key_all$Material<-factor(key_all$Material,levels=c("RNA","DNA","PIT"))
key_all$Depth<-factor(as.character(key_all$Depth),levels=c("15m","DCM","PIT","175m","500m"))
key_all$Sample.type <- factor(as.character(key_all$Sample.type), levels=c( "Water RNA", "Water DNA", "Trap DNA"))

# for plotting
col_depth <-c("15m" = "#fee08b", "DCM" = "#7fbc41", "175m" = "#74add1", "500m" = "#542788")
col_material<- c("Water RNA"="#7fcdbb","Water DNA"="#2c7fb8","Trap DNA"="#fc8d59")

# plot nmds
# a function for creating nmds points using vegan
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
#First need to generate relative abundance
Norm.rel_nometa <-decostand(nometa[,2:60], MARGIN=2, method = "total") #this does it by rows
colSums(Norm.rel_nometa) # check! should all equal 1.
Norm.rel_nometa$ASV.ID<-rownames(Norm.rel_nometa)
melt_norm_nometa<- melt(Norm.rel_nometa)
vars_nometa<-colsplit(melt_norm_nometa$variable, "_", c("Station","Depth","Material"))
new_df_nometa<-data.frame(melt_norm_nometa, vars_nometa)

all_nometa<- nmdsPoints(new_df_nometa, key_all)
nmds_all_nometa<-ggplot(data=all_nometa, aes(x = MDS1, y = MDS2, fill=SLA,shape=Sample.type)) + 
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
nmds_all <- nmds_all_nometa +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Sample.type),linewidth=0.5) + scale_color_manual(values=col_material) #this adds in 95% confidence interval ellipses from ggplot
nmds_all

dna_nometa <- filter(new_df_nometa, Material=="DNA")
dna_df_nometa<-nmdsPoints(dna_nometa, key_all)
dna_nmds_nometa<-ggplot(data=dna_df_nometa, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
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
dna_nmds <- dna_nmds_nometa +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth),linewidth=0.5, linetype="dashed") + scale_color_manual(values=col_depth) #this adds in 95% confidence interval ellipses from ggplot
dna_nmds

rna_nometa <- filter(new_df_nometa, Material=="cDNA")
rna_df_nometa<-nmdsPoints(rna_nometa, key_all)
rna_nmds_nometa<-ggplot(data=rna_df_nometa, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
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
rna_nmds <- rna_nmds_nometa +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth),linewidth=0.5, linetype="dashed") + scale_color_manual(values=col_depth) #this adds in 95% confidence interval ellipses from ggplot
rna_nmds

(nmds_all / dna_nmds / rna_nmds) + plot_layout(guides="collect") # +  plot_annotation(tag_levels = "a")
# ggplot2::ggsave("~/github/mesoscope2017/figures/figure8.pdf", width=3.5, height = 6, units="in", dpi=300)
