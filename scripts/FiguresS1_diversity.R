## Figure S1- Diversity statistics with vegan
## J. L. Beatty
## Last updated: 05/2024

# load packages
library(plyr)
library(tidyverse)
library(reshape2)
library(vegan)
library(patchwork)

# load data
asv_df <- read.csv("~/github/mesoscope2017/processed-data/MS_Clean_No1_Means_ALL_Norm_NewTax_0424.csv", header=T)
asv_df[1]=NULL
row.names(asv_df) = asv_df$ASV1.ID
nometa <- asv_df[asv_df$Level3 != "Metazoa",] # remove metazoan asvs

#metadata
key_all<-read.delim("~/github/mesoscope2017/raw-data/variables_all.txt")
key_all$Sample.type <- factor(as.character(key_all$Sample.type), levels=c( "Water RNA", "Water DNA", "Trap DNA"))

# for plotting
col_material<- c("Water RNA"="#7fcdbb","Water DNA"="#2c7fb8","Trap DNA"="#fc8d59")
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}
# calculate diversity metrics with vegan
shannon_nometa<-diversity(nometa[,1:60], index = "shannon",MARGIN = 2) # shannons diversity
invsimp_nometa<-diversity(nometa[,1:60], index = "invsimpson",MARGIN = 2) # inverse simpson's index
ASV_count_nometa<-colSums(nometa[,1:60]>0) # the number of reads per sample
ASV_number_nometa<-colSums(nometa[,1:60]) # the number of ASVs per sample
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
# group stats by material type
alpha_nometa_means <- alpha.df_nometa %>%
  group_by(Material, variable)%>%
  summarise(mean = mean(value), stdev=sd(value))
# write.csv(alpha_nometa_means, "~/github/mesoscope2017/processed-data/S1_means_sd.csv")

#plot the number of reads
seq_num <-ggplot(alpha.df_nometa[alpha.df_nometa$variable=="ASV_number_nometa",], aes(x=Sample.type, y=value, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material, name="Sample Type")+
  theme_bw()+
  labs(title="", x="", y="Number of reads")+
  scale_y_continuous(labels=scientific_10)+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  theme(text = element_text(size=8),
        axis.text.x=element_blank(),
        legend.position = "right")
seq_stat<-seq_num+ geom_pwc(
  aes(group = Sample.type), tip.length = 0,
  method = "dunn_test", label = "{p.adj.signif}",
  p.adjust.method = "bonferroni",
  hide.ns = TRUE)

# plot the number of ASVs
alpha<-ggplot(alpha.df_nometa[alpha.df_nometa$variable=="ASV_count_nometa",], aes(x=Sample.type, y=value, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material, name="Sample Type")+
  theme_bw()+
  labs(title="", x="", y="Number of ASVs")+
  scale_y_continuous(breaks=seq(0,3000,1000), limits=c(0,3000),labels = scientific_10) +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  theme(text = element_text(size=8),
        axis.text.x=element_blank(),
        legend.position = "right")
alpha_stat<-alpha + geom_pwc(
  aes(group = Sample.type), tip.length = 0,
  method = "dunn_test", label = "{p.adj.signif}",
  p.adjust.method = "bonferroni", 
  hide.ns = TRUE)


# plot shannon's index
shannons<-ggplot(alpha.df_nometa[alpha.df_nometa$variable=="shannon_nometa",], aes(x=Sample.type, y=value, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material, name="Sample Type")+
  theme_bw()+
  labs(title="", x="", y="Shannon's Index")+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  theme(text = element_text(size=8),
        axis.text.x=element_blank(),
        legend.position = "right")+
  scale_y_continuous(labels=scientific_10)
shannons_stat<-shannons+ geom_pwc(
  aes(group = Sample.type), tip.length = 0,
  method = "dunn_test", label = "{p.adj.signif}",
  p.adjust.method = "bonferroni", 
  hide.ns = TRUE)


# plot inverse simpson's
inv_simp<-ggplot(alpha.df_nometa[alpha.df_nometa$variable=="invsimp_nometa",], aes(x=Sample.type, y=value, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material, name="Sample Type")+
  theme_bw()+
  labs(title="", x="", y="Inverse Simpon's")+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  theme(text = element_text(size=8),
        axis.text.x=element_blank(),
        legend.position = "right")+
  scale_y_continuous(labels = scientific_10)
inv_simp_stat<-inv_simp + geom_pwc(
  aes(group = Sample.type), tip.length = 0,
  method = "dunn_test", label = "{p.adj.signif}",
  p.adjust.method = "bonferroni", 
  hide.ns = TRUE)


# combine plots
((seq_stat + alpha_stat) / (shannons_stat + inv_simp_stat)) + plot_layout(guides="collect") + plot_annotation(tag_levels="a")
#ggplot2::ggsave("~/github/mesoscope2017/figures/FigureS1_diversity.pdf",  width=5, height = 5, units="in", dpi=300)
