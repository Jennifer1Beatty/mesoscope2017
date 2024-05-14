## Figure 6- SPC treemaps
## J. L. Beatty
## Last updated: 05/2024

# load required packages
library(plyr)
library(tidyverse)
library(reshape2)
library(treemapify)

# set up variables for plotting
spc_key <- read.csv("~/github/mesoscope2017/raw-data/spc_key.csv",header=T)
spc_colors=c("Rhizaria"="#ca5670","Acantharia"="#ff9ab1", "Foraminifera"="#e33472","Radiolaria"="#ffa9da", "Rad_colony"="#ffd3ff",
             "Mesozooplankton"="#ab62c0","Appendicularia"="#d470de","Chaetognath"="#c47adb","Crustacean"="#9750ae","Gelatinous"="#bb76f6",
             "Diatom"="#72a555","Chain"="#4e8033","Circle"="#3da92d","Single"="#9ced86",
             "Trichodesmium"= "#ffb875", "Trichodesmium"="#ffb875")
spc_order=c("Diatom","Chain","Circle","Single",
            "Trichodesmium",
            "Rhizaria","Acantharia", "Foraminifera","Radiolaria", "Rad_colony",
            "Mesozooplankton","Appendicularia","Chaetognath","Crustacean","Gelatinous")
# load data
spc_100 <- read.csv("~/github/mesoscope2017/raw-data/MS_SPC_100.csv",header=T)
spc_100_m <- melt(filter(spc_100, Group != "Bad")) # remove the bad group
spc_100_m$Type = factor(spc_100_m$Type, levels=spc_order)

# calculate the relative abundance of each type, 
nobad_rel <-spc_100_m %>%
  group_by(variable)%>%
  mutate(relabund = value/sum(value)) %>%
  separate(variable,c("Sample", "Depth"), "_")
# add the metadata
spc_final <- join(nobad_rel, subset(spc_key, select=-c(Depth)), by="Sample", type="left")
head(spc_final)
# take the mean for each eddy type, excluding the sample that was an error
combo_rel <-subset(spc_final, Sample!="SPC.7042017") %>%
  group_by(Group,Type,Eddy)%>%
  summarise(mean_rel = mean(relabund), stdev=sd(relabund))

#### plot ####
ggplot(combo_rel, aes(area=mean_rel, fill=Type, label=Type, subgroup=Group))+
  geom_treemap()+
  geom_treemap_subgroup_border(size=0.75)+
  geom_treemap_text(colour = "white", place = "centre", reflow = T)+
  geom_treemap_text(aes(label = scales::percent(mean_rel, 1)), color = "white")+
  scale_fill_manual(values=spc_colors)+
  theme(text = element_text(size=8),plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.key.size = unit(0.25,'cm'),
        legend.margin=margin(0,0,0,0),
        strip.background = element_blank())+ 
  facet_wrap(~Eddy)
#ggplot2::ggsave("~/github/mesoscope2017/figures/fig9.pdf", width=4, height = 3.5, units="in", dpi=300)


#### stats for groups and types per eddy type ####
nobad_mean_group <- nobad_rel %>%
  group_by(variable, Group) %>%
  mutate(sum_group = sum(relabund)) %>%
  group_by(Group) %>%
  summarise(mean_group = mean(sum_group), sd_group=sd(sum_group))
nobad_mean_type <- nobad_rel %>%
  group_by(variable, Type) %>%
  mutate(sum_type = sum(relabund)) %>%
  group_by(Type) %>%
  summarise(mean_type = mean(sum_type), sd_type=sd(sum_type))


