##Jennifer Beatty
# 03/2023
# OTU Analysis for Mesoscope 2017 Cruise



# Load necessary packages

library(plyr)
library(tidyverse)
library(reshape2)
library(decontam)  # for removing possible contaminated ASVs, doesn't work with updated R
library(vegan)  # for statistical analysis
library(RColorBrewer)
library(edgeR)  # for statistical normalization, doesn't work with udpated R
library(ggpubr)
library(scico) # for scientifically approved colors
library(ggupset)  # for upset plots
library(stringi)  # for maniuplating strings

#### Input all the extra variables and factor them #####
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
key_all$Shape.mat<-as.numeric(key_all$Shape.mat)
key_all$Line.mat<-factor(as.numeric(key_all$Line.mat),levels=c("6", "1", "3"))
key_all$Station.density <- paste(key_all$Station, key_all$Potential_Density, sep="_")
key_all$Stat.depth.density <- paste(key_all$Station.density, key_all$Depth, sep="_")
summary(key_all)


#### Read in OTU table and taxonomy file, and combine ####
OTU_table<-read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/OTU_table_97.csv",header=TRUE)
rownames(OTU_table)= OTU_table$OTU.ID
head(OTU_table)
OTU_taxonomy <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/taxonomy_OTU_97_0323.csv", header=TRUE)
head(OTU_taxonomy)
OTU_df <- merge(x=OTU_table,y=OTU_taxonomy,by="OTU.ID",all.x=TRUE)
head(OTU_df)


######## 2. GET RID OF WEIRD COLUMNS ##########
Table_Cor<-OTU_df
#To remove the sample that is all unknowns
names(Table_Cor[83]) # S6_15m_DNA_B
Table_Cor[,83]=NULL
#To remove 0 column, "S14_500_cDNA_A"
sum(Table_Cor[60]) #To check if the sum is 0
names(Table_Cor)[60] #To check the name of the column- S14_500m_cDNA_A
Table_Cor[,60]=NULL #To remove the column
#To remove the messed up blank (Blank_C)
names(Table_Cor[4]) #- Blank_C
sum(Table_Cor[4])
Table_Cor[,4]=NULL
dim(Table_Cor) # 67973   111

########## 3. USE BLANKS TO REMOVE POSSIBLE CONTAMINATED SEQUENCES #########
#Need to separate datasets into "batches" that match with the specific blanks because the PITs were processed at a different time than the water column
#Water-Column samples 
names(Table_Cor[4:16]) #check where the pits are located
t_OTUs_WC <- Table_Cor %>%
  select(.,-c(1, 4:16, 110:111)) %>% #to remove the pit columns and other non-numeric columns
  t %>% # transpose
  as.matrix
dim(t_OTUs_WC) #95 rows and 67973 columns
#Need to create a vector True if blank, False if other
rownames(t_OTUs_WC)
blanks_WC <- c(TRUE,TRUE, FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
               FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
               FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)
contam_wc_vector <- isContaminant(t_OTUs_WC,neg=blanks_WC, detailed=FALSE) #this returns a logical vector where true is contimate
length(which(contam_wc_vector == TRUE)) #Shows how many are true aka contaminants We got 10
Contam_wc<-which(contam_wc_vector == TRUE) #Shows the location of which are true contaminants
write(Contam_wc, "~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/Contaminated_WC_OTU.txt")
#PIT samples
t_OTUs_PIT <- Table_Cor %>%
  select(.,c(4:16)) %>% # selecting only PIT columns
  t %>% # transpose
  as.matrix
dim(t_OTUs_PIT) #13 rows amd 67973 columns
#Need to create a vector when whether or not it's a blank sample
row.names(t_OTUs_PIT)
blanks_pits <- c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE)
contam_pits_vector <- isContaminant(t_OTUs_PIT,neg=blanks_pits, detailed=FALSE) #this returns a logical vector where true is contimate
length(which(contam_pits_vector == TRUE)) #Shows how many are true aka contaminants We got 9
contam_pit<-which(contam_pits_vector == TRUE) #Shows the location of which are true contaminants: 23445 25910 27070 32385 
write(contam_pit, "~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/Contaminated_PIT_OTU.txt")
#then add this to the document created by the wc code, and rename the file
#Let's remove the contaminants...From Gerid
#Now I want to remove them from a complete transposed data frame. 
#They should be the same row index as column index, because we didn't change the OTU orderever
#First, we need indexes instead of OTU.ID as rownames
Clean_table<-unrowname(Table_Cor)
#Next, I read in the contaminant vector. This must be read in with the scan fx.
contaminates<-scan("~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/Contaminated_ALL_OTU.txt") #This is the vector that contains the index for every contaminant 
# Filter the table based on the contaminate vector
Filtered_table<-Clean_table[!rownames(Clean_table) %in% contaminates, ]
#Verification steps. the number of rows in the Filtered table should equal nrows OTU_table - length of contaminants.
nrow(Table_Cor) - nrow(Filtered_table) # 19, so it worked!
#To add back the ASV IDs into the header as row names
rownames(Filtered_table) = Filtered_table$OTU.ID #Here I want to give the rownames the OTU IDs
head(Filtered_table) 
dim(Filtered_table) #  36371 rows and  111 columns, 108 columns of samples, first column is OTU names and last column is taxon list
#I'm going to remove the Blank samples because they have served their purpose
names(Filtered_table[c(2,3,16)])
Filtered_table[,c(2,3,16)]=NULL
dim(Filtered_table) # 67954 rows and 108 columns 
write.csv(Filtered_table,"~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/MS_Clean_OTU_0323.csv")


######## 03/03/23- Calculate ASV stats ###########
# clean_counts <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/MS_Clean_OTU_0323.csv", header=T) 
# clean_counts[1]= NULL
# rownames(clean_counts) = clean_counts$ASV.ID
# head(clean_counts)
# first need to take the average of duplicates
clean_mean <- Filtered_table %>%
  select(., -108) %>%  # remove confidence of ASV so as not to confuse
  melt
head(clean_mean)
var<-colsplit(clean_mean$variable, '_',c("Station","Depth","Type","Rep")) #Separate the variables in sample name
clean_mean_combo<-data.frame(clean_mean,var) #combine dataframes
sample_mean <- clean_mean_combo %>%
  dcast(., OTU.ID+Station+Depth+Type~Rep, value.var="value")  %>%  #recast data frame by reps
  mutate(mean = ((A+B)/2))  %>%  #Take the mean of rep A and B
  unite(.,"Sample",c("Station","Depth","Type")) %>%  #regroup the sample names
  dcast(.,OTU.ID~Sample, value.var="mean") #recast the data frame with the mean as the value
head(sample_mean)
dim(sample_mean) 

write.csv(sample_mean, "~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/MS_Clean_mean_0323_OTU.csv")
#downloaded it and then combined it with the numbers for the samples without blanks
all_counts<-read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/MS_Clean_mean_0323_OTU_ALL.csv", header=T)
rownames(all_counts)=all_counts$OTU.ID
head(all_counts)
all_counts.t<-t(all_counts)
# Number of sequences per sample
colsum<-apply(all_counts[2:61],2,sum); colsum
mean(colsum)
colsumm<-data.frame(all_counts.t[2:61,1],colsum)
colsumm[,1]=rownames(colsumm)
names(colsumm)=c("Sample","colsum")
# I want to be able to color and rename them 
colsum_df <- join(colsumm, key_all, by="Sample", type="left", match="first") %>%
  arrange(., Material, Station.sla, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))
col_mat <- unique(key_all$Color.mat)

all.plot<-ggplot(colsum_df,aes(x=Sample,y=colsum, fill=Material))+ geom_bar(stat="identity",color="black")+
  labs(title="OTUs Sequence Depth", x="",y="Total Reads")+theme_bw()+
  theme(axis.text.x = element_blank(),axis.text.y=element_text(size=10,color = "black"), strip.background = element_blank())+ theme(legend.title = element_blank())+
  scale_fill_manual(values=c("#7fcdbb","#2c7fb8","#fc8d59"))
all.plot + scale_y_continuous(trans = 'log10')
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/All_reads_mean_0323_OTU.pdf", width=10, height = 10)

# Reports total number of unique OTUs
length(unique(all_counts$OTU.ID)) #67954

# ASV counts and information:
counts_only<-all_counts[2:61]
seq_total<-apply(counts_only,2,sum)
min(seq_total); max(seq_total); mean(seq_total)
OTU_count<-colSums(counts_only>0); OTU_count
min(OTU_count); max(OTU_count); mean(OTU_count)
OTU_single<-colSums(counts_only==1) #by sample
OTU_double<-colSums(counts_only==2) #by sample
OTU_true<-colSums(counts_only>2)
#
#Compile sample information
sample_info<-data.frame(seq_total,OTU_count,OTU_single,OTU_double,OTU_true);sample_info
#Visual representation of ASVs:

sample_info$Sample<-row.names(sample_info)
head(sample_info)
counts.melt<-melt(sample_info[c(6,1:4)])
counts.plot <- join(counts.melt, key_all, by="Sample", type="left", match="first") %>%
  arrange(., Material, Station.sla, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))

head(counts.plot)
means_type <- tapply(counts.plot$value, list(counts.plot$Material, counts.plot$variable), mean)
write.csv(means_type, "~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/OTU_means_0323.csv")

OTUs<-c("OTU_count", "OTU_single", "OTU_double")

All_bar_stats<- ggplot(counts.plot, aes(x=Sample, y=value, fill=variable))+
  geom_bar(stat="identity",position="stack",color="black")+
  labs(title="OTU Counts", x="",y="Total OTUs")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=1,size=12, color = "black"),axis.text.y=element_text(size=12,color = "black"), strip.background = element_blank())+ theme(legend.title = element_blank())
All_bar_stats %+% subset(counts.plot, variable %in% OTUs) + scale_fill_manual("",values=c("#e41a1c","#fee08b","#4393c3"),labels = c("OTUs > 2 seqs", "Singletons", "Doubletons"))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/All_OTU_Stats_0323.pdf", width=10, height = 10)

## could break this up if need be, but the main point is that they are mostly more than doubletons
DNA<- All_bar_stats %+% subset(counts.plot, (variable %in% OTUs & Material=="DNA")) + scale_fill_manual("",values=c("#e41a1c","#fee08b","#4393c3"),labels = c("OTUs > 2 seqs", "Singletons", "Doubletons"))
RNA<- All_bar_stats %+% subset(counts.plot, (variable %in% OTUs & Material=="RNA")) + scale_fill_manual("",values=c("#e41a1c","#fee08b","#4393c3"),labels = c("OTUs > 2 seqs", "Singletons", "Doubletons"))
PIT<- All_bar_stats %+% subset(counts.plot, (variable %in% OTUs & Material=="PIT")) + scale_fill_manual("",values=c("#e41a1c","#fee08b","#4393c3"),labels = c("OTUs > 2 seqs", "Singletons", "Doubletons"))

ggarrange(DNA, RNA, PIT, common.legend=TRUE, legend="right")
# may want to make the y axis the same for them all
# may want to change the titles to the material type, and edit the sample labels


############ 4. Remove global singletons ##########
##Filter out ASVs with only 1 sequence in the whole dataset (global singletons)
OTU_sum<-apply(Filtered_table[2:106],1,sum, na.rm =TRUE) #remove global singletons by rows- summing total of reads an ASV is found across all samples
OTU_no1 = Filtered_table[ OTU_sum>1, ]  #count.no1 = OTU table without global singletons
removed = dim(Filtered_table)[1] - dim(OTU_no1)[1] #Outputs the number of OTUss (total) lost in this step, 30,084
names(OTU_no1)[1] = "OTU.ID"
dim(OTU_no1) #Check dimensions- 37870 rows and 108 columns
write.csv(OTU_no1, "~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/MS_Clean_No1_0323_OTU.csv")

############ 5. Take the mean of duplicate samples #########
No1_m <- OTU_no1 %>%
  select(., -108) %>%  # remove confidence of ASV so as not to confuse
  melt
head(No1_m)
var<-colsplit(No1_m$variable, '_',c("Station","Depth","Type","Rep")) #Separate the variables in sample name
No1_m_var<-data.frame(No1_m,var) #combine dataframes
sample_mean <- No1_m_var %>%
  dcast(., OTU.ID+Station+Depth+Type~Rep) %>%  #recast data frame by reps
  mutate(mean = (A+B)/2) %>%  #Take the mean of rep A and B
  unite(.,"Sample",c("Station","Depth","Type")) %>%  #regroup the sample names
  dcast(.,OTU.ID~Sample, value.var="mean") #recast the data frame with the mean as the value
head(sample_mean)
dim(sample_mean) # 37870    61
write.csv(sample_mean, "~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/MS_Clean_No1_Means_0323_OTU.csv")

###### 6. Now that we have a clean, mean, no global singletons table, we can rarefy to see sequencing depth ######
#I took the means calculated above, and manually added back the rows that didn't have duplicates
OTU_mean<-read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/MS_Clean_No1_Means_0323_OTU_All.csv",header=T)
rownames(OTU_mean)=OTU_mean$OTU.ID
dim(OTU_mean) #34943 rows and 63 columns 
t_mean <- OTU_mean[2:61] %>%  # selecting only numeric columns
  as.matrix %>%
  t %>%  # transpose
  floor  # convert to whole numbers

##### 03/3/23 Trying a new method to visualization the rarefaction curve #######
# from Patt Schloss
options(scipen = 999) # this prevents scientific notation
rarecurve_data <- rarecurve(t_mean, step =20)
col_mat <- unique(key_all$Color.mat)
# c("#7fcdbb","#2c7fb8","#fc8d59")

# there are NA's in the materials???
map_dfr(rarecurve_data, bind_rows) %>%
  bind_cols(Sample = rownames(t_mean),.) %>%
  pivot_longer(-Sample) %>%
  drop_na() %>%
  mutate(n_seqs= as.numeric(str_replace(name, "N", "")))%>%
  select(-name) %>%
  join(.,key_all, by= "Sample") %>%
  ggplot(aes(x=n_seqs, y=value, group=Sample)) + 
  geom_line(aes(color=Material)) +
  scale_color_manual(values=c("#7fcdbb","#2c7fb8","#fc8d59"))+
  labs(title="Rarefaction Curve of OTUs", x="Number of sequences",y="Number of OTUs")+
  theme_bw()
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/Rarefaction_0323_OTU.pdf", width=6, height = 6)

############# 7. Normalize ##################
#From Gerid
ListDGE = DGEList(OTU_mean[,2:61]) 
ListDGE
ListDGE = calcNormFactors(ListDGE, method = "TMM")
ListDGE
TMMNorm_OTU = cpm(ListDGE)
head(TMMNorm_OTU)
TMMNorm_OTU = as.data.frame(TMMNorm_OTU)
TMMNorm_OTU$OTU.ID = row.names(TMMNorm_OTU)
Joined<-join(TMMNorm_OTU, OTU_mean[c(1)], by="OTU.ID", type="left", match="first") 
dim(Joined) # 37870 rows and 61 columns
write.csv(Joined, "~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/MS_Clean_No1_Means_Norm_0223.csv")
############# 8. Rename the taxa from PR2, into nicer names ##########
##### Function Updated 11/17/22 to reflect the PR2, V14
new_tax <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/taxonomy_OTU_97_0323.csv", header=T)
names(new_tax)[1] = "OTU.ID"

pr2_rename_taxa_w2<-function(df){
  library(reshape2)
  split<-colsplit(df$Taxon, ";", c("Level1","Level2","Level3","Level4","Level5","Level6", "Level7","Level8"))
  split[ is.na(split) ] = "XXX"
  split[ split == "" ] = "XXX"
  split$Taxa<-"Other/unknown"
  split$Taxa[split$Level1 == "Eukaryota"]="Other/unknown"
  split$Taxa[split$Level2=="Alveolata"]="Alveolates-Other"
  split$Taxa[split$Level2=="Opisthokonta"]="Opisthokonts-Other"
  split$Taxa[split$Level2=="Rhizaria"]="Rhizaria-Other"
  split$Taxa[split$Level2=="Stramenopiles"]="Stramenopiles-Other"
  split$Taxa[split$Level2=="Hacrobia"]="Hacrobia-Other"
  split$Taxa[split$Level2=="Archaeplastida"]="Archaeplastids-Other"
  split$Taxa[split$Level2=="Amoebozoa"]="Amoebozoa"
  split$Taxa[split$Level2=="Excavata"]="Excavates"
  split$Taxa[split$Level2=="Eukaryota_X"]="Other/unknown"
  split$Taxa[split$Level2=="Apusozoa"]="Other/unknown"
  split$Taxa[split$Level3=="Dinoflagellata"]="Alveolates-Dinoflagellates"
  split$Taxa[split$Level3=="Metazoa"]="Opisthokont-Metazoa"
  split$Taxa[split$Level3=="Radiolaria"]="Rhizaria-Radiolaria"
  split$Taxa[split$Level3=="Ciliophora"]="Alveolates-Ciliates"
  split$Taxa[split$Level3=="Ochrophyta"]="Stramenopiles-Ochrophyta"
  split$Taxa[split$Level3=="Haptophyta"]="Hacrobia-Haptophytes"
  split$Taxa[split$Level3=="Opisthokonta_X"]="Opisthokonts-Other"
  split$Taxa[split$Level3=="Choanoflagellida"]="Opisthokont-Choanoflagellida"
  split$Taxa[split$Level3=="Chlorophyta"]="Archaeplastids-Chlorophytes"
  split$Taxa[split$Level3=="Cercozoa"]="Rhizaria-Cercozoa"
  split$Taxa[split$Level3=="Cryptophyta"]="Hacrobia-Cryptophytes"
  split$Taxa[split$Level3=="Fungi"]="Opisthokont-Fungi"
  split$Taxa[split$Level3=="Discoba"]="Excavates-Discoba"
  split$Taxa[split$Level3=="Foraminifera"]="Rhizaria-Foraminifera"
  split$Taxa[split$Level4=="Dinophyceae"]="Alveolates-Dinophyceae"
  split$Taxa[split$Level4=="Syndiniales"]="Alveolates-Syndiniales"
  split$Taxa[split$Level4=="Polycystinea"]="Rhizaria-Polycystines"
  split$Taxa[split$Level4=="Acantharea"]="Rhizaria-Acantharia"
  split$Taxa[split$Level4=="Pelagophyceae"]="Stramenopiles-Pelagophytes"
  split$Taxa[split$Level4=="RAD-A"]="Rhizaria-RAD (A,B,C)"
  split$Taxa[split$Level4=="RAD-B"]="Rhizaria-RAD (A,B,C)"
  split$Taxa[split$Level4=="RAD-C"]="Rhizaria-RAD (A,B,C)"
  split$Taxa[split$Level4=="Bacillariophyta"]="Stramenopiles-Diatoms"
  split$Taxa[split$Level4=="Chrysophyceae"]="Stramenopiles-Chrysophytes"
  split$Taxa[split$Level4=="MAST-1"]="Stramenopiles-MAST"
  split$Taxa[split$Level4=="MAST-2"]="Stramenopiles-MAST"
  split$Taxa[split$Level4=="MAST-3"]="Stramenopiles-MAST"
  split$Taxa[split$Level4=="MAST-4"]="Stramenopiles-MAST"
  split$Taxa[split$Level4=="MAST-7"]="Stramenopiles-MAST"
  split$Taxa[split$Level4=="MAST-8"]="Stramenopiles-MAST"
  split$Taxa[split$Level4=="MAST-9"]="Stramenopiles-MAST"
  split$Taxa[split$Level4=="MAST-10"]="Stramenopiles-MAST"
  split$Taxa[split$Level4=="MAST-11"]="Stramenopiles-MAST"
  split$Taxa[split$Level4=="MAST-12"]="Stramenopiles-MAST"
  split$Taxa[split$Level4=="MAST-23"]="Stramenopiles-MAST"
  split$Taxa[split$Level4=="MAST-25"]="Stramenopiles-MAST"
  split$Taxa2<-"XXX"
  four<-c("Alveolates-Ciliates","Archaeplastids-Chlorophytes")
  split$Taxa2<-with(split, ifelse(Taxa %in% four, Level4, Taxa2)) #Take taxa named in "four" and place in the level4 name, this is the next level of tax resolution I would like to show.
  five<-c("Alveolates-Syndiniales","Stramenopiles-MAST")
  split$Taxa2<-with(split, ifelse(Taxa %in% five, Level5, Taxa2))
  two<-c("Alveolates-Other", "Other/unknown")
  split$Taxa2<-with(split, ifelse(Taxa %in% two, Level3, Taxa2))
  six<-c("Hacrobia-Haptophytes", "Hacrobia-Crytophytes")
  split$Taxa2<-with(split, ifelse(Taxa %in% six, Level6, Taxa2))
  four<-c("Opisthokont-Metazoa", "Opisthokonts-Other")
  split$Taxa2<-with(split, ifelse(Taxa %in% four, Level4, Taxa2))
  seven<-c("Alveolates-Dinoflagellates","Archaeplastids-Other","Excavates","Rhizaria-Acantharia","Rhizaria-Cercozoa","Rhizaria-Polycystines","Rhizaria-RAD (A,B,C)","Stramenopiles-Diatoms","Stramenopiles-Pelagophytes","Stramenopiles-Chrysophytes","Stramenopiles-Other")
  split$Taxa2<-with(split, ifelse(Taxa %in% seven, Level7, Taxa2))
  # head(split)
  # unique(split$Taxa2)
  split$Taxa2<-gsub("_XXX", "", split$Taxa2)
  split$Taxa2<-gsub("_XX", "", split$Taxa2)
  split$Taxa2<-gsub("_X", "", split$Taxa2)
  #If taxa2 = "XXX" replace with Taxa, and other
  split$Taxa2<-gsub("XXX","Other/unknown",split$Taxa2)
  # want to make another column to match Isha naming
  split$TaxaPlot<-"XXX"
  okay<-c("Alveolates-Ciliates","Alveolates-Dinophyceae","Alveolates-Syndiniales","Archaeplastids-Chlorophytes",
          "Excavates-Discoba", "Hacrobia-Cryptophytes","Hacrobia-Haptophytes",
          "Opisthokont-Choanoflagellida","Opisthokont-Metazoa","Other/unknown",
          "Rhizaria-Acantharia","Rhizaria-Cercozoa","Rhizaria-Polycystines","Rhizaria-Other",
          "Stramenopiles-Chrysophytes", "Stramenopiles-Diatoms","Stramenopiles-MAST","Stramenopiles-Ochrophyta",
          "Stramenopiles-Pelagophytes","Stramenopiles-Other")
  split$TaxaPlot<-with(split, ifelse(Taxa %in% okay, Taxa, TaxaPlot)) #Take taxa named in "okay" and place in the Taxa
  split$TaxaPlot[split$Taxa=="Rhizaria-RAD (A,B,C)"]="Rhizaria-Other"
  # head(split)
  # unique(split$TaxaPlot)
  split$TaxaPlot<-gsub("_XXX", "", split$TaxaPlot)
  split$TaxaPlot<-gsub("_XX", "", split$TaxaPlot)
  split$TaxaPlot<-gsub("_X", "", split$TaxaPlot)
  #If taxaPlot = "XXX" replace with Taxa, and other
  split$TaxaPlot<-gsub("XXX","Other/unknown",split$TaxaPlot)
  head(split)
  return(split)
} 

table <- join(Joined, new_tax, by="OTU.ID")
NewNameTax = pr2_rename_taxa_w2(table)
combined_newtax = data.frame(Joined, NewNameTax)
combined_newtax[1:5,]
write.csv(combined_newtax, "~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/MS_Clean_No1_Means_Norm_NewTax_0323.csv")
###### START HERE FOR OTU TABLE WITH UPDATED TAXONOMY ##########
otu_df <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/MS_Clean_No1_Means_Norm_NewTax_0323.csv", header=T)
otu_df[1]=NULL
row.names(otu_df) = otu_df$OTU.ID

#### To add in the ASV data
asv_df <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Output/MS_Clean_No1_Means_Complete_Norm_NewTax_1122.csv", header=T)
asv_df[1]=NULL
row.names(asv_df) = asv_df$ASV1.ID

######## 9. Shannon's Diveristy #########
col_material<- c("#7fcdbb","#2c7fb8","#fc8d59")

# Needs to transposeee
shannons_df_otu<- otu_df[,2:61] %>%
  t(.)%>%
  diversity(., index = "shannon",MARGIN = 1) %>%
  melt() %>%
  mutate(Sample = rownames(.))%>%
  join(.,key_all, by="Sample") %>%
  arrange(., Material, Station.sla, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))
names(shannons_df_otu)[1]="H"
rownames(shannons_df_otu)=shannons_df_otu$Sample
head(shannons_df_otu)

ggplot(shannons_df_otu, aes(x = Sample, y = H,shape=Material,fill=SLA)) + 
  geom_point(color="black",size=6) + 
  theme_bw() + 
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_fill_gradientn(colours=c("blue", "red"))+
  scale_shape_manual(values = c(21,25,22)) + 
  labs(title="Shannon's Diversity Index- Metazoa Included-OTUs", x="",y="Shannon's Diversity Index")+
  theme(legend.position="right",axis.text.x = element_text(angle=90, hjust=1,vjust=1,size=10,color="black"))

ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/shannons_All_0323_OTU.pdf", width=8, height = 8)

# Now without metazoa
shannons_df_nometa_otu<- otu_df %>%
  filter(.,Taxa!="Opisthokont-Metazoa")%>%
  .[,2:61]%>%
  t(.)%>%
  diversity(., index = "shannon",MARGIN = 1) %>%
  melt() %>%
  mutate(Sample = rownames(.))%>%
  join(.,key_all, by="Sample") %>%
  arrange(., Material, Station.sla, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))
names(shannons_df_nometa_otu)[1]="H"
rownames(shannons_df_nometa_otu)=shannons_df_nometa_otu$Sample
head(shannons_df_nometa_otu)

ggplot(shannons_df_nometa_otu, aes(x = Sample, y = H,shape=Material,fill=SLA)) + 
  geom_point(color="black",size=6) + 
  theme_bw() + 
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_fill_gradientn(colours=c("blue", "red"))+
  scale_shape_manual(values = c(21,25,22)) + 
  labs(title="Shannon's Diversity Index- No Metazoa-OTUs", x="",y="Shannon's Diversity Index")+
  theme(legend.position="right",axis.text.x = element_text(angle=90, hjust=1,vjust=1,size=10,color="black"))

ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/shannons_nometa_0323_OTU.pdf", width=8, height = 8)

####### 10. Taxa Bar Plots #######
# ordered by Station and SLA
data.m_otu<-melt(otu_df) #melt
head(data.m_otu)
data.agg_otu<-aggregate(data.m_otu$value, by=list(Taxa=data.m_otu$TaxaPlot,Sample=data.m_otu$variable),sum) #sum sequences by taxonomic group
data.agg_otu$Sample<-factor(data.agg_otu$Sample,levels=names(otu_df[2:61]))
tax_order<-c("Alveolates-Ciliates","Alveolates-Dinophyceae","Alveolates-Syndiniales", "Archaeplastids-Chlorophytes",
             "Excavates-Discoba", "Hacrobia-Cryptophytes","Hacrobia-Haptophytes",
             "Opisthokont-Choanoflagellida","Opisthokonts-Metazoa","Other/unknown",
             "Rhizaria-Acantharia","Rhizaria-Cercozoa","Rhizaria-Polycystines", "Rhizaria-Other",
             "Stramenopiles-Chrysophytes","Stramenopiles-Diatoms","Stramenopiles-MAST","Stramenopiles-Ochrophyta","Stramenopiles-Pelagophytes","Stramenopiles-Other")
tax_color_phyla2<-c('firebrick4','indianred1','tomato3','forestgreen','yellowgreen','darkblue',
                    'lightblue','moccasin','gold1','grey',"#C291A4",'mediumvioletred','magenta','lightpink',"#DECDBE",'#DDAD4B','tan2','tan3','tan4',"#5C4033")
names(tax_color_phyla2)<-tax_order
bar_plot_df_otu <- join(data.agg_otu, key_all, by="Sample", type="left", match="first") %>%
  arrange(., Material, Station.sla, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))

#Bar plot of community composition ordered by Material, Station and Depth
bars_otu<-ggplot(bar_plot_df_otu, aes(y=x, fill=tax, x=Station.sla))+
  geom_bar(position = "fill",stat="identity",color="black",aes(fill=Taxa))+
  scale_fill_manual(values=tax_color_phyla2)+
  theme_bw()+
  labs(title="Relative abundance of reads-OTUs", x="Sample",y="")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5, vjust = .5))
WC <- c("DNA","RNA")
bars_otu %+% subset(bar_plot_df_otu, Material %in% WC) + coord_flip() +facet_wrap(.~ Material + Depth, scales = "free_y", ncol=2, dir="v")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/Taxa_All_Meta_WC_DepthSLA_OTU.pdf", width=17, height = 10)
bars_otu %+% coord_flip() +facet_wrap(.~ Material + Depth, scales = "free_y", ncol=3, nrow=4, dir="v")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/Taxa_All_Meta_DepthSLA_OTU.pdf", width=17, height = 10)


#now with no meta
nometa_otu <-filter(otu_df, Taxa!="Opisthokont-Metazoa")
data_nometa_otu<-melt(nometa_otu) #melt
head(data_nometa_otu)
data.agg_nometa_otu<-aggregate(data_nometa_otu$value, by=list(Taxa=data_nometa_otu$TaxaPlot,Sample=data_nometa_otu$variable),sum) #sum sequences by taxonomic group
data.agg_nometa_otu$Sample<-factor(data.agg_nometa_otu$Sample,levels=names(nometa_otu[2:61]))
tax_color_phyla2_nometa_otu<-c('firebrick4','indianred1','tomato3','forestgreen','yellowgreen','darkblue',
                           'lightblue','moccasin','grey',"#C291A4",'mediumvioletred','magenta','lightpink',"#DECDBE",'#DDAD4B','tan2','tan3','tan4',"#5C4033")

bar_plot_df_nometa_otu <- join(data.agg_nometa_otu, key_all, by="Sample", type="left", match="first") %>%
  arrange(., Material, Station.sla, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))

#Bar plot of community composition ordered by Material, Station and Depth
bars_nometa_otu<-ggplot(bar_plot_df_otu, aes(y=x, fill=tax, x=Station.sla))+
  geom_bar(position = "fill",stat="identity",color="black",aes(fill=Taxa))+
  scale_fill_manual(values=tax_color_phyla2_nometa_otu)+
  theme_bw()+
  labs(title="Relative abundance of reads-OTUs", x="Sample",y="")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5, vjust = .5))

bars_nometa_otu %+% subset(bar_plot_df_nometa_otu, Material %in% WC) + coord_flip() +facet_wrap(.~ Material + Depth, scales = "free_y", ncol=2, dir="v")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/Taxa_NOMeta_0323_WC_DepthSLA_OTU.pdf", width=17, height = 10)
bars_nometa_otu %+% subset(bar_plot_df_nometa_otu, !(Material %in% WC)) + coord_flip() +facet_wrap(.~ Material + Depth, scales = "free_y", ncol=3, nrow=4, dir="v")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/Taxa_NOMeta_0323_DepthSLA_OTU_PIT.pdf", width=10, height = 10)


######### 11. NMDS of community #######
Norm.rel_otu <-decostand(otu_df[,2:61], MARGIN=2, method = "total") #this does it by rows
colSums(Norm.rel_otu) # check! should all equal 1.
Norm.rel_otu$OTU.ID<-rownames(Norm.rel_otu)
melt_norm_otu<- melt(Norm.rel_otu)
vars1<-colsplit(melt_norm_otu$variable, "_", c("Station","Depth","Material"))
new_df_otu<-data.frame(melt_norm_otu, vars1)

nmdsPoints_otu <- function(df,key_df){
  new_df<-dcast(df, variable~OTU.ID, fill=0) #restructure dataframe to only include ASV1.ID and samples
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

#all with meta
all_otu<- nmdsPoints_otu(new_df_otu, key_all)
nmds_all_otu<-ggplot(data=all_otu, aes(x = MDS1, y = MDS2, fill=Loc,shape=Material)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = c(21,25,22))+ scale_fill_manual(values=(levels(all_otu$Color.sla)))+
  labs(title="NMDS-Bray-Curtis-OTUs", x="NMDS1",y="NMDS2",size=12)+
  guides(fill = guide_legend(override.aes=list(color=levels(all_otu$Color.sla))))+
  theme(legend.title = element_blank())
nmds_all_otu +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Material)) #this adds in 95% confidence interval ellipses from ggplot
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/meso_nmds_all_OTU.pdf", width=6, height = 6)

# all without meta
nometa_otu <-filter(otu_df, Taxa!="Opisthokont-Metazoa")
Norm.rel_nometa_otu <-decostand(nometa_otu[,2:61], MARGIN=2, method = "total") #this does it by rows
colSums(Norm.rel_nometa_otu) # check! should all equal 1.
Norm.rel_nometa_otu$OTU.ID<-rownames(Norm.rel_nometa_otu)
melt_norm_nometa_otu<- melt(Norm.rel_nometa_otu)
vars_nometa_otu<-colsplit(melt_norm_nometa_otu$variable, "_", c("Station","Depth","Material"))
new_df_nometa_otu<-data.frame(melt_norm_nometa_otu, vars_nometa_otu)

all_nometa_otu<- nmdsPoints_otu(new_df_nometa_otu, key_all)
nmds_all_nometa_otu<-ggplot(data=all_nometa_otu, aes(x = MDS1, y = MDS2, fill=Loc,shape=Material)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = c(21,25,22))+ scale_fill_manual(values=(levels(all_nometa_otu$Color.sla)))+
  labs(title="NMDS-Bray-Curtis- No Metazoa-OTUs", x="NMDS1",y="NMDS2",size=12)+
  guides(fill = guide_legend(override.aes=list(color=levels(all_nometa_otu$Color.sla))))+
  theme(legend.title = element_blank())
nmds_all_nometa_otu +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Material)) #this adds in 95% confidence interval ellipses from ggplot
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/meso_nmds_all_nometa_OTU.pdf", width=6, height = 6)

col_material<- c("#7fcdbb","#2c7fb8","#fc8d59")

all_nometa_otu<- nmdsPoints_otu(new_df_nometa_otu, key_all)
all_nometa_otu$Material<- factor(all_nometa_otu$Material, levels=c("DNA","RNA","PIT"))
nmds_all_nometa_otu_2<-ggplot(data=all_nometa_otu, aes(x = MDS1, y = MDS2, fill=SLA,shape=Sample.type)) + 
  geom_point(color="black",size=5,position="jitter") + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = unique(all_nometa_otu$Shape.mat), name="Sample Type")+ 
  labs(title="", x="NMDS1",y="NMDS2",size=12)+
  scale_y_continuous(breaks = seq(-3, 3, by = 1), trans="reverse")+
  scale_x_continuous(breaks = seq(-4, 2, by = 1))+
  scale_fill_gradient2(low = "#0818A8",
                       mid = "white",
                       high = "#880808",
                       midpoint = 0,
                       limits=c(-28,28))

nmds_all_nometa_otu_2 +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Sample.type),linewidth=1.5) + scale_color_manual(values=col_material, name="Sample Type") #this adds in 95% confidence interval ellipses from ggplot
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/meso_nmds_all_nometa_OTU_sla_depth.pdf", width=6, height = 6)


# no meta, dna only
dna_nometa_otu <- filter(new_df_nometa_otu, Material=="DNA")
dna_df_nometa_otu<-nmdsPoints_otu(dna_nometa_otu, key_all)

nmds_dna_nometa_otu<-ggplot(data=dna_df_nometa_otu, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = dna_df_nometa_otu$Shape.depth)+ # scale_fill_manual(values=levels(dna_df_nometa$Color.sla))+
  labs(title="NMDS-Bray-Curtis: DNA-No Metazoa-OTUs", x="NMDS1",y="NMDS2",size=12)+
  # scico::scale_fill_scico(palette = "batlow") +
  scale_fill_gradientn(colours=c("blue", "red"))+
  # guides(fill = guide_legend(override.aes=list(color=levels(dna_df_nometa$Color.sla))))+
  theme(legend.title = element_blank())
nmds_dna_nometa_otu + stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/meso_nmds_dna_nometa_sla0323_OTU.pdf", width=8, height = 8)


# no meta, rna only
rna_nometa_otu <- filter(new_df_nometa_otu, Material=="cDNA")
rna_df_nometa_otu<-nmdsPoints_otu(rna_nometa_otu, key_all)

nmds_rna_nometa_otu<-ggplot(data=rna_df_nometa_otu, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_df_nometa_otu$Shape.depth)+ # scale_fill_manual(values=levels(dna_df_nometa$Color.sla))+
  labs(title="NMDS-Bray-Curtis: RNA-No Metazoa-OTUs", x="NMDS1",y="NMDS2",size=12)+
  # scico::scale_fill_scico(palette = "batlow") +
  scale_fill_gradientn(colours=c("blue", "red"))+
  # guides(fill = guide_legend(override.aes=list(color=levels(dna_df_nometa$Color.sla))))+
  theme(legend.title = element_blank())
nmds_rna_nometa_otu + stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/OTUs/meso_nmds_rna_nometa_sla0323_OTU.pdf", width=8, height = 8)


