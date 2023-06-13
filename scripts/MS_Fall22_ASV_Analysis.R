## Analysis of the MESO-SCOPE 2017 Cruise 18S data
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
library(patchwork)

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
key_all$Sample.type <- factor(as.character(key_all$Sample.type), levels=c("Water column DNA", "Water column RNA", "Trap DNA"))
summary(key_all)
col_material<- c("#7fcdbb","#2c7fb8","#fc8d59")
#### Read in ASV table and taxonomy file, and combine ####
ASV_table<-read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/allmeso_table.csv",header=TRUE)
rownames(ASV_table)= ASV_table$ASV.ID
head(ASV_table)
ASV_taxonomy <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/taxonomy_asv.csv", header=TRUE)
head(ASV_taxonomy)
ASV_df <- merge(x=ASV_table,y=ASV_taxonomy,by="ASV.ID",all.x=TRUE)
head(ASV_df)


######## 2. GET RID OF WEIRD COLUMNS ##########
Table_Cor<-ASV_df
#To remove the sample that is all unknowns
names(Table_Cor[112]) # S6_15m_DNA_B
Table_Cor[,112]=NULL
#To remove 0 column, "S14_500_cDNA_A"
sum(Table_Cor[101]) #To check if the sum is 0
names(Table_Cor)[101] #To check the name of the column- S14_500m_cDNA_A
Table_Cor[,101]=NULL #To remove the column
#To remove the messed up blank (Blank_C)
names(Table_Cor[82]) #- Blank_C
Table_Cor[,82]=NULL
dim(Table_Cor) #36381   111

########## 3. USE BLANKS TO REMOVE POSSIBLE CONTAMINATED SEQUENCES #########
#Need to separate datasets into "batches" that match with the specific blanks because the PITs were processed at a different time than the water column
#Water-Column samples 
names(Table_Cor[82:94]) #check where the pits are located
t_ASVs_WC <- Table_Cor %>%
              select(.,-c(1, 82:94, 110:111)) %>% #to remove the pit columns and other non-numeric columns
              t %>% # transpose
              as.matrix
dim(t_ASVs_WC) #95 rows and 36381 columns
#Need to create a vector True if blank, False if other
rownames(t_ASVs_WC)
blanks_WC <- c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
              TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
              FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)
contam_wc_vector <- isContaminant(t_ASVs_WC,neg=blanks_WC, detailed=FALSE) #this returns a logical vector where true is contimate
length(which(contam_wc_vector == TRUE)) #Shows how many are true aka contaminants We got 6
Contam_wc<-which(contam_wc_vector == TRUE) #Shows the location of which are true contaminants
write(Contam_wc, "~/Desktop/Chapter1/MESOSCOPE2017/Output/Contaminated_WC_ASV.txt")
#PIT samples
t_ASVs_PIT <- Table_Cor %>%
            select(.,c(82:94)) %>% # selecting only PIT columns
            t %>% # transpose
            as.matrix
dim(t_ASVs_PIT) #13 rows amd 36381 columns
#Need to create a vector when whether or not it's a blank sample
row.names(t_ASVs_PIT)
blanks_pits <- c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE)
contam_pits_vector <- isContaminant(t_ASVs_PIT,neg=blanks_pits, detailed=FALSE) #this returns a logical vector where true is contimate
length(which(contam_pits_vector == TRUE)) #Shows how many are true aka contaminants We got 4
which(contam_pits_vector == TRUE) #Shows the location of which are true contaminants: 23445 25910 27070 32385 
#then add this to the document created by the wc code, and rename the file
#Let's remove the contaminants...From Gerid
#Now I want to remove them from a complete transposed data frame. 
#They should be the same row index as column index, because we didn't change the ASV orderever
#First, we need indexes instead of ASV.ID as rownames
Clean_table<-unrowname(Table_Cor)
#Next, I read in the contaminant vector. This must be read in with the scan fx.
contaminates<-scan("~/Desktop/Chapter1/MESOSCOPE2017/Output/Contaminated_WC_PIT_ASV.txt") #This is the vector that contains the index for every contaminant 
# Filter the table based on the contaminate vector
Filtered_table<-Clean_table[!rownames(Clean_table) %in% contaminates, ]
#Verification steps. the number of rows in the Filtered table should equal nrows OTU_table - length of contaminants.
nrow(Table_Cor) - nrow(Filtered_table) # 10, so it worked!
#To add back the ASV IDs into the header as row names
rownames(Filtered_table) = Filtered_table$ASV.ID #Here I want to give the rownames the ASV IDs
head(Filtered_table) 
dim(Filtered_table) #  36371 rows and  111 columns, 108 columns of samples, first column is ASV names and last column is taxon list
#I'm going to remove the Blank samples because they have served their purpose
names(Filtered_table[c(2,42,94)])
Filtered_table[,c(2, 42, 94)]=NULL
dim(Filtered_table) # 36371 rows and 108 columns 
write.csv(Filtered_table,"~/Desktop/Chapter1/MESOSCOPE2017/Output/MS_Clean_1022.csv")

######## 03/03/23- Calculate ASV stats ###########
clean_counts <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Output/MS_Clean_1022.csv", header=T) 
clean_counts[1]= NULL
rownames(clean_counts) = clean_counts$ASV.ID
head(clean_counts)
# first need to take the average of duplicates
clean_mean <- clean_counts %>%
  select(., -108) %>%  # remove confidence of ASV so as not to confuse
  melt
head(clean_mean)
var<-colsplit(clean_mean$variable, '_',c("Station","Depth","Type","Rep")) #Separate the variables in sample name
clean_mean_combo<-data.frame(clean_mean,var) #combine dataframes
sample_mean <- clean_mean_combo %>%
  dcast(., ASV.ID+Station+Depth+Type~Rep, value.var="value")  %>%  #recast data frame by reps
  mutate(mean = ((A+B)/2))  %>%  #Take the mean of rep A and B
  unite(.,"Sample",c("Station","Depth","Type")) %>%  #regroup the sample names
  dcast(.,ASV.ID~Sample, value.var="mean") #recast the data frame with the mean as the value
head(sample_mean)
dim(sample_mean) 

write.csv(sample_mean, "~/Desktop/Chapter1/MESOSCOPE2017/Output/MS_Clean_mean_0323.csv")
#downloaded it and then combined it with the numbers for the samples without blanks
all_counts<-read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Output/MS_Clean_mean_0323_ALL.csv", header=T)
rownames(all_counts)=all_counts$ASV.ID
head(all_counts)
all_counts.t <- t(all_counts)
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
  labs(title="", x="",y="Total Reads")+theme_bw()+
  theme(axis.text.x = element_blank(),axis.text.y=element_text(size=10,color = "black"), strip.background = element_blank())+ theme(legend.title = element_blank())+
  scale_fill_manual(values=col_mat)
all.plot + scale_y_continuous(trans = 'log10')
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/All_reads_mean_0323.pdf", width=10, height = 10)

# Reports total number of unique ASVs
length(unique(all_counts$ASV.ID)) #36371

# ASV counts and information:
counts_only<-all_counts[2:61]
seq_total<-apply(counts_only,2,sum)
min(seq_total); max(seq_total); mean(seq_total)
ASV_count<-colSums(counts_only>0); ASV_count
min(ASV_count); max(ASV_count); mean(ASV_count)
ASV_single<-colSums(counts_only==1) #by sample
ASV_double<-colSums(counts_only==2) #by sample
ASV_true<-colSums(counts_only>2)
#
#Compile sample information
sample_info<-data.frame(seq_total,ASV_count,ASV_single,ASV_double,ASV_true);sample_info
#Visual representation of ASVs:

sample_info$Sample<-row.names(sample_info)
head(sample_info)
counts.melt<-melt(sample_info[c(6,1:4)])
counts.plot <- join(counts.melt, key_all, by="Sample", type="left", match="first") %>%
  arrange(., Material, Station.sla, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))

head(counts.plot)
means_type <- tapply(counts.plot$value, list(counts.plot$Material, counts.plot$variable), mean)
min_type <- tapply(counts.plot$value, list(counts.plot$Material, counts.plot$variable), min)
max_type <- tapply(counts.plot$value, list(counts.plot$Material, counts.plot$variable), max)
type_stats <- data.frame(means_type, min_type, max_type)
write.csv(means_type, "~/Desktop/Chapter1/MESOSCOPE2017/Output/ASV_means_0323.csv")
write.csv(type_stats, "~/Desktop/Chapter1/MESOSCOPE2017/Output/ASV_stats_0323.csv")

ASVs<-c("ASV_count", "ASV_single", "ASV_double")

All_bar_stats<- ggplot(counts.plot, aes(x=Sample, y=value, fill=variable))+
  geom_bar(stat="identity",position="stack",color="black")+
  labs(title="ASV Counts", x="",y="Total ASVs")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=1,size=12, color = "black"),axis.text.y=element_text(size=12,color = "black"), strip.background = element_blank())+ theme(legend.title = element_blank())
All_bar_stats %+% subset(counts.plot, variable %in% ASVs) + scale_fill_manual("",values=c("#e41a1c","#fee08b","#4393c3"),labels = c("ASVs > 2 seqs", "Singletons", "Doubletons"))
  ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/All_ASV_Stats_0323.pdf", width=10, height = 10)
  
## could break this up if need be, but the main point is that they are mostly more than doubletons
DNA<- All_bar_stats %+% subset(counts.plot, (variable %in% ASVs & Material=="DNA")) + scale_fill_manual("",values=c("#e41a1c","#fee08b","#4393c3"),labels = c("ASVs > 2 seqs", "Singletons", "Doubletons"))
RNA<- All_bar_stats %+% subset(counts.plot, (variable %in% ASVs & Material=="RNA")) + scale_fill_manual("",values=c("#e41a1c","#fee08b","#4393c3"),labels = c("ASVs > 2 seqs", "Singletons", "Doubletons"))
PIT<- All_bar_stats %+% subset(counts.plot, (variable %in% ASVs & Material=="PIT")) + scale_fill_manual("",values=c("#e41a1c","#fee08b","#4393c3"),labels = c("ASVs > 2 seqs", "Singletons", "Doubletons"))

ggarrange(DNA, RNA, PIT, common.legend=TRUE, legend="right")
# may want to make the y axis the same for them all
# may want to change the titles to the material type, and edit the sample labels

############ 4. Remove global singletons ##########
##Filter out ASVs with only 1 sequence in the whole dataset (global singletons)
ASV_sum<-apply(Filtered_table[2:106],1,sum, na.rm =TRUE) #remove global singletons by rows- summing total of reads an ASV is found across all samples
ASV_no1 = Filtered_table[ ASV_sum>1, ]  #count.no1 = ASV table without global singletons
removed = dim(Filtered_table)[1] - dim(ASV_no1)[1] #Outputs the number of ASVs (total) lost in this step, 1,428
names(ASV_no1)[1] = "ASV1.ID"
dim(ASV_no1) #Check dimensions- 34943 rows and 108 columns
write.csv(ASV_no1, "~/Desktop/Chapter1/MESOSCOPE2017/Output/MS_Clean_No1_1022.csv")

############ 5. Take the mean of duplicate samples #########
No1_m <- ASV_no1 %>%
        select(., -108) %>%  # remove confidence of ASV so as not to confuse
        melt
head(No1_m)
var<-colsplit(No1_m$variable, '_',c("Station","Depth","Type","Rep")) #Separate the variables in sample name
No1_m_var<-data.frame(No1_m,var) #combine dataframes
sample_mean <- No1_m_var %>%
  dcast(., ASV1.ID+Station+Depth+Type~Rep) %>%  #recast data frame by reps
  mutate(mean = (No1_mean$A+No1_mean$B)/2) %>%  #Take the mean of rep A and B
  unite(.,"Sample",c("Station","Depth","Type")) %>%  #regroup the sample names
  dcast(.,ASV1.ID~Sample, value.var="mean") #recast the data frame with the mean as the value
head(sample_mean)
dim(sample_mean) # 34943 rows and 62 columns
write.csv(sample_mean, "~/Desktop/Chapter1/MESOSCOPE2017/Output/MS_Clean_No1_Means_1022.csv")

###### 6. Now that we have a clean, mean, no global singletons table, we can rarefy to see sequencing depth ######
#I took the means calculated above, and manually added back the rows that didn't have duplicates
ASV_mean<-read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Output/MS_Clean_No1_Means_Complete_1022.csv",header=T)
rownames(ASV_mean)=ASV_mean$ASV1.ID
dim(ASV_mean) #34943 rows and 63 columns 
t_mean <- ASV_mean[2:61] %>%  # selecting only numeric columns
          as.matrix %>%
          t %>%  # transpose
          floor  # convert to whole numbers
# to get color order to match matrix order        
sample_order <- t_mean %>%
          as.data.frame  %>% # revert back to dataframe
          mutate(Sample = rownames(.)) %>%
          join(., key_all, by="Sample", type ="left", match="first")
          
col_rar <- sample_order$Color.depth
rarefac<-rarecurve(t_mean, step =20, col="blue", lty=1, label=F, 
                   xlab="Number of reads sampled", ylab="Number of ASVs",xlim=c(0,250000))
legend("topright",legend = key_all$Depth ,col=key_all$Color.depth,lty=1,lwd=5)

##### 03/3/23 Trying a new method to visualization the rarefaction curve #######
# from Patt Schloss
options(scipen = 999) # this prevents scientific notation
rarecurve_data <- rarecurve(t_mean, step =20)
col_mat <- unique(key_all$Color.mat)
# c("#7fcdbb","#2c7fb8","#fc8d59")
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
  labs(title="Rarefaction Curve of ASVs", x="Number of sequences",y="ASVs")+
  theme_bw()
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Rarefaction_0323.pdf", width=6, height = 6)

############# 7. Normalize ##################
#From Gerid
ListDGE = DGEList(ASV_mean[,2:61]) 
ListDGE
ListDGE = calcNormFactors(ListDGE, method = "TMM")
ListDGE
TMMNorm_ASV = cpm(ListDGE)
head(TMMNorm_ASV)
TMMNorm_ASV = as.data.frame(TMMNorm_ASV)
TMMNorm_ASV$ASV1.ID = row.names(TMMNorm_ASV)
Joined<-join(TMMNorm_ASV, ASV_mean[c(1,62:63)], by="ASV1.ID", type="left", match="first") 
dim(Joined) # 34943 rows and 61 columns
write.csv(Joined, "~/Desktop/Chapter1/MESOSCOPE2017/Output/MS_Clean_No1_Means_Complete_Norm_1022.csv")
############# 8. Rename the taxa from PR2, into nicer names ##########
##### Updated 07/02/20 to reflect the PR2, V12
##### Updated 11/17/22 to reflect the PR2, V14
table <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Output/MS_Clean_No1_Means_Complete_Norm_1022.csv", header=T)
new_tax <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/taxonomy_Nov22.csv", header=T)
table[1] = NULL
table[62:63] = NULL
names(new_tax)[1] = "ASV1.ID"

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
  # split$Taxa[split$Level3=="Foraminifera"]="Rhizaria-Foraminifera"
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

NewNameTax = pr2_rename_taxa_w2(Joined)
combined_newtax = data.frame(Joined, NewNameTax)
combined_newtax[1:5,]
write.csv(combined_newtax, "~/Desktop/Chapter1/MESOSCOPE2017/Output/MS_Clean_No1_Means_Complete_Norm_Renamed_1022.csv")
asv_df <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Output/MS_Clean_No1_Means_Complete_Norm_Renamed_1022.csv", header=T)
asv_df[1]=NULL
row.names(asv_df) = asv_df$ASV1.ID

###### START HERE FOR ASV TABLE WITH UPDATED TAXONOMY ##########
combined_newnametax = join(table, new_tax, by="ASV1.ID", type="left")
NewNameNewTax = pr2_rename_taxa_w2(combined_newnametax)
combined_newnamenewtax = data.frame(combined_newnametax,NewNameNewTax)
write.csv(combined_newnamenewtax, "~/Desktop/Chapter1/MESOSCOPE2017/Output/MS_Clean_No1_Means_Complete_Norm_NewTax_1122.csv")
asv_df <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Output/MS_Clean_No1_Means_Complete_Norm_NewTax_1122.csv", header=T)
asv_df[1]=NULL
row.names(asv_df) = asv_df$ASV1.ID
no_meta <- asv_df[asv_df$Level3 != "Metazoa",]

############ 9. Plot the abundance #########
### with new names and new colors 12/1/22
data.m<-melt(asv_df) #melt
head(data.m)
data.agg<-aggregate(data.m$value, by=list(Taxa=data.m$TaxaPlot,Sample=data.m$variable),sum) #sum sequences by taxonomic group
data.agg$Sample<-factor(data.agg$Sample,levels=names(asv_df[2:61]))
tax_order<-c("Alveolates-Ciliates","Alveolates-Dinophyceae","Alveolates-Syndiniales", "Archaeplastids-Chlorophytes",
             "Excavates-Discoba", "Hacrobia-Cryptophytes","Hacrobia-Haptophytes",
             "Opisthokont-Choanoflagellida","Opisthokonts-Metazoa","Other/unknown",
             "Rhizaria-Acantharia","Rhizaria-Cercozoa","Rhizaria-Polycystines", "Rhizaria-Other",
             "Stramenopiles-Chrysophytes","Stramenopiles-Diatoms","Stramenopiles-MAST","Stramenopiles-Ochrophyta","Stramenopiles-Pelagophytes","Stramenopiles-Other")
tax_color_phyla2<-c('firebrick4','indianred1','tomato3','forestgreen','yellowgreen','darkblue',
                    'lightblue','moccasin','gold1',"light grey","#C291A4",'mediumvioletred','magenta','lightpink',"#DECDBE",'#DDAD4B','tan2','tan3','tan4',"#5C4033")
names(tax_color_phyla2)<-tax_order
bar_plot_df <- join(data.agg, key_all, by="Sample", type="left", match="first") # %>%
arrange(., Material, Station, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))

#Bar plot of community composition ordered by Material, Station and Depth
bars<-ggplot(bar_plot_df, aes(y=x, fill=tax, x=Station.sla))+
  geom_bar(position = "fill",stat="identity",color="black",aes(fill=Taxa))+
  scale_fill_manual(values=tax_color_phyla2)+
  theme_bw()+
  labs(title="Relative abundance of reads", x="Sample",y="")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5, vjust = .5))

bars
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/Taxa_All_Meta_1122Taxa.pdf", width=15, height = 7)
bars %+% facet_wrap(~Material, ncol = 2, scales = "free_x")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/Taxa_All_Meta_1122Taxa_Material.pdf", width=10, height = 10)
bars %+% subset(bar_plot_df, Material == "DNA") + facet_wrap(~Depth , scales = "free_x", ncol = 1)
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/Taxa_All_Meta_1122Taxa_DNA.pdf", width=8, height = 20)
# barplot %+% subset(bar_plot_df, Material == "RNA") + facet_wrap(~Station , scales = "free_x", ncol = 1)
# ggplot2::ggsave("Taxa_ASV_RNA_Station.pdf", width=8, height = 20)

WC <- c("DNA","RNA")
bars %+% subset(bar_plot_df, Material %in% WC) + coord_flip() +facet_wrap(.~ Material + Depth, scales = "free_y", ncol=2, dir="v")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/Taxa_All_Meta_1122Taxa_WC_Depth.pdf", width=15, height = 10)

###### NEW TAXA BARPLOT, first round with outdated names#####
data.m<-melt(asv_df) #melt
head(data.m)
data.agg<-aggregate(data.m$value, by=list(Taxa=data.m$Taxa,Sample=data.m$variable),sum) #sum sequences by taxonomic group
data.agg$Sample<-factor(data.agg$Sample,levels=names(asv_df[2:61]))
# tax_order<-c("Alveolates-Dinophyceae","Alveolates-Ciliates","Alveolates-Syndiniales","Alveolates-Other",
  #           "Stramenopiles-Diatoms","Stramenopiles-Pelagophytes","Stramenopiles-MAST","Stramenopiles-Chrysophytes","Stramenopiles-Ochrophyta","Stramenopiles-Other",
   #          "Archaeplastids-Chlorophytes","Cryptophytes","Haptophytes",
    #         "Rhizaria-Acantharia","Rhizaria-Cercozoa","Rhizaria-Foraminifera","Rhizaria-Polycystines","Rhizaria-RAD (A,B,C)","Rhizaria-Other","Rhizaria-Radiolaria-Other",
     #        "Opisthokont-Choanoflagellida","Opisthokonts-Other","Other/unknown")
tax_order<-c("Alveolates-Ciliates","Alveolates-Dinophyceae","Alveolates-Syndiniales", "Archaeplastids-Chlorophytes",
             "Excavates-Discoba", "Hacrobia-Cryptophytes","Hacrobia-Haptophytes",
             "Opisthokont-Choanoflagellida","Opisthokonts-Metazoa","Other/unknown",
             "Rhizaria-Acantharia","Rhizaria-Cercozoa","Rhizaria-Polycystines", "Rhizaria-Other",
             "Stramenopiles-Chrysophytes","Stramenopiles-Diatoms","Stramenopiles-MAST","Stramenopiles-Ochrophyta","Stramenopiles-Pelagophytes","Stramenopiles-Other")

# tax_color=c("#A6CEE3", "#1F78B4","#daa520", "#33A02C", "#FFFFB3", "#BEBADA","#FB9A99", "#E31A1C", "#01665E", "#FF7F00", "#DE77AE", "#6A3D9A", "chocolate2", "#5E4FA2",
 #           "indianred1",  "mediumpurple2","#E5C494","#ea7e5d", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#A65628", "#CCEBC5", "plum",
  #          "#66C2A5", "#8DA0CB", "#e52b50","#E78AC3", "#A6D854","#FFD92F")

# color_tax<-c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99",
     #        "#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5")
# cols<-c("#6a3d9a","#9970ab","#cab2d6","#bebada","#dfc27d","#fa9fb5","#dd3497","#a6cee3","#1f78b4","#d9d9d9","#FFD92F","#cc4c02",
     #   "#fee391","#ff7f00","#fdb462","#ffbba8","#a1d99b","#b3de69","#33a02c","#8dd3c7","#b2df8a","#d9ef8b")
tax_color_phyla2<-c('firebrick4','indianred1','tomato3','forestgreen','yellowgreen','darkblue',
                    'lightblue','moccasin','gold1','grey',"#C291A4",'mediumvioletred','magenta','lightpink',"#DECDBE",'#DDAD4B','tan2','tan3','tan4',"#5C4033")

names(tax_color)<-tax_order
data.agg$tax<-factor(data.agg$Taxa, levels=rev(tax_order)) #factoring
bar_plot_df <- join(data.agg, key_all, by="Sample", type="left", match="first") %>%
  arrange(., Material, Station, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))

#Bar plot of community composition ordered by Material, Station and Depth
ggplot(bar_plot_df, aes(y=x, fill=tax, x=Sample))+
  geom_bar(position = "fill",stat="identity",color="black",aes(fill=Taxa))+
  scale_fill_manual(values=tax_color_phyla2)+
  theme_bw()+
  labs(title="Relative abundance of reads", x="Sample",y="")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5, vjust = .5))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/New_Taxa_All_Meta.pdf", width=15, height = 7)

# Now with only metazoans
meta_only <- combined_newnamenewtax[combined_newnamenewtax$Level3 == "Metazoa",]
meta_only[63]= NULL
data.meta.m<-melt(meta_only) #melt
head(data.meta.m)
data.meta.agg<-aggregate(data.meta.m$value, by=list(Taxa=data.meta.m$Taxa2,Sample=data.meta.m$variable),sum) #sum sequences by taxonomic group
data.meta.agg$Sample<-factor(data.meta.agg$Sample,levels=names(meta_only[2:61]))
meta_bar_plot_df <- join(data.meta.agg, key_all, by="Sample", type="left", match="first") %>%
  arrange(., Material, Station, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))

color_meta<-c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a",
              "#ffda03", "#ffff99", "#8dd3c7","#fb8072")
ggplot(meta_bar_plot_df, aes(y=x, fill=Taxa, x=Sample))+
  geom_bar(position = "fill",stat="identity",color="black",aes(fill=Taxa))+
  scale_fill_manual(values=color_meta)+
  theme_bw()+
  labs(title="Relative abundance of reads", x="Sample",y="")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5, vjust = .5))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/New_Taxa_Only_Meta.pdf", width=15, height = 7)

## now without metazoans
no_meta <- asv_df[asv_df$Level3 != "Metazoa",]
data.nometa.m<-melt(no_meta) #melt
head(data.nometa.m)
data.nometa.agg<-aggregate(data.nometa.m$value, by=list(Taxa=data.nometa.m$Taxa2,Sample=data.nometa.m$variable),sum) #sum sequences by taxonomic group
data.nometa.agg$Sample<-factor(data.nometa.agg$Sample,levels=names(no_meta[2:61]))
nometa_bar_plot_df <- join(data.nometa.agg, key_all, by="Sample", type="left", match="first") %>%
  arrange(., Material, Station, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))


ggplot(nometa_bar_plot_df, aes(y=x, fill=tax, x=Sample))+
  geom_bar(position = "fill",stat="identity",color="black",aes(fill=Taxa))+
  scale_fill_manual(values=tax_color)+
  theme_bw()+
  labs(title="Relative abundance of reads", x="Sample",y="")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5, vjust = .5))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/New_Taxa_All_NoMeta.pdf", width=15, height = 7)


######### Now with the samples ordered by SLA level not numeric station #######
data.m<-melt(asv_df) #melt
head(data.m)
data.agg<-aggregate(data.m$value, by=list(Taxa=data.m$TaxaPlot,Sample=data.m$variable),sum) #sum sequences by taxonomic group
data.agg$Sample<-factor(data.agg$Sample,levels=names(asv_df[2:61]))
tax_order<-c("Alveolates-Ciliates","Alveolates-Dinophyceae","Alveolates-Syndiniales", "Archaeplastids-Chlorophytes",
             "Excavates-Discoba", "Hacrobia-Cryptophytes","Hacrobia-Haptophytes",
             "Opisthokont-Choanoflagellida","Opisthokonts-Metazoa","Other/unknown",
             "Rhizaria-Acantharia","Rhizaria-Cercozoa","Rhizaria-Polycystines", "Rhizaria-Other",
             "Stramenopiles-Chrysophytes","Stramenopiles-Diatoms","Stramenopiles-MAST","Stramenopiles-Ochrophyta","Stramenopiles-Pelagophytes","Stramenopiles-Other")
tax_color_phyla2<-c('firebrick4','indianred1','tomato3','forestgreen','yellowgreen','darkblue',
                    'lightblue','moccasin','gold1','grey',"#C291A4",'mediumvioletred','magenta','lightpink',"#DECDBE",'#DDAD4B','tan2','tan3','tan4',"#5C4033")
names(tax_color_phyla2)<-tax_order
bar_plot_df <- join(data.agg, key_all, by="Sample", type="left", match="first") %>%
  arrange(., Material, Station.sla, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))

#Bar plot of community composition ordered by Material, Station and Depth
bars<-ggplot(bar_plot_df, aes(y=x, fill=tax, x=Station.sla))+
  geom_bar(position = "fill",stat="identity",color="black",aes(fill=Taxa))+
  scale_fill_manual(values=tax_color_phyla2)+
  theme_bw()+
  labs(title="Relative abundance of reads", x="Sample",y="")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5, vjust = .5))
WC <- c("DNA","RNA")
bars %+% subset(bar_plot_df, Material %in% WC) + coord_flip() +facet_wrap(.~ Material + Depth, scales = "free_y", ncol=2, dir="v")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/Taxa_All_Meta_1222Taxa_WC_DepthSLA.pdf", width=17, height = 10)

#now with no meta
nometa <-filter(asv_df, Taxa!="Opisthokont-Metazoa")
data_nometa<-melt(nometa) #melt
head(data_nometa)
data.agg_nometa<-aggregate(data_nometa$value, by=list(Taxa=data_nometa$TaxaPlot,Sample=data_nometa$variable),sum) #sum sequences by taxonomic group
data.agg_nometa$Sample<-factor(data.agg_nometa$Sample,levels=names(nometa[2:61]))
tax_color_phyla2_nometa<-c('firebrick4','indianred1','tomato3','forestgreen','yellowgreen','darkblue',
                    'lightblue','moccasin','grey',"#C291A4",'mediumvioletred','magenta','lightpink',"#DECDBE",'#DDAD4B','tan2','tan3','tan4',"#5C4033")

bar_plot_df_nometa <- join(data.agg_nometa, key_all, by="Sample", type="left", match="first") %>%
  arrange(., Material, Station.sla, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))

#Bar plot of community composition ordered by Material, Station and Depth
bars_nometa<-ggplot(bar_plot_df, aes(y=x, fill=tax, x=Station.sla))+
  geom_bar(position = "fill",stat="identity",color="black",aes(fill=Taxa))+
  scale_fill_manual(values=tax_color_phyla2_nometa)+
  theme_bw()+
  labs(title="Relative abundance of reads", x="Sample",y="")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5, vjust = .5))

bars_nometa %+% subset(bar_plot_df_nometa, Material %in% WC) + coord_flip() +facet_wrap(.~ Material + Depth, scales = "free_y", ncol=2, dir="v")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/Taxa_All_NOMeta_1222Taxa_WC_DepthSLA.pdf", width=17, height = 10)

########## 02/13/23 Now ordered by density #######
data.m<-melt(asv_df) #melt
head(data.m)
data.agg<-aggregate(data.m$value, by=list(Taxa=data.m$TaxaPlot,Sample=data.m$variable),sum) #sum sequences by taxonomic group
data.agg$Sample<-factor(data.agg$Sample,levels=names(asv_df[2:61]))
tax_order<-c("Alveolates-Ciliates","Alveolates-Dinophyceae","Alveolates-Syndiniales", "Archaeplastids-Chlorophytes",
             "Excavates-Discoba", "Hacrobia-Cryptophytes","Hacrobia-Haptophytes",
             "Opisthokont-Choanoflagellida","Opisthokonts-Metazoa","Other/unknown",
             "Rhizaria-Acantharia","Rhizaria-Cercozoa","Rhizaria-Polycystines", "Rhizaria-Other",
             "Stramenopiles-Chrysophytes","Stramenopiles-Diatoms","Stramenopiles-MAST","Stramenopiles-Ochrophyta","Stramenopiles-Pelagophytes","Stramenopiles-Other")
tax_color_phyla2<-c('firebrick4','indianred1','tomato3','forestgreen','yellowgreen','darkblue',
                    'lightblue','moccasin','gold1','grey',"#C291A4",'mediumvioletred','magenta','lightpink',"#DECDBE",'#DDAD4B','tan2','tan3','tan4',"#5C4033")
names(tax_color_phyla2)<-tax_order
bar_plot_df <- join(data.agg, key_all, by="Sample", type="left", match="first") %>%
  arrange(., Material, Potential_Density, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)), Sample.order = factor(Stat.depth.density, levels=rev(unique(Stat.depth.density))))


bars<-ggplot(bar_plot_df, aes(y=x, fill=tax, x=Sample.order))+
  geom_bar(position = "fill",stat="identity",color="black",aes(fill=Taxa))+
  scale_fill_manual(values=tax_color_phyla2)+
  theme_bw()+
  labs(title="Relative abundance of reads", x="Sample",y="")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5, vjust = .5))
  
WC <- c("DNA","RNA")
bars %+% subset(bar_plot_df, Material %in% WC) + coord_flip() +facet_wrap(.~ Material, scales = "free_y", ncol=2, dir="v")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/Taxa_All_Meta_0223Taxa_WC_Density.pdf", width=20, height = 10)

#nometa
nometa <-filter(asv_df, Taxa!="Opisthokont-Metazoa")
data_nometa<-melt(nometa) #melt
head(data_nometa)
data.agg_nometa<-aggregate(data_nometa$value, by=list(Taxa=data_nometa$TaxaPlot,Sample=data_nometa$variable),sum) #sum sequences by taxonomic group
data.agg_nometa$Sample<-factor(data.agg_nometa$Sample,levels=names(nometa[2:61]))
tax_color_phyla2_nometa<-c('firebrick4','indianred1','tomato3','forestgreen','yellowgreen','darkblue',
                           'lightblue','moccasin','grey',"#C291A4",'mediumvioletred','magenta','lightpink',"#DECDBE",'#DDAD4B','tan2','tan3','tan4',"#5C4033")

bar_plot_df_nometa <- join(data.agg_nometa, key_all, by="Sample", type="left", match="first") %>%
  arrange(., Material, Potential_Density, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)), Sample.order = factor(Stat.depth.density, levels=rev(unique(Stat.depth.density))))

#Bar plot of community composition ordered by Material, Station and Depth
bars_nometa<-ggplot(bar_plot_df_nometa, aes(y=x, fill=tax, x=Sample.order))+
  geom_bar(position = "fill",stat="identity",color="black",aes(fill=Taxa))+
  scale_fill_manual(values=tax_color_phyla2_nometa)+
  theme_bw()+
  labs(title="Relative abundance of reads", x="Sample",y="")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5, vjust = .5))

bars_nometa %+% subset(bar_plot_df_nometa, Material %in% WC) + coord_flip() +facet_wrap(.~ Material, scales = "free_y", ncol=2, dir="v")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/Taxa_All_NOMeta_0223Taxa_WC_density.pdf", width=20, height = 10)

# want it arranged by station across the transect normally

bars_nometa_station<-ggplot(bar_plot_df_nometa, aes(y=x, fill=tax, x=Station))+
  geom_bar(position = "fill",stat="identity",color="black",aes(fill=Taxa))+
  scale_fill_manual(values=tax_color_phyla2_nometa)+
  theme_bw()+
  labs(title="Relative abundance of reads", x="Sample",y="")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5, vjust = .5))

bars_nometa_station %+% subset(bar_plot_df_nometa, Material %in% WC) + coord_flip() +facet_wrap(.~ Material+Depth, scales = "free_y", ncol=2, dir="v")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/Taxa_All_NOMeta_0223Taxa_WC_station.pdf", width=17, height = 10)

### 10. NMDS of community composition color by depth #######
#Hybdrid of code from Sam and Sarah
# Blanks removed, Samples averaged, normalized
# first with meta included 
#First need to generate relative abundance
head(asv_df)
Norm.rel <-decostand(asv_df[,2:61], MARGIN=2, method = "total") #this does it by rows
colSums(Norm.rel) # check! should all equal 1.
Norm.rel$ASV.ID<-rownames(Norm.rel)
melt_norm<- melt(Norm.rel)
vars1<-colsplit(melt_norm$variable, "_", c("Station","Depth","Material"))
new_df<-data.frame(melt_norm, vars1)

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
###### 12/22 Updated Tax NMDS with Metazoan ########
all<- nmdsPoints(new_df, key_all)
nmds_all<-ggplot(data=all, aes(x = MDS1, y = MDS2, fill=Loc,shape=Material)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = c(21,25,22))+ scale_fill_manual(values=(levels(all$Color.sla)))+
  labs(title="NMDS-Bray-Curtis -No Metazoa", x="NMDS1",y="NMDS2",size=12)+
  guides(fill = guide_legend(override.aes=list(color=levels(all$Color.sla))))+
  theme(legend.title = element_blank())
nmds_all_nometa +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Material)) #this adds in 95% confidence interval ellipses from ggplot
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_all.pdf", width=6, height = 6)

pits <- filter(new_df, Material=="PIT")
pit_df<-nmdsPoints(pits, key_all)
color_pits <- unique(factor(as.character(pit_df_nometa$Color.sla), levels=c("#de2d26","#cfcdcd","#2c7fb8")))
nmds_pit<-ggplot(data=pit_df, aes(x = MDS1, y = MDS2, fill=Loc,shape=Material)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = c(25))+ scale_fill_manual(values=levels(color_pits))+
  labs(title="NMDS-Bray-Curtis: PITs-No Metazoa", x="NMDS1",y="NMDS2",size=12)+
  guides(fill = guide_legend(override.aes=list(color=levels(color_pits))))+
  theme(legend.title = element_blank())
nmds_pit
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_pit.pdf", width=6, height = 6)

rna <- filter(new_df, Material=="cDNA")
rna_df<-nmdsPoints(rna, key_all)
nmds_rna<-ggplot(data=rna_df, aes(x = MDS1, y = MDS2, fill=Loc,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_df$Shape.depth)+ scale_fill_manual(values=levels(rna_df$Color.sla))+
  labs(title="NMDS-Bray-Curtis: RNA-No Metazoa", x="NMDS1",y="NMDS2",size=12)+
  guides(fill = guide_legend(override.aes=list(color=levels(rna_df$Color.sla))))+
  theme(legend.title = element_blank())
nmds_rna + stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna.pdf", width=8, height = 8)

dna <- filter(new_df, Material=="DNA")
dna_df<-nmdsPoints(dna, key_all)
nmds_dna<-ggplot(data=dna_df, aes(x = MDS1, y = MDS2, fill=Loc,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = dna_df$Shape.depth)+ scale_fill_manual(values=levels(dna_df$Color.sla))+
  labs(title="NMDS-Bray-Curtis: DNA-No Metazoa", x="NMDS1",y="NMDS2",size=12)+
  guides(fill = guide_legend(override.aes=list(color=levels(dna_df$Color.sla))))+
  theme(legend.title = element_blank())
nmds_dna + stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_dna.pdf", width=8, height = 8)


dna_15 <-filter(dna, Depth=="15m")
dna_15_df<-nmdsPoints(dna_15, key_all)
nmds_dna_15<-ggplot(data=dna_15_df, aes(x = MDS1, y = MDS2, fill=Loc, shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  scale_y_continuous(limits = c(-.2, .2))+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = dna_15_df$Shape.depth)+ scale_fill_manual(values=levels(dna_15_df$Color.sla))+
  labs(title="NMDS-Bray-Curtis: DNA-15m", x="NMDS1",y="NMDS2",size=12)+
  guides(fill = guide_legend(override.aes=list(color=dna_15_df$Color.sla)))+
  theme(legend.title = element_blank())
nmds_dna_15 
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_dna_15m.pdf", width=8, height = 8)

DCM_dna<-filter(dna, Depth=="DCM")
dcm_dna_df<-nmdsPoints(DCM_dna, key_all)
nmds_dcm_dna<-ggplot(data=dcm_dna_df, aes(x = MDS1, y = MDS2, fill=Loc, shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  scale_y_continuous(limits = c(-.4, .4))+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = dcm_dna_df$Shape.depth)+ scale_fill_manual(values=levels(dcm_dna_df$Color.sla))+
  labs(title="NMDS-Bray-Curtis: DNA-DCM", x="NMDS1",y="NMDS2",size=12)+
  guides(fill = guide_legend(override.aes=list(color=levels(dcm_dna_df$Color.sla))))+
  theme(legend.title = element_blank())
nmds_dcm_dna 
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_dcm_dna.pdf", width=8, height = 8)

dna_175<-filter(dna, Depth=="175m")
dna_175_df<-nmdsPoints(dna_175, key_all)
nmds_dna_175<-ggplot(data=dna_175_df, aes(x = MDS1, y = MDS2, fill=Loc, shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  scale_y_continuous(limits = c(-.6, .6))+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = dna_175_df$Shape.depth)+ scale_fill_manual(values=levels(dna_175_df$Color.sla))+
  labs(title="NMDS-Bray-Curtis: DNA-175m", x="NMDS1",y="NMDS2",size=12)+
  guides(fill = guide_legend(override.aes=list(color=levels(dna_175_df$Color.sla))))+
  theme(legend.title = element_blank())
nmds_dna_175
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_dna_175.pdf", width=8, height = 8)

dna_500<-filter(dna, Depth=="500m")
dna_500_df<-nmdsPoints(dna_500, key_all)
nmds_dna_500<-ggplot(data=dna_500_df, aes(x = MDS1, y = MDS2, fill=Loc, shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  scale_y_continuous(limits = c(-.4, .4))+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = dna_500_df$Shape.depth)+ scale_fill_manual(values=levels(dna_500_df$Color.sla))+
  labs(title="NMDS-Bray-Curtis: DNA-500m", x="NMDS1",y="NMDS2",size=12)+
  guides(fill = guide_legend(override.aes=list(color=levels(dna_500_df$Color.sla))))+
  theme(legend.title = element_blank())
nmds_dna_500
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_dna_500m.pdf", width=8, height = 8)

ggarrange(nmds_dna_15,nmds_dcm_dna,nmds_dna_175,nmds_dna_500)
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_dna_depths.pdf", width=10, height = 10)


rna_dcm <- filter(rna, Depth=="DCM")
rna_dcm_df<-nmdsPoints(rna_dcm, key_all)
nmds_rna_dcm <-ggplot(data=rna_dcm_df, aes(x = MDS1, y = MDS2, fill=Loc,shape=Depth)) + 
  geom_point(color="black",size=6, position="jitter") + 
  theme_bw() +
  coord_equal() +
  scale_y_continuous(limits = c(-.4, .4))+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_dcm_df$Shape.depth)+ scale_fill_manual(values=levels(rna_dcm_df$Color.sla))+
  labs(title="NMDS-Bray-Curtis: RNA-DCM", x="NMDS1",y="NMDS2",size=12)+
  guides(fill = guide_legend(override.aes=list(color=levels(rna_dcm_df$Color.sla))))+
  theme(legend.title = element_blank())
nmds_rna_dcm
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_dcm.pdf", width=8, height = 8)

rna_15 <- filter(rna, Depth=="15m")
rna_15_df<-nmdsPoints(rna_15, key_all)
nmds_rna_15 <-ggplot(data=rna_15_df, aes(x = MDS1, y = MDS2, fill=Loc,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  scale_y_continuous(limits = c(-.5, .4))+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_15_df$Shape.depth)+ scale_fill_manual(values=levels(rna_15_df$Color.sla))+
  labs(title="NMDS-Bray-Curtis: RNA-15m", x="NMDS1",y="NMDS2",size=12)+
  guides(fill = guide_legend(override.aes=list(color=levels(rna_15_df$Color.sla))))+
  theme(legend.title = element_blank())
nmds_rna_15
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_15m.pdf", width=8, height = 8)

rna_175 <- filter(rna, Depth=="175m")
rna_175_df<-nmdsPoints(rna_175, key_all)
nmds_rna_175 <-ggplot(data=rna_175_df, aes(x = MDS1, y = MDS2, fill=Loc,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  scale_y_continuous(limits = c(-.4, .4))+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_175_df$Shape.depth)+ scale_fill_manual(values=levels(rna_175_df$Color.sla))+
  labs(title="NMDS-Bray-Curtis: RNA-175m", x="NMDS1",y="NMDS2",size=12)+
  guides(fill = guide_legend(override.aes=list(color=levels(rna_175_df$Color.sla))))+
  theme(legend.title = element_blank())
nmds_rna_175
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_175.pdf", width=8, height = 8)

rna_500 <- filter(rna, Depth=="500m")
rna_500_df<-nmdsPoints(rna_500, key_all)
nmds_rna_500 <-ggplot(data=rna_500_df, aes(x = MDS1, y = MDS2, fill=Loc,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  scale_y_continuous(limits = c(-.4, .4))+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_500_df$Shape.depth)+ scale_fill_manual(values=levels(rna_500_df$Color.sla))+
  labs(title="NMDS-Bray-Curtis: RNA-500m", x="NMDS1",y="NMDS2",size=12)+
  guides(fill = guide_legend(override.aes=list(color=levels(rna_500_df$Color.sla))))+
  theme(legend.title = element_blank())
nmds_rna_500
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_500m.pdf", width=8, height = 8)

ggarrange(nmds_rna_15,nmds_rna_dcm,nmds_rna_175,nmds_rna_500)
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_depths.pdf", width=10, height = 10)

####### 12/22 Updated Tax NMDS without Metazoan #########
head(asv_df)
nometa <-filter(asv_df, Taxa!="Opisthokont-Metazoa")
Norm.rel_nometa <-decostand(nometa[,2:61], MARGIN=2, method = "total") #this does it by rows
colSums(Norm.rel_nometa) # check! should all equal 1.
Norm.rel_nometa$ASV.ID<-rownames(Norm.rel_nometa)
melt_norm_nometa<- melt(Norm.rel_nometa)
vars_nometa<-colsplit(melt_norm_nometa$variable, "_", c("Station","Depth","Material"))
new_df_nometa<-data.frame(melt_norm_nometa, vars_nometa)

all_nometa<- nmdsPoints(new_df_nometa, key_all)
nmds_all_nometa<-ggplot(data=all_nometa, aes(x = MDS1, y = MDS2, fill=Loc,shape=Material)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = c(21,25,22))+ scale_fill_manual(values=(levels(all_nometa$Color.sla)))+
  labs(title="NMDS-Bray-Curtis -No Metazoa", x="NMDS1",y="NMDS2",size=12)+
  guides(fill = guide_legend(override.aes=list(color=levels(all_nometa$Color.sla))))+
  theme(legend.title = element_blank())
nmds_all_nometa +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Material)) #this adds in 95% confidence interval ellipses from ggplot
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_all_nometa.pdf", width=6, height = 6)

pits_nometa <- filter(new_df_nometa, Material=="PIT")
pit_df_nometa<-nmdsPoints(pits_nometa, key_all)
color_pits_nometa <- unique(factor(as.character(pit_df_nometa$Color.sla), levels=c("#de2d26","#cfcdcd","#2c7fb8")))
nmds_pit_nometa<-ggplot(data=pit_df_nometa, aes(x = MDS1, y = MDS2, fill=Loc,shape=Material)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = c(25))+ scale_fill_manual(values=levels(color_pits_nometa))+
  labs(title="NMDS-Bray-Curtis: PITs-No Metazoa", x="NMDS1",y="NMDS2",size=12)+
  guides(fill = guide_legend(override.aes=list(color=levels(color_pits_nometa))))+
  theme(legend.title = element_blank())
nmds_pit_nometa
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_pit_nometa.pdf", width=6, height = 6)

rna_nometa <- filter(new_df_nometa, Material=="cDNA")
rna_df_nometa<-nmdsPoints(rna_nometa, key_all)
nmds_rna_nometa<-ggplot(data=rna_df_nometa, aes(x = MDS1, y = MDS2, fill=Loc,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_df_nometa$Shape.depth)+ scale_fill_manual(values=levels(rna_df_nometa$Color.sla))+
  labs(title="NMDS-Bray-Curtis: RNA-No Metazoa", x="NMDS1",y="NMDS2",size=12)+
  guides(fill = guide_legend(override.aes=list(color=levels(rna_df_nometa$Color.sla))))+
  theme(legend.title = element_blank())
nmds_rna_nometa + stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_nometa.pdf", width=8, height = 8)


nmds_rna_nometa_sla<-ggplot(data=rna_df_nometa, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
  geom_point(color="black",size=5,position="jitter") + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_df_nometa$Shape.depth)+ 
  # scale_fill_manual(values=all_nometa$Col.sla)+
  labs(title="NMDS-Bray-Curtis: RNA- No Metazoa", x="NMDS1",y="NMDS2",size=12)+
  scale_y_continuous(breaks = seq(-3, 3, by = 1))+
  scale_x_continuous(breaks = seq(-4, 2, by = 1))+
  scale_fill_gradient2(low = "#0818A8",
                       mid = "white",
                       high = "#880808",
                       midpoint = 0,
                       limits=c(-28,28))
# guides(fill = guide_legend(override.aes=list(values=all_nometa$Col.sla)))
nmds_rna_nometa_sla +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth),linewidth=1,linetype="dashed") + scale_color_manual(values=unique(rna_df_nometa$Color.depth)) #this adds in 95% confidence interval ellipses from ggplot
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_nometa_sla.pdf", width=6, height = 6)

dna_nometa <- filter(new_df_nometa, Material=="DNA")
dna_df_nometa<-nmdsPoints(dna_nometa, key_all)
nmds_dna_nometa<-ggplot(data=dna_df_nometa, aes(x = MDS1, y = MDS2, fill=Loc,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = dna_df_nometa$Shape.depth)+ scale_fill_manual(values=levels(dna_df_nometa$Color.sla))+
  labs(title="NMDS-Bray-Curtis: DNA-No Metazoa", x="NMDS1",y="NMDS2",size=12)+
  guides(fill = guide_legend(override.aes=list(color=levels(dna_df_nometa$Color.sla))))+
  theme(legend.title = element_blank())
nmds_dna_nometa + stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_dna_nometa.pdf", width=8, height = 8)

nmds_dna_nometa_sla<-ggplot(data=dna_df_nometa, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
  geom_point(color="black",size=5,position="jitter") + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = dna_df_nometa$Shape.depth)+ 
  # scale_fill_manual(values=all_nometa$Col.sla)+
  labs(title="NMDS-Bray-Curtis: DNA- No Metazoa", x="NMDS1",y="NMDS2",size=12)+
  scale_y_continuous(breaks = seq(-3, 3, by = 1))+
  scale_x_continuous(breaks = seq(-4, 2, by = 1))+
  scale_fill_gradient2(low = "#0818A8",
                       mid = "white",
                       high = "#880808",
                       midpoint = 0,
                       limits=c(-28,28))
# guides(fill = guide_legend(override.aes=list(values=all_nometa$Col.sla)))
nmds_dna_nometa_sla +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth),linewidth=1,linetype="dashed") + scale_color_manual(values=unique(dna_df_nometa$Color.depth)) #this adds in 95% confidence interval ellipses from ggplot
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_dna_nometa_sla.pdf", width=6, height = 6)





####### 12/16/22- Adding environmental arrows #####
# first plot way from: https://ourcodingclub.github.io/tutorials/ordination/ 
library(caret)
library(grid)
metadata <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/metadata_edited.csv", header=T)
metadata$identifier <- paste(metadata$Station, metadata$Depth_var, sep="_")
row.names(metadata) = metadata$identifier
meta <- select(metadata, -c("Station","Depth_var","identifier")) # selecting only the numeric variables
dna_df_nometa$identifer = paste(dna_df_nometa$Station, dna_df_nometa$Depth.num, sep="_")
row.names(dna_df_nometa) = dna_df_nometa$identifer
ef <- envfit(select(dna_df_nometa, c("MDS1","MDS2")), meta, permu = 999, na.rm=TRUE)
plot(select(dna_df_nometa, c("MDS1","MDS2")))
plot(ef, p.max = 0.05)

# Normalizing the environmental variables to see which might be impacting it
# from: https://www.digitalocean.com/community/tutorials/normalize-data-in-r 
process <- preProcess(as.data.frame(meta), method=c("range"))
norm_scale <- predict(process, as.data.frame(meta))
ef_2 <- envfit(select(dna_df_nometa, c("MDS1","MDS2")), norm_scale, permu = 999, na.rm=TRUE)
plot(select(dna_df_nometa, c("MDS1","MDS2")))
plot(ef_2, p.max = 0.05)

#trying another way to plot them
# from : https://stackoverflow.com/questions/14711470/plotting-envfit-vectors-vegan-package-in-ggplot2
spp.scrs <- as.data.frame(scores(ef_2, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Samples = rownames(spp.scrs))

p <- ggplot(dna_df_nometa) +
  geom_point(mapping = aes(x = MDS1, y = MDS2, colour = Depth)) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data = spp.scrs, aes(x = MDS1, y = MDS2, label = Samples),
            size = 3)

# to select only the the things are significant
meta_exclusive <- select(norm_scale, -c("Salinity","Oxygen","NO3_NO2","TDN","sd","TDP","PC","Heterotrophic_Bacteria","Si","SLA"))
ef_3 <- envfit(select(dna_df_nometa, c("MDS1","MDS2")), meta_exclusive, permu = 999, na.rm=TRUE)
spp.scrs_2 <- as.data.frame(scores(ef_3, display = "vectors"))
spp.scrs_2<- cbind(spp.scrs_2, Samples = rownames(spp.scrs_2))

p_2 <- ggplot(dna_df_nometa) +
  geom_point(mapping = aes(x = MDS1, y = MDS2, colour = Depth)) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = spp.scrs_2,
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data = spp.scrs_2, aes(x = MDS1, y = MDS2, label = Samples),
            size = 3)

## Trying with just one depth
DCM_dna_nometa<-filter(dna_nometa, Depth=="DCM")
norm_scale$Sample = row.names(norm_scale)
norm_scale_DCM_dna <- norm_scale
norm_scale_DCM_dna <- colsplit(row.names(norm_scale_DCM_dna),"_", names=c("Station","Depth"))
norm_scale_DCM_dna <- join(norm_scale, norm_scale_DCM_dna, join_by="Sample" ,type="left")
norm_scale_DCM_dna <- filter(norm_scale,)
dcm_dna_nometa_df<-nmdsPoints(DCM_dna_nometa, key_all)
ef_dcm_dna <- envfit(select(dcm_dna_nometa_df, c("MDS1","MDS2")), norm_scale, permu = 999, na.rm=TRUE)

spp.scrs_dcm_dna <- as.data.frame(scores(ef_dcm_dna, display = "vectors"))
spp.scrs_dcm_dna <- cbind(spp.scrs_dcm_dna, Samples = rownames(spp.scrs_dcm_dna))

p_dcm_dna <- ggplot(dcm_dna_nometa_df) +
  geom_point(mapping = aes(x = MDS1, y = MDS2, colour = Station.sla)) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = spp.scrs_dcm_dna,
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data = spp.scrs_dcm_dna, aes(x = MDS1, y = MDS2, label = Samples),
            size = 3)

dna_175_nometa<-filter(dna_nometa, Depth=="175m")
dna_175_nometa_df<-nmdsPoints(dna_175_nometa, key_all)
ef_dna_175 <- envfit(select(dna_175_nometa_df, c("MDS1","MDS2")), norm_scale, permu = 999, na.rm=TRUE)

spp.scrs_dna_175 <- as.data.frame(scores(ef_dna_175, display = "vectors"))
spp.scrs_dna_175 <- cbind(spp.scrs_dna_175, Samples = rownames(spp.scrs_dna_175))

p_dna_175 <- ggplot(dna_175_nometa_df) +
  geom_point(mapping = aes(x = MDS1, y = MDS2, colour = Station.sla)) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = spp.scrs_dna_175,
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data = spp.scrs_dna_175, aes(x = MDS1, y = MDS2, label = Samples),
            size = 3)

######## 02/13/23 for density ########

dna_nometa <- filter(new_df_nometa, Material=="DNA")
dna_df_nometa<-nmdsPoints(dna_nometa, key_all)


nmds_dna_nometa<-ggplot(data=dna_df_nometa, aes(x = MDS1, y = MDS2, fill=Potential_Density,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = dna_df_nometa$Shape.depth)+ # scale_fill_manual(values=levels(dna_df_nometa$Color.sla))+
  labs(title="NMDS-Bray-Curtis: DNA-No Metazoa", x="NMDS1",y="NMDS2",size=12)+
  scico::scale_fill_scico(palette = "batlow") +
  # scale_fill_gradientn(colours=c("green", "purple"))+
  # guides(fill = guide_legend(override.aes=list(color=levels(dna_df_nometa$Color.sla))))+
  theme(legend.title = element_blank())
nmds_dna_nometa + stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_dna_nometa_density.pdf", width=8, height = 8)


rna_nometa <- filter(new_df_nometa, Material=="cDNA")
rna_df_nometa<-nmdsPoints(rna_nometa, key_all)

nmds_rna_nometa<-ggplot(data=rna_df_nometa, aes(x = MDS1, y = MDS2, fill=Potential_Density,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_df_nometa$Shape.depth)+ # scale_fill_manual(values=levels(dna_df_nometa$Color.sla))+
  labs(title="NMDS-Bray-Curtis: RNA-No Metazoa", x="NMDS1",y="NMDS2",size=12)+
  scico::scale_fill_scico(palette = "batlow") +
  # scale_fill_gradientn(colours=c("green", "purple"))+
  # guides(fill = guide_legend(override.aes=list(color=levels(dna_df_nometa$Color.sla))))+
  theme(legend.title = element_blank())
nmds_rna_nometa + stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_nometa_density0223.pdf", width=8, height = 8)


### updating with the SLA
dna_nometa <- filter(new_df_nometa, Material=="DNA")
dna_df_nometa<-nmdsPoints(dna_nometa, key_all)


nmds_dna_nometa<-ggplot(data=dna_df_nometa, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = dna_df_nometa$Shape.depth)+ # scale_fill_manual(values=levels(dna_df_nometa$Color.sla))+
  labs(title="NMDS-Bray-Curtis: DNA-No Metazoa", x="NMDS1",y="NMDS2",size=12)+
  # scico::scale_fill_scico(palette = "batlow") +
  scale_fill_gradientn(colours=c("blue", "red"))+
  # guides(fill = guide_legend(override.aes=list(color=levels(dna_df_nometa$Color.sla))))+
  theme(legend.title = element_blank())
nmds_dna_nometa + stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_dna_nometa_sla0223.pdf", width=8, height = 8)


rna_nometa <- filter(new_df_nometa, Material=="cDNA")
rna_df_nometa<-nmdsPoints(rna_nometa, key_all)

nmds_rna_nometa<-ggplot(data=rna_df_nometa, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_df_nometa$Shape.depth)+ # scale_fill_manual(values=levels(dna_df_nometa$Color.sla))+
  labs(title="NMDS-Bray-Curtis: RNA-No Metazoa", x="NMDS1",y="NMDS2",size=12)+
  # scico::scale_fill_scico(palette = "batlow") +
  scale_fill_gradientn(colours=c("blue", "red"))+
  # guides(fill = guide_legend(override.aes=list(color=levels(dna_df_nometa$Color.sla))))+
  theme(legend.title = element_blank())
nmds_rna_nometa + stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_nometa_sla0223.pdf", width=8, height = 8)

########## 2/21/23 second round RNA nmds as Dave requested #######
nometa <-filter(asv_df, Taxa!="Opisthokont-Metazoa")
Norm.rel_nometa <-decostand(nometa[,2:61], MARGIN=2, method = "total") #this does it by rows
colSums(Norm.rel_nometa) # check! should all equal 1.
Norm.rel_nometa$ASV.ID<-rownames(Norm.rel_nometa)
melt_norm_nometa<- melt(Norm.rel_nometa)
vars_nometa<-colsplit(melt_norm_nometa$variable, "_", c("Station","Depth","Material"))
new_df_nometa<-data.frame(melt_norm_nometa, vars_nometa)

rna_nometa <- filter(new_df_nometa, Material=="cDNA")
rna_df_nometa<-nmdsPoints(rna_nometa, key_all)

# 1- just the pit centers
rna_nometa_centers <- filter(rna_nometa, (Loc=="C-Center"|Loc=="AC-Center"))
names(new_df_nometa)[2] = "Sample"
rna_nometa_centers <- new_df_nometa %>%
  join(., select(key_all, -c("Station","Depth","Material")), by="Sample", type="left") %>%
  filter(., (Material=="cDNA" & (Loc=="C-Center"|Loc=="AC-Center")))
  
names(rna_nometa_centers)[2] = "variable"
rna_centers_nometa<-nmdsPoints(select(rna_nometa_centers, c("ASV.ID","variable","value","Station","Depth","Material")), key_all)

nmds_rna_nometa_centers<-ggplot(data=rna_centers_nometa, aes(x = MDS1, y = MDS2, fill=Potential_Density,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_centers_nometa$Shape.depth)+ # scale_fill_manual(values=levels(dna_df_nometa$Color.sla))+
  labs(title="NMDS-Bray-Curtis: RNA-No Metazoa", x="NMDS1",y="NMDS2",size=12)+
  scico::scale_fill_scico(palette = "batlow")
nmds_rna_nometa_centers + stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/nmds_rna_nometa_density_centers.pdf", width=8, height = 8)

nmds_rna_nometa_centers_sla<-ggplot(data=rna_centers_nometa, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_centers_nometa$Shape.depth)+ # scale_fill_manual(values=levels(dna_df_nometa$Color.sla))+
  labs(title="NMDS-Bray-Curtis: RNA-No Metazoa", x="NMDS1",y="NMDS2",size=12)+
  scale_fill_gradientn(colours=c("blue", "red"))
nmds_rna_nometa_centers_sla + stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/nmds_rna_nometa_sla_centers.pdf", width=8, height = 8)


# 2- circles around surface/500m and 175m/DCM in RNA
rna_nometa <- new_df_nometa %>%
  join(., select(key_all, -c("Station","Depth","Material")), by="Sample", type="left") %>%
  filter(., (Material=="cDNA"))
         
rna_nometa$depth_group = case_when(
  (rna_nometa$Depth == "15m" | rna_nometa$Depth == "500m") ~ "15m+500m",
  (rna_nometa$Depth == "DCM" | rna_nometa$Depth == "175m") ~ "DCM+175m",
)

names(rna_nometa)[2] = "variable"
rna_nometa_df<-nmdsPoints(select(rna_nometa, c("ASV.ID","variable","value","Station","Depth","Material")), key_all)

rna_nometa_df$depth_group = case_when(
  (rna_nometa_df$Depth == "15m" | rna_nometa_df$Depth == "500m") ~ "15m+500m",
  (rna_nometa_df$Depth == "DCM" | rna_nometa_df$Depth == "175m") ~ "DCM+175m",
)

nmds_rna_nometa_density<-ggplot(data=rna_nometa_df, aes(x = MDS1, y = MDS2, fill=Potential_Density,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_centers_nometa$Shape.depth)+ # scale_fill_manual(values=levels(dna_df_nometa$Color.sla))+
  labs(title="NMDS-Bray-Curtis: RNA-No Metazoa", x="NMDS1",y="NMDS2",size=12)+
  scico::scale_fill_scico(palette = "batlow")

nmds_rna_nometa_density + stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=depth_group))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/nmds_rna_nometa_density_othergrouping.pdf", width=8, height = 8)

nmds_rna_nometa_sla<-ggplot(data=rna_nometa_df, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_centers_nometa$Shape.depth)+ # scale_fill_manual(values=levels(dna_df_nometa$Color.sla))+
  labs(title="NMDS-Bray-Curtis: RNA-No Metazoa", x="NMDS1",y="NMDS2",size=12)+
  scale_fill_gradientn(colours=c("blue", "red"))
 
nmds_rna_nometa_sla + stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=depth_group))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/nmds_rna_nometa_sla_othergrouping.pdf", width=8, height = 8)

# 3- depths separate
rna_dcm <- filter(rna_nometa, Depth=="DCM")
rna_dcm_df<-nmdsPoints(rna_dcm, key_all)
nmds_rna_dcm_sla <-ggplot(data=rna_dcm_df, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_dcm_df$Shape.depth)+ # scale_fill_manual(values=levels(dna_df_nometa$Color.sla))+
  labs(title="DCM", x="NMDS1",y="NMDS2",size=12)+
  # scale_fill_gradientn(colours=c("blue", "red"))+
  scale_fill_gradient2(low = "#0818A8",
                       mid = "white",
                       high = "#880808",
                       midpoint = 0,
                       limits=c(-28,28))+
  scale_x_continuous(limits = c(-1.2, 1.2))+
  scale_y_continuous(limits = c(-.5, .5))
nmds_rna_dcm_sla
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_dcm_sla_2.pdf", width=8, height = 8)
nmds_rna_dcm_dense <-ggplot(data=rna_dcm_df, aes(x = MDS1, y = MDS2, fill=Potential_Density,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_dcm_df$Shape.depth)+ # scale_fill_manual(values=levels(dna_df_nometa$Color.sla))+
  labs(title="DCM", x="NMDS1",y="NMDS2",size=12)+
  scico::scale_fill_scico(limits=c(22.5,27), palette = "batlow") +
  scale_x_continuous(limits = c(-1.2, 1.2))+
  scale_y_continuous(limits = c(-.5, .5))
nmds_rna_dcm_dense
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_dcm_density.pdf", width=8, height = 8)


rna_15 <- filter(rna_nometa, Depth=="15m")
rna_15_df<-nmdsPoints(rna_15, key_all)
nmds_rna_15_sla <-ggplot(data=rna_15_df, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
  geom_point(color="black",size=6, position="jitter") + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_15_df$Shape.depth)+ # scale_fill_manual(values=levels(dna_df_nometa$Color.sla))+
  labs(title="15m", x="NMDS1",y="NMDS2",size=12)+
  scale_y_continuous(limits = c(-.5, .5))+
  scale_x_continuous(limits = c(-1.2, 1.2))+
  #scale_fill_gradientn(colours=c("blue", "red"))
  scale_fill_gradient2(low = "#0818A8",
                       mid = "white",
                       high = "#880808",
                       midpoint = 0,
                       limits=c(-28,28))
nmds_rna_15_sla
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_15m_sla_2.pdf", width=8, height = 8)
nmds_rna_15_dense <-ggplot(data=rna_15_df, aes(x = MDS1, y = MDS2, fill=Potential_Density,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_15_df$Shape.depth)+ # scale_fill_manual(values=levels(dna_df_nometa$Color.sla))+
  labs(title="15m", x="NMDS1",y="NMDS2",size=12)+
  scale_y_continuous(limits = c(-.5, .5))+
  scale_x_continuous(limits = c(-1.2, 1.2))+
  scico::scale_fill_scico(limits=c(22.5,27), palette = "batlow")
nmds_rna_15_dense
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_15m_dense.pdf", width=8, height = 8)


rna_175 <- filter(rna_nometa, Depth=="175m")
rna_175_df<-nmdsPoints(rna_175, key_all)
nmds_rna_175_sla <-ggplot(data=rna_175_df, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_175_df$Shape.depth)+ # scale_fill_manual(values=levels(dna_df_nometa$Color.sla))+
  labs(title="175m", x="NMDS1",y="NMDS2",size=12)+
  scale_y_continuous(limits = c(-.5, .5))+
  scale_x_continuous(limits = c(-1.2, 1.2))+
  scale_fill_gradient2(low = "#0818A8",
                       mid = "white",
                       high = "#880808",
                       midpoint = 0,
                       limits=c(-28,28))
  # scale_fill_gradientn(colours=c("blue", "red"))
nmds_rna_175_sla
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_175m_sla_2.pdf", width=8, height = 8)
nmds_rna_175_dense <-ggplot(data=rna_175_df, aes(x = MDS1, y = MDS2, fill=Potential_Density,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_175_df$Shape.depth)+ # scale_fill_manual(values=levels(dna_df_nometa$Color.sla))+
  labs(title="175m", x="NMDS1",y="NMDS2",size=12)+
  scale_y_continuous(limits = c(-.5, .5))+
  scale_x_continuous(limits = c(-1.2, 1.2))+
  scico::scale_fill_scico(limits=c(22.5,27), palette = "batlow")
nmds_rna_175_dense
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_175m_dense.pdf", width=8, height = 8)

rna_500 <- filter(rna_nometa, Depth=="500m")
rna_500_df<-nmdsPoints(rna_500, key_all)
nmds_rna_500_sla <-ggplot(data=rna_500_df, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_500_df$Shape.depth)+ # scale_fill_manual(values=levels(dna_df_nometa$Color.sla))+
  labs(title="500m", x="NMDS1",y="NMDS2",size=12)+
  scale_y_continuous(limits = c(-.5, .5))+
  scale_x_continuous(limits = c(-1.2, 1.2))+
  scale_fill_gradient2(low = "#0818A8",
                       mid = "white",
                       high = "#880808",
                       midpoint = 0,
                       limits=c(-28,28))
  # scale_fill_gradientn(colours=c("blue", "red"))
nmds_rna_500_sla
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_500m_sla_2.pdf", width=8, height = 8)
nmds_rna_500_dense <-ggplot(data=rna_500_df, aes(x = MDS1, y = MDS2, fill=Potential_Density,shape=Depth)) + 
  geom_point(color="black",size=6) + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = rna_500_df$Shape.depth)+ # scale_fill_manual(values=levels(dna_df_nometa$Color.sla))+
  labs(title="500m", x="NMDS1",y="NMDS2",size=12)+ 
  scale_y_continuous(limits = c(-.5, .5))+
  scale_x_continuous(limits = c(-1.2, 1.2))+
  scico::scale_fill_scico(limits=c(22.5,27),palette = "batlow")
nmds_rna_500_dense
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_500m_dense.pdf", width=8, height = 8)

ggarrange(nmds_rna_15_sla,nmds_rna_dcm_sla,nmds_rna_175_sla,nmds_rna_500_sla)
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_depths_sla.pdf", width=12, height = 12)

ggarrange(nmds_rna_15_dense,nmds_rna_dcm_dense,nmds_rna_175_dense,nmds_rna_500_dense)
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_depths_dense.pdf", width=10, height = 10)
library(patchwork)

nmds_rna_15_dense+nmds_rna_dcm_dense +nmds_rna_175_dense + nmds_rna_500_dense + plot_layout(guides="collect") +  plot_annotation(title="NMDS-Bray-Curtis: RNA-No Metazoa")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_depths_dense_2.pdf", width=10, height = 10)

nmds_rna_15_sla + nmds_rna_dcm_sla + nmds_rna_175_sla + nmds_rna_500_sla + plot_layout(guides="collect") +  plot_annotation(title="NMDS-Bray-Curtis: RNA-No Metazoa")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_rna_depths_sla_3.pdf", width=10, height = 10)

######## 11. Alpha + Shannon's Diveristy #########
col_material<- c("#7fcdbb","#2c7fb8","#fc8d59")
# Needs to transposeee
shannons_df<- asv_df[,2:61] %>%
  t(.)%>%
  diversity(., index = "shannon",MARGIN = 1) %>%
  melt() %>%
  mutate(Sample = rownames(.))%>%
  join(.,key_all, by="Sample") %>%
  arrange(., Material, Station.sla, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))
names(shannons_df)[1]="H"
rownames(shannons_df)=shannons_df$Sample
head(shannons_df)

ggplot(shannons_df, aes(x = Sample, y = H,shape=Material,fill=SLA)) + 
  geom_point(color="black",size=6) + 
  theme_bw() + 
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_fill_gradientn(colours=c("blue", "red"))+
  scale_shape_manual(values = c(21,25,22)) + 
  labs(title="Shannon's Diversity Index- Metazoa Included", x="",y="Shannon's Diversity Index")+
  theme(legend.position="right",axis.text.x = element_text(angle=90, hjust=1,vjust=1,size=10,color="black"))

ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/shannons_All_0323.pdf", width=8, height = 8)


# box plots
ggplot(shannons_df, aes(x=Material, y=H, fill=Material))+
  geom_boxplot(position = position_dodge(preserve = "single"))+
  scale_fill_manual(values=col_material)+
  # facet_grid(variable~.,scales="free")+
  theme_bw()+
  labs(title="Shannon's Diversity Index", y="Shannon's Index")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size=8),axis.text.y=element_text(size=12),legend.position = "right")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/shannons_boxplot_0323.pdf", width=6, height = 6)
# centers only
shannons_centers= filter(shannons_df, (Station=="S6" | Station=="S12"))
ggplot(shannons_centers, aes(x=Station, y=H, fill=Material))+
  geom_boxplot(position = position_dodge(preserve = "single"))+
  scale_fill_manual(values=col_material)+
  facet_grid(Material~.,scales="free")+
  theme_bw()+
  ylim(5.4,6.5)+
  labs(title="Shannon's Diversity Index", y="Shannon's Index")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size=8),axis.text.y=element_text(size=12),legend.position = "right")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/shannons_boxplot_0323_centers.pdf", width=6, height = 6)

# Now without metazoa
shannons_df_nometa<- asv_df %>%
  filter(.,Taxa!="Opisthokont-Metazoa")%>%
  .[,2:61]%>%
  t(.)%>%
  diversity(., index = "shannon",MARGIN = 1) %>%
  melt() %>%
  mutate(Sample = rownames(.))%>%
  join(.,key_all, by="Sample") %>%
  arrange(., Material, Station.sla, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))
names(shannons_df_nometa)[1]="H"
rownames(shannons_df_nometa)=shannons_df_nometa$Sample
head(shannons_df_nometa)

ggplot(shannons_df_nometa, aes(x = Sample, y = H,shape=Material,fill=SLA)) + 
  geom_point(color="black",size=6) + 
  theme_bw() + 
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_fill_gradientn(colours=c("blue", "red"))+
  scale_shape_manual(values = c(21,25,22)) + 
  labs(title="Shannon's Diversity Index- No Metazoa", x="",y="Shannon's Diversity Index")+
  theme(legend.position="right",axis.text.x = element_text(angle=90, hjust=1,vjust=1,size=10,color="black"))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/shannons_nometa_0323.pdf", width=8, height = 8)

ggplot(shannons_df_nometa, aes(x = Sample, y = H,shape=Depth,fill=Material)) + 
  geom_point(color="black",size=6) + 
  theme_bw() + 
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_fill_manual(values=c("#7fcdbb","#2c7fb8","#fc8d59"))+
  scale_shape_manual(values = c(21,25,22,23,24)) + 
  labs(title="Shannon's Diversity Index- No Metazoa", x="",y="Shannon's Diversity Index")+
  theme(legend.position="right",axis.text.x = element_text(angle=90, hjust=1,vjust=1,size=10,color="black"))
# box plots
ggplot(shannons_df_nometa, aes(x=Material, y=H, fill=Material))+
  geom_boxplot(position = position_dodge(preserve = "single"))+
  scale_fill_manual(values=col_material)+
  # facet_grid(variable~.,scales="free")+
  theme_bw()+
  labs(title="Shannon's Diversity Index- No metazoa", y="Shannon's Index")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size=8),axis.text.y=element_text(size=12),legend.position = "right")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/shannons_boxplot_0323_nometa.pdf", width=6, height = 6)

# centers only
shannons_centers_nometa= filter(shannons_df_nometa, (Station=="S6" | Station=="S12"))
ggplot(shannons_centers_nometa, aes(x=Station, y=H, fill=Material))+
  geom_boxplot(position = position_dodge(preserve = "single"))+
  scale_fill_manual(values=col_material)+
  facet_grid(Material~.,scales="free")+
  theme_bw()+
  ylim(5.4,6.5)+
  labs(title="Shannon's Diversity Index-No metazoa", y="Shannon's Index")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size=8),axis.text.y=element_text(size=12),legend.position = "right")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/shannons_boxplot_0323_centers_nometa.pdf", width=6, height = 6)


#compare alpha diversity (from Sarah Hu:

shannon<-diversity(asv_df[,2:61], index = "shannon",MARGIN = 2)
invsimp<-diversity(asv_df[,2:61], index = "invsimpson",MARGIN = 2)
ASV_count<-colSums(asv_df[,2:61]>0)

alpha<-data.frame(shannon, invsimp, ASV_count)
head(alpha)

alpha$Sample<-row.names(alpha)
alpha.df<- alpha %>%
  melt(.) %>%
  join(.,key_all, by="Sample", type="left", match="first") %>%
  arrange(., Material, Station.sla, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))
head(alpha.df)

ggplot(alpha.df, aes(x=Material, y=value, fill=Material))+
  geom_boxplot(position = position_dodge(preserve = "single"))+
  scale_fill_manual(values=col_material)+
  facet_grid(variable~.,scales="free")+
  theme_bw()+
  labs(title="Alpha Diversity Metrics")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size=8),axis.text.y=element_text(size=12),legend.position = "right")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/aloha_boxplot_0323.pdf", width=6, height = 6)

#now no meta
asv_df_nometa <-asv_df %>%
  filter(.,Taxa!="Opisthokont-Metazoa")
shannon_nometa<-diversity(asv_df_nometa[,2:61], index = "shannon",MARGIN = 2)
invsimp_nometa<-diversity(asv_df_nometa[,2:61], index = "invsimpson",MARGIN = 2)
ASV_count_nometa<-colSums(asv_df_nometa[,2:61]>0)
ASV_number_nometa<-colSums(asv_df_nometa[,2:61])
Num_Seqs_stats <- summary(ASV_number_nometa)
num_seqs_sd <- sd(ASV_number_nometa)

alpha_nometa<-data.frame(ASV_number_nometa, ASV_count_nometa, shannon_nometa, invsimp_nometa)
head(alpha_nometa)

alpha_nometa$Sample<-row.names(alpha_nometa)
alpha.df_nometa<- alpha_nometa %>%
  melt(.) %>%
  join(.,key_all, by="Sample", type="left", match="first")
head(alpha.df_nometa)

ggplot(alpha.df_nometa, aes(x=Sample.type, y=value, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material)+
  facet_grid(variable~.,scales="free")+
  theme_bw()+
  # labs(title="Alpha Diversity Metrics-No Metazoa")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size=8),axis.text.y=element_text(size=12),legend.position = "right")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/alpha_boxplot_0323_nometa.pdf", width=4, height = 6)

seq_num <-ggplot(alpha.df_nometa[alpha.df_nometa$variable=="ASV_number_nometa",], aes(x=Sample.type, y=value, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material, name="Sample Type")+
  theme_bw()+
  labs(title="", x="", y="Number of Sequences")+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size=8),axis.text.y=element_text(size=12),legend.position = "right")

alpha<-ggplot(alpha.df_nometa[alpha.df_nometa$variable=="ASV_count_nometa",], aes(x=Sample.type, y=value, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material, name="Sample Type")+
  theme_bw()+
  labs(title="", x="", y="Number of ASVs")+
  scale_y_continuous(breaks=seq(0,3000,1000), limits=c(0,3000)) +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size=8),axis.text.y=element_text(size=12),legend.position = "right")

shannons<-ggplot(alpha.df_nometa[alpha.df_nometa$variable=="shannon_nometa",], aes(x=Sample.type, y=value, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material, name="Sample Type")+
  theme_bw()+
  labs(title="", x="", y="Shannon's Index")+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size=8),axis.text.y=element_text(size=12),legend.position = "right")

inv_simp<-ggplot(alpha.df_nometa[alpha.df_nometa$variable=="invsimp_nometa",], aes(x=Sample.type, y=value, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material, name="Sample Type")+
  theme_bw()+
  labs(title="", x="", y="Inverse Simpon's")+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size=8),axis.text.y=element_text(size=12),legend.position = "right")

(seq_num / alpha / shannons / inv_simp) + plot_layout(guides="collect") + plot_annotation(tag_levels="a")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/diversity_boxplot_0323_nometa.pdf", width=6, height = 10)

######### 12. Make the PITs and microscopy figs look alike #######
# (Edited from Full_PIT_Microscope_0521.R)
counts<-read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/CellsPit_Final.csv",header=T)
rownames(counts)=counts$Type
# combine ciliates and tintinnid to the same category
counts$Plotting = counts$Type
counts$Plotting[counts$Type == "Tintinnid"]="Ciliate"
# combine acanths and rads to one category
counts$Plotting[counts$Type == "Rhizaria-Acantharia"]="Rhizaria-Radiolaria"
# make the color match the other type
tax_order.m<-c("Ciliate","Dinoflagellate", "Diatom","Dictyocha", "Rhizaria-Foraminifera","Rhizaria-Radiolaria","Rhizaria-other")
cols.m<-c('firebrick4','indianred1','yellowgreen','lightblue','mediumvioletred','magenta','lightpink')

counts.m<-melt(counts, id.vars = c("Type","Plotting"))
head(counts.m)
df<-colsplit(counts.m$variable, "_", c("Station","Rep"))
counts.ms<-data.frame(counts.m,df)
head(counts.ms)

#To calculate average values of taxa by PIT
counts.mean<- counts.ms %>% 
  group_by(Type, Station) %>%
  summarise(mean=mean(value),stdev=sd(value)) %>%
  as.data.frame
counts.mean

# bring back my plotting column
counts.mean.df <- join(counts.mean, select(counts, c("Type","Plotting")), by="Type", type="left")

# need to add the columns
counts_plot <- counts.mean.df %>%
  group_by(Plotting, Station) %>%
  summarise(sum = sum(mean)) %>%
  as.data.frame
counts_plot
counts_plot$tax.order<-factor(counts_plot$Plotting,levels=(tax_order.m), labels=(tax_order.m))
counts_plot$Station <- factor(counts_plot$Station, levels=c("PIT12","PIT11","PIT10","PIT9","PIT8","PIT7","PIT6","PIT5","PIT4","PIT3","PIT2","PIT1"))

#for relative abundance based on mean
data.agg.count<-aggregate(counts_plot$sum, by=list(Taxa=counts_plot$Plotting,Sample=counts_plot$Station),sum) #sum sequences by taxonomic group
data.agg.count$Sample<-factor(data.agg.count$Sample,levels=c("PIT12","PIT11","PIT10","PIT9","PIT8","PIT7","PIT6","PIT5","PIT4","PIT3","PIT2","PIT1"))
data.agg.count$tax.order<-factor(data.agg.count$Taxa,levels=(tax_order.m), labels=(tax_order.m))

barplot<-ggplot(counts_plot, aes(y=sum, fill=tax.order, x=Station))+
  geom_bar(position = "fill",stat="identity",color="black",aes(fill=tax.order))+
  scale_fill_manual(values=cols.m)+
  theme_bw()+
  labs(title="Relative abundance of mean counts per PIT", x="",y="")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5))
barplot
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/MicroscopyCellMean_PIT_relav.pdf", width=8, height = 6)

#Now to restructure the ASV PIT data to look like the microscopy data
#no meta sequences
nometa<-filter(asv_df, Level3 != "Metazoa")
data.m<-melt(nometa) #melt
head(data.m)
names(data.m)[14]="Sample"
names(data.m)
data.m.key <- join(data.m, key_all, by="Sample", type="left")
head(data.m.key)
data.m.pit <- subset(data.m.key, Depth=="PIT")
# Rename so it matches the microscopy taxonomy
data.m.pit$Plotting<-"Other"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Alveolates-Ciliates"]="Ciliate"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Alveolates-Dinophyceae"]="Dinoflagellate"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Stramenopiles-Diatoms"]="Diatom"
data.m.pit$Plotting[data.m.pit$Level4 == "Dictyochophyceae"]="Dictyocha"
# combine acanths and rads to one category
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Rhizaria-Polycystines"]="Rhizaria-Radiolaria"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Rhizaria-Acantharia"]="Rhizaria-Radiolaria"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Rhizaria-Cercozoa"] = "Rhizaria-other"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Rhizaria-Other"] = "Rhizaria-other"
tax_order.pit<-c("Ciliate","Dinoflagellate", "Diatom","Dictyocha", "Rhizaria-Foraminifera","Rhizaria-Radiolaria","Rhizaria-other", "Other")
col.pit <-c('firebrick4','indianred1','yellowgreen','lightblue','magenta','lightpink','lightgrey')

data.agg.pit<-aggregate(data.m.pit$value, by=list(Taxa=data.m.pit$Plotting,Sample=data.m.pit$Station),sum) #sum sequences by taxonomic group
data.agg.pit$tax.order<-factor(data.agg.pit$Taxa,levels=(tax_order.pit), labels=(tax_order.pit))


#Bar plot of community composition ordered by Material, Station and Depth
asv_barplot<-ggplot(data.agg.pit, aes(y=x, fill=tax.order, x=Sample))+
  geom_bar(position = "fill",stat="identity",color="black",aes(fill=tax.order))+
  scale_fill_manual(values=col.pit)+
  theme_bw()+
  labs(title="Relative abundance of ASVs", x="",y="")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5))
asv_barplot
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/ASVPIT_relav.pdf", width=8, height = 6)
asv_pit_bars <- asv_barplot + coord_flip()

mico_pit_bars <- barplot +  coord_flip()

ggarrange(asv_pit_bars, mico_pit_bars)
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/PIT_relav_combo.pdf", width=10, height = 4)

# looks weird so will plot everything
data.agg.pit2<-aggregate(data.m.pit$value, by=list(Taxa=data.m.pit$TaxaPlot,Sample=data.m.pit$Sample),sum) #sum sequences by taxonomic group

tax_order<-c("Alveolates-Ciliates","Alveolates-Dinophyceae","Alveolates-Syndiniales", "Archaeplastids-Chlorophytes",
             "Excavates-Discoba", "Hacrobia-Cryptophytes","Hacrobia-Haptophytes",
             "Opisthokont-Choanoflagellida","Opisthokonts-Metazoa","Other/unknown",
             "Rhizaria-Acantharia","Rhizaria-Cercozoa","Rhizaria-Polycystines", "Rhizaria-Other",
             "Stramenopiles-Chrysophytes","Stramenopiles-Diatoms","Stramenopiles-MAST","Stramenopiles-Ochrophyta","Stramenopiles-Pelagophytes","Stramenopiles-Other")
tax_color_phyla2<-c('firebrick4','indianred1','tomato3','forestgreen','yellowgreen','darkblue',
                    'lightblue','moccasin','gold1',"light grey","#C291A4",'mediumvioletred','magenta','lightpink',"#DECDBE",'#DDAD4B','tan2','tan3','tan4',"#5C4033")
names(tax_color_phyla2)<-tax_order
bar_plot_df <- join(data.agg, key_all, by="Sample", type="left", match="first") # %>%
arrange(., Material, Station, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))

#Bar plot of community composition ordered by Material, Station and Depth
bars.pit<-ggplot(data.agg.pit2, aes(y=x, fill=tax, x=Sample))+
  geom_bar(position = "fill",stat="identity",color="black",aes(fill=Taxa))+
  scale_fill_manual(values=tax_color_phyla2)+
  theme_bw()+
  labs(title="Relative abundance of reads", x="Sample",y="")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5, vjust = .5))

bars.pit

########## 13. Microscopy flux plots ########

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
flux$Plotting[flux$Type == "Diatom"]="Stramenopiles-Diatom"
flux$Plotting[flux$Type == "Dictyocha"]="Stramenopiles-Dictyocha"

flux.m<-melt(select(flux, -c("Type")), id.vars = "Plotting")
df.f<-colsplit(flux.m$variable, "_", c("Station","Rep"))
flux.ms<-data.frame(flux.m,df.f)
head(flux.ms)


#To calculate mean values of taxa by PIT
flux.mean<- flux.ms %>% 
  group_by(Plotting, Station) %>%
  summarise(mean=mean(as.numeric(value)),stdev=sd(as.numeric(value))) %>%
  as.data.frame
flux.mean


# make the color match the other type
tax_order.pit<-c("Alveolates-Ciliates","Alveolates-Dinoflagellates","Other","Rhizaria-Foraminifera","Rhizaria-Radiolaria","Rhizaria-Other", "Stramenopiles-Diatom","Stramenopiles-Dictyocha")
col.pit <-c('firebrick4','indianred1', 'grey','mediumvioletred',"#C291A4",'lightpink','#DDAD4B',"#5C4033")


flux.mean$tax.order<-factor(flux.mean$Plotting,levels=(tax_order.pit), labels=(tax_order.pit))
flux.mean$Station<-factor(flux.mean$Station, levels=c("PIT12","PIT11","PIT10","PIT9","PIT8","PIT7","PIT6","PIT5","PIT4","PIT3","PIT2","PIT1"))

pit_colors <- c("Alveolates-Ciliates"='firebrick4',"Alveolates-Dinoflagellates"='indianred1',"Other"= 'grey',
                  "Rhizaria-Foraminifera"='mediumvioletred',"Rhizaria-Radiolaria"='lightpink',"Rhizaria-Other"='#C291A4',
                  "Stramenopiles-Diatom"='#DDAD4B',"Stramenopiles-Dictyocha"="#5C4033")


#To plot flux values
scatter<-ggplot(flux.mean, aes(y=mean,ymin= mean-stdev, ymax= mean + stdev, x=Station,fill=tax.order))+
  geom_errorbar(width = 0.2) +
  geom_point(size=3.5, shape=21,color="black")+
  scale_fill_manual(values=pit_colors, name="Taxa")+
  scale_y_continuous(trans = 'log10',labels = function(x) format(x, scientific = TRUE))+
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
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5))
scatter %+% facet_wrap(~tax.order, ncol = 2) 
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/mean_flux_facet_stdev_samescale_0323.pdf", width=6.5, height = 7.5)
scatter %+% facet_wrap(~tax.order, ncol = 2, scales="free_y") 
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/mean_flux_facet_stdev_difscale_0323.pdf", width=6.5, height = 7.5)

scatter
ggplot2::ggsave("mean_flux_stdev.pdf", width=5, height = 5)
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/mean_flux_stdev.pdf", width=10, height = 4)

######### 14. UpSetR plots D:< ####

######## Format for UpSetR plots:
# upsetR
melt_ASV<-melt(asv_df)
names(melt_ASV)[14]="Sample"
melt_all<-join(melt_ASV,key_all,by="Sample", type="left", match="first")
head(melt_all)
toUpset<-melt_all[c(1:2,14:18)] # We want the ASV ID, Sample, TaxaPlot, value, Station, Depth, Material
head(toUpset) # Use below to add annotation
toUpset$Uniq<-paste(toUpset$ASV1.ID, toUpset$TaxaPlot, sep="_")

# Change to binary
summed_material<-aggregate(toUpset$value, by=list(Material=toUpset$Material, Uniq=toUpset$Uniq),sum)
summed_material$bin<-ifelse(summed_material$x > 0, 1, 0)
binary_wide_material<-dcast(summed_material[c(1,2,4)], Uniq~Material, fill=0)
row.names(binary_wide_material)<-binary_wide_material$Uniq; binary_wide_material$Uniq<-NULL
head(binary_wide_material)

# Let's make "Intersect" column to tell us which months a given ASV shows up in 
binary_wide_material$Intersect <- apply(binary_wide_material > 0, 1, function(x){toString(names(binary_wide_material)[x])})
head(binary_wide_material)

# Get rid of spaces after common in "Intersect" column
binary_wide_material$Intersect <- stri_replace_all_fixed(binary_wide_material$Intersect, " ", "")
head(binary_wide_material)
# Convert "Intersect" column to list format
binary_wide_material$Intersect <- strsplit(binary_wide_material$Intersect, ",")
head(binary_wide_material)




binary_tax_material <- binary_wide_material %>%
  rownames_to_column(var = "uniq") %>% 
  separate(uniq, c("ASV.ID", "Taxonomy"), sep = "_", remove = FALSE) %>% 
  # inner_join(tmp, by = "Taxonomy") %>% 
  column_to_rownames(var = "uniq") %>% 
  data.frame
head(binary_tax_material[1:2,])
# unique(binary_tax_material$tax_compiled)

tax_order<-c("Alveolates-Ciliates","Alveolates-Dinophyceae","Alveolates-Syndiniales", "Archaeplastids-Chlorophytes",
             "Excavates-Discoba", "Hacrobia-Cryptophytes","Hacrobia-Haptophytes",
             "Opisthokont-Choanoflagellida","Opisthokonts-Metazoa","Other/unknown",
             "Rhizaria-Acantharia","Rhizaria-Cercozoa","Rhizaria-Polycystines", "Rhizaria-Other",
             "Stramenopiles-Chrysophytes","Stramenopiles-Diatoms","Stramenopiles-MAST","Stramenopiles-Ochrophyta","Stramenopiles-Pelagophytes","Stramenopiles-Other")
tax_color_phyla2<-c('firebrick4','indianred1','tomato3','forestgreen','yellowgreen','darkblue',
                    'lightblue','moccasin','gold1','grey',"#C291A4",'mediumvioletred','magenta','lightpink',"#DECDBE",'#DDAD4B','tan2','tan3','tan4',"#5C4033")
names(tax_color_phyla2)<-tax_order
binary_tax_material$TAX_ORDER<-factor(binary_tax_material$Taxonomy, levels = (tax_order), labels = (tax_order))

totals_material<- binary_tax_material %>%
  group_by(Intersect) %>%
  summarise(total=sum(length(ASV.ID))) %>%
  data.frame
write.csv(as.data.frame(totals_material), "~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/totals_material.csv")
capture.output(summary(totals_material), file = "~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/totals_material.txt")
lapply(totals_material, function(x) write.table( data.frame(x), "~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/totals_material.csv"  , append= T, sep=',' ))

totals.num_material<-order(totals_material$total,decreasing=TRUE)
total.data_material<-data.frame(as.list(totals_material$Intersect),totals_material$total)

binary_tax_material%>% ggplot(aes(x=Intersect, fill = Taxonomy)) +
  geom_bar(stat = "count", position="stack", color="black") + 
  # geom_text(aes(y=sum(stat(count)), group=Intersect), vjust=0)+
  geom_text(data=totals_material ,aes(x=Intersect,y=total,label=total,fill=NULL),
            nudge_y = 10)+
  scale_x_upset()+scale_fill_manual(name="Taxonomic Group",values=c(tax_color_phyla2))+
  theme_classic(base_size = 12)+xlab("")+ylab("Number of Shared ASVs")+
  theme(plot.margin = margin(10, 10, 10, 100))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/UpSetR_all_material_0323.pdf", width=12, height = 10)

### Now doing it by depth
# Change to binary
summed_depth<-aggregate(toUpset$value, by=list(Depth=toUpset$Depth, Uniq=toUpset$Uniq),sum)
summed_depth$bin<-ifelse(summed_depth$x > 0, 1, 0)
binary_wide_depth<-dcast(summed_depth[c(1,2,4)], Uniq~Depth, fill=0)
row.names(binary_wide_depth)<-binary_wide_depth$Uniq; binary_wide_depth$Uniq<-NULL
head(binary_wide_depth)

# Let's make "Intersect" column to tell us which months a given ASV shows up in 
binary_wide_depth$Intersect <- apply(binary_wide_depth > 0, 1, function(x){toString(names(binary_wide_depth)[x])})
head(binary_wide_depth)
# Get rid of spaces after common in "Intersect" column
binary_wide_depth$Intersect <- stri_replace_all_fixed(binary_wide_depth$Intersect, " ", "")
head(binary_wide_depth)
# Convert "Intersect" column to list format
binary_wide_depth$Intersect <- strsplit(binary_wide_depth$Intersect, ",")
head(binary_wide_depth)

binary_tax_depth <- binary_wide_depth %>%
  rownames_to_column(var = "uniq") %>% 
  separate(uniq, c("ASV.ID", "Taxonomy"), sep = "_", remove = FALSE) %>% 
  column_to_rownames(var = "uniq") %>% 
  data.frame
head(binary_tax_depth[1:2,])

totals_depth<- binary_tax_depth %>%
  group_by(Intersect) %>%
  summarise(total=sum(length(ASV.ID))) 
totals.num_depth<-order(totals_depth$total,decreasing=TRUE)
total.data_material<-data.frame(as.list(totals_material$Intersect),totals_material$total)

binary_tax_depth%>% ggplot(aes(x=Intersect, fill = Taxonomy)) +
  geom_bar(stat = "count", position="stack", color="black") + 
  # geom_text(aes(y=sum(stat(count)), group=Intersect), vjust=0)+
  geom_text(data=totals_depth ,aes(x=Intersect,y=total,label=total,fill=NULL),
            nudge_y = 10)+
  scale_x_upset()+scale_fill_manual(name="Taxonomic Group",values=c(tax_color_phyla2))+
  theme_classic(base_size = 12)+xlab("")+ylab("Number of Shared ASVs")+
  theme(plot.margin = margin(10, 10, 10, 100))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/UpSetR_all_depth_0323.pdf", width=12, height = 10)

#### lets do depth and material
# Change to binary
summed_depth_mat<-aggregate(toUpset$value, by=list(Material=toUpset$Material,Depth=toUpset$Depth, Uniq=toUpset$Uniq),sum)
summed_depth_mat$bin<-ifelse(summed_depth_mat$x > 0, 1, 0)
binary_wide_depth_mat<-dcast(summed_depth_mat[c(1,2,3,5)], Uniq~Material+Depth, fill=0)
row.names(binary_wide_depth_mat)<-binary_wide_depth_mat$Uniq; binary_wide_depth_mat$Uniq<-NULL
head(binary_wide_depth_mat)

# Let's make "Intersect" column to tell us which months a given ASV shows up in 
binary_wide_depth_mat$Intersect <- apply(binary_wide_depth_mat > 0, 1, function(x){toString(names(binary_wide_depth_mat)[x])})
head(binary_wide_depth_mat)
# Get rid of spaces after common in "Intersect" column
binary_wide_depth_mat$Intersect <- stri_replace_all_fixed(binary_wide_depth_mat$Intersect, " ", "")
head(binary_wide_depth_mat)
# Convert "Intersect" column to list format
binary_wide_depth_mat$Intersect <- strsplit(binary_wide_depth_mat$Intersect, ",")
head(binary_wide_depth_mat)

binary_tax_depth_mat <- binary_wide_depth_mat %>%
  rownames_to_column(var = "uniq") %>% 
  separate(uniq, c("ASV.ID", "Taxonomy"), sep = "_", remove = FALSE) %>% 
  column_to_rownames(var = "uniq") %>% 
  data.frame
head(binary_tax_depth_mat[1:2,])

totals_depth_mat<- binary_tax_depth_mat %>%
  group_by(Intersect) %>%
  summarise(total=sum(length(ASV.ID))) 
totals.num_depth_mat<-order(totals_depth_mat$total,decreasing=TRUE)
total.data_material<-data.frame(as.list(totals_material$Intersect),totals_material$total)

binary_tax_depth_mat%>% ggplot(aes(x=Intersect, fill = Taxonomy)) +
  geom_bar(stat = "count", position="stack", color="black") + 
  # geom_text(aes(y=sum(stat(count)), group=Intersect), vjust=0)+
  geom_text(data=totals_depth_mat ,aes(x=Intersect,y=total,label=total,fill=NULL),
            nudge_y = 10)+
  scale_x_upset()+scale_fill_manual(name="Taxonomic Group",values=c(tax_color_phyla2))+
  theme_classic(base_size = 12)+xlab("")+ylab("Number of Shared ASVs")+
  theme(plot.margin = margin(10, 10, 10, 100))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/UpSetR_all_depth_material0323.pdf", width=12, height = 10)

### Now doing it by station
# Change to binary
summed_stn<-aggregate(toUpset$value, by=list(Station=toUpset$Station, Uniq=toUpset$Uniq),sum)
summed_stn$bin<-ifelse(summed_stn$x > 0, 1, 0)
binary_wide_stn<-dcast(summed_stn[c(1,2,4)], Uniq~Station, fill=0)
row.names(binary_wide_stn)<-binary_wide_stn$Uniq; binary_wide_stn$Uniq<-NULL
head(binary_wide_stn)

# Let's make "Intersect" column to tell us which months a given ASV shows up in 
binary_wide_stn$Intersect <- apply(binary_wide_stn > 0, 1, function(x){toString(names(binary_wide_stn)[x])})
head(binary_wide_stn)
# Get rid of spaces after common in "Intersect" column
binary_wide_stn$Intersect <- stri_replace_all_fixed(binary_wide_stn$Intersect, " ", "")
head(binary_wide_stn)
# Convert "Intersect" column to list format
binary_wide_stn$Intersect <- strsplit(binary_wide_stn$Intersect, ",")
head(binary_wide_stn)

binary_tax_stn <- binary_wide_stn %>%
  rownames_to_column(var = "uniq") %>% 
  separate(uniq, c("ASV.ID", "Taxonomy"), sep = "_", remove = FALSE) %>% 
  column_to_rownames(var = "uniq") %>% 
  data.frame
head(binary_tax_stn[1:2,])


totals_stn<- binary_tax_stn %>%
  group_by(Intersect) %>%
  summarise(total=sum(length(ASV.ID))) 
totals.num_stn<-order(totals_stn$total,decreasing=TRUE)
total.data_material<-data.frame(as.list(totals_stn$Intersect),totals_stn$total)

binary_tax_stn%>% ggplot(aes(x=Intersect, fill = Taxonomy)) +
  geom_bar(stat = "count", position="stack", color="black") + 
  # geom_text(aes(y=sum(stat(count)), group=Intersect), vjust=0)+
  geom_text(data=totals_stn ,aes(x=Intersect,y=total,label=total,fill=NULL),
            nudge_y = 10)+
  scale_x_upset()+scale_fill_manual(name="Taxonomic Group",values=c(tax_color_phyla2))+
  theme_classic(base_size = 12)+xlab("")+ylab("Number of Shared ASVs")+
  theme(plot.margin = margin(10, 10, 10, 100))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/UpSetR_all_station_0323.pdf", width=12, height = 10)
# that didn't work, so lets look at just the wc 

### Now doing it just water column
# Change to binary
WC = c("DNA", "RNA")
toUpset_wc <- filter(toUpset, Material %in% WC)

summed_stn_wc<-aggregate(toUpset_wc$value, by=list(Station=toUpset_wc$Station, Uniq=toUpset_wc$Uniq),sum)
summed_stn_wc$bin<-ifelse(summed_stn_wc$x > 0, 1, 0)
binary_wide_stn_wc<-dcast(summed_stn_wc[c(1,2,4)], Uniq~Station, fill=0)
row.names(binary_wide_stn_wc)<-binary_wide_stn_wc$Uniq; binary_wide_stn_wc$Uniq<-NULL
head(binary_wide_stn_wc)

# Let's make "Intersect" column to tell us which months a given ASV shows up in 
binary_wide_stn_wc$Intersect <- apply(binary_wide_stn_wc > 0, 1, function(x){toString(names(binary_wide_stn_wc)[x])})
head(binary_wide_stn_wc)
# Get rid of spaces after common in "Intersect" column
binary_wide_stn_wc$Intersect <- stri_replace_all_fixed(binary_wide_stn_wc$Intersect, " ", "")
head(binary_wide_stn_wc)
# Convert "Intersect" column to list format
binary_wide_stn_wc$Intersect <- strsplit(binary_wide_stn_wc$Intersect, ",")
head(binary_wide_stn_wc)

binary_tax_stn_wc <- binary_wide_stn_wc %>%
  rownames_to_column(var = "uniq") %>% 
  separate(uniq, c("ASV.ID", "Taxonomy"), sep = "_", remove = FALSE) %>% 
  column_to_rownames(var = "uniq") %>% 
  data.frame
head(binary_tax_stn_wc[1:2,])


totals_stn_wc<- binary_tax_stn_wc %>%
  group_by(Intersect) %>%
  summarise(total=sum(length(ASV.ID))) 
totals.num_stn_wc<-order(totals_stn_wc$total,decreasing=TRUE)
total.data_material<-data.frame(as.list(totals_stn$Intersect),totals_stn$total)

binary_tax_stn_wc%>% ggplot(aes(x=Intersect, fill = Taxonomy)) +
  geom_bar(stat = "count", position="stack", color="black") + 
  # geom_text(aes(y=sum(stat(count)), group=Intersect), vjust=0)+
  geom_text(data=totals_stn_wc ,aes(x=Intersect,y=total,label=total,fill=NULL),
            nudge_y = 10)+
  scale_x_upset()+scale_fill_manual(name="Taxonomic Group",values=c(tax_color_phyla2))+
  theme_classic(base_size = 12)+xlab("")+ylab("Number of Shared ASVs")+
  theme(plot.margin = margin(10, 10, 10, 100))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/UpSetR_WC_station_0323.pdf", width=12, height = 10)

# depth and water column
# Change to binary

summed_depth_wc<-aggregate(toUpset_wc$value, by=list(Depth=toUpset_wc$Depth, Uniq=toUpset_wc$Uniq),sum)
summed_depth_wc$bin<-ifelse(summed_depth_wc$x > 0, 1, 0)
binary_wide_depth_wc<-dcast(summed_depth_wc[c(1,2,4)], Uniq~Depth, fill=0)
row.names(binary_wide_depth_wc)<-binary_wide_depth_wc$Uniq; binary_wide_depth_wc$Uniq<-NULL
head(binary_wide_depth_wc)

# Let's make "Intersect" column to tell us which months a given ASV shows up in 
binary_wide_depth_wc$Intersect <- apply(binary_wide_depth_wc > 0, 1, function(x){toString(names(binary_wide_depth_wc)[x])})
head(binary_wide_depth_wc)
# Get rid of spaces after common in "Intersect" column
binary_wide_depth_wc$Intersect <- stri_replace_all_fixed(binary_wide_depth_wc$Intersect, " ", "")
head(binary_wide_depth_wc)
# Convert "Intersect" column to list format
binary_wide_depth_wc$Intersect <- strsplit(binary_wide_depth_wc$Intersect, ",")
head(binary_wide_depth_wc)

binary_tax_depth_wc <- binary_wide_depth_wc %>%
  rownames_to_column(var = "uniq") %>% 
  separate(uniq, c("ASV.ID", "Taxonomy"), sep = "_", remove = FALSE) %>% 
  column_to_rownames(var = "uniq") %>% 
  data.frame
head(binary_tax_depth_wc[1:2,])


totals_depth_wc<- binary_tax_depth_wc %>%
  group_by(Intersect) %>%
  summarise(total=sum(length(ASV.ID))) 
totals.num_depth_wc<-order(totals_depth_wc$total,decreasing=TRUE)
total.data_material<-data.frame(as.list(totals_stn$Intersect),totals_stn$total)

binary_tax_depth_wc%>% ggplot(aes(x=Intersect, fill = Taxonomy)) +
  geom_bar(stat = "count", position="stack", color="black") + 
  # geom_text(aes(y=sum(stat(count)), group=Intersect), vjust=0)+
  geom_text(data=totals_depth_wc ,aes(x=Intersect,y=total,label=total,fill=NULL),
            nudge_y = 10)+
  scale_x_upset()+scale_fill_manual(name="Taxonomic Group",values=c(tax_color_phyla2))+
  theme_classic(base_size = 12)+xlab("")+ylab("Number of Shared ASVs")+
  theme(plot.margin = margin(10, 10, 10, 100))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/UpSetR_WC_station_0323.pdf", width=12, height = 10)

#material and depth and water column
# Change to binary
summed_depth_mat_wc<-aggregate(toUpset_wc$value, by=list(Material=toUpset_wc$Material,Depth=toUpset_wc$Depth, Uniq=toUpset_wc$Uniq),sum)
summed_depth_mat_wc$bin<-ifelse(summed_depth_mat_wc$x > 0, 1, 0)
binary_wide_depth_mat_wc<-dcast(summed_depth_mat_wc[c(1,2,3,5)], Uniq~Material+Depth, fill=0)
row.names(binary_wide_depth_mat_wc)<-binary_wide_depth_mat_wc$Uniq; binary_wide_depth_mat_wc$Uniq<-NULL
head(binary_wide_depth_mat_wc)

# Let's make "Intersect" column to tell us which months a given ASV shows up in 
binary_wide_depth_mat_wc$Intersect <- apply(binary_wide_depth_mat_wc > 0, 1, function(x){toString(names(binary_wide_depth_mat_wc)[x])})
head(binary_wide_depth_mat_wc)
# Get rid of spaces after common in "Intersect" column
binary_wide_depth_mat_wc$Intersect <- stri_replace_all_fixed(binary_wide_depth_mat_wc$Intersect, " ", "")
head(binary_wide_depth_mat_wc)
# Convert "Intersect" column to list format
binary_wide_depth_mat_wc$Intersect <- strsplit(binary_wide_depth_mat_wc$Intersect, ",")
head(binary_wide_depth_mat_wc)

binary_tax_depth_mat_wc <- binary_wide_depth_mat_wc %>%
  rownames_to_column(var = "uniq") %>% 
  separate(uniq, c("ASV.ID", "Taxonomy"), sep = "_", remove = FALSE) %>% 
  column_to_rownames(var = "uniq") %>% 
  data.frame
head(binary_tax_depth_mat_wc[1:2,])


totals_depth_mat_wc<- binary_tax_depth_mat_wc %>%
  group_by(Intersect) %>%
  summarise(total=sum(length(ASV.ID))) 
totals.num_depth_mat_wc<-order(totals_depth_mat_wc$total,decreasing=TRUE)
total.data_material<-data.frame(as.list(totals_stn$Intersect),totals_stn$total)

binary_tax_depth_mat_wc%>% ggplot(aes(x=Intersect, fill = Taxonomy)) +
  geom_bar(stat = "count", position="stack", color="black") + 
  # geom_text(aes(y=sum(stat(count)), group=Intersect), vjust=0)+
  geom_text(data=totals_depth_mat_wc ,aes(x=Intersect,y=total,label=total,fill=NULL),
            nudge_y = 10)+
  scale_x_upset()+scale_fill_manual(name="Taxonomic Group",values=c(tax_color_phyla2))+
  theme_classic(base_size = 12)+xlab("")+ylab("Number of Shared ASVs")+
  theme(plot.margin = margin(10, 10, 10, 100))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/UpSetR_WC_station_0323.pdf", width=12, height = 10)

######## PICK UP HERE: then centers of WC#########
asv_df <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Output/MS_Clean_No1_Means_Complete_Norm_NewTax_1122.csv", header=T)
asv_df[1]=NULL
row.names(asv_df) = asv_df$ASV1.ID
# upsetR
melt_ASV<-melt(asv_df)
names(melt_ASV)[14]="Sample"
melt_all<-join(melt_ASV,key_all,by="Sample", type="left", match="first")
head(melt_all)
toUpset<-melt_all[c(1:2,14:18)] # We want the ASV ID, Sample, TaxaPlot, value, Station, Depth, Material
head(toUpset) # Use below to add annotation
toUpset$Uniq<-paste(toUpset$ASV1.ID, toUpset$TaxaPlot, sep="_")

WC = c("DNA", "RNA")
toUpset_wc <- filter(toUpset, Material %in% WC)

toUpset_wc_centers<-filter(toUpset_wc, (Station=="S6" | Station=="S12"))

tax_order<-c("Alveolates-Ciliates","Alveolates-Dinophyceae","Alveolates-Syndiniales", "Archaeplastids-Chlorophytes",
             "Excavates-Discoba", "Hacrobia-Cryptophytes","Hacrobia-Haptophytes",
             "Opisthokont-Choanoflagellida","Opisthokonts-Metazoa","Other/unknown",
             "Rhizaria-Acantharia","Rhizaria-Cercozoa","Rhizaria-Polycystines", "Rhizaria-Other",
             "Stramenopiles-Chrysophytes","Stramenopiles-Diatoms","Stramenopiles-MAST","Stramenopiles-Ochrophyta","Stramenopiles-Pelagophytes","Stramenopiles-Other")
tax_color_phyla2<-c('firebrick4','indianred1','tomato3','forestgreen','yellowgreen','darkblue',
                    'lightblue','moccasin','gold1','grey',"#C291A4",'mediumvioletred','magenta','lightpink',"#DECDBE",'#DDAD4B','tan2','tan3','tan4',"#5C4033")
names(tax_color_phyla2)<-tax_order
binary_tax_material$TAX_ORDER<-factor(binary_tax_material$Taxonomy, levels = (tax_order), labels = (tax_order))

# Change to binary
summed_material_wc_centers<-aggregate(toUpset_wc_centers$value, by=list(Material=toUpset_wc_centers$Material, Uniq=toUpset_wc_centers$Uniq),sum)
summed_material_wc_centers$bin<-ifelse(summed_material_wc_centers$x > 0, 1, 0)
binary_wide_material_wc_centers<-dcast(summed_material_wc_centers[c(1,2,4)], Uniq~Material, fill=0)
row.names(binary_wide_material_wc_centers)<-binary_wide_material_wc_centers$Uniq; binary_wide_material_wc_centers$Uniq<-NULL
head(binary_wide_material_wc_centers)

# Let's make "Intersect" column to tell us which months a given ASV shows up in 
binary_wide_material_wc_centers$Intersect <- apply(binary_wide_material_wc_centers > 0, 1, function(x){toString(names(binary_wide_material_wc_centers)[x])})
head(binary_wide_material_wc_centers)
# Get rid of spaces after common in "Intersect" column
binary_wide_material_wc_centers$Intersect <- stri_replace_all_fixed(binary_wide_material_wc_centers$Intersect, " ", "")
head(binary_wide_material_wc_centers)
# Convert "Intersect" column to list format
binary_wide_material_wc_centers$Intersect <- strsplit(binary_wide_material_wc_centers$Intersect, ",")
head(binary_wide_material_wc_centers)

binary_tax_material_wc_centers <- binary_wide_material_wc_centers %>%
  rownames_to_column(var = "uniq") %>% 
  separate(uniq, c("ASV.ID", "Taxonomy"), sep = "_", remove = FALSE) %>% 
  column_to_rownames(var = "uniq") %>% 
  data.frame
head(binary_tax_material_wc_centers[1:2,])


totals_material_wc_centers<- binary_tax_material_wc_centers %>%
  group_by(Intersect) %>%
  summarise(total=sum(length(ASV.ID))) 
totals.num_material_wc_centers<-order(totals_material_wc_centers$total,decreasing=TRUE)
total.data_material<-data.frame(as.list(totals_stn$Intersect),totals_stn$total)

binary_tax_material_wc_centers%>% ggplot(aes(x=Intersect, fill = Taxonomy)) +
  geom_bar(stat = "count", position="stack", color="black") + 
  # geom_text(aes(y=sum(stat(count)), group=Intersect), vjust=0)+
  geom_text(data=totals_material_wc_centers ,aes(x=Intersect,y=total,label=total,fill=NULL),
            nudge_y = 10)+
  scale_x_upset()+scale_fill_manual(name="Taxonomic Group",values=c(tax_color_phyla2))+
  theme_classic(base_size = 12)+xlab("")+ylab("Number of Shared ASVs")+
  theme(plot.margin = margin(10, 10, 10, 100))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/UpSetR_WC_station_0323.pdf", width=12, height = 10)

######## 03/24/23 noMETA #########
melt_ASV_nometa<-melt(no_meta)
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

########## unique vs shared asv plotting. trying something new   ###########
# filter out the ASVs that are only found in all three samples
shared_asvs <- filter(binary_wide_material_nometa, Intersect=="Water column DNA, Water column RNA, Trap DNA")
shared_asvs$ASV1.ID <- row.names(shared_asvs) # add back the ASVs names
shared_asvs_names <- colsplit(shared_asvs$ASV1.ID, "_", names=c("ASV1.ID_real", "Taxon")) # split the ASV names and Taxon group
# filter the whole dataset by the shared ASVs
shared_asvs_only <- filter(no_meta, ASV1.ID %in% shared_asvs_names$ASV1.ID_real)
sum_shared_each <- colSums(shared_asvs_only[,2:61]) # this gets the number of shared reads in all the samples
sum_total_each <-colSums(no_meta[,2:61]) # this gets the total number of reads in all the samples
percent_shared <- sum_shared_each/sum_total_each
comparison_all <- data.frame(names(percent_shared),percent_shared)
names(comparison_all)[1]<-"Sample"
head(comparison_all)
comparison_plot <- join(comparison_all, key_all, by="Sample")
head(comparison_plot)
col_material<- c("#7fcdbb","#2c7fb8","#fc8d59")
kruskal.test(comparison_plot$percent_shared,comparison_plot$Sample.type)
pairwise.wilcox.test(comparison_plot$percent_shared, comparison_plot$Sample.type, p.adjust="none")

percent_shared<-ggplot(comparison_plot, aes(x=Sample.type, y=percent_shared, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material)+
  # facet_grid(variable~.,scales="free")+
  theme_bw()+
  ylab("Percent of total reads")+
  xlab("")+
  labs(title="Percent of all reads of ASVs shared between all sample types")+
  guides(fill=guide_legend(title="Sample Type"))+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size=8),axis.text.y=element_text(size=12),legend.position = "right")
percent_shared + stat_compare_means(method = "kruskal.test")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/percent_shared_byall_nometa.pdf", width=4, height = 6)

pairwise.wilcox.test(comparison_plot$easy_comparison, comparison_plot$Depth, p.adjust="none")
# sig dif between 500m and all other depths
percent_shared_depth<- ggplot(comparison_plot, aes(x=Depth, y=percent_shared, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material)+
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2))+
  theme_bw()+
  # labs(title="", ylabs="Percent shared reads of total reads")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size=8),axis.text.y=element_text(size=12),legend.position = "right")# +
  facet_grid(Sample.type~.)
percent_shared_depth + stat_compare_means(method = "wilcox.test", paired=TRUE)
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/percent_shared_byall_nometa_depth_together.pdf", width=6, height = 6) 

percent_shared_station <- ggplot(comparison_plot, aes(x=Station, y=percent_shared, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material)+
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2))+
  theme_bw()+
  # labs(title="", ylabs="Percent shared reads of total reads")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size=8),axis.text.y=element_text(size=12),legend.position = "right") # +
  facet_wrap(.~Station)
percent_shared_station + stat_compare_means(method = "wilcox.test", paired=TRUE)
kruskal.test(comparison_plot$percent_shared ~ comparison_plot$Station + comparison_plot$Sample.type)
pairwise.wilcox.test(comparison_plot$percent_shared, comparison_plot$Station, p.adjust="none")


#plot of the relative abundance of each group within the reads 
data.nometa.m<-melt(shared_asvs_only) #melt
head(data.nometa.m)
data.nometa.agg<-aggregate(data.nometa.m$value, by=list(Taxa=data.nometa.m$TaxaPlot),sum) #sum sequences by taxonomic group
data.nometa.agg$Sample<-factor(data.nometa.agg$Sample,levels=names(nometa[2:61]))
shared_nometa_bar_plot_df <- join(data.nometa.agg, key_all, by="Sample", type="left", match="first") %>%
  arrange(., Material, Station.sla, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))

ggplot(data.nometa.agg, aes(y=x, fill=tax, x=""))+
  geom_bar(position = position_fill(),stat="identity",color="black",aes(fill=Taxa))+
  scale_fill_manual(values=tax_color_2, name="Taxa")+
  theme_bw()+
  labs(title="", x="",y="Relative abundance of shared reads")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.text = element_text(size = 8),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 10, angle = 315, hjust = 1, vjust = 0),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5, vjust = .5),
        legend.key.size = unit(0.5,'cm'))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/barplot_shared_byall_nometa.pdf", width=6, height = 6) 

tax_order<-c("Alveolates-Ciliates","Alveolates-Dinophyceae","Alveolates-Syndiniales", "Archaeplastids-Chlorophytes",
             "Excavates-Discoba", "Hacrobia-Cryptophytes","Hacrobia-Haptophytes",
             "Opisthokont-Choanoflagellida","Opisthokonts-Metazoa","Other/unknown",
             "Rhizaria-Acantharia","Rhizaria-Cercozoa","Rhizaria-Polycystines", "Rhizaria-Other",
             "Stramenopiles-Chrysophytes","Stramenopiles-Diatoms","Stramenopiles-MAST","Stramenopiles-Ochrophyta","Stramenopiles-Pelagophytes","Stramenopiles-Other")
tax_color_phyla2<-c('firebrick4','indianred1','tomato3','forestgreen','yellowgreen','darkblue',
                    'lightblue','moccasin','gold1','grey',"#C291A4",'mediumvioletred','magenta','lightpink',"#DECDBE",'#DDAD4B','tan2','tan3','tan4',"#5C4033")
names(tax_color_phyla2)<-tax_order
shared_asvs_only$TAX_ORDER<-factor(shared_asvs_only$Taxonomy, levels = (tax_order), labels = (tax_order))

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
DNA_only <- filter(binary_wide_material_nometa, Intersect=="Water column DNA")
RNA_only <- filter(binary_wide_material_nometa, Intersect=="Water column RNA")
PIT_only <- filter(binary_wide_material_nometa, Intersect=="Trap DNA")
all_unique<- rbind(DNA_only, RNA_only, PIT_only)
all_unique$ASV1.ID <- row.names(all_unique) # add back the ASVs names
unique_names <- colsplit(all_unique$ASV1.ID, "_", names=c("ASV1.ID_real", "Taxon")) # split the ASV names and Taxon group
# filter the whole dataset by the unique ASVs
unique_ASVs <- filter(no_meta, ASV1.ID %in% unique_names$ASV1.ID_real)
sum_unique_each <- colSums(unique_ASVs[,2:61]) # this gets the number of shared reads in all the samples
percent_unique_each <- sum_unique_each/sum_total_each
comparison_unique <- data.frame(names(percent_unique_each),percent_unique_each)
names(comparison_unique)[1]<-"Sample"
head(comparison_unique)
comparison_unique_plot <- join(comparison_unique, key_all, by="Sample")
head(comparison_unique_plot)
col_material<- c("#7fcdbb","#2c7fb8","#fc8d59")
ggplot(comparison_unique_plot, aes(x=Sample.type, y=percent_unique_each, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material)+
  # facet_grid(variable~.,scales="free")+
  theme_bw()+
  # labs(title="Alpha Diversity Metrics-No Metazoa")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size=8),axis.text.y=element_text(size=12),legend.position = "right")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/percent_unique_byall_nometa.pdf", width=4, height = 6)
kruskal.test(comparison_unique_plot$percent_unique_each,comparison_unique_plot$Sample.type)
pairwise.wilcox.test(comparison_unique_plot$percent_unique_each, comparison_unique_plot$Sample.type, p.adjust="none")

ggplot(comparison_unique_plot, aes(x=Depth, y=percent_unique_each, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material)+
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2))+
  theme_bw()+
  # labs(title="", ylabs="Percent shared reads of total reads")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size=8),axis.text.y=element_text(size=12),legend.position = "right")+
  facet_grid(Sample.type~.)
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/percent_unique_byall_nometa_depth.pdf", width=6, height = 6) 
pairwise.wilcox.test(comparison_unique_plot$percent_unique_each, comparison_unique_plot$Depth, p.adjust="none")

percent_unique_station <- ggplot(comparison_unique_plot, aes(x=Station, y=percent_unique_each, fill=Sample.type))+
  geom_boxplot(position = position_dodge(preserve = "single"), width=0.5)+
  scale_fill_manual(values=col_material)+
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2))+
  theme_bw()+
  # labs(title="", ylabs="Percent shared reads of total reads")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size=8),axis.text.y=element_text(size=12),legend.position = "right") # +
facet_wrap(.~Station)
percent_unique_station + stat_compare_means(method = "wilcox.test", paired=TRUE)
kruskal.test(comparison_unique_plot$percent_unique_each ~ comparison_unique_plot$Station)
pairwise.wilcox.test(comparison_unique_plot$percent_unique_each, comparison_unique_plot$Station, p.adjust="none")

sum_all <- sum(trial_2[,2]) # sum the number of shared reads in a sample

pit_1 <- sum(nometa[,2])
trial_percent <- trial_2[,1:2] %>%
  mutate(percent_reads = PIT1_150m_PIT/pit_1)
ggplot(trial_percent, aes(x=ASV1.ID, y=percent_reads))+
  geom_point()
sort_trial <- trial_percent[order(-trial_percent$percent_reads),]
######## ##########
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

tax_order<-c("Alveolates-Ciliates","Alveolates-Dinophyceae","Alveolates-Syndiniales", "Archaeplastids-Chlorophytes",
             "Excavates-Discoba", "Hacrobia-Cryptophytes","Hacrobia-Haptophytes",
             "Opisthokont-Choanoflagellida","Opisthokonts-Metazoa","Other/unknown",
             "Rhizaria-Acantharia","Rhizaria-Cercozoa","Rhizaria-Polycystines", "Rhizaria-Other",
             "Stramenopiles-Chrysophytes","Stramenopiles-Diatoms","Stramenopiles-MAST","Stramenopiles-Ochrophyta","Stramenopiles-Pelagophytes","Stramenopiles-Other")
tax_color_phyla2<-c('firebrick4','indianred1','tomato3','forestgreen','yellowgreen','darkblue',
                    'lightblue','moccasin','gold1','grey',"#C291A4",'mediumvioletred','magenta','lightpink',"#DECDBE",'#DDAD4B','tan2','tan3','tan4',"#5C4033")
names(tax_color_phyla2)<-tax_order
binary_tax_material_nometa$TAX_ORDER<-factor(binary_tax_material_nometa$Taxonomy, levels = (tax_order), labels = (tax_order))

totals_material_nometa<- binary_tax_material_nometa %>%
  group_by(Intersect) %>%
  summarise(total=sum(length(ASV.ID))) %>%
  data.frame
write.csv(as.data.frame(totals_material_nometa), "~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/totals_material.csv")
capture.output(summary(totals_material_nometa), file = "~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/totals_material.txt")
lapply(totals_material_nometa, function(x) write.table( data.frame(x), "~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/totals_material.csv"  , append= T, sep=',' ))

totals.num_material_nometa<-order(totals_material_nometa$total,decreasing=TRUE)
total.data_material_nometa<-data.frame(as.list(totals_material_nometa$Intersect),totals_material_nometa$total)

binary_tax_material_nometa%>% ggplot(aes(x=Intersect, fill = Taxonomy)) +
  geom_bar(position="stack",stat = "count", color="black") + 
  #geom_text(aes(y=sum(stat(count)), group=Intersect), vjust=0)+
  geom_text(data=totals_material_nometa ,aes(x=Intersect,y=total,label=total,fill=NULL),
            nudge_y = 10)+
  scale_x_upset()+scale_fill_manual(name="Taxa",values=c(tax_color_2))+
  theme_classic(base_size = 12)+xlab("")+ylab("Number of Shared ASVs")+
  theme(legend.position="right",
        legend.text = element_text(size = 8),
        plot.margin = margin(10, 10, 10, 100))
# ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/UpsetR/UpSetR_nometa_sample_0323.pdf", width=12, height = 10)

######## balloon plots of specific species #######
## level 3 fungi
fungi_only<-asv_df[asv_df$Level3=="Fungi",]
fungi_sum <- data.frame(apply(fungi_only[2:61], 2, sum)) # sum number of fungi in each sample
fungi_df = cbind(rownames(fungi_sum), fungi_sum) # add the sample column to the dataset
names(fungi_df) = c("Sample", "sum")

fungi_df_m = melt(fungi_df) #reformat to long format for plotting
fungi_df_plotting = join(fungi_df_m,key_all, by="Sample", type ="left", match="first")
fungi_df_plotting$SamplePlot = factor(fungi_df_plotting$Sample, levels=fungi_df_plotting[order(fungi_df_plotting$Material),]$Sample)
fungi_balloon = ggplot(fungi_df_plotting, aes(x=SamplePlot, y=value)) +
  geom_point(aes(size=value, fill=Material), shape=21) +
  scale_size_area(max_size=15) + theme_bw() + theme(axis.text.y = element_text(angle = 45)) + theme(axis.text.x = element_text(angle=90, hjust=1,vjust=1,color="red"))
fungi_balloon
ggplot2::ggsave("MS_Fungi.pdf", width=10, height = 7)

# level 7 diatom
diatom_df = filter(asv_df, TaxaPlot=="Stramenopiles-Diatoms")
diatom_sum <- data.frame(apply(diatom_df[2:61], 2, sum)) # sum number of fungi in each sample
diatom_sum_df = cbind(rownames(diatom_sum), diatom_sum) # add the sample column to the dataset
names(diatom_sum_df) = c("Sample", "sum")
diatom_df_m = melt(diatom_sum_df) #reformat to long format for plotting
diatom_df_plotting = join(diatom_df_m,key_all, by="Sample", type ="left", match="first")
diatom_df_plotting$SamplePlot = factor(diatom_df_plotting$Sample, levels=diatom_df_plotting[order(diatom_df_plotting$Material),]$Sample)
diatom_balloon_material = ggplot(diatom_df_plotting, aes(x=SamplePlot, y=value)) +
  geom_point(aes(size=value, fill=Material), shape=21) +
  scale_size_area(max_size=15) + theme_bw() + theme(axis.text.y = element_text(angle = 45)) + theme(axis.text.x = element_text(angle=90, hjust=1,vjust=1,color="red"))
diatom_balloon_material
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/MS_diatom.pdf", width=10, height = 7)

############# Final Figs 03/09/23 ######################
asv_df <- read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Output/MS_Clean_No1_Means_Complete_Norm_NewTax_1122.csv", header=T)
asv_df[1]=NULL
row.names(asv_df) = asv_df$ASV1.ID
nometa <-filter(asv_df, Taxa!="Opisthokont-Metazoa")

tax_order<-c("Alveolates-Ciliates","Alveolates-Dinophyceae","Alveolates-Syndiniales", "Archaeplastids-Chlorophytes",
             "Excavates-Discoba", "Hacrobia-Cryptophytes","Hacrobia-Haptophytes",
             "Opisthokont-Choanoflagellida","Opisthokonts-Metazoa","Other/unknown",
             "Rhizaria-Acantharia","Rhizaria-Cercozoa","Rhizaria-Polycystines", "Rhizaria-Other",
             "Stramenopiles-Chrysophytes","Stramenopiles-Diatoms","Stramenopiles-MAST","Stramenopiles-Ochrophyta","Stramenopiles-Pelagophytes","Stramenopiles-Other")
tax_color_phyla2<-c('firebrick4','indianred1','tomato3','forestgreen','yellowgreen','darkblue',
                    'lightblue','moccasin','gold1','grey',"#C291A4",'mediumvioletred','#db6ea0','lightpink',"#DECDBE",'#DDAD4B','tan2','tan3','tan4',"#5C4033")
tax_color_phyla2_nometa<-c('firebrick4','indianred1','tomato3','forestgreen','yellowgreen','darkblue',
                           'lightblue','moccasin','grey',"#C291A4",'mediumvioletred','#db6ea0','lightpink',"#DECDBE",'#DDAD4B','tan2','tan3','tan4',"#5C4033")

col_material<- c("#7fcdbb","#2c7fb8","#fc8d59")
WC <- c("DNA","RNA")
col_depth <-c("15m" = "#fee08b", "DCM" = "#7fbc41", "175m" = "#74add1", "500m" = "#542788")

### Creating a new color scheme
color.plot <- ggplot(key_all, aes(x=Station, y=SLA, color=SLA)) +
  geom_point(size=6)+
  scale_color_gradient2(low = "#0818A8",
                        mid = "light gray",
                        high = "#880808",
                        midpoint = 0,
                        limits=c(-28,28))+
  theme_bw()+
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white"))
color.plot
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/SLA_Colors_gray.pdf", width=10, height = 4)

### NMDS All
Norm.rel_nometa <-decostand(nometa[,2:61], MARGIN=2, method = "total") #this does it by rows
colSums(Norm.rel_nometa) # check! should all equal 1.
Norm.rel_nometa$ASV.ID<-rownames(Norm.rel_nometa)
melt_norm_nometa<- melt(Norm.rel_nometa)
vars_nometa<-colsplit(melt_norm_nometa$variable, "_", c("Station","Depth","Material"))
new_df_nometa<-data.frame(melt_norm_nometa, vars_nometa)

all_nometa<- nmdsPoints(new_df_nometa, key_all)
nmds_all_nometa<-ggplot(data=all_nometa, aes(x = MDS1, y = MDS2, fill=SLA,shape=Material)) + 
  geom_point(color="black",size=3,position="jitter") + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = unique(all_nometa$Shape.mat))+ 
  # scale_fill_manual(values=all_nometa$Col.sla)+
  labs(title="", x="NMDS1",y="NMDS2",size=12)+
  scale_y_continuous(breaks = seq(-3, 3, by = 1))+
  scale_x_continuous(breaks = seq(-4, 2, by = 1))+
  scale_fill_gradient2(low = "#0818A8",
                        mid = "light gray",
                        high = "#880808",
                        midpoint = 0,
                        limits=c(-28,28))
  # guides(fill = guide_legend(override.aes=list(values=all_nometa$Col.sla)))
nmds_all <- nmds_all_nometa +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Material),linewidth=1.5) + scale_color_manual(values=col_material) #this adds in 95% confidence interval ellipses from ggplot
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_all_nometa_sla_depth.pdf", width=6, height = 6)

nmds_all_nometa_shape<-ggplot(data=all_nometa, aes(x = MDS1, y = MDS2, fill=SLA,shape=Material)) + 
  geom_point(color="black",size=6,position="jitter") + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = c(21,25,22))+
  scale_fill_gradientn(colours=c("blue", "red"), limits=c(-20,20))+
  labs(title="NMDS-Bray-Curtis- No Metazoa", x="NMDS1",y="NMDS2",size=12)+
  scale_y_continuous(breaks = seq(-3, 3, by = 1))+
  scale_x_continuous(breaks = seq(-4, 2, by = 1))
nmds_all_nometa_shape +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Material),linewidth=1) + scale_color_manual(values=col_material) #this adds in 95% confidence interval ellipses from ggplot
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_all_nometa_sla_material.pdf", width=6, height = 6)

dna_nometa <- filter(new_df_nometa, Material=="DNA")
dna_df_nometa<-nmdsPoints(dna_nometa, key_all)
dna_nmds_nometa<-ggplot(data=dna_df_nometa, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
  geom_point(color="black",size=3,position="jitter") + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = unique(dna_df_nometa$Shape.depth))+ 
  # scale_fill_manual(values=all_nometa$Col.sla)+
  labs(title="Water column 18S rDNA", x="NMDS1",y="NMDS2",size=12)+
  scale_y_continuous(breaks = seq(-3, 3, by = 1))+
  scale_x_continuous(breaks = seq(-4, 2, by = 1))+
  scale_fill_gradient2(low = "#0818A8",
                       mid = "light gray",
                       high = "#880808",
                       midpoint = 0,
                       limits=c(-28,28))
# guides(fill = guide_legend(override.aes=list(values=all_nometa$Col.sla)))
dna_nmds <- dna_nmds_nometa +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth),linewidth=1, linetype="dashed") + scale_color_manual(values=col_depth) #this adds in 95% confidence interval ellipses from ggplot
dna_nmds
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_nometa_sla_v2_dna.pdf", width=6, height =6)


rna_nometa <- filter(new_df_nometa, Material=="cDNA")
rna_df_nometa<-nmdsPoints(rna_nometa, key_all)
rna_nmds_nometa<-ggplot(data=rna_df_nometa, aes(x = MDS1, y = MDS2, fill=SLA,shape=Depth)) + 
  geom_point(color="black",size=3,position="jitter") + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = unique(rna_df_nometa$Shape.depth))+ 
  # scale_fill_manual(values=all_nometa$Col.sla)+
  labs(title="", x="NMDS1",y="NMDS2",size=12)+
  scale_y_continuous(breaks = seq(-3, 3, by = 1))+
  scale_x_continuous(breaks = seq(-4, 2, by = 1))+
  scale_fill_gradient2(low = "#0818A8",
                       mid = "light gray",
                       high = "#880808",
                       midpoint = 0,
                       limits=c(-28,28))
# guides(fill = guide_legend(override.aes=list(values=all_nometa$Col.sla)))
rna_nmds <- rna_nmds_nometa +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Depth),linewidth=1, linetype="dashed") + scale_color_manual(values=col_depth) #this adds in 95% confidence interval ellipses from ggplot

(nmds_all / dna_nmds / rna_nmds) + plot_layout(guides="collect") +  plot_annotation(tag_levels = "a")
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/nmds/meso_nmds_nometa_sla_v2.pdf", width=8, height = 10)


#### Taxa barplots for water column
data.nometa.m<-melt(nometa) #melt
head(data.nometa.m)
data.nometa.agg<-aggregate(data.nometa.m$value, by=list(Taxa=data.nometa.m$TaxaPlot,Sample=data.nometa.m$variable),sum) #sum sequences by taxonomic group
data.nometa.agg$Sample<-factor(data.nometa.agg$Sample,levels=names(nometa[2:61]))
nometa_bar_plot_df <- join(data.nometa.agg, key_all, by="Sample", type="left", match="first") %>%
  arrange(., Material, Station.sla, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))
#Bar plot of community composition ordered by Material, Station and Depth
bars_nometa<-ggplot(nometa_bar_plot_df, aes(y=x, fill=tax, x=Station.sla))+
  geom_bar(position = position_fill(reverse=TRUE),stat="identity",color="black",aes(fill=Taxa))+
  scale_fill_manual(values=tax_color_phyla2_nometa, name="Taxa")+
  theme_bw()+
  labs(title="", x="Station",y="Relative abundance of water column reads")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.text = element_text(size = 8),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 10, angle = 315, hjust = 1, vjust = 0),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5, vjust = .5),
        legend.key.size = unit(0.5,'cm'))
bars_2 <-bars_nometa %+% subset(nometa_bar_plot_df, Material %in% WC) + coord_flip() +facet_grid(Depth~Material) 
tag_facet_outside(bars_2)
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/Taxa_All_NOMeta_0323WC_DepthSLA.pdf", width=8, height = 6)

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
                          total_ciliates$total_ciliates, total_dino$total_dinos, total_synd$total_synd) %>%
  mutate(percent_alv=sum_alv.total_alveolates/total_reads,
         percent_chloro=sum_chloro.total_chlorophytes/total_reads,
         percent_pelago=sum_pelago.total_pelago/total_reads,
         percent_rhizaria=sum_rhizaria.total_rhizaria/total_reads,
         percent_stram=sum_stram.total_stramenopiles/total_reads,
         percent_ciliates = total_ciliates.total_ciliates/total_reads,
         percent_dino = total_dino.total_dinos/total_reads,
         percent_synd = total_synd.total_synd/total_reads)

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
            stdev_synd = sd(percent_synd)) 
percent_rnanodepth_stats <- filter(percent_all, (Material=="RNA")) %>%
  summarise(mean_ciliates = mean(percent_ciliates),
            stdev_ciliates = sd(percent_ciliates),
            mean_dino = mean(percent_dino),
            stdev_dino = sd(percent_dino),
            mean_synd = mean(percent_synd),
            stdev_synd = sd(percent_synd))

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
            stdev_synd = sd(percent_synd))

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
            stdev_synd = sd(percent_synd))

percent_pit_stats <- filter(percent_all, (Material=="PIT")) %>%
  summarise(mean_stram = mean(percent_stram),
            stdev_stram = sd(percent_stram),
            mean_rhiz = mean(percent_rhizaria),
            stdev_rhiz = sd(percent_rhizaria))

# Taxa barplots for the pits
tax_order.pit<-c("Alveolates-Ciliates","Alveolates-Dinoflagellates","Other","Rhizaria-Foraminifera","Rhizaria-Radiolaria","Rhizaria-Other", "Stramenopiles-Diatom","Stramenopiles-Dictyocha")
col.pit <-c('firebrick4','indianred1', 'grey','mediumvioletred',"#C291A4",'lightpink','#DDAD4B',"#5C4033")


counts<-read.csv("~/Desktop/Chapter1/MESOSCOPE2017/Data/CellsPit_Final.csv",header=T)
rownames(counts)=counts$Type
# combine ciliates and tintinnid to the same category
counts$Plotting = counts$Type
counts$Plotting[counts$Type == "Tintinnid"]="Alveolates-Ciliates"
counts$Plotting[counts$Type == "Ciliate"]="Alveolates-Ciliates"
counts$Plotting[counts$Type == "Dinoflagellate"]="Alveolates-Dinoflagellates"
# combine acanths and rads to one category
counts$Plotting[counts$Type == "Rhizaria-Acantharia"]="Rhizaria-Radiolaria"
counts$Plotting[counts$Type == "Rhizaria-other"]="Rhizaria-Other"
counts$Plotting[counts$Type == "Diatom"]="Stramenopiles-Diatom"
counts$Plotting[counts$Type == "Dictyocha"]="Stramenopiles-Dictyocha"
# make the color match the other type

counts.m<-melt(counts, id.vars = c("Type","Plotting"))
head(counts.m)
df<-colsplit(counts.m$variable, "_", c("Station","Rep"))
counts.ms<-data.frame(counts.m,df)
head(counts.ms)

#To calculate average values of taxa by PIT
counts.mean<- counts.ms %>% 
  group_by(Type, Station) %>%
  summarise(mean=mean(value),stdev=sd(value)) %>%
  as.data.frame
counts.mean

# bring back my plotting column
counts.mean.df <- join(counts.mean, select(counts, c("Type","Plotting")), by="Type", type="left")

# need to add the columns
counts_plot <- counts.mean.df %>%
  group_by(Plotting, Station) %>%
  summarise(sum = sum(mean)) %>%
  as.data.frame
counts_plot
tax_order.counts<-c("Alveolates-Ciliates","Alveolates-Dinoflagellates","Rhizaria-Foraminifera","Rhizaria-Radiolaria","Rhizaria-Other", "Stramenopiles-Diatom","Stramenopiles-Dictyocha")
col.pit.counts <-c('firebrick4','indianred1','mediumvioletred','lightpink','#db6ea0','#DDAD4B',"#5C4033")
counts_plot$tax.order<-factor(counts_plot$Plotting,levels=(tax_order.counts), labels=(tax_order.counts))
counts_plot$Station <- factor(counts_plot$Station, levels=c("PIT12","PIT11","PIT10","PIT9","PIT8","PIT7","PIT6","PIT5","PIT4","PIT3","PIT2","PIT1"))

#for relative abundance based on mean
data.agg.count<-aggregate(counts_plot$sum, by=list(Taxa=counts_plot$Plotting,Sample=counts_plot$Station),sum) #sum sequences by taxonomic group
data.agg.count$Sample<-factor(data.agg.count$Sample,levels=c("PIT12","PIT11","PIT10","PIT9","PIT8","PIT7","PIT6","PIT5","PIT4","PIT3","PIT2","PIT1"))
data.agg.count$Sample<-factor(data.agg.count$Sample,levels=c("PIT12","PIT11","PIT10","PIT9","PIT8","PIT7","PIT6","PIT5","PIT4","PIT3","PIT2","PIT1"))

data.agg.count$tax.order<-factor(data.agg.count$Taxa,levels=(tax_order.counts), labels=(tax_order.counts))

barplot<-ggplot(data.agg.count, aes(y=x, fill=tax.order, x=Sample))+
  geom_bar(position = position_fill(reverse=TRUE),stat="identity",color="black",aes(fill=tax.order))+
  scale_fill_manual(values=trial_colors,name="Taxa")+
  theme_bw()+
  labs(title="", x="",y="Relative abundance of trap microscopy counts")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 10, angle = 315, hjust = 1, vjust = 0),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5))
barplot_plot <-barplot + coord_flip() + scale_x_discrete(limits=rev)
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/MicroscopyCellMean_PIT_relav.pdf", width=8, height = 6)

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
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Stramenopiles-Diatoms"]="Stramenopiles-Diatom"
data.m.pit$Plotting[data.m.pit$Level4 == "Dictyochophyceae"]="Stramenopiles-Dictyocha"
# combine acanths and rads to one category
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Rhizaria-Polycystines"]="Rhizaria-Radiolaria"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Rhizaria-Acantharia"]="Rhizaria-Radiolaria"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Rhizaria-Cercozoa"] = "Rhizaria-Other"
data.m.pit$Plotting[data.m.pit$TaxaPlot == "Rhizaria-Other"] = "Rhizaria-Other"

tax_order.pit<-c("Alveolates-Ciliates","Alveolates-Dinoflagellates","Other","Rhizaria-Radiolaria","Rhizaria-Other", "Stramenopiles-Diatom","Stramenopiles-Dictyocha")
col.pit <-c('firebrick4','indianred1', 'grey','lightpink','magenta','#DDAD4B',"#5C4033")

data.agg.pit<-aggregate(data.m.pit$value, by=list(Taxa=data.m.pit$Plotting,Sample=data.m.pit$Station),sum) #sum sequences by taxonomic group
data.agg.pit$tax.order<-factor(data.agg.pit$Taxa,levels=(tax_order.pit), labels=(tax_order.pit))


#Bar plot of community composition ordered by Material, Station and Depth
data.agg.pit$Sample<-factor(data.agg.pit$Sample,levels=c("PIT12","PIT11","PIT10","PIT9","PIT8","PIT7","PIT6","PIT5","PIT4","PIT3","PIT2","PIT1"))

asv_barplot<-ggplot(data.agg.pit, aes(y=x, fill=tax.order, x=Sample))+
  geom_bar(position = position_fill(reverse=TRUE),stat="identity",color="black",aes(fill=tax.order))+
  scale_fill_manual(values=trial_colors, name="Taxa")+
  theme_bw()+
  labs(title="", x="",y="Relative abundance of trap DNA reads")+
  theme(axis.ticks = element_blank())+
  theme(legend.position="right",
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 12, angle = 0, hjust = .5, vjust = .5),
        axis.text.x = element_text(size = 10, angle = 315, hjust = 1, vjust = 0),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),
        axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5))
asv_plot <-asv_barplot + coord_flip() + scale_x_discrete(limits=rev) # +  guides(fill = guide_legend(override.aes=list(color=col_full.pit, labels=tax_full.pit)))
# ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/ASVPIT_relav.pdf", width=8, height = 6)


ggarrange(asv_plot, barplot_plot, common.legend = TRUE, legend="right") +  guides(fill = guide_legend(override.aes=list(color=col_full.pit, labels=tax_full.pit)))
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/PIT_relav_combo.pdf", width=10, height = 6)
tax_full.pit<-c("Alveolates-Ciliates","Alveolates-Dinoflagellates","Other","Rhizaria-Foraminifera","Rhizaria-Radiolaria","Rhizaria-Other", "Stramenopiles-Diatom","Stramenopiles-Dictyocha")
col_full.pit <-c('firebrick4','indianred1', 'grey','mediumvioletred','magenta','lightpink','#DDAD4B',"#5C4033")
trial_colors <- c("Alveolates-Ciliates"='firebrick4',"Alveolates-Dinoflagellates"='indianred1',"Other"= 'grey',
                  "Rhizaria-Foraminifera"='mediumvioletred',"Rhizaria-Radiolaria"='lightpink',"Rhizaria-Other"='#C291A4',
                  "Stramenopiles-Diatom"='#DDAD4B',"Stramenopiles-Dictyocha"="#5C4033")
asv_plot + barplot_plot + plot_layout(guides="collect") 
ggplot2::ggsave("~/Desktop/Chapter1/MESOSCOPE2017/Output/Taxa/PIT_relav_combo_newcolor.pdf", width=10, height = 6)

##### Venn diagram
### Venn diagram of ASVs found in different material type
#Using https://www.r-graph-gallery.com/14-venn-diagramm.html


#First I'm going to sum each ASV by material type
ASV_m<-melt(asv_df)
names(ASV_m)[14]="Sample"

ASV_all<-join(ASV_m,key_all,by="Sample", type="left", match="first")

ASV_pit<-subset(ASV_all,Material=="PIT") #get PIT samples separate from WC
dim(ASV_pit) #419316     35
ASV_pit.2<-dcast(ASV_pit[c(1:2,14,15)], ASV1.ID+TaxaPlot~Sample) #restructure data fame by ASV ID
rownames(ASV_pit.2)=ASV_pit.2$ASV1.ID
rowsum<-apply(ASV_pit.2[,3:14],1,sum) #This sums the rows- so sums the amount of each ASV in all PITs
rowsumm<-data.frame(ASV_pit.2[1:2],rowsum) #this creates a dataframe of the ASVs found
pit.no0<-subset(rowsumm,rowsum>0) #I only want ASVs that have more than 0 sequences in all pits
dim(pit.no0) #3502    3, so 3502 different ASVs in the PITs

ASV_DNA<-subset(ASV_all,Material=="DNA")
dim(ASV_DNA) #838632     20
ASV_dna.2<-dcast(ASV_DNA[c(1:2,14,15)], ASV1.ID+TaxaPlot~Sample) #restructure data fame by ASV ID
rownames(ASV_dna.2)=ASV_dna.2$ASV1.ID
rowsum.dna<-apply(ASV_dna.2[3:26],1,sum) #This sums the rows- so sums the amount of each ASV in all PITs
rowsumm.dna<-data.frame(ASV_dna.2[1:2],rowsum.dna) #this creates a dataframe of the ASVs found
dna.no0<-subset(rowsumm.dna,rowsum.dna>0) #I only want ASVs that have more than 0 sequences in all pits
dim(dna.no0) #21993    3, so 21993 different ASVs in the WC DNA


ASV_RNA<-subset(ASV_all,Material=="RNA")
dim(ASV_RNA) #838632     20
ASV_rna.2<-dcast(ASV_RNA[c(1:2,14,15)], ASV1.ID+TaxaPlot~Sample) #restructure data fame by ASV ID
rownames(ASV_rna.2)=ASV_rna.2$ASV1.ID
rowsum.rna<-apply(ASV_rna.2[3:26],1,sum) #This sums the rows- so sums the amount of each ASV in all PITs
rowsumm.rna<-data.frame(ASV_rna.2[1:2],rowsum.rna) #this creates a dataframe of the ASVs found
rna.no0<-subset(rowsumm.rna,rowsum.rna>0) #I only want ASVs that have more than 0 sequences in all pits
dim(rna.no0) #17739    3, so 17739 different ASVs in the WC RNA

trial<-
##Now trying for the ven diragram 
venn.diagram(
  x=list(dna.no0[,1],rna.no0[,1],pit.no0[,1]),
  category.names=c("DNA","RNA","PITs"),
  
  filename="trial2_venn_diagram.png",
  output=TRUE,
    fill=col_material
)

venn.details<-get.venn.partitions(
  x=list(rna.no0[,1],dna.no0[,1],pit.no0[,1])
)
#column 4 show how they are related
#column 5 shows the ASV IDs included in that set
#column 6 shows the # of counts for that set
head(venn.details[,4])
All<-venn.details[1,5]
All.df<-as.data.frame(All)

#trying with no metazoa
nometa <-filter(asv_df, Taxa!="Opisthokont-Metazoa")
ASV_m_nometa<-melt(nometa)
names(ASV_m_nometa)[14]="Sample"

ASV_all_nometa<-join(ASV_m_nometa,key_all,by="Sample", type="left", match="first")

ASV_pit_nometa<-subset(ASV_all_nometa,Material=="PIT") #get PIT samples separate from WC
dim(ASV_pit_nometa) #400404     35
ASV_pit.2_nometa<-dcast(ASV_pit_nometa[c(1:2,14,15)], ASV1.ID+TaxaPlot~Sample) #restructure data fame by ASV ID
rownames(ASV_pit.2_nometa)=ASV_pit.2_nometa$ASV1.ID
rowsum_nometa<-apply(ASV_pit.2_nometa[,3:14],1,sum) #This sums the rows- so sums the amount of each ASV in all PITs
rowsumm_nometa<-data.frame(ASV_pit.2_nometa[1:2],rowsum_nometa) #this creates a dataframe of the ASVs found
pit.no0_nometa<-subset(rowsumm_nometa,rowsum_nometa>0) #I only want ASVs that have more than 0 sequences in all pits
dim(pit.no0_nometa) #2874    3, so 2874 different ASVs in the PITs

ASV_DNA_nometa<-subset(ASV_all_nometa,Material=="DNA")
dim(ASV_DNA_nometa) #800808     20
ASV_dna.2_nometa<-dcast(ASV_DNA_nometa[c(1:2,14,15)], ASV1.ID+TaxaPlot~Sample) #restructure data fame by ASV ID
rownames(ASV_dna.2_nometa)=ASV_dna.2_nometa$ASV1.ID
rowsum.dna_nometa<-apply(ASV_dna.2_nometa[3:26],1,sum) #This sums the rows- so sums the amount of each ASV in all PITs
rowsumm.dna_nometa<-data.frame(ASV_dna.2_nometa[1:2],rowsum.dna_nometa) #this creates a dataframe of the ASVs found
dna.no0_nometa<-subset(rowsumm.dna_nometa,rowsum.dna_nometa>0) #I only want ASVs that have more than 0 sequences in all pits
dim(dna.no0_nometa) #21055    3, so 21993 different ASVs in the WC DNA


ASV_RNA_nometa<-subset(ASV_all_nometa,Material=="RNA")
dim(ASV_RNA_nometa) #800808     20
ASV_rna.2_nometa<-dcast(ASV_RNA_nometa[c(1:2,14,15)], ASV1.ID+TaxaPlot~Sample) #restructure data fame by ASV ID
rownames(ASV_rna.2_nometa)=ASV_rna.2_nometa$ASV1.ID
rowsum.rna_nometa<-apply(ASV_rna.2_nometa[3:26],1,sum) #This sums the rows- so sums the amount of each ASV in all PITs
rowsumm.rna_nometa<-data.frame(ASV_rna.2_nometa[1:2],rowsum.rna_nometa) #this creates a dataframe of the ASVs found
rna.no0_nometa<-subset(rowsumm.rna_nometa,rowsum.rna_nometa>0) #I only want ASVs that have more than 0 sequences in all pits
dim(rna.no0_nometa) #17397    3, so 17739 different ASVs in the WC RNA


  ##Now trying for the ven diragram 
  venn.diagram(
    x=list(dna.no0_nometa[,1],rna.no0_nometa[,1],pit.no0_nometa[,1]),
    category.names=c("Water column DNA","Water column RNA","Trap DNA"),
    col=c("#7fcdbb","#2c7fb8","#fc8d59"),
    fill = c(alpha("#7fcdbb",0.3), alpha("#2c7fb8",0.3), alpha("#fc8d59",0.3)),
    #format text inside circles
    cex = 1.5,
    fontfamily = "sans",
    #format label text
    cat.cex = 1.2,
    cat.fontfamily = "sans",
    cat.fontface = "bold", 
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.055),
    # file output name
    filename="trial2_venn_diagram_nometa_names.png",
    output=TRUE,
    #fill=col_material
)
venn.details<-get.venn.partitions(
  x=list(rna.no0[,1],dna.no0[,1],pit.no0[,1]))

#column 4 show how they are related
#column 5 shows the ASV IDs included in that set
#column 6 shows the # of counts for that set
head(venn.details[,4])
All<-venn.details[1,5]
All.df<-as.data.frame(All)


########### for real for real ##########

# legend title ("Sample Type")
# c("Water Column DNA", "Water Column RNA", "Trap DNA")
# depth colors c("15m" = "#fee08b", "DCM" = "#7fbc41", "175m" = "#74add1", "500m" = "#542788")
