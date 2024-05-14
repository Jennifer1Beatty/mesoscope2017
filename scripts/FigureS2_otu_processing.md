MESO-SCOPE OTU Bioinformatic Pipeline

By: Jennifer Beatty

Most of the code is modified from Dr. Gerid Ollison (<https://github.com/theOlligist/Daily_dynamics-SMP>)

Last modified: 05/13/2024

Load required libraries

```{r}
library(plyr)
library(tidyverse)
library(reshape2)
library(decontam)  # for removing possible contaminated ASVs
library(edgeR)
library(vegan)
```

Read in ASV table, taxonomy and key file, and combine

```{r}
# OTU table generated in DADA2 using PR2
OTU_table<-read.csv("~/github/mesoscope2017/raw-data/OTU_table_97.csv",header=TRUE)
rownames(OTU_table)= OTU_table$OTU.ID # give the rownames OTU ids
OTU_taxonomy <- read.csv("~/github/mesoscope2017/raw-data/taxonomy_OTU_97_0323.csv", header=TRUE)
# merge the two tables
OTU_df <- merge(x=OTU_table,y=OTU_taxonomy,by="OTU.ID",all.x=TRUE)

# key file with important identifying information for each sample
#### Input all the extra variables and factor them #####
key_all<-read.delim("~/github/mesoscope2017/raw-data/variables_all.txt")
key_all$Station<-factor(as.character(key_all$Station),levels=c("S4","S6","S8","S10","S12","S14",
                                                 "PIT12","PIT11","PIT10","PIT9","PIT8","PIT7","PIT6","PIT5","PIT4","PIT3","PIT2","PIT1"))
key_all$Material<-factor(key_all$Material,levels=c("RNA","DNA","PIT"))
key_all$Depth<-factor(as.character(key_all$Depth),levels=c("15m","DCM","PIT","175m","500m"))
key_all$Sample.type <- factor(as.character(key_all$Sample.type), levels=c( "Water RNA", "Water DNA", "Trap DNA"))

```

Remove error columns

-   One column had zero reads, another column was all unknowns and one sample was mislabeled as "blank"

```{r}
Table_Cor<-OTU_df
Table_Cor[,83]=NULL # S6_15m_DNA_B, remove the sample that is all unknowns
Table_Cor[,60]=NULL  # S14_500m_cDNA_A, sample has 10 reads total
Table_Cor[,4]=NULL # Remove mislabeled blank- Blank_C
dim(Table_Cor) # 67973   111
```

Blank treatment to remove possible contaminants

-   Use the decontam package

```{r}
#Need to separate datasets into "batches" that match with the specific blanks because the PITs were processed at a different time than the water column
# Water RNA and DNA samples 
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
contam_wc<-which(contam_wc_vector == TRUE) #Shows the location of which are true contaminants

### Now for the PIT samples
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

#Remove the contaminants...From Gerid
#First, we need indexes instead of ASV.ID as rownames
Clean_table<-unrowname(Table_Cor)
#Next, combine the pit and wc vector
contaminates<-c(contam_pit,contam_wc)

Filtered_table<-Clean_table[!rownames(Clean_table) %in% contaminates, ]# Filter the table based on the contaminate vector

#Verification steps. the number of rows in the Filtered table should equal nrows ASV_table - length of contaminants.
nrow(Table_Cor) - nrow(Filtered_table) # 19, so it worked!
#To add back the ASV IDs into the header as row names
rownames(Filtered_table) = Filtered_table$OTU.ID #Here I want to give the rownames the OTU IDs
head(Filtered_table) 
dim(Filtered_table) #  67954 rows and  111 columns, 108 columns of samples, first column is OTU names and last column is taxon list
#I'm going to remove the Blank samples because they have served their purpose
names(Filtered_table[c(2,3,16)])
Filtered_table[,c(2,3,16)]=NULL
dim(Filtered_table) # 67954 rows and 108 columns 
# write.csv(Filtered_table,"~/github/mesoscope2017/processed-data/OTU/MS_Clean_OTU_0823.csv")

```

Remove global singletons

```{r}
##Filter out OTUs with only 1 sequence in the whole dataset (global singletons)
OTU_sum<-apply(Filtered_table[2:106],1,sum, na.rm =TRUE) #remove global singletons by rows- summing total of reads an OTU is found across all samples
OTU_no1 = Filtered_table[ OTU_sum>1, ]  #count.no1 = OTU table without global singletons
removed = dim(Filtered_table)[1] - dim(OTU_no1)[1] #Outputs the number of OTUss (total) lost in this step, 30,084
names(OTU_no1)[1] = "OTU.ID"
dim(OTU_no1) #Check dimensions- 37870 rows and 108 columns
#write.csv(OTU_no1, "~/github/mesoscope2017/processed-data/OTU/MS_Clean_No1_0823_OTU.csv")

```

Take mean of duplicate samples

```{r}

No1_m <- OTU_no1 %>%
  select(., -108) %>%  # remove confidence of ASV so as not to confuse
  melt %>%
  separate(variable, c("Station","Depth","Type","Rep")) %>%
  dcast(., OTU.ID+Station+Depth+Type~Rep) %>%  #recast data frame by reps
  mutate(mean = (A+B)/2) %>%  #Take the mean of rep A and B
  unite(.,"Sample",c("Station","Depth","Type")) %>%  #regroup the sample names
  dcast(.,OTU.ID~Sample, value.var="mean") #recast the data frame with the mean as the value
head(No1_mean)
dim(No1_mean) # 37870    61
# write.csv(sample_mean, "~/github/mesoscope2017/processed-data/OTU/MS_Clean_No1_Means_0823_OTU.csv")
#I took the means calculated above, and manually added back the rows that didn't have duplicates
OTU_mean<-read.csv("~/github/mesoscope2017/processed-data/OTU/MS_Clean_No1_Means_0823_OTU_All.csv",header=T)
rownames(OTU_mean)=OTU_mean$OTU.ID
dim(OTU_mean) #37870 rows and 61 columns
```

Normalize using edgeR

```{r}
ListDGE = DGEList(OTU_mean[,2:61]) 
#ListDGE
ListDGE = calcNormFactors(ListDGE, method = "TMM")
#ListDGE
TMMNorm_OTU = cpm(ListDGE)
# head(TMMNorm_OTU)
TMMNorm_OTU = as.data.frame(TMMNorm_OTU)
TMMNorm_OTU$OTU.ID = row.names(TMMNorm_OTU)
Joined<-join(TMMNorm_OTU, OTU_mean[c(1)], by="OTU.ID", type="left", match="first") 
dim(Joined) # 37870 rows and 61 columns
# write.csv(Joined, "~/github/mesoscope2017/processed-data/OTU/MS_Clean_No1_Means_Norm_0823.csv")
```

Rename function to customize PR2 output

```{r}
pr2_rename_taxa_w2<-function(df){
  library(reshape2)
  split<-colsplit(df$Taxon, ";", c("Level1","Level2","Level3","Level4","Level5","Level6", "Level7","Level8"))
  split[ is.na(split) ] = "XXX"
  split[ split == "" ] = "XXX"
  split$Taxa<-"Other/unknown"
  split$Taxa[split$Level1 == "Eukaryota"]="Other/unknown"
  split$Taxa[split$Level2=="Alveolata"]="Alveolata-Other"
  split$Taxa[split$Level2=="Opisthokonta"]="Opisthokonta-Other"
  split$Taxa[split$Level2=="Rhizaria"]="Rhizaria-Other"
  split$Taxa[split$Level2=="Stramenopiles"]="Stramenopila-Other"
  split$Taxa[split$Level2=="Hacrobia"]="Hacrobia-Other"
  split$Taxa[split$Level2=="Archaeplastida"]="Archaeplastida-Other"
  split$Taxa[split$Level2=="Amoebozoa"]="Amoebozoa"
  split$Taxa[split$Level2=="Excavata"]="Excavata"
  split$Taxa[split$Level2=="Eukaryota_X"]="Other/unknown"
  split$Taxa[split$Level2=="Apusozoa"]="Other/unknown"
  split$Taxa[split$Level3=="Dinoflagellata"]="Alveolata-Dinoflagellates"
  split$Taxa[split$Level3=="Metazoa"]="Opisthokonta-Metazoa"
  split$Taxa[split$Level3=="Radiolaria"]="Rhizaria-Radiolaria"
  split$Taxa[split$Level3=="Ciliophora"]="Alveolata-Ciliates"
  split$Taxa[split$Level3=="Ochrophyta"]="Stramenopila-Ochrophytes"
  split$Taxa[split$Level3=="Haptophyta"]="Hacrobia-Haptophytes"
  split$Taxa[split$Level3=="Opisthokonta_X"]="Opisthokonta-Other"
  split$Taxa[split$Level3=="Choanoflagellida"]="Opisthokonta-Choanoflagellates"
  split$Taxa[split$Level3=="Chlorophyta"]="Archaeplastida-Chlorophytes"
  split$Taxa[split$Level3=="Cercozoa"]="Rhizaria-Cercozoans"
  split$Taxa[split$Level3=="Cryptophyta"]="Hacrobia-Cryptophytes"
  split$Taxa[split$Level3=="Fungi"]="Opisthokonta-Fungi"
  split$Taxa[split$Level3=="Discoba"]="Excavata-Discobids"
  split$Taxa[split$Level3=="Foraminifera"]="Rhizaria-Foraminiferans"
  split$Taxa[split$Level4=="Dinophyceae"]="Alveolata-Dinoflagellates"
  split$Taxa[split$Level4=="Syndiniales"]="Alveolata-Syndiniales"
  split$Taxa[split$Level4=="Polycystinea"]="Rhizaria-Polycystines"
  split$Taxa[split$Level4=="Acantharea"]="Rhizaria-Acantharians"
  split$Taxa[split$Level4=="Pelagophyceae"]="Stramenopila-Pelagophytes"
  split$Taxa[split$Level4=="RAD-A"]="Rhizaria-RAD (A,B,C)"
  split$Taxa[split$Level4=="RAD-B"]="Rhizaria-RAD (A,B,C)"
  split$Taxa[split$Level4=="RAD-C"]="Rhizaria-RAD (A,B,C)"
  split$Taxa[split$Level4=="Bacillariophyta"]="Stramenopila-Diatoms"
  split$Taxa[split$Level4=="Chrysophyceae"]="Stramenopila-Chrysophytes"
  split$Taxa[split$Level4=="MAST-1"]="Stramenopila-MAST"
  split$Taxa[split$Level4=="MAST-2"]="Stramenopila-MAST"
  split$Taxa[split$Level4=="MAST-3"]="Stramenopila-MAST"
  split$Taxa[split$Level4=="MAST-4"]="Stramenopila-MAST"
  split$Taxa[split$Level4=="MAST-7"]="Stramenopila-MAST"
  split$Taxa[split$Level4=="MAST-8"]="Stramenopila-MAST"
  split$Taxa[split$Level4=="MAST-9"]="Stramenopila-MAST"
  split$Taxa[split$Level4=="MAST-10"]="Stramenopila-MAST"
  split$Taxa[split$Level4=="MAST-11"]="Stramenopila-MAST"
  split$Taxa[split$Level4=="MAST-12"]="Stramenopila-MAST"
  split$Taxa[split$Level4=="MAST-23"]="Stramenopila-MAST"
  split$Taxa[split$Level4=="MAST-25"]="Stramenopila-MAST"
  split$Taxa2<-"XXX"
  four<-c("Alveolata-Ciliates","Archaeplastida-Chlorophytes")
  split$Taxa2<-with(split, ifelse(Taxa %in% four, Level4, Taxa2)) #Take taxa named in "four" and place in the level4 name, this is the next level of tax resolution I would like to show.
  five<-c("Alveolata-Syndiniales","Stramenopila-MAST")
  split$Taxa2<-with(split, ifelse(Taxa %in% five, Level5, Taxa2))
  two<-c("Alveolata-Other", "Other/unknown")
  split$Taxa2<-with(split, ifelse(Taxa %in% two, Level3, Taxa2))
  six<-c("Hacrobia-Haptophytes", "Hacrobia-Crytophytes")
  split$Taxa2<-with(split, ifelse(Taxa %in% six, Level6, Taxa2))
  four<-c("Opisthokonta-Metazoa", "Opisthokonta-Other")
  split$Taxa2<-with(split, ifelse(Taxa %in% four, Level4, Taxa2))
  seven<-c("Alveolata-Dinoflagellates","Archaeplastida-Other","Excavata","Rhizaria-Acantharians","Rhizaria-Cercozoans","Rhizaria-Polycystines","Rhizaria-RAD (A,B,C)","Stramenopila-Diatoms","Stramenopila-Pelagophytes","Stramenopila-Chrysophytes","Stramenopila-Other")
  split$Taxa2<-with(split, ifelse(Taxa %in% seven, Level7, Taxa2))
  split$Taxa2<-gsub("_XXX", "", split$Taxa2)
  split$Taxa2<-gsub("_XX", "", split$Taxa2)
  split$Taxa2<-gsub("_X", "", split$Taxa2)
  #If taxa2 = "XXX" replace with Taxa, and other
  split$Taxa2<-gsub("XXX","Other/unknown",split$Taxa2)
  # For consistent naming
  split$TaxaPlot<-"XXX"
  selected<-c("Alveolata-Ciliates","Alveolata-Dinoflagellates","Alveolata-Syndiniales","Archaeplastida-Chlorophytes",
          "Excavata-Discobids", "Hacrobia-Cryptophytes","Hacrobia-Haptophytes",
          "Opisthokonta-Choanoflagellates","Opisthokonta-Metazoa","Other/unknown",
          "Rhizaria-Acantharians","Rhizaria-Cercozoans","Rhizaria-Polycystines","Rhizaria-Other",
          "Stramenopila-Chrysophytes", "Stramenopila-Diatoms","Stramenopila-MAST","Stramenopila-Ochrophytes",
          "Stramenopila-Pelagophytes","Stramenopila-Other")
  split$TaxaPlot<-with(split, ifelse(Taxa %in% selected, Taxa, TaxaPlot)) #Take taxa named in "okay" and place in the Taxa
  split$TaxaPlot[split$Taxa=="Rhizaria-RAD (A,B,C)"]="Rhizaria-Other"
  split$TaxaPlot<-gsub("_XXX", "", split$TaxaPlot)
  split$TaxaPlot<-gsub("_XX", "", split$TaxaPlot)
  split$TaxaPlot<-gsub("_X", "", split$TaxaPlot)
  #If taxaPlot = "XXX" replace with Taxa, and other
  split$TaxaPlot<-gsub("XXX","Other/unknown",split$TaxaPlot)
  head(split)
  return(split)
}
```

Apply the renaming function to the ASV taxonomy

```{r}
tax <- read.csv("~/github/mesoscope2017/raw-data/taxonomy_OTU_97_0323.csv", header=T)
names(tax)[1] = "OTU.ID"
table <- join(Joined, tax, by="OTU.ID")
NameTax = pr2_rename_taxa_w2(table)
combined_tax = data.frame(Joined, NameTax)
combined_tax[1:5,]
#write.csv(combined_tax, "~/github/mesoscope2017/processed-data/OTU/MS_Clean_No1_Means_Norm_NewTax_0823.csv")
otu_df<-combined_tax
```

Create figure S2- nMDS of OTUs

```{r}
# filter out metazoa
nometa_otu <-filter(otu_df, Taxa!="Opisthokont-Metazoa")

# function that creates an nmds using vegan package 
nmdsPoints_otu <- function(df,key_df){
  new_df<-dcast(df, variable~OTU.ID, fill=0) #restructure dataframe to only include OTU.ID and samples
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

# all without meta

Norm.rel_nometa_otu <-decostand(nometa_otu[,2:61], MARGIN=2, method = "total") # calculate the relative abundance of each OTU by rows
colSums(Norm.rel_nometa_otu) # check! should all equal 1.
Norm.rel_nometa_otu$OTU.ID<-rownames(Norm.rel_nometa_otu)
# restructure dataframe
new_df_nometa_otu<- Norm.rel_nometa_otu %>% 
  melt() %>%
  separate(variable, c("Station","Depth","Material"))

#run restructured dataframe through nmds function
all_nometa_otu<- nmdsPoints_otu(new_df_nometa_otu, key_all)

#plot results
# color scheme for plot
col_material<- c("Water RNA"="#7fcdbb","Water DNA"="#2c7fb8","Trap DNA"="#fc8d59")

nmds_all_nometa<-ggplot(data=all_nometa_otu, aes(x = MDS1, y = MDS2, fill=SLA,shape=Sample.type)) + 
  geom_point(color="black",size=2,position="jitter") + 
  theme_bw() +
  coord_equal() +
  theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +
  scale_shape_manual(values = c(22,25,21))+ 
  labs(title="", x="NMDS1",y="NMDS2",size=10)+
  scale_fill_gradient2(low = "#0818A8",
                       mid = "light gray",
                       high = "#880808",
                       midpoint = 0,
                       limits=c(-28,28))+
  theme(text = element_text(size=8))
nmds_all <- nmds_all_nometa +  stat_ellipse(inherit.aes = F, mapping=aes(x = MDS1, y = MDS2,col=Sample.type),linewidth=0.5) + scale_color_manual(values=col_material) #this adds in 95% confidence interval ellipses from ggplot
nmds_all

# ggplot2::ggsave("~/github/mesoscope2017/figures/Supplemental1.pdf",  width=3.5, height = 3.5, units="in", dpi=300)

```
