MESO-SCOPE ASV Bioinformatic Pipeline

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
# ASV table generated in DADA2 using PR2
ASV_table<-read.csv("~/github/mesoscope2017/raw-data/allmeso_table.csv",header=TRUE)
rownames(ASV_table)= ASV_table$ASV.ID # give the rownames ASVs
ASV_taxonomy <- read.csv("~/github/mesoscope2017/raw-data/taxonomy_asv.csv", header=TRUE)
# merge the two tables
ASV_df <- merge(x=ASV_table,y=ASV_taxonomy,by="ASV.ID",all.x=TRUE)

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
Table_Cor<-ASV_df
Table_Cor[,112]=NULL # S6_15m_DNA_B, sample that is all unknowns
sum(Table_Cor[101]) #Check if the sum is 0
Table_Cor[,101]=NULL #S14_500_cDNA_A, sample that has zero reads
Table_Cor[,82]=NULL # Blank_C that is mislabeled
dim(Table_Cor) #36381   111

```

Blank treatment to remove possible contaminants

-   Use the decontam package

```{r}
#Need to separate datasets into "batches" that match with the specific blanks because the PITs were processed at a different time than the water column
# Water RNA and DNA samples 
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
contam_wc<-which(contam_wc_vector == TRUE) #Shows the location of which are true contaminants

### Now for the PIT samples
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
contam_pit<-which(contam_pits_vector == TRUE) #Shows the location of which are true contaminants: 23445 25910 27070 32385 

#Remove the contaminants...From Gerid
#First, we need indexes instead of ASV.ID as rownames
Clean_table<-unrowname(Table_Cor)
#Next, combine the pit and wc vector
contaminates<-c(contam_pit,contam_wc)

Filtered_table<-Clean_table[!rownames(Clean_table) %in% contaminates, ]# Filter the table based on the contaminate vector

#Verification steps. the number of rows in the Filtered table should equal nrows ASV_table - length of contaminants.
nrow(Table_Cor) - nrow(Filtered_table) # 10, so it worked!
#To add back the ASV IDs into the header as row names
rownames(Filtered_table) = Filtered_table$ASV.ID #Here I want to give the rownames the ASV IDs
head(Filtered_table) 
dim(Filtered_table) #  36371 rows and  111 columns, 108 columns of samples, first column is ASV names and last column is taxon list
#I'm going to remove the Blank samples because they have served their purpose
names(Filtered_table[c(2,42,94)])
Filtered_table[,c(2, 42, 94)]=NULL
dim(Filtered_table) # 36371 rows and 108 columns 
#write.csv(Filtered_table,"~/github/mesoscope2017/processed-data/MS_Clean_0723.csv")

```

Calculate mean and standard error on raw reads

```{r}
clean_counts <- read.csv("~/github/mesoscope2017/processed-data/MS_Clean_0723.csv", header=T) 
clean_counts[1]= NULL # remove index column
rownames(clean_counts) = clean_counts$ASV.ID # make rownames the ASV IDs
# first need to take the average of duplicates
clean_mean <- clean_counts %>%
  select(., -108) %>%  # remove confidence of ASV so as not to confuse
  melt %>%
  separate(variable, c("Station","Depth","Type","Rep")) %>% # separate the sample name into components
  dcast(., ASV.ID+Station+Depth+Type~Rep, value.var="value")  %>%  #recast data frame by reps
  mutate(mean = ((A+B)/2))  %>%  #Take the mean of rep A and B
  unite(.,"Sample",c("Station","Depth","Type")) %>%  #regroup the sample names
  dcast(.,ASV.ID~Sample, value.var="mean") #recast the data frame with the mean as the value
dim(clean_mean) # 36371 rows and 61 columns
# write.csv(sample_mean, "~/github/mesoscope2017/processed-data/MS_Clean_mean_0723.csv")

#downloaded it and then combined it with the numbers for the samples without blanks
all_counts<-read.csv("~/github/mesoscope2017/processed-data/MS_Clean_mean_ALL_0723.csv", header=T)
rownames(all_counts)=all_counts$ASV.ID # make rownames ASV IDs

# Reports total number of unique ASVs
length(unique(all_counts$ASV.ID)) #36371

# ASV counts and information:
counts_only<-all_counts[2:61] # remove the taxonomy information
seq_total<-apply(counts_only,2,sum) # sum of all reads
min(seq_total); max(seq_total); mean(seq_total)
ASV_count<-colSums(counts_only>0); ASV_count # sum of ASVs 
min(ASV_count); max(ASV_count); mean(ASV_count)
ASV_single<-colSums(counts_only==1) # sum of ASVs that are only found once by sample 
ASV_double<-colSums(counts_only==2) # sum of ASVs that are only observed twice by sample
ASV_true<-colSums(counts_only>2) # sum of ASVs that are observed more than twice by sample
#
#Compile sample information
sample_info<-data.frame(seq_total,ASV_count,ASV_single,ASV_double,ASV_true);sample_info
#Visual representation of ASV stats:
sample_info$Sample<-row.names(sample_info)
head(sample_info)
counts.melt<-melt(sample_info[c(6,1:4)])
# Need to add in the information about material type of sample
counts.plot <- join(counts.melt, key_all, by="Sample", type="left", match="first") %>%
  arrange(., Material, Station, Depth) %>%  ## Decide here how to order them
  mutate(Sample = factor(Sample, levels=unique(Sample)))

means_type <- tapply(counts.plot$value, list(counts.plot$Material, counts.plot$variable), mean)
min_type <- tapply(counts.plot$value, list(counts.plot$Material, counts.plot$variable), min) # min per type of sample
max_type <- tapply(counts.plot$value, list(counts.plot$Material, counts.plot$variable), max) # max per type of sample
sterror_type <- tapply(counts.plot$value, list(counts.plot$Material, counts.plot$variable), function(x) sd(x) / sqrt(length(x))) # standard error per sample type
type_stats <- data.frame(means_type, min_type, max_type) # combine information
# write.csv(means_type, "~/github/mesoscope2017/processed-data/ASV_means_0524.csv")
# write.csv(type_stats, "~/github/mesoscope2017/processed-data/ASV_stats_0524.csv")

```

Remove global singletons

```{r}
##Filter out ASVs with only 1 sequence in the whole dataset (global singletons)
ASV_sum<-apply(Filtered_table[2:106],1,sum, na.rm =TRUE) #remove global singletons by rows- summing total of reads an ASV is found across all samples
ASV_no1 = Filtered_table[ ASV_sum>1, ]  #count.no1 = ASV table without global singletons
removed = dim(Filtered_table)[1] - dim(ASV_no1)[1] #Outputs the number of ASVs (total) lost in this step, 1,428
names(ASV_no1)[1] = "ASV1.ID"
dim(ASV_no1) #Check dimensions- 34943 rows and 108 columns
#write.csv(ASV_no1, "~/github/mesoscope2017/processed-data/MS_Clean_No1_0723.csv")

```

Take mean of duplicate samples

```{r}
No1_mean <- ASV_no1 %>%
  select(., -108) %>%  # remove confidence of ASV so as not to confuse
  melt %>%
  separate(variable, c("Station","Depth","Type","Rep")) %>% # separate the sample name into components
  dcast(., ASV1.ID+Station+Depth+Type~Rep) %>%  #recast data frame by reps
  mutate(mean = (A+B)/2) %>%  #Take the mean of rep A and B
  unite(.,"Sample",c("Station","Depth","Type")) %>%  #regroup the sample names
  dcast(.,ASV1.ID~Sample, value.var="mean") #recast the data frame with the mean as the value
head(No1_mean)
dim(No1_mean) # 34943 rows and 61 columns
#write.csv(No1_mean, "~/github/mesoscope2017/processed-data/MS_Clean_No1_Means_0723.csv")
#I took the means calculated above, and manually added back the rows that didn't have duplicates
```

Rarefy to check sequencing depth

```{r}
# using the rarecurve from vegan package
# Took the means calculated above, and manually added back the rows that didn't have duplicates
ASV_mean<-read.csv("~/github/mesoscope2017/processed-data/MS_Clean_No1_Means_ALL_0723.csv",header=T)
rownames(ASV_mean)=ASV_mean$ASV1.ID # make row names ASV IDs

t_mean <- ASV_mean[2:61] %>%  # selecting only numeric columns
  as.matrix %>%
  t %>%  # transpose
  floor  # convert to whole numbers
# calculate the rarecurve data
rarecurve_data <- rarecurve(t_mean, step =20)
# color scheme for plot
col_material<- c("Water RNA"="#7fcdbb","Water DNA"="#2c7fb8","Trap DNA"="#fc8d59")
# visualize rarefaction curves from Patt Schloss
map_dfr(rarecurve_data, bind_rows) %>%
  bind_cols(Sample = rownames(t_mean),.) %>%
  pivot_longer(-Sample) %>%
  drop_na() %>%
  mutate(n_seqs= as.numeric(str_replace(name, "N", "")))%>%
  select(-name) %>%
  join(.,key_all, by= "Sample") %>%
  ggplot(aes(x=n_seqs, y=value, group=Sample)) + 
  geom_line(aes(color=Sample.type)) +
  scale_color_manual(values=col_material)+
  labs(title="Rarefaction Curve of ASVs", x="Number of sequences",y="ASVs")+
  theme_bw()

```

Normalize using edgeR

```{r}
ListDGE = DGEList(ASV_mean[,2:61]) 
#ListDGE
ListDGE = calcNormFactors(ListDGE, method = "TMM")
#ListDGE
TMMNorm_ASV = cpm(ListDGE)
#head(TMMNorm_ASV)
TMMNorm_ASV = as.data.frame(TMMNorm_ASV)
TMMNorm_ASV$ASV1.ID = row.names(TMMNorm_ASV)

#write.csv(TMMNorm_ASV, "~/github/mesoscope2017/processed-data/MS_Clean_No1_Means_ALL_Norm_0723.csv")

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
table <- read.csv("~/github/mesoscope2017/processed-data/MS_Clean_No1_Means_ALL_Norm_0723.csv", header=T)
tax <- read.csv("~/github/mesoscope2017/raw-data/taxonomy_asv.csv", header=T)
table[1] = NULL
table[62:63] = NULL
names(tax)[1] = "ASV1.ID"
combined_tax = join(table, tax, by="ASV1.ID", type="left")
NewTax = pr2_rename_taxa_w2(combined_tax)
combined_newnamenewtax = data.frame(combined_tax,NewTax)
#write.csv(combined_newnamenewtax, "~/github/mesoscope2017/processed-data/MS_Clean_No1_Means_ALL_Norm_NewTax_0424.csv")

```
