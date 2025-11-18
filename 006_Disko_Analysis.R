#Jason Bosch
#Analysis partly adapted from Gabriele Tosadori's scripts

#####SET UP#####

# ###Reproducibility###
# 
# library(groundhog)
# 
# groundhog_day <- "2024-06-01"
# 
# meta.groundhog(groundhog_day)

###Libraries###

#CRAN libraries
#CHECK WHAT IS NEEDED AND WHAT IS EXCESS
# groundhog.library(c("ggplot2", #Plotting
#                     "scales", #Provide transparency for backgrounds in environmental variables plot
#                     "plyr",
#                     "cowplot", #Combine multiple plots
#                     "dplyr",
#                     "reshape2",
#                     "RColorBrewer",
#                     "ellipse",
#                     "ggpubr",
#                     "ggsignif",
#                     "ggrepel",
#                     "ggtext",
#                     "vegan",
#                     "nlme", #Linear models for soil properties
#                     "multcompView", #Grouping stats
#                     #"stringr", #Extract strings when preparing clustering
#                     #"cluster", #Clustering samples
#                     #"NbClust", #Clustering samples
#                     #"ggdendro", #Colouring the clusters
#                     #"shipunov", #For bootstrapping
#                     "DescTools", #Rounding up for graphs
#                     "lubridate", #For dates in temp measurements
#                     "MicrobiomeStat", #Linda differential abundance
#                     "ggh4x"), #Combining facets
#                   groundhog_day)

library("ggplot2") #Plotting
library("scales") #Provide transparency for backgrounds in environmental variables plot
library("plyr")
library("cowplot") #Combine multiple plots
library("dplyr")
library("reshape2")
library("RColorBrewer")
library("ellipse")
library("ggpubr")
library("ggsignif")
library("ggrepel")
library("ggtext")
library("vegan")
library("nlme")#Linear models for soil properties
library("multcompView") #Grouping stats
library("DescTools") #Rounding up for graphs
library("lubridate") #For dates in temp measurements
library("MicrobiomeStat") #Linda differential abundance 
library("ggh4x") #Combining facets
library("stringr") #Need str_trim()
library("ggnewscale") #Add two scales for the DA plot

#BioConductor libraries
library("phyloseq") #For phyloseq objects

#Github libraries
library("mirlyn") #For rarefaction (Technically Github but installed through BioConductor)
library("ggpubfigs") #For friendly palette
#library("ggdendroplot") #Help in plotting the dendrograms
library("metagMisc") #Get average OTU abundance over multiple rarefactions

###Directory###

setwd("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/04_Analysis_Results/")

#####FUNCTIONS#####

impute_missing_metadata <- function(physeq_object) {
  
  # make a copy of centered_metadata to perform imputation
  imputed_metadata <- as.data.frame(sample_data(physeq_object))
  
  # for the remaining columns (if any!), fill the NA with
  # the median value computed using the available values
  # that are collected from the biological replicates
  for (COL in colnames(imputed_metadata)) {
    # if yes, compute average/median and substitute
    if (any(is.na(imputed_metadata[, COL]))) {
      # find NAs
      is_na <- which(is.na(imputed_metadata[, COL]))
      # for each NA that is found, compute the average using
      # samples from same Month and Treatment
      for (MISSING in is_na) {
        # find sample Treatment
        treat <- imputed_metadata$Treatment[MISSING]
        # find sample Month
        timep <- imputed_metadata$Month[MISSING]
        # get all the datapoints with that Treatment and that Month
        interesting_data <- which(imputed_metadata$Month==timep & imputed_metadata$Treatment==treat)
        all_data <- imputed_metadata[interesting_data, COL]
        # compute median with the available values
        average <- median(unlist(all_data[!is.na(all_data)]))
        # substitute the NAs into the original column
        imputed_metadata[MISSING, COL] <- average
      }
    }
  }
  
  # Set factors
  imputed_metadata$Treatment <- factor(imputed_metadata$Treatment, levels = c("C","SW"),ordered = F)
  imputed_metadata$Treatment <- relevel(imputed_metadata$Treatment,ref = "C")
  imputed_metadata$Month <- factor(imputed_metadata$Month, levels = c("June","July","August"),ordered = T)
  imputed_metadata$Site <- factor(imputed_metadata$Site, levels = unique(imputed_metadata$Site),ordered = F)
  imputed_metadata$Community <- factor(imputed_metadata$Community, levels = c("Total","Active"),ordered = F)
  
  #Replace the metadata
  sample_data(physeq_object) <- imputed_metadata
  
  return(physeq_object)
}

#To pull p-values out of a linear model
#https://stackoverflow.com/questions/5587676/pull-out-p-values-and-r-squared-from-a-linear-regression
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#####LOAD DATA#####

# get metadata
metadata <- read.csv("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/06_Miscellaneous/Disko2013-control_SW_metadata_for_Jason.csv")

# Clean up the names for better plotting
colnames(metadata)[8:37] <- c("SM", "pH", "NH4_soil", "NO3_soil", "PO4_soil","percent_N", "percent_C", "C:N", "TOC_soil", "TN_soil", "DON", "microbial_C", "microbial_N", "microCN", "d15N", "d13C", "percent_SOM", "LOI", "mcrA", "mxaF", "nifH", "nirS", "nosZ", "Deciduous_shrub", "Evergreen_shrub", "Lichen", "Forb", "Litter", "plant_diversity", "plant_eveness")
#Fix the micro_CN to use the same units as everything else
metadata$microCN <- metadata$microbial_C/metadata$microbial_N
#Fix the LOI to use percentages
metadata$LOI <- metadata$LOI*100
# Reshape the metadata to be used by all
metadata_final <- Reduce(rbind, list(setNames(metadata[,c(1,5:37)],c("ID",colnames(metadata[5:37]))),
                                     setNames(metadata[,c(2,5:37)],c("ID",colnames(metadata[5:37]))),
                                     setNames(metadata[,c(3,5:37)],c("ID",colnames(metadata[5:37]))),
                                     setNames(metadata[,c(4,5:37)],c("ID",colnames(metadata[5:37])))))
# Add in DNA/RNA column
metadata_final$Community <- rep(c(rep("Total",36),rep("Active",36)),2)
# Readable names and set factors
#metadata_final$Treatment <- gsub("C","Control",gsub("SW","Treatment",metadata_final$Treatment))
#metadata_final$Treatment <- gsub("SW","T",metadata_final$Treatment)
metadata_final$Treatment <- factor(metadata_final$Treatment, levels = c("C","SW"),ordered = F)
metadata_final$Treatment <- relevel(metadata_final$Treatment,ref = "C")
metadata_final$TimePoint <- factor(metadata_final$TimePoint, levels = c("June","July","August"),ordered = T)
metadata_final$CollectionSite <- factor(metadata_final$CollectionSite, levels = unique(metadata_final$CollectionSite),ordered = F)
metadata_final$Community <- factor(metadata_final$Community, levels = c("Total","Active"),ordered = F)
metadata_final$Description <- paste(metadata_final$Treatment,metadata_final$CollectionSite,metadata_final$TimePoint,metadata_final$Community,sep = "_")
#New column names
#colnames(metadata_final)[2] <- "Treatment"
colnames(metadata_final)[3] <- "Site"
colnames(metadata_final)[4] <- "Month"
# Remove the unused samples and prepare for phyloseq
metadata_final <- subset(metadata_final, !metadata_final$ID=="missing sample")
#metadata_final <- sample_data(metadata_final)
rownames(metadata_final) <- metadata_final$ID

# Bacterial (16S) data
# import counts
bact_counts <- readRDS("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/02_Cleaned_Data/Bacteria_ASV/counts_table_ASV.Rds")
# rename rownames to match metadata sample names
rownames(bact_counts) <- sapply(strsplit(rownames(bact_counts), "_"), "[", 1)
# transpose dada_counts for phyloseq
bact_counts <- t(bact_counts)
#Filtering
#Remove any ASV with fewer 0.1% of the mean number of reads per sample per library size
bact_counts_final <- otu_table(bact_counts[(rowSums(bact_counts)>round(mean(colSums(bact_counts))*0.001)),],taxa_are_rows = TRUE)
cat("Removed ",nrow(bact_counts) - nrow(bact_counts_final)," taxa.") #Removed  890  taxa.
# taxonomy created through naive bayes classifier
bact_taxonomy <- readRDS("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/02_Cleaned_Data/Bacteria_ASV/taxonomy_table_boot_70_NBC_ASV.Rds")
bact_taxonomy_final <- tax_table(bact_taxonomy)
#Create phyloseq object
phyloseq_bact <- phyloseq(bact_counts_final, bact_taxonomy_final, sample_data(metadata_final))

# Fungal (ITS) data
# import counts
#fung_counts <- readRDS("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/02_Cleaned_Data/Fungi_ASV_OTU/counts_table_ASV.Rds")
fung_counts <- readRDS("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/02_Cleaned_Data/Fungi_ASV_OTU/counts_table_OTU_ITSx.Rds")
# rename rownames to match metadata sample names
rownames(fung_counts) <- sapply(strsplit(rownames(fung_counts), "_"), "[", 1)
# transpose dada_counts for phyloseq
fung_counts <- t(fung_counts)
#Filtering
#Remove any ASV with fewer 0.1% of the mean number of reads per sample per library size
fung_counts_final <- otu_table(fung_counts[(rowSums(fung_counts)>round(mean(colSums(fung_counts))*0.001)),],taxa_are_rows = TRUE)
cat("Removed ",nrow(fung_counts) - nrow(fung_counts_final)," taxa.") #Removed  324  taxa.
# taxonomy created through naive bayes classifier
fung_taxonomy <- readRDS("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/02_Cleaned_Data/Fungi_ASV_OTU/taxonomy_table_boot_70_NBC_OTU.Rds")
#Clean the taxonomy names
for (COL in 1:ncol(fung_taxonomy)) {
  fung_taxonomy[,COL] <- gsub("[a-z]__","",fung_taxonomy[,COL])
}
#Phyloseq format
fung_taxonomy_final <- tax_table(fung_taxonomy)
#Create phyloseq object
phyloseq_fung <- phyloseq(fung_counts_final, fung_taxonomy_final, sample_data(metadata_final))

#Filter sequences with incorrect taxonomy
#Bacteria
table(phyloseq_bact@tax_table[,"Kingdom"])
#2 Archaea, 2176 Bacteria
table(as.vector(phyloseq_bact@tax_table=="Mitochondria"))
table(as.vector(phyloseq_bact@tax_table=="Chloroplast"))
phyloseq_bact <- subset_taxa(phyloseq_bact, tax_table(phyloseq_bact)[,"Family"] != "Mitochondria" & tax_table(phyloseq_bact)[,"Order"] != "Chloroplast")
#Removed 71 Mitochondria and 2 Chloroplast
#Fungi
table(phyloseq_fung@tax_table[,"Kingdom"])
# 25 Alveolata, 640 Fungi, 5 Rhizaria, 18 Viridiplantae
phyloseq_fung <- subset_taxa(phyloseq_fung, tax_table(phyloseq_fung)[,"Kingdom"] == "Fungi")

#Make sample names more useful for downstream work
sample_names(phyloseq_bact) <- phyloseq_bact@sam_data$Description
sample_names(phyloseq_fung) <- phyloseq_fung@sam_data$Description

#Library sizes
#Bacteria
summary(colSums(phyloseq_bact@otu_table))
#Fungi
summary(colSums(phyloseq_fung@otu_table))

# #Import FungalTraits database
# FungalTraits <- read.csv("../06_Miscellaneous/FungalTraits 1.2_ver_16Dec_2020 - V.1.2.csv", header = T)
#Import Camelia's altered version of FungalTraits with better identification of ericoid mycorrhiza.
CameliaTraits <- read.csv("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/06_Miscellaneous/guilds_from_seed_changed_by_camelia_for_Disko2013Paper.csv",header = F, sep = ",")
for (COL in 1:6) {
  CameliaTraits[,COL] <- gsub("[a-z]__","",CameliaTraits[,COL])
}

# Load raw temperature data
raw_temps <- read.csv("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/06_Miscellaneous/Top soil temperature Disko data summary summer 2012-summer 2014_raw_data.csv",header = F, sep = ",")
#Have to completely restructure the table into a useable form!
raw_temps_reshaped <- data.frame()
#Cycle over the "C" and "SW" columns
for (COL in which(raw_temps[6,]%in%c("C","SW"))) {
  Sensor <- rep(raw_temps[1,COL],length(10:nrow(raw_temps)))
  SN <- rep(raw_temps[2,COL],length(10:nrow(raw_temps)))
  Type <- rep(raw_temps[3,COL],length(10:nrow(raw_temps)))
  Logger <- rep(raw_temps[4,COL],length(10:nrow(raw_temps)))
  Block <- rep(raw_temps[5,COL],length(10:nrow(raw_temps)))
  Treatment <- rep(raw_temps[6,COL],length(10:nrow(raw_temps)))
  Depth <- rep(raw_temps[7,COL],length(10:nrow(raw_temps)))
  Plot <- rep(raw_temps[8,COL],length(10:nrow(raw_temps)))
  Temps <- raw_temps[10:nrow(raw_temps),COL]
  current_measure <- cbind(raw_temps[10:nrow(raw_temps),1:7],Sensor,SN,Type,Logger,Block,Treatment,Depth,Plot,Temps)
  raw_temps_reshaped <- rbind(raw_temps_reshaped,current_measure)
}
colnames(raw_temps_reshaped) <- c("Copenhagen_time","Disko_date","Disko_time", "Year", "Month", "Day", "DOY","Sensor","SN","Type","Logger","Block","Treatment","Depth","Plot","Temps")
#Only need data from June 2012 - August 2013 (establishment to end of sampling)
raw_temps_reshaped$Time <- as.POSIXct(paste(raw_temps_reshaped[,"Year"],raw_temps_reshaped[,"Month"],raw_temps_reshaped[,"Day"],raw_temps_reshaped[,"Disko_time"],sep="-"), format="%Y-%m-%d-%H:%M", tz = "UTC")
raw_temps_reshaped <- raw_temps_reshaped[!is.na(raw_temps_reshaped$Time),]
raw_temps_reshaped <- raw_temps_reshaped[raw_temps_reshaped$Time <= date("2013-08-31"),]
#Clean up
raw_temps_reshaped <- raw_temps_reshaped[,c("Time","Sensor","SN","Type","Logger","Block","Treatment","Depth","Plot","Temps")]
raw_temps_reshaped$Temps <- as.numeric(raw_temps_reshaped$Temps)
raw_temps_reshaped <- raw_temps_reshaped[!is.na(raw_temps_reshaped$Temps),]
#raw_temps_reshaped$Treatment <- gsub("SW","T",raw_temps_reshaped$Treatment)
raw_temps_reshaped$Treatment <- factor(raw_temps_reshaped$Treatment, levels = c("C","SW"),ordered = T)
raw_temps_reshaped$Depth <- gsub("T_surf","0 cm",gsub("Ts_10cm","10 cm",gsub("Ts_20cm","20 cm",gsub("Ts_5cm","5 cm",raw_temps_reshaped$Depth))))
raw_temps_reshaped <- raw_temps_reshaped[raw_temps_reshaped$Depth%in%c("0 cm","5 cm"),]
raw_temps_reshaped$Depth <- factor(raw_temps_reshaped$Depth, levels = c("0 cm","5 cm"),ordered = T)

#####ANALYSIS#####

###Constant colours for figures###

site_colours <- friendly_pal("vibrant_seven",7)[2:7]
names(site_colours) <- unique(metadata_final$Site)
treatment_colours <- friendly_pal("ito_seven",2)[c(1,2)]
names(treatment_colours) <- unique(metadata_final$Treatment)
timepoint_colours <- friendly_pal("contrast_three",3)[1:3]
names(timepoint_colours) <- unique(metadata_final$Month)
group_colours <- friendly_pal("zesty_four",4)[1:4]
names(group_colours) <- c("C-Total","C-Active","SW-Total","SW-Active")
environmental_colours <- c("#a6611a","#dfc27d","#80cdc1","#018571")
names(environmental_colours) <- c("Soil","Microbial","Functional","Plant")

###Temperatures###

#Get daily averages
raw_temps_reshaped_daily <- aggregate(.~as.Date(Time)*Sensor*SN*Type*Logger*Block*Treatment*Depth*Plot, raw_temps_reshaped,mean)
#Clean up
raw_temps_reshaped_daily <- subset(raw_temps_reshaped_daily,select = -Time)
colnames(raw_temps_reshaped_daily)[1] <- "Time"
raw_temps_reshaped_daily$Time <- as.POSIXct(raw_temps_reshaped_daily$Time, format="%Y-%m-%d", tz = "UTC")
raw_temps_reshaped_daily$Treatment <- factor(raw_temps_reshaped_daily$Treatment, levels = c("C","SW"),ordered = F)
raw_temps_reshaped_daily$Depth <- factor(raw_temps_reshaped_daily$Depth, levels = c("0 cm","5 cm"),ordered = T)

#Plot
soil_temp_full <- ggplot(raw_temps_reshaped_daily, aes(x = Time)) +
  stat_summary(aes(y = Temps, colour = Treatment), fun = mean, geom="line") +
  scale_x_datetime(breaks = "1 month", date_labels = "%b %Y") +
  annotate('rect', xmin=as.POSIXct(date("2013-06-01")), xmax=as.POSIXct(date("2013-07-01")), ymin=-17, ymax=23, alpha=.2, fill=timepoint_colours["June"]) +
  annotate('rect', xmin=as.POSIXct(date("2013-07-01")), xmax=as.POSIXct(date("2013-08-01")), ymin=-17, ymax=23, alpha=.2, fill=timepoint_colours["July"]) +
  annotate('rect', xmin=as.POSIXct(date("2013-08-01")), xmax=as.POSIXct(date("2013-08-31")), ymin=-17, ymax=23, alpha=.2, fill=timepoint_colours["August"]) +
  annotate('rect', xmin=as.POSIXct(date("2013-06-13")), xmax=as.POSIXct(date("2013-06-14")), ymin=-15, ymax=-13, fill="black") +
  annotate('rect', xmin=as.POSIXct(date("2013-07-01")), xmax=as.POSIXct(date("2013-07-05")), ymin=-15, ymax=-13, fill="black") +
  annotate('rect', xmin=as.POSIXct(date("2013-08-22")), xmax=as.POSIXct(date("2013-08-25")), ymin=-15, ymax=-13, fill="black") +
  scale_color_manual(values=treatment_colours) +
  scale_y_continuous(breaks = c(-15,-10,-5,0,5,10,15,20)) +
  coord_cartesian(ylim = c(-17,23),expand = FALSE) +
  facet_wrap(~Depth) +
  labs(x = NULL, y = "Temperature (°C)") +
  theme_bw() +
  theme(panel.grid.major = element_line(color = "gray80", linewidth = 0.5)) +
  theme(text=element_text(size=16), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = NULL, axis.text = element_text(size=16))

soil_temp_zoom <- ggplot(raw_temps_reshaped_daily, aes(x = Time)) +
  stat_summary(aes(y = Temps, colour = Treatment), fun = mean, geom="line") +
  scale_x_datetime(breaks = "1 week", date_labels = "%d %b") +
  annotate('rect', xmin=as.POSIXct(date("2013-06-01")), xmax=as.POSIXct(date("2013-07-01")), ymin=-17, ymax=23, alpha=.2, fill=timepoint_colours["June"]) +
  annotate('rect', xmin=as.POSIXct(date("2013-07-01")), xmax=as.POSIXct(date("2013-08-01")), ymin=-17, ymax=23, alpha=.2, fill=timepoint_colours["July"]) +
  annotate('rect', xmin=as.POSIXct(date("2013-08-01")), xmax=as.POSIXct(date("2013-08-31")), ymin=-17, ymax=23, alpha=.2, fill=timepoint_colours["August"]) +
  annotate('rect', xmin=as.POSIXct(date("2013-06-13")), xmax=as.POSIXct(date("2013-06-14")), ymin=0, ymax=2, fill="black") +
  annotate('rect', xmin=as.POSIXct(date("2013-07-01")), xmax=as.POSIXct(date("2013-07-05")), ymin=0, ymax=2, fill="black") +
  annotate('rect', xmin=as.POSIXct(date("2013-08-22")), xmax=as.POSIXct(date("2013-08-25")), ymin=0, ymax=2, fill="black") +
  scale_color_manual(values=treatment_colours) +
  scale_y_continuous(breaks = c(-5,0,5,10,15,20)) +
  coord_cartesian(ylim = c(-5,23),xlim = c(as.POSIXct(date("2013-06-01")),as.POSIXct(date("2013-08-31"))),expand = FALSE) +
  facet_wrap(~Depth) +
  labs(x = NULL, y = "Temperature (°C)") +
  theme_bw() +
  theme(panel.grid.major = element_line(color = "gray80", linewidth = 0.5)) +
  theme(text=element_text(size=16), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "bottom", axis.text = element_text(size=16))

#Monthly comparison
soil_temp_compare_monthly <- ggplot(raw_temps_reshaped_daily, aes(x = Treatment, y = Temps)) +
  geom_rect(data=data.frame(Time=as.POSIXct(date("2013-06-15")),Treatment="C",Depth=c("0 cm","5 cm"),Temps=0), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = timepoint_colours["June"]) +
  geom_rect(data=data.frame(Time=as.POSIXct(date("2013-07-15")),Treatment="C",Depth=c("0 cm","5 cm"),Temps=0), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = timepoint_colours["July"]) +
  geom_rect(data=data.frame(Time=as.POSIXct(date("2013-08-15")),Treatment="C",Depth=c("0 cm","5 cm"),Temps=0), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = timepoint_colours["August"]) +
  geom_boxplot(aes( colour = Treatment), fill = NA) +
  geom_point(aes(colour = Treatment),position = position_jitterdodge(jitter.width = 0.3),size = 0.3) +
  scale_color_manual(values=treatment_colours) +
  scale_y_continuous(limits = c(-19,26.5), breaks = c(-14,-7,0,7,14,21)) +
  geom_signif(comparisons = list(c("C", "SW")), colour="black", y_position = 21, map_signif_level = T) +
  facet_wrap(~Depth*year(Time)+month(Time,label = T, abbr = F), ncol = 7) +
  labs(x = NULL, y = "Temperature (°C)") +
  theme_bw() +
  theme(panel.grid.major.y = element_line(color = "gray80", linewidth = 0.5)) +
  theme(text=element_text(size=16), legend.position = "none", axis.text = element_text(size=16))

#Get information for results text.
for (DEPTH in c("0 cm","5 cm")) {
  for (MONTH in c("Jun","Jul","Aug")) {
    pvalue <- wilcox.test(raw_temps_reshaped_daily[year(raw_temps_reshaped_daily$Time)==2013 & month(raw_temps_reshaped_daily$Time,label = TRUE)==MONTH & raw_temps_reshaped_daily$Depth==DEPTH & raw_temps_reshaped_daily$Treatment=="SW","Temps"],
                          raw_temps_reshaped_daily[year(raw_temps_reshaped_daily$Time)==2013 & month(raw_temps_reshaped_daily$Time,label = TRUE)==MONTH & raw_temps_reshaped_daily$Depth==DEPTH & raw_temps_reshaped_daily$Treatment=="C","Temps"])$p.value
    difference <- median(raw_temps_reshaped_daily[year(raw_temps_reshaped_daily$Time)==2013 & month(raw_temps_reshaped_daily$Time,label = TRUE)==MONTH & raw_temps_reshaped_daily$Depth==DEPTH & raw_temps_reshaped_daily$Treatment=="SW","Temps"]) -
      median(raw_temps_reshaped_daily[year(raw_temps_reshaped_daily$Time)==2013 & month(raw_temps_reshaped_daily$Time,label = TRUE)==MONTH & raw_temps_reshaped_daily$Depth==DEPTH & raw_temps_reshaped_daily$Treatment=="C","Temps"])
    cat("Median difference between Treatment and control at",DEPTH,"in",MONTH,"is",signif(difference,digits = 3),"and the p-value is",formatC(pvalue, format = "e", digits = 2),"\n",sep = " ")
  }
}

for (DEPTH in c("0 cm","5 cm")) {
  difference <- median(raw_temps_reshaped_daily[raw_temps_reshaped_daily$Time >= as.POSIXct("2012-11-01") & raw_temps_reshaped_daily$Time < as.POSIXct("2013-06-01") & raw_temps_reshaped_daily$Depth==DEPTH & raw_temps_reshaped_daily$Treatment=="SW","Temps"]) -
    median(raw_temps_reshaped_daily[raw_temps_reshaped_daily$Time >= as.POSIXct("2012-11-01") & raw_temps_reshaped_daily$Time < as.POSIXct("2013-06-01") & raw_temps_reshaped_daily$Depth==DEPTH & raw_temps_reshaped_daily$Treatment=="C","Temps"])
  pvalue <- wilcox.test(raw_temps_reshaped_daily[raw_temps_reshaped_daily$Time >= as.POSIXct("2012-11-01") & raw_temps_reshaped_daily$Time < as.POSIXct("2013-06-01") & raw_temps_reshaped_daily$Depth==DEPTH & raw_temps_reshaped_daily$Treatment=="SW","Temps"],
                        raw_temps_reshaped_daily[raw_temps_reshaped_daily$Time >= as.POSIXct("2012-11-01") & raw_temps_reshaped_daily$Time < as.POSIXct("2013-06-01") & raw_temps_reshaped_daily$Depth==DEPTH & raw_temps_reshaped_daily$Treatment=="C","Temps"])$p.value
  cat("Median difference between Treatment and control at",DEPTH,"in between November and June is",signif(difference,digits = 3),"and the p-value is",formatC(pvalue, format = "e", digits = 2),"\n",sep = " ")
}

###Graph the soil chem data###

#Only first 34 should be unique (2 samples are missing)
soil_chem <- metadata_final[1:34,]
#soil_chem <- soil_chem[!soil_chem$Month=="June",]

#Create a table to give human readable titles and units for the graphs
soil_chem_annot <- as.data.frame(matrix(nrow = 30,ncol = 2, dimnames = list(NULL,c("Variable","Unit"))))
soil_chem_annot[,"Variable"] <- c("SM","pH","NH4_soil","NO3_soil","PO4_soil","percent_N","percent_C","C:N","TOC_soil","TN_soil","DON","microbial_C","microbial_N","microCN","d15N","d13C","percent_SOM","LOI","mcrA","mxaF","nifH","nirS","nosZ","Deciduous_shrub","Evergreen_shrub","Lichen","Forb","Litter","plant_diversity","plant_eveness")
soil_chem_annot[,"Name"] <- c("Soil moisture","pH","Soil NH4","Soil NO3","Soil PO4","Nitrogen (%)","Carbon (%)","Carbon:nitrogen","DOC","Dissolved nitrogen","DON","Mic. carbon","Mic. nitrogen","Mic. carbon:nitrogen","d15N","d13C","Soil organic matter","Loss on ignition","mcrA","mxaF","nifH","nirS","nosZ","Dec. shrub cover","EG shrub cover","Lichen cover","Herb cover","Litter cover","Plant diversity","Plant evenness")
soil_chem_annot[,"Unit"] <- c("g/g","","μg/g","μg/g","μg/g","%","%","","mg/g","mg/g","mg/g","mg/g","mg/g","","‰","‰","%","%","copies/mg","copies/mg","copies/mg","copies/mg","copies/mg","%","%","%","%","%","","")
soil_chem_annot[,"Category"] <- c("Soil","Soil","Soil","Soil","Soil","Soil","Soil","Soil","Soil","Soil","Soil","Microbial","Microbial","Microbial","Soil","Soil","Soil","Soil","Functional","Functional","Functional","Functional","Functional","Plant","Plant","Plant","Plant","Plant","Plant","Plant")

#Combined figure showing the different soil properties
nbchem <- list()
nbchem_stat_linear <- list()
nbchem_stat_anova <- list()
nbchem_normality <- list()
for (VAR in soil_chem_annot$Variable) {
  
  #Site is treated as a separate block with random effect
  #Some variables were only measured a year later and so are not split by month
  if (VAR%in%c("Deciduous_shrub","Evergreen_shrub","Lichen","Forb","Litter","plant_diversity","plant_eveness")) {
    #Check for normality
    current_values <- soil_chem[!is.na(soil_chem[,VAR]) & soil_chem$Month=="July", VAR]
    current_n <- length(current_values)
    #Perform statistical tests
    current_ks <- ks.test(current_values,y = "pnorm")$`p.value`
    current_sw <- shapiro.test(current_values)$`p.value`
    #Create a normal distribution over the same distance
    norm_table <- data.frame(X = seq(min(current_values),max(current_values),length.out=100),
                                Y = dnorm(seq(min(current_values),max(current_values),length.out=100),
                                          mean = mean(current_values), 
                                          sd = sd(current_values)))
    #Adjust to the height of the graph
    scaling_factor <- max(hist(current_values, breaks=seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)/max(norm_table$Y)
    norm_table$Y <- norm_table$Y*scaling_factor
    #norm_table$Y <- norm_table$Y*max(hist(current_values, breaks=seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)/dnorm(mean(current_values),mean = mean(current_values), sd = sd(current_values))
    #Plot histogram of the data
    hist_counts <- round(seq(range(hist(current_values, breaks = seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)[1],
                                 range(hist(current_values, breaks = seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)[2],
                                 length.out = 3))
    graph_normality <- ggplot(soil_chem[!is.na(soil_chem[,VAR]) & soil_chem$Month=="July", c("Treatment", "Month", "Site", VAR)],aes(!! sym(VAR))) +
      geom_histogram(breaks = seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),closed="right",colour="black") +
      #stat_function(fun = dnorm, args = list(mean = mean(current_values), sd = sd(current_values)), colour = "blue", 
      #              aes(y = after_stat(y)*max(hist(current_values, breaks = seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)/dnorm(mean(current_values),mean = mean(current_values), sd = sd(current_values)))) +
      #stat_function(fun = dnorm, args = list(mean = mean(soil_chem[!is.na(soil_chem[,VAR]) & soil_chem$Month=="July", VAR]), sd = sd(soil_chem[!is.na(soil_chem[,VAR]) & soil_chem$Month=="July", VAR])), colour = "blue",
      #             aes(y = after_stat(y)*max(hist(soil_chem[!is.na(soil_chem[,VAR]) & soil_chem$Month=="July", VAR], breaks = seq(min(soil_chem[!is.na(soil_chem[,VAR]) & soil_chem$Month=="July", VAR]),max(soil_chem[!is.na(soil_chem[,VAR]) & soil_chem$Month=="July", VAR]),length.out=nclass.Sturges(soil_chem[!is.na(soil_chem[,VAR]) & soil_chem$Month=="July", VAR])+1),plot=F)$counts)/dnorm(mean(soil_chem[!is.na(soil_chem[,VAR]) & soil_chem$Month=="July", VAR]),mean = mean(soil_chem[!is.na(soil_chem[,VAR]) & soil_chem$Month=="July", VAR]), sd = sd(soil_chem[!is.na(soil_chem[,VAR]) & soil_chem$Month=="July", VAR])))) +
      geom_point(data = norm_table, aes(x = X, y = Y)) +
      annotate("text", label = paste("n = ",current_n,sep = ""), x = Inf, y = Inf,hjust=1.1,vjust=1.1, colour = ifelse(current_n < 30,'red','blue')) +
      annotate("text", label = paste("KS = ",formatC(current_ks, format = "e", digits = 1),sep = ""), x = Inf, y = Inf,hjust=1.1,vjust=2.2, colour = ifelse(current_ks < 0.05,'red','blue')) +
      annotate("text", label = paste("SW = ",formatC(current_sw, format = "e", digits = 1),sep = ""), x = Inf, y = Inf,hjust=1.1,vjust=3.3, colour = ifelse(current_sw < 0.05,'red','blue')) +
      labs(title = paste(soil_chem_annot[soil_chem_annot[,"Variable"]==VAR,"Name"])) +
      scale_y_continuous(breaks = hist_counts) +
      ylab(NULL) +
      xlab(label = paste(soil_chem_annot[soil_chem_annot[,"Variable"]==VAR,"Unit"])) +
      theme_bw() +
      theme(text=element_text(size=16, family = "ArialMT"), plot.title = element_text(hjust = 0.5, size=16), axis.text = element_text(size=16),
            panel.background = element_rect(fill = alpha(environmental_colours[soil_chem_annot[soil_chem_annot$Variable==VAR,"Category"]],alpha = 0.4)))

    #Perform linear modelling
    current_data <- soil_chem[soil_chem$Month=="July", c("Treatment", "Month", "Site", VAR)]
    colnames(current_data) <- c("Treatment", "Month", "Site", "an_env")
    set.seed(123456)
    an_env_model <- lme(an_env~Treatment, random=~1|Site, data=current_data, na.action=na.omit)
    current_summary <- summary(an_env_model)
    current_res <- as.data.frame(current_summary$tTable)
    nbchem_stat_linear[VAR] <- list(current_res)

    #Perform anova
    set.seed(123456)
    current_aov <- aov(get(VAR) ~ Site+Treatment, data=soil_chem[soil_chem$Month=="July",])
    current_aov_summary <- summary(current_aov)
    nbchem_stat_anova[VAR] <- list(current_aov_summary)
    current_aov_tukey <- TukeyHSD(current_aov,ordered = TRUE)
    current_tukey_table <- as.data.frame(current_aov_tukey$Treatment)
    current_res <- current_tukey_table$`p adj`
    names(current_res) <- row.names(current_tukey_table)
    current_letters <- multcompLetters(current_res)
    current_letters_df <- data.frame("comp" = names(current_letters$Letters), "Letter"=as.vector(current_letters$Letters))
    current_letters_df$Treatment <- gsub(":.+","",current_letters_df$comp)
    current_letters_df$Treatment <- factor(current_letters_df$Treatment, levels = c("C","SW"),ordered = T)
    current_letters_df[,VAR] <- max(soil_chem[,VAR],na.rm = T)+(0.05*(max(abs(soil_chem[,VAR]),na.rm = T)-min(abs(soil_chem[,VAR]),na.rm = T)))
    
    #Perform the plotting
    graph <- ggplot(soil_chem[soil_chem$Month=="July",], aes(x = Treatment, y = !! sym(VAR))) +
      geom_boxplot(aes(fill=NULL)) +
      geom_point(aes(colour=Site),position = position_dodge(width = 0.2)) +
      geom_text(data = current_letters_df, aes(label=Letter)) +
      labs(title = paste(soil_chem_annot[soil_chem_annot[,"Variable"]==VAR,"Name"])) +
      ylab(label = paste(soil_chem_annot[soil_chem_annot[,"Variable"]==VAR,"Unit"])) +
      scale_color_manual(values = site_colours) +
      xlab(label = NULL) +
      theme_bw() +
      theme(text=element_text(size=16, family = "ArialMT"), legend.position = "right", plot.title = element_text(hjust = 0.5, size=16),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), panel.spacing=unit(0,"lines"),axis.text = element_text(size=16),
            panel.background = element_rect(fill = alpha(environmental_colours[soil_chem_annot[soil_chem_annot$Variable==VAR,"Category"]],alpha = 0.4)))

  }
  else {
    #Check for normality
    current_values <- soil_chem[!is.na(soil_chem[,VAR]),VAR]
    current_n <- length(current_values)
    #Perform statistical tests
    current_ks <- ks.test(current_values,y = "pnorm")$`p.value`
    current_sw <- shapiro.test(current_values)$`p.value`
    #Create a normal distribution over the same distance
    norm_table <- data.frame(X = seq(min(current_values),max(current_values),length.out=100),
                             Y = dnorm(seq(min(current_values),max(current_values),length.out=100),
                                       mean = mean(current_values), 
                                       sd = sd(current_values)))
    #Adjust to the height of the graph
    scaling_factor <- max(hist(current_values, breaks=seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)/max(norm_table$Y)
    norm_table$Y <- norm_table$Y*scaling_factor
    #Plot histogram of the data
    hist_counts <- round(seq(range(hist(current_values, breaks = seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)[1],
                             range(hist(current_values, breaks = seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)[2],
                             length.out = 3))
    graph_normality <- ggplot(soil_chem[!is.na(soil_chem[,VAR]), c("Treatment", "Month", "Site", VAR)],aes(!! sym(VAR))) +
      geom_histogram(breaks = seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),closed="right",colour="black") +
      geom_point(data = norm_table, aes(x = X, y = Y)) +
      annotate("text", label = paste("n = ",current_n,sep = ""), x = Inf, y = Inf,hjust=1.1,vjust=1.1, colour = ifelse(current_n < 30,'red','blue')) +
      annotate("text", label = paste("KS = ",formatC(current_ks, format = "e", digits = 1),sep = ""), x = Inf, y = Inf,hjust=1.1,vjust=2.2, colour = ifelse(current_ks < 0.05,'red','blue')) +
      annotate("text", label = paste("SW = ",formatC(current_sw, format = "e", digits = 1),sep = ""), x = Inf, y = Inf,hjust=1.1,vjust=3.3, colour = ifelse(current_sw < 0.05,'red','blue')) +
      labs(title = paste(soil_chem_annot[soil_chem_annot[,"Variable"]==VAR,"Name"])) +
      scale_y_continuous(breaks = hist_counts) +
      ylab(NULL) +
      xlab(label = paste(soil_chem_annot[soil_chem_annot[,"Variable"]==VAR,"Unit"])) +
      theme_bw() +
      theme(text=element_text(size=16, family = "ArialMT"), plot.title = element_text(hjust = 0.5, size=16), axis.text = element_text(size=16),
            panel.background = element_rect(fill = alpha(environmental_colours[soil_chem_annot[soil_chem_annot$Variable==VAR,"Category"]],alpha = 0.4)))
 
    #Perform linear modelling
    current_data <- soil_chem[, c("Treatment", "Month", "Site", VAR)]
    colnames(current_data) <- c("Treatment", "Month", "Site", "an_env")
    set.seed(123456)
    an_env_model <- lme(an_env~Treatment*Month, random=~1|Site, data=current_data, na.action=na.omit)
    current_summary <- summary(an_env_model)
    current_res <- as.data.frame(current_summary$tTable)
    nbchem_stat_linear[VAR] <- list(current_res)
    
    #Perform anova
    set.seed(123456)
    current_aov <- aov(get(VAR) ~ Site+Treatment*Month, data=soil_chem)
    current_aov_summary <- summary(current_aov)
    nbchem_stat_anova[VAR] <- list(current_aov_summary)
    current_aov_tukey <- TukeyHSD(current_aov,ordered = TRUE)
    current_tukey_table <- as.data.frame(current_aov_tukey$`Treatment:Month`)
    current_res <- current_tukey_table$`p adj`
    names(current_res) <- row.names(current_tukey_table)
    current_letters <- multcompLetters(current_res)
    current_letters_df <- data.frame("comp" = names(current_letters$Letters), "Letter"=as.vector(current_letters$Letters))
    current_letters_df$Treatment <- gsub(":.+","",current_letters_df$comp)
    current_letters_df$Month <- gsub("^.+:{1}","",current_letters_df$comp)
    current_letters_df$Treatment <- factor(current_letters_df$Treatment, levels = c("C","SW"),ordered = T)
    current_letters_df$Month <- factor(current_letters_df$Month, levels = c("June","July","August"),ordered = T)
    current_letters_df[,VAR] <- max(soil_chem[,VAR],na.rm = T)+(0.05*(max(abs(soil_chem[,VAR]),na.rm = T)-min(abs(soil_chem[,VAR]),na.rm = T)))
    
    #Perform the plotting
    graph <- ggplot(soil_chem, aes(x = Month, y = !! sym(VAR))) +
      geom_boxplot(aes(fill=NULL)) +
      stat_summary(geom = "point", fun = mean, shape = 18, size = 3) +
      geom_point(aes(colour=Site),position = position_dodge(width = 0.2)) +
      geom_text(data = current_letters_df, aes(label=Letter)) +
      labs(title = paste(soil_chem_annot[soil_chem_annot[,"Variable"]==VAR,"Name"])) +
      ylab(label = paste(soil_chem_annot[soil_chem_annot[,"Variable"]==VAR,"Unit"])) +
      facet_wrap(~Treatment) +
      scale_color_manual(values = site_colours) +
      xlab(label = NULL) +
      theme_bw() +
      theme(text=element_text(size=16, family = "ArialMT"), legend.position = "right", plot.title = element_text(hjust = 0.5, size = 16),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), panel.spacing=unit(0,"lines"),axis.text = element_text(size=16),
              panel.background = element_rect(fill = alpha(environmental_colours[soil_chem_annot[soil_chem_annot$Variable==VAR,"Category"]],alpha = 0.2)))
  }
  
  #Add to list
  nbchem[VAR] <- list(graph)
  nbchem_normality[VAR] <- list(graph_normality)
  
}

#Make a table of VARs with significant differences according to linear models
Soil_chem_table <- data.frame()
signif_chems <- NULL
for (VAR in names(nbchem_stat_linear)) {
  
  if (VAR%in%c("Deciduous_shrub","Evergreen_shrub","Lichen","Forb","Litter","plant_diversity","plant_eveness")) {
    current <- nbchem_stat_linear[[VAR]]
    current_df <- as.data.frame(matrix(nrow = 2,ncol = 7, dimnames = list(NULL,c("Variable","Test","Value","Std.Error","DF","t value","p value"))))
    current_df$Variable <- VAR
    current_df$Test <- str_trim(row.names(current))
    current_df$Value <- current$Value
    current_df$Std.Error <- current$Std.Error
    current_df$DF <- current$DF
    current_df$`t value` <- current$`t-value`
    current_df$`p value` <- current$`p-value`
    #current_df$Test <- gsub(".L","",current_df$Test)
    #if (current_df[current_df$Test=="Treatment","p value"] < 0.05) {signif_chems <- c(signif_chems,VAR)}
  }
  
  else {
    current <- nbchem_stat_linear[[VAR]]
    current_df <- as.data.frame(matrix(nrow = 6,ncol = 7, dimnames = list(NULL,c("Variable","Test","Value","Std.Error","DF","t value","p value"))))
    current_df$Variable <- VAR
    current_df$Test <- str_trim(row.names(current))
    current_df$Value <- current$Value
    current_df$Std.Error <- current$Std.Error
    current_df$DF <- current$DF
    current_df$`t value` <- current$`t-value`
    current_df$`p value` <- current$`p-value`
    #current_df$Test <- gsub(".L","",current_df$Test)
    #if (current_df[current_df$Test=="Treatment","p value"] < 0.05 | current_df[current_df$Test=="Treatment:Month","p value"] < 0.05) {signif_chems <- c(signif_chems,VAR)}
  }
  
  Soil_chem_table <- rbind(Soil_chem_table,current_df)
}
# Soil_chem_table_signif <- Soil_chem_table[Soil_chem_table$Variable%in%signif_chems,]
# Soil_chem_table_signif <- Soil_chem_table_signif[!Soil_chem_table_signif$Test=="(Intercept)",]
# Soil_chem_table_signif <- Soil_chem_table_signif[!Soil_chem_table_signif$Test=="Month",]
Soil_chem_table_signif <- Soil_chem_table
Soil_chem_table_signif <- Soil_chem_table_signif[!Soil_chem_table_signif$Test=="(Intercept)",]
Soil_chem_table_signif <- Soil_chem_table_signif[Soil_chem_table_signif$`p value` < 0.05,]
#Soil_chem_table_signif <- Soil_chem_table_signif[grepl("Treatment",Soil_chem_table_signif$Test),]
Soil_chem_table_signif <- Soil_chem_table_signif[!grepl(".Q",Soil_chem_table_signif$Test),]
Soil_chem_table_signif$Test <- gsub(".L","",Soil_chem_table_signif$Test)
Soil_chem_table_signif$Test <- gsub("TreatmentSW","Treatment",Soil_chem_table_signif$Test)
#For interpretation: https://stackoverflow.com/questions/57297771/interpretation-of-l-q-c-4-for-logistic-regression

#Clean up the names
for (VAR in Soil_chem_table$Variable) {
  Soil_chem_table[which(Soil_chem_table$Variable==VAR),"Variable"] <- soil_chem_annot[soil_chem_annot$Variable==VAR,"Name"]
}
for (VAR in Soil_chem_table_signif$Variable) {
  Soil_chem_table_signif[which(Soil_chem_table_signif$Variable==VAR),"Variable"] <- soil_chem_annot[soil_chem_annot$Variable==VAR,"Name"]
}

# #Export tables
# write.table(Soil_chem_table, file = "Soil_chem_table.csv",quote = T,row.names = F, col.names = T, sep = ",")
# write.table(Soil_chem_table_signif, file = "Soil_chem_table_signif.csv",quote = T,row.names = F, col.names = T, sep = ",")

# #Plot the variables which are important for the RDA in order of how often they come up
# soil_chem_plot <- plot_grid(nbchem[[names(sort(table(RDA_arrows_signif$envvars),decreasing = T))[1]]] + theme(legend.position = "none"),
#                             nbchem[[names(sort(table(RDA_arrows_signif$envvars),decreasing = T))[2]]] + theme(legend.position = "none"),
#                             nbchem[[names(sort(table(RDA_arrows_signif$envvars),decreasing = T))[3]]] + theme(legend.position = "none"),
#                             nbchem[[names(sort(table(RDA_arrows_signif$envvars),decreasing = T))[4]]] + theme(legend.position = "none"),
#                             nbchem[[names(sort(table(RDA_arrows_signif$envvars),decreasing = T))[5]]] + theme(legend.position = "none"),
#                             nbchem[[names(sort(table(RDA_arrows_signif$envvars),decreasing = T))[6]]] + theme(legend.position = "none"),
#                             nbchem[[names(sort(table(RDA_arrows_signif$envvars),decreasing = T))[7]]] + theme(legend.position = "none"),
#                             nbchem[[names(sort(table(RDA_arrows_signif$envvars),decreasing = T))[8]]] + theme(legend.position = "none"),
#                             nbchem[[names(sort(table(RDA_arrows_signif$envvars),decreasing = T))[9]]] + theme(legend.position = "none"),
#                             nbchem[[names(sort(table(RDA_arrows_signif$envvars),decreasing = T))[10]]] + theme(legend.position = "none"),
#                             nbchem[[names(sort(table(RDA_arrows_signif$envvars),decreasing = T))[11]]] + theme(legend.position = "none"),
#                             nbchem[[names(sort(table(RDA_arrows_signif$envvars),decreasing = T))[12]]] + theme(legend.position = "none"),
#                             ncol = 4)
# 
# ggsave("RDA_Soil_Chem.svg",soil_chem_plot,width = 3*168, height = 2*168, units = "mm")
# ggsave("RDA_Soil_Chem.jpg",soil_chem_plot,width = 3*168, height = 2*168, units = "mm")

#Save normality plot for internal use
Environmental_Variables_Normality <- plot_grid(nbchem_normality[[1]] + theme(legend.position = "none"), nbchem_normality[[2]] + theme(legend.position = "none"),
                                               nbchem_normality[[3]] + theme(legend.position = "none"), nbchem_normality[[4]] + theme(legend.position = "none"),
                                               nbchem_normality[[5]] + theme(legend.position = "none"), nbchem_normality[[6]] + theme(legend.position = "none"),
                                               nbchem_normality[[7]] + theme(legend.position = "none"), nbchem_normality[[8]] + theme(legend.position = "none"),
                                               nbchem_normality[[9]] + theme(legend.position = "none"), nbchem_normality[[10]] + theme(legend.position = "none"),
                                               nbchem_normality[[11]] + theme(legend.position = "none"), nbchem_normality[[15]] + theme(legend.position = "none"),
                                               nbchem_normality[[16]] + theme(legend.position = "none"), nbchem_normality[[17]] + theme(legend.position = "none"),
                                               nbchem_normality[[18]] + theme(legend.position = "none"), nbchem_normality[[12]] + theme(legend.position = "none"),
                                               nbchem_normality[[13]] + theme(legend.position = "none"), nbchem_normality[[14]] + theme(legend.position = "none"),
                                               nbchem_normality[[19]] + theme(legend.position = "none"), nbchem_normality[[20]] + theme(legend.position = "none"),
                                            nbchem_normality[[21]] + theme(legend.position = "none"), nbchem_normality[[22]] + theme(legend.position = "none"),
                                            nbchem_normality[[23]] + theme(legend.position = "none"), nbchem_normality[[24]] + theme(legend.position = "none"),
                                            nbchem_normality[[25]] + theme(legend.position = "none"), nbchem_normality[[26]] + theme(legend.position = "none"),
                                            nbchem_normality[[27]] + theme(legend.position = "none"), nbchem_normality[[28]] + theme(legend.position = "none"),
                                            nbchem_normality[[29]] + theme(legend.position = "none"), nbchem_normality[[30]] + theme(legend.position = "none"),
                                            ncol = 5, align="vh", axis = "tblr")
#SVG
svg("Environmental_Variables_Normality.svg", width=2*210*0.03937008, height=2*297*0.03937008)
Environmental_Variables_Normality
dev.off()
#JPG
jpeg("Environmental_Variables_Normality.jpg", width=2*210, height=2*297, units = "mm", res = 300)
Environmental_Variables_Normality
dev.off()

###Alpha diversity###

#Rarefy communities
#Using mirl to repeatedly rarefy to specified library size 1000 times.
mirl_object_bact <- mirl(phyloseq_bact, libsize = min(colSums(phyloseq_bact@otu_table)), rep = 1000, mc.cores = 6L, set.seed = 123456)
mirl_object_fung <- mirl(phyloseq_fung, libsize = min(colSums(phyloseq_fung@otu_table)), rep = 1000, mc.cores = 6L, set.seed = 123456)

# Generates dataframe of alpha-diversity metric from mirl_object
# Can't use the default because we want different metrics
alphadiv_bact <- data.frame()
for (REP in 1:length(mirl_object_bact)) {
  cat(paste("Calculating rep ",REP,"\n",sep = ""))
  md <- sample_data(mirl_object_bact[[REP]])
  div <- estimate_richness(mirl_object_bact[[REP]],measures = c("Observed","Simpson","Shannon"))
  final <- cbind(md, div)
  alphadiv_bact <- rbind(alphadiv_bact,final)
}

alphadiv_fung <- data.frame()
for (REP in 1:length(mirl_object_fung)) {
  cat(paste("Calculating rep ",REP,"\n",sep = ""))
  md <- sample_data(mirl_object_fung[[REP]])
  div <- estimate_richness(mirl_object_fung[[REP]],measures = c("Observed","Simpson","Shannon"))
  final <- cbind(md, div)
  alphadiv_fung <- rbind(alphadiv_fung,final)
}

#Delete the mirl objects as they are quite large and unneeded now
rm(mirl_object_bact,mirl_object_fung)

#Visualise the variation in the rarefactions
observed_bact_rares <- ggplot(alphadiv_bact, aes(x=Treatment, y=Observed)) +
  geom_point (aes(colour = Site),position = position_jitterdodge(jitter.width = 0.1),size = 0.3)  +
  geom_boxplot(fill=NA,size=0.5) +
  theme_bw() +
  scale_y_continuous(limits = c(0,RoundTo(max(alphadiv_bact[,"Observed"]),multiple = 50,ceiling)),n.breaks = 6) +
  scale_color_manual(values=site_colours) +
  facet_nested(~Month+Community) +
  labs(title = "Bacterial richness with 1000X rarefactions", x = NULL, y = "Richness") +
  theme(text=element_text(size=16, family = "ArialMT"), legend.position = "right", plot.title = element_text(hjust = 0.5), 
        panel.spacing=unit(0,"lines"))

observed_fung_rares <- ggplot(alphadiv_fung, aes(x=Treatment, y=Observed)) +
  geom_point (aes(colour = Site),position = position_jitterdodge(jitter.width = 0.1),size = 0.3)  +
  geom_boxplot(fill=NA,size=0.5) +
  theme_bw() +
  scale_y_continuous(limits = c(0,RoundTo(max(alphadiv_fung[,"Observed"]),multiple = 20,ceiling)),n.breaks = 6) +
  scale_color_manual(values=site_colours) +
  facet_nested(~Month+Community) +
  labs(title = "Fungal richness with 1000X rarefactions", x = NULL, y = "Richness") +
  theme(text=element_text(size=16, family = "ArialMT"), legend.position = "right", plot.title = element_text(hjust = 0.5), 
        panel.spacing=unit(0,"lines"))

#Mean of the 1000 rarefactions
alphadiv_bact_mean <- aggregate(alphadiv_bact[,c("ID","Observed","Simpson","Shannon")],.~ID,mean)
alphadiv_bact_mean <- merge.data.frame(alphadiv_bact_mean,metadata_final,by="ID")
alphadiv_bact_mean$library_size <-  colSums(phyloseq_bact@otu_table)[alphadiv_bact_mean$Description]
alphadiv_fung_mean <- aggregate(alphadiv_fung[,c("ID","Observed","Simpson","Shannon")],.~ID,mean)
alphadiv_fung_mean <- merge.data.frame(alphadiv_fung_mean,metadata_final,by="ID")
alphadiv_fung_mean$library_size <-  colSums(phyloseq_fung@otu_table)[alphadiv_fung_mean$Description]

#Check the correlation with sequencing depth
plot_bact_depth_corr <- ggplot(alphadiv_bact_mean, aes(x=library_size, y=Observed)) +
  geom_point(aes(colour = Site)) +
  geom_smooth(method = "lm", colour = "black") +
  annotate("text", label = paste("Adj. R^2 = ",formatC(summary(lm(Observed~library_size,alphadiv_bact_mean))$adj.r.squared, format = "f", digits = 3),"\n",
                                 "p = ",formatC(lmp(lm(Observed~library_size,alphadiv_bact_mean)), format = "e", digits = 1),sep = ""),x = Inf, y = -Inf,hjust=1,vjust=-0.3) +
  theme_bw() +
  scale_y_continuous(n.breaks = 6) +
  scale_x_continuous(n.breaks = 6) +
  labs(x = "Library size", y = "Richness") +
  scale_color_manual(values=site_colours) +
  theme(text=element_text(size=16), legend.position = "right", plot.title = element_text(hjust = 0.5, size=16), 
        panel.spacing=unit(0,"lines"),axis.text = element_text(size=16))

plot_fung_depth_corr <- ggplot(alphadiv_fung_mean, aes(x=library_size, y=Observed)) +
  geom_point(aes(colour = Site)) +
  geom_smooth(method = "lm", colour = "black")  +
  annotate("text", label = paste("Adj. R^2 = ",formatC(summary(lm(Observed~library_size,alphadiv_fung_mean))$adj.r.squared, format = "f", digits = 3),"\n",
                                 "p = ",formatC(lmp(lm(Observed~library_size,alphadiv_fung_mean)), format = "f", digits = 3),sep = ""),x = Inf, y = -Inf,hjust=1,vjust=-0.3) +
  theme_bw() +
  scale_y_continuous(n.breaks = 6) +
  scale_x_continuous(n.breaks = 6) +
  labs(x = "Library size", y = "Richness") +
  scale_color_manual(values=site_colours) +
  theme(text=element_text(size=16), legend.position = "right", plot.title = element_text(hjust = 0.5, size=16), 
        panel.spacing=unit(0,"lines"),axis.text = element_text(size=16))

#Plots and stats
alpha_plots <- NULL
alpha_plots_normality <- NULL
alpha_anova <- NULL
alpha_Tukey <- NULL
for (ALPHA in c("alphadiv_bact_mean","alphadiv_fung_mean")) {
  if (grepl("bact",ALPHA)) {taxonomic_group <- "Bacterial" }
  if (grepl("fung",ALPHA)) {taxonomic_group <- "Fungal" }
  ##Richness##
  #Check for normality
  current_values <- get(ALPHA)[,"Observed"]
  current_n <- length(current_values)
  #Perform statistical tests
  current_ks <- ks.test(current_values,y = "pnorm")$`p.value`
  current_sw <- shapiro.test(current_values)$`p.value`
  #Create a normal distribution over the same distance
  norm_table <- data.frame(X = seq(min(current_values),max(current_values),length.out=100),
                           Y = dnorm(seq(min(current_values),max(current_values),length.out=100),
                                     mean = mean(current_values), 
                                     sd = sd(current_values)))
  #Adjust to the height of the graph
  scaling_factor <- max(hist(current_values, breaks=seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)/max(norm_table$Y)
  norm_table$Y <- norm_table$Y*scaling_factor
  #Plot histogram of the data
  hist_counts <- round(seq(range(hist(current_values, breaks = seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)[1],
                           range(hist(current_values, breaks = seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)[2],
                           length.out = 3))
  observed_normality <- ggplot(get(ALPHA),aes(Observed)) +
    geom_histogram(breaks = seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),closed="right",colour="black") +
    geom_point(data = norm_table, aes(x = X, y = Y), colour = "blue") +
    annotate("text", label = paste("n = ",current_n,sep = ""), x = Inf, y = Inf,hjust=1.1,vjust=1.1, colour = ifelse(current_n < 30,'red','blue')) +
    annotate("text", label = paste("KS = ",formatC(current_ks, format = "e", digits = 1),sep = ""), x = Inf, y = Inf,hjust=1.1,vjust=2.2, colour = ifelse(current_ks < 0.05,'red','blue')) +
    annotate("text", label = paste("SW = ",formatC(current_sw, format = "e", digits = 1),sep = ""), x = Inf, y = Inf,hjust=1.1,vjust=3.3, colour = ifelse(current_sw < 0.05,'red','blue')) +
    labs(title = paste(taxonomic_group," richness",sep = "")) +
    scale_y_continuous(breaks = hist_counts) +
    ylab(label = "Count") +
    xlab(NULL) +
    theme_bw() +
    theme(text=element_text(size=16, family = "ArialMT"), plot.title = element_text(hjust = 0.5, size=16), axis.text = element_text(size=16))
  #Perform anova
  #Site is treated as a separate block with random effect
  set.seed(123456)
  observed_aov <- aov(Observed ~ Site+Treatment*Month*Community, data=get(ALPHA))
  observed_aov_summary <- summary(observed_aov)
  observed_aov_tukey <- TukeyHSD(observed_aov,ordered = TRUE)
  # get a table for tukey results to plot multiple comparisons
  observed_tukey_table <- as.data.frame(observed_aov_tukey$`Treatment:Month:Community`)
  observed_res <- observed_tukey_table$`p adj`
  names(observed_res) <- row.names(observed_tukey_table)
  observed_letters <- multcompLetters(observed_res)
  observed_letters_df <- data.frame("comp" = names(observed_letters$Letters), "Letter"=as.vector(observed_letters$Letters))
  observed_letters_df$Treatment <- gsub(":.+","",observed_letters_df$comp)
  observed_letters_df$Month <- gsub("^.+:{1}","",gsub(":{1}[A-z]+$","",observed_letters_df$comp))
  observed_letters_df$Community <- gsub(".+:","",observed_letters_df$comp)
  observed_letters_df$Treatment <- factor(observed_letters_df$Treatment, levels = c("C","SW"),ordered = F)
  observed_letters_df$Month <- factor(observed_letters_df$Month, levels = c("June","July","August"),ordered = T)
  observed_letters_df$Community <- factor(observed_letters_df$Community, levels = c("Total","Active"),ordered = F)
  observed_letters_df$Observed <- max(get(ALPHA)$Observed)*1.05
  #Create plot
  observed_plot <- ggplot(get(ALPHA), aes(x=Treatment, y=Observed)) +
    geom_point (aes(colour = Site),position = position_jitterdodge(jitter.width = 0.1))  +
    geom_boxplot(fill=NA,size=0.5) +
    geom_text(data = observed_letters_df, aes(label=Letter)) +
    theme_bw() +
    scale_y_continuous(limits = c(0,RoundTo(max(get(ALPHA)[,"Observed"])*1.05,multiple = 25,ceiling)),n.breaks = 6) +
    scale_color_manual(values=site_colours) +
    facet_nested(~Month+Community) +
    labs(title = paste(taxonomic_group," richness",sep = ""), x = NULL, y = "Richness" ,colour="Site") +
    theme(text=element_text(size=16, family = "ArialMT"), legend.position = "right", plot.title = element_text(hjust = 0.5, size=16), 
          panel.spacing=unit(0,"lines"),axis.text = element_text(size=16))
  #Export plot
  # ggsave(paste("Alpha_observed_",ALPHA,".svg",sep = ""),observed_plot,width = 2*84, height = 2*84, units = "mm")
  #Save results
  alpha_plots[paste("observed_",ALPHA,sep = "")] <- list(observed_plot)
  alpha_plots_normality[paste("observed_",ALPHA,sep = "")] <- list(observed_normality)
  alpha_anova[paste("observed_",ALPHA,sep = "")] <- list(observed_aov_summary)
  alpha_Tukey[paste("observed_",ALPHA,sep = "")] <- list(observed_aov_tukey)
  
  ##Simpson##
  #Check for normality
  current_values <- get(ALPHA)[,"Simpson"]
  current_n <- length(current_values)
  #Perform statistical tests
  current_ks <- ks.test(current_values,y = "pnorm")$`p.value`
  current_sw <- shapiro.test(current_values)$`p.value`
  #Create a normal distribution over the same distance
  norm_table <- data.frame(X = seq(min(current_values),max(current_values),length.out=100),
                           Y = dnorm(seq(min(current_values),max(current_values),length.out=100),
                                     mean = mean(current_values), 
                                     sd = sd(current_values)))
  #Adjust to the height of the graph
  scaling_factor <- max(hist(current_values, breaks=seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)/max(norm_table$Y)
  norm_table$Y <- norm_table$Y*scaling_factor
  #Plot histogram of the data
  hist_counts <- round(seq(range(hist(current_values, breaks = seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)[1],
                           range(hist(current_values, breaks = seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)[2],
                           length.out = 3))
  simpson_normality <- ggplot(get(ALPHA),aes(Simpson)) +
    geom_histogram(breaks = seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),closed="right",colour="black") +
    geom_point(data = norm_table, aes(x = X, y = Y), colour = "blue") +
    annotate("text", label = paste("n = ",current_n,sep = ""), x = Inf, y = Inf,hjust=1.1,vjust=1.1, colour = ifelse(current_n < 30,'red','blue')) +
    annotate("text", label = paste("KS = ",formatC(current_ks, format = "e", digits = 1),sep = ""), x = Inf, y = Inf,hjust=1.1,vjust=2.2, colour = ifelse(current_ks < 0.05,'red','blue')) +
    annotate("text", label = paste("SW = ",formatC(current_sw, format = "e", digits = 1),sep = ""), x = Inf, y = Inf,hjust=1.1,vjust=3.3, colour = ifelse(current_sw < 0.05,'red','blue')) +
    labs(title = paste(taxonomic_group,"  Simpson (1-D)",sep = "")) +
    scale_y_continuous(breaks = hist_counts) +
    ylab(label = "Count") +
    xlab(NULL) +
    theme_bw() +
    theme(text=element_text(size=16, family = "ArialMT"), plot.title = element_text(hjust = 0.5, size=16), axis.text = element_text(size=16))
  #Perform anova
  set.seed(123456)
  simpson_aov <- aov(Simpson ~ Site+Treatment*Month*Community, data=get(ALPHA))
  simpson_aov_summary <- summary(simpson_aov)
  simpson_aov_tukey <- TukeyHSD(simpson_aov,ordered = TRUE)
  # get a table for tukey results to plot multiple comparisons
  simpson_tukey_table <- as.data.frame(simpson_aov_tukey$`Treatment:Month:Community`)
  simpson_res <- simpson_tukey_table$`p adj`
  names(simpson_res) <- row.names(simpson_tukey_table)
  simpson_letters <- multcompLetters(simpson_res)
  simpson_letters_df <- data.frame("comp" = names(simpson_letters$Letters), "Letter"=as.vector(simpson_letters$Letters))
  simpson_letters_df$Treatment <- gsub(":.+","",simpson_letters_df$comp)
  simpson_letters_df$Month <- gsub("^.+:{1}","",gsub(":{1}[A-z]+$","",simpson_letters_df$comp))
  simpson_letters_df$Community <- gsub(".+:","",simpson_letters_df$comp)
  simpson_letters_df$Treatment <- factor(simpson_letters_df$Treatment, levels = c("C","SW"),ordered = F)
  simpson_letters_df$Month <- factor(simpson_letters_df$Month, levels = c("June","July","August"),ordered = T)
  simpson_letters_df$Community <- factor(simpson_letters_df$Community, levels = c("Total","Active"),ordered = F)
  simpson_letters_df$Simpson <- 1.05
  #Create plot
  simpson_plot <- ggplot(get(ALPHA), aes(x=Treatment, y=Simpson)) +
    geom_point (aes(colour = Site),position = position_jitterdodge(jitter.width = 0.1))  +
    geom_boxplot(fill=NA,size=0.5) +
    geom_text(data = simpson_letters_df, aes(label=Letter)) +
    theme_bw() +
    scale_y_continuous(limits = c(0,1.1),breaks = c(0.0,0.2,0.4,0.6,0.8,1.0)) +
    scale_color_manual(values=site_colours) +
    facet_nested(~Month+Community) +
    labs(title = paste(taxonomic_group," Simpson (1-D)",sep = ""), x = NULL, colour="Site") +
    theme(text=element_text(size=16, family = "ArialMT"), legend.position = "right", plot.title = element_text(hjust = 0.5, size=16), 
          panel.spacing=unit(0,"lines"),axis.text = element_text(size=16))
  #Export plot
  # ggsave(paste("Alpha_Simpson_",ALPHA,".svg",sep = ""),simpson_plot,width = 2*84, height = 2*84, units = "mm")
  #Save results
  alpha_plots[paste("Simpson_",ALPHA,sep = "")] <- list(simpson_plot)
  alpha_plots_normality[paste("Simpson_",ALPHA,sep = "")] <- list(simpson_normality)
  alpha_anova[paste("Simpson_",ALPHA,sep = "")] <- list(simpson_aov_summary)
  alpha_Tukey[paste("Simpson_",ALPHA,sep = "")] <- list(simpson_aov_tukey)
  
  ##Shannon##
  #Check for normality
  current_values <- get(ALPHA)[,"Shannon"]
  current_n <- length(current_values)
  #Perform statistical tests
  current_ks <- ks.test(current_values,y = "pnorm")$`p.value`
  current_sw <- shapiro.test(current_values)$`p.value`
  #Create a normal distribution over the same distance
  norm_table <- data.frame(X = seq(min(current_values),max(current_values),length.out=100),
                           Y = dnorm(seq(min(current_values),max(current_values),length.out=100),
                                     mean = mean(current_values), 
                                     sd = sd(current_values)))
  #Adjust to the height of the graph
  scaling_factor <- max(hist(current_values, breaks=seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)/max(norm_table$Y)
  norm_table$Y <- norm_table$Y*scaling_factor
  #Plot histogram of the data
  hist_counts <- round(seq(range(hist(current_values, breaks = seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)[1],
                           range(hist(current_values, breaks = seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),plot=F)$counts)[2],
                           length.out = 3))
  shannon_normality <- ggplot(get(ALPHA),aes(Shannon)) +
    geom_histogram(breaks = seq(min(current_values),max(current_values),length.out=nclass.Sturges(current_values)+1),closed="right",colour="black") +
    geom_point(data = norm_table, aes(x = X, y = Y), colour = "blue") +
    annotate("text", label = paste("n = ",current_n,sep = ""), x = Inf, y = Inf,hjust=1.1,vjust=1.1, colour = ifelse(current_n < 30,'red','blue')) +
    annotate("text", label = paste("KS = ",formatC(current_ks, format = "e", digits = 1),sep = ""), x = Inf, y = Inf,hjust=1.1,vjust=2.2, colour = ifelse(current_ks < 0.05,'red','blue')) +
    annotate("text", label = paste("SW = ",formatC(current_sw, format = "e", digits = 1),sep = ""), x = Inf, y = Inf,hjust=1.1,vjust=3.3, colour = ifelse(current_sw < 0.05,'red','blue')) +
    labs(title = paste(taxonomic_group,"  Shannon",sep = "")) +
    scale_y_continuous(breaks = hist_counts) +
    ylab(label = "Count") +
    xlab(NULL) +
    theme_bw() +
    theme(text=element_text(size=16, family = "ArialMT"), plot.title = element_text(hjust = 0.5, size=16), axis.text = element_text(size=16))
  #Perform anova
  set.seed(123456)
  shannon_aov <- aov(Shannon ~ Site+Treatment*Month*Community, data=get(ALPHA))
  shannon_aov_summary <- summary(shannon_aov)
  shannon_aov_tukey <- TukeyHSD(shannon_aov,ordered = TRUE)
  # get a table for tukey results to plot multiple comparisons
  shannon_tukey_table <- as.data.frame(shannon_aov_tukey$`Treatment:Month:Community`)
  shannon_res <- shannon_tukey_table$`p adj`
  names(shannon_res) <- row.names(shannon_tukey_table)
  shannon_letters <- multcompLetters(shannon_res)
  shannon_letters_df <- data.frame("comp" = names(shannon_letters$Letters), "Letter"=as.vector(shannon_letters$Letters))
  shannon_letters_df$Treatment <- gsub(":.+","",shannon_letters_df$comp)
  shannon_letters_df$Month <- gsub("^.+:{1}","",gsub(":{1}[A-z]+$","",shannon_letters_df$comp))
  shannon_letters_df$Community <- gsub(".+:","",shannon_letters_df$comp)
  shannon_letters_df$Treatment <- factor(shannon_letters_df$Treatment, levels = c("C","SW"),ordered = F)
  shannon_letters_df$Month <- factor(shannon_letters_df$Month, levels = c("June","July","August"),ordered = T)
  shannon_letters_df$Community <- factor(shannon_letters_df$Community, levels = c("Total","Active"),ordered = F)
  shannon_letters_df$Shannon <- max(get(ALPHA)$Shannon)*1.05
  #Create plot
  shannon_plot <- ggplot(get(ALPHA), aes(x=Treatment, y=Shannon)) +
    geom_point (aes(colour = Site),position = position_jitterdodge(jitter.width = 0.1))  +
    geom_boxplot(fill=NA,size=0.5) +
    geom_text(data = shannon_letters_df, aes(label=Letter)) +
    theme_bw() +
    scale_y_continuous(limits = c(0,RoundTo(max(get(ALPHA)[,"Shannon"]),multiple = 2,ceiling)),n.breaks = 6) +
    scale_color_manual(values=site_colours) +
    facet_nested(~Month+Community) +
    labs(title = paste(taxonomic_group," Shannon",sep = ""), x = NULL, colour="Site") +
    theme(text=element_text(size=16, family = "ArialMT"), legend.position = "right", plot.title = element_text(hjust = 0.5, size=16), 
          panel.spacing=unit(0,"lines"),axis.text = element_text(size=16))
  #Export plot
  # ggsave(paste("Alpha_shannon_",ALPHA,".svg",sep = ""),shannon_plot,width = 2*84, height = 2*84, units = "mm")
  #Save results
  alpha_plots[paste("shannon_",ALPHA,sep = "")] <- list(shannon_plot)
  alpha_plots_normality[paste("shannon_",ALPHA,sep = "")] <- list(shannon_normality)
  alpha_anova[paste("shannon_",ALPHA,sep = "")] <- list(shannon_aov_summary)
  alpha_Tukey[paste("shannon_",ALPHA,sep = "")] <- list(shannon_aov_tukey)
}

#Create a table of the anova test results
alpha_anova_table <- cbind(alpha_anova$observed_alphadiv_bact_mean[1][[1]],alpha_anova$observed_alphadiv_fung_mean[1][[1]])
alpha_anova_table[,"Comparison"] <- row.names(alpha_anova_table)
alpha_anova_table[nrow(alpha_anova_table)+1,] <- c(rep("Bacteria",5),rep("Fungi",5),"Comparison")
alpha_anova_table[nrow(alpha_anova_table)+1,] <- c(rep(colnames(alpha_anova_table[1:5]),2),"Comparison")
alpha_anova_table <- alpha_anova_table[c(10,11,1:9),c(11,1:10)]

# #Plots for report
# Site_legend <- get_legend(alpha_plots[["observed_alphadiv_bact_mean"]] + theme(legend.position = "bottom")+ guides(colour = guide_legend(nrow = 1)))
# report_observed <- plot_grid(
#   plot_grid(alpha_plots[["observed_alphadiv_bact_mean"]] + theme(legend.position = "none"),
#             alpha_plots[["observed_alphadiv_fung_mean"]] + theme(legend.position = "none"),
#             ncol = 2),
#   Site_legend,
#   ncol = 1, rel_heights = c(1,0.1))
# report_shannon <- plot_grid(
#   plot_grid(alpha_plots[["shannon_alphadiv_bact_mean"]] + theme(legend.position = "none"),
#             alpha_plots[["shannon_alphadiv_fung_mean"]] + theme(legend.position = "none"),
#             ncol = 2),
#   Site_legend,
#   ncol = 1, rel_heights = c(1,0.1))
# report_simpson <- plot_grid(
#   plot_grid(alpha_plots[["Simpson_alphadiv_bact_mean"]] + theme(legend.position = "none"),
#             alpha_plots[["Simpson_alphadiv_fung_mean"]] + theme(legend.position = "none"),
#             ncol = 2),
#   Site_legend,
#   ncol = 1, rel_heights = c(1,0.1))
# 
# ggsave("report_observed.svg",report_observed,width = 2*168, height = 2*84, units = "mm")
# ggsave("report_shannon.svg",report_shannon,width = 2*168, height = 2*84, units = "mm")
# ggsave("report_simpson.svg",report_simpson,width = 2*168, height = 2*84, units = "mm")
# ggsave("report_observed.jpg",report_observed,width = 2*168, height = 2*84, units = "mm")
# ggsave("report_shannon.jpg",report_shannon,width = 2*168, height = 2*84, units = "mm")
# ggsave("report_simpson.jpg",report_simpson,width = 2*168, height = 2*84, units = "mm")


#Want to pull out the significant variables

#pull out a table with the details of the significant variables
alpha_significant_variables <- data.frame()
for (STAT in names(alpha_anova)) {
  current <- as.data.frame((alpha_anova[[STAT]][[1]])[alpha_anova[[STAT]][[1]]$'Pr(>F)'<0.05,])
  current$STAT <- STAT
  current$Variable <- row.names(current)
  alpha_significant_variables <- rbind(alpha_significant_variables,current)
}
#Clean up the table to make it readable
alpha_significant_variables <- alpha_significant_variables[!alpha_significant_variables$Variable=="NA",]
row.names(alpha_significant_variables) <- 1:nrow(alpha_significant_variables)
alpha_significant_variables$Measurement <- gsub("_.+","",alpha_significant_variables$STAT)
alpha_significant_variables$Taxon <- gsub("_mean","", gsub(".+_alphadiv_","",alpha_significant_variables$STAT))
alpha_significant_variables[which(alpha_significant_variables$Taxon == "bact"),"Taxon"] <- "Bacteria"
alpha_significant_variables[which(alpha_significant_variables$Taxon == "fung"),"Taxon"] <- "Fungi"
alpha_significant_variables <- alpha_significant_variables[,c(9,8,7,1:6)]
#Exclude collection site
alpha_significant_variables$Variable <- trimws(alpha_significant_variables$Variable)
alpha_significant_variables <- alpha_significant_variables[!alpha_significant_variables$Variable=="Site",]
# #Export table
# write.table(subset(alpha_significant_variables,select = -STAT),"alpha_significant_variables.csv",quote = T,sep = ",",row.names = F,col.names = T)

#Save normality plot for internal use
#Observed
Alpha_normality <- plot_grid(alpha_plots_normality$observed_alphadiv_bact_mean + labs(tag="A"),
                                      alpha_plots_normality$observed_alphadiv_fung_mean + labs(tag="B"),
                                      alpha_plots_normality$Simpson_alphadiv_bact_mean + labs(tag="C"),
                                      alpha_plots_normality$Simpson_alphadiv_fung_mean + labs(tag="D"),
                                      alpha_plots_normality$shannon_alphadiv_bact_mean + labs(tag="E"),
                                      alpha_plots_normality$shannon_alphadiv_fung_mean + labs(tag="F"),
                                      ncol = 2)
#SVG
svg("Alpha_normality.svg", width=2*168*0.03937008, height=2*180*0.03937008)
Alpha_normality
dev.off()
#JPG
jpeg("Alpha_normality.jpg", width=2*168, height=2*180, units = "mm", res = 300)
Alpha_normality
dev.off()

###Beta diversity###

#Distances
#Use recommended Aitchison and Hellinger
# set.seed(123456)
# Dist_Aitchison_bact <- as.dist(avgdist(x = as.data.frame(t(phyloseq_bact@otu_table)), sample = min(colSums(phyloseq_bact@otu_table)), distfun = vegdist, dmethod = "robust.aitchison", iterations = 1000))
set.seed(123456)
Dist_Hellinger_bact <- as.dist(avgdist(x = as.data.frame(t(phyloseq_bact@otu_table)), sample = min(colSums(phyloseq_bact@otu_table)), distfun = vegdist, dmethod = "hellinger", iterations = 1000))
# set.seed(123456)
# Dist_Aitchison_fung <- as.dist(avgdist(x = as.data.frame(t(phyloseq_fung@otu_table)), sample = min(colSums(phyloseq_fung@otu_table)), distfun = vegdist, dmethod = "robust.aitchison", iterations = 1000))
set.seed(123456)
Dist_Hellinger_fung <- as.dist(avgdist(x = as.data.frame(t(phyloseq_fung@otu_table)), sample = min(colSums(phyloseq_fung@otu_table)), distfun = vegdist, dmethod = "hellinger", iterations = 1000))

#Calculate the ordinations
# set.seed(123456)
# ordination_Aitchison_bact <- ordinate(phyloseq_bact, method="NMDS", distance=Dist_Aitchison_bact)
set.seed(123456)
ordination_Hellinger_bact <- ordinate(phyloseq_bact, method="NMDS", distance=Dist_Hellinger_bact)
# set.seed(123456)
# ordination_Aitchison_fung <- ordinate(phyloseq_fung, method="NMDS", distance=Dist_Aitchison_fung)
set.seed(123456)
ordination_Hellinger_fung <- ordinate(phyloseq_fung, method="NMDS", distance=Dist_Hellinger_fung)

#Plotting
# NMDS_Aitchison_bact <- 
#   plot_ordination(phyloseq_bact, ordination_Aitchison_bact, color = "Treatment", shape = "Month") +
#   geom_point(size=3) + 
#   scale_colour_manual(values = c(treatment_colours)) + 
#   annotate(geom="text", x = max(ordination_Aitchison_bact$points[,"MDS1"]), y = max(ordination_Aitchison_bact$points[,"MDS2"])+0.1, label = paste("Stress: ", round(ordination_Aitchison_bact$stress, digits = 3)), vjust = "inward", hjust = "inward") +
#   facet_wrap(~Community) +
#   theme_bw() +
#   labs(title = "Bacteria: NMDS (Aitchison)") +
#   theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=16),axis.text = element_text(size=16))

NMDS_Hellinger_bact <- 
  plot_ordination(phyloseq_bact, ordination_Hellinger_bact, color = "Treatment", shape = "Month") +
  geom_point(size=3) + 
  scale_colour_manual(values = c(treatment_colours)) + 
  annotate(geom="text", x = max(ordination_Hellinger_bact$points[,"MDS1"]), y = max(ordination_Hellinger_bact$points[,"MDS2"])+0.1, label = paste("Stress: ", round(ordination_Hellinger_bact$stress, digits = 3)), vjust = "inward", hjust = "inward") +
  facet_wrap(~Community) +
  theme_bw() +
  labs(title = "Bacteria: NMDS (Hellinger)") +
  theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=16),axis.text = element_text(size=16))

# NMDS_Aitchison_fung <- 
#   plot_ordination(phyloseq_fung, ordination_Aitchison_fung, color = "Treatment", shape = "Month") +
#   geom_point(size=3) + 
#   scale_colour_manual(values = c(treatment_colours)) + 
#   annotate(geom="text", x = max(ordination_Aitchison_fung$points[,"MDS1"]), y = max(ordination_Aitchison_fung$points[,"MDS2"])+0.1, label = paste("Stress: ", round(ordination_Aitchison_fung$stress, digits = 3)), vjust = "inward", hjust = "inward") +
#   facet_wrap(~Community) +
#   theme_bw() +
#   labs(title = "Fungi: NMDS (Aitchison)") +
#   theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=16),axis.text = element_text(size=16))

NMDS_Hellinger_fung <- 
  plot_ordination(phyloseq_fung, ordination_Hellinger_fung, color = "Treatment", shape = "Month") +
  geom_point(size=3) + 
  scale_colour_manual(values = c(treatment_colours)) + 
  annotate(geom="text", x = max(ordination_Hellinger_fung$points[,"MDS1"]), y = max(ordination_Hellinger_fung$points[,"MDS2"])+0.1, label = paste("Stress: ", round(ordination_Hellinger_fung$stress, digits = 3)), vjust = "inward", hjust = "inward") +
  facet_wrap(~Community) +
  theme_bw() +
  labs(title = "Fungi: NMDS (Hellinger)") +
  theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=16),axis.text = element_text(size=16))

#Figures showing the ordination of sites instead of months
NMDS_Hellinger_bact_site <- 
  plot_ordination(phyloseq_bact, ordination_Hellinger_bact, color = "Site", shape = "Treatment") +
  geom_point(size=3) + 
  scale_colour_manual(values = c(site_colours)) + 
  annotate(geom="text", x = max(ordination_Hellinger_bact$points[,"MDS1"]), y = max(ordination_Hellinger_bact$points[,"MDS2"])+0.1, label = paste("Stress: ", round(ordination_Hellinger_bact$stress, digits = 3)), vjust = "inward", hjust = "inward") +
  facet_wrap(~Community) +
  theme_bw() +
  labs(title = "Bacteria: NMDS (Hellinger)") +
  theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=16),axis.text = element_text(size=16))

NMDS_Hellinger_fung_site <- 
  plot_ordination(phyloseq_fung, ordination_Hellinger_fung, color = "Site", shape = "Treatment") +
  geom_point(size=3) + 
  scale_colour_manual(values = c(site_colours)) + 
  annotate(geom="text", x = max(ordination_Hellinger_fung$points[,"MDS1"]), y = max(ordination_Hellinger_fung$points[,"MDS2"])+0.1, label = paste("Stress: ", round(ordination_Hellinger_fung$stress, digits = 3)), vjust = "inward", hjust = "inward") +
  facet_wrap(~Community) +
  theme_bw() +
  labs(title = "Fungi: NMDS (Hellinger)") +
  theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=16),axis.text = element_text(size=16))

# #Export figures
# ggsave("NMDS_Aitchison_bact.svg",NMDS_Aitchison_bact,width = 2*84, height = 2*84, units = "mm")
# ggsave("NMDS_Hellinger_bact.svg",NMDS_Hellinger_bact,width = 2*84, height = 2*84, units = "mm")
# ggsave("NMDS_Aitchison_fung.svg",NMDS_Aitchison_fung,width = 2*84, height = 2*84, units = "mm")
# ggsave("NMDS_Hellinger_fung.svg",NMDS_Hellinger_fung,width = 2*84, height = 2*84, units = "mm")

# #Plots for report
# Beta_legend <- get_legend(NMDS_Hellinger_bact + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 1)))
# report_beta_hellinger <- plot_grid(
#   plot_grid(NMDS_Hellinger_bact + theme(legend.position = "none"),
#             NMDS_Hellinger_fung + theme(legend.position = "none"),
#             ncol = 2),
#   Beta_legend,
#   ncol = 1, rel_heights = c(1,0.1))
# report_beta_aitchison <- plot_grid(
#   plot_grid(NMDS_Aitchison_bact + theme(legend.position = "none"),
#             NMDS_Aitchison_fung + theme(legend.position = "none"),
#             ncol = 2),
#   Beta_legend,
#   ncol = 1, rel_heights = c(1,0.1))
# ggsave("report_beta_hellinger.jpg",report_beta_hellinger,width = 2*168, height = 2*84, units = "mm")
# ggsave("report_beta_aitchison.jpg",report_beta_aitchison,width = 2*168, height = 2*84, units = "mm")

#Check for differences in composition
#Using Hellinger as it had the lowest stress in the NMDS

#Checks whether the samples have the same centroid. Null hypothesis is same centroid.
#https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
#R2 shows how much is explained, residuals are unexplained variation.
sampledf_bact <- data.frame(sample_data(phyloseq_bact))
#Create interaction columns
sampledf_bact$`Treatment:Month` <- interaction(sampledf_bact$Treatment, sampledf_bact$Month)
sampledf_bact$`Treatment:Community` <- interaction(sampledf_bact$Treatment, sampledf_bact$Community)
sampledf_bact$`Month:Community` <- interaction(sampledf_bact$Month, sampledf_bact$Community)
sampledf_bact$`Treatment:Month:Community` <- interaction(sampledf_bact$Treatment, sampledf_bact$Month, sampledf_bact$Community)
#Perform permanova
set.seed("123456")
adonis_hellinger_bact <- adonis2(Dist_Hellinger_bact ~ Site+Treatment*Month*Community, data = sampledf_bact, permutations = 1000, by = "terms")

sampledf_fung <- data.frame(sample_data(phyloseq_fung))
#Create interaction columns
sampledf_fung$`Treatment:Month` <- interaction(sampledf_fung$Treatment, sampledf_fung$Month)
sampledf_fung$`Treatment:Community` <- interaction(sampledf_fung$Treatment, sampledf_fung$Community)
sampledf_fung$`Month:Community` <- interaction(sampledf_fung$Month, sampledf_fung$Community)
sampledf_fung$`Treatment:Month:Community` <- interaction(sampledf_fung$Treatment, sampledf_fung$Month, sampledf_fung$Community)
#Perform permanova
set.seed("123456")
adonis_hellinger_fung <- adonis2(Dist_Hellinger_fung ~ Site+Treatment*Month*Community, data = sampledf_fung, permutations = 1000, by = "terms")

#Create a table of the permanova test results
beta_perm_table <- cbind(as.data.frame(adonis_hellinger_bact),as.data.frame(adonis_hellinger_fung))
beta_perm_table[,"Comparison"] <- row.names(beta_perm_table)
beta_perm_table[,"Beta Dispersion (p)"] <- row.names(beta_perm_table)
beta_perm_table[,"Beta Dispersion (p).1"] <- row.names(beta_perm_table)
beta_perm_table[nrow(beta_perm_table)+1,] <- c(rep("Bacteria",5),rep("Fungi",5),"Comparison","Bacteria","Fungi")
beta_perm_table[nrow(beta_perm_table)+1,] <- c(rep(colnames(beta_perm_table[1:5]),2),"Comparison", "Beta Dispersion (p)","Beta Dispersion (p)")
beta_perm_table <- beta_perm_table[c(11,12,1:10),c(11,1:5,12,6:10,13)]
beta_perm_table[c(11,12),c(7,13)] <- NA

#Homogeneity of dispersion tests
for (GROUP in beta_perm_table[3:10,"Comparison"]) {
  beta_b <- betadisper(Dist_Hellinger_bact, sampledf_bact[,GROUP])
  beta_f <- betadisper(Dist_Hellinger_fung, sampledf_fung[,GROUP])
  set.seed("123456")
  beta_perm_table[beta_perm_table$Comparison==GROUP,7] <- permutest(beta_b)$tab$`Pr(>F)`[1]
  set.seed("123456")
  beta_perm_table[beta_perm_table$Comparison==GROUP,13] <- permutest(beta_f)$tab$`Pr(>F)`[1]
}

#Adjust betadispers results for multiple testing
beta_perm_table[3:10,7] <- p.adjust(beta_perm_table[3:10,7], method = "fdr")
beta_perm_table[3:10,13] <- p.adjust(beta_perm_table[3:10,13], method = "fdr")

# #export
# write.table(beta_perm_table,"beta_perm_table.csv",quote = T,row.names = F,sep = ",")

# ###Cluster the samples###
# 
# #Split by DNA and RNA
# phyloseq_bact_DNA <- subset_samples(phyloseq_bact, Type=="DNA")
# phyloseq_bact_RNA <- subset_samples(phyloseq_bact, Type=="RNA")
# phyloseq_fung_DNA <- subset_samples(phyloseq_fung, Type=="DNA")
# phyloseq_fung_RNA <- subset_samples(phyloseq_fung, Type=="RNA")
# 
# #Distances
# #Hellinger performed best in ordinations
# Dist_Hellinger_bact_DNA <- as.dist(avgdist(x = as.data.frame(t(phyloseq_bact_DNA@otu_table)), sample = min(colSums(phyloseq_bact_DNA@otu_table)), distfun = vegdist, dmethod = "hellinger", iterations = 100))
# Dist_Hellinger_bact_RNA <- as.dist(avgdist(x = as.data.frame(t(phyloseq_bact_RNA@otu_table)), sample = min(colSums(phyloseq_bact_RNA@otu_table)), distfun = vegdist, dmethod = "hellinger", iterations = 100))
# Dist_Hellinger_fung_DNA <- as.dist(avgdist(x = as.data.frame(t(phyloseq_fung_DNA@otu_table)), sample = min(colSums(phyloseq_fung_DNA@otu_table)), distfun = vegdist, dmethod = "hellinger", iterations = 100))
# Dist_Hellinger_fung_RNA <- as.dist(avgdist(x = as.data.frame(t(phyloseq_fung_RNA@otu_table)), sample = min(colSums(phyloseq_fung_RNA@otu_table)), distfun = vegdist, dmethod = "hellinger", iterations = 100))
# 
# #Find the best Linkage Method to use for clustering 
# #define linkage methods
# m <- c( "average", "single", "complete", "ward")
# names(m) <- c( "average", "single", "complete", "ward")
# #function to compute agglomerative coefficient
# ac_bact_DNA <- function(x) {
#   agnes(Dist_Hellinger_bact_DNA, diss = TRUE, method = x)$ac
# }
# ac_bact_RNA <- function(x) {
#   agnes(Dist_Hellinger_bact_RNA, diss = TRUE, method = x)$ac
# }
# ac_fung_DNA <- function(x) {
#   agnes(Dist_Hellinger_fung_DNA, diss = TRUE, method = x)$ac
# }
# ac_fung_RNA <- function(x) {
#   agnes(Dist_Hellinger_fung_RNA, diss = TRUE, method = x)$ac
# }
# #calculate agglomerative coefficient for each clustering linkage method
# sapply(m, ac_bact_DNA)
# sapply(m, ac_bact_RNA)
# sapply(m, ac_fung_DNA)
# sapply(m, ac_fung_RNA)
# #Ward has highest ac in all cases, so best method to use
# 
# #Cluster samples
# final_clust_bact_DNA <- hclust(Dist_Hellinger_bact_DNA, method = "ward.D2")
# final_clust_bact_RNA <- hclust(Dist_Hellinger_bact_RNA, method = "ward.D2")
# final_clust_fung_DNA <- hclust(Dist_Hellinger_fung_DNA, method = "ward.D2")
# final_clust_fung_RNA <- hclust(Dist_Hellinger_fung_RNA, method = "ward.D2")
# 
# #Determine the Optimal Number of Clusters by choosing the one returned most often by a list of 20 indices
# inds <- c("kl","ch","hartigan","cindex","db","silhouette","duda","pseudot2","ratkowsky","ball","ptbiserial","gap","frey","mcclain","gamma","gplus","tau","dunn","sdindex","sdbw")
# #In previous agnes function, the "ward" option is equivalent to "ward.D2" in the hclust function. So will use the "ward.D2" method in NbClust as well for equivalence.
# #https://link.springer.com/article/10.1007/s00357-014-9161-z
# restable_bact_DNA <- data.frame()
# for (i in (1:length(inds))) {
#   restable1 <- data.frame()
#   restable1[1,1] <- inds[i]
#   restable1[1,2] <- NbClust(as.data.frame(t(otu_table(phyloseq_bact_DNA@otu_table,taxa_are_rows = TRUE))),diss=Dist_Hellinger_bact_DNA,distance = NULL,method="ward.D2",index=inds[i])$Best.nc[1]
#   restable_bact_DNA <- rbind(restable_bact_DNA,restable1)
# }
# restable_bact_RNA <- data.frame()
# for (i in (1:length(inds))) {
#   restable1 <- data.frame()
#   restable1[1,1] <- inds[i]
#   restable1[1,2] <- NbClust(as.data.frame(t(otu_table(phyloseq_bact_RNA@otu_table,taxa_are_rows = TRUE))),diss=Dist_Hellinger_bact_RNA,distance = NULL,method="ward.D2",index=inds[i])$Best.nc[1]
#   restable_bact_RNA <- rbind(restable_bact_RNA,restable1)
# }
# restable_fung_DNA <- data.frame()
# for (i in (c(1:11,13:length(inds)))) { #Trouble with "gap" - can't allocate 170GB vector.
#   restable1 <- data.frame()
#   restable1[1,1] <- inds[i]
#   restable1[1,2] <- NbClust(as.data.frame(t(otu_table(phyloseq_fung_DNA@otu_table,taxa_are_rows = TRUE))),diss=Dist_Hellinger_fung_DNA,distance = NULL,method="ward.D2",index=inds[i])$Best.nc[1]
#   restable_fung_DNA <- rbind(restable_fung_DNA,restable1)
# }
# restable_fung_RNA <- data.frame()
# for (i in (c(1:11,13:length(inds)))) { #Trouble with "gap" - can't allocate 170GB vector.
#   restable1 <- data.frame()
#   restable1[1,1] <- inds[i]
#   restable1[1,2] <- NbClust(as.data.frame(t(otu_table(phyloseq_fung_RNA@otu_table,taxa_are_rows = TRUE))),diss=Dist_Hellinger_fung_RNA,distance = NULL,method="ward.D2",index=inds[i])$Best.nc[1]
#   restable_fung_RNA <- rbind(restable_fung_RNA,restable1)
# }
# 
# #plot histogram
# Clustering_inds_bact_DNA <- ggplot(restable_bact_DNA, aes(x=V2)) +
#   geom_histogram(bins = max(restable_bact_DNA$V2),col="black",fill="dark grey",binwidth = 1) +
#   scale_x_continuous(n.breaks = max(restable_bact_DNA$V2)) +
#   scale_y_continuous(breaks = c(0,2,4,6,8,10)) +
#   theme_classic() +
#   labs(title = "Optimal number of clusters: Bacteria (DNA)", y="Count", x = "Clusters") +
#   theme(plot.title = element_text(hjust = 0.5))
# Clustering_inds_bact_RNA <- ggplot(restable_bact_RNA, aes(x=V2)) +
#   geom_histogram(bins = max(restable_bact_RNA$V2),col="black",fill="dark grey",binwidth = 1) +
#   scale_x_continuous(n.breaks = max(restable_bact_RNA$V2)) +
#   scale_y_continuous(breaks = c(0,2,4,6,8,10)) +
#   theme_classic() +
#   labs(title = "Optimal number of clusters: Bacteria (RNA)", y="Count", x = "Clusters") +
#   theme(plot.title = element_text(hjust = 0.5))
# Clustering_inds_fung_DNA <- ggplot(restable_fung_DNA, aes(x=V2)) +
#   geom_histogram(bins = max(restable_fung_DNA$V2),col="black",fill="dark grey",binwidth = 1) +
#   scale_x_continuous(n.breaks = max(restable_fung_DNA$V2)) +
#   scale_y_continuous(breaks = c(0,2,4,6,8,10)) +
#   theme_classic() +
#   labs(title = "Optimal number of clusters: Fungi (DNA)", y="Count", x = "Clusters") +
#   theme(plot.title = element_text(hjust = 0.5))
# Clustering_inds_fung_RNA <- ggplot(restable_fung_RNA, aes(x=V2)) +
#   geom_histogram(bins = max(restable_fung_RNA$V2),col="black",fill="dark grey",binwidth = 1) +
#   scale_x_continuous(n.breaks = max(restable_fung_RNA$V2)) +
#   scale_y_continuous(breaks = c(0,2,4,6,8,10)) +
#   theme_classic() +
#   labs(title = "Optimal number of clusters: Fungi (RNA)", y="Count", x = "Clusters") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# #Colour clusters by vector is not officially supported but seems to work for now.
# bact_clust_colours_DNA <- treatment_colours[gsub("_.+","",final_clust_bact_DNA$labels[final_clust_bact_DNA$order])]
# bact_clust_colours_RNA <- treatment_colours[gsub("_.+","",final_clust_bact_RNA$labels[final_clust_bact_RNA$order])]
# fung_clust_colours_DNA <- treatment_colours[gsub("_.+","",final_clust_fung_DNA$labels[final_clust_fung_DNA$order])]
# fung_clust_colours_RNA <- treatment_colours[gsub("_.+","",final_clust_fung_RNA$labels[final_clust_fung_RNA$order])]
# names(bact_clust_colours_DNA) <- final_clust_bact_DNA$labels[final_clust_bact_DNA$order]
# names(bact_clust_colours_RNA) <- final_clust_bact_RNA$labels[final_clust_bact_RNA$order]
# names(fung_clust_colours_DNA) <- final_clust_fung_DNA$labels[final_clust_fung_DNA$order]
# names(fung_clust_colours_RNA) <- final_clust_fung_RNA$labels[final_clust_fung_RNA$order]
# 
# #Perform bootstrapping of clusters
# set.seed(123456)
# Bootstrap_clust_bact_DNA <- Bclust(t(as.data.frame(phyloseq_bact_DNA@otu_table)), FUN = function(.x) hclust(as.dist(avgdist(x = .x, sample = min(rowSums(.x)), distfun = vegdist, dmethod = "hellinger", iterations = 100)), method = "ward.D2"), iter=1000,mc.cores = 7)
# set.seed(123456)
# Bootstrap_clust_bact_RNA <- Bclust(t(as.data.frame(phyloseq_bact_RNA@otu_table)), FUN = function(.x) hclust(as.dist(avgdist(x = .x, sample = min(rowSums(.x)), distfun = vegdist, dmethod = "hellinger", iterations = 100)), method = "ward.D2"), iter=1000,mc.cores = 7)
# set.seed(123456)
# Bootstrap_clust_fung_DNA <- Bclust(t(as.data.frame(phyloseq_fung_DNA@otu_table)), FUN = function(.x) hclust(as.dist(avgdist(x = .x, sample = min(rowSums(.x)), distfun = vegdist, dmethod = "hellinger", iterations = 100)), method = "ward.D2"), iter=1000,mc.cores = 7)
# set.seed(123456)
# Bootstrap_clust_fung_RNA <- Bclust(t(as.data.frame(phyloseq_fung_RNA@otu_table)), FUN = function(.x) hclust(as.dist(avgdist(x = .x, sample = min(rowSums(.x)), distfun = vegdist, dmethod = "hellinger", iterations = 100)), method = "ward.D2"), iter=1000,mc.cores = 7)
# 
# #Get labels and co-ordinates 
# #Not how this is meant to be used but it seems to work
# plot(Bootstrap_clust_bact_DNA)
# Bootstrap_labels_bact_DNA <- Bclabels(hcl = final_clust_bact_DNA,values = Bootstrap_clust_bact_DNA$values)
# plot(Bootstrap_clust_bact_RNA)
# Bootstrap_labels_bact_RNA <- Bclabels(hcl = final_clust_bact_RNA,values = Bootstrap_clust_bact_RNA$values)
# plot(Bootstrap_clust_fung_DNA)
# Bootstrap_labels_fung_DNA <- Bclabels(hcl = final_clust_fung_DNA,values = Bootstrap_clust_fung_DNA$values)
# plot(Bootstrap_clust_fung_RNA)
# Bootstrap_labels_fung_RNA <- Bclabels(hcl = final_clust_fung_RNA,values = Bootstrap_clust_fung_RNA$values)
# 
# #Final set up
# dendro_clust_bact_DNA <- dendro_data(final_clust_bact_DNA)
# dendro_clust_bact_RNA <- dendro_data(final_clust_bact_RNA)
# dendro_clust_fung_DNA <- dendro_data(final_clust_fung_DNA)
# dendro_clust_fung_RNA <- dendro_data(final_clust_fung_RNA)
# dendro_clust_bact_labels_DNA <- as.data.frame(cbind(Bootstrap_labels_bact_DNA$coords[,"x"],cbind(Bootstrap_labels_bact_DNA$coords[,"y"],Bootstrap_labels_bact_DNA$labels)))
# dendro_clust_bact_labels_RNA <- as.data.frame(cbind(Bootstrap_labels_bact_RNA$coords[,"x"],cbind(Bootstrap_labels_bact_RNA$coords[,"y"],Bootstrap_labels_bact_RNA$labels)))
# dendro_clust_fung_labels_DNA <- as.data.frame(cbind(Bootstrap_labels_fung_DNA$coords[,"x"],cbind(Bootstrap_labels_fung_DNA$coords[,"y"],Bootstrap_labels_fung_DNA$labels)))
# dendro_clust_fung_labels_RNA <- as.data.frame(cbind(Bootstrap_labels_fung_RNA$coords[,"x"],cbind(Bootstrap_labels_fung_RNA$coords[,"y"],Bootstrap_labels_fung_RNA$labels)))
# colnames(dendro_clust_bact_labels_DNA) <- c("x","y","label")
# colnames(dendro_clust_bact_labels_RNA) <- c("x","y","label")
# colnames(dendro_clust_fung_labels_DNA) <- c("x","y","label")
# colnames(dendro_clust_fung_labels_RNA) <- c("x","y","label")
# 
# #Plot clustering
# #Dendrocut needs to be tried manually
# Cluster_graph_bact_DNA <- 
#   ggdendrogram(dendro_clust_bact_DNA) + 
#   geom_dendro(final_clust_bact_DNA, dendrocut=32, groupCols=c(colorRampPalette(brewer.pal(8, "Dark2"))(14),"black")) +
#   geom_text(data = dendro_clust_bact_labels_DNA, aes(x = x, y = y-0.05, label = paste(round(label*100),"%",sep="")), size = 2, vjust = 0) +
#   labs(y = "", x = "") +
#   scale_y_continuous(breaks = NULL) +
#   coord_cartesian(xlim = c(0, 41)) +
#   theme_classic() +
#   theme(text=element_text(size=16, family = "ArialMT"), axis.text.x = element_text(angle = 90, vjust = 0.5, colour = bact_clust_colours_DNA))
# Cluster_graph_bact_RNA <- 
#   ggdendrogram(dendro_clust_bact_RNA) + 
#   geom_dendro(final_clust_bact_RNA, dendrocut=32, groupCols=c(colorRampPalette(brewer.pal(8, "Dark2"))(14),"black")) +
#   geom_text(data = dendro_clust_bact_labels_RNA, aes(x = x, y = y-0.05, label = paste(round(label*100),"%",sep="")), size = 2, vjust = 0) +
#   labs(y = "", x = "") +
#   scale_y_continuous(breaks = NULL) +
#   coord_cartesian(xlim = c(0, 41)) +
#   theme_classic() +
#   theme(text=element_text(size=16, family = "ArialMT"), axis.text.x = element_text(angle = 90, vjust = 0.5, colour = bact_clust_colours_RNA))
# Cluster_graph_fung_DNA <- 
#   ggdendrogram(dendro_clust_fung_DNA) + 
#   geom_dendro(final_clust_fung_DNA, dendrocut=19, groupCols=c(colorRampPalette(brewer.pal(8, "Dark2"))(14),"black")) +
#   geom_text(data = dendro_clust_fung_labels_DNA, aes(x = x, y = y-0.05, label = paste(round(label*100),"%",sep="")), size = 2, vjust = 0) +
#   labs(y = "", x = "") +
#   scale_y_continuous(breaks = NULL) +
#   coord_cartesian(xlim = c(0, 41)) +
#   theme_classic() +
#   theme(text=element_text(size=16, family = "ArialMT"), axis.text.x = element_text(angle = 90, vjust = 0.5, colour = fung_clust_colours_DNA))
# Cluster_graph_fung_RNA <- 
#   ggdendrogram(dendro_clust_fung_RNA) + 
#   geom_dendro(final_clust_fung_RNA, dendrocut=18, groupCols=c(colorRampPalette(brewer.pal(8, "Dark2"))(14),"black")) +
#   geom_text(data = dendro_clust_fung_labels_RNA, aes(x = x, y = y-0.05, label = paste(round(label*100),"%",sep="")), size = 2, vjust = 0) +
#   labs(y = "", x = "") +
#   scale_y_continuous(breaks = NULL) +
#   coord_cartesian(xlim = c(0, 41)) +
#   theme_classic() +
#   theme(text=element_text(size=16, family = "ArialMT"), axis.text.x = element_text(angle = 90, vjust = 0.5, colour = fung_clust_colours_RNA))
# 
# #Export Clusters
# 
# ggsave("Cluster_graph_bact_DNA.svg",Cluster_graph_bact_DNA,width = 2*84, height = 2*84, units = "mm")
# ggsave("Cluster_graph_bact_RNA.svg",Cluster_graph_bact_RNA,width = 2*84, height = 2*84, units = "mm")
# ggsave("Cluster_graph_fung_DNA.svg",Cluster_graph_fung_DNA,width = 2*84, height = 2*84, units = "mm")
# ggsave("Cluster_graph_fung_RNA.svg",Cluster_graph_fung_RNA,width = 2*84, height = 2*84, units = "mm")

###Taxonomic Composition###

#Get average relative abundance over 1000 rarefactions
set.seed(123456)
phyloseq_bact_raref_RA <- phyloseq_mult_raref_avg(phyloseq_bact,min(colSums(phyloseq_bact@otu_table)),iter = 1000)
set.seed(123456)
phyloseq_fung_raref_RA <- phyloseq_mult_raref_avg(phyloseq_fung,min(colSums(phyloseq_fung@otu_table)),iter = 1000)

#Create plots and tables of community composition
Taxonomy_RA_plots <- list()
Taxonomy_RA_tables <- list()
Taxonomy_RA_tables_sum <- list()
Taxonomy_RA_plot_tables <- list()
Taxonomy_RA_plots_gabriele <- list()
for (BACFUN in c("bact","fung")) {
  for (LEVEL in c("Phylum","Class","Order","Family","Genus")) {
    
    #Plot abundance (merge <1%)
    physeq_current <- tax_glom(get(paste("phyloseq_",BACFUN,"_raref_RA",sep = "")), taxrank = LEVEL,NArm = FALSE)

    physeq_current@tax_table[taxa_names(filter_taxa(physeq_current, function(x) mean(x) < 0.01, TRUE)),] <- "Rare taxa"
    physeq_current <- merge_taxa(physeq_current,taxa_names(filter_taxa(physeq_current, function(x) mean(x) < 0.01, TRUE)))
    colours_current <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tax_table(physeq_current)[,LEVEL])))
    names(colours_current) <- sort(unique(physeq_current@tax_table[,LEVEL]),decreasing = F)
    plot_current <-
      plot_bar(physeq_current, fill = LEVEL, x = "Site") +
      geom_bar(stat = "identity", position = "fill", colour = "black") +
      labs(x= element_blank(), y = "Relative Abundance") +
      facet_wrap(~Treatment*Community*Month,ncol = 12, scales = "free_x") +
      scale_fill_manual(values = colours_current, labels = names(colours_current), guide = guide_legend(ncol = 1)) +
      theme_bw() +
      theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=16), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), ,axis.text = element_text(size=16))
    
    #Create RA table
    RA_table_current <- as.data.frame(physeq_current@otu_table)
    RA_table_current$Taxon <- as.character(physeq_current@tax_table[row.names(RA_table_current),LEVEL])
    RA_table_current <- RA_table_current[,c(ncol(RA_table_current),1:(ncol(RA_table_current)-1))]
    
    #Create summary RA table for groupings'
    RA_table_current_sum <- as.data.frame(matrix(nrow = nrow(RA_table_current),ncol = 18, dimnames = list(NULL,c("Taxon","Overall",unique(paste(physeq_current@sam_data$Treatment,physeq_current@sam_data$Community,sep = "-")),unique(paste(physeq_current@sam_data$Treatment,physeq_current@sam_data$Community,physeq_current@sam_data$Month,sep = "-"))))))
    RA_table_current_sum$Taxon <- RA_table_current$Taxon
    RA_table_current_sum$Overall <- rowMeans(physeq_current@otu_table)
    for (COMBINATION in unique(paste(physeq_current@sam_data$Treatment,physeq_current@sam_data$Community,sep = "-"))) {
      RA_table_current_sum[,COMBINATION] <- rowMeans(physeq_current@otu_table[,physeq_current@sam_data$Treatment==strsplit(COMBINATION,"-")[[1]][1] &
                                                                                physeq_current@sam_data$Community==strsplit(COMBINATION,"-")[[1]][2]])

    }
    for (COMBINATION in unique(paste(physeq_current@sam_data$Treatment,physeq_current@sam_data$Community,physeq_current@sam_data$Month,sep = "-"))) {
      RA_table_current_sum[,COMBINATION] <- rowMeans(physeq_current@otu_table[,physeq_current@sam_data$Treatment==strsplit(COMBINATION,"-")[[1]][1] &
                                                                                physeq_current@sam_data$Community==strsplit(COMBINATION,"-")[[1]][2] &
                                                                                physeq_current@sam_data$Month==strsplit(COMBINATION,"-")[[1]][3]])
    }
    
    #Create summary RA table for groupings that is suitable for ggplot
    RA_plot_table_current <- data.frame()
    for (COMBINATION in unique(paste(physeq_current@sam_data$Treatment,physeq_current@sam_data$Community,sep = "-"))) {
      RA_plot_table_current_partial <- as.data.frame(matrix(nrow = nrow(RA_table_current),ncol = 6,dimnames = list(NULL,c("Taxon","Sample","Treatment","Community","RA","sd"))))
      RA_plot_table_current_partial$Taxon <- RA_table_current$Taxon
      Taxon_NAs <- which(is.na(RA_plot_table_current_partial$Taxon))
      if (length(Taxon_NAs) >= 1) {
        for (UNKNOWNS in 1:length(Taxon_NAs)) {
          RA_plot_table_current_partial[Taxon_NAs[UNKNOWNS],"Taxon"] <- paste("Unidentified_",UNKNOWNS,sep = "")
        }
      }
      RA_plot_table_current_partial$Sample <- COMBINATION
      RA_plot_table_current_partial$Treatment <- strsplit(COMBINATION,"-")[[1]][1]
      RA_plot_table_current_partial$Community <- strsplit(COMBINATION,"-")[[1]][2]
      RA_plot_table_current_partial$RA <- rowMeans(physeq_current@otu_table[,physeq_current@sam_data$Treatment==strsplit(COMBINATION,"-")[[1]][1] &
                                                                           physeq_current@sam_data$Community==strsplit(COMBINATION,"-")[[1]][2]])
      RA_plot_table_current_partial$sd <- apply(physeq_current@otu_table[,physeq_current@sam_data$Treatment==strsplit(COMBINATION,"-")[[1]][1] &
                                                                           physeq_current@sam_data$Community==strsplit(COMBINATION,"-")[[1]][2]],
                                                ,MARGIN = 1, FUN = sd)
      RA_plot_table_current <- rbind(RA_plot_table_current,RA_plot_table_current_partial)
    }
    #Neaten it up
    colnames(RA_plot_table_current) <- c("Taxon","Sample","Treatment","Community","RA","sd")
    RA_plot_table_current$Taxon <- factor(RA_plot_table_current$Taxon, levels = c(sort(unique(RA_plot_table_current$Taxon))[-c(which(sort(unique(RA_plot_table_current$Taxon))=="Rare taxa"),which(grepl("Unknown",sort(unique(RA_plot_table_current$Taxon)))))], "Rare taxa", sort(unique(RA_plot_table_current$Taxon))[which(grepl("Unknown",sort(unique(RA_plot_table_current$Taxon))))]))
    
    #Find the top 15 taxa by mean RA
    #In order to clean up some visualisations
    Top_taxa <- aggregate(RA_plot_table_current[,c("Taxon","RA")],.~Taxon,FUN = mean)
    Top_taxa <- Top_taxa[order(-Top_taxa$RA),]
    if (nrow(Top_taxa)>15) {
      Top_taxa <- Top_taxa[1:15,]
    }
    
    #Try out Gabriele's style of presentation
    plot_gabriele_current <- ggplot(RA_plot_table_current[RA_plot_table_current$Taxon%in%Top_taxa$Taxon,]) +
      geom_errorbar(aes(group=Sample, y = Taxon, xmin=RA, xmax = RA+sd), position=position_dodge2(reverse = TRUE),color="black", width=0.6) +
      geom_col(aes(fill=Sample, y=Taxon, x=RA), position=position_dodge2(reverse = TRUE),color="black", width=0.6) +
      scale_y_discrete(limits = rev) +
      labs(x= "Relative abundance", y = element_blank()) +
      scale_fill_manual(values = group_colours) +
      theme_bw() +
      theme(text=element_text(size=16),axis.text = element_text(size=16),legend.position = "bottom", panel.grid.major = element_line(color = "gray80", linewidth = 0.5))

    #Export results
    # ggsave(paste("taxonomy_RA_",BACFUN,"_",LEVEL,".eps",sep = ""),plot_current,width = 2*168, height = 2*168, units = "mm")
    # write.table(RA_table_current, file = paste("RA_table_",BACFUN,"_",LEVEL,".csv",sep = ""),quote = T,row.names = F, col.names = T, sep = ",")
    write.table(RA_table_current_sum, file = paste("RA_table_",BACFUN,"_",LEVEL,"_sum.csv",sep = ""),quote = T,row.names = F, col.names = T, sep = ",")
    
    #Save to list
    Taxonomy_RA_plots[paste(BACFUN,"_",LEVEL,sep = "")] <- list(plot_current)
    Taxonomy_RA_tables[paste(BACFUN,"_",LEVEL,sep = "")] <- list(RA_table_current)
    Taxonomy_RA_tables_sum[paste(BACFUN,"_",LEVEL,sep = "")] <- list(RA_table_current_sum)
    Taxonomy_RA_plot_tables[paste(BACFUN,"_",LEVEL,sep = "")] <- list(RA_plot_table_current)
    Taxonomy_RA_plots_gabriele[paste(BACFUN,"_",LEVEL,sep = "")] <- list(plot_gabriele_current)
  }
}

#How high can the rare taxa be in other samples?
rare_RAs <- data.frame()
for (BACFUN in c("bact","fung")) {
  physeq_current <- tax_glom(get(paste("phyloseq_",BACFUN,"_raref_RA",sep = "")), taxrank = "Genus",NArm = FALSE)
  
  #Select only the rare genera
  physeq_rare <- filter_taxa(physeq_current, function(x) mean(x) < 0.01, TRUE)
  
  #Table of abundances
  rare_breaks <- seq(0,0.3,by=0.05)
  rare_ranges <- paste(head(rare_breaks,-1), rare_breaks[-1], sep=" - ")
  rare_freq <- hist(as.vector(physeq_rare@otu_table), breaks=rare_breaks, include.lowest=TRUE, plot=FALSE)
  
  rare_RA_current <- data.frame(RA = rare_ranges, frequency = rare_freq$counts)
  colnames(rare_RA_current)[2] <- ifelse(BACFUN=="bact","Bacteria","Fungi")
  
  if (nrow(rare_RAs)==0) {
    rare_RAs <- rare_RA_current
  } else {
    rare_RAs <- merge.data.frame(rare_RAs,rare_RA_current, by = "RA")
  }
}

# #Export a figure for publication
# Taxonomy_figure <- plot_grid(Taxonomy_RA_plots$bact_Phylum + labs(tag = "A"),
#                       Taxonomy_RA_plots$fung_Phylum + labs(tag = "B"),
#                       ncol = 1)
# ggsave("Taxonomy_figure.svg",Taxonomy_figure,width = 2*230, height = 2*114, units = "mm")


# #Export a figure for publication
# Taxonomy_legend <- get_legend(Taxonomy_RA_plots_gabriele$bact_Phylum + guides(colour = guide_legend(nrow = 1)))
# Taxonomy_figure_bact <- plot_grid(
#   plot_grid(Taxonomy_RA_plots_gabriele$bact_Phylum + labs(tag = "A") + theme(legend.position = "none"),
#             Taxonomy_RA_plots_gabriele$bact_Family + labs(tag = "B") + theme(legend.position = "none"),
#             ncol = 2, rel_widths = c(1,1.5)),
#   Taxonomy_legend,
#   ncol = 1, rel_heights = c(1,0.1))
# Taxonomy_figure_fung <- plot_grid(
#   plot_grid(Taxonomy_RA_plots_gabriele$fung_Phylum + labs(tag = "A") + theme(legend.position = "none"),
#             Taxonomy_RA_plots_gabriele$fung_Family + labs(tag = "B") + theme(legend.position = "none"),
#             ncol = 2, rel_widths = c(1,1.5)),
#   Taxonomy_legend,
#   ncol = 1, rel_heights = c(1,0.1))
# ggsave("Taxonomy_figure_bact.svg",Taxonomy_figure_bact,width = 2*168, height = 2*86, units = "mm")
# ggsave("Taxonomy_figure_fung.svg",Taxonomy_figure_fung,width = 2*168, height = 2*86, units = "mm")
# ggsave("Taxonomy_figure_bact.jpg",Taxonomy_figure_bact,width = 2*168, height = 2*86, units = "mm")
# ggsave("Taxonomy_figure_fung.jpg",Taxonomy_figure_fung,width = 2*168, height = 2*86, units = "mm")

###Differential abundance###

#Split by DNA and RNA
phyloseq_bact_DNA <- subset_samples(phyloseq_bact, Community=="Total")
phyloseq_bact_RNA <- subset_samples(phyloseq_bact, Community=="Active")
phyloseq_fung_DNA <- subset_samples(phyloseq_fung, Community=="Total")
phyloseq_fung_RNA <- subset_samples(phyloseq_fung, Community=="Active")

#Rarefy to RA and do Linda with proportions
set.seed(123456)
phyloseq_bact_DNA_raref_RA <- phyloseq_mult_raref_avg(phyloseq_bact_DNA,min(colSums(phyloseq_bact_DNA@otu_table)),iter = 1000)
set.seed(123456)
phyloseq_bact_RNA_raref_RA <- phyloseq_mult_raref_avg(phyloseq_bact_RNA,min(colSums(phyloseq_bact_RNA@otu_table)),iter = 1000)
set.seed(123456)
phyloseq_fung_DNA_raref_RA <- phyloseq_mult_raref_avg(phyloseq_fung_DNA,min(colSums(phyloseq_fung_DNA@otu_table)),iter = 1000)
set.seed(123456)
phyloseq_fung_RNA_raref_RA <- phyloseq_mult_raref_avg(phyloseq_fung_RNA,min(colSums(phyloseq_fung_RNA@otu_table)),iter = 1000)

#Perform differential abundance analysis with Linda
linda_outs <- list()
linda_volcanos <- data.frame()
linda_plots <- list()
for (PHYSEQ in c("phyloseq_bact_DNA_raref_RA","phyloseq_bact_RNA_raref_RA","phyloseq_fung_DNA_raref_RA","phyloseq_fung_RNA_raref_RA")) {
  for (LEVEL in c("Phylum", "Class", "Order", "Family", "Genus")) {
    
    #Perform analysis at genus level
    physeq_current <- tax_glom(get(PHYSEQ), taxrank = LEVEL,NArm = FALSE)
    
    #Set factors with a reference for Linda
    physeq_current@sam_data$Treatment <- factor(physeq_current@sam_data$Treatment, levels = unique(physeq_current@sam_data$Treatment), ordered = F)
    physeq_current@sam_data$Treatment <- relevel(physeq_current@sam_data$Treatment,ref = "C")
    physeq_current@sam_data$Month <- factor(physeq_current@sam_data$Month, levels = unique(physeq_current@sam_data$Month), ordered = F)
    physeq_current@sam_data$Month <- relevel(physeq_current@sam_data$Month,ref = "July")
    
    # set seed
    set.seed(123456)
    # run linda
    # Using prev.filter that removes anything with less than 3 non-zero values. This is because Linda warns they have no statistical power.
    # linda_current <- linda(as(otu_table(physeq_current), "matrix"), as(sample_data(physeq_current), "data.frame"), 
    #                        formula = '~Treatment+Month+(1|Site)', alpha = 0.05, prev.filter=3/34, feature.dat.type="count", 
    #                        zero.handling="pseudo-count", n.cores=6)
    linda_current <- linda(as(otu_table(physeq_current), "matrix"), as(sample_data(physeq_current), "data.frame"), 
                           formula = '~Treatment+Month+(1|Site)', alpha = 0.05, prev.filter=3/34, feature.dat.type="proportion", 
                           n.cores=6)
    
    #Save result
    linda_outs[paste(PHYSEQ,LEVEL,sep = "_")] <- list(linda_current)
    
    # get results from all the tables in linda_out
    for (which_tab in names(linda_current$output)) {
      
      # get table
      volcano_data <- linda_current$output[[which_tab]]
      
      #Add taxon names and other useful information
      volcano_data[,c("Phylum", "Class", "Order", "Family", "Genus")] <- physeq_current@tax_table[rownames(volcano_data),c("Phylum", "Class", "Order", "Family", "Genus")]
      volcano_data$Group <- ifelse(grepl("bact", PHYSEQ),"Bacteria","Fungi")
      volcano_data$Community <- ifelse(grepl("DNA", PHYSEQ),"Total","Active")
      if (grepl("Treatment", which_tab)) { volcano_data$Comparison <- "C vs SW"}
      if (grepl("June", which_tab)) { volcano_data$Comparison <- "July vs June"}
      if (grepl("August", which_tab)) { volcano_data$Comparison <- "July vs August"}
      volcano_data$ASV <- row.names(volcano_data)
      volcano_data$Level <- LEVEL
      # apply threshold
      volcano_data <- volcano_data %>% mutate(threshold = ifelse(log2FoldChange >= 1 & padj <= 0.05,"up", ifelse(log2FoldChange <= -1  & padj <= 0.05, "down", "other")))
      volcano_data <- volcano_data %>% mutate(dsize=ifelse(log2FoldChange >= 1 & padj <= 0.05, 4, ifelse(log2FoldChange <= -1  & padj <= 0.05, 4, 2.5)))
      # compute log of padj
      volcano_data$log_p <- abs(log10(volcano_data$padj))
      #Add in RA of the base and comparison
      for (ASV in row.names(volcano_data)) {
        if (grepl("Treatment", which_tab)) { 
          volcano_data[ASV,"Base"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Treatment=="C"])
          volcano_data[ASV,"Comp"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Treatment=="SW"])
        }
        if (grepl("June", which_tab)) { 
          volcano_data[ASV,"Base"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Month=="July"])
          volcano_data[ASV,"Comp"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Month=="June"])
        }
        if (grepl("August", which_tab)) { 
          volcano_data[ASV,"Base"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Month=="July"])
          volcano_data[ASV,"Comp"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Month=="August"])
        }
      }
      
      #Rearrange columns
      volcano_data <- volcano_data[,c("Group","Community","Level","Comparison","threshold","log2FoldChange","pvalue","padj","Phylum","Class","Order","Family","Genus",
                                      "Base","Comp","baseMean","lfcSE","stat","reject","df","dsize","log_p","ASV")]
      
      #Save Volcano information
      linda_volcanos <- rbind(linda_volcanos,volcano_data)
      
      #Plot using p adjusted
      #Unlabel non significant dots, if any exist
      volcano_data[which(volcano_data$threshold=="other"),LEVEL] <- ""
      volcano_data[grep("Unidentified", volcano_data[,LEVEL]),LEVEL] <- ""
      
      # plot
      volcano <- ggplot(volcano_data, aes(x=log2FoldChange, y=log_p, label=!! sym(LEVEL))) +
        geom_point(aes(fill=threshold), shape=21, colour="black", size=volcano_data$dsize) +
        scale_fill_manual(values = c("up"="firebrick1", "down"="steelblue1",  "other"= "gray50")) +
        geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "gray30", linewidth=1) +
        geom_vline(xintercept=1, linetype="dashed", color = "gray30", linewidth=1) +
        geom_vline(xintercept=-1, linetype="dashed", color = "gray30", linewidth=1) +
        geom_text_repel(nudge_x=-0.2, nudge_y=0.05, max.overlaps=25) +
        ylab("-log10(p-value)") +
        xlab("log2FoldChange") +
        labs(title = paste(volcano_data$Group[1],": ",volcano_data$Type[1]," ",volcano_data$Comparison[1]," (",LEVEL,")",sep = "")) +
        theme_bw() +
        theme(text=element_text(size=16), legend.position = "none", plot.title = element_text(hjust = 0.5, size=16),axis.text = element_text(size=16))
      
      # save and export
      linda_plots[paste(volcano_data$Group[1],volcano_data$Type[1],gsub(" ","_",volcano_data$Comparison[1]),LEVEL,sep = "_")] <- list(volcano)
      # ggsave(paste("DA_",volcano_data$Group[1],"_",volcano_data$Type[1],"_",gsub(" ","_",volcano_data$Comparison[1]),".svg",sep = ""),volcano,width = 2*84, height = 2*84, units = "mm")
    }
  }
}

#Also test differential abundance between DNA and RNA samples

#Rarefy to RA and do Linda with proportions
set.seed(123456)
phyloseq_bact_raref_RA <- phyloseq_mult_raref_avg(phyloseq_bact,min(colSums(phyloseq_bact@otu_table)),iter = 1000)
set.seed(123456)
phyloseq_fung_raref_RA <- phyloseq_mult_raref_avg(phyloseq_fung,min(colSums(phyloseq_fung@otu_table)),iter = 1000)

#Perform differential abundance analysis with Linda
for (PHYSEQ in c("phyloseq_bact_raref_RA","phyloseq_fung_raref_RA")) {
  for (LEVEL in c("Phylum", "Class", "Order", "Family", "Genus")) {
    
    #Perform analysis at genus level
    physeq_current <- tax_glom(get(PHYSEQ), taxrank = LEVEL,NArm = FALSE)
    
    #Set factors with a reference for Linda
    physeq_current@sam_data$Community <- factor(physeq_current@sam_data$Community, levels = unique(physeq_current@sam_data$Community), ordered = F)
    physeq_current@sam_data$Community <- relevel(physeq_current@sam_data$Community,ref = "Total")
    
    # set seed
    set.seed(123456)
    # run linda
    # Using prev.filter that removes anything with less than 3 non-zero values. This is because Linda warns they have no statistical power.
    linda_current <- linda(as(otu_table(physeq_current), "matrix"), as(sample_data(physeq_current), "data.frame"), 
                           formula = '~Community+(1|Site)', alpha = 0.05, prev.filter=3/68, feature.dat.type="proportion", 
                           n.cores=6)
    
    #Save result
    linda_outs[paste(PHYSEQ,LEVEL,sep = "_")] <- list(linda_current)
    
    # get table
    volcano_data <- linda_current$output[["CommunityActive"]]
    
    #Add taxon names and other useful information
    volcano_data[,c("Phylum", "Class", "Order", "Family", "Genus")] <- physeq_current@tax_table[rownames(volcano_data),c("Phylum", "Class", "Order", "Family", "Genus")]
    volcano_data$Group <- ifelse(grepl("bact", PHYSEQ),"Bacteria","Fungi")
    volcano_data$Comparison <- "Total vs Active"
    volcano_data$Community <- "-"
    volcano_data$Level <- LEVEL
    volcano_data$ASV <- row.names(volcano_data)
    # apply threshold
    volcano_data <- volcano_data %>% mutate(threshold = ifelse(log2FoldChange >= 1 & padj <= 0.05,"up", ifelse(log2FoldChange <= -1  & padj <= 0.05, "down", "other")))
    volcano_data <- volcano_data %>% mutate(dsize=ifelse(log2FoldChange >= 1 & padj <= 0.05, 4, ifelse(log2FoldChange <= -1  & padj <= 0.05, 4, 2.5)))
    # compute log of padj
    volcano_data$log_p <- abs(log10(volcano_data$padj))
    #Add in RA of the base and comparison
    for (ASV in row.names(volcano_data)) {
      volcano_data[ASV,"Base"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Community=="Total"])
      volcano_data[ASV,"Comp"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Community=="Active"])
    }
    
    #Rearrange columns
    volcano_data <- volcano_data[,c("Group","Community","Level","Comparison","threshold","log2FoldChange","pvalue","padj","Phylum","Class","Order","Family","Genus",
                                    "Base","Comp","baseMean","lfcSE","stat","reject","df","dsize","log_p","ASV")]
    
    #Save Volcano information
    linda_volcanos <- rbind(linda_volcanos,volcano_data)
    
    #Plot using p adjusted
    #Unlabel non significant dots, if any exist
    volcano_data[which(volcano_data$threshold=="other"),LEVEL] <- ""
    volcano_data[grep("Unidentified", volcano_data[,LEVEL]),LEVEL] <- ""
    
    # plot
    volcano <- ggplot(volcano_data, aes(x=log2FoldChange, y=log_p, label=!! sym(LEVEL))) +
      geom_point(aes(fill=threshold), shape=21, colour="black", size=volcano_data$dsize) +
      scale_fill_manual(values = c("up"="firebrick1", "down"="steelblue1",  "other"= "gray50")) +
      geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "gray30", linewidth=1) +
      geom_vline(xintercept=1, linetype="dashed", color = "gray30", linewidth=1) +
      geom_vline(xintercept=-1, linetype="dashed", color = "gray30", linewidth=1) +
      geom_text_repel(nudge_x=-0.2, nudge_y=0.05, max.overlaps=25) +
      ylab("-log10(p-value)") +
      xlab("log2FoldChange") +
      labs(title = paste(volcano_data$Group[1],": ",volcano_data$Comparison[1]," (",LEVEL,")",sep = "")) +
      theme_bw() +
      theme(text=element_text(size=16), legend.position = "none", plot.title = element_text(hjust = 0.5, size=16),axis.text = element_text(size=16))
    
    # save and export
    linda_plots[paste(volcano_data$Group[1],gsub(" ","_",volcano_data$Comparison[1]),LEVEL,sep = "_")] <- list(volcano)
    
  }
}

#Save a table of the significantly differentially abundant taxa
linda_volcanos_significant <- linda_volcanos[linda_volcanos$padj<=0.05 & abs(linda_volcanos$log2FoldChange) >= 1,]
# write.table(linda_volcanos_significant,"linda_volcanos_significant.csv",quote = T,row.names = F, sep = ",")

# #Plot of the differentially abundant taxa for Control vs Treatment
# DA_combined <- plot_grid(linda_plots$Bacteria_DNA_Control_vs_Treatment + labs(tag = "A"),
#                          linda_plots$Bacteria_RNA_Control_vs_Treatment + labs(tag = "B"),
#                          linda_plots$Fungi_DNA_Control_vs_Treatment + labs(tag = "C"),
#                          linda_plots$Fungi_RNA_Control_vs_Treatment + labs(tag = "D"),
#                          ncol = 2)
# ggsave("DA_combined.svg",DA_combined,width = 2*114, height = 2*114, units = "mm")
# ggsave("DA_combined.jpg",DA_combined,width = 2*114*1.15, height = 2*114*1.15, units = "mm")

###Redundancy Analysis (RDA)###

#Need to fill in the missing metadata
#Removing June samples due to an excess of missing metadata in the treatments.
phyloseq_bact_Total_imputed <- impute_missing_metadata(subset_samples(phyloseq_bact_DNA, !Month=="June"))
phyloseq_bact_Active_imputed <- impute_missing_metadata(subset_samples(phyloseq_bact_RNA, !Month=="June"))
phyloseq_fung_Total_imputed <- impute_missing_metadata(subset_samples(phyloseq_fung_DNA, !Month=="June"))
phyloseq_fung_Active_imputed <- impute_missing_metadata(subset_samples(phyloseq_fung_RNA, !Month=="June"))

#Get average relative abundance over 1000 rarefactions
set.seed(123456)
phyloseq_bact_Total_imputed_raref_RA <- phyloseq_mult_raref_avg(phyloseq_bact_Total_imputed,min(colSums(phyloseq_bact_Total_imputed@otu_table)),iter = 1000)
set.seed(123456)
phyloseq_bact_Active_imputed_raref_RA <- phyloseq_mult_raref_avg(phyloseq_bact_Active_imputed,min(colSums(phyloseq_bact_Active_imputed@otu_table)),iter = 1000)
set.seed(123456)
phyloseq_fung_Total_imputed_raref_RA <- phyloseq_mult_raref_avg(phyloseq_fung_Total_imputed,min(colSums(phyloseq_fung_Total_imputed@otu_table)),iter = 1000)
set.seed(123456)
phyloseq_fung_Active_imputed_raref_RA <- phyloseq_mult_raref_avg(phyloseq_fung_Active_imputed,min(colSums(phyloseq_fung_Active_imputed@otu_table)),iter = 1000)

#Loop through and RDA for both genus and phylum level using DNA and RNA. Now only running at ASV/OTU level
# define lists to store the results from the ordinations
RDA_r_square_partial <- list()
RDA_signif_test_model <- list()
RDA_signif_test_terms <- list()
RDA_Figures <- list()
RDA_arrows <- data.frame()
for (BACFUN in c("bact","fung")) {
  for (TYPE in c("Total","Active")) {
    #for (LEVEL in c("Phylum","Genus")) {
      if (grepl("bact",BACFUN)) {taxonomic_group <- "Bacteria" }
      if (grepl("fung",BACFUN)) {taxonomic_group <- "Fungi" }
      
      #physeq_current <- tax_glom(get(paste("phyloseq_",BACFUN,"_",TYPE,"_imputed_raref_RA",sep = "")), taxrank = LEVEL,NArm = FALSE)
      physeq_current <- get(paste("phyloseq_",BACFUN,"_",TYPE,"_imputed_raref_RA",sep = ""))
      
      #Remove any taxa where there is no information
      physeq_current <- prune_taxa(apply(tax_table(physeq_current), 1, function(x) length(which(x==0))!=length(x)), physeq_current)
      
      #Perform the RDA with the set factors
      set.seed(123456)
      RDA_current <- ordinate(physeq_current, "RDA", formula=as.formula(".~Treatment+Month+Condition(Site)"), scale=1)
      
      # get adjusted R2
      r_square_partial_current <- RsquareAdj(RDA_current)$adj.r.squared
      
      #Test model significance
      set.seed(123456)
      RDA_signif_test_model_current <- anova.cca(RDA_current, step=1000)
      
      # get p-value
      set.seed(123456)
      RDA_signif_test_terms_current <- anova.cca(RDA_current, step=1000, by="terms")
      
      # create data frame for ggplot plotting
      # get scores for sites
      RDA_current_df <- as.data.frame(scores(RDA_current, display="sites", scaling = 1)) #Altered to use the same scores settings as plot_ordination
      RDA_current_df$label <- rownames(RDA_current_df)
      RDA_current_df[grep("SW", RDA_current_df$label),"Treatment"] <- "SW"
      RDA_current_df[grep("C", RDA_current_df$label),"Treatment"] <- "C"
      RDA_current_df[grep("June", RDA_current_df$label),"Month"] <- "June"
      RDA_current_df[grep("July", RDA_current_df$label),"Month"] <- "July"
      RDA_current_df[grep("August", RDA_current_df$label),"Month"] <- "August"
      
      # set factors
      RDA_current_df$Treatment <- factor(RDA_current_df$Treatment, levels=c("C", "SW"))
      RDA_current_df$Month <- factor(RDA_current_df$Month, levels=c("June", "July", "August"))
      
      #Variance explained by axes
      RDA_var_current <- as.data.frame(summary(RDA_current)$cont$importance)

      # compute fitting
      set.seed(123456)
      var_fitting_enviro_current <- envfit(RDA_current, sample_data(physeq_current)[, colnames(physeq_current@sam_data)[5:34]], perm=1000, na.rm=TRUE, scaling=1)
      # get multiplying factors for the arrows
      #Needs a plot to get the right coordinates
      plot(RDA_current, type="none", scaling=1, main="", display ="sites")
      arrow_factor <- ordiArrowMul(var_fitting_enviro_current, choices = c(1,2))
      # get arrows and scale them
      arrows <- as.data.frame(scores(var_fitting_enviro_current, display="vectors")) * arrow_factor #Scaling means nothing here
      arrows$envvars <- rownames(arrows)
      arrows$r2 <- var_fitting_enviro_current$vectors$r
      arrows$Taxon <- BACFUN
      arrows$Community <- TYPE
      #arrows$Level <- LEVEL
      arrows$p <- var_fitting_enviro_current$vectors$pvals
      # select arrows which have p-value <= 0.05
      sig_arrows <- arrows[names(which(var_fitting_enviro_current$vectors$pvals<=0.05)), ]
      
      #Fix the labels
      sig_arrows$envvars <- soil_chem_annot[soil_chem_annot$Variable%in%sig_arrows$envvars,"Name"]
      
      # define shapes for Month. both shapes must be in the same order
      timepoint_shapes <- c(21, 23, 24)
      names(timepoint_shapes) <- unique(RDA_current_df$Month)
      
      # ggplotting
      RDA_current_plot <- ggplot(RDA_current_df) + #Now matches
        geom_point(mapping=aes(x=RDA1, y=RDA2, fill=Treatment, shape=Month), size = 4) +
        coord_fixed() +
        scale_shape_manual(values=timepoint_shapes) +
        scale_fill_manual(values=treatment_colours) +
        geom_segment(data=sig_arrows, aes(x=0, xend=RDA1, y=0, yend=RDA2), arrow=arrow(length=unit(0.25, "cm")), colour="black", lwd=0.5) +
        geom_text_repel(data=sig_arrows, aes(x=RDA1, y=RDA2, label=envvars), col="black", alpha = 70,max.overlaps = 13)+
        guides(fill=guide_legend(override.aes=list(shape=22, size=4)), shape=guide_legend(override.aes=list(size=4))) +
        #labs(title = paste(taxonomic_group,": ",TYPE," - ",LEVEL,sep = ""),
        labs(title = paste(taxonomic_group,": ",TYPE,sep = ""),
             x = paste("RDA1 (",round(RDA_var_current[2,"RDA1"]*100,digits = 2),"%)",sep = ""), 
             y = paste("RDA2 (",round(RDA_var_current[2,"RDA2"]*100,digits = 2),"%)",sep = "")) +
        theme_bw() +
        theme(text=element_text(size=16), aspect.ratio = 1, axis.text = element_text(size=16), plot.title = element_text(hjust = 0.5, size=16))

      # # export ggplot
      # ggsave(paste("RDA_",BACFUN,"_",TYPE,"_",LEVEL,".eps",sep = ""),RDA_current_plot,width = 2*84, height = 2*84, units = "mm")
      
      #Save results to lists
      # RDA_r_square_partial[paste(BACFUN,"_",TYPE,"_",LEVEL,sep = "")] <- list(r_square_partial_current)
      # RDA_signif_test_terms[paste(BACFUN,"_",TYPE,"_",LEVEL,sep = "")] <- list(RDA_signif_test_terms_current)
      # RDA_Figures[paste(BACFUN,"_",TYPE,"_",LEVEL,sep = "")] <- list(RDA_current_plot)
      RDA_r_square_partial[paste(BACFUN,"_",TYPE,sep = "")] <- list(r_square_partial_current)
      RDA_signif_test_model[paste(BACFUN,"_",TYPE,sep = "")] <- list(RDA_signif_test_model_current)
      RDA_signif_test_terms[paste(BACFUN,"_",TYPE,sep = "")] <- list(RDA_signif_test_terms_current)
      RDA_Figures[paste(BACFUN,"_",TYPE,sep = "")] <- list(RDA_current_plot)
      RDA_arrows <- rbind(RDA_arrows, arrows)
      
   # }
  }
}

# #Plots for report
# report_RDA_Phylum <- plot_grid(RDA_Figures[["bact_DNA_Phylum"]],
#                                RDA_Figures[["bact_RNA_Phylum"]],
#                                RDA_Figures[["fung_DNA_Phylum"]],
#                                RDA_Figures[["fung_RNA_Phylum"]],
#                              ncol = 2)
# report_RDA_Genus <- plot_grid(RDA_Figures[["bact_DNA_Genus"]],
#                               RDA_Figures[["bact_RNA_Genus"]],
#                               RDA_Figures[["fung_DNA_Genus"]],
#                               RDA_Figures[["fung_RNA_Genus"]],
#                               ncol = 2)
# ggsave("report_RDA_Phylum.svg",report_RDA_Phylum,width = 2*168, height = 2*168, units = "mm")
# ggsave("report_RDA_Phylum.jpg",report_RDA_Phylum,width = 2*168, height = 2*168, units = "mm")
# ggsave("report_RDA_Genus.svg",report_RDA_Genus,width = 2*168, height = 2*168, units = "mm")
# ggsave("report_RDA_Genus.jpg",report_RDA_Genus,width = 2*168, height = 2*168, units = "mm")

#Find the variables which are brought up in the RDA

#Make a table of the RDA r2 values
RDA_arrows_signif <- RDA_arrows[RDA_arrows$p < 0.05,]
#RDA_arrows_signif <- RDA_arrows_signif[,c("Taxon","Level","Type","p","r2","envvars")]
RDA_arrows_signif <- RDA_arrows_signif[,c("Taxon","Community","p","r2","envvars")]
RDA_arrows_signif[grep("bact", RDA_arrows_signif$Taxon),"Taxon"] <- "Bacteria"
RDA_arrows_signif[grep("fung", RDA_arrows_signif$Taxon),"Taxon"] <- "Fungi"
# #export table
# write.table(RDA_arrows_signif, file = "RDA_significant_values.csv",quote = T,row.names = F, col.names = T, sep = ",")

###Fungal Guilds###

#Agglomerate to genus level because that's the level of the database
phyloseq_fung_guilds <- tax_glom(phyloseq_fung, taxrank="Genus", NArm = FALSE)

#Insert Camelia's FungalTraits data into the phyloseq object using the most common lifestyle for a given genus.
Guilds <- as.data.frame(tax_table(phyloseq_fung_guilds))
Guilds[,"Lifestyle"] <- "Unidentified"
#Search by genus
#For genus, the lifestyle taken is the most common lifestyle for that genus.
for (GEN in unique(Guilds$Genus)) {
  if (GEN%in%unique(CameliaTraits$V6)) {
    Guilds[which(Guilds$Genus==GEN,"Lifestyle"),"Lifestyle"] <- names(sort(table((CameliaTraits[CameliaTraits$V6 == GEN,"V8"])),decreasing = T)[1])
  }
}
#Search by family if missing genus
#For family, the lifestyle is taken as the most common lifestyle for that family
for (GEN in unique(Guilds[is.na(Guilds$Genus),"Family"])) {
  if (GEN%in%unique(CameliaTraits$V5)) {
    Guilds[which(Guilds$Family==GEN & is.na(Guilds$Genus)),"Lifestyle"] <- names(sort(table((CameliaTraits[CameliaTraits$V5 == GEN,"V8"])),decreasing = T)[1])
  }
}
# #Search by order if missing family
# #For order, the lifestyle is taken as the most common lifestyle for that family
# for (GEN in unique(Guilds[is.na(Guilds$Family),"Order"])) {
#   if (GEN%in%unique(CameliaTraits$V4)) {
#     Guilds[which(Guilds$Order==GEN & is.na(Guilds$Family)),"Lifestyle"] <- names(sort(table((CameliaTraits[CameliaTraits$V4 == GEN,"V8"])),decreasing = T)[1])
#   }
# }
Guilds$Lifestyle <- str_to_sentence(Guilds$Lifestyle)

##Lichen aside
#pull out lichenised fungi
lichenised <- Guilds[Guilds$Lifestyle=="Lichenized",]
phyloseq_fung_DNA_raref_RA_lichen <- subset_taxa(phyloseq_fung_DNA_raref_RA, Genus%in%lichenised$Genus)
phyloseq_fung_DNA_raref_RA_lichen_table <- psmelt(phyloseq_fung_DNA_raref_RA_lichen)
plot_lichen <- ggplot(phyloseq_fung_DNA_raref_RA_lichen_table, aes(fill = Genus, x = Site, y = Abundance)) +
  geom_col(color="black") +
  labs(x= element_blank(), y = "Relative Abundance") +
  facet_wrap(~Treatment*Month,ncol = 12, scales = "free_x") +
  #scale_fill_manual(values = colours_guilds, labels = gsub("_"," ",names(colours_guilds)), guide = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(text=element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),legend.position = "bottom",axis.text = element_text(size=16)) +
  guides(fill = guide_legend(nrow = 2))
#lichen_labeller <- as_labeller(c(`C` = "Control", `T` = "Treatment"))
lichen_plot <- ggplot(soil_chem[soil_chem$Month=="July",], aes(x = Treatment, y = Lichen)) + 
  #geom_col(aes(fill=Site)) + 
  geom_boxplot(aes(fill=NULL)) + 
  geom_point(aes(colour=Site),position = position_dodge(width = 0.2)) +
  ylab(label = "Lichen cover (%)") + 
  #facet_wrap(~Treatment, labeller = lichen_labeller) +
  scale_colour_manual(values = site_colours) +
  scale_x_discrete(labels = c("Control", "Treatment")) +
  scale_y_continuous(limits = c(0,52)) +
  xlab(label = NULL) +
  theme_bw() +
  theme(text=element_text(size=16, family = "ArialMT"), legend.position = "right", plot.title = element_text(hjust = 0.5, size = 16), 
        panel.spacing=unit(0,"lines"),axis.text = element_text(size=16))
jpeg("Lichen_Abundance.jpg", width=2*62, height=2*62, units = "mm", res = 300)
lichen_plot
dev.off()
##

Guilds <- Guilds[,c("Kingdom","Lifestyle")]
phyloseq_fung_guilds <- phyloseq(otu_table(phyloseq_fung_guilds), tax_table(as.matrix(Guilds)), sample_data(phyloseq_fung_guilds))

#Rarefy for plotting
set.seed(123456)
phyloseq_fung_guilds_raref_RA <- phyloseq_mult_raref_avg(phyloseq_fung_guilds,min(colSums(phyloseq_fung_guilds@otu_table)),iter = 1000)

#Plot abundance (merge <1%)
phyloseq_fung_guilds_raref_RA_glom <- tax_glom(phyloseq_fung_guilds_raref_RA, taxrank = "Lifestyle",NArm = FALSE)
phyloseq_fung_guilds_raref_RA_glom_dominant <- merge_taxa(phyloseq_fung_guilds_raref_RA_glom,taxa_names(filter_taxa(phyloseq_fung_guilds_raref_RA_glom, function(x) mean(x) < 0.01, TRUE)))
phyloseq_fung_guilds_raref_RA_glom_dominant@tax_table[is.na(phyloseq_fung_guilds_raref_RA_glom_dominant@tax_table[,"Lifestyle"]),"Lifestyle"] <- "Rare Lifestyles"

phyloseq_fung_guilds_raref_RA_glom_dominant_plot <- psmelt(phyloseq_fung_guilds_raref_RA_glom_dominant)
phyloseq_fung_guilds_raref_RA_glom_dominant_plot$Lifestyle <- factor(phyloseq_fung_guilds_raref_RA_glom_dominant_plot$Lifestyle,levels = c("Ectomycorrhizal","Ericoid_mycorrhizal","Lichen_parasite","Litter_saprotroph","Soil_saprotroph","Wood_saprotroph","Rare Lifestyles","Unidentified"))

colours_guilds <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(phyloseq_fung_guilds_raref_RA_glom_dominant_plot$Lifestyle)))
names(colours_guilds) <- c("Ectomycorrhizal","Ericoid_mycorrhizal","Lichen_parasite","Litter_saprotroph","Soil_saprotroph","Wood_saprotroph","Rare Lifestyles","Unidentified")

plot_guilds <- ggplot(phyloseq_fung_guilds_raref_RA_glom_dominant_plot, aes(fill = Lifestyle, x = Site, y = Abundance)) +
  geom_bar(position="fill", stat="identity", color="black") +
  labs(x= element_blank(), y = "Relative Abundance") +
  facet_wrap(~Treatment*Community*Month,ncol = 12, scales = "free_x") +
  scale_fill_manual(values = colours_guilds, labels = gsub("_"," ",names(colours_guilds)), guide = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(text=element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),legend.position = "bottom",axis.text = element_text(size=16)) +
  guides(fill = guide_legend(nrow = 2))

#Create RA table
RA_table_guilds <- as.data.frame(phyloseq_fung_guilds_raref_RA_glom_dominant@otu_table)
RA_table_guilds$Lifestyle <- as.character(phyloseq_fung_guilds_raref_RA_glom_dominant@tax_table[row.names(RA_table_guilds),"Lifestyle"])
RA_table_guilds <- RA_table_guilds[,c(ncol(RA_table_guilds),1:(ncol(RA_table_guilds)-1))]

#Create a tidy table without merging taxa
RA_table_guilds_full <- as.data.frame(phyloseq_fung_guilds_raref_RA_glom@otu_table)
RA_table_guilds_full$Lifestyle <- as.character(phyloseq_fung_guilds_raref_RA_glom@tax_table[row.names(RA_table_guilds_full),"Lifestyle"])
RA_table_guilds_full <- RA_table_guilds_full[,c(ncol(RA_table_guilds_full),1:(ncol(RA_table_guilds_full)-1))]

Fungal_Guilds <- data.frame()
for (DESC in colnames(RA_table_guilds_full)[2:ncol(RA_table_guilds_full)]) {
  current <- as.data.frame(matrix(nrow = nrow(RA_table_guilds_full), ncol = 7, dimnames = list(NULL,c("Description","Lifestyle","Abundance","Treatment","Site","Month","Community"))))
  current$Description <- DESC
  current$Lifestyle <- RA_table_guilds_full$Lifestyle
  current$Abundance <- RA_table_guilds_full[,DESC]
  Fungal_Guilds <- rbind(Fungal_Guilds,current)
}
for (ROW in 1:nrow(Fungal_Guilds)) {
  Fungal_Guilds[ROW,4:7] <- as.vector(unlist(sample_data(phyloseq_fung_guilds_raref_RA_glom_dominant)[Fungal_Guilds[ROW,"Description"],c("Treatment","Site","Month","Community")]))
}
Fungal_Guilds$Treatment <- factor(Fungal_Guilds$Treatment, levels = c("C","SW"),ordered = F)
Fungal_Guilds$Month <- factor(Fungal_Guilds$Month, levels = c("June","July","August"),ordered = T)
Fungal_Guilds$Site <- factor(Fungal_Guilds$Site, levels = unique(Fungal_Guilds$Site),ordered = F)
Fungal_Guilds$Community <- factor(Fungal_Guilds$Community, levels = c("Total","Active"),ordered = F)

#Create RA table with summary data
Fungal_Guilds_sum <- Fungal_Guilds[,c("Description","Lifestyle","Abundance")]
Fungal_Guilds_sum$Description <- gsub("s[0-9]_","",Fungal_Guilds_sum$Description)
Fungal_Guilds_sum <- aggregate(Abundance~Lifestyle+Description,mean,data = Fungal_Guilds_sum)
for (DESC in Fungal_Guilds_sum$Description) {
  Fungal_Guilds_sum[,DESC] <- Fungal_Guilds_sum[Fungal_Guilds_sum$Description==DESC,"Abundance"]
}
Fungal_Guilds_sum <- Fungal_Guilds_sum[1:17,c(-2,-3)]
Fungal_Guilds_sum[,"Overall"] <- rowMeans(Fungal_Guilds_sum[,2:13])

#Perform differential abundance of guilds with Linda

#Split by DNA and RNA
phyloseq_fung_guilds_raref_RA_DNA <- subset_samples(phyloseq_fung_guilds_raref_RA, Community=="Total")
phyloseq_fung_guilds_raref_RA_RNA <- subset_samples(phyloseq_fung_guilds_raref_RA, Community=="Active")

linda_volcanos_guild <- data.frame()
linda_plots_guild <- list()
linda_outs_guild <- list()
for (PHYSEQ in c("phyloseq_fung_guilds_raref_RA_DNA","phyloseq_fung_guilds_raref_RA_RNA")) {

  #Perform analysis at Lifestyle level
  physeq_current <- tax_glom(get(PHYSEQ), taxrank = "Lifestyle",NArm = FALSE)

  #Set factors with a reference for Linda
  physeq_current@sam_data$Treatment <- factor(physeq_current@sam_data$Treatment, levels = unique(physeq_current@sam_data$Treatment), ordered = F)
  physeq_current@sam_data$Treatment <- relevel(physeq_current@sam_data$Treatment,ref = "C")
  physeq_current@sam_data$Month <- factor(physeq_current@sam_data$Month, levels = unique(physeq_current@sam_data$Month), ordered = F)
  physeq_current@sam_data$Month <- relevel(physeq_current@sam_data$Month,ref = "July")

  # set seed
  set.seed(123456)
  # run linda
  linda_current <- linda(as(otu_table(physeq_current), "matrix"), as(sample_data(physeq_current), "data.frame"),
                           formula = '~Treatment+Month+(1|Site)', alpha = 0.05, prev.filter=3/34, feature.dat.type="proportion",
                           n.cores=6)
  

  #Save result
  linda_outs_guild[PHYSEQ] <- list(linda_current)

  # get results from all the tables in linda_out_guild
  for (which_tab in names(linda_current$output)) {

    # get table
    volcano_data <- linda_current$output[[which_tab]]

    #Add taxon names and other useful information
    volcano_data[,"Lifestyle"] <- as.vector(physeq_current@tax_table[rownames(volcano_data),"Lifestyle"])
    volcano_data$Group <- "Fungal Guilds"
    volcano_data$Community <- ifelse(grepl("DNA", PHYSEQ),"Total","Active")
    if (grepl("Treatment", which_tab)) { volcano_data$Comparison <- "Control vs Treatment"}
    if (grepl("June", which_tab)) { volcano_data$Comparison <- "July vs June"}
    if (grepl("August", which_tab)) { volcano_data$Comparison <- "July vs August"}
    volcano_data$ASV <- row.names(volcano_data)
    # apply threshold
    volcano_data <- volcano_data %>% mutate(threshold = ifelse(log2FoldChange >= 1 & padj <= 0.05,"up", ifelse(log2FoldChange <= -1  & padj <= 0.05, "down", "other")))
    volcano_data <- volcano_data %>% mutate(dsize=ifelse(log2FoldChange >= 1 & padj <= 0.05, 4, ifelse(log2FoldChange <= -1  & padj <= 0.05, 4, 2.5)))
    # compute log of padj
    volcano_data$log_p <- abs(log10(volcano_data$padj))
    #Add in RA of the base and comparison
    for (ASV in row.names(volcano_data)) {
      if (grepl("Treatment", which_tab)) {
        volcano_data[ASV,"Base"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Treatment=="C"])
        volcano_data[ASV,"Comp"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Treatment=="SW"])
      }
      if (grepl("June", which_tab)) {
        volcano_data[ASV,"Base"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Month=="July"])
        volcano_data[ASV,"Comp"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Month=="June"])
      }
      if (grepl("August", which_tab)) {
        volcano_data[ASV,"Base"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Month=="July"])
        volcano_data[ASV,"Comp"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Month=="August"])
      }
    }

    #Rearrange columns
    volcano_data <- volcano_data[,c("Group","Community","Comparison","threshold","log2FoldChange","pvalue","padj","Lifestyle",
                                    "Base","Comp","baseMean","lfcSE","stat","reject","df","dsize","log_p","ASV")]

    #Save Volcano information
    linda_volcanos_guild <- rbind(linda_volcanos_guild,volcano_data)

    #Plot using p adjusted
    #Unlabel non significant dots, if any exist
    volcano_data$Lifestyle[which(volcano_data$threshold=="other")] <- ""

    # plot
    volcano <- ggplot(volcano_data, aes(x=log2FoldChange, y=log_p, label=Lifestyle)) +
      geom_point(aes(fill=threshold), shape=21, colour="black", size=volcano_data$dsize) +
      scale_fill_manual(values = c("up"="firebrick1", "down"="steelblue1",  "other"= "gray50")) +
      geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "gray30", linewidth=1) +
      geom_vline(xintercept=1, linetype="dashed", color = "gray30", linewidth=1) +
      geom_vline(xintercept=-1, linetype="dashed", color = "gray30", linewidth=1) +
      geom_text_repel(nudge_x=-0.2, nudge_y=0.05, max.overlaps=25) +
      ylab("-log10(p-value)") +
      xlab("log2FoldChange") +
      labs(title = paste(volcano_data$Group[1],": ",volcano_data$Comparison[1],sep = "")) +
      theme_bw() +
      theme(text=element_text(size=16), legend.position = "none", plot.title = element_text(hjust = 0.5, size=16),axis.text = element_text(size=16))

    # save and export
    linda_plots_guild[paste(volcano_data$Group[1],volcano_data$Community[1],gsub(" ","_",volcano_data$Comparison[1]),sep = "_")] <- list(volcano)
  }
}

# #Perform analysis at Lifestyle level
# physeq_current <- tax_glom(phyloseq_fung_guilds_raref_RA, taxrank = "Lifestyle",NArm = FALSE)
# 
# #Set factors with a reference for Linda
# physeq_current@sam_data$Treatment <- factor(physeq_current@sam_data$Treatment, levels = unique(physeq_current@sam_data$Treatment), ordered = F)
# physeq_current@sam_data$Treatment <- relevel(physeq_current@sam_data$Treatment,ref = "C")
# physeq_current@sam_data$Month <- factor(physeq_current@sam_data$Month, levels = unique(physeq_current@sam_data$Month), ordered = F)
# physeq_current@sam_data$Month <- relevel(physeq_current@sam_data$Month,ref = "July")
# 
# # set seed
# set.seed(123456)
# # run linda
# linda_current <- linda(as(otu_table(physeq_current), "matrix"), as(sample_data(physeq_current), "data.frame"),
#                        formula = '~Treatment+Month+(1|Site)', alpha = 0.05, prev.filter=3/68, feature.dat.type="proportion",
#                        n.cores=6)
# 
# #Save result
# linda_outs_guild["phyloseq_fung_guilds_raref_RA"] <- list(linda_current)
# 
# # get results from all the tables in linda_out_guild
# for (which_tab in names(linda_current$output)) {
#   
#   # get table
#   volcano_data <- linda_current$output[[which_tab]]
#   
#   #Add taxon names and other useful information
#   volcano_data[,"Lifestyle"] <- as.vector(physeq_current@tax_table[rownames(volcano_data),"Lifestyle"])
#   volcano_data$Group <- "Fungal Guilds"
#   if (grepl("Treatment", which_tab)) { volcano_data$Comparison <- "Control vs Treatment"}
#   if (grepl("June", which_tab)) { volcano_data$Comparison <- "July vs June"}
#   if (grepl("August", which_tab)) { volcano_data$Comparison <- "July vs August"}
#   volcano_data$ASV <- row.names(volcano_data)
#   # apply threshold
#   volcano_data <- volcano_data %>% mutate(threshold = ifelse(log2FoldChange >= 1 & padj <= 0.05,"up", ifelse(log2FoldChange <= -1  & padj <= 0.05, "down", "other")))
#   volcano_data <- volcano_data %>% mutate(dsize=ifelse(log2FoldChange >= 1 & padj <= 0.05, 4, ifelse(log2FoldChange <= -1  & padj <= 0.05, 4, 2.5)))
#   # compute log of padj
#   volcano_data$log_p <- abs(log10(volcano_data$padj))
#   #Add in RA of the base and comparison
#   for (ASV in row.names(volcano_data)) {
#     if (grepl("Treatment", which_tab)) {
#       volcano_data[ASV,"Base"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Treatment=="C"])
#       volcano_data[ASV,"Comp"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Treatment=="SW"])
#     }
#     if (grepl("June", which_tab)) {
#       volcano_data[ASV,"Base"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Month=="July"])
#       volcano_data[ASV,"Comp"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Month=="June"])
#     }
#     if (grepl("August", which_tab)) {
#       volcano_data[ASV,"Base"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Month=="July"])
#       volcano_data[ASV,"Comp"] <- mean(physeq_current@otu_table[ASV,physeq_current@sam_data$Month=="August"])
#     }
#   }
#   
#   #Rearrange columns
#   volcano_data <- volcano_data[,c("Group","Comparison","threshold","log2FoldChange","pvalue","padj","Lifestyle",
#                                   "Base","Comp","baseMean","lfcSE","stat","reject","df","dsize","log_p","ASV")]
#   
#   #Save Volcano information
#   linda_volcanos_guild <- rbind(linda_volcanos_guild,volcano_data)
#   
#   #Plot using p adjusted
#   #Unlabel non significant dots, if any exist
#   volcano_data$Lifestyle[which(volcano_data$threshold=="other")] <- ""
#   
#   # plot
#   volcano <- ggplot(volcano_data, aes(x=log2FoldChange, y=log_p, label=Lifestyle)) +
#     geom_point(aes(fill=threshold), shape=21, colour="black", size=volcano_data$dsize) +
#     scale_fill_manual(values = c("up"="firebrick1", "down"="steelblue1",  "other"= "gray50")) +
#     geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "gray30", linewidth=1) +
#     geom_vline(xintercept=1, linetype="dashed", color = "gray30", linewidth=1) +
#     geom_vline(xintercept=-1, linetype="dashed", color = "gray30", linewidth=1) +
#     geom_text_repel(nudge_x=-0.2, nudge_y=0.05, max.overlaps=25) +
#     ylab("-log10(p-value)") +
#     xlab("log2FoldChange") +
#     labs(title = paste(volcano_data$Group[1],": ",volcano_data$Comparison[1],sep = "")) +
#     theme_bw() +
#     theme(text=element_text(size=16), legend.position = "none", plot.title = element_text(hjust = 0.5, size=16),axis.text = element_text(size=16))
#   
#   # save and export
#   linda_plots_guild[paste(volcano_data$Group[1],volcano_data$Community[1],gsub(" ","_",volcano_data$Comparison[1]),sep = "_")] <- list(volcano)
# }

#Save a table of the significantly differentially abundant taxa
linda_volcanos_guild_significant <- linda_volcanos_guild[linda_volcanos_guild$padj<=0.05 & abs(linda_volcanos_guild$log2FoldChange) >= 1,]

#Only want to look at mychorrizal and saprotrophic
#Want to group the saprotrophs
saps <- aggregate(Fungal_Guilds[grep("saprotroph", Fungal_Guilds$Lifestyle),c("Description","Abundance")],.~Description,FUN = sum )
saps$Lifestyle <- "Saprotroph"
for (ROW in 1:nrow(saps)) {
  saps[ROW,4:7] <- as.vector(unlist(sample_data(phyloseq_fung_guilds_raref_RA_glom_dominant)[saps[ROW,"Description"],c("Treatment","Site","Month","Community")]))
}
colnames(saps)[4:7] <- c("Treatment","Site","Month","Community")
saps <- saps[,colnames(Fungal_Guilds)]
Fungal_Guilds2 <- Fungal_Guilds[-grep("saprotroph", Fungal_Guilds$Lifestyle),]
Fungal_Guilds2 <- rbind(Fungal_Guilds2, saps)

# #Abundance over time
# Lifestyle_colours <- friendly_pal("wong_eight",3)
# names(Lifestyle_colours) <- c("Saprotroph","Ectomycorrhizal","Ericoid_mycorrhizal")
# plot_guilds_time <-
#   ggplot(Fungal_Guilds2[Fungal_Guilds2$Lifestyle=="Saprotroph" | Fungal_Guilds2$Lifestyle=="Ectomycorrhizal" | Fungal_Guilds2$Lifestyle=="Ericoid_mycorrhizal" ,], aes(x = Month, y = Abundance, colour = Lifestyle)) +
#   stat_summary(aes(y = Abundance, color = Lifestyle, group = Lifestyle), fun=mean, geom="line") +
#   geom_point(size = 1) +
#   facet_wrap(~Treatment,ncol = 2, scales = "free_x") +
#   labs(x= element_blank(), y = "Relative Abundance") +
#   scale_color_manual(values = Lifestyle_colours) +
#   theme_bw() +
#   theme(text=element_text(size=16), legend.position = "bottom",axis.text = element_text(size=16)) +
#   guides(colour = guide_legend(nrow = 2))

#Statistics
# set.seed(123456)
# guilds_aov_ecto <- aov(Abundance ~ Site+Treatment*Month, data=Fungal_Guilds2[Fungal_Guilds2$Lifestyle=="Ectomycorrhizal" ,])
# guilds_aov_summary_ecto <- summary(guilds_aov_ecto)
set.seed(123456)
guilds_aov_sapro <- aov(Abundance ~ Site+Treatment*Month*Community, data=Fungal_Guilds2[Fungal_Guilds2$Lifestyle=="Saprotroph" ,])
guilds_aov_summary_sapro <- summary(guilds_aov_sapro)
guilds_aov_tukey_sapro <- TukeyHSD(guilds_aov_sapro,ordered = TRUE)
# set.seed(123456)
# guilds_aov_erico <- aov(Abundance ~ Site+Treatment*Month, data=Fungal_Guilds2[Fungal_Guilds2$Lifestyle=="Ericoid_mycorrhizal" ,])
# guilds_aov_summary_erico <- summary(guilds_aov_erico)
# guilds_aov_tukey_erico <- TukeyHSD(guilds_aov_erico,ordered = TRUE)

#Export results
# ggsave("Fungal_guilds.svg",plot_guilds,width = 2*168, height = 2*84, units = "mm")
# ggsave("Fungal_guilds_time.svg",plot_guilds_time,width = 2*80, height = 2*60, units = "mm")
# ggsave("Fungal_guilds.jpg",plot_guilds,width = 2*168, height = 2*84, units = "mm")
# ggsave("Fungal_guilds_time.jpg",plot_guilds_time,width = 2*80, height = 2*60, units = "mm")
# write.table(RA_table_guilds, "RA_table_guilds.csv",quote = T,row.names = F, col.names = T, sep = ",")
# write.table(lichenised, file = "lichenised_arctic_tundra.csv",quote = T,row.names = F, col.names = T, sep = ",")
# write.table(Fungal_Guilds_sum, "RA_table_guilds_sum.csv",quote = T,row.names = F, col.names = T, sep = ",")

###Publication Exports###

#Main figures
Fig_X_Environmental_Variables <- plot_grid(nbchem[[1]] + theme(legend.position = "none"), nbchem[[2]] + theme(legend.position = "none"),
                                           nbchem[[3]] + theme(legend.position = "none"), nbchem[[4]] + theme(legend.position = "none"),
                                           nbchem[[5]] + theme(legend.position = "none"), nbchem[[6]] + theme(legend.position = "none"),
                                           nbchem[[7]] + theme(legend.position = "none"), nbchem[[8]] + theme(legend.position = "none"),
                                           nbchem[[9]] + theme(legend.position = "none"), nbchem[[10]] + theme(legend.position = "none"),
                                           nbchem[[11]] + theme(legend.position = "none"), nbchem[[15]] + theme(legend.position = "none"),
                                           nbchem[[16]] + theme(legend.position = "none"), nbchem[[17]] + theme(legend.position = "none"),
                                           nbchem[[18]] + theme(legend.position = "none"), nbchem[[12]] + theme(legend.position = "none"),
                                           nbchem[[13]] + theme(legend.position = "none"), nbchem[[14]] + theme(legend.position = "none"),
                                           nbchem[[19]] + theme(legend.position = "none"), nbchem[[20]] + theme(legend.position = "none"),
                                           nbchem[[21]] + theme(legend.position = "none"), nbchem[[22]] + theme(legend.position = "none"),
                                           nbchem[[23]] + theme(legend.position = "none"), nbchem[[24]] + theme(legend.position = "none"),
                                           nbchem[[25]] + theme(legend.position = "none"), nbchem[[26]] + theme(legend.position = "none"),
                                           nbchem[[27]] + theme(legend.position = "none"), nbchem[[28]] + theme(legend.position = "none"),
                                           nbchem[[29]] + theme(legend.position = "none"), nbchem[[30]] + theme(legend.position = "none"),
                                           ncol = 5, align = "hv", axis = "tblr")
#SVG
svg("Fig_1.svg", width=2*210*0.03937008, height=2*297*0.03937008)
Fig_X_Environmental_Variables
dev.off()
#JPG
jpeg("Fig_1.jpg", width=2*210, height=2*297, units = "mm", res = 300)
Fig_X_Environmental_Variables
dev.off()


Site_legend <- get_legend(alpha_plots[["observed_alphadiv_bact_mean"]] + theme(legend.position = "bottom")+ guides(colour = guide_legend(nrow = 1)))
Beta_legend <- get_legend(NMDS_Hellinger_bact + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 1)))
Fig_X_Alpha_Beta_Diversity <- plot_grid(
  plot_grid(
    plot_grid(alpha_plots[["observed_alphadiv_bact_mean"]] + theme(legend.position = "none") + labs(tag="A"),
              alpha_plots[["observed_alphadiv_fung_mean"]] + theme(legend.position = "none") + labs(tag="B"),
              ncol = 2),
    Site_legend, ncol = 1, rel_heights = c(1,0.1)),
  plot_grid(
    plot_grid(NMDS_Hellinger_bact + theme(legend.position = "none") + labs(tag="C"),
              NMDS_Hellinger_fung + theme(legend.position = "none") + labs(tag="D"),
              ncol = 2),
    Beta_legend, ncol = 1, rel_heights = c(1,0.1)),
  ncol = 1)
#SVG
svg("Fig_2.svg", width=2*168*0.03937008, height=2*120*0.03937008)
Fig_X_Alpha_Beta_Diversity
dev.off()
#JPG
jpeg("Fig_2.jpg", width=2*168, height=2*120, units = "mm", res = 300)
Fig_X_Alpha_Beta_Diversity
dev.off()






PresentationTable <- linda_volcanos_significant[linda_volcanos_significant$Level=="Order" & linda_volcanos_significant$Comparison == "C vs SW",]
PresentationTable["BaseLabel"] <- "RA-C"
PresentationTable["CompLabel"] <- "RA-SW"
PresentationTable[!PresentationTable$Order=="Russulales","BaseLabel"] <- ""
PresentationTable[!PresentationTable$Order=="Russulales","CompLabel"] <- ""
PresentationTable$Order <- factor(PresentationTable$Order,levels = sort(unique(PresentationTable$Order),decreasing = T),ordered = T)
Presentation_DA <- ggplot(PresentationTable, aes(x = Order, y = log2FoldChange, fill = log2FoldChange>1)) +
  geom_col(colour="black") +
  scale_fill_manual(values = c(treatment_colours["C"][[1]],treatment_colours["SW"][[1]])) +
  new_scale_colour() +
  geom_text(aes(label=round(Base, digits = 3)), colour="#0072B2", y=3, hjust=0, size=20, size.unit = "pt", vjust=-0.5) +
  geom_text(aes(label=round(Comp, digits = 3)), colour="#D55E00", y=4, hjust=0, size=20, size.unit = "pt", vjust=-0.5) +
  geom_text(aes(label=BaseLabel, y=3, hjust=0), colour="black", size=20, size.unit = "pt", vjust=1.5) +
  geom_text(aes(label=CompLabel, y=4, hjust=0), colour="black", size=20, size.unit = "pt", vjust=1.5) +
  scale_y_continuous(breaks = c(-5,-4,-3,-2,-1,0,1,2), minor_breaks = c(-4.5,-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5) , limits = c(-5,5)) +
  facet_grid(rows=vars(Group, Community), scales = "free_y", switch = "y", space = "free_y") +
  labs(x = NULL, y = "Log2Fold Change") +
  theme_minimal() +
  theme(strip.placement = "outside", strip.background = element_rect(fill = "grey85", colour = "grey20")) +
  theme(panel.grid.major = element_line(color = "gray80", linewidth = 0.5)) +
  theme(text=element_text(size=20), axis.text.x = element_text(hjust=1), legend.position = "none", axis.text = element_text(size=20)) +
  coord_flip()
#SVG
svg("Fig_4.svg", width=2*168*0.03937008, height=2*86*0.03937008)
Presentation_DA
dev.off()
#JPG
jpeg("Fig_4.jpg", width=2*168, height=2*86, units = "mm", res = 300)
Presentation_DA
dev.off()
#For graphical abstract
svg("Graphical_Abstract_panel.svg", width=2*135*0.03937008, height=2*86*0.03937008)
Presentation_DA + scale_y_continuous(breaks = c(-5,-4,-3,-2,-1,0,1,2), minor_breaks = c(-4.5,-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5) , limits = c(-5,2.5))
dev.off()


RDA_legend <- get_legend(RDA_Figures$bact_Total + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 1)))
Fig_X_RDA <- plot_grid(
  plot_grid(RDA_Figures[["bact_Total"]] + theme(legend.position = "none") + labs(tag = "A"),
                              RDA_Figures[["bact_Active"]] + theme(legend.position = "none") + labs(tag = "B"),
                              RDA_Figures[["fung_Total"]] + theme(legend.position = "none") + labs(tag = "C"),
                              RDA_Figures[["fung_Active"]] + theme(legend.position = "none") + labs(tag = "D"),
                              ncol = 2),
  RDA_legend, ncol = 1, rel_heights = c(1,0.1))
#SVG
svg("Fig_3.svg", width=2*114*0.03937008, height=2*114*0.03937008)
Fig_X_RDA
dev.off()
#JPG
jpeg("Fig_3.jpg", width=2*114, height=2*114, units = "mm", res = 300)
Fig_X_RDA
dev.off()

Guild_DA_Table <- linda_volcanos_guild_significant
Guild_DA_Table$Lifestyle <- gsub("_"," ",Guild_DA_Table$Lifestyle)
Guild_DA_Table$Lifestyle <- factor(Guild_DA_Table$Lifestyle,levels = sort(unique(Guild_DA_Table$Lifestyle),decreasing = T),ordered = T)
Guild_DA_Table$Comparison <- ifelse(Guild_DA_Table$Comparison == "July vs June", "JlvJn", "JlvAu")
Guild_DA_Table$Comparison <- factor(Guild_DA_Table$Comparison,levels = c("JlvJn","July","JlvAu"), ordered = T)
Guild_DA_Table$MonthColour <- ifelse(Guild_DA_Table$log2FoldChange < 0, "July", ifelse(Guild_DA_Table$Comparison == "JlvJn", "June", "August"))
Guild_DA_Plot <- ggplot(Guild_DA_Table, aes(x = Lifestyle, y = log2FoldChange, fill = MonthColour)) +
  geom_col(colour="black", show.legend = T) +
  scale_fill_manual(name = "Month", 
                    values = timepoint_colours) +
  new_scale_colour() +
  geom_text(aes(label=round(Base, digits = 3)), colour=timepoint_colours["July"], y=-3.5, size=16, size.unit = "pt",vjust=-0.2) +
  geom_text(aes(label=round(Comp, digits = 3),colour=Comparison), y=4.5, size=16, size.unit = "pt",show.legend = F,vjust=-0.2) +
  scale_color_manual(values = list("JlvJn" = timepoint_colours["June"],"JlvAu" = timepoint_colours["August"]), guide="none") +
  scale_y_continuous(breaks = c(-3,-2,-2,-1,0,1,2,3,4), minor_breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5) , limits = c(-4,4.5)) +
  facet_grid(rows=vars(Community,Comparison), scales = "free_y", switch = "y", space = "free_y") +
  #facet_grid(rows=vars(Group, Community), scales = "free_y", switch = "y", space = "free_y") +
  labs(x = NULL, y = "Log2Fold Change") +
  theme_minimal() +
  theme(strip.placement = "outside", strip.background = element_rect(fill = "grey85", colour = "grey20")) +
  theme(panel.grid.major = element_line(color = "gray80", linewidth = 0.5)) +
  theme(text=element_text(size=20), axis.text.x = element_text(hjust=1), legend.position = "bottom", axis.text = element_text(size=20)) +
  coord_flip()

Fig_X_Guilds <- plot_grid(plot_guilds + labs(tag = "A"),
                          Guild_DA_Plot + labs(tag = "B"),
                          ncol = 1)
#SVG
svg("Fig_5.svg", width=2*168*0.03937008, height=2*150*0.03937008)
Fig_X_Guilds
dev.off()
#JPG
jpeg("Fig_5.jpg", width=2*168, height=2*150, units = "mm", res = 300)
Fig_X_Guilds
dev.off()

#Supplementary figures
Fig_SX_Soil_Temp <- plot_grid(soil_temp_full + theme(legend.position = "none") + labs(tag="A"),
                             soil_temp_zoom + labs(tag="B"),
                             ncol = 1, align = "v", axis = "lr")
#SVG
svg("Fig_S2.svg", width=2*168*0.03937008, height=2*84*0.03937008)
Fig_SX_Soil_Temp
dev.off()
#JPG
jpeg("Fig_S2.jpg", width=2*168, height=2*84, units = "mm", res = 300)
Fig_SX_Soil_Temp
dev.off()

#SVG
svg("Fig_S3.svg", width=2*168*0.03937008, height=2*168*0.03937008)
soil_temp_compare_monthly
dev.off()
#JPG
jpeg("Fig_S3.jpg", width=2*168, height=2*168, units = "mm", res = 300)
soil_temp_compare_monthly
dev.off()
# 
# Fig_SX_Aitch_ordination <- plot_grid(
#   plot_grid(NMDS_Aitchison_bact + theme(legend.position = "none") + labs(tag="A"),
#             NMDS_Aitchison_fung + theme(legend.position = "none") + labs(tag="B"),
#             ncol = 2),
#   Beta_legend, ncol = 1, rel_heights = c(1,0.1))
# #SVG
# svg("Fig_SX_Aitch_ordination.svg", width=2*168*0.03937008, height=2*60*0.03937008)
# Fig_SX_Aitch_ordination
# dev.off()
# #JPG
# jpeg("Fig_SX_Aitch_ordination.jpg", width=2*168, height=2*60, units = "mm", res = 300)
# Fig_SX_Aitch_ordination
# dev.off()

Beta_site_legend <- get_legend(NMDS_Hellinger_bact_site + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 1)))
Fig_SX_Hell_ordination_sites <- plot_grid(
  plot_grid(NMDS_Hellinger_bact_site + theme(legend.position = "none") + labs(tag="A"),
            NMDS_Hellinger_fung_site + theme(legend.position = "none") + labs(tag="B"),
            ncol = 2),
  Beta_site_legend, ncol = 1, rel_heights = c(1,0.1))
#SVG
svg("Fig_S6.svg", width=2*168*0.03937008, height=2*60*0.03937008)
Fig_SX_Hell_ordination_sites
dev.off()
#JPG
jpeg("Fig_S6.jpg", width=2*168, height=2*60, units = "mm", res = 300)
Fig_SX_Hell_ordination_sites
dev.off()

Site_legend_rare <- get_legend(observed_bact_rares + theme(legend.position = "bottom")+ guides(colour = guide_legend(nrow = 1)))
Fig_SX_Rich_Rare <-   plot_grid(
  plot_grid(observed_bact_rares + theme(legend.position = "none") + labs(tag="A"),
            observed_fung_rares + theme(legend.position = "none") + labs(tag="B"),
            plot_bact_depth_corr + theme(legend.position = "none") + labs(tag="C"),
            plot_fung_depth_corr + theme(legend.position = "none") + labs(tag="D"),
            ncol = 2),
  Site_legend_rare, ncol = 1, rel_heights = c(1,0.05))
#SVG
svg("Fig_S1.svg", width=2*168*0.03937008, height=2*168*0.03937008)
Fig_SX_Rich_Rare
dev.off()
#JPG
jpeg("Fig_S1.jpg", width=2*168, height=2*168, units = "mm", res = 300)
Fig_SX_Rich_Rare
dev.off()

Fig_SX_Shannon <- plot_grid(
  plot_grid(alpha_plots[["shannon_alphadiv_bact_mean"]] + theme(legend.position = "none") + labs(tag="A"),
            alpha_plots[["shannon_alphadiv_fung_mean"]] + theme(legend.position = "none") + labs(tag="B"),
            ncol = 2),
  Site_legend, ncol = 1, rel_heights = c(1,0.1))
#SVG
svg("Fig_S4.svg", width=2*168*0.03937008, height=2*60*0.03937008)
Fig_SX_Shannon
dev.off()
#JPG
jpeg("Fig_S4.jpg", width=2*168, height=2*60, units = "mm", res = 300)
Fig_SX_Shannon
dev.off()

Fig_SX_Simpson <- plot_grid(
  plot_grid(alpha_plots[["Simpson_alphadiv_bact_mean"]] + theme(legend.position = "none") + labs(tag="A"),
            alpha_plots[["Simpson_alphadiv_fung_mean"]] + theme(legend.position = "none") + labs(tag="B"),
            ncol = 2),
  Site_legend, ncol = 1, rel_heights = c(1,0.1))
#SVG
svg("Fig_S5.svg", width=2*168*0.03937008, height=2*60*0.03937008)
Fig_SX_Simpson
dev.off()
#JPG
jpeg("Fig_S5.jpg", width=2*168, height=2*60, units = "mm", res = 300)
Fig_SX_Simpson
dev.off()

Taxonomy_legend <- get_legend(Taxonomy_RA_plots_gabriele$bact_Phylum + guides(colour = guide_legend(nrow = 1)))
Fig_X_RA <- plot_grid(plot_grid(
  plot_grid(Taxonomy_RA_plots_gabriele$bact_Phylum + labs(tag = "A") + theme(legend.position = "none"),
            Taxonomy_RA_plots_gabriele$bact_Order + labs(tag = "B") + theme(legend.position = "none"),
            Taxonomy_RA_plots_gabriele$fung_Phylum + labs(tag = "C") + theme(legend.position = "none"),
            Taxonomy_RA_plots_gabriele$fung_Order + labs(tag = "D") + theme(legend.position = "none"),
            ncol = 2, rel_widths = c(1,1.5)),
  Taxonomy_legend, ncol = 1, rel_heights = c(1,0.05)),
  ncol = 1)
#SVG
svg("Fig_S7.svg", width=2*168*0.03937008, height=2*172*0.03937008)
Fig_X_RA
dev.off()
#JPG
jpeg("Fig_S7.jpg", width=2*168, height=2*172, units = "mm", res = 300)
Fig_X_RA
dev.off()

# Fig_SX_DA <- plot_grid(linda_plots$Bacteria_DNA_vs_RNA_Genus + labs(tag = "A"),
#                        linda_plots$Fungi_DNA_vs_RNA_Genus + labs(tag = "B"),
#                        linda_plots$Bacteria_DNA_July_vs_June_Genus + labs(tag = "C"),
#                       linda_plots$Bacteria_DNA_July_vs_August_Genus + labs(tag = "D"),
#                       linda_plots$Bacteria_RNA_July_vs_June_Genus + labs(tag = "E"),
#                       linda_plots$Bacteria_RNA_July_vs_August_Genus + labs(tag = "F"),
#                       linda_plots$Fungi_DNA_July_vs_June_Genus + labs(tag = "G"),
#                       linda_plots$Fungi_DNA_July_vs_August_Genus + labs(tag = "H"),
#                       linda_plots$Fungi_RNA_July_vs_June_Genus + labs(tag = "I"),
#                       linda_plots$Fungi_RNA_July_vs_August_Genus + labs(tag = "J"),
#                       ncol = 2)
# ggsave("Fig_SX_DA.svg",Fig_SX_DA,width = 2*114, height = 2*285, units = "mm")
# ggsave("Fig_SX_DA.jpg",Fig_SX_DA,width = 2*114*1.15, height = 2*285*1.15, units = "mm")

#Main tables

write.table(Soil_chem_table_signif, file = "Table 1.csv",quote = T,row.names = F, col.names = T, sep = ",")
write.table(beta_perm_table,"Table 2.csv",quote = T,row.names = F,col.names = F ,sep = ",")

# Table_X_DA <- linda_volcanos_significant[linda_volcanos_significant$Comparison == "Control vs Treatment" & linda_volcanos_significant$Level%in%c("Phylum","Order"),
#                                          c("Group","Community","Level","Comparison","threshold","log2FoldChange","padj","Phylum","Class","Order","Base","Comp")]
# write.table(Table_X_DA,"Table_X_DA.csv",quote = T,row.names = F ,sep = ",")

#Supplementary tables

write.table(alpha_anova_table, file = "Table S1.csv",quote = T,row.names = F, col.names = F, sep = ",")
write.table(linda_volcanos_significant,"Table S2.csv",quote = T,row.names = F, sep = ",")
write.table(linda_volcanos_guild_significant,"Table S3.csv",quote = T,row.names = F, sep = ",")

# sessionInfo()

# R version 4.5.0 (2025-04-11)
# Platform: x86_64-pc-linux-gnu
# Running under: Arch Linux
# 
# Matrix products: default
# BLAS/LAPACK: /usr/lib/libopenblas.so.0.3;  LAPACK version 3.12.0
# 
# locale:
# [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8   
# [6] LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Prague
# tzcode source: system (glibc)
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] metagMisc_0.5.0     ggpubfigs_0.0.1     mirlyn_1.4.2        phyloseq_1.52.0     stringr_1.5.1       ggh4x_0.3.0         MicrobiomeStat_1.2 
# [8] lubridate_1.9.4     DescTools_0.99.60   multcompView_0.1-10 nlme_3.1-168        vegan_2.6-10        lattice_0.22-6      permute_0.9-7      
# [15] ggtext_0.1.2        ggrepel_0.9.6       ggsignif_0.6.4      ggpubr_0.6.0        ellipse_0.5.0       RColorBrewer_1.1-3  reshape2_1.4.4     
# [22] dplyr_1.1.4         cowplot_1.1.3       plyr_1.8.9          scales_1.3.0        ggplot2_3.5.2      
# 
# loaded via a namespace (and not attached):
# [1] rstudioapi_0.17.1       jsonlite_2.0.0          magrittr_2.0.3          farver_2.1.2            nloptr_2.2.1            fs_1.6.6               
# [7] vctrs_0.6.5             multtest_2.64.0         minqa_1.2.8             rstatix_0.7.2           forcats_1.0.0           haven_2.5.4            
# [13] broom_1.0.8             cellranger_1.1.0        Rhdf5lib_1.30.0         Formula_1.2-5           rhdf5_2.52.0            rootSolve_1.8.2.4      
# [19] igraph_2.1.4            lifecycle_1.0.4         iterators_1.0.14        pkgconfig_2.0.3         Matrix_1.7-3            R6_2.6.1               
# [25] GenomeInfoDbData_1.2.14 rbibutils_2.3           clue_0.3-66             digest_0.6.37           Exact_3.3               numDeriv_2016.8-1.1    
# [31] colorspace_2.1-1        spatial_7.3-18          S4Vectors_0.46.0        labeling_0.4.3          timechange_0.3.0        httr_1.4.7             
# [37] abind_1.4-8             mgcv_1.9-1              compiler_4.5.0          proxy_0.4-27            withr_3.0.2             backports_1.5.0        
# [43] carData_3.0-5           MASS_7.3-65             biomformat_1.36.0       fBasics_4041.97         gld_2.6.7               tools_4.5.0            
# [49] ape_5.8-1               glue_1.8.0              stabledist_0.7-2        rhdf5filters_1.20.0     gridtext_0.1.5          grid_4.5.0             
# [55] cluster_2.1.8.1         ade4_1.7-23             generics_0.1.3          gtable_0.3.6            tzdb_0.5.0              class_7.3-23           
# [61] tidyr_1.3.1             data.table_1.17.0       lmom_3.2                hms_1.1.3               xml2_1.3.8              car_3.1-3              
# [67] XVector_0.48.0          rmutil_1.1.10           BiocGenerics_0.54.0     foreach_1.5.2           pillar_1.10.2           splines_4.5.0          
# [73] survival_3.8-3          tidyselect_1.2.1        Biostrings_2.76.0       reformulas_0.4.0        IRanges_2.42.0          stats4_4.5.0           
# [79] expm_1.0-0              Biobase_2.68.0          statmod_1.5.0           factoextra_1.0.7        timeDate_4041.110       matrixStats_1.5.0      
# [85] stringi_1.8.7           UCSC.utils_1.4.0        boot_1.3-31             codetools_0.2-20        timeSeries_4041.111     tibble_3.2.1           
# [91] cli_3.6.4               rpart_4.1.24            Rdpack_2.6.4            munsell_0.5.1           Rcpp_1.0.14             GenomeInfoDb_1.44.0    
# [97] readxl_1.4.5            stable_1.1.6            parallel_4.5.0          modeest_2.4.0           readr_2.1.5             lme4_1.1-37            
# [103] mvtnorm_1.3-3           lmerTest_3.1-3          e1071_1.7-16            statip_0.2.3            purrr_1.0.4             crayon_1.5.3           
# [109] rlang_1.1.6 