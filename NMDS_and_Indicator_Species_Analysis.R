#Jornada_DxD_Fungal_Metabarcoding

#Protocol for processing Metabarcoding of soil
#File name: NMDS_and_Indicator_Species_Analysis

#In this script, you will do an: 
# 1. Non-metric Multi-dimensional Scaling (NMDS) Analysis
# 2. Indicator Species Analysis by treatment and location.


#Clear environment
rm(list = ls(all.names = TRUE))

#Load packages
library(vegan) #For the NMDS
library(funrar) #For calculating relative abundance
library(ggplot2) #For plotting
library(indicspecies) #For indicator species analysis

#1. Non-metric Multi-dimensional Scaling (NMDS) Analysis

#Import the 'otu.tables.clean.csv' file. You can find it in this repository.  
otu.clean <- read.csv("/project/egcc/stats/csv/otu.tables.clean.csv", row.names = 1)

#Select only the columns with the ASV data.
com <- otu.clean[,5:ncol(otu.clean)]

#Make it into a matrix.
m_com <- as.matrix(com)

#Calculate relative abundances
rel_com <- make_relative(m_com)

#Run the NMDS
set.seed(100)
nmds <- metaMDS(rel_com, distance = "bray", trymax = 2000, k=3, autotransform = F, noshare = TRUE)

nmds

#Call:
#metaMDS(comm = rel_com, distance = "bray", k = 3, trymax = 2000, autotransform = F, noshare = TRUE) 

#global Multidimensional Scaling using monoMDS

#Data:     rel_com 
#Distance: bray shortest 

#Dimensions: 3 
#Stress:     0.1449432 
#Stress type 1, weak ties
#Two convergent solutions found after 29 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on 'rel_com' 

plot(nmds) #for quick visualization.

#Extract NMDS scores (x and y coordinates)
data.scores <- as.data.frame(scores(nmds))

#Add columns with 'Treatment' and 'Location' information to the 'data.scores' data frame 
data.scores$Treatment <- otu.clean$treatment
data.scores$Location <- otu.clean$location

head(data.scores)

#Create a plot using ggplot2
xx <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( shape = Location, colour = Treatment))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "MDS1", colour = "Treatment", y = "MDS2", shape = "Location")  + 
  scale_colour_manual(values=c("#FFC300","#CD7F32","#9E0142", "#F46D43")) 
xx

ggsave(width = 8, height = 5, "/project/egcc/stats/NMDS/NMDS.jpg")


# 2. Indicator Species Analysis by treatment and location.

#For this step, use the 'rel_com' matrix created above since need relative abundances.

head(rel_com)

#Store the values for the factors treatment and location.
treatment <- otu.clean$treatment
location <- otu.clean$location

#Run analysis by treatment.
inv <- multipatt(rel_com, treatment, func = "r.g", control = how(nperm = 9999))
summary(inv)

#Multilevel pattern analysis
#---------------------------
  
#  Association function: r.g
#Significance level (alpha): 0.05

#Total number of species: 705
#Selected number of species: 3 
#Number of species associated to 1 group: 3 
#Number of species associated to 2 groups: 0 
#Number of species associated to 3 groups: 0 

#List of species associated to each combination: 
  
#  Group control   #sps.  1 
#       stat p.value   
#ASV33 0.363  0.0054 **
  
#  Group DxD  #sps.  2 
#       stat p.value    
#ASV39 0.373  0.0007 ***
#ASV24 0.332  0.0207 *  
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#Run analysis by location.
inv2 <- multipatt(rel_com, location, func = "r.g", control = how(nperm = 9999))
summary(inv2)

#Multilevel pattern analysis
#---------------------------
  
#  Association function: r.g
#Significance level (alpha): 0.05

#Total number of species: 705
#Selected number of species: 14 
#Number of species associated to 1 group: 14 

#List of species associated to each combination: 
  
#  Group interspace  #sps.  9 
#        stat p.value    
#ASV6   0.360  0.0003 ***
#ASV28  0.315  0.0002 ***
#ASV24  0.278  0.0031 ** 
#ASV4   0.274  0.0072 ** 
#ASV38  0.257  0.0124 *  
#ASV110 0.241  0.0261 *  
#ASV21  0.229  0.0177 *  
#ASV51  0.218  0.0118 *  
#ASV36  0.203  0.0327 *  
  
#  Group under  #sps.  5 
#stat p.value    
#ASV1  0.590  0.0001 ***
#ASV5  0.461  0.0001 ***
#ASV15 0.280  0.0005 ***
#ASV52 0.184  0.0052 ** 
#ASV57 0.148  0.0261 *  
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
