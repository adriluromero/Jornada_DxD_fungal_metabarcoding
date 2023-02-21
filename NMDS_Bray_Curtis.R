#Jornada_DxD_Fungal_Metabarcoding

#Protocol for processing Metabarcoding of soil
#File name: NMDS_Bray_Curtis


#Clear environment
rm(list = ls(all.names = TRUE))

#Load packages
library(vegan)
library(funrar) #For relative abundance
library(ggplot2)

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
