######################################################################
######### Sula leucogaster's  Biogeographical Patterns in Diet #######
######################################################################

rm(list=ls()) #clean dataset
graphics.off()

setwd("C:/Users/Outra pessoa/Desktop/Sula_leucogaster/R/final_leuco_script")

#### OBJECTS GUIDE ####

#abund_prey - complete prey matrix

##FAMILIES###
#reg.fam - complete matrix
#data.fam - richness
#fam.col - resume matrix
#fam.stats - metrics

##SPECIES###
#leuco - complete clean matrix (whithout non-ID preys)
#data.spp - richness, shannon, shelf, biogeographic zonings, clusters
#spp.col - resumed matrix
#spp.stats - only metrics (occurrence, FO, abund, abund per colony)
#spp.bin - complete binary matrix (all studies)
#col.bin - resumed binary matrix

#### 1. REGURGITATE ANALYSIS ###########################################################################
#
#### 1.1 Prey abundance ####
abund_prey <- read.table(file="species_prey_leuco.csv", sep=';', h=T)
names(abund_prey)
abund_prey$total_abund <- rowSums(leuco[ , 14:152])
sum(abund_prey$total_abund)

lat.order <- factor(abund_prey$colony, levels=c("MS", "CU", "CA", "CF", "ST", "AB", "FN",
                                            "RO", "SPSP"))
tiff("prey_abund.tiff", height=20, width=30, 
     units='cm', compression="lzw", res=300)
ggplot(abund_prey, aes(x=as.factor(lat.order), y=total_abund)) + 
  geom_boxplot(fill="cornflowerblue", alpha=0.5) + xlab("Colony") + ylab("Total prey") +
  theme(text = element_text(size=15))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                            size= 15, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                            size= 15, colour="black", face="bold")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1135))
dev.off()   



##### 1.2 FAMILIES ####
###### 1.2.1 General metrics  ####

reg.fam <- read.table(file="families_prey_leuco.csv", sep=';', h=T)


#dataset structure
head(reg.fam)
dim(reg.fam)
str(reg.fam)
names(reg.fam)

############ richness #####
library (ggplot2)
library(vegan)

reg.fam$colony <- c("MS", "MS", "CU", "CA", "CF", "ST", "AB","AB","AB", "FN", "FN", 
                   "RO", "RO", "SPSP","SPSP","SPSP")

specnumber_f <- specnumber(reg.fam [,12:53])
site <- reg.fam$colony
shelf <-reg.fam$shelf
data.fam <- data.frame(specnumber_f,site,shelf)

lat.order.levels <- c("MS", "CU", "CA", "CF", "ST", "AB", "FN",
                      "RO", "SPSP")

site.data <- data.frame(site = lat.order.levels, 
                        site_num = 1:length(lat.order.levels)) #make another dataframe to map the continuous column

df <- dplyr::left_join(data.fam, site.data, by= "site") #merge the two dataframes

tiff("n_families.tiff", height=15, width=20, 
     units='cm', compression="lzw", res=300)
ggplot(df, aes(x=site_num, y=specnumber_f))+
  geom_point()+
  geom_smooth()+ scale_x_continuous(
    breaks=1:length(lat.order.levels),
    labels=lat.order.levels) +
  xlab("Colonies") + ylab("Number of families") +
  theme(text = element_text(size=15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                            size= 15, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                            size= 15, colour="black", face="bold")) +
  theme(legend.text = element_text(colour="black", size=15)) +
  theme(legend.title = element_blank()) 
dev.off()

############ frequency of occurence ####
library(plyr)


fam.col <- sapply(reg.fam[, 12:53], tapply, INDEX=reg.fam$site, sum, na.rm=T)
fam.col <- data.frame(fam.col)
fam.col$site <- c("AB", "CF", "CA", "CU", "MS", "FN", "RO", 
                  "ST", "SPSP")
names(fam.col)
fam.oc <- ldply(fam.col[,1:42], function(c) sum(c >0)) 
fam.stats <- data.frame (fam.oc)
colnames(fam.stats) <- c("family", "occ")

fam.stats$freq_occ <- (fam.stats$occ*100)/9 
sort(fam.stats$freq_occ) 
fam.stats[order(fam.stats[,"freq_occ"]), , drop=FALSE] #rank

############ abundance ####
names(fam.col)
fam.stats$fam.abund <- apply(fam.col[,1:42], 2, sum, na.rm=TRUE)
fam.stats[order(fam.stats[,"fam.abund"]), , drop=FALSE] #rank
 
############ relative abundance ####
sum(fam.stats$fam.abund) # total prey considered 
fam.stats$tot.rel.abu <- ((fam.stats$fam.abund*100)/sum(fam.stats$fam.abund)) 
head(fam.stats)
fam.stats[order(fam.stats[,"tot.rel.abu"]), , drop=FALSE] #rank
sum(fam.stats$tot.rel.abu) #sum 100%

############ abundance per colony ####

fam.stats$abund.Abrolhos <- apply(fam.col[1,1:42], 2, sum, na.rm=TRUE)
fam.stats$abund.Cabo_Frio <- apply(fam.col[2,1:42], 2, sum, na.rm=TRUE)
fam.stats$abund.Cagarras <- apply(fam.col[3,1:42], 2, sum, na.rm=TRUE)
fam.stats$abund.Currais <- apply(fam.col[4,1:42], 2, sum, na.rm=TRUE)
fam.stats$abund.Moleques <- apply(fam.col[5,1:42], 2, sum, na.rm=TRUE)
fam.stats$abund.Noronha <- apply(fam.col[6,1:42], 2, sum, na.rm=TRUE)
fam.stats$abund.Rocas <- apply(fam.col[7,1:42], 2, sum, na.rm=TRUE)
fam.stats$abund.Santana <- apply(fam.col[8,1:42], 2, sum, na.rm=TRUE)
fam.stats$abund.SPSP <- apply(fam.col[9,1:42], 2, sum, na.rm=TRUE)

############ relative abundance per colony ####

fam.stats$relabund.Abrolhos <- ((fam.stats$abund.Abrolhos*100)/sum(fam.stats$abund.Abrolhos)) 
fam.stats$relabund.Cabo_Frio <- ((fam.stats$abund.Cabo_Frio*100)/sum(fam.stats$abund.Cabo_Frio))
fam.stats$relabund.Cagarras<- ((fam.stats$abund.Cagarras*100)/sum(fam.stats$abund.Cagarras)) 
fam.stats$relabund.Currais <- ((fam.stats$abund.Currais*100)/sum(fam.stats$abund.Currais)) 
fam.stats$relabund.Moleques <- ((fam.stats$abund.Moleques*100)/sum(fam.stats$abund.Moleques)) 
fam.stats$relabund.Noronha <- ((fam.stats$abund.Noronha*100)/sum(fam.stats$abund.Noronha)) 
fam.stats$relabund.Rocas <- ((fam.stats$abund.Rocas*100)/sum(fam.stats$abund.Rocas)) 
fam.stats$relabund.Santana <- ((fam.stats$abund.Santana *100)/sum(fam.stats$abund.Santana)) 
fam.stats$relabund.SPSP <- ((fam.stats$abund.SPSP*100)/sum(fam.stats$abund.SPSP)) 
names(fam.stats)

#PLOT of prey families consumed by 
#Sula leucogaster with relative abundance (RAC) above 10% in each colony.

names(fam.stats)

MS <- data.frame (fam.stats[which(fam.stats$relabund.Moleques>10), c(1,19)])
MS$site <- rep("MS")
colnames(MS) <- c("family","relabund", "site")
CU <- data.frame (fam.stats[which(fam.stats$relabund.Currais>10), c(1,18)])
CU$site <- rep("CU")
colnames(CU) <- c("family","relabund", "site")
CA <- data.frame (fam.stats[which(fam.stats$relabund.Cagarras>10), c(1,17)])
CA$site <- rep("CA")
colnames(CA) <- c("family","relabund", "site")
CF <- data.frame (fam.stats[which(fam.stats$relabund.Cabo_Frio>10), c(1,16) ])
CF$site <- rep("CF")
colnames(CF) <- c("family","relabund", "site")
ST <- data.frame (fam.stats[which(fam.stats$relabund.Santana>10), c(1,22)])
ST$site <- rep("ST")
colnames(ST) <- c("family","relabund", "site")
AB <- data.frame (fam.stats[which(fam.stats$relabund.Abrolhos>10), c(1,15)])
AB$site <- rep("AB")
colnames(AB) <- c("family","relabund", "site")
FN <- data.frame (fam.stats[which(fam.stats$relabund.Noronha>10),c(1,20)])
FN$site <- rep("FN")
colnames(FN) <- c("family","relabund", "site")
RO <- data.frame (fam.stats[which(fam.stats$relabund.Rocas>10), c(1,21)])
RO$site <- rep("RO")
colnames(RO) <- c("family","relabund", "site")
SPSP <- data.frame (fam.stats[which(fam.stats$relabund.SPSP>10), c(1,23)])
SPSP$site <- rep("SPSP")
colnames(SPSP) <- c("family","relabund", "site")

abund_10 <- rbind(MS, CU, CA, CF, ST, AB, FN, RO, SPSP)

lat.order <- factor(abund_10$site, levels=c("MS", "CU", "CA", "CF", "ST", "AB", "FN",
                                              "RO", "SPSP"))
tiff("abund_10.tiff", height=20, width=25, 
     units='cm', compression="lzw", res=300)
par(mfrow = c(1,1))
ggplot(abund_10, aes(y=relabund, x=lat.order, fill=family)) + xlab("Colony") + ylab("Relative Abundance (%)") +
  geom_bar(position="stack", stat="identity") + scale_fill_manual("family", values = c("Batrachoididae" = "#8dd3c7", 
                                                                                       "Carangidae" =  "#b3de69",
                                                                                       "Clupeidae" = "steelblue",
                                                                                       "Engraulidae" = "lightgoldenrod1",
                                                                                       "Exocoetidae" = "#fb8072", 
                                                                                       "Hemiramphidae" = "#fdb462",
                                                                                       "Malacanthidae" = "violet", 
                                                                                       "Pristigasteridae" = "#fccde5",
                                                                                       "Sciaenidae" = "seagreen3",
                                                                                       "Synodontidae" = "#bc80bd")) +
  theme(text = element_text(size=15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 15), 
                                    size= 15, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                    size= 15, colour="black", face="bold")) +
  theme(legend.text = element_text(colour="black", size=15)) + theme(legend.title = element_blank())
dev.off()

##### 1.3 SPECIES ####
leuco <- read.table(file="species_prey_leuco.csv", sep=';', h=T)

## Removing columns that we didn't want (non-ID prey)

library(tidyverse)
leuco <- select(leuco, -(ends_with(c("nID","_sp","_spp", "Sula_leucogaster"))))
names(leuco)
leuco$colony <- c("MS", "MS", "CU", "CA", "CF", "ST", "AB","AB","AB", "FN", "FN", 
                    "RO", "RO", "SPSP","SPSP","SPSP")


write.csv(leuco,"species_prey_leuco_CLEAN.csv", row.names= F)

head(leuco)
dim(leuco)
str(leuco)
###### 1.3.1 General metrics ####
############ abundance in each study and total ####

names(leuco)

leuco <- read.table(file="species_prey_leuco_CLEAN.csv", sep=',', h=T)
leuco$total_abund <- rowSums(leuco[ , 14:113])

sum(leuco$total_abund)
############ frequency of occurence ####
library(plyr)

#join all the species prey in each colony, making a resumed matrix with one colony per row

names(leuco)
spp.col <- sapply(leuco[, 14:113], tapply, INDEX=leuco$colony, sum, na.rm=T)
spp.col <- data.frame(spp.col)

names(spp.col)

spp.oc <- ldply(spp.col[,1:100], function(c) sum(c >0)) 
spp.stats <- data.frame (spp.oc)
colnames(spp.stats) <- c("species", "occ")

spp.stats$freq_occ <- (spp.stats$occ*100)/9 
sort(spp.stats$freq_occ) 
spp.stats[order(spp.stats[,"freq_occ"]), , drop=FALSE] #rank

############ abundance ####
names(spp.col)
spp.stats$spp.abund <- apply(spp.col[,1:100], 2, sum, na.rm=TRUE)
spp.stats[order(spp.stats[,"spp.abund"]), , drop=FALSE] #rank

############ relative abundance ####
sum(spp.stats$spp.abund) # total prey considered 
spp.stats$tot.rel.abu <- ((spp.stats$spp.abund*100)/sum(spp.stats$spp.abund)) 
head(spp.stats)
spp.stats[order(spp.stats[,"tot.rel.abu"]), , drop=FALSE] #rank
sum(spp.stats$tot.rel.abu) #sum 100%

############ abundance per colony ####
spp.stats$abund.Abrolhos <- apply(spp.col[1,1:100], 2, sum, na.rm=TRUE)
spp.stats$abund.Cabo_Frio <- apply(spp.col[2,1:100], 2, sum, na.rm=TRUE)
spp.stats$abund.Cagarras <- apply(spp.col[3,1:100], 2, sum, na.rm=TRUE)
spp.stats$abund.Currais <- apply(spp.col[4,1:100], 2, sum, na.rm=TRUE)
spp.stats$abund.Moleques <- apply(spp.col[5,1:100], 2, sum, na.rm=TRUE)
spp.stats$abund.Noronha <- apply(spp.col[6,1:100], 2, sum, na.rm=TRUE)
spp.stats$abund.Rocas <- apply(spp.col[7,1:100], 2, sum, na.rm=TRUE)
spp.stats$abund.Santana <- apply(spp.col[8,1:100], 2, sum, na.rm=TRUE)
spp.stats$abund.SPSP <- apply(spp.col[9,1:100], 2, sum, na.rm=TRUE)

############ relative abundance per colony ####
spp.stats$relabund.Abrolhos <- ((spp.stats$abund.Abrolhos*100)/sum(spp.stats$abund.Abrolhos)) 
spp.stats$relabund.Cabo_Frio <- ((spp.stats$abund.Cabo_Frio *100)/sum(spp.stats$Cabo_Frio))
spp.stats$relabund.Cagarras<- ((spp.stats$abund.Cagarras*100)/sum(spp.stats$abund.Cagarras)) 
spp.stats$relabund.Currais <- ((spp.stats$abund.Currais*100)/sum(spp.stats$abund.Currais)) 
spp.stats$relabund.Moleques <- ((spp.stats$abund.Moleques*100)/sum(spp.stats$abund.Moleques)) 
spp.stats$relabund.Noronha <- ((spp.stats$abund.Noronha*100)/sum(spp.stats$abund.Noronha)) 
spp.stats$relabund.Rocas <- ((spp.stats$abund.Rocas*100)/sum(spp.stats$abund.Rocas)) 
spp.stats$relabund.Santana <- ((spp.stats$abund.Santana *100)/sum(spp.stats$abund.Santana)) 
spp.stats$relabund.SPSP <- ((spp.stats$abund.SPSP*100)/sum(spp.stats$abund.SPSP)) 
names(spp.stats)

############ richness ####
library (ggplot2)
library(vegan)
names(leuco)
specnumber_s <- specnumber(leuco[,14:113])
site <- leuco$colony
shelf <-leuco$shelf
data.spp <- data.frame(specnumber_s,site,shelf)

lat.order.levels <- c("MS", "CU", "CA", "CF", "ST", "AB", "FN",
                      "RO", "SPSP")

site.spp <- data.frame(site = lat.order.levels, 
                        site_num = 1:length(lat.order.levels)) #make another dataframe to map the continuous column

df.spp <- dplyr::left_join(data.spp, site.spp, by= "site") #merge the two dataframes

tiff("richness.tiff", height=15, width=20, 
     units='cm', compression="lzw", res=300)
ggplot(df.spp, aes(x=site_num, y=specnumber_s))+
  geom_point()+
  geom_smooth()+ scale_x_continuous(
    breaks=1:length(lat.order.levels),
    labels=lat.order.levels) +
  xlab("Colonies") + ylab("Richness") +
  theme(text = element_text(size=15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                    size= 15, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                    size= 15, colour="black", face="bold")) +
  theme(legend.text = element_text(colour="black", size=15)) +
  theme(legend.title = element_blank()) 
dev.off() 

###### 1.3.2 Shannon Index ####
# Shannon

data.spp$shannon.div <- diversity(leuco[,14:113], index = "shannon", 
                                   MARGIN = 1, base = exp(1)) 

plot1 <- ggplot(data.spp, aes(x=as.factor(lat.order), y=shannon.div)) + 
  geom_boxplot(fill="cornflowerblue", alpha=0.5) + xlab("Colony") + ylab("Shannon Index") +
  theme(text = element_text(size=15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                    size= 15, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                    size= 15, colour="black", face="bold")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
  theme(legend.text = element_text(colour="black", size=15)) +
  theme(legend.title = element_blank()) 

# Shannon X Shelf

plot2 <- ggplot(data.spp, aes(x=as.factor(shelf), y=shannon.div)) + 
  geom_boxplot(fill="chartreuse3", alpha=0.5) + xlab("Shelf") + ylab("") +
  theme(text = element_text(size=15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                    size= 15, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                    size= 15, colour="black", face="bold")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
  theme(legend.text = element_text(colour="black", size=15)) +
  theme(legend.title = element_blank()) 
  
#Shannon X cluster K=2
data.spp$cluster_k2 <-leuco$cluster_k2

plot3 <-ggplot(data.spp, aes(x=as.factor(cluster_k2), y=shannon.div)) + 
  geom_boxplot(fill="darkorchid2", alpha=0.5) + xlab("Cluster K=2") + ylab("Shannon Index") +
  theme(text = element_text(size=15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                    size= 15, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                    size= 15, colour="black", face="bold")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
  theme(legend.text = element_text(colour="black", size=15)) +
  theme(legend.title = element_blank())

# Shannon X cluster K=3

data.spp$cluster_k3 <-leuco$cluster_k3

cluster.order <- factor(data.spp$cluster_k3, levels=c("Coastal", "FN_RO", "SPSP")) 

plot4 <-ggplot(data.spp, aes(x=as.factor(cluster.order), y=shannon.div)) + 
  geom_boxplot(fill="darkorchid2", alpha=0.5) + xlab("Cluster K=3") + ylab("") +
  theme(text = element_text(size=15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                    size= 15, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                    size= 15, colour="black", face="bold")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
  theme(legend.text = element_text(colour="black", size=15)) +
  theme(legend.title = element_blank())

library(gridExtra)
tiff("shannon_shelf_cluster.tiff", height=30, width=35, 
     units='cm', compression="lzw", res=300)
par(mfrow = c(1,1)) 
grid.arrange(plot1, plot2, plot3, plot4, nrow=2)
dev.off()

#Shannon X LME

data.spp$LME <-leuco$LME

LME.order <- factor(data.spp$LME, levels=c("SBS15", "EBS16"))

plot5 <- ggplot(data.spp, aes(x=as.factor(LME.order), y=shannon.div)) + 
  geom_boxplot(fill="cornflowerblue", alpha=0.5) + xlab("LME") + ylab("Shannon Index") +
  theme(text = element_text(size=15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                    size= 15, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                    size= 15, colour="black", face="bold")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
  theme(legend.text = element_text(colour="black", size=15)) +
  theme(legend.title = element_blank())

#Shannon X Spalding's Ecorregions

data.spp$MEOW_ECO <-leuco$MEOW_ECO

MEOWECO.order <- factor(data.spp$MEOW_ECO, levels=c("WTSA180", "TSA76", "TSA74", "TSA73"))

plot8 <- ggplot(data.spp, aes(x=as.factor(MEOWECO.order), y=shannon.div)) + 
  geom_boxplot(fill="darkorchid2", alpha=0.5) + xlab("Spalding's Ecorregion") + ylab("") +
  theme(text = element_text(size=15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                    size= 15, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                    size= 15, colour="black", face="bold")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
  theme(legend.text = element_text(colour="black", size=15)) +
  theme(legend.title = element_blank())

#Shannon X Spaldin's Provinces

data.spp$MEOW_PROV <-leuco$MEOW_PROV

MEOWPROV.order <- factor(data.spp$MEOW_PROV, levels=c("WTSA", "TSA"))

plot7 <-ggplot(data.spp, aes(x=as.factor(MEOWPROV.order), y=shannon.div)) + 
  geom_boxplot(fill="darkorchid2", alpha=0.5) + xlab("Spalding's Province") + ylab("Shannon Index") +
  theme(text = element_text(size=15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                    size= 15, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                    size= 15, colour="black", face="bold")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
  theme(legend.text = element_text(colour="black", size=15)) +
  theme(legend.title = element_blank())

#Shannon X Longhurst's Provinces

data.spp$longhurst <-leuco$longhurst

long.order <- factor(data.spp$longhurst, levels=c("BRAZ", "SATL", "WTRA"))

plot6 <-ggplot(data.spp, aes(x=as.factor(long.order), y=shannon.div)) + 
  geom_boxplot(fill="chartreuse3", alpha=0.5) + xlab("Longhurst's Province") + ylab("") +
  theme(text = element_text(size=15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                    size= 15, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                    size= 15, colour="black", face="bold")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
  theme(legend.text = element_text(colour="black", size=15)) +
  theme(legend.title = element_blank())

tiff("shannon_biogeo.tiff", height=30, width=35, 
     units='cm', compression="lzw", res=300)
par(mfrow = c(1,1)) 
grid.arrange(plot5, plot6, plot7, plot8, nrow=2)
dev.off()

###### 1.3.3 accumulation curve ####

sac <- specaccum(leuco[,14:113], method = "random")

tiff("accumulation_curve.tiff", height=30, width=35, 
     units='cm', compression="lzw", res=300)
plot(sac, ci.type="poly", col="lightseagreen", lwd=2, ci.lty=0, ci.col="lightskyblue", xlab = "Datasets", ylab="Number of Species")
boxplot(sac, col="lightcoral", add=TRUE, pch="+")
dev.off()

###### 1.3.4 PERMANOVA ####

#creating a binary matrix
spp.bin <- ifelse(leuco[,14:113]>0, 1, 0)
spp.bin <- data.frame(spp.bin)
spp.bin$shelf <- leuco$shelf
spp.bin$cluster_k2 <- leuco$cluster_k2
spp.bin$cluster_k3 <- leuco$cluster_k3
spp.bin$site <- leuco$site
spp.bin$LME <- leuco$LME
spp.bin$MEOW_ECO <- leuco$MEOW_ECO
spp.bin$MEOW_PROV <- leuco$MEOW_PROV
spp.bin$LONG <- leuco$longhurst

permanova.LME <- adonis(spp.bin[,1:100] ~ spp.bin$LME, permutations=999, 
                        distance='bray')
permanova.LME

permanova.MEOW_ECO <- adonis(spp.bin[,1:100] ~ spp.bin$MEOW_ECO, permutations=999, 
                             distance='bray')
permanova.MEOW_ECO

permanova.MEOW_PROV <- adonis(spp.bin[,1:100] ~ spp.bin$MEOW_PROV, permutations=999, 
                              distance='bray')
permanova.MEOW_PROV

permanova.LONG <- adonis(spp.bin[,1:100] ~ spp.bin$LONG, permutations=999, 
                         distance='bray')
permanova.LONG

permanova.SHELF <- adonis(spp.bin[,1:100] ~ spp.bin$shelf, permutations=999, 
                          distance='bray')
permanova.SHELF

permanova.K2 <- adonis(spp.bin[,1:100] ~ spp.bin$cluster_k2, permutations=999, 
                       distance='bray')
permanova.K2

permanova.K3 <- adonis(spp.bin[,1:100] ~ spp.bin$cluster_k3, permutations=999, 
                       distance='bray')
permanova.K3

###### 1.3.5 nMDS ####
#creating a binary matrix from the complete matrix (all studies)

names(leuco)

spp.bin <- ifelse(leuco[,14:113]>0, 1, 0)
head(spp.bin)
spp.bin <- data.frame(spp.bin)
spp.bin$shelf <- leuco$shelf
spp.bin$cluster_k2 <- leuco$cluster_k2
spp.bin$cluster_k3 <- leuco$cluster_k3
spp.bin$site <- leuco$site
spp.bin$colony <- leuco$colony
spp.bin$LME <- leuco$LME
spp.bin$MEOW_ECO <- leuco$MEOW_ECO
spp.bin$MEOW_PROV <- leuco$MEOW_PROV
spp.bin$LONG <- leuco$longhurst
names(spp.bin)
head(spp.bin)

#nMDS
library (vegan)
set.seed(123)
initial_nMDS1 <- metaMDS(spp.bin[,1:100], distance="bray", k=2, trymax=1000)
nMDS1 <- metaMDS(spp.bin[,1:100], previous.best = initial_nMDS1, k=2, trymax=1000)
nMDS1
#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nMDS1))
data.scores$shelf <- spp.bin$shelf
data.scores$cluster_k2 <- spp.bin$cluster_k2
data.scores$cluster_k3 <- spp.bin$cluster_k3
data.scores$site <- spp.bin$site
data.scores$colony <- spp.bin$colony
data.scores$LME <- spp.bin$LME
data.scores$MEOW_ECO <- spp.bin$MEOW_ECO
data.scores$MEOW_PROV <- spp.bin$MEOW_PROV
data.scores$LONG <- spp.bin$LONG
head(data.scores)

# Plot the stress
stressplot(nMDS1)

#Shelf
library(dplyr)
library(stringr)
library(tidyverse)

par(mfrow = c(2,2)) 

hull <- data.scores %>% group_by(shelf) %>%
  slice(c(chull(NMDS1,NMDS2),
          chull(NMDS1,NMDS2)))
plot1 <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = shelf), size = 2, alpha = 5) + 
  scale_colour_manual(values = c("indianred1", "dodgerblue")) + xlab("") + ylab("NMDS2")+
  theme(text = element_text(size=15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                    size= 12, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                    size= 12, colour="black", face="bold")) +
  theme(legend.text = element_text(colour="black", size=15)) + 
  stat_ellipse(geom = "polygon", aes(fill = shelf),alpha = 0.25, level=0.75) +
  theme(legend.title = element_blank())

#Cluster K=2
hull <- data.scores %>% group_by(cluster_k2) %>%
  slice(c(chull(NMDS1,NMDS2),
          chull(NMDS1,NMDS2)))
plot2 <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = cluster_k2), size = 2, alpha = 5) + 
  xlab("NMDS1") + ylab("")+
  scale_colour_manual(values = c("indianred1", "dodgerblue")) + 
  theme(text = element_text(size=15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                    size= 12, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                    size= 12, colour="black", face="bold")) +
  theme(legend.text = element_text(colour="black", size=15)) + 
  stat_ellipse(geom = "polygon", aes(fill = cluster_k2),alpha = 0.25, level=0.75) +
  theme(legend.title = element_blank())

#Cluster K=3

hull <- data.scores %>% group_by(cluster_k3) %>%
  slice(c(chull(NMDS1,NMDS2),
          chull(NMDS1,NMDS2)))
plot3 <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = cluster_k3), size = 2, alpha = 5) +
  xlab("NMDS1") + ylab("NMDS2")+
  scale_colour_manual(values = c("indianred1", "dodgerblue", "sienna1")) + 
  theme(text = element_text(size=15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                    size= 12, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                    size= 12, colour="black", face="bold")) +
  theme(legend.text = element_text(colour="black", size=15))  +
  stat_ellipse(geom = "polygon", aes(fill = cluster_k3),alpha = 0.25, level=0.75) +
  theme(legend.title = element_blank())


library(gridExtra)

tiff("plot_ndms_1.tiff", height=30, width=35, 
     units='cm', compression="lzw", res=300)
grid.arrange(plot1, plot2, plot3, nrow=2)
dev.off() 

#LME
hull <- data.scores %>% group_by(LME) %>%
  slice(c(chull(NMDS1,NMDS2),
          chull(NMDS1,NMDS2)))
plot4 <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = LME), size = 2, alpha = 5) + xlab("") + ylab("NMDS2") +
  scale_colour_manual(values = c("indianred1", "dodgerblue", "sienna1")) + 
  theme(text = element_text(size=15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                    size= 12, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                    size= 12, colour="black", face="bold")) +
  theme(legend.text = element_text(colour="black", size=15)) + 
  stat_ellipse(geom = "polygon", aes(fill = LME),alpha = 0.25, level=0.75) +
  theme(legend.title = element_blank())

#Spalding's Provinces

hull <- data.scores %>% group_by(MEOW_PROV) %>%
  slice(c(chull(NMDS1,NMDS2),
          chull(NMDS1,NMDS2)))
plot6 <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = MEOW_PROV), size = 2, alpha = 5) + xlab("NMDS1") + ylab("NMDS2") +
  scale_colour_manual(values = c("indianred1", "dodgerblue", "sienna1")) + 
  theme(text = element_text(size=15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                    size= 12, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                    size= 12, colour="black", face="bold")) +
  theme(legend.text = element_text(colour="black", size=15)) +
  stat_ellipse(geom = "polygon", aes(fill = MEOW_PROV),alpha = 0.25, level=0.75) +
  theme(legend.title = element_blank())


#Spalding's Ecoregions

hull <- data.scores %>% group_by(MEOW_ECO) %>%
  slice(c(chull(NMDS1,NMDS2),
          chull(NMDS1,NMDS2)))

plot7 <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = MEOW_ECO), size = 2, alpha = 5) + xlab("NMDS1") + ylab("") +
  scale_colour_manual(values = c("sienna1","mediumseagreen", "dodgerblue","indianred1")) + 
  theme(text = element_text(size=15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                    size= 12, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                    size= 12, colour="black", face="bold")) +
  theme(legend.text = element_text(colour="black", size=15)) +
  stat_ellipse(geom = "polygon", aes(fill = MEOW_ECO),alpha = 0.25, level=0.75) +
  theme(legend.title = element_blank())

#Longhurst's Provinces

hull <- data.scores %>% group_by(LONG) %>%
  slice(c(chull(NMDS1,NMDS2),
          chull(NMDS1,NMDS2)))

plot5 <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = LONG), size = 2, alpha = 5) + xlab("") + ylab("") +
  scale_colour_manual(values = c("indianred1", "mediumseagreen", "dodgerblue")) + 
  theme(text = element_text(size=15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                    size= 12, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                    size= 12, colour="black", face="bold")) +
  theme(legend.text = element_text(colour="black", size=15)) +
  stat_ellipse(geom = "polygon", aes(fill = LONG),alpha = 0.25, level=0.75) +
  theme(legend.title = element_blank())

library(gridExtra)


tiff("grafico_ndms_2.tiff", height=30, width=35, 
     units='cm', compression="lzw", res=300)
grid.arrange(plot4, plot5, plot6, plot7, nrow=2)
dev.off() 

###### 1.3.6 Dendogram ####

#making a binary matrix from the resumed matrix (one colony per row)
head(spp.col) #resumed matrix
names(spp.col)

col.bin <- ifelse(spp.col[,1:100]>0, 1, 0)

col.bin <- data.frame(col.bin)

head(col.bin)

names(col.bin)

library(ggdendro)
library(dendextend)
library(betapart)

# Beta diversity indexes - betapart package

betapair <- beta.pair(col.bin[,1:100], index.family = "sorensen")
betapair
beta.sim <- betapair$beta.sim ###  Simpson pair-wise dissimilarity
beta.sor <- betapair$beta.sor ###  Sorensen pair-wise dissimilarity

# Hierarchical Clustering with hclust

hc <- hclust(beta.sim)
dd <- as.dendrogram(hc)
dendr <- dendro_data(dd, type="rectangle")
ggdendrogram(dd)
pp <- click_rotate(dd, continue=TRUE)
dd <- ggdendrogram(rev(pp), rotate = TRUE, theme_dendro = TRUE)
dendro <- dd + theme(text = element_text(size = 17, colour="black", face="bold")) 

tiff("dendrogram.tiff", height=15, width=20, 
     units='cm', compression="lzw", res=300)
dendro
dev.off() 
##### 2. STABLE ISOTOPES ANALYSIS ############################################################################
###### 2.1 General metrics ####

library(FSA)

data.iso <- read.csv(file="isotopes_siber.csv", sep=";", h=T)
attach(data.iso)

C_met <-Summarize(iso1 ~ group,
                  data = data.iso, digits=2)
C_met

N_met <-Summarize(iso2 ~ group,
                  data = data.iso, digits=2)
N_met

###### 2.2 Kruskall-Wallis/Mann-Witney ####

data.iso2 <- read.csv(file="isotopes_cgd.csv", sep=";", h=T)

attach(data.iso2)

# carbon x colonies
KWC <- kruskal.test(d_carbon~group)
KWC
## pairwise
pairwise.wilcox.test(d_carbon, group, p.adjust.method = "fdr")

# nitrogen x colonies
KWN <- kruskal.test(d_nitrogen~group)
KWN
## pairwise
pairwise.wilcox.test(d_nitrogen, group, p.adjust.method = "fdr")

#nitrogen x criteria for definig groups

##shelf
pairwise.wilcox.test(d_nitrogen, shelf, p.adjust.method = "fdr")

## clusters
pairwise.wilcox.test(d_nitrogen, cluster_k2, p.adjust.method = "fdr")
pairwise.wilcox.test(d_nitrogen, cluster_k3, p.adjust.method = "fdr")

## biogeographical zonings
### Longhusts's Biogeochemical Zones
pairwise.wilcox.test(d_nitrogen, LONG, p.adjust.method = "fdr")
### Large Marine Ecosystems
pairwise.wilcox.test(d_nitrogen, LME, p.adjust.method = "fdr")
### Spalding's Provinces
pairwise.wilcox.test(d_nitrogen, MEOW_PROV, p.adjust.method = "fdr")
### Spalding's Ecorregions
pairwise.wilcox.test(d_nitrogen, MEOW_ECO, p.adjust.method = "fdr")

###### 2.3 SIBER ####
library(ggplot2)
library(hrbrthemes)
library(SIBER)
library(ellipse)
library(dplyr)

leuco_SIBER<- read.table("isotopes_siber.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE)
siber.leuco <- createSiberObject(leuco_SIBER)

# C-N points plot
ggplot(leuco_SIBER, aes(x=leuco_SIBER$iso1,leuco_SIBER$iso2, color=leuco_SIBER$group)) + 
  geom_point(size=2) +
  theme_ipsum() + ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) 

community.hulls.args <- list(col = 2, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

group.ML <- groupMetricsML(siber.leuco)
print(group.ML)

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

ellipses.posterior <- siberMVN(siber.leuco, parms, priors)

# calculate the SEA.B for each colony.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group")

#----create-ellipse-df-
# how many of the posterior draws do you want 
n.posts <- 5

# how big an ellipse you want to draw
p.ell <- 0.95

# a list to store the results
all_ellipses <- list()

# loop over groups
for (i in 1:length(ellipses.posterior)){
  
  # a dummy variable to build in the loop
  ell <- NULL
  post.id <- NULL
  
  for ( j in 1:n.posts){
    
    # covariance matrix
    Sigma  <- matrix(ellipses.posterior[[i]][j,1:4], 2, 2)
    
    # mean
    mu     <- ellipses.posterior[[i]][j,5:6]
    
    # ellipse points
    out <- ellipse::ellipse(Sigma, centre = mu , level = p.ell)
    ell <- rbind(ell, out)
    post.id <- c(post.id, rep(j, nrow(out)))}
  
  ell <- as.data.frame(ell)
  ell$rep <- post.id
  all_ellipses[[i]] <- ell}

ellipse_df <- bind_rows(all_ellipses, .id = "id")

# now we need the group and community names
# extract them from the ellipses.posterior list
group_comm_names <- names(ellipses.posterior)[as.numeric(ellipse_df$id)]

# split them and conver to a matrix, NB byrow = T
split_group_comm <- matrix(unlist(strsplit(group_comm_names, "[.]")),
                           nrow(ellipse_df), 2, byrow = TRUE)

ellipse_df$community <- split_group_comm[,1]
ellipse_df$group     <- split_group_comm[,2]

ellipse_df <- dplyr::rename(ellipse_df, iso1 = x, iso2 = y)

## plot-data-
first.plot <- ggplot(data = ellipse_df, aes(iso1, iso2)) +
  geom_point(aes(color = factor(group)), size = 1)+ ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + labs(colour="Colony") +
  theme(text = element_text(size=15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                    size= 15, colour="black", face="bold")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0),
                                    size= 15, colour="black", face="bold")) +
  theme(legend.text = element_text(colour="black", size=15)) + 
          theme(legend.title= element_text(colour="black", size=15, face="bold")) 

tiff("siber1.tiff", height=15, width=20, 
     units='cm', compression="lzw", res=300)
first.plot
dev.off() 



second.plot <- first.plot + facet_wrap(~factor(group))
print(second.plot)

# rename columns of ellipse_df to match the aesthetics

third.plot <- second.plot + 
  geom_polygon(data = ellipse_df, mapping = aes(iso1, iso2, group = group, color = factor(group),
                             fill = NULL), fill = NA, alpha = 0.2)
print(third.plot)

tiff("siber2.tiff", height=15, width=20, 
     units='cm', compression="lzw", res=300)
third.plot
dev.off() 

# to calculate the overlap area between the ellipses
ellipse1 <- "1.Cagarras" 
ellipse2 <- "1.Abrolhos"
# the overlap betweeen the corresponding 95% prediction ellipses is given by: 
ellipse95.overlap <- maxLikOverlap(ellipse1, ellipse2, siber.leuco,
                                   p.interval = 0.95, n = 100)
ellipse95.overlap


#################################### the end ##########################################

