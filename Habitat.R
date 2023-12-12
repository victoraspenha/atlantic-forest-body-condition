##################################################################################################################################################################################################
##################################################################################################################################################################################################
## Removing leftovers
rm(list=ls())

## Needed packages
library(Rmisc)
library(ggExtra)
library(lubridate)
library(rstanarm)
library(data.table)
library(forecast)
library(nlme)
library(ape)
library(picante)
library(phytools)
library(wesanderson)
library(diversitree)
library(RColorBrewer)
library(grDevices)
library(BBmisc)
library(phyr)
library(INLA)
library(devtools)
library(rr2)
library(MuMIn)
library(lme4)
library(phyloch)
library(sf)
library(rgdal)
library(raster)  
library(CGPfunctions)
library(precrec)
library(tidyverse)
library(regclass)
library(ggpubr)
library(lmerTest)
library(smatr)
library(magrittr)
library(MASS)
library(cluster)
library(factoextra)
library(ggplot2)
library(rlang)
library(scales)
library(FactoMineR)
library(heatmaply)
library(gridExtra)
library(effectsize)
library(car)
library(rstatix)
library(Hmisc)
library(PerformanceAnalytics)
library(corrplot)
library(caret)
library(sjstats)
library(knitr)
library(dplyr)
library(tidyr)
library(clustertend)
library(hopkins)
library(flexclust)
library(fpc)
library(ClusterR)
library(clusterSim)
library(smacof)
library(paran)
library(stats)
library(rworldmap)
library(ClustGeo)
library(readxl)	
library(rio)
library(caTools)
library(jtools) 
library(interactions)
library(effects)
library(rgl)
library(geiger)
library(phangorn)
library(tibble)
library(magrittr)
library(roperators) 
library(rgbif)
library(gdalUtils)
library(ncdf4)
library(elevatr)
library(phylobase)
library(grid)
library(phylosignal)
library(pwr)
require(subplex)

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
## Getting data -- you can find the data here: https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2647
data <- read.csv("ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv", h=T)
str(data)

## Species correction:
data <- data[-which(data$Binomial == "Elaenia sp."),]
data[which(data$Binomial == "Stephanoxis loddigesii"),"Binomial"] <- "Stephanoxis lalandi"
data[which(data$Binomial == "Campylopterus diamantinensis"),"Binomial"] <- "Campylopterus largipennis"
data[which(data$Binomial == "Campylopterus largipennis"),"Binomial"] <- "Campylopterus calcirupicola"
data[which(data$Binomial == "Pyriglena leuconota"),"Binomial"] <- "Pyriglena pernambucensis"

## Remove rows without longitude and latitude
data <- data[-which(is.na(data$Longitude_decimal_degrees)),]
which(is.na(data$Latitude_decimal_degrees))

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
## Habitat: https://www.earthenv.org/texture: 30 arc-second (1km)
mapa <- raster("./Dissimilarity_01_05_1km_uint32.tif")

## Getting habitat cover
occ <- unique(as.data.frame.matrix(data[,c("Longitude_decimal_degrees", "Latitude_decimal_degrees")]))
global.veg <- raster::extract(mapa, occ)
forest_cover <- cbind(occ, global.veg)

data$ForestCover <- NA
## Getting data
for(i in 1:dim(data)[1]){
  filtro <- data[i,c("Longitude_decimal_degrees", "Latitude_decimal_degrees")]
  result <- c(filtro[1,1], filtro[1,2])
  value <- forest_cover[which(forest_cover[,c("Longitude_decimal_degrees")] == result[1] & forest_cover[,c("Latitude_decimal_degrees")] == result[2]),"global.veg"]
  data[i,"ForestCover"] <- value
}

length(which(!complete.cases(data$ForestCover))) / dim(data)[1]
data <- data[-which(!complete.cases(data$ForestCover)),]
hist(log(data$ForestCover))

## CHecking species
test <- as.data.frame(sort(table(data$Binomial)))
test <- test[-which(test$Freq <= 3),]

## Testing frequency
data2 <- data[which(data$Binomial %in% test$Var1),]
sort(table(data2$Binomial))

## Body conditiom
length(which(is.na(data2$Body_mass.g.)))
length(which(is.na(data2$Body_length.mm.)))
length(which(is.na(data2$Wing_length_right.mm.)))
length(which(is.na(data2$Tarsus_length_right.mm.)))
length(which(is.na(data2$Bill_length.mm.)))

## Removing NAs data
data3 <- data2[-which(is.na(data2$Body_mass.g.)),]
data3 <- data3[-which(is.na(data3$Body_length.mm.)),]
dim(data3)

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
## Body condition: available from http://apansharing.blogspot.com/2018/05/an-r-function-olsrobust-caled-mass-index.html
## Credits to Chen-Pan Liao
## Testing again
test2 <- as.data.frame(sort(table(data3$Binomial)))
test2 <- test2[-which(test2$Freq <= 3),]
data4 <- data3[which(data3$Binomial %in% test2$Var1),]
sort(table(data4$Binomial))

## Body condition
data4$BodyCond <- NA
especies <- unique(data4$Binomial)
for(i in 1:length(especies)){
  print(especies[i])
  print("################################")
  dentro <- data4[which(data4$Binomial == especies[i]),]
  estados <- unique(dentro$State)
  
  for(j in 1:length(estados)){
    massa <- dentro[which(dentro$Binomial == especies[i] & dentro$State == estados[j]),"Body_mass.g."] 
    comprimento <- dentro[which(dentro$Binomial == especies[i] & dentro$State == estados[j]),"Body_length.mm."] 
    
    if((length(massa) <= 10) == TRUE){
      data4[which(data4$Binomial == especies[i] & data4$State == estados[j]),"BodyCond"] <- "NA"
      print("NOT ENOUGH DATA")
    }
    else{
      residuos <- round(scaledMassIndex(massa, comprimento)[,1],2)
      data4[which(data4$Binomial == especies[i] & data4$State == estados[j]),"BodyCond"] <- residuos
      print("########")
    }
    
  }
}

## Removing the ones that we did not have enough data
data4 <- data4[-which(data4$BodyCond == "NA"),]
data4$BodyCond <- as.numeric(data4$BodyCond)

## Variable distribution
hist(log(data4$BodyCond))

## Checking cofounders
any(is.na(data4$Year))
data4[which(is.na(data4$Year)),"Year"] <- "Unknown"
table(data4$Year)
any(is.na(data4$State))
table(data4$State)
any(is.na(data4$Sex))
data4[which(is.na(data4$Sex)),"Sex"] <- "Unknown"
table(data4$Sex)
any(is.na(data4$Age))
data4[which(is.na(data4$Age)),"Age"] <- "Unknown"
table(data4$Age)
data4[which(is.na(data4$Molt)),"Molt"] <- "Unknown"
table(data4$Molt)

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
dim(data4)
## Phylogeny: you can get phylogeny here: https://birdtree.org/subsets/
phylogeny <- read.nexus("/output.nex")
phylogeny <- maxCladeCred(phylogeny)

## Checking for similarities between phylogeny and dataset
data4$Binomial <- gsub(" ", "_", data4$Binomial)
phylogeny$tip.label[which(!phylogeny$tip.label %in% unique(data4$Binomial))]
unique(data4$Binomial)[which(!unique(data4$Binomial) %in% phylogeny$tip.label)]

## Correcting phylogeny
phylogeny$tip.label[which(phylogeny$tip.label == "Poospiza_lateralis")] <- "Microspingus_lateralis"
phylogeny$tip.label[which(phylogeny$tip.label == "Oreophylax_moreirae")] <- "Asthenes_moreirae"
phylogeny$tip.label[which(phylogeny$tip.label == "Thraupis_palmarum")] <- "Tangara_palmarum"
phylogeny$tip.label[which(phylogeny$tip.label == "Phaethornis_malaris")] <- "Phaethornis_margarettae"
phylogeny$tip.label[which(phylogeny$tip.label == "Aratinga_aurea")] <- "Eupsittula_aurea"
phylogeny$tip.label[which(phylogeny$tip.label == "Thraupis_sayaca")] <- "Tangara_sayaca"
phylogeny$tip.label[which(phylogeny$tip.label == "Tangara_velia")] <- "Tangara_cyanomelas"
phylogeny$tip.label[which(phylogeny$tip.label == "Pipra_pipra")] <- "Dixiphia_pipra"
phylogeny$tip.label[which(phylogeny$tip.label == "Pipra_rubrocapilla")] <- "Ceratopipra_rubrocapilla"
phylogeny$tip.label[which(phylogeny$tip.label == "Thryothorus_genibarbis")] <- "Pheugopedius_genibarbis"
phylogeny$tip.label[which(phylogeny$tip.label == "Picumnus_exilis")] <- "Picumnus_pernambucensis"
phylogeny$tip.label[which(phylogeny$tip.label == "Myrmeciza_loricata")] <- "Myrmoderus_loricatus"
phylogeny$tip.label[which(phylogeny$tip.label == "Basileuterus_flaveolus")] <- "Myiothlypis_flaveola"
phylogeny$tip.label[which(phylogeny$tip.label == "Oryzoborus_angolensis")] <- "Sporophila_angolensis"
phylogeny$tip.label[which(phylogeny$tip.label == "Parula_pitiayumi")] <- "Setophaga_pitiayumi"
phylogeny$tip.label[which(phylogeny$tip.label == "Basileuterus_leucoblepharus")] <- "Myiothlypis_leucoblephara"
phylogeny$tip.label[which(phylogeny$tip.label == "Poospiza_cabanisi")] <- "Microspingus_cabanisi"
phylogeny$tip.label[which(phylogeny$tip.label == "Carduelis_magellanica")] <- "Spinus_magellanicus"
phylogeny$tip.label[which(phylogeny$tip.label == "Basileuterus_leucophrys")] <- "Myiothlypis_leucophrys"
phylogeny$tip.label[which(phylogeny$tip.label == "Gallinula_chloropus")] <- "Gallinula_galeata"
phylogeny$tip.label[which(phylogeny$tip.label == "Buteo_magnirostris")] <- "Rupornis_magnirostris"
phylogeny$tip.label[which(phylogeny$tip.label == "Phaeothlypis_rivularis")] <- "Myiothlypis_rivularis"
phylogeny$tip.label[which(phylogeny$tip.label == "Philydor_lichtensteini")] <- "Anabacerthia_lichtensteini"
phylogeny$tip.label[which(phylogeny$tip.label == "Pyriglena_leuconota")] <- "Pyriglena_pernambucensis"
phylogeny$tip.label[which(phylogeny$tip.label == "Dendrocincla_fuliginosa")] <- "Dendrocincla_taunayi"
phylogeny$tip.label[which(phylogeny$tip.label == "Myrmeciza_ruficauda")] <- "Myrmoderus_ruficauda"
phylogeny$tip.label[which(phylogeny$tip.label == "Cercomacra_laeta")] <- "Cercomacroides_laeta"
phylogeny$tip.label[which(phylogeny$tip.label == "Ortalis_guttata")] <- "Ortalis_squamata"
phylogeny$tip.label[which(phylogeny$tip.label == "Clytolaema_rubricauda")] <- "Heliodoxa_rubricauda"
phylogeny$tip.label[which(phylogeny$tip.label == "Myrmotherula_gularis")] <- "Rhopias_gularis"
phylogeny$tip.label[which(phylogeny$tip.label == "Thraupis_ornata")] <- "Tangara_ornata"
phylogeny$tip.label[which(phylogeny$tip.label == "Tyto_alba")] <- "Tyto_furcata"
phylogeny$tip.label[which(phylogeny$tip.label == "Saltator_atricollis")] <- "Saltatricula_atricollis"
phylogeny$tip.label[which(phylogeny$tip.label == "Hylocryptus_rectirostris")] <- "Clibanornis_rectirostris"
phylogeny$tip.label[which(phylogeny$tip.label == "Pseudoscops_clamator")] <- "Asio_clamator"
phylogeny$tip.label[which(phylogeny$tip.label == "Campylopterus_largipennis")] <- "Campylopterus_calcirupicola"
phylogeny$tip.label[which(phylogeny$tip.label == "Cinclodes_pabsti")] <- "Cinclodes_espinhacensis"
phylogeny$tip.label[which(phylogeny$tip.label == "Cyanocompsa_brissonii")] <- "Cyanoloxia_brissonii"
phylogeny$tip.label[which(phylogeny$tip.label == "Aramides_cajanea")] <- "Aramides_cajaneus"
phylogeny$tip.label[which(phylogeny$tip.label == "Scytalopus_pachecoi")] <- "Scytalopus_petrophilus"
phylogeny$tip.label[which(phylogeny$tip.label == "Myrmeciza_squamosa")] <- "Myrmoderus_squamosus"
phylogeny$tip.label[which(phylogeny$tip.label == "Aratinga_leucophthalma")] <- "Psittacara_leucophthalmus"
phylogeny$tip.label[which(phylogeny$tip.label == "Agelaioides_badius")] <- "Agelaioides_fringillarius"
phylogeny$tip.label[which(phylogeny$tip.label == "Thryothorus_longirostris")] <- "Cantorchilus_longirostris"
phylogeny$tip.label[which(phylogeny$tip.label == "Caprimulgus_longirostris")] <- "Hydropsalis_longirostris"
phylogeny$tip.label[which(phylogeny$tip.label == "Caprimulgus_parvulus")] <- "Hydropsalis_parvula"

## Final check
phylogeny$tip.label[which(!phylogeny$tip.label %in% unique(data4$Binomial))]
cort <- phylogeny$tip.label[which(!phylogeny$tip.label %in% unique(data4$Binomial))]
phylogeny <- drop.tip(phylogeny, cort)
phylogeny$tip.label[which(!phylogeny$tip.label %in% unique(data4$Binomial))]
unique(data4$Binomial)[which(!unique(data4$Binomial) %in% phylogeny$tip.label)]
all(sort(unique(data4$Binomial)) ==  sort(phylogeny$tip.label))

## By correcting the variables, do we have any infinite values?
any(is.infinite(log(data4$BodyCond)))
any(is.infinite(log(data4$ForestCover)))

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
## habitat breadth: you can get here: https://esajournals.onlinelibrary.wiley.com/doi/10.1890/ES14-00332.1
hb_trau <- read.csv("./Habitat_breadth_index.csv", sep =";")
hb_trau$Species <- gsub("-", "_", hb_trau$Species)
hb_trau <- hb_trau[which(hb_trau$Species %in% data4$Binomial),]

## Checking presence
all(hb_trau$Species %in% unique(data4$Binomial))
all(unique(data4$Binomial) %in% hb_trau$Species)

## Removing absences
data5 <- data4[-which(!data4$Binomial %in% hb_trau$Species),]

## Checking again
all(hb_trau$Species %in% unique(data5$Binomial))
all(unique(data5$Binomial) %in% hb_trau$Species)

## PCA
## Variables
especies <- unique(data5$Binomial)
data5$HB1 <- NA
data5$HB2 <- NA
data5$HB3 <- NA
data5$HB4 <- NA

for(i in 1:length(especies)){
  data5[which(data5$Binomial == especies[i]),"HB1"] <- hb_trau[which(hb_trau$Species == especies[i]),"Habitat.breadth..co.occurrence.among.species.from.all.four.classes.multiplicative.beta."]
  data5[which(data5$Binomial == especies[i]),"HB2"] <- hb_trau[which(hb_trau$Species == especies[i]),"Habitat.breadth..Within.class.co.occurrence.index.additive.beta."]
  data5[which(data5$Binomial == especies[i]),"HB3"] <- hb_trau[which(hb_trau$Species == especies[i]),"Habitat.breadth..co.occurrence.among.species.from.all.four.classes.multiplicative.beta..1"]
  data5[which(data5$Binomial == especies[i]),"HB4"] <- hb_trau[which(hb_trau$Species == especies[i]),"Habitat.breadth..co.occurrence.among.species.from.all.four.classes.additive.beta."]
  
}

## PCA for habitat breadth
pca.hb <- prcomp(data5[,c("HB1","HB2","HB3","HB4")], scale = T)
summary(pca.hb)
biplot(pca.hb)
fviz_eig(pca.hb)
cor.clima.eixos <- pca.hb$rotation
eixos.pca <- pca.hb$x
data5$PCHB1 <-  eixos.pca[,1]

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
## Cluster analysis
## Creating single row for each species
matriz <- as.data.frame(matrix(nrow = length(unique(data5$Binomial)), ncol = 7, NA))
row.names(matriz) <- unique(data5$Binomial)
colnames(matriz) <- c("Order","Family","Forest_cover","Body_condition", "Habitat_breadth", "Body_mass", "Body_length")

unicas <- unique(data5$Binomial)
for(i in 1:dim(matriz)[1]){
  ## ID
  matriz[which(row.names(matriz) == unicas[i]),"Order"] <- unique(data5[which(data5$Binomial == unicas[i]),"Order"])
  matriz[which(row.names(matriz) == unicas[i]),"Family"] <- unique(data5[which(data5$Binomial == unicas[i]),"Family"])
  
  ## Data
  matriz[which(row.names(matriz) == unicas[i]),"Forest_cover"] <- round(mean(data5[which(data5$Binomial == unicas[i]),"ForestCover"]),2)
  matriz[which(row.names(matriz) == unicas[i]),"Body_condition"] <- round(mean(sqrt(data5[which(data5$Binomial == unicas[i]),"BodyCond"])),2)
  matriz[which(row.names(matriz) == unicas[i]),"Body_mass"] <- round(mean(sqrt(data5[which(data5$Binomial == unicas[i]),"Body_mass.g."])),2)
  matriz[which(row.names(matriz) == unicas[i]),"Body_length"] <- round(mean(sqrt(data5[which(data5$Binomial == unicas[i]),"Body_length.mm."])),2)
  matriz[which(row.names(matriz) == unicas[i]),"Habitat_breadth"] <- unique(data5[which(data5$Binomial == unicas[i]),"PCHB1"])
  
}

#write.csv(matriz, "./Dados.csv")
matriz2 <- matriz
matriz2 <- matriz2[,-c(1,2,3,5)]
## How many clusters shoulw we use?
fviz_nbclust(matriz2, kmeans, method = "silhouette")

## Cluster analysis
k2 <- kmeans(matriz2, centers = 2, nstart = 25)
str(k2)
fviz_cluster(k2, data = matriz2, repel = T, ggtheme = theme_bw())

## Accurary
regular_scaling = train(log(Body_condition) ~ log(Body_mass)  + normalize(Body_length), data = matriz, method = "knn", 
                        trControl = trainControl(method = "cv", number = 10), 
                        preProcess = c("center", "scale"), 
                        tuneGrid = expand.grid(k = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)))
plot(regular_scaling)
summary(regular_scaling)
max(regular_scaling$results$RMSE)

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
## Getting data from cluster
data5$Cluster <- NA
unicas <- unique(data5$Binomial)
grupos <- as.data.frame(k2$cluster)
grupos$Species <- row.names(grupos)
row.names(grupos) <- NULL

## looping through
for(i in 1:dim(data5)[1]){
  data5[which(data5$Binomial == unicas[i]),"Cluster"] <- grupos[which(grupos$Species == unicas[i]),1]
}

## Final check
any(is.na(data5$Cluster))

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
## Climate
df_final <- data5[,c(5,6,9,35,37,44,59,60,62,73,76,84,85,90,91)]
colnames(df_final)
str(df_final)

## Get climate data
r <- raster::getData("worldclim",var="bio",res=10)

## Temp and rainfall
occurrences <- unique(as.data.frame.matrix(df_final[,c("Longitude_decimal_degrees", "Latitude_decimal_degrees")]))
coords <- data.frame(x=occurrences$Longitude_decimal_degrees, y=occurrences$Latitude_decimal_degrees)
points <- SpatialPoints(coords, proj4string = r@crs)
values <- extract(r,points)
clima <- as.data.frame(cbind(occurrences$Longitude_decimal_degrees,occurrences$Latitude_decimal_degrees,values))
colnames(clima)[c(1,2)] <- c("Longitude", "Latitude")

## Transfering data
df_final$Temp <- NA
df_final$Rainfall <- NA

## Getting data
for(i in 1:dim(df_final)[1]){
  ## Filter
  filtro <- df_final[i,c("Longitude_decimal_degrees", "Latitude_decimal_degrees")]
  result <- c(filtro[1,1], filtro[1,2])
  
  ## Getting data
  value <- clima[which(clima[,c("Longitude")] == result[1] & clima[,c("Latitude")] == result[2]),"bio1"]
  df_final[i,"Temp"] <- value
  
  value2 <- clima[which(clima[,c("Longitude")] == result[1] & clima[,c("Latitude")] == result[2]),"bio12"]
  df_final[i,"Rainfall"] <- value2

}

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
## Bird diversity: you can get the data here: https://www.globio.info/globio-data-downloads: MSA tif birds and mammals 2015
diversity_map <- raster("./TerrestrialMSA_2015_World_wbvert.tif")
plot(diversity_map)

## Getting habitat cover
occurrence_3 <- unique(as.data.frame.matrix(df_final[,c("Longitude_decimal_degrees", "Latitude_decimal_degrees")]))
local_diversity <- raster::extract(diversity_map, occurrence_3)
diversity_species <- cbind(occurrence_3, local_diversity)
diversity_species[which(is.na(diversity_species$local_diversity)),]

## Some didn't appear
problems_2 <- diversity_species[which(is.na(diversity_species$local_diversity)),]
problems_2$Longitude_decimal_degrees <- as.character(problems_2$Longitude_decimal_degrees)
problems_2$Latitude_decimal_degrees <- as.character(problems_2$Latitude_decimal_degrees)

## Check why
diversity_species$Longitude_decimal_degrees <- as.character(diversity_species$Longitude_decimal_degrees)
diversity_species$Latitude_decimal_degrees <- as.character(diversity_species$Latitude_decimal_degrees)

## Getting closest value
for(i in 1:dim(problems_2)[1]){
  correct___2 <- diversity_species[grep(substr(problems_2$Longitude_decimal_degrees[i], 1,3), diversity_species$Longitude_decimal_degrees),]
  correct___2_2 <- correct___2[-which(is.na(correct___2$local_diversity)),]
  avg_2 <- mean(as.numeric(correct___2_2$local_diversity))
  
  problems_2$local_diversity[i] <- avg_2
}

## Bringing to final data frame
diversity_species[which(diversity_species$Longitude_decimal_degrees %in% problems_2$Longitude_decimal_degrees & 
                          diversity_species$Latitude_decimal_degrees %in% problems_2$Latitude_decimal_degrees & 
                          is.na(diversity_species$local_diversity)==T ),"local_diversity"] <- problems_2$local_diversity

## Bringing data
df_final$Speciesdiversity <- NA
## Getting data
for(i in 1:dim(df_final)[1]){
  filtro <- df_final[i,c("Longitude_decimal_degrees", "Latitude_decimal_degrees")]
  result <- c(filtro[1,1], filtro[1,2])
  value <- diversity_species[which(diversity_species[,c("Longitude_decimal_degrees")] == result[1] & diversity_species[,c("Latitude_decimal_degrees")] == result[2]),"local_diversity"]
  df_final[i,"Speciesdiversity"] <- value
}

str(df_final)
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
df_final <- read.csv("./Planilha_BD_final3.csv")
colnames(df_final)

## NA presence
any(is.na(df_final$Binomial))
any(is.na(df_final$Age))
any(is.na(df_final$Sex))
any(is.na(df_final$Molt))
any(is.na(df_final$Latitude_decimal_degrees))
any(is.na(df_final$State))
any(is.na(df_final$Elevation))
any(is.na(df_final$Temp))
any(is.na(df_final$Rainfall))
any(is.na(df_final$Year))
any(is.na(df_final$ForestCover))
any(is.na(df_final$BodyCond))
any(is.na(df_final$PCHB1))
any(is.na(df_final$Cluster))
any(is.na(df_final$Speciesdiversity))

## Distributions
hist(df_final$Temp)
hist(log(df_final$Rainfall))
hist(normalize(df_final$Elevation))
hist(normalize(df_final$Speciesdiversity))
hist(sqrt(df_final$ForestCover))
hist(log(df_final$BodyCond))
hist(scale(df_final$PCHB1))

colnames(df_final)
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
## Subtropical x tropical
table(round(df_final$Latitude_decimal_degrees,0))

## Creating variable
df_final$Region <- NA
df_final[which(round(df_final$Latitude_decimal_degrees,2) > -23),"Region"] <- "Tropical"
df_final[which(is.na(df_final$Region)),"Region"] <- "Sub-tropical"

## Checking
any(is.na(df_final$Region))
table(df_final$Region)
df_final$Region <- as.factor(df_final$Region)

## Final datasheet
unicas <- unique(df_final$Binomial)
for(i in 1:length(unicas)){
  if(length(table(df_final[which(df_final$Binomial == unicas[i]),"Region"])) == 2){
    print(unicas[i])
    print(table(df_final[which(df_final$Binomial == unicas[i]),"Region"]))
    print("##################")
  }
}

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
## Cliuster analysis
## Getting final data
df_final <- read.csv("./Planilha_BD_final2.csv")
str(df_final)

## Phylogeny
phylogeny$tip.label[which(!phylogeny$tip.label %in% unique((df_final)$Binomial))]
cort <- phylogeny$tip.label[which(!phylogeny$tip.label %in% unique((df_final)$Binomial))]
phylogeny <- drop.tip(phylogeny, cort)

## Creating new data frame
matriz <- as.data.frame(matrix(nrow = length(unique(df_final$Binomial)), ncol = 10, NA))
colnames(matriz) <- c("Species","BodyCond", "ForestCover", "PCHB1", "Rainfall", "SpeciesDiversity", "cluster", "Temp", "Elevation", "Region")
matriz$Species <- unique(df_final$Binomial)

## Making averages per species
colnames(df_final)
unicas <- unique(df_final$Binomial)
for(i in 1:length(unicas)){
  matriz[which(matriz$Species == unicas[i]),"BodyCond"] <- mean(df_final[which(df_final$Binomial == unicas[i]),"BodyCond"])
  matriz[which(matriz$Species == unicas[i]),"ForestCover"] <- mean(df_final[which(df_final$Binomial == unicas[i]),"ForestCover"])
  matriz[which(matriz$Species == unicas[i]),"PCHB1"] <- mean(df_final[which(df_final$Binomial == unicas[i]),"PCHB1"])
  matriz[which(matriz$Species == unicas[i]),"Rainfall"] <- mean(df_final[which(df_final$Binomial == unicas[i]),"Rainfall"])
  matriz[which(matriz$Species == unicas[i]),"SpeciesDiversity"] <- mean(df_final[which(df_final$Binomial == unicas[i]),"Speciesdiversity"])
  matriz[which(matriz$Species == unicas[i]),"Temp"] <- mean(df_final[which(df_final$Binomial == unicas[i]),"Temp"])
  matriz[which(matriz$Species == unicas[i]),"Elevation"] <- mean(df_final[which(df_final$Binomial == unicas[i]),"Elevation"])
  matriz[which(matriz$Species == unicas[i]),"cluster"] <- mean(df_final[which(df_final$Binomial == unicas[i]),"Cluster"])
  
  ## Region
  if(length(unique(df_final[which(df_final$Binomial == unicas[i]),"Region"])) > 1){
    matriz[which(matriz$Species == unicas[i]),"Region"] <- "Both"
  }
  else{
    matriz[which(matriz$Species == unicas[i]),"Region"] <- unique(df_final[which(df_final$Binomial == unicas[i]),"Region"])
  }
}

## Checking
table(matriz$Region)

## Checking
all(matriz$Species %in% phylogeny$tip.label)
all(phylogeny$tip.label %in% matriz$Species)

## Species in cluster 1
as.data.frame(matriz[which(matriz$cluster == 1), "Species"])

## GLS for cluster: models
model_cluster_BM <- gls(scale(sqrt(BodyCond)) ~ as.factor(cluster), data = matriz,
                     correlation = corBrownian(phy = phylogeny, form = ~Species), method = "ML")
model_cluster_OU <- gls(scale(sqrt(BodyCond)) ~ as.factor(cluster), data = matriz,
                     correlation = corPagel(1,phy = phylogeny, fixed = F,form = ~Species), method = "ML")

## AIC comparision
AIC(model_cluster_BM, model_cluster_OU)

## Null mmodel
model_cluster_OU_null <- gls(scale(sqrt(BodyCond)) ~ 1, data = matriz,
                        correlation = corPagel(1,phy = phylogeny, fixed = F,form = ~Species), method = "ML")

## variable significance
drop1(model_cluster_OU, test = "Chisq")

## P-value
summary(model_cluster)

## Plotting results
jpeg("./Cluster.jpg", width = 4, height = 4, units = 'in', res = 800) 
ggplot(matriz, aes(x = as.factor(as.character(cluster)), y = BodyCond)) + geom_boxplot() + xlab("Cluster identity") + ylab("SMI") + theme_bw()
dev.off()

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
## Phylogenetic models of body condition without the cluster 1
df_final_2 <- df_final[-which(df_final$Cluster == 1), ]
cortar_2 <- which(!phylogeny$tip.label %in% df_final_2$Binomial)
phylogeny_2 <- drop.tip(phylogeny, cortar_2)

## Histogram
hist(log(df_final_2$BodyCond))
hist(log(df_final_2$ForestCover))
hist(normalize(df_final_2$PCHB1))
hist(df_final_2$Rainfall)
hist(df_final_2$Speciesdiversity)
hist(df_final_2$Temp)
hist(sqrt(df_final_2$Elevation))

## Models
## Multicolinearity 
colnames(df_final)
model_full_collinearity <- lm(scale(log(BodyCond)) ~ scale(sqrt(ForestCover)) + scale(PCHB1) + scale(Temp) + Region + 
                                scale(log(Rainfall)) + Region + scale(Speciesdiversity), data = df_final_2, na.action = na.fail)
VIF(model_full_collinearity)

## Full model
model_full_2 <- phyr::pglmm(scale(log(BodyCond)) ~ scale(sqrt(ForestCover)) + scale(normalize(PCHB1)) + scale(Temp) + 
                            scale(log(Rainfall)) + scale(normalize(Speciesdiversity)) + 
                              (1|Age) + (1|Sex) + (1|Molt) + (1|Year) + (1|State), 
                            data = df_final_2, cov_ranef = list(sp = phylogeny_2), REML = TRUE, verbose = TRUE)
## Null model
model_null <- phyr::pglmm(scale(log(BodyCond)) ~ 1 +  
                            (1|Age) + (1|Sex) + (1|Molt) + (1|Year) + (1|State), data = df_final_2,
                    cov_ranef = list(sp = phylogeny_2), REML = TRUE, verbose = TRUE)


## Model comparison
pvalue <- pchisq(2*(model_full$logLik - model_null$logLik), df = 4, lower.tail = FALSE)
round(pvalue,10)

## Full was significantly different
summary(model_full_2)

## Predicted values
df_final_2$Predicted <- NA
df_final_2$Predicted <- as.data.frame(predict(model_full_2))[,1]

## Plotting results
ggplot(df_final_2, aes(x = scale(normalize(Speciesdiversity)), y =Predicted)) + geom_point() +  geom_smooth(method = "glm", color = "black", se = T) +  
  theme_bw() + ylab("Predicted values of SMI") + xlab("Local species diversity")
ggplot(df_final_2, aes(x =  scale(normalize(PCHB1)), y = Predicted)) + geom_point() +  geom_smooth(method = "glm", color = "black", se = T) +  
  theme_bw() + ylab("Predicted values of SMI") + xlab("Habitat breadht")

## End of code
##############################################################################################################################################
##############################################################################################################################################