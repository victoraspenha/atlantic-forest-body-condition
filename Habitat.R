##################################################################################################################################################################################################
##################################################################################################################################################################################################
## Removing leftovers
rm(list=ls())

## Needed packages
library(arm)
library(fossil)
library(dismo)
library(rgeos)
library(geosphere)
library(sf)
library(rgdal)
library(raster)
library(sp)
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
library(geosphere)
library(rnaturalearth)
library(ggspatial)
library(forestplot)
require(smatr)
require(magrittr)
require(MASS)
require(data.table)

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
## Getting data -- you can find the data here: https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2647
data <- read.csv("./Documents/Manuscripts/0 - Habitat/Planilhas/ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv", h=T)
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
mapa <- raster("./Documents/Manuscripts/0 - Habitat/Dissimilaridade/Dissimilarity_01_05_1km_uint32.tif")

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

scaledMassIndex <-
  function(x, y, x.0 = mean(x)) {
    logM.ols <- lm(log(y) ~ log(x))
    logM.rob <- rlm(log(y) ~ log(x), method = "M")
    b.msa.ols <- coef(sma(log(y) ~ log(x)))[2]
    b.msa.rob <- coef(sma(log(y) ~ log(x), robust = T))[2]
    SMI.ols <- y * (x.0 / x) ^ b.msa.ols
    SMI.rob <- y * (x.0 / x) ^ b.msa.rob
    res <- data.frame(SMI.ols, SMI.rob, x, y)
    pred.DT <-
      data.table(x = seq(min(x), max(x), length = 100)) %>%
      .[, y.ols := predict(logM.ols, newdata = .) %>% exp] %>%
      .[, y.rob := predict(logM.rob, newdata = .) %>% exp]
    attr(res, "b.msa") <- c(ols = b.msa.ols, rob = b.msa.rob)
    return(res)
  }

## Body condition
data4$BodyCond <- NA
especies <- unique(data4$Binomial)
for(i in 1:length(especies)){
  print(especies[i])
  print("################################")
  dentro <- data4[which(data4$Binomial == especies[i]),]
  locais <- unique(cbind(dentro$Latitude_decimal_degrees, dentro$Longitude_decimal_degrees))
  
  for(j in 1:dim(locais)[1]){
    body_mass <- dentro[which(dentro$Latitude_decimal_degrees == locais[j,1] & dentro$Longitude_decimal_degrees == locais[j,2]),"Body_mass.g."]
    body_ltngh <- dentro[which(dentro$Latitude_decimal_degrees == locais[j,1] & dentro$Longitude_decimal_degrees == locais[j,2]),"Body_length.mm."]

    if((length(body_mass) <= 5) == TRUE){
      data4[data4$Binomial == especies[i] & data4$Latitude_decimal_degrees == locais[j,1] & data4$Longitude_decimal_degrees == locais[j,2],"BodyCond"] <- "NA"
      print("NOT ENOUGH DATA")
    }
    else{
      tryCatch(
        expr = {
          residuos <- round(scaledMassIndex(body_mass, body_ltngh)[,1],2)
          data4[which(data4$Binomial == especies[i] & data4$Latitude_decimal_degrees == locais[j,1] & data4$Longitude_decimal_degrees == locais[j,2]),"BodyCond"] <- residuos
          message("Successfully executed..")
        },
        error = function(e){
          message('Caught an error!')
          data4[data4$Binomial == especies[i] & data4$Latitude_decimal_degrees == locais[j,1] & data4$Longitude_decimal_degrees == locais[j,2],"BodyCond"] <- "Error"
          print(e)
        },
        warning = function(w){
          data4[data4$Binomial == especies[i] & data4$Latitude_decimal_degrees == locais[j,1] & data4$Longitude_decimal_degrees == locais[j,2],"BodyCond"] <- "Error"
          message('Caught an warning!')
          print(w)
        },
        finally = {
          message('All done, quitting.')
        }
      )
    }
  }
}

which(data4$BodyCond == "Error")
## None
which(is.na(data4$BodyCond))
which(data4$BodyCond == "NA")

## Removing the ones that we did not have enough data
data5 <- data4[-which(data4$BodyCond == "NA"),]
data5 <- data5[-which(is.na(data5$BodyCond)),]
data5$BodyCond <- as.numeric(data5$BodyCond)

## Variable distribution
hist(log(data5$BodyCond))

## Checking cofounders
any(is.na(data5$Year))
any(is.na(data5$Sex))
any(is.na(data4$Age))
any(is.na(data4$Molt))

data5 <- data5[-which(is.na(data5$Year)),]
data5 <- data5[-which(is.na(data5$Sex)),]
data5 <- data5[-which(is.na(data5$Age)),]
data5 <- data5[-which(is.na(data5$Molt)),]

## Remove unknown
table(data5$Year)
table(data5$Sex)
table(data5$Age)
table(data4$Molt)
data5 <- data5[-which(data5$Age == "Unknown"),]
data5 <- data5[-which(data5$Sex == "Unknown"),]

## Final Check
table(data5$Year)
table(data5$Sex)
table(data5$Age)
table(data5$Molt)

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
dim(data5)
## Phylogeny: you can get phylogeny here: https://birdtree.org/subsets/
phylogeny <- read.nexus("./Documents/Manuscripts/0 - Habitat/Filogenias/output.nex")
phylogeny <- maxCladeCred(phylogeny)

## Checking for similarities between phylogeny and dataset
data5$Binomial <- gsub(" ", "_", data5$Binomial)
phylogeny$tip.label[which(!phylogeny$tip.label %in% unique(data5$Binomial))]
unique(data5$Binomial)[which(!unique(data5$Binomial) %in% phylogeny$tip.label)]

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
phylogeny$tip.label[which(!phylogeny$tip.label %in% unique(data5$Binomial))]
cort <- phylogeny$tip.label[which(!phylogeny$tip.label %in% unique(data5$Binomial))]
phylogeny <- drop.tip(phylogeny, cort)
phylogeny$tip.label[which(!phylogeny$tip.label %in% unique(data5$Binomial))]
unique(data5$Binomial)[which(!unique(data5$Binomial) %in% phylogeny$tip.label)]
all(sort(unique(data5$Binomial)) ==  sort(phylogeny$tip.label))

## By correcting the variables, do we have any infinite values?
any(is.infinite(log(data5$BodyCond)))
any(is.infinite(log(data5$ForestCover)))

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
## habitat breadth: you can get here: https://esajournals.onlinelibrary.wiley.com/doi/10.1890/ES14-00332.1
hb_trau <- read.csv("./Documents/Manuscripts/0 - Habitat/Planilhas//Habitat_breadth_index.csv", sep =";")
hb_trau$Species <- gsub("-", "_", hb_trau$Species)
hb_trau <- hb_trau[which(hb_trau$Species %in% data5$Binomial),]

## Checking presence
all(hb_trau$Species %in% unique(data5$Binomial))
all(unique(data5$Binomial) %in% hb_trau$Species)

## Removing absences
data6 <- data5[-which(!data5$Binomial %in% hb_trau$Species),]

## Checking again
all(hb_trau$Species %in% unique(data6$Binomial))
all(unique(data6$Binomial) %in% hb_trau$Species)

## PCA
## Variables
especies <- unique(data6$Binomial)
data6$HB1 <- NA
data6$HB2 <- NA
data6$HB3 <- NA
data6$HB4 <- NA

for(i in 1:length(especies)){
  data6[which(data6$Binomial == especies[i]),"HB1"] <- hb_trau[which(hb_trau$Species == especies[i]),"Habitat.breadth..co.occurrence.among.species.from.all.four.classes.multiplicative.beta."]
  data6[which(data6$Binomial == especies[i]),"HB2"] <- hb_trau[which(hb_trau$Species == especies[i]),"Habitat.breadth..Within.class.co.occurrence.index.additive.beta."]
  data6[which(data6$Binomial == especies[i]),"HB3"] <- hb_trau[which(hb_trau$Species == especies[i]),"Habitat.breadth..co.occurrence.among.species.from.all.four.classes.multiplicative.beta..1"]
  data6[which(data6$Binomial == especies[i]),"HB4"] <- hb_trau[which(hb_trau$Species == especies[i]),"Habitat.breadth..co.occurrence.among.species.from.all.four.classes.additive.beta."]
  
}

## PCA for habitat breadth
pca.hb <- prcomp(data6[,c("HB1","HB2","HB3","HB4")], scale = T)
summary(pca.hb)
biplot(pca.hb)
fviz_eig(pca.hb)
cor.clima.eixos <- pca.hb$rotation
eixos.pca <- pca.hb$x
data6$PCHB1 <-  eixos.pca[,1]

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
matriz <- as.data.frame(matrix(nrow = length(unique(data6$Binomial)), ncol = 7, NA))
row.names(matriz) <- unique(data6$Binomial)
colnames(matriz) <- c("Order","Family","Forest_cover","Body_condition", "Habitat_breadth", "Body_mass", "Body_length")

unicas <- unique(data6$Binomial)
for(i in 1:dim(matriz)[1]){
  ## ID
  matriz[which(row.names(matriz) == unicas[i]),"Order"] <- unique(data6[which(data6$Binomial == unicas[i]),"Order"])
  matriz[which(row.names(matriz) == unicas[i]),"Family"] <- unique(data6[which(data6$Binomial == unicas[i]),"Family"])
  
  ## Data
  matriz[which(row.names(matriz) == unicas[i]),"Forest_cover"] <- round(mean(data6[which(data6$Binomial == unicas[i]),"ForestCover"]),2)
  matriz[which(row.names(matriz) == unicas[i]),"Body_condition"] <- round(mean(sqrt(data6[which(data6$Binomial == unicas[i]),"BodyCond"])),2)
  matriz[which(row.names(matriz) == unicas[i]),"Body_mass"] <- round(mean(sqrt(data6[which(data6$Binomial == unicas[i]),"Body_mass.g."])),2)
  matriz[which(row.names(matriz) == unicas[i]),"Body_length"] <- round(mean(sqrt(data6[which(data6$Binomial == unicas[i]),"Body_length.mm."])),2)
  matriz[which(row.names(matriz) == unicas[i]),"Habitat_breadth"] <- unique(data6[which(data6$Binomial == unicas[i]),"PCHB1"])
  
}

#write.csv(matriz, "./Dados.csv")
matriz2 <- matriz
matriz2 <- matriz2[,-c(1,2,3,5,7)]
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
data6$Cluster <- NA
unicas <- unique(data6$Binomial)
grupos <- as.data.frame(k2$cluster)
grupos$Species <- row.names(grupos)
row.names(grupos) <- NULL

## looping through
for(i in 1:dim(data6)[1]){
  data6[which(data6$Binomial == unicas[i]),"Cluster"] <- grupos[which(grupos$Species == unicas[i]),1]
}

## Final check
any(is.na(data6$Cluster))
table(data6$Cluster)
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
## Climate
colnames(data6)
df_final <- data6
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
diversity_map <- raster("./Documents/Manuscripts/0 - Habitat/TerrestrialMSA_2015_World_wbvert.tif")
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
df_final$Speciesdiversity
str(df_final)
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
## df_final <- read.csv("./Documents/Manuscripts/0 - Habitat/Data/Planilha_BD_final3.csv")
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
## Getting final data
## df_final <- read.csv("./Documents/Manuscripts/0 - Habitat/Data/Planilha_BD_final2.csv")
str(df_final)

## Remove RS
df_final <- df_final[-which(df_final$State == "RS"),]
## df_final <- df_final[-which(df_final$State == "PE" & df_final$Longitude_decimal_degrees == -39.656944),]

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
## Habitat
## Getting shapefile
## Source: https://ecoregions.appspot.com/
shape1 <- raster("./Documents/Manuscripts/0 - Habitat/lvl1_frac_1km_ver003/lvl1_frac_1km_ver003/iucn_habitatclassification_fraction_lvl1__100_Forest__ver003.tif")
shape2 <- raster("./Documents/Manuscripts/0 - Habitat/lvl1_frac_1km_ver003/lvl1_frac_1km_ver003/iucn_habitatclassification_fraction_lvl1__1400_Artificial - Terrestrial__ver003.tif")
shape3 <- raster("./Documents/Manuscripts/0 - Habitat/lvl1_frac_1km_ver003/lvl1_frac_1km_ver003/iucn_habitatclassification_fraction_lvl1__200_Savanna__ver003.tif")
shape4 <- raster("./Documents/Manuscripts/0 - Habitat/lvl1_frac_1km_ver003/lvl1_frac_1km_ver003/iucn_habitatclassification_fraction_lvl1__300_Shrubland__ver003.tif")
shape5 <- raster("./Documents/Manuscripts/0 - Habitat/lvl1_frac_1km_ver003/lvl1_frac_1km_ver003/iucn_habitatclassification_fraction_lvl1__400_Grassland__ver003.tif")
shape6 <- raster("./Documents/Manuscripts/0 - Habitat/lvl1_frac_1km_ver003/lvl1_frac_1km_ver003/iucn_habitatclassification_fraction_lvl1__500_Wetlands inland__ver003.tif")
shape7 <- raster("./Documents/Manuscripts/0 - Habitat/lvl1_frac_1km_ver003/lvl1_frac_1km_ver003/iucn_habitatclassification_fraction_lvl1__600_Rocky Areas__ver003.tif")
shape8 <- raster("./Documents/Manuscripts/0 - Habitat/lvl1_frac_1km_ver003/lvl1_frac_1km_ver003/iucn_habitatclassification_fraction_lvl1__800_Desert__ver003.tif")

occurrence_4 <- as.data.frame.matrix(df_final[,c("Longitude_decimal_degrees", "Latitude_decimal_degrees")])

local_habitat1 <- raster::extract(shape1, occurrence_4)
local_habitat2 <- raster::extract(shape2, occurrence_4)
local_habitat3 <- raster::extract(shape3, occurrence_4)
local_habitat4 <- raster::extract(shape4, occurrence_4)
local_habitat5 <- raster::extract(shape5, occurrence_4)
local_habitat6 <- raster::extract(shape6, occurrence_4)
local_habitat7 <- raster::extract(shape7, occurrence_4)
local_habitat8 <- raster::extract(shape8, occurrence_4)

habitat_species <- cbind(occurrence_4, local_habitat1, local_habitat2, local_habitat3, local_habitat4, local_habitat5, local_habitat6, local_habitat7, local_habitat8)
habitat_species$Habitat <- NA

## Loop through all habitat types and choose the best one
for(i in 1:dim(habitat_species)[1]){
  filtro <- habitat_species[i,-c(1,2)]
  maior <- names(which.max(filtro))
  
  if(is.null(maior)){
    habitat_species$Habitat[i] <- NA
  }
  else{
    if(maior == "local_habitat1"){
      habitat_species$Habitat[i] <- "Forest"
    }
    else{
      if(maior == "local_habitat2"){
        habitat_species$Habitat[i] <- "Artifical_land"
      }
      else{
        if(maior == "local_habitat3"){
          habitat_species$Habitat[i] <- "Savanna"
        }
        else{
          if(maior == "local_habitat4"){
            habitat_species$Habitat[i] <- "Shrubland"
          }
          else{
            if(maior == "local_habitat5"){
              habitat_species$Habitat[i] <- "Grassland"
            }
            else{
              if(maior == "local_habitat6"){
                habitat_species$Habitat[i] <- "Wetlands"
              }
              else{
                if(maior == "local_habitat7"){
                  habitat_species$Habitat[i] <- "Rocky_areas"
                }
                else{
                  if(aior == "local_habitat8"){
                    habitat_species$Habitat[i] <- "Desert"
                  }
                  else{
                    print("NADA")
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
}

df_final$Habitat <- NA

for(i in 1:dim(df_final)[1]){
  df_final[which(df_final$Latitude_decimal_degrees[i] == habitat_species$Latitude_decimal_degrees & 
                   df_final$Longitude_decimal_degrees[i] == habitat_species$Longitude_decimal_degrees),"Habitat"] <- 
    habitat_species[which(habitat_species$Latitude_decimal_degrees == df_final$Latitude_decimal_degrees[i] & 
                            habitat_species$Longitude_decimal_degrees == df_final$Longitude_decimal_degrees[i]),"Habitat"]
}

table(df_final$Habitat)

#############################################################################################################################################
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
## Diet
dieta <- read.csv("./Documents/Manuscripts/0 - Habitat/BirdFuncDat.csv")
dieta$Scientific <- gsub(" ", "_", dieta$Scientific)
dieta_2 <- dieta[which(dieta$Scientific %in% df_final$Binomial),]

df_final$Diet <- NA
unicas <- unique(df_final$Binomial)
for(i in 1:length(unicas)){
  if(length(which(dieta_2$Scientific == unicas[i])) == 0){
    df_final[which(df_final$Binomial == unicas[i]),"Diet"] <- NA
  }
  else{
    df_final[which(df_final$Binomial == unicas[i]),"Diet"] <- dieta_2[which(dieta_2$Scientific == unicas[i]),"Diet.5Cat"]
  }
  
}

df_final[which(!complete.cases(df_final$Diet)),"Diet"] <- "Omnivore"
table(df_final$Diet)

## Export dataframe so far
#write.csv(df_final, "./Documents/Manuscripts/0 - Habitat/Final.csv")
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
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
## https://gis.stackexchange.com/questions/17638/clustering-spatial-data-in-r
## Density
## df_final <- read.csv("./Documents/Manuscripts/0 - Habitat/Final.csv")
colnames(df_final)

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
## Final check
phylogeny$tip.label[which(!phylogeny$tip.label %in% unique(df_final$Binomial))]
cort <- phylogeny$tip.label[which(!phylogeny$tip.label %in% unique(df_final$Binomial))]
phylogeny <- drop.tip(phylogeny, cort)
phylogeny$tip.label[which(!phylogeny$tip.label %in% unique(df_final$Binomial))]
unique(df_final$Binomial)[which(!unique(df_final$Binomial) %in% phylogeny$tip.label)]
all(sort(unique(df_final$Binomial)) ==  sort(phylogeny$tip.label))

## Clster analysis
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

table(df_final$Cluster)
## Checking
table(matriz$Region)

## Checking
all(matriz$Species %in% phylogeny$tip.label)
all(phylogeny$tip.label %in% matriz$Species)
table(matriz$cluster)

## Species in cluster 1
as.data.frame(matriz[which(matriz$cluster == 1),])
dim(data6[which(data6$Binomial == "Geotrygon_montana"),])
dim(data6[which(data6$Binomial == "Vanellus_chilensis"),])
dim(data6[which(data6$Binomial == "Celeus_flavescens"),])
dim(data6[which(data6$Binomial == "Colaptes_melanochloros"),])
dim(data6[which(data6$Binomial == "Selenidera_maculirostris"),])
dim(data6[which(data6$Binomial == "Trogon_viridis"),])

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
summary(model_cluster_OU)

## Plotting results
jpeg("./Cluster.jpg", width = 4, height = 4, units = 'in', res = 800) 
ggplot(matriz, aes(x = as.factor(as.character(cluster)), y = BodyCond)) + geom_boxplot() + xlab("Cluster identity") + ylab("SMI") + theme_bw()
dev.off()


######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
## Removing Cluster 1
df_final <- df_final[-which(df_final$Cluster == 1), ]
cortar_2 <- which(!phylogeny$tip.label %in% df_final$Binomial)
phylogeny <- drop.tip(phylogeny, cortar_2)

## Histogram
hist(log(df_final$BodyCond))
hist(log(df_final$ForestCover))
hist(normalize(df_final$PCHB1))
hist(df_final$Rainfall)
hist(df_final$Speciesdiversity)
hist(df_final$Temp)

## Remove nesltigls
table(df_final$Age)
table(df_final$Sex)
table(df_final$Molt)
table(df_final$Year)
table(df_final$Latitude_decimal_degrees)
table(df_final$Diet)

## write.csv(df_final_3, "./Documents/Final.csv")
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
## Maps
## Shapefiles: https://tapiquen-sig.jimdo.com/english-version/free-downloads/europe/
nordic <- ne_countries(scale = "medium", returnclass = "sf", continent = "South America", country = "Brazil")
## nordic <- ne_states(country = "Brazil", returnclass = "sf")
world_points<- st_centroid(nordic)
world_points <- cbind(nordic, st_coordinates(st_centroid(nordic$geometry)))

jpeg("Map.jpeg", width = 8, height = 8, units = 'in', res = 800) 
ggplot() +
  geom_sf(data = nordic, fill= "white") +
  #coord_sf(xlim = c(10, 30), ylim = c(60, 71)) + 
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "black", fontface = "bold", check_overlap = TRUE) +
  ylab("Latitude") + xlab("Longitude") + 
  annotation_scale(location = "bl", width_hint = 0.5) + 
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "white")) + 
  geom_point(aes(x = df_final$Longitude_decimal_degrees, y = df_final$Latitude_decimal_degrees), color="black", size = 2)
dev.off()

table(df_final$Binomial)
## Models
## Multicolinearity 
colnames(df_final)
model_full_collinearity <- lm(scale(log(BodyCond)) ~ scale(sqrt(ForestCover)) + scale(PCHB1) + scale(Temp) + 
                                scale(log(Rainfall)) + scale(Speciesdiversity)+Region, data = df_final, na.action = na.fail)
VIF(model_full_collinearity)

## Remove region
model_full_collinearity_2 <- lm(scale(log(BodyCond)) ~ scale(sqrt(ForestCover)) + scale(PCHB1) + scale(Temp) + 
                                scale(log(Rainfall)) + scale(Speciesdiversity), data = df_final, na.action = na.fail)
VIF(model_full_collinearity_2)

## Full model
model_full_2 <- phyr::pglmm(scale(log(BodyCond)) ~ scale(sqrt(ForestCover)) + scale(normalize(PCHB1)) + scale(Temp) + 
                            scale(log(Rainfall)) + scale(normalize(Speciesdiversity)) + 
                              (1|Age) + (1|Sex) + (1|Molt) + (1|Year) + (1|Latitude_decimal_degrees)+
                              (1|Diet), 
                            data = df_final, cov_ranef = list(sp = phylogeny), REML = TRUE, verbose = TRUE)

## Null model
model_null <- phyr::pglmm(scale(log(BodyCond)) ~ 1 +  
                            (1|Age) + (1|Sex) + (1|Molt) + (1|Year) + (1|Latitude_decimal_degrees)+
                            (1|Diet), data = df_final,
                    cov_ranef = list(sp = phylogeny), REML = TRUE, verbose = TRUE)


## Model comparison
pvalue <- pchisq(2*(model_full_2$logLik - model_null$logLik), df = 4, lower.tail = FALSE)
round(pvalue,10)

## Full was significantly different
summary(model_full_2)

## Predicted values
df_final$Predicted <- NA
df_final$Predicted <- as.data.frame(predict(model_full_2))[,1]

## Model performance
res <- residuals(model_full_2)
binnedplot(df_final$Predicted, res)

# Getting confidence intervals through bootstrapping for the categoreis
n_bootstrap <- 5
n_coefs <- length(fixef(model_full_2)$Value)
bootstrap_coefs <- matrix(NA, nrow = n_bootstrap, ncol = n_coefs)

## Confidence intervals
for (i in 1:n_bootstrap) {
  indices <- sample(nrow(df_final), replace = TRUE)
  sampled_data <- df_final[indices, ]
  
  tryCatch({
    bootstrap_model <- phyr::pglmm(scale(log(BodyCond)) ~ scale(sqrt(ForestCover)) + scale(normalize(PCHB1)) + scale(Temp) + 
                                     scale(log(Rainfall)) + scale(normalize(Speciesdiversity)) + 
                                     (1|Age) + (1|Sex) + (1|Molt) + (1|Year) + (1|Latitude_decimal_degrees)+
                                     (1|Diet), 
                                   data = sampled_data, cov_ranef = list(sp = phylogeny))
  }, error = function(e) {
    cat("Error occurred in iteration", i, ":", conditionMessage(e), "\n")
    next
  })
  bootstrap_coefs[i, ] <- fixef(bootstrap_model)$Value
  print(i)
}

##write.csv(bootstrap_coefs, "./A.csv")
conf_intervals <- apply(bootstrap_coefs, 2, quantile, c(0.025, 0.975))
print(conf_intervals)
conf_intervals <- conf_intervals[,-1]
bootstrap_coefs <- bootstrap_coefs[,-1]
conf_intervals <- conf_intervals[,-c(1,3)]
bootstrap_coefs <- bootstrap_coefs[,-c(1,3)]

## Plotting results
colnames(conf_intervals) <- c("Forest_cover", "Habitat_breadth", "Temperature", "Rainfall", "Species_diversity")

## Plotting results
# Create a data frame for forest plot
forest_data <- data.frame(
  Variables = colnames(conf_intervals),
  Estimate = colMeans(bootstrap_coefs),
  CI.Lower = conf_intervals[1, ],
  CI.Upper = conf_intervals[2, ]
)

forest_data <- data.frame(
  Variables = c("Habitat_breadth", "Rainfall", "Species_diversity"),
  Estimate = c(-0.05718933, 0.19410780, -0.25800196),
  CI.Lower = c(-0.06812262, 0.17388757, -0.28496741),
  CI.Upper = c(-0.05045376, 0.21546152, -0.22851038)
)

# Create the forest plot
forest_plot <- ggplot(forest_data, aes(x = Estimate, y = Variables)) +
  geom_point(size = 3, color = "black") +
  geom_errorbarh(aes(xmin = CI.Lower, xmax = CI.Upper), height = 0, size = 1, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Confidance interval",
       y = "Variables") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 18), 
    panel.grid.major = element_line(color = "gray", size = 0.2),  
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white")
  )

# Print the forest plot
jpeg("./Significant_variables.jpg", width = 10, height = 8, units = 'in', res = 800)
print(forest_plot)
dev.off()

## End of code
## save.image("./Documents/Habitat.RData")
## write.csv(df_final, "./Documents/Habitat_final.csv")
##############################################################################################################################################
##############################################################################################################################################
