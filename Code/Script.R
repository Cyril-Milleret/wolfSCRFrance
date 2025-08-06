rm(list=ls())
#LOAD PACKAGES
library(nimble)
library(nimbleSCR)
library(raster)
library(sf)
library(stars)
library(terra)
library(R.utils)
library(basicMCMCplots)
library(coda)
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
### ==== 1. GENERAL VARIABLES DECLARATION ====
myVars <- list(
  ## WORKING DIRECTORY & MODEL NAME
  # WD = "C:/My_documents/NIMBLE/WOLVERINE",
  WD = "C:/Personal_Cloud/OneDrive/Work/CNRS/Papers/SCRWolfFrance/SCRWolfFrance",
  # modelName = "s31_[W,Hum,Gras,F]sigma[sex]p0[reg[dpt1],geaco,snow,rds.indcov[prev,sex]Loc2024FinalToPublish",
  
  # HABITAT SPECIFICATIONS
  HABITAT = list( habResolution = 10000,
                  habBuffer = 20000),
  
  # NGS DATA SPECIFICATIONS
  DATA = list( years = 2024,
               sex = c("F","M"),                   ## "Hann","Hunn","Ukjent"
               samplingMonths = list(11,12,1:3)), ## list(10:12,1:4), list(1:XXX), list(XX:XX,YY:YY)
  
  # DETECTORS SPECIFICATIONS
  DETECTORS = list( detSubResolution = 1000,
                    detResolution = 5000,
                    detDeadResolution = 15000),
  
  # DATA GENERATION
  DETECTIONS = list( maxDetDist = 40000,
                     resizeFactor = 0,
                     aug.factor = 4),

  ## MISCELLANEOUS
  plot.check = TRUE)

years <- myVars$DATA$years
nYears <- length(years)
YEARS <- lapply(years, function(x)c(x,x+1))


##source functions
sourceDirectory(file.path(myVars$WD,"Functions"), modifiedOnly = FALSE)

## ==== 1. LOAD DATA ####
## ====  1.1 DNA DATA ####
## ====   1.1.1 ALIVE DATA ####
load(file="C:/Personal_Cloud/OneDrive/Work/CNRS/Papers/SCRWolfFrance/SCRWolfFrance/Data/DNA.RData")

## ====  1.2 GIS DATA ####
## ====   1.2.1 10*10 km SAMPLING GRID #####
grid1010 <- read_sf(file.path(myVars$WD,"Data","Indices_biologiques_recolte_Hiver2021_2022_grid_10x10km_03102023.shp"))
st_crs(DNA) <- st_crs(grid1010)
## ====   1.2.2 EUROPE ####
Europe <- read_sf(file.path(myVars$WD,"Data","EuropeCountries"))
##subset to neighboring countries 
Europe <- Europe[Europe$NAME %in%c("France","Spain","Italy","Switzerland","Germany","Luxembourg","Belgium","Andorra","Monaco"),]
Europe <- st_transform(Europe,crs = st_crs(DNA))
## ====   1.2.3 FRANCE ####
France <- Europe[Europe$NAME %in%c("France"),]
France <- st_transform(France,crs = st_crs(DNA))
plot(Europe$geometry)
plot(DNA$geometry,add=T,col="red")

## ====   1.2.4 DEPARTEMENT #####
Departement <- read_sf(file.path(myVars$WD,"Data","departements-20180101-shp","departements-20180101.shp"))
Departement <- st_transform(Departement,crs = st_crs(DNA))

plot(Departement$geometry,add=T)

## ====   1.2.5 REGION #####
Region <- read_sf(file.path(myVars$WD,"Data","Region","regions_2015_metropole_region.shp"))
Region <- st_transform(Region,crs = st_crs(DNA))

plot(Region$geometry,add=T,col="red")



## ====   1.2.8 LCIE WOLF DISTRIBUTION #####
#2018 DISTRIBUTION MAP 
# LCIEDistribution <- read_sf(file.path(myVars$WD,"GIS","LCIE_wolfDistribution","2018_06_06_Wolf_IUCN_RedList.shp"))
# str(LCIEDistribution)
# LCIEDistribution$SPOIS <- ifelse(LCIEDistribution$SPOIS %in% "Sporadic",1,3)
#2022 DISTRIBUTION MAP 
# LCIEDistribution2012permanent <- read_sf(file.path(myVars$WD,"GIS","LCIE_wolfDistribution","Wolf_LCIE_2012","Clip_2012_12_01_Wolves_permanent.shp"))
# LCIEDistribution2012permanent$SPOIS <- 3
# LCIEDistribution2012sporadic <- read_sf(file.path(myVars$WD,"GIS","LCIE_wolfDistribution","Wolf_LCIE_2012","Clip_2012_12_01_Wolves_sporadic.shp"))
# LCIEDistribution2012sporadic$SPOIS <- 1
# LCIEDistribution2012 <- rbind(LCIEDistribution2012permanent, LCIEDistribution2012sporadic)
# plot(LCIEDistribution2012["SPOIS"])
# 
# LCIEDistribution2012 <- st_transform(LCIEDistribution2012, crs=st_crs(DNA))
# LCIEDistribution <- st_transform(LCIEDistribution, crs=st_crs(DNA))
# 
# # Both grids classified wolf presence as sporadic or permanent (we scored 1 for sporadic and 3 for permanent). In a second step,
# #we summed the two grids (weighted by the area overlap between the 10km LCIE grid and the 5km habitat grid used in the model) 
# #to create an index of the historical presence of wolves in each 
# #habitat grid cell (from 0, no wolf presence registered in the two periods, to 6, permanent presence during the two periods).
# plot(France$geometry)
# plot(LCIEDistribution2012["SPOIS"],add=T)
# 
# str(LCIEDistribution)
# plot(st_transform(LCIEDistribution,st_crs(Europe))$geometry,add=T)
# LCIEDistribution <- st_transform(LCIEDistribution,crs = st_crs(DNA))

# DATA AND SCRIPT FROM MARRUCCO ET AL  
iucn_2012_1 <- read_sf(file.path(myVars$WD,"Data","LCIE_wolfDistribution","Wolf_LCIE_2012","Clip_2012_12_01_Wolves_permanent.shp"))
iucn_2012_1$SPOIS <- 3
iucn_2012_1 <- st_transform(iucn_2012_1, st_crs(DNA))
plot(st_geometry(DNA))
plot(st_geometry(iucn_2012_1), add = T, col = "lightskyblue4")

iucn_2012_2 <- read_sf(file.path(myVars$WD,"Data","LCIE_wolfDistribution","Wolf_LCIE_2012","Clip_2012_12_01_Wolves_sporadic.shp"))
iucn_2012_2$SPOIS <- 1
iucn_2012_2 <- st_transform(iucn_2012_2, st_crs(DNA))
plot(st_geometry(iucn_2012_2), add = T, col = "lightskyblue2")

iucn_2018 <- read_sf(file.path(myVars$WD,"Data","LCIE_wolfDistribution","2018_06_06_Wolf_IUCN_RedList.shp"))
iucn_2018 <- st_transform(iucn_2018, st_crs(DNA))
plot(st_geometry(DNA))
plot(iucn_2018[ ,"SPOIS"], add = T)
## Code back to numeric
iucn_2018$SPOIS <- ifelse(iucn_2018$SPOIS == "Sporadic", 1, 3)

## ====   1.2.9 EFFORT #####
## ====     1.2.9.1 POTENTIAL Baduin et al 2023 #####
#LOAD GRID
load(file.path(myVars$WD, "Data","Effort","gridFr.RData"))
EffortGrid <- st_as_sf(gridFr)
EffortGrid <- st_transform(EffortGrid,crs = st_crs(DNA))
#LOAD EFFORT VALUES FROM SARAH
load(file.path(myVars$WD, "Data","Effort","effort.RData"))
#FILL IN EFFORT VALUE
EffortGrid$value <- effort$eff_2020_2021
#plot
plot(EffortGrid["value"]$geometry)
plot(EffortGrid["value"])
plot(DNA$geometry,add=T,col="black")
# ## ====     1.2.9.2 GEACO #####
# effort <- read.csv(file.path(myVars$WD, "GIS","EffortOFB","Bilan2021mars2024.csv"), 
#                    fileEncoding="latin1")
# ##CONVERT TIME AND SEASON ##
# effort$date1 <- as.POSIXct(strptime(effort$Date, "%d/%m/%Y"))
# effort <- effort[!is.na(effort$date1 ),]
# effort$Year <- as.numeric(format(effort$date1,"%Y"))
# effort$mois <- as.numeric(format(effort$date1,"%m"))
# effort$saisonyear <- ifelse(effort$mois < 11, effort$Year, effort$Year+1) #--- need to change for other species
# #MAKE IT NUMERIC
# effort$Heures <- as.numeric(effort$Heures)
# #SUBSET TO SEASON OF INTEREST
# #MONTHS
# effort <- effort[effort$mois %in% unlist(myVars$DATA$samplingMonths),]
# #YEARS
# effort <- effort[effort$saisonyear %in% myVars$DATA$years ,]
# #CHECK 
# effort$date1
# # subset unknown commune
# effort <- effort[!effort$Code.Insee %in% "-",]
# effort$insee <- effort$Code.Insee
# 
# #summarize per communes 
# df <- 
#   effort %>%
#   group_by(saisonyear,insee) %>%
#   summarize(#across(saisonyear,list(
#     n = n(),
#     sum = sum(Heures),
#     logSum = log(sum(Heures))
#     
#   )
# 
# 
# #mergeCommunes
# Communes <- merge(Communes,df,by="insee")

# save(Communes,file="C:/Personal_Cloud/OneDrive/Work/CNRS/Papers/SCRWolfFrance/SCRWolfFrance/Data/Geaco.RData")
load("C:/Personal_Cloud/OneDrive/Work/CNRS/Papers/SCRWolfFrance/SCRWolfFrance/Data/Geaco.RData")
mapview::mapview(Communes["n"],border=NA)

## ====   1.2.10 SNOW #####
# dropbox <- "C:/Users/cymi/Dropbox (Old)/AQEG Dropbox/AQEG Team Folder/RovQuant/DATA/"
# SNOW <- stack(paste(dropbox,"/GISData/SNOW/ModisSnowCover0.1degrees/AverageSnowCoverModisSeason2008_2024_Wolf.tif", sep=""))
# SNOW <- stack(file.path(myVars$WD,"GIS","SNOW","AverageSnowCoverModisSeason2008_2024_Wolf.tif"))
# 
# ## RENAME THE LAYERS
# names(SNOW) <- paste(2008:2023,(2008:2023)+1, sep="_")
# ## SELECT SNOW DATA CORRESPONDING TO THE MONITORING PERIOD
# SNOW <- SNOW[[paste("X", myVars$DATA$years-1, "_", myVars$DATA$years, sep="")]]
# SNOW <- raster::crop(SNOW, st_transform(Europe,crs = st_crs(SNOW)))

# save(SNOW,file="C:/Personal_Cloud/OneDrive/Work/CNRS/Papers/SCRWolfFrance/SCRWolfFrance/Data/SNOW.RData")
load("C:/Personal_Cloud/OneDrive/Work/CNRS/Papers/SCRWolfFrance/SCRWolfFrance/Data/SNOW.RData")

## ====   1.2.11 HUMAN DENSITY #####
# HumanDensity1 <- raster(file.path(myVars$WD,"GIS","HumanDensity","GHS_POP_E2020_GLOBE_R2023A_54009_1000_V1_0_R4_C19.tif"))
# HumanDensity2 <- raster(file.path(myVars$WD,"GIS","HumanDensity","GHS_POP_E2020_GLOBE_R2023A_54009_1000_V1_0_R4_C18.tif"))
# HumanDensity <- merge(HumanDensity2, HumanDensity1)

# save(HumanDensity,file="C:/Personal_Cloud/OneDrive/Work/CNRS/Papers/SCRWolfFrance/SCRWolfFrance/Data/HumanDensity.RData")
load("C:/Personal_Cloud/OneDrive/Work/CNRS/Papers/SCRWolfFrance/SCRWolfFrance/Data/HumanDensity.RData")

plot(HumanDensity)
#use focal to get an approx 10km resolution 
plot(log(HumanDensity))
plot(France$geometry,add=T)

## ====   1.2.12 ROAD DENSITY #####
#LENGTH OF ROADS 
# load(file.path(myVars$WD,"GIS","Roads","Roads1.RData"))#ROADS1 IS THE GOOD ONE
# roads <- aggregate(roads1,fact=10)
# plot(roads)

# save(roads,file="C:/Personal_Cloud/OneDrive/Work/CNRS/Papers/SCRWolfFrance/SCRWolfFrance/Data/roads.RData")
load("C:/Personal_Cloud/OneDrive/Work/CNRS/Papers/SCRWolfFrance/SCRWolfFrance/Data/roads.RData")

## ====   1.2.13 FOREST #####
#Corinne land cover data should be downloaded https://land.copernicus.eu/en/products/corine-land-cover
#CLC CLIPPED TO FRANCE IN QGIS
CLC_FR <- st_read(file.path("C:/Personal_Cloud/OneDrive/Work/CNRS/SCR","GIS","CLC","u2018_clc2018_v2020_20u1_fgdb","CLC_FR.shp"))

CLCForest <- CLC_FR[CLC_FR$Code_18 %in% c("311", "312", "313"),]
## ====   1.2.14 GRASSLAND #####
CLCGrassland <- CLC_FR[CLC_FR$Code_18 %in% c("321", "322", "324"),]


## ==== 2. Check NGS DATA ####
barplot(table(DNA$Year))
table(DNA$mois,DNA$Year)
NDet <- tapply(DNA$Id, DNA$Id, length)
mean(NDet)


## ====     2.3.2 PLOT AND SUMMARY TO CHECK EVERYTHING IS BEING USED ####
nrow(DNA)

plot(France$geometry)
plot(DNA$geometry,pch=16, cex=0.5, col="red",add=T)

#MAKE A TABLE WITH SUMMARY OF DETECTIONS PER MONTHS AND SEXE
matDets <- matrix(NA,nrow=4,ncol=6)
colnames(matDets) <- c("Nov","Dec","Jan","Feb","Mar","Total")
rownames(matDets) <- c("F","M","NI","Total")
tab <- table(DNA$SEXE,DNA$mois)
matDets[1:3,1:5]<- tab[c(2,3,1),c(4,5,1,2,3)]
matDets[,"Total"] <- rowSums(matDets,na.rm=T)
matDets["Total",] <- colSums(matDets,na.rm=T)

apply(table(DNA$Id,DNA$SEXE),2,function(x) sum(x>0))


###
# DNA <- DNA[,c(11,17:21)]
# save(DNA, file="C:/Personal_Cloud/OneDrive/Work/CNRS/Papers/SCRWolfFrance/SCRWolfFrance/Data/DNA.RData")


## ==== 3. GENERATE HABITAT ====
## ====   3.1 DEFINE AREA SEARCHED ====
#BUFFER THE DETECTIONS TO FIND THE HABITAT
areaSearched <- st_buffer(DNA,dist = 100000)
areaSearched <- st_union(areaSearched)
#REMOVE ALL AREAS OUTSIDE OF FRANCE (NOT SEARCHED)
areaSearched <- st_intersection(areaSearched, France)
areaSearched <- st_union(areaSearched)
areaSearched <- st_as_sf(areaSearched)

plot(areaSearched)


## ====   3.2 DEFINE HABTIAT EXTENT ====
#BUFFER AREA SEARCHED
habitatExtent <- st_buffer(areaSearched,dist = myVars$HABITAT$habBuffer)
habitatExtent <- st_union(habitatExtent)
#REMOVE THE SEA FROM THE BUFFER
habitatExtent <- st_intersection(habitatExtent, Europe)
habitatExtent <- st_union(habitatExtent)
habitatExtent <- st_as_sf(habitatExtent)
plot(habitatExtent)
plot(France$geometry,add=T)

## ====   3.3 RASTERIZE HABITAT ====
#BASE THE HABITAT ON THE 10*10 KM GRID USED FOR THE MONITORING 
mat <- as.data.frame(cbind(st_coordinates(grid1010)[,c("X","Y")], rep(1,nrow(grid1010))))
colnames(mat) <- c("x","y","z")
mat[,c(1:2)] <- mat[,c(1:2)] - 5000## for some reasons the cells dont align; do it manually. 
template.rSpa  <- rast(mat, type="xyz",crs=crs(DNA))# rasterFromXYZ(mat)
template.rSpa <- crop(mask(template.rSpa,habitatExtent),habitatExtent)
st_crs(template.rSpa) == st_crs(grid1010)

#SOME SEARCHED GRID CELLS ARE OUTSIDE THE HABITAT (VERY LOW OVERLAP WITH HABITAT; DO NOT INCLUDE THEM)
#keep habitat prop>50%
template.rSpa <- rasterize(habitatExtent, template.rSpa, getCover=TRUE)        ## Add field=1 so its 1 for study area and NA for the rest.
template.rSpa[template.rSpa < 50] <- 0
template.rSpa[template.rSpa >= 50] <- 1

template.r <- r <- raster(template.rSpa)

#PLOT CHECK 
plot(template.r)
plot(grid1010$geometry,add=T)
plot(France$geometry,add=T)
plot(habitatExtent,add=T)
plot(st_centroid(grid1010)$geometry,add=T,pch=16,cex=0.5)

# STEP REQUIRED WHEN AGGREGATING THE HABITAT? 
habitat.r <- template.r# habitat.r <- raster::rasterize(habitatExtent, template.r,crs=st_crs(DNA))
##AGGREGATE HABITAT RESOLUTION TO 20KM ?
habitat.r1 <- aggregate(habitat.r,fact=1)

# GET THE HABITAT OBJECTS NECESSARY FOR NIMBLESCR
habitatxy <- data.frame(raster::coordinates(habitat.r)[!is.na(habitat.r)[],])#,df=T)#[is.na(habitat.r[][,1]),]
habitatxysf <- st_as_sf(habitatxy,coords = c("x", "y"), crs = st_crs(DNA))#double check the projection system
habitatSF.pol <- sf::st_as_sf(stars::st_as_stars(habitat.r), 
                              as_points = FALSE, merge = F)

#PLOT CHECK 
plot(habitat.r)
plot(habitatxysf$geometry,add=T)
plot(France$geometry,add=T)
plot(DNA$geometry,add=T,pch=16,col="red")
plot(grid1010$geometry,add=T,col=NA)
#subset DNA samples falling outside habitat
# idremoved <- raster::extract(habitat.r,DNA)
# which(is.na(idremoved))

## ====   3.4 DEFINE DETECTOR GRID ====
detector.r <- rasterize(areaSearched, template.rSpa)
#BASE THE DETECTOR GRID ON SUBDETECTORS
subdetector.r <- disagg(detector.r, fact=10)#1km resolution for the PAB
myDetectors <- MakeSearchGrid( subdetector.r=subdetector.r,
                                 detResolution=myVars$DETECTORS$detResolution,
                                 plot = F
)
myDetectors$detector.xy <- st_coordinates(myDetectors$main.detector.sf)
colnames(myDetectors$detector.xy) <- c("x","y")

#PLOT CHECK 
plot(habitat.r)
plot(France$geometry,add=T)
plot(myDetectors$main.detector.sf$geometry,add=T,col="red",pch=16,cex=0.1)


#CHECK THAT SUBDECTORS ALIGNS WITH MAIN DETECTORS 
plot(habitat.r)
myDetectors$detector.xy
tmp <- myDetectors$main.detector.sf[myDetectors$main.detector.sf$main.cell.id%in%4835,]
plot(tmp$geometry,add=T,pch=16)
tmpsub <- myDetectors$detector.sf[myDetectors$detector.sf$main.cell.id%in%4835,]
plot(tmpsub$geometry,add=T,col="red")
## ==== 4. GENERATE HABITAT-LEVEL COVARIATES ====
## ====   4.2 WOLF PRESENCE ====
## EXTRACT WOLF PERMANENT PRESENCE  IN EACH HABITAT CELL
habitatGrid  <- sf::st_as_sf(stars::st_as_stars(template.rSpa), 
                    as_points = FALSE, merge = F)

sizehabitatCell <- st_area(habitatGrid)[1]
habitatGrid$id <- 1:nrow(habitatGrid)
intersection <- st_intersection(habitatGrid, iucn_2012_1) %>%
  mutate(iucn = st_area(.)) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(IUCN = sum(iucn*SPOIS)/sizehabitatCell)
  #summarise(IUCN = sum(SPOIS))

habitatGrid <- habitatGrid %>%
  left_join(intersection, by = "id")
habitatGrid$IUCN[is.na(habitatGrid$IUCN)] <- 0
plot(habitatGrid[,"IUCN"])

## Extract LCIE wolf sporadic presence in each habitat grid cell
intersection <- st_intersection(habitatGrid, iucn_2012_2) %>%
  mutate(iucn = st_area(.)) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(iucn_2 = sum(iucn*SPOIS)/sizehabitatCell)
  #summarise(iucn_2 = sum(SPOIS))

tmp <- habitatGrid %>%
  left_join(intersection, by = "id")
tmp$iucn_2[is.na(tmp$iucn_2)] <- 0
habitatGrid$IUCN <- habitatGrid$IUCN + tmp$iucn_2
plot(habitatGrid[,"IUCN"])

## Extract LCIE wolf presence in each habitat grid cell
intersection <- st_intersection(habitatGrid, iucn_2018) %>%
  mutate(iucn = st_area(.)) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  #summarise(iucn_2 = sum(SPOIS))
  summarise(iucn_2 = sum(iucn*SPOIS)/sizehabitatCell)

tmp <- habitatGrid %>%
  left_join(intersection, by = "id")
tmp$iucn_2[is.na(tmp$iucn_2)] <- 0

habitatGrid$IUCN <- habitatGrid$IUCN + tmp$iucn_2
plot(habitatGrid[,"IUCN"])

habitatGrid$IUCN <- scale(habitatGrid$IUCN)

## ====   4.3 HUMAN DENSITY ====
habitatxysf1 <- st_transform(habitatxysf, st_crs(HumanDensity))
HumanCov1 <- raster::extract( HumanDensity, 
                             habitatxysf1,
                             buffer=5800)#radius of 10km2 circle

HumanCov <- unlist(lapply(HumanCov1,function(x) sum(x,na.rm=T)))
plot(HumanDensity)
#check if nas?
sum(is.na(HumanCov))


## ====   4.4 PROPORTION OF FOREST ====
habitatSF.pol1 <- st_transform(habitatSF.pol, st_crs(CLCForest))
habitatSF.pol1$ID <- 1:nrow(habitatSF.pol1) 

cellForest <- st_intersection(habitatSF.pol1 , CLCForest) %>%
  mutate(areaInter = st_area(.)) %>%
  group_by(ID) %>%
  summarise(areaforestCell = sum(areaInter)) %>%
  mutate(propForestCell = areaforestCell / sizehabitatCell) %>%
  as_tibble() %>%
  select(ID, propForestCell)

#PLOT CHECK  
forest.r <- habitat.r
forest.r[] <- 0

forest.r[which(!is.na(habitat.r[]))[cellForest$ID]] <- cellForest$propForestCell
plot(forest.r)

## ====   4.5 PROPORTION OF GRASSLAND ====
cellGrassland <- st_intersection(habitatSF.pol1 , CLCGrassland) %>%
  mutate(areaInter = st_area(.)) %>%
  group_by(ID) %>%
  summarise(areaGrasslandCell = sum(areaInter)) %>%
  mutate(propGrasslandCell = areaGrasslandCell / sizehabitatCell) %>%
  as_tibble() %>%
  select(ID, propGrasslandCell)
#plot#check 
grassland.r <- habitat.r
grassland.r[] <- 0

grassland.r[which(!is.na(habitat.r[]))[cellGrassland$ID]] <- cellGrassland$propGrasslandCell
plot(grassland.r)


## ====   4.4 BIND COVARIATES ====
habCovs <-  cbind(habitatGrid$IUCN,
                  log(HumanCov+0.0001),
                  grassland.r[!is.na(habitat.r[])],
                  forest.r[!is.na(habitat.r[])])
#SCALE COVARIATES 
habCovs <- scale(habCovs)
colnames(habCovs) <- c(#"DEM",
                       "Wolf","Human","Grassland","Forest")

##PlotCheck
tmp <- habitat.r
for(i in 1:dim(habCovs)[2]){
  tmp[!is.na(habitat.r[])] <- habCovs[,i]
  plot(tmp,main=colnames(habCovs)[i])
  plot(France$geometry,add=T)
}

## ==== 5. GENERATE DECTOR-LEVEL COVARIATES ====
## ====   5.1 EFFORT ====
## ====         5.1.1 SARAH
id <- as.numeric(st_intersects(myDetectors$main.detector.sf, EffortGrid["value"]))
myDetectors$main.detector.sf$effortCov <- EffortGrid$value[id]

tmpr <- myDetectors$maindetector.r
# id <- terra::extract(tmpr, myDetectors$main.detector.sf,cells=T)


## if NA returns the average value of the cells within 20000m 
isna <- which(is.na(myDetectors$main.detector.sf$effortCov))#, 1, function(x)any(is.na(x))))
if(length(isna)>0){
whichClose <- apply(st_distance(EffortGrid, myDetectors$main.detector.sf[isna, ]),2,function(x){order(x) [1:3]})
for(i in 1:length(isna)){
  myDetectors$main.detector.sf$effortCov[isna[i]] <- mean(EffortGrid$value[whichClose[,i]],na.rm=T)
}
}
##check 
tmpr[terra::extract(tmpr, myDetectors$main.detector.sf, cells=T)[,3]]<- myDetectors$main.detector.sf$effortCov
 plot(tmpr)

 ## ====         5.1.2 GEACO
nbVisits <- 
   Communes %>%  
   # filter over the second condition
   st_intersects(myDetectors$grid.poly, .)
 
 
myDetectors$grid.poly$nbVisits <- 0
 for(i in 1:nrow(myDetectors$grid.poly)){
   myDetectors$grid.poly$nbVisits[i] <-  sum(Communes$n[nbVisits[[i]]])
 }
 
 # grid$nbHours <- 0
 # for(i in 1:nrow(grid)){
 #   grid$nbHours[i] <-  sum(Communes2021$sum[nbVisits[[i]]])
 # }
myDetectors$grid.poly$nbVisits[myDetectors$grid.poly$nbVisits>15] <- 15
 #quick ckeck 
mapview::mapview(list(myDetectors$grid.poly["nbVisits"],Communes["n"]))

## ====   5.2 REGION ====
Departement$id <- 1:nrow(Departement)
tmpDis <- st_distance(myDetectors$main.detector.sf, st_simplify(Departement,dTolerance = 1000,preserveTopology = T))
# tmpDis[500,]
DetRegion <- 0
for(i in 1:nrow(tmpDis)){
  DetRegion[i] <- which.min(tmpDis[i,])
}
DetRegion <- apply(tmpDis, 1, which.min)

table(DetRegion)
centr <- st_centroid(Departement)$geometry
text(st_coordinates(centr))

Departement$code_insee[c(53,11,17,16,42,100,99,79,78,87,60)]
Departement$code_insee[c(88,89,95,97,4,52,18,37)]
Departement$code_insee[c(27,34,2,19,57,10,15,24,96,3,8,14,7,101,102,71,64,65,59)]

inseRegion1 <- c("01","39","71","25","21","70","90","68","69D","69M","74","58","03","89","42")

inseRegion2 <- c("07","63","43","26")

#inseRegion2 <- c("26","05","84","04","13","83","06")
inseRegion3 <- c("66","09","11","34","31","76","81","12","30","32","47",
                 "48","43","15","63","19","46","82","71","84","13")


DetRegion[DetRegion %in% which(Departement$code_insee %in%c(inseRegion1)) ] <- "A"
DetRegion[DetRegion %in% which(Departement$code_insee %in%c(inseRegion2)) ] <- "Z"
DetRegion[DetRegion %in% which(Departement$code_insee %in%c(inseRegion3)) ] <- "C"

# DetRegion[DetRegion %in% which(Departement$code_insee %in%c("84","13")) ] <- "84"
# DetRegion[DetRegion %in% which(Departement$code_insee %in%inseRegion2) ] <- "B"
# DetRegion[DetRegion %in% which(Departement$code_insee %in%inseRegion3) ] <- "C"

trapIntercept <- as.numeric(as.factor(DetRegion))
plot(France$geometry)
col <- rainbow(max(as.numeric(as.factor(DetRegion)))+5)
plot(myDetectors$main.detector.sf$geometry,col=col[as.numeric(as.factor(DetRegion))+2])
plot(Departement$geometry,add=T)
text(st_coordinates(centr))
mapview::mapview(Departement)

## ====   5.3 SNOW ====
snowCov <- raster::extract(SNOW, st_transform(myDetectors$main.detector.sf,crs = st_crs(SNOW)))

## if NA returns the average value of the cells within 20000m 
snowXY <- coordinates(SNOW) 
snowXY <- st_as_sf(data.frame(snowXY), coords=c("x","y"), crs=st_crs(SNOW))

isna <- which(is.na(snowCov))#, 1, function(x)any(is.na(x))))
if(length(isna)>0){
  whichClose <- apply(st_distance(snowXY, st_transform(myDetectors$main.detector.sf[isna, ],crs=st_crs(SNOW))),2,function(x){order(x) [1:8]})
  for(i in 1:length(isna)){
    snowCov[isna[i]] <- mean(SNOW[whichClose[,i]],na.rm=T)
  }
}

##check 
plot(tmpr)
tmpr[!is.na(tmpr[])]<- snowCov
plot(tmpr)



## ====   5.4 ROADS ====
roads <- focal(roads, matrix(1,3,3),mean,na.rm=T)
# plot(Road.r)
RoadCov <- terra::extract(roads, st_transform(myDetectors$main.detector.sf,crs = st_crs(roads)))
isna <- which(is.na(RoadCov))#, 1, function(x)any(is.na(x))))


roadXY <- coordinates(roads) 
roadXY <- st_as_sf(data.frame(roadXY), coords=c("x","y"), crs=st_crs(roads))

isna <- which(is.na(RoadCov))#, 1, function(x)any(is.na(x))))
if(length(isna)>0){
  whichClose <- apply(st_distance(roadXY, st_transform(myDetectors$main.detector.sf[isna, ],crs=st_crs(roads))),2,function(x){order(x) [1:8]})
  for(i in 1:length(isna)){
    RoadCov[isna[i]] <- mean(roads[whichClose[,i]],na.rm=T)
  }
}
tmpr[!is.na(tmpr[])]<- RoadCov
plot(tmpr)
## ====   5.2 BIND COVARIATES ====
trapCov <- scale(cbind(myDetectors$grid.poly$nbVisits , snowCov, log(RoadCov)))
colnames(trapCov) <- c("Effort","Snow","Roads")

cor(trapCov)

##PlotCheck
  tmp <- myDetectors$maindetector.r
  for(i in 1:dim(trapCov)[2]){
    tmp[!is.na(myDetectors$maindetector.r[])] <- trapCov[,i]
    plot(tmp,main=colnames(trapCov)[i])
    plot(France$geometry,add=T)
    
  }
  
  tmp[!is.na(myDetectors$maindetector.r[])] <- trapIntercept
  plot(tmp,main="Region")
  plot(France$geometry,add=T)



## ==== 6. ASSIGN DETECTIONS TO DETECTORS ====
DNA$Year <- DNA$saisonyear
myData.alive <- AssignDetectors_v3sf( myData = DNA
                                      ,                
                                      myDetectors = myDetectors$main.detector.sf
                                      ,
                                      mysubDetectors = myDetectors$detector.sf
                                      ,
                                      radius = myVars$DETECTORS$detResolution)

## ==== 7. MAKE Y ==== 
y.ar <- MakeYsf( myData = myData.alive$myData.sp,
                 myDetectors = myDetectors$main.detector.sf,
                 method = "Binomial",
                 returnIdvector = TRUE)

## ==== 8. INDIVIDUAL COVARIATES ====
## ====   8.1 SEX ====
sex <- 0
i=183
for(i in 1:length(y.ar$Id.vector)){
  tmp <- DNA[DNA$Id %in% y.ar$Id.vector[i],]
  tmpsexe <- unique(tmp$SEXE)
  #Check if id assign to all 4 categories, if yes, print id
  if(length(tmpsexe)>2){print(i)}
  #if NI get the sex from other samples
  if(sum(tmpsexe%in% "NI")>0){
  tmpsexe <- tmpsexe[!names(tmpsexe)%in% "NI"]
  }
  
  if(length(tmpsexe)>1){
   tabSexe <- table(tmp$SEXE)
   #find most comon sexe
   whichMostComon <- which(tabSexe>1)#[tabSexe>0]
    if(length(whichMostComon)>0){
      sex[i] <- tabSexe[whichMostComon]
      #print(i)
    }else{#assign sex to highest quality sample
      sex[i] <-  tmp$SEXE[which.max(tmp$INDEX_CAL)]
    }
    
  }
   sex[i] <- tmpsexe[1]
}
sex[sex%in% "XX"] <- "F"
sex[sex%in% "XY"] <- "M"
#LET MODEL ASSIGN THE SEX OF THE INDIVIDUAL
sex[sex%in% "NI"] <- NA
sex[sex%in% "#N"] <- NA

sex[sex%in% "F"] <- "0"
sex[sex%in% "M"] <- "1"
sex <- as.numeric(sex)
#CHECK IT 
table(sex,useNA="always")

## ====   8.2 PREVIOUS DETECTION ====
#check if individuals was detected during the previous years 3 years
load(file.path(myVars$WD,"Data/DNAAllPev.RData"))

#SOME CHECKS 
tab1 <- table(DNAAllPev$Id,DNAAllPev$saisonyear)
whichDets <- apply(tab1,1,function(x) sum(x)>1)#ids detected at least during two winters previous to 2021-2023
idGen <- names(whichDets)[whichDets]#unique(DNA1$GENOTYPE_TRANS)
PrevDets <- as.numeric(y.ar$Id.vector %in% idGen)
sum(PrevDets)/length(PrevDets)#PROPORTION OF IDS DETECTED THIS YEAR VS LAST YEAR

#mean number of detections 
mean(rowSums(y.ar$y.ar[y.ar$Id.vector %in% idGen,]))


## ==== 9. CHECK PATTERNS OF INDIVIDUAL DETECTIONS ====
distances <- CheckDistanceDetectionsV2sf( y = y.ar$y.ar, 
                                         detector.xy = myDetectors$main.detector.sf, 
                                         max.distance = 40000,
                                         method = "pairwise",
                                         plot.check = F)

#EXLCUDE DETECTIONS OF INDIVIDUALS DOING "LARGE MOVEMENT" DURING THE SAMPLING. 
  par(mfrow=c(1,1),mar=c(1,1,1,1))
  if(sum(distances$y.flagged) > 0){
    affected.ids <- which(apply(distances$y.flagged,1,sum)>0)
    for(i in affected.ids){
      plot(st_geometry(myDetectors$main.detector.sf),pch=16, cex=0.1)
      plot(st_geometry(France), add = T,border="red")
      
      tmp <- DNA[DNA$Id == y.ar$Id.vector[i], ]
      tmp <- tmp[order(tmp$date), ]
      tmp.xy <- st_coordinates(tmp)
      n.det <- nrow(tmp.xy)
      
      plot(st_geometry(tmp),add=T, col = "orange", pch = 16, cex = 1)
      
      arrows(x0 = tmp.xy[1:(n.det-1),1], y0 = tmp.xy[1:(n.det-1),2],
             x1 = tmp.xy[2:n.det,1], y1 = tmp.xy[2:n.det,2], length = 0.1, lwd = 1,col="orange")
      
      plot(st_geometry(myDetectors$main.detector.sf[which(y.ar$y.ar[i,] > 0), ]), pch = 16, col = "red",add=T)
      
      tmp2 <- myDetectors$main.detector.sf[which(y.ar$y.ar[i,] > 0 & distances$y.flagged[i,] == 1), ]
      plot(st_geometry(tmp2), col = "blue", pch = 13, cex = 1.5, lwd = 1,add=T)
    }#i
  }#if
#if plot.check

## ==== 10. DATA AUGMENTATION ====
y.alive <- MakeAugmentation(y = y.ar$y.ar, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
SEX <- MakeAugmentation(y = sex, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = NA)
PrevDets <- MakeAugmentation(y = PrevDets, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = NA)

## ==== 11. PREPARE NIMBLE SCR OBJECTS ====
## ====   11.1 SCALED DETECTORS ====
ScaledDetectors <- scaleCoordsToHabitatGrid(coordsData = myDetectors$detector.xy,
                                            coordsHabitatGridCenter = habitatxy,
                                            scaleToGrid =T )

## ====   11.2 GET WINDOW COORDINATES ====
ScaledLowUpCoords <- getWindowCoords(scaledHabGridCenter = ScaledDetectors$coordsHabitatGridCenterScaled,
                                     scaledObsGridCenter = ScaledDetectors$coordsDataScaled)


habitat.mx <- ScaledLowUpCoords$habitatGrid 
habitat.mx[habitat.mx[]>0]<-1

## ====   11.3 GET LOCAL OBJECTS ====
LocalDetectors <- getLocalObjects(habitatMask = habitat.mx,
                                     coords = ScaledDetectors$coordsDataScaled,
                                     dmax =  7,
                                     resizeFactor = 1,
                                     plot.check = TRUE
)
## ====   11.4 GET SPARSE OBJECTS ====
ySparse <- getSparseY(y.alive)


## ==== 12. DEFINE MODEL CODE ====
modelCode <- nimbleCode({
  ##---- SPATIAL PROCESS 
  ## Prior for AC distribution parameter
  for(i in 1:nhabCov){
    habCoeffSlope[i] ~ dunif(-10,10)
  }

  habIntensity[1:numHabWindows] <- exp(habCovs[1:numHabWindows,1:nhabCov] %*% habCoeffSlope[1:nhabCov])
  
  sumHabIntensity <- sum(habIntensity[1:numHabWindows])
  logHabIntensity[1:numHabWindows] <- log(habIntensity[1:numHabWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  ## AC distribution
  for(i in 1:M){
    sxy[i, 1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:numGridRows,1:numGridCols],
      numGridRows =  numGridRows,
      numGridCols = numGridCols
    )
  }
  
  ##---- DEMOGRAPHIC PROCESS
  ## Prior for data augmentation
  psi ~ dunif(0,1)
  probMale ~ dunif(0,1)
  teta ~ dunif(0,1)
  ## Data augmentation
  for (i in 1:M){
    z[i] ~ dbern(psi)
     indCov[i,1] ~ dbern(teta)#prev dets
     indCov[i,2] ~ dbern(probMale)#sex
     

  }
  
  ##---- DETECTION PROCESS
  ## Priors for detection parameters
  for(c in 1:nCovTraps){
    betaTraps[c] ~ dunif(-5,5)
  }
  
   sigma ~ dunif(0, 50)
  
  for(c in 1:nCounties){
    p0[c] ~ dunif(0,0.50)
  }
  

   for(c in 1:2){
   indBetas[c] ~ dunif(-5,5)
  }
  
  ## Detection process
  for (i in 1:M){
    y[i, 1:lengthYCombined] ~  dbinomLocal_normalCovsis(size = trials[1:n.traps],
                                                      p0Traps = p0[1:nCounties],
                                                      sigma = sigma,
                                                      s = sxy[i,1:2],
                                                      trapCoords = trapCoords[1:n.traps,1:2],
                                                      localTrapsIndices = trapIndex[1:n.cells,1:maxNBDets],
                                                      localTrapsNum = nTraps[1:n.cells],
                                                      resizeFactor = ResizeFactor,
                                                      habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                      indicator = z[i],
                                                      lengthYCombined = lengthYCombined,
                                                      allowNoLocal = 0,
                                                      trapCovs =  trapCov[1:n.traps,1:nCovTraps],
                                                      trapCovsIntercept =  trapIntercept[1:n.traps],
                                                      trapBetas = betaTraps[1:nCovTraps],
                                                       indCov = indCov[i,1:2],
                                                       indBetas = indBetas[1:2]
                                                      
    )
    
  }
  
  ##---- DERIVED QUANTITIES
  ## Number of individuals in the population
  N <- sum(z[1:M])
  
  #density
  dens[1:numHabWindows] <- calculateDensity(s=sxy[1:M,1:2],
                                            habitatGrid = habitatGrid[1:numGridRows,1:numGridCols],
                                            indicator = z[1:M],
                                            numWindows = numHabWindows,
                                            nIndividuals = M )
  #get N within the french border
  for(i in 1:numHabWindowsFR){
    densFr[i] <-  dens[HabWindowsFr[i]]
  }
  NFrance <- sum(densFr[1:numHabWindowsFR])
  
  
})
## ==== 13. INITIAL VALUES ====
## ====   13.1 Z ====
z <- ifelse(apply(y.alive, 1,sum)>0,1,NA)
z <- z
zinits <- ifelse(!is.na(z),NA,rbinom(sum(is.na(z)),1,0.5))

## ====   13.2 SEX ====
sexInits <- ifelse(!is.na(SEX),NA,rbinom(sum(is.na(SEX)),1,0.5))
PrevDetsInits <- ifelse(!is.na(PrevDets),NA,rbinom(sum(is.na(PrevDets)),1,0.5))

## ====   13.3 S ====
#IDENTIFY THE AUGMENTED IDS
idAugmented <- which(names(z) %in% "Augmented")

#sxy 
xy <- st_coordinates(myData.alive$myData.sp)# cbind(myData.alive$myData.sp$x, myData.alive$myData.sp$y)
colnames(xy) <- c("x","y")
ScaledCoords <- scaleCoordsToHabitatGrid(coordsData =  xy,
                                            coordsHabitatGridCenter = habitatxy,
                                            scaleToGrid =T )
myData.alivetmp <- myData.alive$myData.sp
myData.alivetmp$x <- ScaledCoords$coordsDataScaled[,1]
myData.alivetmp$y <- ScaledCoords$coordsDataScaled[,2]
 
 
#ADD A FAKE YEAR TO GET getSInits TO WORK
tmp1 <- myData.alivetmp[1:10,]
tmp1$Year <- as.numeric(tmp1$Year) +1

# GET INITIAL SXY VALUES 
sxy.init <- getSInits( AllDetections = rbind(myData.alivetmp,tmp1),
                       Id.vector = y.ar$Id.vector,
                       idAugmented = idAugmented,
                       lowerCoords = as.matrix(ScaledLowUpCoords$lowerHabCoords),
                       upperCoords = as.matrix(ScaledLowUpCoords$upperHabCoords),
                       habitatGrid = ScaledLowUpCoords$habitatGrid,
                       intensity = NULL,
                       sd = 4,
                       movementMethod = "dbernppACmovement_normal"
                       
)


# get the cells that in France to extract abundance. 
load(file.path(myVars$WD,"Data",
               paste("cellInFrance", ".RData", sep = "")))

## ====   13.4 SOTHER INITS VALUES ====
nimInits <- list( "sxy" = sxy.init[,,1],
                  "z" = zinits,
                  "sigma" = runif(1,1,1.1),
                  "habCoeffSlope" = runif(dim(habCovs)[2],-0.1,0.1),#[CM]#0,
                  "p0" = runif(max(trapIntercept),0.2,0.4),#array(runif(max(trapIntercept),0.2,0.4),c(max(trapIntercept),2)),#[CM]rep(0,dim(detCovs)[3]),
                  "betaTraps"  = runif(dim(trapCov)[2], 0.4, 0.5),
                  "psi" = runif(1,0.4,0.5),#,
                  "probMale"= runif(1,0.4,0.5),
                  "teta" = runif(1,0.4,0.5),
                  "indCov" = cbind(PrevDetsInits,sexInits),
                  "indBetas" = runif(2,0.2,0.3)
)

## ==== 14. NIMDATA ====
nimData <- list(y = ySparse$yCombined[,,1],
                  trapCov = trapCov,
                  habCovs = habCovs,
                  z = z,
                  indCov = cbind(PrevDets,SEX),
                  trapIntercept = trapIntercept
                  )

## ==== 15. NIMCONSTANTS ====
nimConstants <- list( M = dim(ySparse$yCombined)[1],
                      nhabCov= dim(nimData$habCovs)[2],
                      numHabWindows = dim(ScaledLowUpCoords$upperHabCoords)[1],
                      nCounties = max(trapIntercept),
                      habitatGrid = ScaledLowUpCoords$habitatGrid,
                      numGridRows = dim(ScaledLowUpCoords$habitatGrid)[1],
                      numGridCols = dim(ScaledLowUpCoords$habitatGrid)[2],
                      lowerHabCoords = as.matrix(ScaledLowUpCoords$lowerHabCoords),
                      upperHabCoords = as.matrix(ScaledLowUpCoords$upperHabCoords),
                      y.maxDet = dim(LocalDetectors$habitatGrid)[1],
                      x.maxDet = dim(LocalDetectors$habitatGrid)[2],
                      ResizeFactor = LocalDetectors$resizeFactor,
                      n.cells = dim(LocalDetectors$localIndices)[1],
                      maxNBDets = LocalDetectors$numLocalIndicesMax,
                      trapIndex = LocalDetectors$localIndices,
                      nTraps = LocalDetectors$numLocalIndices,
                      habitatIDDet = LocalDetectors$habitatGrid,
                      lengthYCombined = ySparse$lengthYCombined,
                      trials = myData.alive$n.trials[[1]],
                      trapCoords = ScaledDetectors$coordsDataScaled,
                      n.traps=dim(ScaledDetectors$coordsDataScaled)[1],
                      nCovTraps=dim(trapCov)[2],
                      HabWindowsFr = cellInFrance,
                      numHabWindowsFR = length(cellInFrance)
)
  
## ==== 16. BUILD AND FIT MODEL  ====
## ==== 17. NIMPARAMETERS  ====
nimParams <- c("N","sigma","psi",
               "p0", "habCoeffSlope","betaTraps","teta","probMale","indBetas")
# SAVE LESS ITERATIONS FOR AC LOCAITON AND INDIVIDUAL STATES
nimParams2 <- c("sxy","z")

## ==== 18. SAVE MODEL  ====
save(nimData,
     nimConstants,
     nimParams,
     nimParams2,
     modelCode,
     nimInits,
     file = file.path(myVars$WD,"Data","model.RData", sep = ""))


## ==== 19. BUILD AND FIT MODEL  ====
## ====   19.1 BUILD AND CACULATE  ====
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,
                      calculate = F)
model$calculate()
model$initializeInfo()
## ====   19.2 SOME CHECKS AND THE REST  ====
cmodel <- compileNimble(model)
MCMCconf <- configureMCMC(model = model,
                          monitors  = c("N","sigma","psi","p0",
                                        "habCoeffSlope","betaTraps","probMale","deltap0","teta","NFrance"),
                          control = list(reflective = TRUE),
                          monitors2= c("sxy","z","sex","indDet"),
                          thin2 = 10,
                          thin = 1)
MCMC <- buildMCMC(MCMCconf)

cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)

## Run MCMC
MCMCRuntime <- system.time(samples <- runMCMC( mcmc = cMCMC,
                                               nburnin = 100,
                                               niter = 500,
                                               nchains = 2,
                                               samplesAsCodaMCMC = TRUE))

###

myResults <- ProcessCodaOutput(samples$samples,params.omit = c("sxy","z"))

## abundance estimates
myResults$mean$NFrance
myResults$q2.5$NFrance
myResults$q97.5$NFrance


## reproduce figures ##
###map covariates###
r1 <- r2 <- r3 <- r4 <- hab.r 
r1[!is.na(hab.r)] <- nimData$habCovs[,1] 
r2[!is.na(hab.r)] <- nimData$habCovs[,2] 
r3[!is.na(hab.r)] <- nimData$habCovs[,3] 
r4[!is.na(hab.r)] <- nimData$habCovs[,4] 


box <- st_as_sfc(st_bbox(r1))

#### HABITAT 
pdf(file="C:/Personal_Cloud/OneDrive/Work/CNRS/Papers/SCRWolfFrance/Figure/FigureS1HabCov.pdf",width = 9,height = 8)
par(mfrow=c(2,2),mar=c(1,1,1,1.5))
plot(box,border=NA)
plot(gg[1,]$geometry,col=grey(0.8),add=T)
plot(r1,axes=F,box=F,col=map.pal("viridis",100),add=T)
mtext("Historical presence")

plot(box,border=NA)
plot(gg[1,]$geometry,col=grey(0.8),add=T)
plot(r2,axes=F,box=F,col=map.pal("viridis",100),add=T)
mtext("Human")

plot(box,border=NA)
plot(gg[1,]$geometry,col=grey(0.8),add=T)
plot(r3,axes=F,box=F,col=map.pal("viridis",100),add=T)
mtext("Low natural vegetation")

plot(box,border=NA)
plot(gg[1,]$geometry,col=grey(0.8),add=T)
plot(r4,axes=F,box=F,col=map.pal("viridis",100),add=T)
mtext("Forest")

dev.off()


####DETECTORS 
det.r <- myDetectors$maindetector.r
habitat.pol <- sf::st_as_sf(stars::st_as_stars(hab.r), 
                            as_points = FALSE, merge = F)
rdet1 <- det.rp0 <- rdet2 <- rdet3  <- raster(det.r) 
rdet1[!is.na(rdet1[])] <- nimData$trapCov[,1] 
rdet2[!is.na(rdet1[])]  <- nimData$trapCov[,2] 
rdet3[!is.na(rdet1[])] <- nimData$trapCov[,3] 

pdf(file="C:/Personal_Cloud/OneDrive/Work/CNRS/Papers/SCRWolfFrance/Figure/FigureS1DetCov.pdf",width = 9,height = 8)
par(mfrow=c(2,2),mar=c(1,1,1,1.5))
plot(box,border=NA)
plot(gg[1,]$geometry,col=grey(0.8),add=T)
plot(rdet1,axes=F,box=F,col=map.pal("viridis",100),add=T)
mtext("Effort")

plot(box,border=NA)
plot(gg[1,]$geometry,col=grey(0.8),add=T)
plot(rdet2,axes=F,box=F,col=map.pal("viridis",100),add=T)
mtext("Snow")

plot(box,border=NA)
plot(gg[1,]$geometry,col=grey(0.8),add=T)
plot(rdet3,axes=F,box=F,col=map.pal("viridis",100),add=T)
mtext("Roads")
dev.off()


###########p0####
## ====     20.5.5 P0 #### 
pdf(file="C:/Personal_Cloud/OneDrive/Work/CNRS/Papers/SCRWolfFrance/Figure/FigureS3P0.pdf",width = 9,height = 6)

det.rp0[!is.na(rdet1[])]  <- nimData$trapIntercept
det.pol <- sf::st_as_sf(stars::st_as_stars(det.rp0), 
                        as_points = FALSE, merge = F)
det.pol <- det.pol %>%    group_by(layer) %>%summarize()

nregions <- dim(myResults$sims.list$p0)[2]
par(mfrow=c(1,2),mar=c(5,5,0.5,1))
plot(-1000,xlim=c(0, nregions+1), ylim=c(0,0.02),ylab="p0",xlab="Region",xaxt="n")
regNames <- 0

for(c in 1:nregions){
  regNames[c] <- DetRegion[which(trapIntercept %in% c)[1]]
  
  #axis(1, at=c(1:dim(habCovs)[2]), labels = colnames(habCovs) )
  tmp <- myResults$sims.list$p0[,]
  col <- viridis::viridis(nregions)# c("red","blue","green4","purple")
  plotBars(x=tmp[,c],
           at = c,
           quantile = c(0.0275, 0.975),
           quantile1 = c(0.25, 0.75),
           widthBar = 0.25,
           col = col[c],alpha= 0.5,
           alpha1= 0.8
  )
}

labels <- Departement$code_insee[as.numeric(regNames)]
labels[labels%in% NA] <- c("A","B","C")
# axis(1, at=c(1:nregions), labels = labels )
axis(1, at=c(1:nregions), labels = c(1:nregions) )

#PLOT REGIONS 
par(mar=c(0,1,1,0))

plot(box,border=NA)
plot(gg[1,]$geometry,col=grey(0.8),add=T)
plot(det.pol$geometry,col=adjustcolor(col[det.pol$layer],alpha.f = 0.8),cex=0.3,add=T,pch=16,border=col[det.pol$layer])
text(x=st_coordinates(st_centroid(det.pol))[,1],
     y=st_coordinates(st_centroid(det.pol))[,2],det.pol$layer,col=grey(1),font=2)

dev.off()


##
####DETECTORS 
det.r <- myDetectors$maindetector.r
habitat.pol <- sf::st_as_sf(stars::st_as_stars(hab.r), 
                            as_points = FALSE, merge = F)
rdet1 <- det.rp0 <- rdet2 <- rdet3  <- raster(det.r) 
rdet1[!is.na(rdet1[])] <- nimData$trapCov[,1] 
rdet2[!is.na(rdet1[])]  <- nimData$trapCov[,2] 
rdet3[!is.na(rdet1[])] <- nimData$trapCov[,3] 

pdf(file="C:/Personal_Cloud/OneDrive/Work/CNRS/Papers/SCRWolfFrance/Figure/FigureS2DetCov.pdf",width = 9,height = 8)
par(mar=c(1,1,1,1.5))
plot(box,border=NA)
plot(gg[1,]$geometry,col=grey(0.8),add=T)
plot(rdet1,axes=F,box=F,col=map.pal("viridis",100),add=T)
mtext("Effort")

dev.off()
## ==== 20. PLOT  ====
## ====   20.1 GET MCMC/POSTERIORS   ====
c=1
load(file.path(myVars$WD,"Output", myVars$modelName,
          paste(unlist(strsplit(myVars$modelName, '_'))[1],"Chain", c, ".RData", sep = "")))

outDirectories <- list.files(file.path(myVars$WD,"Output", myVars$modelName))[grep(paste("NimbleOutFOR",unlist(strsplit(myVars$modelName, '_'))[1],sep=""), list.files(file.path(myVars$WD,"Output", myVars$modelName)))]
path.list <- file.path(myVars$WD,"Output", myVars$modelName, outDirectories)

# Retrieve the minimum number of bites per chain
numBites <- unlist(lapply(path.list, function(x){
  files <- list.files(x)
  files <- files[grep(".RData", files)]
  length(files)
}))
minBites <- min(numBites)/2
#DEFINE NUMBER OF ITERS TO REMOVE (BITE SIZE IS 250 ITERATIONS) 
NSkipBites <- 30
nthin <- 1
nimOutput <- RUNTIME <-nimOutputsxy<- list()

c=1

for(p in length(path.list):1){
  print(path.list[p])
  outfiles <- list.files(path.list[p])
  setwd(path.list[p])
  out <- outsxy <- runtime <- list()#[CM]
  for(x in NSkipBites:minBites){
    print(x)
    load(file.path( paste("bite_", x, ".RData", sep = "")))
    runtime[[x]] <- RunTime[3] 
    params.simple <- sapply(strsplit(colnames(this.sample), "\\["), "[", 1)
    parmIndex <- which(! params.simple %in% c("sxy","z"))
    nthins <- seq(1,dim(this.sample)[1], by=nthin)
    out[[x]] <- this.sample[nthins,parmIndex]#[ ,parmIndex] 
    if(sum(is.na(out[[x]]))>0){
      out[[x]] <- out[[x]][-unique(which(is.na(out[[x]]),arr.ind = T)[,1]),]
    }
    load(file.path(paste("biteSxyZ_", x, ".RData", sep = "")))
    outsxy[[x]] <- this.sampleSxyZ#[ ,parmIndex] 
    
    
  }#x
  RUNTIME[[p]] <- unlist(runtime)#[CM]
  out.mx <- do.call(rbind, out)
  nimOutput[[p]] <- as.mcmc(out.mx)
  
  out.mxsxy <- do.call(rbind, outsxy)
  nimOutputsxy[[p]] <- as.mcmc(out.mxsxy)
}#p

# COMPILE AND PROCESS THE RESULTS 
#MAIN PARAMETERS
nimOutput <- as.mcmc.list(nimOutput)
myResults <- ProcessCodaOutput(nimOutput,params.omit = c("sxy","z"))

#SXY AND Z PARAMETERS
nimOutputsxy <- as.mcmc.list(nimOutputsxy)
myResultsSZ <- ProcessCodaOutput(nimOutputsxy,params.omit = c("sxy","z","sex","indDet"))
dimnames(myResultsSZ$sims.list$sxy)[[3]] <- c("x","y")


## ====   20.2 PLOT CHAINS AND CONVERGENCE   ====
chainsPlot(nimOutput,
           file=file.path(myVars$WD, "Output",myVars$modelName,paste("Chains_",unlist(strsplit(myVars$modelName, '_'))[1],".pdf",sep=""))
           ,height = 15)

#CHECK RHAT
myResults$Rhat

## ====   20.3 UNSCALE COORDINATES FOR PLOTTING  ====
sRescaled <- scaleCoordsToHabitatGrid(coordsData = myResultsSZ$sims.list$sxy,
                                      coordsHabitatGridCenter = habitatxy,
                                      scaleToGrid = F)$coordsDataScaled

#PLOT CHECK HTAT POSTERIORS FALLS IN HABITAT
plot(habitat.r)
i=1
points(sRescaled[i,myResultsSZ$sims.list$z[i,]%in%1,2]~sRescaled[1,myResultsSZ$sims.list$z[i,]%in%1,1])


## ====   20.4 ESTIMATE DENSITY  ====
## ====      20.4.1 GET THE OBJECTS TO RUN THE DENSITY FUNCTION ====
## IDENTIFY BUFFER AREA 
plot(habitat.r)
plot(areaSearched,add=T)
habitatBuffer <- habitat.r
id <- raster::extract(habitatBuffer, areaSearched, cellnumbers=T)[[1]][,1]
habitatBuffer[id] <- 2
plot(habitatBuffer)

## SET RESOLUTION FOR EXTRACTION 
habDensity.r <- raster::disaggregate(habitatBuffer, fact=2)
habDensity.r[habDensity.r%in% 0] <- 1

## GET OBJECTS FOR DENSITY EXTRACTION 
densityInputCountries <- getDensityInput( regions = habDensity.r
                                          , 
                                          habitat = habDensity.r
                                          ,
                                          s = sRescaled
                                          ,
                                          plot.check = TRUE
)
dimnames(densityInputCountries$regions.rgmx)[[1]] <- c("Buffer","France")

## ====      20.4.2 GET AC DENSITY USING C++ ====
alive.states <- c(1)
DensityCountriesRegions <- GetDensity_PD(
      sx = densityInputCountries$sx,
      sy =  densityInputCountries$sy,
      z = myResultsSZ$sims.list$z,
      IDmx = densityInputCountries$habitat.id,
      aliveStates = alive.states,
      regionID = densityInputCountries$regions.rgmx,
      returnPosteriorCells = F)
##
##get density using the nimbleSCR function 
## this will do the same thing but with nimbleSCR
# ite <- 2
# which(habDensity.r[]>1)
# #
# dens <- calculateDensity(s=myResultsSZ$sims.list$sxy[ite,,],
#                          habitatGrid = nimConstants$habitatGrid,
#                          indicator = myResultsSZ$sims.list$z[ite,],
#                          numWindows = nimConstants$numHabWindows,
#                          nIndividuals =nimConstants$M )
# 
# plot(habDensity.r)
# plot(habitatxysf$geometry,add=T)
#cellInFrance <- which(raster::extract(habDensity.r,habitatxysf)>1)

save(cellInFrance,
     file = file.path(myVars$WD,"Output", myVars$modelName,
                      paste("cellInFrance", ".RData", sep = "")))

# sum(dens[cellInFrance])
# DensityCountriesRegions$PosteriorRegions[,ite]

ScaledDetectors$coordsHabitatGridCenterScaled
ScaledDetectors$coordsDataScaled



## CHECK SUMMARY 
DensityCountriesRegions$summary
breaks = seq(range(DensityCountriesRegions$PosteriorRegions["France",])[1],
             range(DensityCountriesRegions$PosteriorRegions["France",])[2],1)
             
quantile(DensityCountriesRegions$PosteriorRegions["France",],probs=c(0.025,0.975))

hist1 <-hist(DensityCountriesRegions$PosteriorRegions["France",],main="Distribution des effectifs hiver 2023/24",
     xlab="Effectifs",breaks	=breaks)
probs <-hist1$counts/sum(hist1$counts)
sum(probs[hist1$breaks<1014],na.rm = T)
sum(probs[hist1$breaks>1021],na.rm = T)

sum(probs[hist1$breaks%in%1022],na.rm = T)
sum(probs[hist1$breaks%in%1013],na.rm = T)


abline(v=DensityCountriesRegions$summary["France","mean"],col="red",lwd=2)
abline(v=DensityCountriesRegions$summary["France","median"],col="blue",lwd=2)
abline(v=DensityCountriesRegions$summary["France","mode"],col="orange",lwd=2)

legend("topright",legend = c("mean","median","mode"),lty=c(1,1,1),col =c("red","blue","orange"))

##prop of detected ids ##
mean(length(unique(DNA$Id))/DensityCountriesRegions$PosteriorRegions["France",])
quantile(length(unique(DNA$Id))/DensityCountriesRegions$PosteriorRegions["France",],probs=c(0.025,0.975))

## PLOT CHECK 
density.r <- habDensity.r

#DensityCountriesRegions$MeanCell[1] <- 0
density.r[!is.na(density.r[])] <- DensityCountriesRegions$MeanCell
plot(density.r)
plot(France$geometry,add=T)
## ====      20.4.3 GET AC DENSITY OF AUGMENTED IDS USING C++ ====
##get AC of augmented ids ##
augID <- which(names(nimData$z)%in% "Augmented")

DensityaugmCountriesRegions <- GetDensity_PD(
  sx = densityInputCountries$sx[,augID],
  sy =  densityInputCountries$sy[,augID],
  z = myResultsSZ$sims.list$z[,augID],
  IDmx = densityInputCountries$habitat.id,
  aliveStates = alive.states,
  regionID = densityInputCountries$regions.rgmx,
  returnPosteriorCells = F)

## CHECK SUMMARY 
DensityaugmCountriesRegions$summary

## PLOT CHECK 
density.rAugm <- habDensity.r
density.rAugm[!is.na(density.rAugm[])] <- DensityaugmCountriesRegions$MeanCell
plot(density.rAugm)
density.rAugmWthBuf <- density.rAugm
density.rAugmWthBuf <- crop(density.rAugmWthBuf, areaSearched )
density.rAugmWthBuf <- mask(density.rAugmWthBuf, France )

# density.rAugmWthBuf[habitatBuffer[]%in% 0] <- NA
plot(density.rAugmWthBuf)
plot(France$geometry,add=T)

## ====      20.4.4 GET DENSITY SPACE USED C++ ====
habitatMask <- densityInputCountries$habitat.id
habitatMask[!is.na(habitatMask)] <- 1
#RESCALE SIGMA TO METERS
sigma <- myResults$sims.list$sigma*res(habitat.r)[1]
#RESCALE SIGMA TO THE HABITAT SCALE
sigmaRescaled <- sigma/res(densityInputCountries$regions.r)[1]
thinnediter <- seq(1,length(sigmaRescaled)[1],
                   by=length(sigmaRescaled)[1]/dim(myResultsSZ$sims.list$z)[1])
sigmaRescaled <- sigmaRescaled[thinnediter]

### COMPUTE THE UD BASED FOR A FEW ITERATIONS
iter <- sample(1:dim(densityInputCountries$sy)[1], 50)

# BECAUSE SIGMA IS ESTIMATED SEPARATELY FOR MALES/FEMALES/STATE, WE HAVE TO DO IT SEPARATELY FOR EACH SEX AND STATE.
spaceUSED <- GetSpaceUse( sx = densityInputCountries$sx[iter,],
                          sy =  densityInputCountries$sy[iter,],
                          z = myResultsSZ$sims.list$z[iter,],
                          sigma = sigmaRescaled[iter],#,"F"],
                          densityInputCountries$habitat.xy,
                          aliveStates = alive.states,
                          regionID = densityInputCountries$regions.rgmx,
                          display_progress = T,
                          returnPosteriorCells = T)
#PLOT CHECK 
SpaceUsed.r <- habDensity.r
spaceUSED$MeanCell[1] <- 0
SpaceUsed.r[!is.na(SpaceUsed.r[])] <- spaceUSED$MeanCell
plot(SpaceUsed.r)
plot(France$geometry,add=T)
mapview::mapview(SpaceUsed.r)
mapview::mapview(list(myData.alive$myData.sp))
mapview::mapview(LCIEDistribution)

#################################


myResults$mean$indBetas
myResults$q2.5$indBetas
myResults$q97.5$indBetas

##############################

## ====   20.5  SUMMARY OF THE RESULTS VIOLINS ####
pdf(file=file.path(myVars$WD,"output",myVars$modelName,paste("SummaryResults_",paste(unlist(strsplit(myVars$modelName, '_')))[1],
                                                             ".pdf",sep="")))
## ====     20.5.1 POP SIZE #### 
plot(-1000,xlim=c(0,2), ylim=c(0,2000),ylab="Population size",xlab="",xaxt="n")
t=1
widthPolygon <- 0.1
tmp <- DensityCountriesRegions$PosteriorRegions["France",]
#N
plotBars(x=tmp,
         at = 1,
         quantile = c(0.0275, 0.975),
         quantile1 = c(0.25, 0.75),
         widthBar = 0.1,
         col = "red",
         alpha= 0.5,
         alpha1= 0.8
         )
#ADD THE VALUE TO THE PLOT 
est <- DensityCountriesRegions$summary["France",]
text(1.2,250,paste(est["95%CILow"],"-",est["95%CIHigh"], "\n",
                   " median=",est["median"],
                   " mean=", round(est["mean"],digits=2)))


## ====     20.5.2 SIGMA #### 
par(mar=c(5,5,5,5))
plot(-1000,xlim=c(0,3), ylim=c(0,5000),ylab="Sigma (m)",xlab="Sex",xaxt="n")
axis(1,at=c(1:2),labels = c("F","M"))
tmp <- myResults$sims.list$sigma*res(habitat.r)[1]
plotBars(x=tmp,
         at = 1,
         quantile = c(0.0275, 0.975),
         quantile1 = c(0.25, 0.75),
         widthBar = 0.1,
         col = "red",alpha= 0.5,
         alpha1= 0.8
)

## ====     20.5.3 DETECTOR COVARIATE #### 
par(mar=c(5,5,5,5))
plot(-1000,xlim=c(0,dim(trapCov)[2]+1), ylim=c(-1,1),ylab="BetaTraps",xlab="",xaxt="n")
abline(h=0)
axis(1,at=c(1:dim(trapCov)[2]),labels = colnames(trapCov))
tmp <- myResults$sims.list$betaTraps
col <- c("red","blue","green4")
for(c in 1:dim(trapCov)[2]){
  plotBars(x=tmp[,c],
           at = c,
           quantile = c(0.0275, 0.975),
           quantile1 = c(0.25, 0.75),
           widthBar = 0.1,
           col = col[c],alpha= 0.5,
           alpha1= 0.8
  )
}

## ====     20.5.4 DENSITY COVARIATE #### 
par(mar=c(5,5,5,5))
plot(-1000,xlim=c(0, dim(habCovs)[2]+1), ylim=c(-2,2),ylab="BetaHab",xlab="",xaxt="n")
abline(h=0)


axis(1, at=c(1:dim(habCovs)[2]), labels = colnames(habCovs) )
tmp <- myResults$sims.list$habCoeffSlope
col <- c("red","blue","green4","purple")
for(c in 1:dim(habCovs)[2]){
  plotBars(x=tmp[,c],
           at = c,
           quantile = c(0.0275, 0.975),
           quantile1 = c(0.25, 0.75),
           widthBar = 0.1,
           col = col[c],alpha= 0.5,
           alpha1= 0.8
  )
}

## ====     20.5.5 P0 #### 
nregions <- dim(myResults$sims.list$p0)[2]
par(mfrow=c(1,2),mar=c(5,5,0.5,1))
plot(-1000,xlim=c(0, nregions+1), ylim=c(0,0.06),ylab="p0",xlab="Region",xaxt="n")
regNames <- 0
for(c in 1:nregions){
  regNames[c] <- DetRegion[which(trapIntercept %in% c)[1]]
  
  #axis(1, at=c(1:dim(habCovs)[2]), labels = colnames(habCovs) )
  tmp <- myResults$sims.list$p0[,]
  col <- viridis::viridis(nregions)# c("red","blue","green4","purple")
  plotBars(x=tmp[,c],
               at = c,
               quantile = c(0.0275, 0.975),
               quantile1 = c(0.25, 0.75),
               widthBar = 0.2,
               col = col[c],alpha= 0.5,
               alpha1= 0.8
  )
}

labels <- Departement$code_insee[as.numeric(regNames)]
labels[labels%in% NA] <- c("A","B","C")
axis(1, at=c(1:nregions), labels = labels )

#PLOT REGIONS 
plot(st_geometry(habitatExtent))
plot(st_geometry(France),add=T)
plot(myDetectors$main.detector.sf[,]$geometry,col=col[trapIntercept],cex=0.3,add=T,pch=16)

  
## ====     20.5.6 PROB DETECTED  #### 
par(mfrow=c(1,1),mar=c(5,5,0.5,1))
plot(-1000,xlim=c(0, 2), ylim=c(0,1),
     ylab="teta",
     xlab="",
     xaxt="n")
tmp <- myResults$sims.list$teta
col <- c("red","blue","green4")
plotBars(x=tmp,
         at = 1,
         quantile = c(0.0275, 0.975),
         quantile1 = c(0.25, 0.75),
         widthBar = 0.2,
         col = col[1],alpha= 0.5,
         alpha1= 0.8
)
## ====     20.5.7 EFFECT PREV DETE ON DETECTION  #### 
plot(-1000,xlim=c(0, 2), ylim=c(-1,1),
     ylab="IndCov(PrevCovs)",
     xlab="",
     xaxt="n")
abline(h=0)
tmp <- myResults$sims.list$indBetas
col <- c("red","blue","green4")
plotBars(x=tmp[,1],
         at = 1,
         quantile = c(0.0275, 0.975),
         quantile1 = c(0.25, 0.75),
         widthBar = 0.2,
         col = col[1],alpha= 0.5
)

## ====     20.5.8 PROB MALE #### 
par(mfrow=c(1,1),mar=c(5,5,0.5,1))
plot(-1000,xlim=c(0, 2), ylim=c(0,1),
     ylab="probMale",
     xlab="",
     xaxt="n")
tmp <- myResults$sims.list$probMale
col <- c("red","blue","green4")
plotBars(x=tmp,
         at = 1,
         quantile = c(0.0275, 0.975),
         quantile1 = c(0.25, 0.75),
         widthBar = 0.2,
         col = col[1],alpha= 0.5,
         alpha1= 0.8
)
## ====     20.5.9 EFFECT OF SEX ON P0#### 
plot(-1000,xlim=c(0, 2), ylim=c(-1,1),
     ylab="IndCov(sex)",
     xlab="",
     xaxt="n")
abline(h=0)
tmp <- myResults$sims.list$indBetas
col <- c("red","blue","green4")
plotBars(x=tmp[,2],
         at = 1,
         quantile = c(0.0275, 0.975),
         quantile1 = c(0.25, 0.75),
         widthBar = 0.2,
         col = col[1],alpha= 0.5
)


#### PLOT DENSITY MAP ###
plot(SpaceUsed.r)
plot(France$geometry,add=T)

dev.off()
######
myResults$q2.5$habCoeffSlope
myResults$q97.5$habCoeffSlope

myResults$q2.5$indBetas
myResults$q97.5$indBetas

myResults$q2.5$betaTraps
myResults$q97.5$betaTraps


myResults$q2.5$p0
myResults$q97.5$p0


mapview::mapview(SpaceUsed.r)
mapview::mapview(density.r)
# mapview::mapview(list(density.r,DNA))
mapview::mapview(list(Departement,DNA))


#chose a region 
reg <- 1
ite <- 1
indCov <- c(1,1)#prev det ; sex
idhab <- 1540
s <- nimConstants$lowerHabCoords[idhab,]
plot(nimConstants$lowerHabCoords[,2]~nimConstants$lowerHabCoords[,1])
points(nimConstants$lowerHabCoords[idhab,2]~nimConstants$lowerHabCoords[idhab,1],col="red",pch=16)

det <- 1
iter <- sample(dim(myResults$sims.list$betaTraps)[1],500) 
s <- myResultsSZ$sims.list$sxy
z <- myResultsSZ$sims.list$z

sAliveList <-lapply(1:length(iter), function(x){
  s[x,z[x,]==1,]
})
maxAlive <-max(unlist(lapply(sAliveList, dim)))

p <- array(NA,c(dim(nimData$trapCov)[1], length(iter)))#,maxAlive ))
 p0Trap <- array(NA,c(dim(nimData$trapCov)[1], length(iter) ))

                     
pdet <- numeric()
#pdet <- array(NA, c(length(iter),maxAlive ))

count <- 0
#average trap conditions 
trapCovMean <- colMeans(nimData$trapCov)

for(ite in iter){
  count <- count +1

    
for(det in 1:dim(nimData$trapCov)[1]){
  p0Traps <- myResults$sims.list$p0[ite,nimData$trapIntercept[det]]
  trapBetas <- myResults$sims.list$betaTraps[ite,]
  indBetas <- myResults$sims.list$indBetas[ite,]
  alpha <- -1/(2 * myResults$sims.list$sigma[ite] * myResults$sims.list$sigma[ite] )
    
  p0Trap[det,count] <- ilogit(logit(p0Traps) + 
                     # inprod(trapBetas, nimData$trapCov[det,])+
                       #asssume "average" spatial conditions 
                     inprod(trapBetas, trapCovMean)+#[det,])+
                     
                    inprod(indBetas, indCov))

}
#   ss <- sAliveList[[count]]
# for(s1 in 1:dim(ss)[1]){
  # d2 <- pow(nimConstants$trapCoords[, 1] - ss[s1,1],2) + 
  #   pow(nimConstants$trapCoords[, 2] - ss[s1,2], 2)
  d2 <- pow(nimConstants$trapCoords[, 1] - s[1],2) +
    pow(nimConstants$trapCoords[, 2] - s[2], 2)
  p[,count] <- p0Trap[,count] * exp(alpha * d2)
  
  pdet[count] <- 1-prod( (1-p[,count])^nimConstants$trials[1] )
# }
}

mean(pdet[],na.rm=T)
#0.5955782
quantile(pdet,probs = c(0.0275,0.975),na.rm=T)
# 2.75%     97.5% 
#   0.4956927 0.7075147 

#0.7397374
# 2.75%     97.5% 
#   0.6094860 0.8628602

mean(pdet)
quantile(pdet,probs = c(0.0275,0.975))
#0.421468
#0.3198269 0.5375026


detPoly <- myDetectors$grid.poly 
detPoly$p0Trap <- rowMeans(p0Trap)
detPoly$p <- rowMeans(p)

plot(detPoly["p0Trap"])
plot(detPoly["p"])

nimData$trapCov


pZero <- ilogit(logit(p0) + betaTraps* covTraps)
for(i in 1:length(myHabitat$habitat.sp)){
  d2 <- (traplocs[,1]-Sgrid[i,1])^2 + (traplocs[,2]-Sgrid[i,2])^2
  pvec <-  pZero*exp(-(1/(2*sigma*sigma))*d2)
  netp[i] <- 1-prod( (1-pvec)^myDetectors$main.detector.sp$count )  # prob of being caught in any trap
}



###   ====      15.GOF ####
iter <- 150
#RECREATE THE MODEL TO SIMULATE DATA FROM THE POSTERIORS 
nimData2 <- nimData
nimData2$sxy <- myResultsSZ$sims.list$sxy[iter,,]
nimData2$z <- myResultsSZ$sims.list$z[iter,]
nimData2$p0 <- myResults$sims.list$p0[iter,,]
nimData2$psi <- myResults$sims.list$psi[iter]
nimData2$sigma <- myResults$sims.list$sigma[iter,]
nimData2$betaTraps <- myResults$sims.list$betaTraps[iter,]
nimData2$habCoeffSlope <- myResults$sims.list$habCoeffSlope[iter,]
nimData2$deltap0 <- myResults$sims.list$deltap0[iter,]
nimData2$teta <- myResults$sims.list$teta[iter,]
nimData2$probMale <- myResults$sims.list$probMale[iter,]

lengthYCombined <- nimConstants$lengthYCombined 
nimConstants$lengthYCombined <-100

Ytrue <-  nimData$y 
nimData2$y <- NULL

model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData2,
                      check = F,
                      calculate = F)
cmodel <- compileNimble(model)

cmodel$simulate()

#COMPILE FUNCTION TO CALCULATE EXPECTED P 
source(file.path("C:/My_documents/rovquant/analyses/Rgit/RovQuant/","Temp/CM/functions/Nimble/calculateP.R"))

calculatePc <- compileNimble(calculateP)

#INITIALIZE OBJECTS
gof_obs_ik <- gof_sim_ik <- gof_obs_i <- gof_sim_i <- 
  gof_obs_j <- gof_sim_j <- 0

# SIMULATE DATA FOR xx ITERATIONS 
niter <- 200
itera <- sample(dim(myResultsSZ$sims.list$sxy)[1], niter)
for(iter in 1:niter){
  print(iter)
  #FETCH THE POSTEIROS 
  cmodel$sxy <- myResultsSZ$sims.list$sxy[itera[iter],,]
  cmodel$z <- myResultsSZ$sims.list$z[itera[iter],]
  cmodel$p0 <- myResults$sims.list$p0[itera[iter],,]
  cmodel$psi <- myResults$sims.list$psi[itera[iter]]
  cmodel$sigma <- myResults$sims.list$sigma[itera[iter],]
  cmodel$betaTraps <- myResults$sims.list$betaTraps[itera[iter],]
  cmodel$habCoeffSlope <- myResults$sims.list$habCoeffSlope[itera[iter],]
  cmodel$deltap0 <- myResults$sims.list$deltap0[itera[iter]]
  cmodel$teta <- myResults$sims.list$teta[itera[iter]]
  cmodel$probMale <- myResults$sims.list$probMale[itera[iter]]

  #Simulate
  cmodel$simulate()
  
  ########################
  ypred <- y <- p <- matrix(0,nrow=dim(cmodel$y)[1],
                            ncol = nimConstants$n.traps)
  #reconstruct y #
  
  cmodel$y
  i=4
  for(i in 1:nimConstants$M){
    detNums <- cmodel$y[i,1]
    if(detNums>0){
      nMaxDetectors <- (nimConstants$lengthYCombined - 1)/2
      #GET DET INDICES
      detIndices1 <- cmodel$y[i,(nMaxDetectors + 2):nimConstants$lengthYCombined]
      detIndices1 <- detIndices1[1:detNums]
      # GET NUMBER OF DETECTIONS
      x1 <- cmodel$y[i,2:(nMaxDetectors + 1)]
      x1 <- x1[1:detNums]
      
      ypred[i,detIndices1] <- x1
      
    }
  }
  
  sum(ypred)
  sum(cmodel$y[,1])
  
  for(i in 1:nimConstants$M){
    detNums <- nimData$y[i,1]
    if(detNums>0){
      nMaxDetectors <- (lengthYCombined - 1)/2
      #GET DET INDICES
      detIndices1 <- nimData$y[i,(nMaxDetectors + 2):lengthYCombined]
      detIndices1 <- detIndices1[1:detNums]
      # GET NUMBER OF DETECTIONS
      x1 <- nimData$y[i,2:(nMaxDetectors + 1)]
      x1 <- x1[1:detNums]
      
      y[i,detIndices1] <- x1
      
    }
  }
  
  sum(y)
  sum(nimData$y[,1])
  p <- calculatePc(  s = cmodel$sxy
                     , trapCoords = nimConstants$trapCoords
                     , betaTraps = cmodel$betaTraps
                     , p0 = cmodel$p0
                     , z = cmodel$z
                     , sigma = cmodel$sigma)
  
  # alpha <- -1/(2 * cmodel$sigma[] * cmodel$sigma[])
  # 
  # for(i in 1:nimConstants$M){
  #   if(cmodel$z[i]==0){next}
  #   s <- cmodel$sxy[i,]
  #   d2 <- pow(nimConstants$trapCoords[, 1] - s[1], 2) + pow(nimConstants$trapCoords[, 2] - s[2], 2)
  #   print(i)
  #     for(j in 1:nimConstants$n.traps){
  #      p0Trap <- ilogit(logit(cmodel$p0[1]) + inprod(cmodel$trapCov[j,], cmodel$betaTraps))
  #    }
  #  p[i,] <- p0Trap * exp(alpha * d2[])
  # 
  # }
  
  ##subset to detected ID 
  pexpected <- p[cmodel$z[]>0,]
  dim(pexpected)
  
  # Expected dets 
  EY <- pexpected#apply(test_stats$expected, c(1,3), sum)
  # Y true observations#
  OYobs <- y[cmodel$z[]>0,] # rbind(apply(tmp_object$capthist,c(1,3),sum),
  #      matrix(0,ntot-nind,nrow(traps(tmp_object$capthist))))
  # Simulated Y 
  OYsim <- ypred[cmodel$z[]>0,]#apply(newCH,c(1,3),sum)
  
  #
  gof_obs_ik[iter] <- FT_stat(OYobs, EY)
  gof_sim_ik[iter] <- FT_stat(OYsim, EY)
  
  ## Test: individual
  EY2 <- apply(pexpected,1,sum)
  OYobs2 <- apply(OYobs,1,sum)
  OYsim2 <- apply(OYsim,1,sum)
  gof_obs_i[iter] <- FT_stat(OYobs2, EY2)
  gof_sim_i[iter] <- FT_stat(OYsim2, EY2)
  
  ## Test: traps
  EY3 <-  apply(pexpected,2,sum)
  OYobs3 <- apply(OYobs,2,sum)
  OYsim3 <- apply(OYsim,2,sum)
  gof_obs_j[iter] <- FT_stat(OYobs3, EY3)
  gof_sim_j[iter] <- FT_stat(OYsim3, EY3)
  
}

pdf(file=file.path(myVars$WD,"output",myVars$modelName,"GOF.pdf"))
par(mfrow = c(2, 2), oma = c(0, 0, 0, 0))
mm <- c(min(c(gof_sim_ik, gof_obs_ik)), max(c(gof_sim_ik, 
                                              gof_obs_ik)))
clr <- ifelse(gof_sim_ik > gof_obs_ik, adjustcolor("navy", 
                                                   0.5), adjustcolor("darkgreen", 0.5))
plot(gof_sim_ik, gof_obs_ik, pch = 16, col = clr, las = 1, 
     asp = 1, xlim = c(mm[1], mm[2]), ylim = c(mm[1], mm[2]), 
     main = paste("pr(sim > obs):", round(mean(gof_sim_ik > 
                                                 gof_obs_ik), 2)))
abline(0, 1, col = 2, lwd = 2)
mm <- c(min(c(gof_sim_i, gof_obs_i)), max(c(gof_sim_i, gof_obs_i)))
clr <- ifelse(gof_sim_i > gof_obs_i, adjustcolor("navy", 
                                                 0.5), adjustcolor("darkgreen", 0.5))
plot(gof_sim_i, gof_obs_i, pch = 16, col = clr, las = 1, 
     asp = 1, xlim = c(mm[1], mm[2]), ylim = c(mm[1], mm[2]), 
     main = paste("pr(sim > obs):", round(mean(gof_sim_i > 
                                                 gof_obs_i), 2)))
abline(0, 1, col = 2, lwd = 2)
mm <- c(min(c(gof_sim_j, gof_obs_j)), max(c(gof_sim_j, gof_obs_j)))
clr <- ifelse(gof_sim_j > gof_obs_j, adjustcolor("navy", 
                                                 0.5), adjustcolor("darkgreen", 0.5))
plot(gof_sim_j, gof_obs_j, pch = 16, col = clr, las = 1, 
     asp = 1, xlim = c(mm[1], mm[2]), ylim = c(mm[1], mm[2]), 
     main = paste("pr(sim > obs):", round(mean(gof_sim_j > 
                                                 gof_obs_j), 2)))
abline(0, 1, col = 2, lwd = 2)
scr_gof1nimbleSCR <- 
  list(scrgof_pval = c(`FT-ind-trap` = mean(gof_sim_ik > 
                                              gof_obs_ik), `FT-individuals` = mean(gof_sim_i > gof_obs_i), 
                       `FT-traps` = mean(gof_sim_j > gof_obs_j)),
       gof_ik = data.frame(simulated = gof_sim_ik, 
                           observed = gof_obs_ik), 
       gof_i = data.frame(simulated = gof_sim_i, observed = gof_obs_i), 
       gof_j = data.frame(simulated = gof_sim_j, observed = gof_obs_j))

dev.off()

####
scr_gof <- scrgof(model.0, nsamp = niter)

scr_gof$gof_ik$observed


###try close SCR
library(Rcapture) # https://www.jstatsoft.org/article/view/v019i05
wolf <- table(DNA$Id,DNA$mois)
wolf <- wolf[,c(4,5,1,2,3)]
wolf[wolf>0] <- 1

closedp.0(wolf)
closedp.t(wolf)

