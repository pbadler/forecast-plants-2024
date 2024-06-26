# PBA March 2016
# modfied by ARK July 2016
# call from removal_analysis_wrapper.r

#########################################
#  1. Import data and calculate W's
#########################################

doSpp <- "POSE"
sppList <- c("ARTR","HECO","POSE","PSSP","allcov","allpts")
dataDir1 <- paste("lib","/driversdata/data/idaho",sep="")
dataDir2 <- paste("lib","/driversdata/data/idaho_modern",sep="")
nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 

# set up distance weights------------------------------------------------
#dists <- read.csv('~/driversdata/data/idaho/speciesData/IdahoDistanceWeights.csv')
dists <- read.csv(paste(dataDir2,"/speciesData/IdahoModDistanceWeights_noExptl.csv",sep=""));
dists$allcov <- rowMeans(dists[,1:4])  # for "other" polygons use average of big 4
dists$allpts <- dists$POSE  # set forb dist wts = smallest grass (POSE)

# import old data--------------------------------------------------------

D1 <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir1,distWts=dists)
D1$Treatment <- "Control"
D1$Period <- "Historical"

# import modern data--------------------------------------------------------

D2 <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir2,distWts=dists)
D2$Period <- "Modern"

# merge in treatment data
tmp <- read.csv(paste(dataDir2,"/quad_info.csv",sep=""))
tmp <- tmp[,c("quad","Treatment")]
D2 <- merge(D2,tmp, all.x=T)

# account for removal in baseline years
if(doSpp!="ARTR"){
  ii <- which(D2$year>=2011 & D2$Treatment=="No_shrub")
  D2$W.ARTR[ii] <- 0
}else{
  ii <- which(D2$year>=2011 & D2$Treatment=="No_grass")
   D2$W.HECO[ii] <- 0 ; D2$W.POSE[ii] <- 0 ; D2$W.PSSP[ii] <- 0
}

# combine old and modern
allD <- rbind(D1,D2)
rm(D1,D2,tmp)

# merge data on removals at individual level
tmp <- read.csv(paste(dataDir2,"/speciesData/",doSpp,"/",doSpp,"_within_ARTRremovals.csv",sep=""))
tmp <- tmp[,c("quad","year","trackID","inARTR")] 
allD<-merge(allD,tmp,all.x=T)
allD$inARTR[is.na(allD$inARTR)] <- 0

# clean up dataset ----------------------------------------------
allD$year[allD$year<2000] <- allD$year[allD$year<2000] + 1900

if(doSpp=="ARTR"){
  keep <- which(is.element(allD$Treatment,c("Control","No_grass", "Irrigation", "Drought")))
}else{
  keep <- which(is.element(allD$Treatment,c("Control","No_shrub", "Irrigation", "Drought")))
}
allD <- allD[keep,]

# remove outliers (large plants that obviously do not turn into tiny plants)
allD <- subset(allD, ! ( quad == 'Q11' & trackID == '26' & year %in% c(1931,1932, 1938, 1939) )) # drop these observations 


tmp <- which( allD$Grazing == 'G')
plot(data = allD[-tmp, ], logarea.t1 ~ logarea.t0)
plot(data = allD[tmp, ], logarea.t1 ~ logarea.t0)

summary(lm(data=allD[-tmp, ], logarea.t1 ~ logarea.t0))
summary(lm(data=allD[tmp, ], logarea.t1 ~ logarea.t0))

plot(resid(lm(data=allD[-tmp, ], logarea.t1 ~ logarea.t0)))
plot(resid(lm(data=allD[tmp, ], logarea.t1 ~ logarea.t0)))




# set up indicator variables
allD$Treatment2 <- allD$Treatment
allD$Treatment2[allD$year>2000] <- "Modern"
allD$Treatment3 <- allD$Treatment
allD$Treatment3[allD$Treatment=="Control" & allD$year>2000] <- "ControlModern"
#allD$Treatment[ allD$year < 2012 & allD$Treatment %in% c('Drought', 'Irrigation') ] <- 'Control'  # set initial treatment to control

# drop removal treatments
tmp <- which(allD$Treatment=="No_shrub" | allD$Treatment=="No_grass")
if(length(tmp)>0) allD <- allD[-tmp,]

# ----------- use this data for prediction ------------------------------------------------------------------------------

write.csv(allD, 'data/temp/POSE_growth.csv') 

# -----------------------------------------------------------------------------------------------------------------------
