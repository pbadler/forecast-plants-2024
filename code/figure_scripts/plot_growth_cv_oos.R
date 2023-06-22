
# import results
cv_dat <- read.csv("output/growth_models/cross_validation_growth_models.csv",header=T)
oos_dat <- read.csv("output/growth_models/oos_growth_validation.csv",header=T)

#choose error metric and filter cross-validation data
my_err <- "RMSE"
tmp <- grep(my_err,names(cv_dat))
cv_dat <- cv_dat[,c(1:5,tmp)]
cv_dat <- cv_dat[,c(1:5,8,7,6)]
cv_dat_large <- cv_dat[cv_dat$size_class=="large",]

#choose error metric and filter OOS data
tmp <- c("size_class","Treatment","name",my_err,"species")
oos_dat <- oos_dat[,tmp]
oos_all_Lg <- subset(oos_dat,size_class=="large") # no small plants, all treatments averaged together
oos_all_Lg <- oos_all_Lg %>% pivot_wider(names_from = name, values_from = all_of(my_err))
oos_all_Lg <- oos_all_Lg[,c(1:3,6,5,4)]
names(oos_all_Lg)[4:6] <- names(cv_dat_large)[6:8]

# compare CV error and OOS error by species 
# focus on overall OOS error across all treatments
mycols <- c("black","lightgrey")
sppList <- sort(unique(cv_dat_large$species))

png("figures/growth_cv_oos.png",height=7, width=7,res=400,units="in")

par(mfrow=c(2,2),mar=c(3,2,2,1),oma=c(0,2,0,0))

for(iSpp in 1:length(sppList)){
  doSpp <- sppList[iSpp]
  tmpD <- rbind(cv_dat_large[cv_dat_large$species==doSpp,6:8],
                oos_all_Lg[oos_all_Lg$species==doSpp & 
                             oos_all_Lg$Treatment=="Overall",4:6])
  barplot(as.matrix(tmpD),beside=T,col=mycols,ylab="",
          names.arg=c("Null","Climate","Clim x Size"))
  if(iSpp==1) legend("topleft",c("Training","Testing"),fill=mycols,bty="n")
  mtext(doSpp,side=3,line=1,adj=0)
}
mtext("RMSE",side=2,line=1, outer=T)


dev.off()


# OOS error by treatment

mycols <- c("black","lightgrey","white")
sppList <- sort(unique(cv_dat_large$species))

png("figures/growth_oos_by_trt.png",height=7, width=7,res=400,units="in")

par(mfrow=c(2,2),mar=c(3,2,2,1),oma=c(0,2,0,0))

for(iSpp in 1:length(sppList)){
  doSpp <- sppList[iSpp]
  tmpD <- subset(oos_all_Lg,species==doSpp & Treatment!="Overall")
  barplot(t(as.matrix(tmpD[,4:6])),beside=T,col=mycols,ylab="",
          names.arg=c("Control","Drought","Irrigation"))
  if(iSpp==1) legend("topleft",c("Null","Climate","Clim x Size"),fill=mycols,bty="n")
  mtext(doSpp,side=3,line=1,adj=0)
}
mtext("RMSE",side=2,line=1, outer=T)

dev.off()


# barplots for cross validation error (all species together)
mycols=c("white","grey","black")
barplot(as.matrix(t(cv_dat_large[1:4,6:8])),beside=T,names.arg=c(cv_dat_large$species),
        ylab="RMSE",col=mycols)
legend("topleft",c("Null","Climate","ClimxSize"),fill=mycols,bty="n")
