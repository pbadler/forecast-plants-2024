
# import results
cv_dat <- read.csv("output/survival_models/cross_validation_survival_models.csv",header=T)
oos_dat <- read.csv("output/survival_models/oos_survival_validation.csv",header=T)

#choose error metric and filter cross-validation data
my_err <- "RMSE"
tmp <- grep(my_err,names(cv_dat))
cv_dat <- cv_dat[,c(1:4,tmp)]
cv_dat <- cv_dat[,c(1:4,7,6,5)]

#choose error metric and filter OOS data
tmp <- c("Treatment","name",my_err,"species")
oos_dat <- oos_dat[,tmp]
oos_all <- oos_dat %>% pivot_wider(names_from = name, values_from = all_of(my_err))
oos_all <- oos_all[,c(1:2,5,4,3)]
names(oos_all)[3:5] <- names(cv_dat)[5:7]

# compare CV error and OOS error by species 
# focus on overall OOS error across all treatments
mycols <- c("black","lightgrey")
sppList <- sort(unique(cv_dat$species))

png("figures/survival_cv_oos.png",height=7, width=7,res=400,units="in")

par(mfrow=c(2,2),mar=c(3,2,2,1),oma=c(0,2,0,0))

for(iSpp in 1:length(sppList)){
  doSpp <- sppList[iSpp]
  tmpD <- rbind(cv_dat[cv_dat$species==doSpp,5:7],
                oos_all[oos_all$species==doSpp & 
                             oos_all$Treatment=="Overall",3:5])
  barplot(as.matrix(tmpD),beside=T,col=mycols,ylab="", ylim=c(0,0.38),
          names.arg=c("Null","Climate","Clim x Size"))
  if(iSpp==1) legend("topleft",c("Training","Testing"),fill=mycols,bty="n")
  mtext(doSpp,side=3,line=1,adj=0)
}
mtext("RMSE",side=2,line=1, outer=T)


dev.off()


# OOS error by treatment

mycols <- c("black","lightgrey","white")
sppList <- sort(unique(cv_dat$species))

png("figures/survival_oos_by_trt.png",height=7, width=7,res=400,units="in")

par(mfrow=c(2,2),mar=c(3,2,2,1),oma=c(0,2,0,0))

for(iSpp in 1:length(sppList)){
  doSpp <- sppList[iSpp]
  tmpD <- subset(oos_all,species==doSpp & Treatment!="Overall")
  barplot(t(as.matrix(tmpD[,3:5])),beside=T,col=mycols,ylab="",
          names.arg=c("Control","Drought","Irrigation"))
  if(iSpp==1) legend("top",c("Null","Climate","Clim x Size"),fill=mycols,bty="n")
  mtext(doSpp,side=3,line=1,adj=0)
}
mtext("RMSE",side=2,line=1, outer=T)

dev.off()



