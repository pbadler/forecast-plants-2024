
# import results
sppList <- c("ARTR","HECO","POSE","PSSP")

oos_dat <- NULL

for(iSpp in sppList){
  tmpD <- read.csv(paste0("output/",iSpp,"_cover_validation.csv"),header=T)
  oos_dat <- rbind(oos_dat,tmpD)
}
# remove duplicate entries
tmp <- which(oos_dat$Period=="Training" & oos_dat$Treatment=="Control")
oos_dat <- oos_dat[-tmp,]

#choose error metric and filter OOS data
my_err <- "RMSE"
tmp <- c("Period","Treatment","name",my_err,"species")
oos_dat <- oos_dat[,tmp]
oos_dat <- oos_dat %>% pivot_wider(names_from = name, values_from = all_of(my_err))
oos_dat <- oos_dat[,c(1:3,6,5,4)]

# compare training and OOS error by species 
mycols <- c("black","grey40","lightgrey","white")

png("figures/cover_validation.png",height=7, width=7,res=400,units="in")

par(mfrow=c(2,2),mar=c(3,2,2,1),oma=c(0,2,0,0))

for(iSpp in 1:length(sppList)){
  doSpp <- sppList[iSpp]
  tmpD <- oos_dat[oos_dat$species==doSpp,]
  tmpD <- tmpD[c(4,1:3),] # put training first
  barplot(as.matrix(tmpD[,4:6]),beside=T,col=mycols,ylab="",
          names.arg=c("Null","Climate","Clim x Size"))
  if(iSpp==1) legend("topleft",c("Training","Control","Drought","Irrigation"),fill=mycols,bty="n")
  mtext(doSpp,side=3,line=1,adj=0)
}
mtext("RMSE",side=2,line=1, outer=T)

dev.off()


