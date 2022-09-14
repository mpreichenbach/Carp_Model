################################################################################################
# To import and model the movement of acoustic-tagged Silver Carp in ponds exposed to various Sound Stimuli
# Date November 26 2019
################################################################################################
#packages
library(momentuHMM)
library(rgdal)
library(animation)
library(lubridate)
library(circular)
library(runjags)
library(tidyverse)
library(StreamMetabolism) 
library(mvtnorm)
library(plyr)
source("list.to.array.R")

# package_path <- "C:/Users/jplumb/Desktop/JAGS-VONMISES-MODULE-0.2.tar"
# install.packages(package_path, repos=NULL, type="source")
# load.module("vonmises")
#library(Rcmdr)
#Map to the animation directory
#ani.options(outdir = getwd(), convert = 'C:/Program Files/ImageMagick-7.0.9-Q16/convert.exe') 

#Set working directory
setwd("C:/Users/jplumb/Desktop/Carp Pond Analysis")
source("list.to.array.R")
#Please specify the trial, pond number, and scenario
TrialNum = 1
PondNum = 26
Senario = "Control" #"ChirpSaw"

#-------------------------------------------------------------------------------------
#assign file names
if (PondNum == 31) {
  filename = 'TagHistoryFor_AllTags-Trial1-System100_ChirpSaw_DroppedTagsRemoved.csv'    #<===Pond 31
} else if (PondNum == 30) {
  filename = 'TagHistoryFor_AllTags-Trial1-System200_BoatM_DroppedTagsRemoved.csv'       #<===Pond 30
} else if (PondNum == 27) {
  filename = 'TagHistoryFor_AllTags-Trial1-System300_ChirpSquare_DroppedTagsRemoved.csv' #<===Pond 27
} else if (PondNum == 26) {
  filename = 'TagHistoryFor_AllTags-Trial1-System400_Control_DroppedTagsRemoved.csv'     #<===Pond 26
}

#-------------------------------------------------------------------------------------
#read in the detection data
realDat <- read.csv(paste0(getwd(),"/Trial ", TrialNum,"/",filename), stringsAsFactors=F)
head(realDat)
realDat$X <- NULL #Beware some data sets have an extra column

#-------------------------------------------------------------------------------------
#Read in the Sound on off time data
SoundDat <- read.csv(paste0(getwd(),"/Trial ", TrialNum,"/AEH_18_CERCSOUND_ALL.csv"), stringsAsFactors=F)
SoundDat <- subset(SoundDat, Pond == PondNum & Trial == TrialNum)
#creates and adds a null off time for the first observation in the sound data
sTemp <- SoundDat[1,]
sTemp$DT <- "6/07/2018 00:00:00"
sTemp$DOY <- "6/07/2018"
sTemp$LocalTime <- "00:00.0"
sTemp$Sound <- "OFF"
SoundDat <- rbind(sTemp, SoundDat)

#create a datetime variable in as POSIXct format
SoundDat$locTimes <- as.POSIXct(SoundDat$DT,format="%m/%d/%Y %H:%M:%S")
SoundDat <- SoundDat[order(SoundDat$locTimes),]
SoundDat$Sindex <- 1:nrow(SoundDat)

#head(SoundDat, n=5)
#dim(SoundDat)

#-------------------------------------------------------------------------------------
#Read in the Pond and attribute location data
LocsDat <- read.csv(paste0(getwd(),"/CERC2018Positions.csv")) #AEH_18_CERCSOUND_01_GPS
LocsDat <- subset(LocsDat, Pond==PondNum)
names(LocsDat)[c(1, 2)] <- c("x", "y")

SPK <- subset(LocsDat, Type=='SPK')
KET <- subset(LocsDat, Type=='KET')
BND <- subset(LocsDat, Type=='BND')
HYD <- subset(LocsDat, Type=='HYD')

TagCodes <- unique(realDat$TagCode)
TagCodes

FishDT <- as.data.frame(matrix(NA, nrow=length(TagCodes), ncol=2), stringsAsFactors=F)
colnames(FishDT) <- c("minDT", "maxDT")
FishN <- as.data.frame(as.matrix(NA, nrow=length(TagCodes), ncol=1), stringsAsFactors=F)
colnames(FishN) <- c("Nobs")

#-------------------------------------------------------------------------------------
# ### Now process the data and format date times
#-------------------------------------------------------------------------------------

for (i in 1:length(TagCodes)) {
    
    TempDat <- subset(realDat, TagCode==TagCodes[i])
    TempDat <- subset(TempDat, DOY!=0)
    
    TempDat$DT <- as.POSIXct(TempDat$UtcTime, format="%m/%d/%Y %H:%M:%S", digits=3)
    TempDat <- TempDat[order(TempDat$DT),]
    
    TempDat$TimeDiff = NA
    TempDat$TimeDiff[2:nrow(TempDat)] <- TempDat$DT[2:nrow(TempDat)] - TempDat$DT[1:(nrow(TempDat)-1)]
    
    minDT = TempDat$DT[1]
    maxDT = TempDat$DT[nrow(TempDat)]
    DateDat <- as.data.frame(seq.POSIXt(from=minDT, to=maxDT, by="2 sec"), stringsAsFactors=F)
    colnames(DateDat) <- "DT"
    
    TempDat <- subset(TempDat, TimeDiff>0)
  
    TempDat <- TempDat[TempDat$DT%in%DateDat$DT,] 
    rownames(TempDat) <- 1:nrow(TempDat)
     
    # TempDat$Steps <- NA
    # TempDat$Steps[2:nrow(TempDat)] <- sqrt(((TempDat$Easting[1:nrow(TempDat)-1] - TempDat$Easting[2:nrow(TempDat)])^2) + ((TempDat$Northing[1:nrow(TempDat)-1] - TempDat$Northing[2:nrow(TempDat)])^2)) 
     
    TempDat <- merge(TempDat, DateDat, by="DT", all=T)
    TempDat <- subset(TempDat, DT >= minDT & DT <= maxDT)
    
TempDat$ID <- TagCodes[i]
TempDat$TagCode <- NULL

FishDT[i,1] <- minDT #minimum date time for the fish
FishDT[i,2] <- maxDT #maximum date time for the fish
FishN[i,1] <- nrow(TempDat)

  if (i==1) {
    AllData <- TempDat
  }else{
    AllData <- rbind(AllData, TempDat)
  }
}

FishDT[,1] <- as.POSIXct(FishDT[,1], origin="1970-01-01")
FishDT[,2] <- as.POSIXct(FishDT[,2], origin="1970-01-01")

TagCodes 
unique(AllData$ID)
head(AllData, n=10)
tail(AllData, n=10)

#-------------------------------------------------------------------------------------
# Now merge all fish positions in the pond by the second to the sound on off times
#-------------------------------------------------------------------------------------
minDT0 <- min(FishDT[,1])
maxDT0 <- max(FishDT[,2])

DateDatall <- as.data.frame(seq.POSIXt(from=minDT0, to=maxDT0, by="1 sec"), stringsAsFactors=F)
colnames(DateDatall) <- "DT"

DateDatall$hour <- as.integer(hour(DateDatall$DT))
DateDatall$Sindex <- findInterval(DateDatall$DT, SoundDat$locTimes)

TempSound <- SoundDat[,c(8,6,5)]
colnames(TempSound) <- c('Sindex','Sound','Soundtime')

DateDatall <- merge(DateDatall, TempSound, by = "Sindex")

DateDatall$Sound[DateDatall$Sound=='OFF'] <- 0
DateDatall$Sound[DateDatall$Sound=='ON'] <- 1

DateDatall$Diel <- NA

#-------------------------------------------------------------------------------------
#Obtain sunrise times
#Here's the function - sunrise.set(lat, long, date, timezone = "UTC", num.days = 1)
#using the location of the Control Pond 
SunRS <- sunrise.set(38.9122061924, -92.2795993947, "2018/06/07",  num.days = 8, timezone='America/North_Dakota/Center')
SunRS[,1] <- as.POSIXct(SunRS[,1], origin="1970-01-01", timezone='PDT')
SunRS[,2] <- as.POSIXct(SunRS[,2], origin="1970-01-01", timezone='PDT')
SunRS$Date <- as.Date(SunRS[,1])

for(j in 1:(nrow(SunRS)-1)){
  DateDatall$Diel[(DateDatall$DT >= SunRS$sunrise[j]) & (DateDatall$DT < SunRS$sunset[j])] <- 0
}
DateDatall$Diel[is.na(DateDatall$Diel)] <- 1

DateDatall$CumSnd <- cumsum(DateDatall$Sound)

AllData <- merge(AllData, DateDatall, by = "DT", sort=T)

head(AllData)
table(AllData$Diel)

#-------------------------------------------------------------------------------------
############### Export the processed data from above ###################
#-------------------------------------------------------------------------------------
# save the output data to the directory
save(AllData, file = paste0('ProcessData_',TrialNum,'_pond_',PondNum,'.RDATA'))


#-------------------------------------------------------------------------------------
#Load the processed data
load(paste0('ProcessData_',TrialNum,'_pond_',PondNum,'.RDATA'))

#Now subset dataset to just time before and after specified time interval
AnalDat <- subset(AllData, DT >= as.POSIXct("2018-06-11 06:00:00") & DT <= as.POSIXct("2018-06-11 18:00:00"))
AnalDat <- AnalDat[order(AnalDat$ID, AnalDat$DT),]
rownames(AnalDat)<- 1:nrow(AnalDat)

#Just checking the amount of time for the subset
as.POSIXct("2018-06-11 12:00:00")-as.POSIXct("2018-06-11 00:00:00")
#-------------------------------------------------------------------------------------
#just some various subsets of the data I use for graphing fish movements at sound off to on times
#Sgrab <- 2 #even numbers result in an off to on viewing of the data where 2 is the 1st sound exposure

#Tspan1 = 600 #These are seconds BEFORE the interval specified by Sgrab
#Tspan2 = 600 #These are seconds AFTER the interval specified by Sgrab
#AnalDat <- AllData
#AnalDat <- subset(AllData, Sindex==Sgrab-1 | Sindex==Sgrab)
#AnalDat <- subset(AnalDat, (as.numeric(AnalDat$DT) >= as.numeric(SoundDat$locTimes[Sgrab] - Tspan1) & as.numeric(AnalDat$DT) <= as.numeric(SoundDat$locTimes[Sgrab] + Tspan2)))
#AnalDat <- subset(AnalDat, !is.na(Sound))

#-------------------------------------------------------------------------------------
#WARNING: THIS SECTION OF CODE REQUIRES SOME PERSONAL ATTENTION
#realign code by removing leading missing values
#This will allow each fish to start with an observed position - a requirement of the model
#Adjustment will be required for each new trial, pond, and start datetime speciefied above
#Data sets for each fish
Dat1 <- subset(AnalDat, ID==TagCodes[1])#
Dat2 <- subset(AnalDat, ID==TagCodes[2])
Dat3 <- subset(AnalDat, ID==TagCodes[3])
Dat4 <- subset(AnalDat, ID==TagCodes[4])
Dat5 <- subset(AnalDat, ID==TagCodes[5])
Dat6 <- subset(AnalDat, ID==TagCodes[6])
Dat7 <- subset(AnalDat, ID==TagCodes[7])
Dat8 <- subset(AnalDat, ID==TagCodes[8])
Dat9 <- subset(AnalDat, ID==TagCodes[9])
Dat10 <- subset(AnalDat, ID==TagCodes[10])

#check the data
   cbind(Dat1[1:10,5],Dat2[1:10,5],Dat3[1:10,5],Dat4[1:10,5],Dat5[1:10,5],
        Dat6[1:10,5],Dat7[1:10,5],
        Dat8[1:10,5],Dat9[1:10,5],Dat10[1:10,5])
   
#Now readjust the data so each fish has a valid first observation
   Dat1 <- subset(AnalDat, ID==TagCodes[1])
   Dat2 <- subset(AnalDat, ID==TagCodes[2])
   Dat3 <- subset(AnalDat, ID==TagCodes[3])[-c(1:2),]
   Dat4 <- subset(AnalDat, ID==TagCodes[4])[-c(1:6),]
   Dat5 <- subset(AnalDat, ID==TagCodes[5])
   Dat6 <- subset(AnalDat, ID==TagCodes[6])
   Dat7 <- subset(AnalDat, ID==TagCodes[7])[-c(1),]
   Dat8 <- subset(AnalDat, ID==TagCodes[8])
   Dat9 <- subset(AnalDat, ID==TagCodes[9])[-c(1:4),]
   Dat10 <- subset(AnalDat, ID==TagCodes[10])[-c(1:2),]
   
#check the data again
   cbind(Dat1[1:10,5],Dat2[1:10,5],Dat3[1:10,5],Dat4[1:10,5],Dat5[1:10,5],
         Dat6[1:10,5],Dat7[1:10,5],
         Dat8[1:10,5],Dat9[1:10,5],Dat10[1:10,5])

#create a library of start times for each fish
StartDTs <- rbind(Dat1$DT[1],Dat2$DT[1],Dat3$DT[1],Dat4$DT[1],Dat5$DT[1],
                  Dat6$DT[1],Dat7$DT[1],
                  Dat8$DT[1],Dat9$DT[1],Dat10$DT[1])
StartDTs <- as.POSIXct(StartDTs, origin="1970-01-01")

#create a library of start postions for each fish
StartXs <- rbind(Dat1[1,5],Dat2[1,5],Dat3[1,5],Dat4[1,5],Dat5[1,5],
                  Dat6[1,5],Dat7[1,5],
                  Dat8[1,5],Dat9[1,5],Dat10[1,5])
StartYs <- rbind(Dat1[1,6],Dat2[1,6],Dat3[1,6],Dat4[1,6],Dat5[1,6],
                 Dat6[1,6],Dat7[1,6],
                 Dat8[1,6],Dat9[1,6],Dat10[1,6])

StartXYs <- cbind(StartXs,StartYs)

#-------------------------------------------------------------------------------------
#!!!!!! This is a Key Step
#Subset data to include only certian fish to make estimation of parameters possible
#Append the desired data sets

#Two fish at a time seems best for most computers
AnalDat <- rbind(Dat1, Dat2)
#AnalDat <- rbind(Dat3, Dat4)
#AnalDat <- rbind(Dat5, Dat6)
#AnalDat <- rbind(Dat7, Dat8)
#AnalDat <- rbind(Dat9, Dat10)

head(AnalDat, n=10)
table(AnalDat$Diel)

#-------------------------------------------------------------------------------------
#create index for the observations of each fish
nFish = length(unique(AnalDat$ID))
FishInd <- rbind(0,as.matrix(cumsum(table(AnalDat$ID))))
FishInd <- cbind(FishInd[1:nFish,1], FishInd[2:(nFish+1),1])
FishInd[,1] <- FishInd[,1] + 1
rownames(FishInd) <- 1:nrow(FishInd)


#-------------------------------------------------------------------------------------
POS <- as.data.frame(ddply(AnalDat, .(ID), summarize, 
                               DT = as.POSIXct(DT, origin="1970-01-01"),
                               Diel = Diel,
                               Easting = round((Easting-mean(Easting, na.rm=T))/sd(Easting, na.rm=T),3),
                               Northing = round((Northing-mean(Northing, na.rm=T))/sd(Northing, na.rm=T),3)), stringsAsFactors = F)
POS <- POS[order(POS$ID, POS$DT),]
POS <- POS[,-c(1:3)]

#POS <- as.matrix(round(scale(AnalDat[,5:6]),3))
rownames(POS) <- 1:nrow(POS)

#-------------------------------------------------------------------------------------
POSstats <- as.matrix(ddply(AnalDat, .(ID), summarize, 
                            MeanX = round(mean(Easting, na.rm=T),3),
                            MeanY = round(mean(Northing, na.rm=T),3),
                            SDX = round(sd(Easting, na.rm=T),3),
                            SDYY = round(sd(Northing, na.rm=T),3)))

#-------------------------------------------------------------------------------------
#initial values for locations (not needed)
initPOS <- matrix(NA, ncol=ncol(POS), nrow=nrow(POS))
initPOS[is.na(POS[,1]),1] <- 0
initPOS[is.na(POS[,2]),2] <- 0

#initial values for position differences
dif <- POS

dif[is.na(POS[,1]),1] <- 0
dif[is.na(POS[,2]),2] <- 0

diffs <-  matrix(NA, ncol=ncol(POS), nrow=nrow(POS))
diffs[2:nrow(dif),1] <- round(dif[2:nrow(dif),1],3) - round(dif[1:nrow(dif)-1,1],3)
diffs[2:nrow(dif),2] <- round(dif[2:nrow(dif),2],3) - round(dif[1:nrow(dif)-1,2],3)

Diffs <-  matrix(NA, ncol=ncol(POS), nrow=nrow(POS))
Diffs[2:nrow(POS),1] <- round(POS[2:nrow(POS),1],3) - round(POS[1:nrow(POS)-1,1],3)
Diffs[2:nrow(POS),2] <- round(POS[2:nrow(POS),2],3) - round(POS[1:nrow(POS)-1,2],3)

diffs[!is.na(Diffs[,1]),1]<-NA
diffs[!is.na(Diffs[,2]),2]<-NA

for(i in 2:nrow(diffs)){
  if (AnalDat$ID[i] != AnalDat$ID[i-1]) diffs[i,] <- NA
  if (AnalDat$ID[i] != AnalDat$ID[i-1]) Diffs[i,] <- NA
}

#-------------------------------------------------------------------------------------
#check the columns
as.matrix(colnames(AnalDat))
###
x <- as.matrix(cbind(as.numeric(AnalDat[,30]), AnalDat[,32]))
X <- matrix(NA, nrow=nrow(x), ncol=1)
#assign treatment condtions (interaction between diel period and Sound)
for (i in 1:nrow(X)){
   # if(x[i,1] == 0 & x[i,2] == 0) X[i] <- 1 # Day no sound
   # if(x[i,1] == 0 & x[i,2] == 1) X[i] <- 2 # Night no Sound
   # if(x[i,1] == 1 & x[i,2] == 0) X[i] <- 3 # Day and sound
   # if(x[i,1] == 1 & x[i,2] == 1) X[i] <- 4 # Night and Sound
  if(x[i,1] == 0) X[i] = 1 
  if(x[i,1] == 1) X[i] = 2 # no sound
  
}
table(x[,1],x[,2])
table(X)  #<---This is the covariate data for OnOff and DayNight
#-------------------------------------------------------------------------------------
#Side note: check out this movement model
#https://sourceforge.net/p/mcmc-jags/discussion/610037/thread/0f030e9a/
#-------------------------------------------------------------------------------------
#Compile the data
data.test<-list(
  nfish = nFish,
  nfact = max(X),
  N = as.matrix(FishInd),
  POS = as.matrix(POS),
  POSdiff = as.matrix(Diffs),
  X = as.vector(X),
  mu = as.vector(c(0, 0))
 # stPOS = as.matrix(StartXYs)
 # muPOS = as.vector(colMeans(AnalDat[,5:6],na.rm=T)),
 # sdPOS = as.vector(sd(AnalDat[,5:6],na.rm=T))
)
#-------------------------------------------------------------------------------------
#clean house to save memory
rm(AllData)
rm(realDat)
#-------------------------------------------------------------------------------------
#Initial values to the model
inits1 <- function() { 
  list(POSdiff=diffs)
}

#-------------------------------------------------------------------------------------
#The JAGS model output as a text file
cat(file = "XYmodbase.txt", "
model{
    
  for(j in 1:nfish){
    for(i in (N[j,1]+1):N[j,2]){

        # Likelihood
          POS[i,] ~ dsum(POS[i-1,],POSdiff[i,])
          POSdiff[i,1:2] ~ dmnorm(mu[],prec[j,X[i],,])

        # Derived quantities        
         # Stps[i] <- pow((pow(POSdiff[i,1],2) + pow(POSdiff[i,2],2)),0.5)  # Step length
         # TA[i] <- atan(POSdiff[i,2] / POSdiff[i,1])                       # Turning Angle
         # TruPOS[i,1] <- sum(POSdiff[(N[j,1]+1):i,1]) + stPOS[j,1]         # Observed value of x position
         # TruPOS[i,2] <- sum(POSdiff[(N[j,1]+1):i,2]) + stPOS[j,2]         # Observed value of y position
    }
  }  
     
  for(j in 1:nfish){
    for (k in 1:nfact){

    # Flat priors
      rho[j,k] ~ dunif(-1,1)
      sigma[j,k,1] ~ dunif(0,5)
      sigma[j,k,2] ~ dunif(0,5)

    # Construct covariance matrix and corresponding precision
      prec[j,k,1:2,1:2] <- inverse(cov[j,k,1:2,1:2])
      cov[j,k,1,1] <- sigma[j,k,1] * sigma[j,k,1]
      cov[j,k,1,2] <- sigma[j,k,1] * sigma[j,k,2] * rho[j,k]
      cov[j,k,2,1] <- sigma[j,k,2] * sigma[j,k,1] * rho[j,k]
      cov[j,k,2,2] <- sigma[j,k,2] * sigma[j,k,2]

    }
  }
}", fill = TRUE)

#-------------------------------------------------------------------------------------
#Parameters to keep in output
parms <- c("rho","sigma","cov","POS") #"POSdiff", "TA", "Stps"

# MCMC settings
ni <- 1000
nt <- 20
nb <- 1000
nc <- 3

#-------------------------------------------------------------------------------------
#Call JAGS from R using the runjags package
rm(jags.out)
start = Sys.time()
jags.out <- run.jags(model = "XYmodbase.txt",  monitor = parms,  inits=inits1,
                     data = data.test, n.chains = nc, burnin = nb, sample = ni, adapt = 1000,
                     thin = nt, method="parallel")
elapsed = Sys.time() - start
elapsed

#-------------------------------------------------------------------------------------
#Summary of output and diagnotic plots of the posteriors
SumStats <-as.data.frame(summary(jags.out, vars = c("rho","sigma")))
SumStats$name <- rownames(SumStats)
rownames(SumStats) <- 1:nrow(SumStats)

#basic diagnostic plots
plot(jags.out, vars = c("rho","sigma"))#,"sigma"

#Save JAGS model output
save(jags.out, file = paste0('Model_Output',TrialNum,'_pond_',PondNum,'fish 1 and 2.RDATA'))


#-------------------------------------------------------------------------------------
#Calculate derived variables from the imputed and observed positions
Stps <- summary(jags.out, vars = c("Stps"))

mcmc <- combine.mcmc(jags.out)

testlist <- list.to.array(mcmcList=mcmc, parNames=c("POS"), NAfill = NA)
str(testlist)


Stps <- array(testlist$POSdiff)
TAs <- array(testlist$POSdiff)

for(Fnum in 1:nFish){
  
  Stps[,Fnum,] <- sqrt(testlist$POSdiff[,FishInd[Fnum,1]:FishInd[Fnum,2],1]^2 - testlist$POSdiff[,FishInd[Fnum,1]:FishInd[Fnum,2],2]^2)
  
  if(testlist$POSdiff[,FishInd[Fnum,1]:FishInd[Fnum,2],1] != 0){
  
    TAs[,Fnum,] <- atan(testlist$POSdiff[,FishInd[Fnum,1]:FishInd[Fnum,2],2] / testlist$POSdiff[,FishInd[Fnum,1]:FishInd[Fnum,2],1])    
  }else{
    TAs[,Fnum,] <- 0
  }

}

##############################################################
rhos <- SumStats[1:(2*nFish),c(12,1:3)]
tiffName = paste0('/Trial_', TrialNum , '_Pond_', PondNum, 'Rho by fish.tif')
tiff(paste0(getwd(), tiffName),
    width = 8, height = 6, units = "in", res= 300)
#windows(8, 6) 
    plot(1:nFish, rhos[1:nFish,3], pch=NA,
         xlim = c(0,(nFish+1)), ylim=c(-1,1),
         xlab = 'Fish number',
         ylab = 'Movement correlation')
    abline(h=0, lty=3)
    polygon(c(1:nFish, rev(1:nFish)), c(rhos[(nFish+1):(2*nFish),2], rev(rhos[(nFish+1):(2*nFish),4])),  lwd=1, lty=1, col="pink", border="pink")
    polygon(c(1:nFish, rev(1:nFish)), c(rhos[1:nFish,2], rev(rhos[1:nFish,4])),  lwd=1, lty=1, col="grey", border="grey")
    
    points(1:nFish, rhos[1:nFish,3], type='o',  col='black', pch=19)
    points(1:nFish, rhos[(nFish+1):(2*nFish),3], type='o', col='red', pch=19)
    legend('topleft', c("Before sound", "After sound"), lty = c(1,1), 
           lwd=c(2,2), col=c('black','red'), pch = c(19,19), cex=1, bty="n", y.intersp = 1)
    legend('topright', c(paste0('Trial ',TrialNum), paste0('Pond ',PondNum), paste0(Senario)), lty = c(NA,NA,NA), 
           lwd=c(NA,NA,NA), col=c('black','black','black'), pch = c(NA,NA,NA), cex=1, bty="n", y.intersp = 1)
dev.off()


Sigs <- SumStats[(2*nFish+1):nrow(SumStats),c(12,1:3)]
rownames(Sigs) <- 1:nrow(Sigs)

tiffName = paste0('/Trial_', TrialNum , '_Pond_', PondNum, 'Easting Variance by fish.tif')
tiff(paste0(getwd(), tiffName),
     width = 8, height = 6, units = "in", res= 300)

#windows(8, 6) 
plot(1:nFish, Sigs[1:nFish,3], pch=NA,
     xlim = c(0,(nFish+1)), ylim=c(0,2),
     xlab = 'Fish number',
     ylab = 'Variance in Easting')
abline(h=0, lty=3)
polygon(c(1:nFish, rev(1:nFish)), c(Sigs[1:nFish,2], rev(Sigs[1:nFish,4])),  lwd=1, lty=1, col="grey", border="grey")
polygon(c(1:nFish, rev(1:nFish)), c(Sigs[(1+nFish):(2*nFish),2], rev(Sigs[(1+nFish):(2*nFish),4])),  lwd=1, lty=1, col="pink", border="pink")
points(1:nFish, Sigs[1:nFish,3], type='o', col='black', pch=19)
points(1:nFish, Sigs[(1+nFish):(2*nFish),3], type='o', col='red', pch=19)
legend('topleft', c("Before sound", "After sound"), lty = c(1,1), 
       lwd=c(2,2), col=c('black','red'), pch = c(19,19), cex=1, bty="n", y.intersp = 1)
legend('topright', c(paste0('Trial ',TrialNum), paste0('Pond ',PondNum), paste0(Senario)), lty = c(NA,NA,NA), 
       lwd=c(NA,NA,NA), col=c('black','black','black'), pch = c(NA,NA,NA), cex=1, bty="n", y.intersp = 1)
dev.off()


tiffName = paste0('/Trial_', TrialNum , '_Pond_', PondNum, 'Norhting Variance by fish.tif')
tiff(paste0(getwd(), tiffName),
     width = 8, height = 6, units = "in", res= 300)
#windows(8, 6) 
plot(1:nFish, Sigs[(2*nFish+1):(nFish*3),3], pch=NA,
     xlim = c(0,(nFish+1)), ylim=c(0,2),
     xlab = 'Fish number',
     ylab = 'Variance in Northing')
abline(h=0, lty=3)
polygon(c(1:nFish, rev(1:nFish)), c(Sigs[(2*nFish+1):(nFish*3),2], rev(Sigs[(2*nFish+1):(nFish*3),4])),  lwd=1, lty=1, col="grey", border="grey")
polygon(c(1:nFish, rev(1:nFish)), c(Sigs[(3*nFish+1):(4*nFish),2], rev(Sigs[(3*nFish+1):(4*nFish),4])),  lwd=1, lty=1, col="pink", border="pink")
points(1:nFish, Sigs[(2*nFish+1):(nFish*3),3], type='o', col='black', pch=19)
points(1:nFish, Sigs[(3*nFish+1):(4*nFish),3], type='o', col='red', pch=19)
legend('topleft', c("Before sound", "After sound"), lty = c(1,1), 
       lwd=c(2,2), col=c('black','red'), pch = c(19,19), cex=1, bty="n", y.intersp = 1)
legend('topright', c(paste0('Trial ',TrialNum), paste0('Pond ',PondNum), paste0(Senario)), lty = c(NA,NA,NA), 
       lwd=c(NA,NA,NA), col=c('black','black','black'), pch = c(NA,NA,NA), cex=1, bty="n", y.intersp = 1)
dev.off()

#############################################################
###Simulate the results
#Fnum = 6
muPOS = as.vector(colMeans(AnalDat[,5:6],na.rm=T))
sdPOS = as.vector(sd(AnalDat[,5:6],na.rm=T))
Nsim = 60
mu <- c(0, 0)
sigmaN <- cbind(SumStats[21:30,2], SumStats[41:50,2])
rhoN <- SumStats[1:10,2]

sigmaS <- cbind(SumStats[31:40,2], SumStats[51:60,2])
rhoS <- SumStats[11:20,2]

Bnd <- matrix(NA, 2, 2)
Bnd[1,1] <- range(BND$x)[1]
Bnd[2,1] <- range(BND$x)[2]
Bnd[1,2] <- range(BND$y)[1]
Bnd[2,2] <- range(BND$y)[2]
colnames(Bnd)<- c("x","y")

for (Fnum in 1:10){
  
#  Fnum = 1

  cov_matN <- as.matrix(rbind(c(sigmaN[Fnum,1]^2, sigmaN[Fnum,1]*sigmaN[Fnum,2]*rhoN[Fnum]),
                   c( sigmaN[Fnum,1]*sigmaN[Fnum,2]*rhoN[Fnum], sigmaN[Fnum,2]^2)))
  
  cov_matS <- as.matrix(rbind(c(sigmaS[Fnum,1]^2, sigmaS[Fnum,1]*sigmaS[Fnum,2]*rhoS[Fnum]),
                    c(sigmaS[Fnum,1]*sigmaS[Fnum,2]*rhoS[Fnum],sigmaS[Fnum,2]^2)))
  
  XN <- as.matrix(rmvnorm(Nsim, mu, cov_matN))
  XS <- as.matrix(rmvnorm(Nsim, mu, cov_matS))
  Xn <- XN*0
  Xs <- XS*0
  
  Xn[1,1] <- muPOS[1] 
  Xn[1,2] <- muPOS[2] 
  Xs[1,1] <- muPOS[1] 
  Xs[1,2] <- muPOS[2] 
  
  for (i in 1:(Nsim-1)){
    Xn[i+1,1] <- Xn[i,1] + XN[i,1]
    Xn[i+1,2] <- Xn[i,2] + XN[i,2]
    
    Xs[i+1,1] <- Xs[i,1] + XS[i,1]
    Xs[i+1,2] <- Xs[i,2] + XS[i,2]
    
    if(Xn[i,1] <= Bnd[1,1]) Xn[i,1] <- Bnd[1,1]+1
    if(Xn[i,1] >= Bnd[2,1]) Xn[i,1] <- Bnd[2,1]-1

    if(Xn[i,2] <= Bnd[1,2]) Xn[i,2] <- Bnd[1,2]+1
    if(Xn[i,2] >= Bnd[2,2]) Xn[i,2] <- Bnd[2,2]-1
    
    if(Xs[i,1] <= Bnd[1,1]) Xs[i,1] <- Bnd[1,1]+1
    if(Xs[i,1] >= Bnd[2,1]) Xs[i,1] <- Bnd[2,1]-1
    
    if(Xs[i,2] <= Bnd[1,2]) Xs[i,2] <- Bnd[1,2]+1
    if(Xs[i,2] >= Bnd[2,2]) Xs[i,2] <- Bnd[2,2]-1
    
  }

  Sim <- cbind(rep(Fnum, Nsim), Xn, Xs)
  colnames(Sim) <- c("Fnum","Xoff","Yoff", "Xon","Yon")

  if (Fnum == 1){
    SimDat <- Sim
  }else{
    SimDat <- rbind(SimDat, Sim)
  }
}

################
split.screen(c(2,1))
screen(1)
hist(SimDat[,2], xlim=range(SimDat[,c(2,4)]))
screen(2)
hist(SimDat[,4], xlim=range(SimDat[,c(2,4)]))

################
Fnum=6

plot(SimDat[SimDat[,1]==Fnum,2], SimDat[SimDat[,1]==Fnum,3], pch=NA, lty=1, type='o', col='black', 
    ylim=c(Bnd[1,2], Bnd[2,2]+10), 
     xlim=c(Bnd[1,1], Bnd[2,1]+10),
     ylab='Northing', xlab='Easting', cex.lab=1.5)
     lines(BND$x, BND$y, col='brown')
     points(SimDat[SimDat[,1]==Fnum,4], SimDat[SimDat[,1]==Fnum,5], 
            pch=NA, col='red', type='o')
     points(AnalDat[FishInd[Fnum,1]:FishInd[Fnum,2],5],AnalDat[FishInd[Fnum,1]:FishInd[Fnum,2],6], 
            pch=NA, col='blue', type='o')

plot(Xn, type='o')#, xlim=c(-2, 2), ylim=c(-2, 2))
points(Xs, col='red', type='o')

##################################
split.screen(c(1,2))
screen(1)
boxplot(SimDat[,2] ~ SimDat[,1],
        xlim = c(0,11),
        ylim = c(min(BND$x),max(BND$x)+10),
        ylab = "Simulated change in easting (Off)",
        xlab = "Fish number")
abline(h=0, lty=3)

screen(2)
boxplot(SimDat[,4] ~ SimDat[,1],
        xlim = c(0,11),
        ylim = c(min(BND$x),max(BND$x)+10),
        ylab = "Simulated change in easting (On)",
        xlab = "Fish number")
abline(h=0, lty=3)

##################################
split.screen(c(1,2))
screen(1)
boxplot(SimDat[,3] ~ SimDat[,1],
        xlim = c(0,11),
        ylim = c(min(BND$y),max(BND$y)+10),
        ylab = "Simulated change in northing (Off)",
        xlab = "Fish number")
abline(h=0, lty=3)

screen(2)
boxplot(SimDat[,5] ~ SimDat[,1],
        xlim = c(0,11),
        ylim = c(min(BND$y),max(BND$y)+10),
        ylab = "Simulated change in northing (On)",
        xlab = "Fish number")
abline(h=0, lty=3)


















# ########################################################################
# #######################################################################
# 
# 
# 
# for (i in 1:length(TagCodes)) {
#   
#   GData <- subset(AllData, ID==TagCodes[i])
#   GData <- GData[,c(1,5,6)]
#   colnames(GData)<-c('DT',paste0('X_',i),paste0('Y_',i))
#   GData <- GData[order(GData$DT),]
#   
#   if (i==1) {
#     FishX = merge(DateDatall, GData,  by="DT", all=T)
#   } else {
#     FishX = merge(FishX, GData, by = 'DT', all=T)
#   } 
#   
# }
# 
# ##########################################################################
# FishX <- subset(FishX, DT >= as.POSIXct("2018-06-10 00:00:00"))
# FishX <- subset(FishX, DT <= as.POSIXct("2018-06-14 00:00:00"))
# 
# #creates a color vector
# # FishX$colindx <- with(FishX, NA)
# # FishX$colindx[FishX$Sound==0] <- 'blue'
# # FishX$colindx[FishX$Sound==1] <- 'red'
# 
# #######Check out the objects
# head(FishX, n=90)
# tail(FishX, n=30)
# dim(FishX)

########################################################################
############### Export the processed data from above ###################
########################################################################
# save the output data to the directory
save(AllData, file = paste0('Time series Trial_',TrialNum,'_pond_',PondNum,'.RDATA'))
save(FishX, file = paste0('FishPositions_Trial_',TrialNum,'_pond_',PondNum,'.RDATA'))


########################################################################
# This code is for graphing the data in .gif movie format
########################################################################
Sgrab <- 2 #even numbers result in an off to on viewing of the data where 2 is the 1st sound burs

#GraphDat<-FishX
GraphDat <- subset(FishX, Sindex==Sgrab-1 | Sindex==Sgrab)
GraphDat <- subset(GraphDat, (as.numeric(GraphDat$DT) >= as.numeric(SoundDat$locTimes[Sgrab] - 300) & as.numeric(GraphDat$DT) <= as.numeric(SoundDat$locTimes[Sgrab] + 300)))
GraphDat <- subset(GraphDat, !is.na(Sound))
nrow(GraphDat)
# unique(GraphDat$Sound)
# head(GraphDat, n=15)

#Read in the Pond and attribute location data
LocsDat <- read.csv(paste0(getwd(),"/CERC2018Positions.csv")) #AEH_18_CERCSOUND_01_GPS
LocsDat <- subset(LocsDat, Pond==PondNum)
names(LocsDat)[c(1, 2)] <- c("x", "y")

SPK <- subset(LocsDat, Type=='SPK')
KET <- subset(LocsDat, Type=='KET')
BND <- subset(LocsDat, Type=='BND')
HYD <- subset(LocsDat, Type=='HYD')

saveGIF(convert = "convert", clean=T, interval=0.1, movie.name = paste0("Trial ",TrialNum," Pond ",PondNum,".gif"), ani.width = 1500, ani.height = 800,
        expr = {   
          
          for(i in 1:nrow(GraphDat)){
            #windows(9, 11) 
            plot(GraphDat[1,2], GraphDat[1,3], pch="", col='white', 
                 ylim=c(min(LocsDat$y),max(LocsDat$y)+10), 
                 xlim=c(min(LocsDat$x),max(LocsDat$x)+10),
                 ylab='Northing', xlab='Easting', cex.lab=1.5)
            
            #Speaker locations
            points(SPK$x,SPK$y, col='black', cex=3.5, pch=12) 
            
            #Pond boundary locations
            lines(BND$x, BND$y, type='o', col='brown', lwd=3, pch="")
            
            #Kettle locations
            lines(KET$x, KET$y, type='o', col='magenta', lwd=3, pch="")
            
            legend('topright', c("Sound off", "Sound on", "Speaker location", "Kettle location","Pond boundary"), lty = c(1,1,NA,1,1), 
                   lwd=c(3,3,NA,3,3), col=c('blue','red','black','magenta','brown'), pch = c(NA,NA,12,NA,NA), cex=1, bty="n", y.intersp = 1)
            
            legend('topleft', c(paste0('Trial ',TrialNum), paste0('Pond ',PondNum), paste0('Fish 6-10')), lty = c(NA,NA,NA), 
                   lwd=c(NA,NA,NA), col=c('black','black','black'), pch = c(NA,NA,NA), cex=1, bty="n", y.intersp = 1)
            
            # #fish 1
            segments(GraphDat[1:i,2], GraphDat[1:i,3], GraphDat[1:i+1,2], GraphDat[1:i+1,3], col=GraphDat$colindx[1:i], lty=1, lwd=0.5)
            # #fish 2
            segments(GraphDat[1:i,4], GraphDat[1:i,5], GraphDat[1:i+1,4], GraphDat[1:i+1,5], col=GraphDat$colindx[1:i], lty=1, lwd=0.75)
            # #fish 3
            segments(GraphDat[1:i,6], GraphDat[1:i,7], GraphDat[1:i+1,6], GraphDat[1:i+1,7], col=GraphDat$colindx[1:i], lty=1, lwd=1)
            # #fish 4
            segments(GraphDat[1:i,8], GraphDat[1:i,9], GraphDat[1:i+1,8], GraphDat[1:i+1,9], col=GraphDat$colindx[1:i], lty=1, lwd=1.25)
            # #fish 5
            segments(GraphDat[1:i,10], GraphDat[1:i,11], GraphDat[1:i+1,10], GraphDat[1:i+1,11], col=GraphDat$colindx[1:i], lty=1, lwd=1.5)
            # #fish 6
            segments(GraphDat[1:i,12], GraphDat[1:i,13], GraphDat[1:i+1,12], GraphDat[1:i+1,13], col=GraphDat$colindx[1:i], lty=1, lwd=1.75)
            # #fish 7
            segments(GraphDat[1:i,14], GraphDat[1:i,15], GraphDat[1:i+1,14], GraphDat[1:i+1,15], col=GraphDat$colindx[1:i], lty=1, lwd=2)
            # #fish 8
            segments(GraphDat[1:i,16], GraphDat[1:i,17], GraphDat[1:i+1,16], GraphDat[1:i+1,17], col=GraphDat$colindx[1:i], lty=1, lwd=2.25)
            # #fish 9
            segments(GraphDat[1:i,18], GraphDat[1:i,19], GraphDat[1:i+1,18], GraphDat[1:i+1,19], col=GraphDat$colindx[1:i], lty=1, lwd=2.5)
            # #fish 10
            segments(GraphDat[1:i,20], GraphDat[1:i,21], GraphDat[1:i+1,20], GraphDat[1:i+1,21], col=GraphDat$colindx[1:i], lty=1, lwd=2.75)
          }
        })
dev.off()
























########################################################################
############## Read in the processed data from above ###################
########################################################################
load(paste0(getwd(), '/Time series Trial_',TrialNum,'_pond_',PondNum,'.RDATA'))
   TmpDat <- subset(AllData, ID==TagCodes[1])
   windows(4, 6) 
    plot(TmpDat$time[-1],TmpDat$step[-1])
load(paste0(getwd(), '/FishPositions_Trial_',TrialNum,'_pond_',PondNum,'.RDATA'))
load(paste0(getwd(), '/FishAngles_Trial_',TrialNum,'_pond_',PondNum,'.RDATA'))
load(paste0(getwd(), '/FishSteps_Trial_',TrialNum,'_pond_',PondNum,'.RDATA'))

range(FishS$DT)
range(FishA$DT)
range(FishX$DT)
dim(FishX)
dim(FishA)
dim(FishS)


#############################################################################
#Plots of X and Y over time
tiffName = paste0('/Trial_', TrialNum , '_Pond_', PondNum, 'Easting by time.tif')
tiff(paste0(getwd(), tiffName),
     width = 9, height = 11, units = "in", res= 300)
windows(9, 11) 
par(mfrow = c(5,2), mar = c(4,4,2,2), oma = c(1.5, 2,0,0))

for (i in seq(2,20, by=2)){
 # i=2
plot(FishX$DT, FishX[,2], 
     xlim = c(min(FishX$DT), max(FishX$DT+(86400*2))),
     ylim = c(min(LocsDat$x), max(LocsDat$x)+1),
     xlab = 'Date',
     ylab = 'Easting (UTM)',
     type = 'l',
     lty = 1,
     cex = 1,
     col = 'white')
  segments(FishX$DT,min(LocsDat$x)*as.numeric(FishX$Diel),FishX$DT, max(LocsDat$x)*as.numeric(FishX$Diel), col='grey')
  segments(FishX$DT,min(LocsDat$x)*as.numeric(FishX$Sound),FishX$DT, max(LocsDat$x)*as.numeric(FishX$Sound), col='red')
  lines(FishX$DT, FishX[,i], lty=1, col='black')
  lines(FishX$DT, rep(max(LocsDat$x),nrow(FishX)), lty=3)
  lines(FishX$DT, rep(min(LocsDat$x),nrow(FishX)), lty=3)
  
  legend('topright', c('Sound On', 'Night', 'Pond boundaries', paste0('Tag = ', substr(colnames(FishX)[i],3,9)), paste0('Trial ', TrialNum) , paste0('Pond ',PondNum)), 
         lty = c(1,1,3,1,NA,NA), 
         lwd = c(3,3,1,1,NA,NA), 
         col=c('red','grey','black','black',NA,NA), 
         pch = c(NA,NA,NA,NA,NA), 
         cex=0.8, 
         bty="n",
        # box.col='white',
        # bg = 'white',
         y.intersp = 0.8)
}
dev.off()

#################################

tiffName = paste0('/Trial_', TrialNum , '_Pond_', PondNum, 'Northing by time.tif')
tiff(paste0(getwd(), tiffName),
     width= 9, height= 11, units= "in", res= 300)

#windows(9, 11) 
par(mfrow = c(5,2), mar = c(4,4,2,2), oma = c(1.5, 2,0,0))

for (i in seq(3,21, by=2)){
  
  plot(FishX$DT, FishX[,3], 
       xlim = c(min(FishX$DT), max(FishX$DT+(86400*2))),
       ylim =  c(min(LocsDat$y), max(LocsDat$y)+1),
       xlab = 'Date',
       ylab = 'Northing (UTM)',
       type = 'l',
       lty = 1,
       cex = 1,
       col = 'white')
  segments(FishX$DT,min(LocsDat$y)*as.numeric(FishX$Diel),FishX$DT, max(LocsDat$y)*as.numeric(FishX$Diel), col='grey')
  segments(FishX$DT,min(LocsDat$y)*as.numeric(FishX$Sound),FishX$DT, max(LocsDat$y)*as.numeric(FishX$Sound), col='red')
  lines(FishX$DT, FishX[,i], lty=1, col='black')
  lines(FishX$DT, rep(max(LocsDat$y),nrow(FishX)), lty=3)
  lines(FishX$DT, rep(min(LocsDat$y),nrow(FishX)), lty=3)
  
  legend('topright', c('Sound On', 'Night', 'Pond boundaries', paste0('Tag = ', substr(colnames(FishX)[i],3,9)), paste0('Trial ', TrialNum) , paste0('Pond ',PondNum)), 
         lty = c(1,1,3,1,NA,NA), 
         lwd = c(3,3,1,1,NA,NA), 
         col=c('red','grey','black','black',NA,NA), 
         pch = c(NA,NA,NA,NA,NA), 
         cex=0.8, 
         bty="n",
         # box.col='white',
         # bg = 'white',
         y.intersp = 0.8)
}
dev.off()

#############################################################################
#Plots of step length for each fish
tiffName = paste0('/Trial_', TrialNum , '_Pond_', PondNum, 'Step length by time.tif')
tiff(paste0(getwd(), tiffName),
     width = 9, height = 11, units = "in", res= 300)

#windows(9, 11) 
par(mfrow = c(5,2), mar = c(4,4,2,2), oma = c(1.5, 2,0,0))

for (i in seq(2,11)){
  
  plot(FishS$DT, FishS[,2], 
       xlim = c(min(FishS$DT), max(FishS$DT+(86400*2.5))),
       ylim = c(0, 15),
       xlab = 'Date',
       ylab = 'Step length (m/s)',
       type = 'l',
       lty = 1,
       cex = 1,
       col = 'white')
  segments(FishS$DT,0*as.numeric(FishS$Diel),FishS$DT, 15*as.numeric(FishS$Diel), col='grey')
  segments(FishS$DT,0*as.numeric(FishS$Sound),FishS$DT, 15*as.numeric(FishS$Sound), col='red')
  lines(FishS$DT, FishS[,i], lty=1, col='black')
 
  
  legend('topright', c('Sound On', 'Night',  paste0('Tag = ', substr(colnames(FishS)[i],3,9)), paste0('Trial ', TrialNum) , paste0('Pond ',PondNum)), 
         lty = c(1,1,1,NA,NA), 
         lwd = c(3,3,1,NA,NA), 
         col=c('red','grey','black',NA,NA), 
         pch = c(NA,NA,NA,NA), 
         cex=1, 
         bty="n",
         # box.col='white',
         # bg = 'white',
         y.intersp = 0.8)
}
dev.off()

#####################################################################
#now angle
# tiffName= paste0("/Truning angles by time.tif")
# tiff(paste0(getwd(), tiffName),
#      width= 9, height= 11, units= "in", res= 300)
# 
# windows(9, 11) 
# par(mfrow = c(5,2), mar = c(4,4,2,2), oma = c(1.5, 2,0,0))
# 
# for (i in seq(2,11)){
#  
#   plot(FishA$DT, FishA[,2], 
#        xlim = c(min(FishA$DT), max(FishA$DT+(86400*3.1))),
#        ylim = c(-3.142, 3.142),
#        xlab = 'Date',
#        ylab = 'Turning angle (radians)',
#        type = 'l',
#        lty = 1,
#        cex = 1,
#        col = 'white')
#   segments(FishA$DT,-3.142*as.numeric(FishA$Diel),FishA$DT, 3.142*as.numeric(FishA$Diel), col='grey')
#   segments(FishA$DT,-3.142*as.numeric(FishA$Sound),FishS$DT, 3.142*as.numeric(FishA$Sound), col='red')
#   lines(FishA$DT, FishA[,i], lty=1, col='black')
#   
#   
#   legend('topright', c('Sound On', 'Night',  paste0('Tag = ', substr(colnames(FishA)[i],3,9)), paste0('Trial ', TrialNum) , paste0('Pond ',PondNum)), 
#          lty = c(1,1,1,NA,NA), 
#          lwd = c(3,3,1,NA,NA), 
#          col=c('red','grey','black',NA,NA), 
#          pch = c(NA,NA,NA,NA), 
#          cex=0.8, 
#          bty="n",
#          # box.col='white',
#          # bg = 'white',
#          y.intersp = 0.8)
# }
# dev.off()

########################################################################
#plots of angles and circular stats
library(circular)

tiffName= paste0('/Trial_', TrialNum , '_Pond_', PondNum, 'Angles by fish.tif')
tiff(paste0(getwd(), tiffName),
     width= 9, height= 11, units= "in", res= 300)

#windows(9, 11) 
par(mfrow = c(5,2), mar = c(1,0,1,1), oma = c(1, 1, 0, 0))

for (i in seq(2,11)){

  Cdir<-circular(FishA[,i])
  
  plot(Cdir, shrink=2, stack=TRUE, pch=1, bins=1, cex=1, col=NA)

  ests <- mle.vonmises(Cdir)

  legend('topleft', c('Density', 
                      'Frequency',  
                      paste0('Tag = ', substr(colnames(FishA)[i],3,9)), 
                      paste0('Trial ', TrialNum) , 
                      paste0('Pond ', PondNum), 
                      paste0('Mean = ', round(ests$mu[[1]],4), ' (',round(ests$se.mu[[1]],4),')'),
                      paste0('kappa = ', round(ests$kappa[[1]],4), ' (',round(ests$se.kappa[[1]],4),')')), 
         lty = c(1,1,NA,NA,NA,NA,NA), 
         lwd = c(3,3,NA,NA,NA,NA,NA), 
         col=c('red','lightblue','black',NA,NA,'black','black'), 
         pch = c(NA,NA,NA,NA,NA,NA), 
         cex=0.9, 
         bty="n",
         # box.col='white',
         # bg = 'white',
         y.intersp = 0.8)
  
  lines(density.circular(Cdir, bw=100), lwd=2, col='red')
  rose.diag(Cdir, bins=36, col='lightblue', units='radians', 
            cex=1, prop=1.5, add=TRUE)

}
dev.off()
##########################################################

Cdir<-circular(FishA[,2])
covs1<-as.matrix(cbind(as.numeric(FishA$Sound),  as.numeric(FishA$Diel)))
Circmod1<-lm.circular(type="c-l", y=Cdir, x=covs1, init=c(0,0))

estimates <- mle.vonmises(Cdir)
estimates$mu[[1]]
estimates$se.mu[[1]]
estimates$kappa[[1]]
estimates$se.kappa[[1]]


####################################################
#Using the momentHMM package and MLE
#fit an HMM with covariate for ON vs Off...
#first let's prep the data...
head(AllData)
tail(AllData)

chckit <- subset(AllData, is.na(step))

chck <- hist(AllData$step)
range(AllData$step, na.rm=T)

hist(AllData$angle, xlim=c(0,10))

##autocorrelation among step length and turning angle
acf(FishS[,2], lag.max = 50)
acf(FishA[,2], lag.max = 10)

AllData <- AllData[-nrow(AllData),]

##############################################################################
#Bayesian model fitting

n_kappas <- 15
kappas <- 1:n_kappas
I0_s   <- I.0(kappas) # from package 'circular' 


#######################################################################
forDat <- merge(FishS,FishA[,1:11], by='DT')
forDat$check <- 1:nrow(forDat)/60

forDat <- forDat[forDat$check%in%seq(1, round(nrow(forDat)/60,0), by=1),]

head(forDat)
dim(forDat)

FishX$check <- 1:nrow(FishX)/60
FishX <- FishX[FishX$check%in%seq(1, round(nrow(FishX)/60,0), by=1),]
#
 FishS$check <- 1:nrow(FishS)/60
 FishS <- FishS[FishS$check%in%seq(1, round(nrow(FishS)/60,0), by=1),]
# 
 FishA$check <- 1:nrow(FishA)/60
 FishA <- FishA[FishA$check%in%seq(1, round(nrow(FishA)/60,0), by=1),]
 
 hist(FishS[,2])


#######################################################################
#for just one fish
jags.dat = list(S = as.vector(forDat[,2]),
                A = as.vector(forDat[,19]),
                N = nrow(forDat),
                I0s = as.vector(I0_s),
                dummy = as.vector(rep(0, nrow(forDat))),
                kappas = as.vector(kappas),
                n_kappas = n_kappas,
                Sound = as.vector(as.numeric(forDat$Sound)),
                Diel = as.vector(as.numeric(forDat$Diel)))#,
                #CumS = as.vector(FishS$CumS))
jags.dat 

# inits1 <- function() { 
#   list(n_arr=arr_inits)
# }

cat(file = "HMMmodel.txt", "
    model {
    # priors and constraints

    #b1 ~ dnorm(0, .0001) I(0, ) #shape parameter for weibull
    b2 ~ dnorm(0, .0001) I(0, ) #scale parameter for weibull

    beta0 ~ dnorm(0, 0.01) # intercept
	  beta1 ~ dnorm(0, 0.01) # slope
    beta2 ~ dnorm(0, 0.01) # slope
    #beta3 ~ dnorm(0, 0.01) # slope
    
    mu  ~ dunif(0,pi2)
    C   <- 1000000   # for the zero's trick
    pi  <- 3.14159
    pi2 <- 2*pi

 # Priors for Autocorrelation
 #rho ~ dunif(-1,1)
 #gamma ~ dnorm(0, 1/100) I(0, )

### Model for step lengths

  # Autocorrelation on step length
    # S[1] ~ dgamma(shape[1], scale[1])
    # shape[1] <- exp(beta0 + beta1*Sound[1] + cor_err[1]) #+ beta2*Diel[1]) #) 
    # 
    # scale[1] <- b2
    # 
    # cor_err[1] ~ dnorm(0, gamma * (1-rho*rho))

  # angle model
    # phi[1] <- - kappa_hat * cos(A[1]-mu) + log(2*pi*IO_hat) + C
    # dummy[1] ~ dpois(phi[1])

 for (i in 1:N){
    
      # pr_err[i] ~ dnorm(0, gamma)
      # cor_err[i] <- pr_err[i] + rho * cor_err[i-1]
      

### Model for step lengths
      S[i] ~ dgamma(shape[i], scale[i])
      shape[i] <- exp(beta0 + beta1*Sound[i]+ beta2*Diel[i]) #+ cor_err[i] # + beta3*CumS 
    
      #shape[i] <- b1
      scale[i] <- b2

### Model for angles
      phi[i] <- - kappa_hat * cos(A[i]-mu) + log(2*pi*IO_hat) + C
      #dummy[i] <- 0
      dummy[i] ~ dpois(phi[i])
    }
    
    k  ~ dcat(p[])
    
    for(j in 1:n_kappas) {
      p[j] <- 1/n_kappas  # creating vector p = {1/n_kappas,...}
    }
    
    kappa_hat <- kappas[k]
    IO_hat    <- I0s[k]

}", fill = TRUE)

#------ Run the JAGS model ---------
parms <- c("beta0","beta1","beta2","shape","scale","mu") #"rho""mean_n",,"rho",

# MCMC settings
ni <- 1000
nt <- 20
nb <- 1000
nc <- 3
#rm(jags.out)
#Call JAGS from R using the runjags package
start = Sys.time()
jags.out <- run.jags(model = "HMMmodel.txt",  monitor = parms,  #inits=inits1,
                     data = jags.dat, n.chains = nc, burnin = nb, sample = ni, adapt = 1000,
                     thin = nt, method="parallel")
elapsed = Sys.time() - start
elapsed

SumStats <-summary(jags.out, vars = c("beta0","beta1","beta2","mu"))#"beta2","shape","scale","rho",
SumStats

plot(jags.out, vars = c("beta0","beta1","beta2","mu")) #"beta2""rho",




#############################################################3
#now for adding random effect of multiple fish
jags.dat2 = list(S = as.matrix(forDat[,2:11]),
                A = as.matrix(forDat[,19:28]),
                N = nrow(forDat),
                I0s = as.vector(I0_s),
                dummy = matrix(0, nrow=nrow(forDat), ncol=10),
                kappas = as.vector(kappas),
                n_kappas = n_kappas,
                Sound = as.vector(as.numeric(forDat$Sound)),
                Diel = as.vector(as.numeric(forDat$Diel)))#,
#CumS = as.vector(FishS$CumS))
jags.dat2


cat(file = "HMMmodelallfish.txt", "
    model {
    # priors and constraints
    
    #b1 ~ dnorm(0, .0001) I(0, ) #shape parameter for weibull
    b2 ~ dnorm(0, .0001) I(0, ) #scale parameter for weibull
    
    beta0 ~ dnorm(0, 0.01) # intercept
    beta1 ~ dnorm(0, 0.01) # slope
    beta2 ~ dnorm(0, 0.01) # slope
    #beta3 ~ dnorm(0, 0.01) # slope
    
    mu  ~ dunif(0,pi2)
    C   <- 1000000   # for the zero's trick
    pi  <- 3.14159
    pi2 <- 2*pi
    
    ## Random effect of individual
    prec.eps <- 1/(sd.eps * sd.eps)
    sd.eps ~ dunif(0, 10)

for (j in 1:10) {
 
  eps[j] ~ dnorm(0, prec.eps) 

  for (i in 1:N){

          ### Model for step lengths
          S[i,j] ~ dweib(shape[i,j], scale[i,j])
          shape[i,j] <- exp(beta0 + beta1*Sound[i] + beta2*Diel[i] + eps[j])    
                                                       
          ### Model for angles
          phi[i,j] <- - kappa_hat * cos(A[i,j]-mu) + log(2*pi*IO_hat) + C
  
          scale[i,j] <- b2
          dummy[i,j] ~ dpois(phi[i,j])
    }
  } 
    
    k ~ dcat(p[])
    
    for(j in 1:n_kappas) {
      p[j] <- 1/n_kappas  # creating vector p = {1/n_kappas,...}
    }
    
    kappa_hat <- kappas[k]
    IO_hat <- I0s[k]
    
}", fill = TRUE)

#------ Run the JAGS model ---------
parms <- c("beta0","beta1","beta2","shape","scale","mu","eps","phi") #"rho""mean_n",,"rho",

# MCMC settings
ni <- 1000
nt <- 20
nb <- 1000
nc <- 3
#rm(jags.out)
#Call JAGS from R using the runjags package
start = Sys.time()
jags.out <- run.jags(model = "HMMmodelallfish.txt",  monitor = parms,  #inits=inits1,
                     data = jags.dat2, n.chains = nc, burnin = nb, sample = ni, adapt = 1000,
                     thin = nt, method="parallel")
elapsed = Sys.time() - start
elapsed

SumStats <-summary(jags.out, vars = c("beta0","beta1","beta2","mu","eps"))#"beta2","shape","scale","rho",
SumStats

plot(jags.out, vars = c("beta0","beta1","beta2","mu","eps")) #"beta2""rho",























