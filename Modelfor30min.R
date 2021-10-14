#---------------------------------------------------------------------------------
library(animation)
library(lubridate)
library(circular)
library(runjags)
library(tidyverse)
library(StreamMetabolism)
library(mvtnorm)
library(plyr)
library(sp)
library(rgdal)
library(maps)
library(momentuHMM)

# Matt added this line 8/3/2021
usgs.path <- 'C:/Users/RDGRLMPR/Documents/Carp_Model/USGS'

#Map to the animation directory
LocsDat <- read.csv(paste0(usgs.path, "/Supplementary Files", "/CERC2018Positions.csv")) #AEH_18_CERCSOUND_01_GPS
names(LocsDat)[c(1, 2)] <- c("x", "y")

SPK <- subset(LocsDat, Type=='SPK')
KET <- subset(LocsDat, Type=='KET')
BND <- subset(LocsDat, Type=='BND')
HYD <- subset(LocsDat, Type=='HYD')

# Bnd <- matrix(NA, 2, 2)
# Bnd[1,1] <- range(BND$x)[1]
# Bnd[2,1] <- range(BND$x)[2]
# Bnd[1,2] <- range(BND$y)[1]
# Bnd[2,2] <- range(BND$y)[2]
# colnames(Bnd)<- c("x","y")i.options(outdir = getwd(), convert = 'C:/Program Files/ImageMagick-7.1.0-Q16-HDRI/convert.exe')

#Set working directory
setwd(paste0(usgs.path, "/Carp Pond Analysis"))

#read in the pond position data


#--------------------------------------------------------------------------------
#Sound data
SoundDat <- read.csv('C:/Users/RDGRLMPR/Documents/Carp/Master_Sound_Tag_20200925.csv', stringsAsFactors=FALSE)
TrialNum <- sort(unique(SoundDat$Trial))
PondNum <- sort(unique(SoundDat$Pond))# Matt added the following two lines 8/3/2021
colnames(SoundDat)[1] <- 'Trial'
SoundDat$DT <- paste(SoundDat$DOY, SoundDat$LocalTime..CT.)
SoundDat$locTimes <- as.POSIXct(SoundDat$DT, format="%m/%d/%Y %H:%M:%S", tz = "America/Chicago")
# Matt added the following line on 8/3/2021
SoundDat <- SoundDat %>% relocate(locTimes, .after = 'Pond')
SoundDat$DT <- NULL
SoundDat$Sound[SoundDat$Sound=="ON "] <- "ON"
SoundDat <- SoundDat[order(SoundDat$locTimes),]
#----------------------------------------------------------------------------------
#Water temperature data
TData <- read.csv(paste0(usgs.path, "/Supplementary Files/2018CERC_WaterTemperature_AllTrialsPonds.csv"), stringsAsFactors=F)
TData$DT <- as.POSIXct(TData$DateTime,format="%Y-%m-%d %H:%M:%S")
TData$DateTime <- NULL

# TrialNum <- c(1, 2, 3, 4, 5)
# PondNum <- c(26, 27, 30, 31)


#-----------------------------------------------------------------------------------------------------------------------
#load the fitted model parameters
for (t in TrialNum) {
  for (p in PondNum) {
    # John Plumb says:
  #Note: because the correlated random walk is finicky about starting values, I had to run each trial/pond combination separately. 
    # Matt agrees, but the issue appears to go away if you make 'attempts = n' for n > 1 in the crawlWrap() function below.
  #These outputs are saved to the directory and then read back into R for fitting the movement model.
  #Set the trail and pond to run
    # t<-1
    # p<-26
    
    #Sound data
    Sdat <- subset(SoundDat, Sound=="ON" & Trial==t & Pond==p)
    Sdat <- Sdat[order(Sdat$locTimes),]
    
    #Load the processed data
    load(paste0(usgs.path, '/Carp Pond Analysis/ProcessData_', t, '_pond_', p, '.RDATA'))
    
    AllData$Trial <- t
    AllData$Pond <- p

    #assign treatments
    if (t == 1) {
      if (p == 31) {
        AllData$trmt <- 'ChirpSaw'
      } else if (p == 30) {
        AllData$trmt <- 'BoatM'
      } else if (p == 27) {
        AllData$trmt <- 'ChirpSquare'
      } else if (p == 26) {
        AllData$trmt <- 'Control'
      }
    }
    if (t == 2) {
      if (p == 31) {
        AllData$trmt <- 'Control'
      } else if (p == 30) {
        AllData$trmt <- 'ChirpSaw'
      } else if (p == 27) {
        AllData$trmt <- 'BoatM'
      } else if (p == 26) {
        AllData$trmt <- 'ChirpSquare'
      }
    }
    if (t == 3) {
      if (p == 31) {
        AllData$trmt <- 'ChirpSquare'
      } else if (p == 30) {
        AllData$trmt <- 'Control'
      } else if (p == 27) {
        AllData$trmt <- 'ChirpSaw'
      } else if (p == 26) {
        AllData$trmt <- 'BoatM'
      }
    }
    if (t == 4) {
      if (p == 31) {
        AllData$trmt <- 'BoatM'
      } else if (p == 30) {
        AllData$trmt <- 'ChirpSquare'
      } else if (p == 27) {
        AllData$trmt <- 'Control'
      } else if (p == 26) {
        AllData$trmt <- 'ChirpSaw'
      }
    }
    if (t == 5) {
      if (p == 31) {
        AllData$trmt <- 'ChirpSquare'
      } else if (p == 30) {
        AllData$trmt <- 'ChirpSaw'
      } else if (p == 27) {
        AllData$trmt <- 'BoatM'
      } else if (p == 26) {
        AllData$trmt <- 'Control'
      }
    }
    
    #Now subset dataset to just time before and after specified time interval
      # the processed data file has lots of NA rows (one row per second, it seems)
    crawldat0 <- subset(AllData, !is.na(Easting))
      # the following line picks out the first time-on time; is this the right thing to do?
    crawldat0$OnDT <- as.POSIXct(Sdat$locTimes[length(Sdat$locTimes)], format="%Y-%m-%d %H:%M:%S")
    # Matt added the following lines 8/4/2021; the subsequent commented out line seems to have a typo (<<- instead of <=)
      #picks out entries in the 30 minutes before the first sound on-time, but do we ever look at later ones?
    crawldat0 <- subset(crawldat0, DT >= OnDT-1800 & DT <= OnDT+1800)
    # crawldat0 <- subset(crawldat0, DT >= OnDT-1800 & DT <<- OnDT+1800)
    
    # do we need this line?
    crawldat0 <- crawldat0[,c("ID","Easting","Northing","DT","Trial","Pond","trmt","Sound")]
    
    TagCodes <- unique(crawldat0$ID)
    
    rawdat <- crawldat0[, c("ID","Easting","Northing","DT")]
    colnames(rawdat)<-c('ID','x','y','time')
    
    # it may be useful to mess with this line
    inits <- c(2, 0.001) #set initial values for CRWM
  
    #Fit the correlated random walk Model
    tempDat0 <- crawlWrap(obsData=rawdat, timeStep="6 sec",
                          theta=inits, fixPar=c(NA, NA), attempts=20)
    
    # Add date-time column called 'DT'
    tempDat0$crwPredict$DT <- tempDat0$crwPredict$time
    
    # only keep the necessary covariates
    covDat <- AllData[,c("ID","DT","Trial","Pond","trmt","Sound")]
    
    # merge the CRW data with covariate infor from telemetry
    tempDat0$crwPredict <- merge(tempDat0$crwPredict, covDat, by=c("ID","DT"))
    
    #prepare data for input into movement model and define covariates
    ModDat <- prepData(data=tempDat0, covNames=c("Trial", "Pond", "trmt", "Sound"))
    
    #output to directory
    save(ModDat, file = paste0("FishDat_Trial_", t,"_Pond_", p,".RDATA"))

    rm(tempDat0)
    rm(AllData)
    rm(covDat)
    rm(ModDat)
#End of run
 }
}
-----------------------------------------------------------
# Read in all correlated random walk outputs files 
# loop through all trials and ponds
TrialNum <- c(1,2,3,4,5)  
PondNum <- c(26,27,30,31) 
IDs <- matrix(NA, nrow=20, ncol=14)
count<-0
for (t in TrialNum) {
  for (p in PondNum) {
    
    # load the fitted CRW file
    load(paste0(usgs.path, "/Carp Pond Analysis/FishDat_Trial_",t,"_Pond_",p,".RDATA"))
    # get only sound-on times for the particular trial and pond
    Sdat <- subset(SoundDat, Sound=="ON" & Trial==t & Pond==p)
    # order from oldest to newest
    Sdat <- Sdat[order(Sdat$locTimes),]
    
    # grabs the first sound on time
    OnDT <- as.POSIXct(Sdat$locTimes[length(Sdat$locTimes)], format="%Y-%m-%d %H:%M:%S")
    
    # gets the relevant temperature data
    tdat <- subset(TData, Trial==t & Pond==p)
    # order from oldest to newest
    tdat <- tdat[order(tdat$DT),]
    
    # change the time frame below to include more or less data (1800 = 60 secs/min * 30 min)
    tdat0 <- subset(tdat, DT >= OnDT-1800 & DT <= OnDT+1800) 
    
    ##### unknown, perhaps take the mean time for the day?
    MuT <- ddply(tdat0, .(Trial, Pond),
                summarize,
                TempC = mean(Temp_C, na.rm=T))
    
    # Label loaded model data with correct sound on/off info
    ModDat$Sound[ModDat$DT<OnDT] <- "0ff"
    ModDat$Sound[ModDat$DT>=OnDT] <- "On"
    # compute distance to speaker
    ModDat$Dist2SPK <- sqrt((SPK$x[SPK$Pond==p] - ModDat$x)^2 + 
                            (SPK$y[SPK$Pond==p] - ModDat$y)^2)
    
    # incorporate temperature
    ModDat$TempC <- MuT$TempC
    
    # just makes a new copy?
    ModDat0 <- ModDat 
 
    # iteratively fill a dataframe with trial/pond tag
    count<-count+1
    IDs[count,1] <- t
    IDs[count,2] <- p
    IDs[count,3:(2+length(unique(ModDat0$ID)))] <- as.numeric(as.character(unique(ModDat0$ID)))[1:length(unique(ModDat0$ID))]
   
    #just keeps track of codes progress
    print(paste0("Trial_",t,"Pond_",p,"Nfish_", length(unique(ModDat0$ID))))

    #Compile all one-hour tracks into one data set for model fitting
    if (p==PondNum[1]) { 
      FishData0 <- ModDat0
    }else{
      FishData0 <- rbind(FishData0, ModDat0)  
    }
  }   
  if (t==TrialNum[1]) { 
    FishData <- FishData0
  }else{
    FishData <- rbind(FishData, FishData0)  
  }
}
#Evaluate the step lengths and angles
# Doesn't this mix all trials/ponds into one metric? How can it be useful (note that it's not used later)
acf(FishData$step[!is.na(FishData$step)],lag.max=300, xlim=c(0,150))

Exdat <- FishData[FishData$Dist2SPK>40,]
Exdat <- Exdat[!is.na(Exdat$step),]
Exdat <- Exdat[Exdat$ID==2047.311,]
# 
 plot(FishData$x[FishData$ID==2047.311], FishData$y[FishData$ID==2047.311],
      xlim=c(min(FishData$x[FishData$ID==2047.311], SPK$x[SPK$Pond==p], na.rm=T),max(FishData$x[FishData$ID==2047.311], SPK$x[SPK$Pond==p], na.rm=T)),
      ylim=c(min(FishData$y[FishData$ID==2047.311], SPK$y[SPK$Pond==p], na.rm=T),max(FishData$y[FishData$ID==2047.311], SPK$y[SPK$Pond==p], na.rm=T)))
 points(Exdat$x, Exdat$y, col='black', cex=2)
 points(SPK$x[SPK$Pond==p], SPK$y[SPK$Pond==p], col='magenta', cex=1.5, pch=19)
 
p<-31
 
 
 plot(AllData$Easting[AllData$ID==2047.311], AllData$Northing[AllData$ID==2047.311],
      xlim=c(min(AllData$Easting[AllData$ID==2047.311],Exdat$x, na.rm=T),max(AllData$Easting[AllData$ID==2047.311],Exdat$x, na.rm=T)),
      ylim=c(min(AllData$Northing[AllData$ID==2047.311],Exdat$y, na.rm=T),max(AllData$Northing[AllData$ID==2047.311],Exdat$y, na.rm=T)))
 points(Exdat$x, Exdat$y, col='red', cex=2)
# 
# plot(FishData[FishData$Trial!=4,]$x,FishData[FishData$Trial!=4,]$Dist2SPK)


tDat<- FishData[FishData$Trial!=4,]
mean(tDat$step, na.rm=T)
sd(tDat$step, na.rm=T)

mean(tDat$angle, na.rm=T)
sd(tDat$angle, na.rm=T)

boxplot(tDat$step~as.factor(tDat$Sound)+as.factor(tDat$trmt), ylim=c(0,10), notch=T)

plot(tDat$step~tDat$ID)

table(FishData$Sound,FishData$ID, useNA="always")
table(FishData$trmt,FishData$Sound, useNA="always")

boxplot(FishData$Dist2SPK~FishData$trmt)


library(plyr)
#summarize output from correlated random walk
Sumfish = ddply(FishData, .(ID, Trial, Pond, trmt, Sound),
                summarize,
                StepMu = mean(step, na.rm=T),
                StepSD = sd(step, na.rm=T),
                AngleMu = mean(angle, na.rm=T),
                AngleSD = sd(angle, na.rm=T),
                DistMu = mean(Dist2SPK, na.rm=T),
                DistSD = sd(Dist2SPK, na.rm=T))



#Remove trial 4
tDat<- Sumfish[Sumfish$Trial!=4,]

windows(15,13)
par(mfrow = c(3,1), mar = c(5, 4, 0, 2), oma = c(2, 5, 1, 5))
boxplot(tDat$StepMu~as.factor(tDat$Sound)+as.factor(tDat$trmt), ylim=c(0,2), cex.lab=1.5,
        xlab = NA,
        ylab = "Mean step length (m)")

boxplot(tDat$AngleMu~as.factor(tDat$Sound)+as.factor(tDat$trmt), ylim=c(-0.5,0.5), cex.lab=1.5,
        xlab = NA,
        ylab = "Mean turning angle (radians)")

boxplot(tDat$DistMu~as.factor(tDat$Sound)+as.factor(tDat$trmt), ylim=c(0,30), cex.lab=1.5,
        xlab = "Sound operation and type",
        ylab = "Distance from speaker (m)")



#Summary stats without trial 4
Sumfish2 = ddply(FishData[FishData$Trial!=4,], .(trmt, Sound),
                summarize,
                StepMu = mean(step, na.rm=T),
                StepSD = sd(step, na.rm=T),
                AngleMu = mean(angle, na.rm=T),
                AngleSD = sd(angle, na.rm=T),
                DistMu = mean(Dist2SPK, na.rm=T),
                DistSD = sd(Dist2SPK, na.rm=T))

Sumfish2
#summary stats
FishData$INDX <- 1
Sumfish3 = ddply(FishData, .(ID, Trial, Pond, trmt),
                summarize,
                Nobs = sum(INDX, na.rm=T),
                MinDT = min(DT, na.rm=T),
                maxDT = max(DT, na.rm=T))

Sumfish3
write.table(Sumfish3,paste0(usgs.path, "/Carp Pond Analysis/SummaryModFish.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA")
head(FishData)
IDs
#-----------------------------------------------------------  


######################################################################################################
#Now fit the movement model
#---------------------------------------------------------------------
#----Step 1
# label states
stateNames <- c("exploratory","encamped") #
# distributions for observation processes
dist = list(step = "gamma", angle = "vm") #"wrpcauchy"
# initial parameters
Par0_m1 <- list(step=c(2,1,2,1,0,0), angle=c(0.004,0.004,0.002,0.002))
# fit model
m1 <- fitHMM(data = FishData[FishData$Trial!=4,], 
             nbStates = 2, 
             dist = dist, 
             Par0 = Par0_m1,
             estAngleMean = list(angle=TRUE), 
             stateNames = stateNames)

#estimates of the state for each observation
states <- viterbi(m1)
#table(states)/nrow(FishData[FishData$Trial!=4,])
m1$data$states <- viterbi(m1)

#--------------------------------------------------------------------
#----Step 2
#First convert two factor sound variable into a continuous variable to accommodate fit function
m1$data$trmt[m1$data$trmt=="Control"] <- "AControl"

m1$data <- within(m1$data, {
  Sound <- as.factor(Sound)
  trmt <- as.factor(trmt)
  Pond <- as.factor(Pond)
  Trial <- as.factor(Trial)
  TempC <- as.numeric(TempC)
})
#----------------------------------------------------------------------
#The Full Model
# formula for transition probabilities
formulaF <- ~ TempC + Trial + Pond + Sound * trmt 
# initial parameters (obtained from nested model m1)
Par0_m2 <- getPar0(model=m1, formula=formulaF)
# first fit the transition probability model to obtain resonable starting values for the full model
m2 <- fitHMM(data = m1$data, nbStates = 2, 
             dist = dist, 
             Par0 = Par0_m2$Par, 
             estAngleMean = list(angle=TRUE),
             stateNames = stateNames, 
             beta0 = Par0_m2$beta,
             formula = formulaF)

DM <- list(step = list(mean = ~ TempC + Trial + Pond + Sound * trmt,
                       sd = ~ 1,
                       zeromass = ~ 1),
           angle = list(mean = ~ TempC + Trial + Pond + Sound * trmt,
                        concentration = ~ 1))

# initial parameters (obtained from nested model m2)
Par0_m3 <- getPar0(model=m2, formula=formulaF, DM=DM)
# fit model
FullMod <- fitHMM(data =  m2$data, nbStates = 2, 
             dist = dist, 
             Par0 = Par0_m3$Par, 
             beta0 = Par0_m3$beta, 
             DM = DM, 
             stateNames = stateNames, 
             estAngleMean = list(angle=TRUE),
             formula = formulaF)

#----------------------------------------------------------------------
#The Null Model
# formula for transition probabilities
formula <- ~ 1 
# initial parameters (obtained from nested model m1)
Par0_m2 <- getPar0(model=m1, formula=formula)
# first fit the transition probability model to obtain resonable starting values 
m2 <- fitHMM(data = m1$data, nbStates = 2, 
             dist = dist, 
             Par0 = Par0_m2$Par, 
             estAngleMean = list(angle=TRUE),
             stateNames = stateNames, 
             beta0 = Par0_m2$beta,
             formula = formula)

DM <- list(step = list(mean = ~ 1,
                       sd = ~ 1,
                       zeromass = ~ 1),
           angle = list(mean = ~ 1,
                       concentration = ~ 1))

# initial parameters (obtained from nested model m2)
Par0_m3 <- getPar0(model=m2, formula=formula, DM=DM)
# fit model
NullMod <- fitHMM(data =  m2$data, nbStates = 2, 
                  dist = dist, 
                  Par0 = Par0_m3$Par, 
                  beta0 = Par0_m3$beta, 
                  DM = DM, 
                  stateNames = stateNames, 
                  estAngleMean = list(angle=TRUE),
                  formula = formula)
#----------------------------------------------------------------------
#The full Additive Model
# formula for transition probabilities
formulaA <- ~ TempC + Trial + Pond + Sound + trmt 

# initial parameters (obtained from nested model m1)
Par0_m2A <- getPar0(model=m1, formula=formulaA)

# first fit the transition probability model to obtain resonable starting values for the full model
m2A <- fitHMM(data = m1$data, nbStates = 2, 
             dist = dist, 
             Par0 = Par0_m2A$Par, 
             estAngleMean = list(angle=TRUE),
             stateNames = stateNames, 
             beta0 = Par0_m2A$beta,
             formula = formulaA)

DMA <- list(step = list(mean = ~ TempC + Trial + Pond + Sound + trmt,
                       sd = ~ 1,
                       zeromass = ~ 1),
           angle = list(mean = ~ TempC + Trial + Pond + Sound + trmt,
                       concentration = ~ 1))

# initial parameters (obtained from nested model m2)
Par0_m3A <- getPar0(model=m2A, formula=formulaA, DM=DMA)

# fit model
AddMod <- fitHMM(data =  m2A$data, nbStates = 2, 
                  dist = dist, 
                  Par0 = Par0_m3A$Par, 
                  beta0 = Par0_m3A$beta, 
                  DM = DMA, 
                  stateNames = stateNames, 
                  estAngleMean = list(angle=TRUE),
                  formula = formulaA)

#----------------------------------------------------------------------
# The full model without temperature 
# formula for transition probabilities
formula3 <- ~ Trial + Pond + Sound * trmt 

# initial parameters (obtained from nested model m1)
Par0_m23 <- getPar0(model=m1, formula=formula3)

# first fit the transition probability model to obtain resonable starting values for the full model
m23 <- fitHMM(data = m1$data, nbStates = 2, 
              dist = dist, 
              Par0 = Par0_m23$Par, 
              estAngleMean = list(angle=TRUE),
              stateNames = stateNames, 
              beta0 = Par0_m23$beta,
              formula = formula3)

DM3 <- list(step = list(mean = ~ Trial + Pond + Sound * trmt,
                        sd = ~ 1,
                        zeromass = ~ 1),
            angle = list(mean = ~ Trial + Pond + Sound * trmt,
                         concentration = ~ 1))

# initial parameters (obtained from nested model m2)
Par0_m33 <- getPar0(model=m23, formula=formula3, DM=DM3)

# fit model
Mod3 <- fitHMM(data =  m23$data, nbStates = 2, 
                 dist = dist, 
                 Par0 = Par0_m33$Par, 
                 beta0 = Par0_m33$beta, 
                 DM = DM3, 
                 stateNames = stateNames, 
                 estAngleMean = list(angle=TRUE),
                 formula = formula3)

#----------------------------------------------------------------
# formula for transition probabilities
formula4 <- ~ TempC + Sound * trmt 
# initial parameters (obtained from nested model m1)
Par0_m24 <- getPar0(model=m1, formula=formula4)
# fit model
m24 <- fitHMM(data = m1$data, nbStates = 2, 
              dist = dist, 
              Par0 = Par0_m24$Par, 
              estAngleMean = list(angle=TRUE),
              stateNames = stateNames, 
              beta0 = Par0_m24$beta,
              formula = formula4)

DM4 <- list(step = list(mean = ~ TempC + Sound * trmt,
                        sd = ~ 1,
                        zeromass = ~ 1),
            angle = list(mean = ~ TempC + Sound * trmt,
                         concentration = ~ 1))

# initial parameters (obtained from nested model m2)
Par0_m34 <- getPar0(model=m24, formula=formula4, DM=DM4)

# fit model
Mod4 <- fitHMM(data =  m24$data, nbStates = 2, 
               dist = dist, 
               Par0 = Par0_m34$Par, 
               beta0 = Par0_m34$beta, 
               DM = DM4, 
               stateNames = stateNames, 
               estAngleMean = list(angle=TRUE),
               formula = formula4)
#----------------------------------------------------------------
# formula for transition probabilities
formula5 <- ~ TempC + Trial + Sound * trmt 
# initial parameters (obtained from nested model m1)
Par0_m25 <- getPar0(model=m1, formula=formula5)
# fit model
m25 <- fitHMM(data = m1$data, nbStates = 2, 
              dist = dist, 
              Par0 = Par0_m25$Par, 
              estAngleMean = list(angle=TRUE),
              stateNames = stateNames, 
              beta0 = Par0_m25$beta,
              formula = formula5)

DM5 <- list(step = list(mean = ~ TempC + Trial + Sound * trmt,
                        sd = ~ 1,
                        zeromass = ~ 1),
            angle = list(mean = ~ TempC + Trial + Sound * trmt,
                         concentration = ~ 1))

# initial parameters (obtained from nested model m2)
Par0_m35 <- getPar0(model=m25, formula=formula5, DM=DM5)

# fit model
Mod5 <- fitHMM(data =  m25$data, nbStates = 2, 
               dist = dist, 
               Par0 = Par0_m35$Par, 
               beta0 = Par0_m35$beta, 
               DM = DM5, 
               stateNames = stateNames, 
               estAngleMean = list(angle=TRUE),
               formula = formula5)
#---------------------------------------------------------------
# formula for transition probabilities
formula6 <- ~ TempC + Pond + Sound * trmt 
# initial parameters (obtained from nested model m1)
Par0_m26 <- getPar0(model=m1, formula=formula6)
# fit model
m26 <- fitHMM(data = m1$data, nbStates = 2, 
              dist = dist, 
              Par0 = Par0_m26$Par, 
              estAngleMean = list(angle=TRUE),
              stateNames = stateNames, 
              beta0 = Par0_m26$beta,
              formula = formula6)

DM6 <- list(step = list(mean = ~ TempC + Pond + Sound * trmt,
                        sd = ~ 1,
                        zeromass = ~ 1),
            angle = list(mean = ~ TempC + Pond + Sound * trmt,
                         concentration = ~ 1))

# initial parameters (obtained from nested model m2)
Par0_m36 <- getPar0(model=m26, formula=formula6, DM=DM6)

# fit model
Mod6 <- fitHMM(data =  m26$data, nbStates = 2, 
               dist = dist, 
               Par0 = Par0_m36$Par, 
               beta0 = Par0_m36$beta, 
               DM = DM6, 
               stateNames = stateNames, 
               estAngleMean = list(angle=TRUE),
               formula = formula6)
#---------------------------------------------------------------
# formula for transition probabilities
formula7 <- ~ Trial 
# initial parameters (obtained from nested model m1)
Par0_m27 <- getPar0(model=m1, formula=formula7)
# fit model
m27 <- fitHMM(data = m1$data, nbStates = 2, 
              dist = dist, 
              Par0 = Par0_m27$Par, 
              estAngleMean = list(angle=TRUE),
              stateNames = stateNames, 
              beta0 = Par0_m27$beta,
              formula = formula7)

DM7 <- list(step = list(mean = ~ Trial,
                        sd = ~ 1,
                        zeromass = ~ 1),
            angle = list(mean = ~ Trial,
                         concentration = ~ 1))

# initial parameters (obtained from nested model m2)
Par0_m37 <- getPar0(model=m27, formula=formula7, DM=DM7)

# fit model
Mod7 <- fitHMM(data =  m27$data, nbStates = 2, 
               dist = dist, 
               Par0 = Par0_m37$Par, 
               beta0 = Par0_m37$beta, 
               DM = DM7, 
               stateNames = stateNames, 
               estAngleMean = list(angle=TRUE),
               formula = formula7)
#---------------------------------------------------------------
# formula for transition probabilities
formula8 <- ~ Pond 
# initial parameters (obtained from nested model m1)
Par0_m28 <- getPar0(model=m1, formula=formula8)
# fit model
m28 <- fitHMM(data = m1$data, nbStates = 2, 
              dist = dist, 
              Par0 = Par0_m28$Par, 
              estAngleMean = list(angle=TRUE),
              stateNames = stateNames, 
              beta0 = Par0_m28$beta,
              formula = formula8)

DM8 <- list(step = list(mean = ~ Pond,
                        sd = ~ 1,
                        zeromass = ~ 1),
            angle = list(mean = ~ Pond,
                         concentration = ~ 1))

# initial parameters (obtained from nested model m2)
Par0_m38 <- getPar0(model=m28, formula=formula8, DM=DM8)

# fit model
Mod8 <- fitHMM(data =  m28$data, nbStates = 2, 
               dist = dist, 
               Par0 = Par0_m38$Par, 
               beta0 = Par0_m38$beta, 
               DM = DM8, 
               stateNames = stateNames, 
               estAngleMean = list(angle=TRUE),
               formula = formula8)
#---------------------------------------------------------------
# formula for transition probabilities
formula9 <- ~ Sound * trmt  
# initial parameters (obtained from nested model m1)
Par0_m29 <- getPar0(model=m1, formula=formula9)
# fit model
m29 <- fitHMM(data = m1$data, nbStates = 2, 
              dist = dist, 
              Par0 = Par0_m29$Par, 
              estAngleMean = list(angle=TRUE),
              stateNames = stateNames, 
              beta0 = Par0_m29$beta,
              formula = formula9)

DM9 <- list(step = list(mean = ~ Sound * trmt,
                        sd = ~ 1,
                        zeromass = ~ 1),
            angle = list(mean = ~ Sound * trmt,
                         concentration = ~ 1))

# initial parameters (obtained from nested model m2)
Par0_m39 <- getPar0(model=m29, formula=formula9, DM=DM9)

# fit model
Mod9 <- fitHMM(data =  m29$data, nbStates = 2, 
               dist = dist, 
               Par0 = Par0_m39$Par, 
               beta0 = Par0_m39$beta, 
               DM = DM9, 
               stateNames = stateNames, 
               estAngleMean = list(angle=TRUE),
               formula = formula9)
#---------------------------------------------------------------
# formula for transition probabilities
formula10 <- ~ TempC 
# initial parameters (obtained from nested model m1)
Par0_m210 <- getPar0(model=m1, formula=formula10)
# fit model
m210 <- fitHMM(data = m1$data, nbStates = 2, 
              dist = dist, 
              Par0 = Par0_m210$Par, 
              estAngleMean = list(angle=TRUE),
              stateNames = stateNames, 
              beta0 = Par0_m210$beta,
              formula = formula10)

DM10 <- list(step = list(mean = ~ TempC,
                        sd = ~ 1,
                        zeromass = ~ 1),
            angle = list(mean = ~ TempC,
                         concentration = ~ 1))

# initial parameters (obtained from nested model m2)
Par0_m310 <- getPar0(model=m210, formula=formula10, DM=DM10)

# fit model
Mod10 <- fitHMM(data =  m210$data, nbStates = 2, 
               dist = dist, 
               Par0 = Par0_m310$Par, 
               beta0 = Par0_m310$beta, 
               DM = DM10, 
               stateNames = stateNames, 
               estAngleMean = list(angle=TRUE),
               formula = formula10)
#--------------------------------------------------------------
#Model selection tables
ModSelect <- data.frame(nrow=11, ncol=4, stringAsFactors=F)
ModSelect[1,2] <- as.character(formulaF)[2]
ModSelect[2,2] <- as.character(formulaA)[2]
ModSelect[3,2] <- as.character(formula3)[2]
ModSelect[4,2] <- as.character(formula4)[2]
ModSelect[5,2] <- as.character(formula5)[2]
ModSelect[6,2] <- as.character(formula6)[2]
ModSelect[7,2] <- as.character(formula7)[2] #Trial only
ModSelect[8,2] <- as.character(formula8)[2] #Pond only
ModSelect[9,2] <- as.character(formula9)[2] #Sound and Trmt interaction
ModSelect[10,2] <- as.character(formula10)[2] #temperature only
ModSelect[11,2] <- as.character(formula)[2]

ModSelect[1,3] <- AIC(FullMod)
ModSelect[2,3] <- AIC(AddMod)
ModSelect[3,3] <- AIC(Mod3)
ModSelect[4,3] <- AIC(Mod4)
ModSelect[5,3] <- AIC(Mod5)
ModSelect[6,3] <- AIC(Mod6)
ModSelect[7,3] <- AIC(Mod7) #Trial only
ModSelect[8,3] <- AIC(Mod8) #Pond only
ModSelect[9,3] <- AIC(Mod9) #Sound and Trmt interaction
ModSelect[10,3] <- AIC(Mod10) #temperature only
ModSelect[11,3] <- AIC(NullMod)

ModSelect[,4] <- min(ModSelect[,3]) - ModSelect[,3]
ModSelect[,1] <- as.numeric(rownames(ModSelect))
colnames(ModSelect)<-c('Number','Model','AIC','DeltaAIC')

ModSelect
#-------------------------------------------------------------------
#pseudo residuals to assess model fit
PResids <- as.data.frame(pseudoRes(FullMod, ncores=1))[-1,]
plot(PResids, ylab="Angle pseudo residual", xlab="Step pseudo residual")

qqnorm(PResids$stepRes[!is.na(PResids$stepRes)], ylim=c(-6,6),xlim=c(-6,6))
abline(0,1, lty=3)
qqnorm(PResids$angleRes, ylim=c(-6,6),xlim=c(-6,6))
abline(0,1, lty=3)

#-------------------------------------------------------------------------------------------------
#Plot the best model
plot(FullMod, plotCI = T, breaks=20)#, covs= data.frame(Sound="On", trmt="ChirpSquare"))#"BoatM" ChirpSquare" "Control"

states <- viterbi(FullMod)
table(states)/nrow(FullMod$data)

FullMod$data$states <- viterbi(FullMod)
# write.table(FullMod$data,paste0(getwd(),"/TopModel_FullMod.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA")
# The following line has caused me untold anguish; I don't know why it exists
# FullMod <-read.table(paste0(getwd(),"/TopModel_FullMod.csv"), header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE) 


# table(m3$data$states, m3$data$Sound, m3$data$trmt)
#--------------------------------------------------------------------

Controlls <- as.matrix(table(FullMod$data$states[FullMod$data$trmt=="AControl"], FullMod$data$Sound[FullMod$data$trmt=="AControl"]))
rownames(Controlls) <- stateNames
colnames(Controlls) <- c("Off", "On")
t(Controlls)

#Conditional on Sound OFF
Controlls/sum(Controlls)
Controlls[,1]/sum(Controlls[,1])
#Conditional on Sound ON
Controlls[,2]/sum(Controlls[,2])

#------
BoatMls <- as.matrix(table(FullMod$data$states[FullMod$data$trmt=="BoatM"], FullMod$data$Sound[FullMod$data$trmt=="BoatM"]))
rownames(BoatMls) <- stateNames
colnames(BoatMls) <- c("Off", "On")
t(BoatMls)
#Conditional on Sound OFF
BoatMls/sum(BoatMls)
BoatMls[,1]/sum(BoatMls[,1])
#Conditional on Sound ON
BoatMls[,2]/sum(BoatMls[,2])

#------
ChirpSquarels <- as.matrix(table(FullMod$data$states[FullMod$data$trmt=="ChirpSquare"], FullMod$data$Sound[FullMod$data$trmt=="ChirpSquare"]))
rownames(ChirpSquarels) <- stateNames
colnames(ChirpSquarels) <- c("Off", "On")
t(ChirpSquarels)

#Conditional on Sound OFF
ChirpSquarels/sum(ChirpSquarels)
ChirpSquarels[,1]/sum(ChirpSquarels[,1])
#Conditional on Sound ON
ChirpSquarels[,2]/sum(ChirpSquarels[,2])

#------
ChirpSawls <- as.matrix(table(FullMod$data$states[FullMod$data$trmt=="ChirpSaw"], FullMod$data$Sound[FullMod$data$trmt=="ChirpSaw"]))
rownames(ChirpSawls) <- stateNames
colnames(ChirpSawls) <- c("Off", "On")
t(ChirpSawls)
#Conditional on Sound OFF
ChirpSawls/sum(ChirpSawls)
ChirpSawls[,1]/sum(ChirpSawls[,1])
#Conditional on Sound ON
ChirpSawls[,2]/sum(ChirpSawls[,2])

#summary for comparison
#off
Controlls[,1]/sum(Controlls[,1])
BoatMls[,1]/sum(BoatMls[,1])
ChirpSawls[,1]/sum(ChirpSawls[,1])
ChirpSquarels[,1]/sum(ChirpSquarels[,1])
#on
Controlls[,2]/sum(Controlls[,2])
BoatMls[,2]/sum(BoatMls[,2])
ChirpSawls[,2]/sum(ChirpSawls[,2])
ChirpSquarels[,2]/sum(ChirpSquarels[,2])

########################################################################
# This code is for graphing the data in .gif movie format
########################################################################
# Tnum=5
# FishX <- FullMod
# FishX$Sound <- as.numeric(as.factor(FishX$Sound))
# FishX$trmt <- as.character(FishX$trmt)
# FishX$trmt[FishX$trmt=="AControl"] <- "Control"
#
# #GraphDat<-FishX
# GraphDat <- subset(FishX, Trial==Tnum)
# Graph26 <- subset(GraphDat, Pond==26)
# Graph27 <- subset(GraphDat, Pond==27)
# Graph30 <- subset(GraphDat, Pond==30)
# Graph31 <- subset(GraphDat, Pond==31)
#
# ID26 <- unique(Graph26$ID)
# ID27 <- unique(Graph27$ID)
# ID30 <- unique(Graph30$ID)
# ID31 <- unique(Graph31$ID)
#
# SPK26 <- subset(SPK, Pond==26)
# SPK27 <- subset(SPK, Pond==27)
# SPK30 <- subset(SPK, Pond==30)
# SPK31 <- subset(SPK, Pond==31)
#
# BND26 <- subset(BND, Pond==26)
# BND27 <- subset(BND, Pond==27)
# BND30 <- subset(BND, Pond==30)
# BND31 <- subset(BND, Pond==31)
#
# KET26 <- subset(KET, Pond==26)
# KET27 <- subset(KET, Pond==27)
# KET27$y[KET27$y<min(BND27$y)]<-min(BND27$y)
# KET30 <- subset(KET, Pond==30)
# KET31 <- subset(KET, Pond==31)
#
# tiffName = paste0('/Trial_', Tnum ,'_fish tracks.tif')
# tiff(paste0(getwd(), tiffName),
#      width = 12, height = 8, units = "in", res= 300)
#
# #windows(12, 8)
# par(mfrow = c(1,1), mar = c(4,4,2,2), oma = c(1.5, 2,0,0))
#         #    windows(12, 8)
#             plot(GraphDat$x, GraphDat$y, pch="", col='white',
#                  ylim=c(min(LocsDat$y),max(LocsDat$y)+7),
#                  xlim=c(min(LocsDat$x),max(LocsDat$x)+7),
#                  ylab='Northing', xlab='Easting', cex.lab=1.5)
#
#             #Speaker locations
#             points(SPK26$x,SPK26$y, col='black', cex=3.5, pch=12)
#             points(SPK27$x,SPK27$y, col='black', cex=3.5, pch=12)
#             points(SPK30$x,SPK30$y, col='black', cex=3.5, pch=12)
#             points(SPK31$x,SPK31$y, col='black', cex=3.5, pch=12)
#
#             #Pond boundary locations
#             lines(BND26$x, BND26$y, type='o', col='brown', lwd=3, pch="")
#             lines(BND27$x, BND27$y, type='o', col='brown', lwd=3, pch="")
#             lines(BND30$x, BND30$y, type='o', col='brown', lwd=3, pch="")
#             lines(BND31$x, BND31$y, type='o', col='brown', lwd=3, pch="")
#
#             #Kettle locations
#             lines(KET26$x, KET26$y, type='o', col='magenta', lwd=3, pch="")
#             lines(KET27$x, KET27$y, type='o', col='magenta', lwd=3, pch="")
#             lines(KET30$x, KET30$y, type='o', col='magenta', lwd=3, pch="")
#             lines(KET31$x, KET31$y, type='o', col='magenta', lwd=3, pch="")
#
#             legend('topright', c(paste0("Trial ", unique(GraphDat$Trial)), "Off", "On", "Speaker location", "Kettle location","Pond boundary"), lty = c(NA,1,1,NA,1,1),
#                    lwd=c(NA,3,3,NA,3,3), col=c(NA,1,2,'black','magenta','brown'), pch = c(NA,NA,NA,12,NA,NA), cex=1, bty="n", y.intersp = 1)
#
#           #--Pond 26
#             lines(Graph26$x[Graph26$ID==ID26[6]][1:300], Graph26$y[Graph26$ID==ID26[6]][1:300], col=Graph26$Sound[Graph26$ID==ID26[6]][1:300], lty=1, pch=Graph26$states[Graph26$ID==ID26[6]][1:300], lwd=0.5, cex=1, type='o')
#             lines(Graph26$x[Graph26$ID==ID26[6]][301:595], Graph26$y[Graph26$ID==ID26[6]][301:595], col=Graph26$Sound[Graph26$ID==ID26[6]][301:595], lty=1, lwd=0.5, cex=1, type='o', pch=Graph26$states[Graph26$ID==ID26[6]][301:595])
#
#             lines(Graph26$x[Graph26$ID==ID26[7]][1:300], Graph26$y[Graph26$ID==ID26[7]][1:300], col=Graph26$Sound[Graph26$ID==ID26[7]][1:300], lty=1, pch=Graph26$states[Graph26$ID==ID26[7]][1:300], lwd=0.5, cex=1, type='o')
#             lines(Graph26$x[Graph26$ID==ID26[7]][301:595], Graph26$y[Graph26$ID==ID26[7]][301:595], col=Graph26$Sound[Graph26$ID==ID26[7]][301:595], lty=1, lwd=0.5, cex=1, type='o', pch=Graph26$states[Graph26$ID==ID26[7]][301:595])
#
#             lines(Graph26$x[Graph26$ID==ID26[4]][1:300], Graph26$y[Graph26$ID==ID26[4]][1:300], col=Graph26$Sound[Graph26$ID==ID26[4]][1:300], lty=1, pch=Graph26$states[Graph26$ID==ID26[4]][1:300], lwd=0.5, cex=1, type='o')
#             lines(Graph26$x[Graph26$ID==ID26[4]][301:595], Graph26$y[Graph26$ID==ID26[4]][301:595], col=Graph26$Sound[Graph26$ID==ID26[4]][301:595], lty=1, lwd=0.5, cex=1, type='o', pch=Graph26$states[Graph26$ID==ID26[4]][301:595])
#             text(562461,4307300, "Pond 26")
#             text(562461,4307298, paste0(unique(Graph26$trmt)))
#
#             #--Pond 27
#             lines(Graph27$x[Graph27$ID==ID27[1]][1:300], Graph27$y[Graph27$ID==ID27[1]][1:300], col=Graph27$Sound[Graph27$ID==ID27[1]][1:300], lty=1, pch=Graph27$states[Graph27$ID==ID27[1]][1:300], lwd=0.5, cex=1, type='o')
#             lines(Graph27$x[Graph27$ID==ID27[1]][301:595], Graph27$y[Graph27$ID==ID27[1]][301:595], col=Graph27$Sound[Graph27$ID==ID27[1]][301:595], lty=1, pch=Graph27$states[Graph27$ID==ID27[1]][301:595], lwd=0.5, cex=1, type='o')
#
#             lines(Graph27$x[Graph27$ID==ID27[2]][1:300], Graph27$y[Graph27$ID==ID27[2]][1:300], col=Graph27$Sound[Graph27$ID==ID27[2]][1:300], lty=1, pch=Graph27$states[Graph27$ID==ID27[2]][1:300], lwd=0.5, cex=1, type='o')
#             lines(Graph27$x[Graph27$ID==ID27[2]][301:595], Graph27$y[Graph27$ID==ID27[2]][301:595], col=Graph27$Sound[Graph27$ID==ID27[2]][301:595], lty=1, pch=Graph27$states[Graph27$ID==ID27[2]][301:595], lwd=0.5, cex=1, type='o')
#
#             lines(Graph27$x[Graph27$ID==ID27[3]][1:300], Graph27$y[Graph27$ID==ID27[3]][1:300], col=Graph27$Sound[Graph27$ID==ID27[3]][1:300], lty=1, pch=Graph27$states[Graph27$ID==ID27[3]][1:300], lwd=0.5, cex=1, type='o')
#             lines(Graph27$x[Graph27$ID==ID27[3]][301:595], Graph27$y[Graph27$ID==ID27[3]][301:595], col=Graph27$Sound[Graph27$ID==ID27[3]][301:595], lty=1, pch=Graph27$states[Graph27$ID==ID27[3]][301:595], lwd=0.5, cex=1, type='o')
#             text(562428,4307300, "Pond 27")
#             text(562428,4307298, paste0(unique(Graph27$trmt)))
#             #--Pond 30
#             lines(Graph30$x[Graph30$ID==ID30[1]][1:300], Graph30$y[Graph30$ID==ID30[1]][1:300], col=Graph30$Sound[Graph30$ID==ID30[1]][1:300], lty=1, pch=Graph30$states[Graph30$ID==ID30[1]][1:300], lwd=0.5, cex=1, type='o')
#             lines(Graph30$x[Graph30$ID==ID30[1]][301:595], Graph30$y[Graph30$ID==ID30[1]][301:595], col=Graph30$Sound[Graph30$ID==ID30[1]][301:595], lty=1, pch=Graph30$states[Graph30$ID==ID30[1]][301:595], lwd=0.5, cex=1, type='o')
#
#             lines(Graph30$x[Graph30$ID==ID30[2]][1:300], Graph30$y[Graph30$ID==ID30[2]][1:300], col=Graph30$Sound[Graph30$ID==ID30[2]][1:300], lty=1, pch=Graph30$states[Graph30$ID==ID30[2]][1:300], lwd=0.5, cex=1, type='o')
#             lines(Graph30$x[Graph30$ID==ID30[2]][301:595], Graph30$y[Graph30$ID==ID30[2]][301:595], col=Graph30$Sound[Graph30$ID==ID30[2]][301:595], lty=1, pch=Graph30$states[Graph30$ID==ID30[2]][301:595], lwd=0.5, cex=1, type='o')
#
#             lines(Graph30$x[Graph30$ID==ID30[3]][1:300], Graph30$y[Graph30$ID==ID30[3]][1:300], col=Graph30$Sound[Graph30$ID==ID30[3]][1:300], lty=1, pch=Graph30$states[Graph30$ID==ID30[3]][1:300], lwd=0.5, cex=1, type='o')
#             lines(Graph30$x[Graph30$ID==ID30[3]][301:595], Graph30$y[Graph30$ID==ID30[3]][301:595], col=Graph30$Sound[Graph30$ID==ID30[3]][301:595], lty=1, pch=Graph30$states[Graph30$ID==ID30[3]][301:595], lwd=0.5, cex=1, type='o')
#             text(562395,4307300, "Pond 30")
#             text(562395,4307298, paste0(unique(Graph30$trmt)))
#             #--Pond 31
#             lines(Graph31$x[Graph31$ID==ID31[1]][1:300], Graph31$y[Graph31$ID==ID31[1]][1:300], col=Graph31$Sound[Graph31$ID==ID31[1]][1:300], lty=1, pch=Graph31$states[Graph31$ID==ID31[1]][1:300], lwd=0.5, cex=1, type='o')
#             lines(Graph31$x[Graph31$ID==ID31[1]][301:595], Graph31$y[Graph31$ID==ID31[1]][301:595], col=Graph31$Sound[Graph31$ID==ID31[1]][301:595], lty=1, pch=Graph31$states[Graph31$ID==ID31[1]][301:595], lwd=0.5, cex=1, type='o')
#
#             lines(Graph31$x[Graph31$ID==ID31[2]][1:300], Graph31$y[Graph31$ID==ID31[2]][1:300], col=Graph31$Sound[Graph31$ID==ID31[2]][1:300], lty=1, pch=Graph31$states[Graph31$ID==ID31[2]][1:300], lwd=0.5, cex=1, type='o')
#             lines(Graph31$x[Graph31$ID==ID31[2]][301:595], Graph31$y[Graph31$ID==ID31[2]][301:595], col=Graph31$Sound[Graph31$ID==ID31[2]][301:595], lty=1, pch=Graph31$states[Graph31$ID==ID31[2]][301:595], lwd=0.5, cex=1, type='o')
#
#             lines(Graph31$x[Graph31$ID==ID31[3]][1:300], Graph31$y[Graph31$ID==ID31[3]][1:300], col=Graph31$Sound[Graph31$ID==ID31[3]][1:300], lty=1, pch=Graph31$states[Graph31$ID==ID31[3]][1:300], lwd=0.5, cex=1, type='o')
#             lines(Graph31$x[Graph31$ID==ID31[3]][301:595], Graph31$y[Graph31$ID==ID31[3]][301:595], col=Graph31$Sound[Graph31$ID==ID31[3]][301:595], lty=1, pch=Graph31$states[Graph31$ID==ID31[3]][301:595], lwd=0.5, cex=1, type='o')
#             text(562365,4307300, "Pond 31")
#             text(562365,4307298, paste0(unique(Graph31$trmt)))
# dev.off()
