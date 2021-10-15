################################################################################################
# To import and model the movement of acoustic-tagged Silver Carp in ponds exposed to various Sound Stimuli
# Date November 26 2019
################################################################################################
#packages
library(lubridate)
library(tidyverse)
library(StreamMetabolism) 
library(mvtnorm)
# Matt added the following line on 8/3/2021
library(dplyr)
library(plyr)

# Matt added the following code
data_path <- "~/Carp Telemetry/Telemetry Data/withDropTags/"

#Set working directory
setwd("~/Carp-Model/Carp Pond Analysis/")
#Please specify the trial, pond number, and scenario
TrialNums <- c(1, 2, 3, 4, 5)
PondNums <- c(26, 27, 30, 31)

#-------------------------------------------------------------------------------------
for (TrialNum in TrialNums) {
  for (PondNum in PondNums) {
    print(paste0("Working on Trial ", TrialNum, " and Pond ", PondNum, "."))
    #assign file names by Trial
    
      if (PondNum == 31){
          system <- 100
      }else if (PondNum ==30){
          system <- 200
      }else if (PondNum == 27){
          system <- 300
      }else if (PondNum == 26){
          system <- 400
      }
      
    realDat <- read.csv(file=paste0(data_path, "TagHistoryFor_AllTags-Trial", TrialNum, "-System",
                                      system, ".csv"), stringsAsFactors=F)
    
    # some datasets have an extra column "X"  
    realDat$X <- NULL

    TagCodes <- unique(realDat$TagCode)
    realDat <- realDat[realDat$Easting!=0,]
    realDat <- realDat[order(realDat$TagCode, realDat$LocalTime),]
    
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
      
      TempDat$DT <- as.POSIXct(TempDat$LocalTime, format="%m/%d/%Y %H:%M:%S", digits=3, tz = "America/Chicago")
      TempDat <- TempDat[order(TempDat$DT),]
      
      TempDat$TimeDiff[2:nrow(TempDat)] <- diff(TempDat$DT)
      
      minDT <- TempDat$DT[1]
      maxDT <- TempDat$DT[nrow(TempDat)]
      
      TempDat <- subset(TempDat, TimeDiff>0)
      rownames(TempDat) <- 1:nrow(TempDat)
      
      TempDat <- subset(TempDat, !is.na(TagCode))
      
      TempDat$ID <- TagCodes[i]
      TempDat$TagCode <- NULL
      
      # minimum datetime for the fish
      FishDT[i,1] <- minDT
      # maximum datetime for the fish
      FishDT[i,2] <- maxDT
      FishN[i,1] <- nrow(TempDat)
      
      if (i==1) {
        AllData0 <- TempDat
      }else{
        AllData0 <- rbind(AllData0, TempDat)
      }
    }
    
    FishDT[,1] <- as.POSIXct(FishDT[,1], origin="1970-01-01", tz = "America/Chicago")
    FishDT[,2] <- as.POSIXct(FishDT[,2], origin="1970-01-01", tz = "America/Chicago")
    
    #-------------------------------------------------------------------------------------
    # Now merge all fish positions in the pond by the second to the sound on off times
    #-------------------------------------------------------------------------------------
    minDT0 <- min(FishDT[,1])
    maxDT0 <- max(FishDT[,2])
    #-------------------------------------------------------------------------------------
    #Read in the Sound on/off times
    SoundDat <- read.csv("~/Carp-Model/Supplementary Files/Master_Sound_Tag_20200925.csv", 
                         stringsAsFactors=FALSE)
    colnames(SoundDat)[1] <- 'Trial'
    SoundDat$DT <- paste(SoundDat$DOY, SoundDat$LocalTime..CT.)
    SoundDat$locTimes <- as.POSIXct(SoundDat$DT, format="%m/%d/%Y %H:%M:%S", tz = "America/Chicago")
    SoundDat <- SoundDat %>% relocate(locTimes, .after = 'Pond')
    SoundDat$DT <- NULL
    SoundDat$Sound[SoundDat$Sound=="ON "] <- "ON"
    SoundDat <- SoundDat[order(SoundDat$locTimes),]
    SoundDat <- subset(SoundDat, Trial==TrialNum & Pond==PondNum)
    
    Offdat <- subset(SoundDat, Sound=="OFF")
    Offdat$Sound  <- NULL
    colnames(Offdat) <- c("Trial", "Pond", "OffDT")
    Offdat <- Offdat[, 1:3]
    
    Ondat <- subset(SoundDat, Sound=="ON")
    Ondat$Sound  <- NULL
    colnames(Ondat) <- c("Trial", "Pond", "OnDT")
    Ondat <- Ondat[, 1:3]
    
    Sdat <- cbind(Ondat, Offdat[,-c(1:2)])
    colnames(Sdat) <- c("Trial", "Pond", "OnDT", "OffDT")
    
    stemp <- Sdat[1,]
    stemp$OffDT <- as.POSIXct("06/01/2018 00:00:00",format="%m/%d/%Y %H:%M:%S", 
                              tz = "America/Chicago")
    stemp$OnDT <- as.POSIXct("06/01/2018 00:00:00",format="%m/%d/%Y %H:%M:%S", 
                             tz = "America/Chicago")
    
    Sdat <- rbind(stemp, Sdat)
    
    #-----------------------------------------------------------------------------------------------
    
    Dates <- as.data.frame(seq.POSIXt(minDT0, maxDT0, by="1 sec"), stringsAsFactors=F)
    colnames(Dates) <- "DT"
    
    Dates$Sound <- NA
    Dates$Sindex <- NA
    Dates$Diel <- NA
    offcounter <- seq(1, nrow(Sdat)*2, by=2)
    oncounter <- seq(2, nrow(Sdat)*2, by=2)
    for(j in 2:nrow(Sdat)){
      Dates$Sound[Dates$DT >= Sdat$OffDT[j-1] & Dates$DT < Sdat$OnDT[j]] <- 0
      Dates$Sound[Dates$DT >= Sdat$OnDT[j] & Dates$DT < Sdat$OffDT[j]]  <- 1
      
      
      Dates$Sindex[Dates$DT >= Sdat$OffDT[j-1] & Dates$DT < Sdat$OnDT[j]] <- offcounter[j-1]
      Dates$Sindex[Dates$DT >= Sdat$OnDT[j] & Dates$DT < Sdat$OffDT[j]]  <- oncounter[j-1]
    } 
    
    Dates$Sound[is.na(Dates$Sound)]  <- 0
    Dates$Sindex[is.na(Dates$Sindex)]  <- offcounter[j]
    
    #-------------------------------------------------------------------------------------
    # Obtain sunrise times
    # Here's the function: sunrise.set(lat, long, date, timezone = "UTC", num.days = 1)
    # Using the location of Pond 26
    SunRS <- sunrise.set(38.9122061924, -92.2795993947,
                         paste0(year(min(AllData0$DT)), '/', month(min(AllData0$DT)), '/', 
                                day(min(AllData0$DT))-1), num.days = 8, 
                                timezone='America/North_Dakota/Center')
    SunRS[,1] <- as.POSIXct(SunRS[,1], origin="1970-01-01", tz = "America/Chicago")
    SunRS[,2] <- as.POSIXct(SunRS[,2], origin="1970-01-01", tz = "America/Chicago")
    SunRS$Date <- as.Date(SunRS[,1])
    
    # Day is 0, night is 1
    for(j in 1:(nrow(SunRS)-1)){
      Dates$Diel[(Dates$DT >= SunRS$sunrise[j]) & (Dates$DT < SunRS$sunset[j])] <- 0
    }
    Dates$Diel[is.na(Dates$Diel)] <- 1
    
    # Merge the environmental conditons
    for(k in unique(AllData0$ID)) {
      
      AllData1 <- subset(AllData0, ID==k)
      AllData1 <- merge(AllData1, Dates, by= c("DT"), all.y=T, sort=T)
      AllData1$ID <- k
      AllData1$Trial <- unique(AllData0$Trial)
      AllData1$Pond <- unique(AllData0$Pond)
      AllData1 <- AllData1[order(AllData1$DT),]
      AllData1 <- subset(AllData1, !is.na(ID))
    
      if (k==unique(AllData0$ID)[1]){
        AllData <- AllData1
    
      }else{
        AllData <- rbind(AllData, AllData1)
      }
    } 
    
    #read in the tags and datetimes to remove
    ToRemv <-  read.csv('~/Carp Model/Supplementary Files/TagsToRemoveForJohn.csv', stringsAsFactors=F)
    ToRemv$DroppedDT <- as.POSIXct(ToRemv$DroppedDT, format="%m/%d/%Y %H:%M:%S", digits=3, tz = "America/Chicago")
    ToRemv <- subset(ToRemv, Trial==TrialNum & Pond==PondNum)
    ToRemv <- ToRemv[,c(1,2,11,12)]

    AllData$indx <- 0
      for (g in 1:nrow(ToRemv)){
       AllData$indx[AllData$ID == ToRemv$TagCode[g] & AllData$DT > ToRemv$DroppedDT[g]] <- 1
      }
    AllData <- subset(AllData, indx==0)
    AllData$indx <- NULL

    save(AllData, file = paste0('ProcessData_',TrialNum,'_pond_',PondNum,'.RDATA'))
  }
}