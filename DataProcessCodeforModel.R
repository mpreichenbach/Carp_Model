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
data.path <- 'C:/Users/RDGRLMPR/Documents/Carp/Aggregated Data'
usgs.path <- 'C:/Users/RDGRLMPR/Documents/Carp_Model/USGS'

#Set working directory
setwd(paste0(usgs.path, '/Carp Pond Analysis'))
#Please specify the trial, pond number, and scenario
TrialNums <- c(1, 2, 3, 4, 5)
PondNums <- c(26, 27, 30, 31)

#-------------------------------------------------------------------------------------
rm(realDat)
for (TrialNum in TrialNums) {
  for (PondNum in PondNums) {
    print(paste0('Working on Trial ', TrialNum, ' and Pond ', PondNum, '.'))
    #assign file names by Trial
    
     # TrialNum <- 3
     # PondNum <- 27
    
    if (TrialNum == 1) {
      if (PondNum == 31) {
        filename <- paste0(data.path, '/Trial1/TagHistoryFor_AllTags-Trial1-System100_Pond31_ChirpSaw.csv')
      } else if (PondNum == 30) {
        filename <- paste0(data.path, '/Trial1/TagHistoryFor_AllTags-Trial1-System200_Pond30_BoatMotor.csv')
      } else if (PondNum == 27) {
        filename <- paste0(data.path, '/Trial1/TagHistoryFor_AllTags-Trial1-System300_Pond27_ChirpSquare.csv')
      } else if (PondNum == 26) {
        filename <- paste0(data.path, '/Trial1/TagHistoryFor_AllTags-Trial1-System400_Pond26_Control.csv')
      }
    }
    if (TrialNum == 2) {
      if (PondNum == 31) {
        filename <- paste0(data.path, '/Trial2/TagHistoryFor_AllTags-Trial2_Pond31_Control_DroppedTagsRemoved.csv')
      } else if (PondNum == 30) {
        filename <- paste0(data.path, '/Trial2/TagHistoryFor_AllTags-Trial2_Pond30_ChirpSaw_DroppedTagsRemoved.csv')
      } else if (PondNum == 27) {
        filename <- paste0(data.path, '/Trial2/TagHistoryFor_AllTags-Trial2_Pond27_BoatMotor_DroppedTagsRemoved.csv')
      } else if (PondNum == 26) {
        filename <- paste0(data.path, '/Trial2/TagHistoryFor_AllTags-Trial2_Pond26_ChirpSquare_DroppedTagsRemoved.csv')
      }
    }
    if (TrialNum == 3) {
      if (PondNum == 31) {
        filename <- paste0(data.path, '/Trial3/TagHistoryFor_AllTags-Trial3_Pond31_ChirpSquare_DroppedTagsRemoved.csv')
      } else if (PondNum == 30) {
        filename <- paste0(data.path, '/Trial3/TagHistoryFor_AllTags-Trial3_Pond30_Control_DroppedTagsRemoved.csv')
      } else if (PondNum == 27) {
        filename <- paste0(data.path, '/Trial3/TagHistoryFor_AllTags-Trial3_Pond27_ChirpSaw_DroppedTagsRemoved.csv')
      } else if (PondNum == 26) {
        filename <- paste0(data.path, '/Trial3/TagHistoryFor_AllTags-Trial3_Pond26_BoatMotor_DroppedTagsRemoved.csv')
      }
    }
    if (TrialNum == 4) {
      if (PondNum == 31) {
        filename <- paste0(data.path, '/Trial4/TagHistoryFor_AllTags-Trial4_Pond31_BoatMotor_DroppedTagsRemoved.csv')
      } else if (PondNum == 30) {
        filename <- paste0(data.path, '/Trial4/TagHistoryFor_AllTags-Trial4_Pond30_ChirpSquare_DroppedTagsRemoved.csv')
      } else if (PondNum == 27) {
        filename <- paste0(data.path, '/Trial4/TagHistoryFor_AllTags-Trial4_Pond27_Control_DroppedTagsRemoved.csv')
      } else if (PondNum == 26) {
        filename <- paste0(data.path, '/Trial4/TagHistoryFor_AllTags-Trial4_Pond26_ChirpSaw_DroppedTagsRemoved.csv')
      }
    }
    if (TrialNum == 5) {
      if (PondNum == 31) {
        filename <- paste0(data.path, '/Trial5/TagHistoryFor_AllTags-Trial5_Pond31_ChirpSquare_DroppedTagsRemoved.csv')
      } else if (PondNum == 30) {
        filename <- paste0(data.path, '/Trial5/TagHistoryFor_AllTags-Trial5_Pond30_ChirpSaw_DroppedTagsRemoved.csv')
      } else if (PondNum == 27) {
        filename <- paste0(data.path, '/Trial5/TagHistoryFor_AllTags-Trial5_Pond27_BoatMotor_DroppedTagsRemoved.csv')
      } else if (PondNum == 26) {
        filename <- paste0(data.path, '/Trial5/TagHistoryFor_AllTags-Trial5_Pond26_Control_DroppedTagsRemoved.csv')
      }
    }
    #-------------------------------------------------------------------------------------
    #read in the detection data
    realDat <- read.csv(filename, stringsAsFactors=F)
    realDat$X <- NULL #Beware some data sets have an extra column

    TagCodes <- unique(realDat$TagCode)
    
    realDat<-realDat[realDat$Easting!=0,]
    
    realDat <- realDat[order(realDat$TagCode, realDat$LocalTime),]
    
    #-------------------------------------------------------------------------------------
    #Read in the Pond and attribute location data
    # LocsDat <- read.csv(paste0(getwd(),"/CERC2018Positions.csv")) #AEH_18_CERCSOUND_01_GPS
    # LocsDat <- subset(LocsDat, Pond==PondNum)
    # names(LocsDat)[c(1, 2)] <- c("x", "y")
    # 
    # SPK <- subset(LocsDat, Type=='SPK')
    # KET <- subset(LocsDat, Type=='KET')
    # BND <- subset(LocsDat, Type=='BND')
    # HYD <- subset(LocsDat, Type=='HYD')
    
    
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
      
      FishDT[i,1] <- minDT #minimum date time for the fish
      FishDT[i,2] <- maxDT #maximum date time for the fish
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
    #Read in the Sound on off time data hard wired to trial 1 folder
    SoundDat <- read.csv(paste0(data.path, '/Master_Sound_Tag_20200925.csv'), stringsAsFactors=FALSE)
    # Matt added the following two lines 8/3/2021
    colnames(SoundDat)[1] <- 'Trial'
    SoundDat$DT <- paste(SoundDat$DOY, SoundDat$LocalTime..CT.)
    SoundDat$locTimes <- as.POSIXct(SoundDat$DT, format="%m/%d/%Y %H:%M:%S", tz = "America/Chicago")
    # Matt added the following line on 8/3/2021
    SoundDat <- SoundDat %>% relocate(locTimes, .after = 'Pond')
    SoundDat$DT <- NULL
    SoundDat$Sound[SoundDat$Sound=="ON "] <- "ON"
    SoundDat <- SoundDat[order(SoundDat$locTimes),]
    SoundDat <- subset(SoundDat, Trial==TrialNum & Pond==PondNum)
    
    Offdat <- subset(SoundDat, Sound=="OFF")
    Offdat$Sound  <- NULL
    colnames(Offdat) <- c("Trial", "Pond", "OffDT")
    # Matt added the following line on 8/3/2021
    Offdat <- Offdat[, 1:3]
    
    Ondat <- subset(SoundDat, Sound=="ON")
    Ondat$Sound  <- NULL
    colnames(Ondat) <- c("Trial", "Pond", "OnDT")
    # Matt added the following line on 8/3/2021
    Ondat <- Ondat[, 1:3]
    
    Sdat <- cbind(Ondat, Offdat[,-c(1:2)])
    colnames(Sdat) <- c("Trial", "Pond", "OnDT", "OffDT")
    
    stemp <- Sdat[1,]
    stemp$OffDT <- as.POSIXct("06/01/2018 00:00:00",format="%m/%d/%Y %H:%M:%S", tz = "America/Chicago")
    stemp$OnDT <- as.POSIXct("06/01/2018 00:00:00",format="%m/%d/%Y %H:%M:%S", tz = "America/Chicago")
    
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
    
    # table(Dates$Sound, Dates$Sindex, useNA = "ifany")
    # tail(Dates)
    #-------------------------------------------------------------------------------------
    #Obtain sunrise times
    #Here's the function - sunrise.set(lat, long, date, timezone = "UTC", num.days = 1)
    #using the location of Pond 26
    SunRS <- sunrise.set(38.9122061924, -92.2795993947, paste0(year(min(AllData0$DT)),'/',month(min(AllData0$DT)),'/',day(min(AllData0$DT))-1), 
                         num.days = 8, timezone='America/North_Dakota/Center')
    SunRS[,1] <- as.POSIXct(SunRS[,1], origin="1970-01-01", tz = "America/Chicago")
    SunRS[,2] <- as.POSIXct(SunRS[,2], origin="1970-01-01", tz = "America/Chicago")
    SunRS$Date <- as.Date(SunRS[,1])
    
    for(j in 1:(nrow(SunRS)-1)){
      Dates$Diel[(Dates$DT >= SunRS$sunrise[j]) & (Dates$DT < SunRS$sunset[j])] <- 0
    }
    Dates$Diel[is.na(Dates$Diel)] <- 1
    
    #and merge the environmental conditons
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
    ToRemv <-  read.csv('D:/Carp Pond Analysis/FTP/withDropTags/TagsToRemoveForJohn.csv', stringsAsFactors=F)
    ToRemv$DroppedDT <- as.POSIXct(ToRemv$DroppedDT, format="%m/%d/%Y %H:%M:%S", digits=3, tz = "America/Chicago")
    ToRemv <- subset(ToRemv, Trial==TrialNum & Pond==PondNum)
    ToRemv <- ToRemv[,c(1,2,11,12)]

    AllData$indx <- 0
      for (g in 1:nrow(ToRemv)){
       AllData$indx[AllData$ID == ToRemv$TagCode[g] & AllData$DT > ToRemv$DroppedDT[g]] <- 1
      }
    AllData <- subset(AllData, indx==0)
    AllData$indx <- NULL

    #-------------------------------------------------------------------------------------
    ###################### Export the processed data from above #########################
    #-------------------------------------------------------------------------------------
    # save the output data to the directory
    save(AllData, file = paste0('ProcessData_',TrialNum,'_pond_',PondNum,'.RDATA'))
    #####################################################################################
    #subset to the time interval currently used
    # if (TrialNum==1) UseDat <- subset(AllData, DT >= as.POSIXct("2018-06-10 06:00:00") & DT <= as.POSIXct("2018-06-13 06:00:00"))
    # if (TrialNum==2) UseDat <- subset(AllData, DT >= as.POSIXct("2018-06-25 06:00:00") & DT <= as.POSIXct("2018-06-28 06:00:00"))
    # if (TrialNum==3) UseDat <- subset(AllData, DT >= as.POSIXct("2018-07-09 06:00:00") & DT <= as.POSIXct("2018-07-12 06:00:00"))
    # if (TrialNum==4) UseDat <- subset(AllData, DT >= as.POSIXct("2018-07-23 06:00:00") & DT <= as.POSIXct("2018-07-26 06:00:00"))
    # if (TrialNum==5) UseDat <- subset(AllData, DT >= as.POSIXct("2018-08-06 06:00:00") & DT <= as.POSIXct("2018-08-09 06:00:00"))
    # 
    # runStats <- matrix(NA, ncol=3, nrow=1)
    # runStats[,1] = TrialNum
    # runStats[,2]= PondNum
    # runStats[,3] = length(unique(UseDat$ID))
    # colnames(runStats) <- c("Trial", "Pond", "NFish")
    #   
    # print(paste0("Used fish Trial ", TrialNum, "; Pond ", PondNum, "; Nfish ", length(unique(UseDat$ID))))
    # print(paste0("Entire data set Trial ", TrialNum, "; Pond ", PondNum, "; Nfish ", length(unique(AllData$ID))))
    # 
    #  rm(realDat)
    #  rm(AllData)
    #  rm(AllData0)
    #  rm(AllData1)
    #  rm(UseDat)
    #  
    #Compile fish stats
      # if (TrialNum==TrialNums[1] & PondNum==PondNums[1]) {
      #   StatRun <- runStats
      # }else{
      #   StatRun <- rbind(StatRun, runStats)   
      # }
  }
}

#write.table(StatRun, paste0(getwd(),"/IncludedTags.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA")

