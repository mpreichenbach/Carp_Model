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


pond.locations <- function(path=file.path(getwd(), "Supplementary Files")){
    # loads GPS coordinates of speakers, kettles, hydrophones in each pond, along
    # with the boundaries of each pond. Returns a list of data frames, one for each
    # type of object.

    LocsDat <- read.csv(file.path(path, "CERC2018Positions.csv"))
    names(LocsDat)[c(1, 2)] <- c("x", "y")
  
    SPK <- subset(LocsDat, Type=="SPK")
    KET <- subset(LocsDat, Type=="KET")
    BND <- subset(LocsDat, Type=="BND")
    HYD <- subset(LocsDat, Type=="HYD")
  
    loc.list = list(speakers=SPK, kettles=KET, boundary=BND, hydrophones=HYD)
    return(loc.list)
}


sound.data <- function(path=file.path(getwd(), "Supplementary Files")){
    # loads the data frame with sound on/off times, pond, and sound type; also processes date/time
    # values to POSIXct format. Returns a dataframe ordered by local times.
    
    SoundDat <- read.csv(file.path(path, "Master_Sound_Tag_20200925.csv"), stringsAsFactors=FALSE)
    
    # remove the following two lines in the future
    TrialNum <- sort(unique(SoundDat$Trial))
    PondNum <- sort(unique(SoundDat$Pond))
    
    colnames(SoundDat)[1] <- 'Trial'
    SoundDat$DT <- paste(SoundDat$DOY, SoundDat$LocalTime..CT.)
    SoundDat$locTimes <- as.POSIXct(SoundDat$DT, format="%m/%d/%Y %H:%M:%S", tz = "America/Chicago")

    SoundDat <- SoundDat %>% relocate(locTimes, .after = 'Pond')
    SoundDat$DT <- NULL
    SoundDat$Sound[SoundDat$Sound=="ON "] <- "ON"
    SoundDat <- SoundDat[order(SoundDat$locTimes),]
    
    return(SoundDat)
}


temperature.data <- function(path=file.path(getwd(), "Supplementary Files")){
    # loads the pond temperature data, converts date/times to POSIXct, and returns a data frame.
    
    TData <- read.csv(file.path(path, "2018CERC_WaterTemperature_AllTrialsPonds.csv"), 
                      stringsAsFactors=F)
    
    TData$DT <- as.POSIXct(TData$DateTime,format="%Y-%m-%d %H:%M:%S")
    TData$DateTime <- NULL
    
    return(TData)
}


treatment.key <- function(trial, pond){
    # this returns the sound treatment type ("Control", "ChirpSquare", "BoatMotor", "ChirpSaw"),
    # given the trial and pond numbers. Note that this function is only valid for the 2018 CERC
    # trials; it will need to change for any future experimental setup.
    
    if (trial == 1) {
        if (pond == 31) {
            treatment <- "ChirpSaw"
        } else if (pond == 30) {
            treatment <- 'BoatMotor'
        } else if (pond == 27) {
            treatment <- 'ChirpSquare'
        } else if (pond == 26) {
            treatment <- "Control"
        }
    }
    if (trial == 2) {
        if (pond == 31) {
            treatment <- "Control"
        } else if (pond == 30) {
            treatment <- "ChirpSaw"
        } else if (pond == 27) {
            treatment <- "BoatMotor"
        } else if (pond == 26) {
            treatment <- 'ChirpSquare'
        }
    }
    if (trial == 3) {
        if (pond == 31) {
            treatment <- 'ChirpSquare'
        } else if (pond == 30) {
            treatment <- "Control"
        } else if (pond == 27) {
            treatment <- "ChirpSaw"
        } else if (pond == 26) {
            treatment <- "BoatMotor"
        }
    }
    if (trial == 4) {
        if (pond == 31) {
            treatment <- "BoatMotor"
        } else if (pond == 30) {
            treatment <- 'ChirpSquare'
        } else if (pond == 27) {
            treatment <- 'Control'
        } else if (pond == 26) {
            treatment <- "ChirpSaw"
        }
    }
    if (trial == 5) {
        if (pond == 31) {
            treatment <- 'ChirpSquare'
        } else if (pond == 30) {
            treatment <- "ChirpSaw"
        } else if (pond == 27) {
            treatment <- "BoatMotor"
        } else if (pond == 26) {
            treatment <- "Control"
        }
    }
    
    return(treatment)
}


fit.crw <- function(trials, ponds, sound_time, seconds_ba, timestep="6 sec", inits=c(2, 0.001), 
                    attempts=20, data_path=file.path(getwd(), "Carp Pond Analysis"),
                    out_path=file.path(getwd(), "Fitted CRWs")){
    # this function loads sound and processed telemetry data, and fits correlated random-walks to
    # the tracks. Instead of returning an object, this function saves .RDATA files to the out_path.

    for (trial in trials){
        for (pond in ponds){
            # s.dat <- subset(sound_data, Sound=="ON" & Trial==trial & Pond==pond)
            # s.dat <- s.dat[order(s.dat$locTimes),]
            treatment <- treatment.key(trial=trial, pond=pond)
            
            # for more flexibility, use RDS instead of RDATA in the processing step; RDATA saves the
            # name of objects too, which must be known in subsequent code.
            load(file.path(data_path, paste0("ProcessData_Trial_", trial, "_pond_", pond, ".RDATA")))
            AllData$Trial <- trial
            AllData$Pond <- pond
            AllData$Treatment <- treatment
            
            #Now subset dataset to just time before and after specified time interval
            crawldat0 <- subset(AllData, !is.na(Easting))
            crawldat0$OnDT <- sound_time
            crawldat0 <- subset(crawldat0, DT>=sound_time-seconds_ba & DT<=sound_time+seconds_ba)

            TagCodes <- unique(crawldat0$ID)
            
            rawdat <- crawldat0[, c("ID","Easting","Northing","DT")]
            colnames(rawdat)<-c('ID','x','y','time')
            print(paste("Data subsetting complete for Trial ", trial, ", Pond ", pond, 
                        ", sound-on time ", as.character(sound_time), "; fitting CRW."))
            
            #Fit the correlated random walk Model
            tempDat0 <- crawlWrap(obsData=rawdat, timeStep=timestep,
                                  theta=inits, fixPar=c(NA, NA), attempts=attempts)
            
            # Add date-time column called 'DT'
            tempDat0$crwPredict$DT <- tempDat0$crwPredict$time
            
            # only keep the necessary covariates
            covDat <- AllData[,c("ID","DT","Trial","Pond","Treatment","Sound")]
            
            # merge the CRW data with covariate info from telemetry
            tempDat0$crwPredict <- merge(tempDat0$crwPredict, covDat, by=c("ID","DT"))
            
            #prepare data for input into movement model and define covariates
            ModDat <- prepData(data=tempDat0, covNames=c("Trial", "Pond", "Treatment", "Sound"))
            
            #output to directory
            sound_time_str <- str_replace_all(as.character(sound_time), ":", "_")
            
            dir.create(file.path(out_path, paste0("Trial ", trial)))
            dir.create(file.path(out_path, paste0("Trial ", trial), sound_time_str))
            file_name <- paste0("CRW_Trial_", trial, "_Pond_", pond, "_SoundOn_", sound_time_str, ".RDATA")
            save(ModDat, file = file.path(out_path, paste0("Trial ", trial), sound_time_str))
            
            print(paste0("Saved fitted random walk for Trial ", trial, 
                         ", Pond ", pond, ", Sound-on time ", sound_time_str, "."))
        }
    }
}


# compile.crw <- function(path="~/Carp Pond Analysis/", sound_data, temperature_data, pond_locations, 
#                         sound_time, trials, ponds, seconds_ba){
#     # this function loads the various saved random-walk outputs and compiles them into one dataframe
#     # per trial.
#     
#     # Read in all correlated random walk outputs files 
#     IDs <- matrix(NA, nrow=length(trials) * length(ponds), ncol=14)
#     count<-0
#     for (trial in trials) {
#         for (pond in ponds) {
#             speaker <- subset(pond_locations$SPK, pond_locations$SPK$Pond==p)
#             
#             # load the fitted CRW file
#             load(paste0("FishDat_Trial_",t,"_Pond_",p,".RDATA"))
#             # get only sound-on times for the particular trial and pond
#             sdat <- subset(sound_data, Sound=="ON" & Trial==t & Pond==p)
#             # order from oldest to newest
#             sdat <- sdat[order(sdat$locTimes),]
#             
#             # gets the relevant temperature data
#             tdat <- subset(TData, Trial==trial & Pond==pond)
#             
#             # order from oldest to newest
#             tdat <- tdat[order(tdat$DT),]
#             
#             # change the time frame below to include more or less data (1800 = 60 secs/min * 30 min)
#             tdat0 <- subset(tdat, DT >= sound_time-seconds_ba & DT <= sound_time+seconds_ba) 
#             
#             # averages the temperatures in a given trial and pond
#             MuT <- ddply(tdat0, .(Trial, Pond),
#                          summarize,
#                          TempC = mean(Temp_C, na.rm=T))
#             
#             # Label loaded model data with correct sound on/off info
#             ModDat$Sound[ModDat$DT<sound_time] <- "0ff"
#             ModDat$Sound[ModDat$DT>=sound_time] <- "On"
#             
#             # compute distance to speaker
#             ModDat$Dist2SPK <- sqrt((speaker$x - ModDat$x)^2 + 
#                                         (speaker$y - ModDat$y)^2)
#             
#             # incorporate temperature
#             ModDat$TempC <- MuT$TempC
#             
#             # just makes a new copy?
#             ModDat0 <- ModDat 
#             
#             # iteratively fill a dataframe with trial/pond tag
#             count<-count+1
#             IDs[count,1] <- trial
#             IDs[count,2] <- pond
#             IDs[count,3:(2+length(unique(ModDat0$ID)))] <- 
#                 as.numeric(as.character(unique(ModDat0$ID)))[1:length(unique(ModDat0$ID))]
#             
#             #just keeps track of codes progress
#             print(paste0("Trial_",t,"Pond_",p,"Nfish_", length(unique(ModDat0$ID))))
#             
#             #Compile all one-hour tracks into one data set for model fitting
#             if (pond==Ponds[1]) { 
#                 FishData0 <- ModDat0
#             }else{
#                 FishData0 <- rbind(FishData0, ModDat0)  
#             }
#         }   
#         if (t==TrialNum[1]) { 
#             FishData <- FishData0
#         }else{
#             FishData <- rbind(FishData, FishData0)  
#         }
#     }
#     
#     return(list(FishData=FishData, IDs=IDs))
# }
# 
# ponds <- c(26, 27, 30, 31)
# for (trial in c(1, 2, 3, 4, 5)){
#     dir.create(paste0("~/Carp-Model/Carp Pond Analysis/Correlated Random Walks/Trial ", trial))
#     time_subset <- sound_data[sound_data$Sound == "ON" & sound_data$Trial == trial, ]
#     times <- sort(unique(round_date(time_subset$locTimes, "hour")))
#     for (i in length(times)){
#         sound_time <- times[i]
#         if (length(hour(sound_time)) == 1){
#             folder_hour <- paste0("0", hour(sound_time), "00")
#         }else{
#             folder_hour <- paste0(hour(sound_time), "00")
#         }
#         folder_name <- paste(date(sound_time), folder_hour)
#         path <- paste0("~/Carp-Model/Carp Pond Analysis/Correlated Random Walks/Trial ", trial, "/",
#                        folder_name, "/")
#         dir.create(path)
#         fit.crw(out_path=path, sound_data=sound_data, trials=c(trial), ponds=ponds, 
#                 sound_time=sound_time, seconds_ba=1800)
#         print(paste0("Finished fitting CRW's for Trial ", trial, " ", folder_name))
#     }
# }