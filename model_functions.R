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


pond.locations <- function(path="~/Supplementary Files/"){
    # loads GPS coordinates of speakers, kettles, hydrophones in each pond, along
    # with the boundaries of each pond. Returns a list of dataframes, one for each
    # type of object.

    LocsDat <- read.csv(paste0(path, "CERC2018Positions.csv"))
    names(LocsDat)[c(1, 2)] <- c("x", "y")
  
    SPK <- subset(LocsDat, Type=="SPK")
    KET <- subset(LocsDat, Type=="KET")
    BND <- subset(LocsDat, Type=="BND")
    HYD <- subset(LocsDat, Type=="HYD")
  
    loc.list = list(speakers=SPK, kettles=KET, boundary=BND, hydrophones=HYD)
    return(loc.list)
}


sound.data <- function(path="~/Supplementary Files/"){
    # loads the dataframe with sound on/off times, pond, and sound type; also processes date/time
    # values to POSIXct format. Returns a dataframe ordered by local times.
    
    SoundDat <- read.csv(paste0(path, "Master_Sound_Tag_20200925.csv", stringsAsFactors=FALSE))
    
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
 

temperature.data <- function(path="~/Supplementary Files/"){
    # loads the pond temperature data, converts date/times to POSIXct, and returns a dataframe.
    
    TData <- read.csv(paste0(path, "2018CERC_WaterTemperature_AllTrialsPonds.csv"), 
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


fit.crw <- function(data_path="~/Carp Pond Analysis/", out_path, sound_data, trials, ponds,
                    sound_time, seconds_ba, timestep="6 sec", inits=c(2, 0.001, attempts=20)){
    # this function loads sound and processed telemetry data, and fits correlated random-walks to
    # the tracks. Instead of returning an object, this function saves .RDATA files to the out_path.
    
    for (trial in trials){
        for (pond in ponds){
            s.dat <- subset(sound_data, Sound=="ON" & Trial==trial & Pond==pond)
            s.dat <- s.dat[order(s.dat$locTimes),]
            treatment <- treatment.key(trial=trial, pond=pond)
            
            # for more flexibility, use RDS instead of RDATA in the processing step; RDATA saves the
            # name of objects too, which must be known in subsequent code.
            load(paste0('~/Carp Pond Analysis/ProcessData_', trial, '_pond_', pond, '.RDATA'))
            AllData$Treatment <- treatment
            
            #Now subset dataset to just time before and after specified time interval
            crawldat0 <- subset(AllData, !is.na(Easting))
            crawldat0$OnDT <- as.POSIXct(Sdat$locTimes[length(Sdat$locTimes)], 
                                         format="%Y-%m-%d %H:%M:%S")
            crawldat0 <- subset(crawldat0, DT >= OnDT-1800 & DT <= OnDT+1800)

            # do we need this line?
            crawldat0 <- crawldat0[,c("ID","Easting","Northing","DT","Trial","Pond","Treatment",
                                      "Sound")]
            
            TagCodes <- unique(crawldat0$ID)
            
            rawdat <- crawldat0[, c("ID","Easting","Northing","DT")]
            colnames(rawdat)<-c('ID','x','y','time')
            
            #Fit the correlated random walk Model
            tempDat0 <- crawlWrap(obsData=rawdat, timeStep=timestep,
                                  theta=inits, fixPar=c(NA, NA), attempts=20)
            
            # Add date-time column called 'DT'
            tempDat0$crwPredict$DT <- tempDat0$crwPredict$time
            
            # only keep the necessary covariates
            covDat <- AllData[,c("ID","DT","Trial","Pond","Treatment","Sound")]
            
            # merge the CRW data with covariate infor from telemetry
            tempDat0$crwPredict <- merge(tempDat0$crwPredict, covDat, by=c("ID","DT"))
            
            #prepare data for input into movement model and define covariates
            ModDat <- prepData(data=tempDat0, covNames=c("Trial", "Pond", "Treatment", "Sound"))
            
            #output to directory
            save(ModDat, file = paste0("FishDat_Trial_", t,"_Pond_", p,".RDATA"))
            
            print(paste0("Saved fitted random walk for Trial ", trial, " Pond ", pond, "."))
        }
    }
}






















