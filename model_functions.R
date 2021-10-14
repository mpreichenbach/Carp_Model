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
    
    SoundDat <- read.csv('C:/Users/RDGRLMPR/Documents/Carp/Master_Sound_Tag_20200925.csv', 
                         stringsAsFactors=FALSE)
    
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
    
}