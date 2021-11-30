library(animation)
library(automap)
library(circular)
library(lubridate)
library(maps)
library(momentuHMM)
library(mvtnorm)
library(parallel)
library(plyr)
library(raster)
library(rgdal)
library(sf)
library(sp)
library(StreamMetabolism)
library(tidyverse)
library(viridis)


compile.crw <- function(on_time, ponds = c(26, 27, 30, 31), path="~/Carp-Model/Fitted CRWs"){
    # compiles various fitted random walks into a single dataset, to be used for HMM fitting.
    
    dt_str <- time.to.str(on_time)
    date_str <- as.character(date(on_time))
    S
    if  (date_str %in% c("2018-06-11", "2018-06-12", "2018-06-13", "2018-06-14")){
        trial <- 1
    }else if (date_str %in% c("2018-06-26", "2018-06-27", "2018-06-28", "2018-06-29")){
        trial <- 2
    }else if (date_str %in% c("2018-07-10", "2018-07-11", "2018-07-12", "2018-07-13")){
        trial <- 3
    }else if (date_str %in% c("2018-07-24", "2018-07-25", "2018-07-26", "2018-07-27")){
        trial <- 4
    }else if (date_str %in% c("2018-08-07", "2018-08-08", "2018-08-09", "2018-08-10")){
        trial <- 5
    }
    
    pond_files <- list.files(path = file.path(path, paste("Trial", trial), dt_str), full.names=TRUE)
    load(file=pond_files[1])
    df <- ModDat
    for (i in 2:length(pond_files)){
        load(pond_files[i])
        df0 <- ModDat
        df <- rbind(df, df0)
    }
    rm(ModDat)
    
    return(df)
}


fit.krig <- function(sound_data, new_data, 
                     crs_string="+proj=utm +zone=15 +ellps=WGS84 +datum=WGS84 +units=m"){
    # this function performs an autoKriging on new_data, and extracts dataframe. Assumes labels of
    # x, y, and dB.
    
    sf_sound <- st_as_sf(sound_data, coords = c("x", "y"), crs = CRS(crs_string))
    sp_new_data <- SpatialPoints(new_data, proj4string = CRS(crs_string))
    
    fit_KRIG <- automap::autoKrige(
        formula = dB ~ 1,
        input_data = as(sf_sound, "Spatial"),
        new_data = sp_new_data
    ) %>%
        .$krige_output %>%
        as.data.frame() %>%
        dplyr::select(x, y, dB = var1.pred)

    return(fit_KRIG)
}


fit.crw <- function(telemetry, trial, pond, on_time, seconds_ba, timestep="6 sec", 
                    inits=c(2, 0.001), retry_fits=100, attempts=100){
    # this function loads sound and processed telemetry data, and fits correlated random-walks to
    # the tracks.
    
    if (hour(on_time) == 0 & minute(on_time) == 0){
        sound_time_str <-paste(as.character(on_time), "00_00_00")
    }else{
        sound_time_str <- str_replace_all(as.character(on_time), ":", "_")
    }
    
    # store sound treatment for the particular trial and pond
    treatment <- treatment.key(trial=trial, pond=pond)
    
    # load processed telemetry data
    telemetry$Trial <- trial
    telemetry$Pond <- pond
    telemetry$Treatment <- treatment
    
    #Now subset dataset to just time before and after specified time interval
    crawldat0 <- subset(telemetry, !is.na(Easting))
    crawldat0$OnDT <- sound_time
    crawldat0 <- subset(crawldat0, DT>=sound_time-seconds_ba & DT<=sound_time+seconds_ba)
    
    rawdat <- crawldat0[, c("ID","Easting","Northing","DT")]
    colnames(rawdat)<-c('ID','x','y','time')
    print(paste0("Data subsetting complete for Trial ", trial, ", Pond ", pond, 
                 ", sound-on time ", sound_time_str, "; fitting CRW."))
    
    #Fit the correlated random walk Model
    tempDat0 <- crawlWrap(obsData=rawdat, timeStep=timestep,
                          theta=inits, fixPar=c(NA, NA), retryFits = retry_fits, attempts=attempts)
    
    # only keep the necessary covariates
    covDat <- telemetry[,c("ID","Time","Trial","Pond","Treatment","Sound", "Diel")]
    
    # merge the CRW data with covariate info from telemetry
    tempDat0$crwPredict <- merge(tempDat0$crwPredict, covDat, by=c("ID","Time"))
    
    #prepare data for input into movement model and define covariates
    ModDat <- prepData(data=tempDat0, covNames=c("Trial", "Pond", "Treatment", "Sound", "Diel"))
    
    #output to directory
    return(ModDat)
}


plot.interp <- function(points, pond, sound, grad_min, grad_max, 
                        colours = c("blue", "red", "yellow")){
    plt <- ggplot(points, aes(x = x, y = y, fill = dB)) +
        geom_raster() +
        ggtitle(label = paste0(sound, ", Pond ", pond)) +
        xlab("Easting") +
        ylab("Northing") +
        scale_fill_gradientn(limits=c(grad_min, grad_max), colours = colours) +
        theme_bw() +
        theme(
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size = 18, hjust = 0.5, margin=margin(t = 15)),
            axis.title.x = element_text(size = 16, margin=margin(b = 15)),
            axis.title.y = element_text(size = 16, margin=margin(l = 15)),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_blank(),
        )
    
    return(plt)
}


pond.locations <- function(path=file.path(getwd(), "Supplementary Files"), bnd_corners_only=TRUE){
    # loads GPS coordinates of speakers, kettles, hydrophones in each pond, along
    # with the boundaries of each pond. Returns a list of data frames, one for each
    # type of object.
    
    # bnd_corners_only=FALSE yields a rounded outline of the pond corners (at least with the
    # CERC 2018 file).

    LocsDat <- read.csv(file.path(path, "CERC2018Positions.csv"), header=TRUE)
    colnames(LocsDat)[c(1, 2)] <- c("x", "y")

    SPK <- LocsDat %>%
        subset(., Type=="SPK") %>%
        select(x, y, Pond)
    KET <- LocsDat %>%
        subset(., Type=="KET") %>%
        select(x, y, Pond)
    HYD <- LocsDat %>%
        subset(., Type=="HYD") %>%
        select(x, y, Pond)
    BND <- LocsDat %>%
        subset(., Type=="BND") %>%
        select(x, y, Pond)

    rownames(SPK) <- rownames(KET) <- rownames(BND) <- rownames(HYD) <- NULL
    
    if (bnd_corners_only){
        ponds = unique(BND$Pond)
        df <- data.frame(matrix(ncol = 3, nrow = 0))
        colnames(df) <- c("x", "y", "Pond")
        for (i in 1:length(ponds)){
            pond = ponds[i]
            bnd_pond <- subset(BND, Pond == pond)
            
            x_limits <- c(min(bnd_pond$x), max(bnd_pond$x))
            y_limits <- c(min(bnd_pond$y), max(bnd_pond$y))
            df_pond <- expand.grid(x_limits, y_limits)
            
            df_pond$Pond <- pond
            colnames(df_pond) <- c("x", "y", "Pond")
            df <- rbind(df, df_pond)
        }
    }
    
    BND <- df
  
    loc.list = list("speakers"=SPK, "kettles"=KET, "hydrophones"=HYD, "boundary"=BND)
    return(loc.list)
}


sound.data <- function(path=file.path(getwd(), "Supplementary Files"), round_str="30 min"){
    # loads the data frame with sound on/off times, pond, and sound type; also processes date/time
    # values to POSIXct format. Rounds values to nearest interval given in round_str; round_str=None
    # keeps unrounded values. Returns a dataframe ordered by local times.
    
    SoundDat <- read.csv(file.path(path, "Master_Sound_Tag_20200925.csv"), stringsAsFactors=FALSE)
    
    colnames(SoundDat)[1] <- 'Trial'
    SoundDat$DT <- paste(SoundDat$DOY, SoundDat$LocalTime..CT.)
    SoundDat$locTimes <- as.POSIXct(SoundDat$DT, format="%m/%d/%Y %H:%M:%S", tz = "America/Chicago")

    SoundDat <- SoundDat %>% relocate(locTimes, .after = 'Pond')
    SoundDat$DT <- NULL
    SoundDat$Sound[SoundDat$Sound=="ON "] <- "ON"
    SoundDat <- SoundDat[order(SoundDat$locTimes),]
    
    if (typeof(round_str) == "character"){
        SoundDat$locTimes <- round_date(SoundDat$locTimes, round_str)
    }
    
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


time.to.str <- function(timeVal, sep = "_", hasDate = TRUE, hasTime = TRUE){
    # turns a POSIXct value into a string
    
    if (!(hasDate & hasTime)) stop("At least one of hasDate or hasTime must be TRUE.")
    
    date_str <- as.character(date(timeVal))
    
    if (length(hour(timeVal)) == 1){
        hr_str <- paste0(0, hour(timeVal))
    }
    if (length(minute(timeVal)) == 1){
        min_str <- paste0(0, minute(timeVal))
    }
    if (length(second(timeVal)) == 1){
        sec_str <- paste0(0, second(timeVal))
    }
    
    time_str <- paste(hr_str, min_str, sec_str, sep=sep)
    
    if (hasDate == FALSE){
        dt_str <- time_str
    }else if (hasTime == FALSE){
        dt_str <- hasDate
    }else{
        dt_str <- paste(date_str, time_str)
    }
    
    return(dt_str)
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

##### this code will add a column of interpolated dB levels to the fitted CRW files
# trials <- c(1, 2, 3, 4, 5)
# for (trial in trials){
#     files0 <- list.files(paste0("~/Carp-Model/Fitted CRWs/Trial ", trial))
#     for (i in 1:length(files0)){
#         file0 <- files0[i]
#         files1 <- list.files(paste0("~/Carp-Model/Fitted CRWs/Trial ", trial, "/", file0))
#         for (j in 1:length(files1)){
#             file1 <- files1[j]
#             print(file1)
#             load(paste0("~/Carp-Model/Fitted CRWs/Trial ", trial, "/", file0, "/", file1))
#             pond <- strtoi(substr(file1, 18, 19))
#             if (unique(ModDat$Treatment) == "Control"){
#                 next
#             }else if (unique(ModDat$Treatment) == "ChirpSaw"){
#                 sound <- "Saw"
#             }else if (unique(ModDat$Treatment) == "ChirpSquare"){
#                 sound <- "Square"
#             }else if (unique(ModDat$Treatment) == "BoatMotor"){
#                 sound <- "BoatMotor"
#             }
#             sound_data <- read.csv(paste0("~/Carp-Model/Supplementary Files/Sound Mapping/UTM, Zone 15/Pond", pond, sound, ".csv"))
#             krig <- fit.krig(sound_data=sound_data, new_data=as.data.frame(ModDat[, c("x", "y")]))
#             ModDat$dB <- krig$dB
#             ModDat$dB[ModDat$Sound == 0] <- 0
#             save(ModDat, file=paste0("~/Carp-Model/Fitted CRWs/Trial ", trial, "/", file0, "/", file1))
#         }
#     }
#}

# trials <- c(1, 2, 3, 4, 5)
# path <- '~/Carp-Model/Fitted CRWs'
# for (trial in trials){
#     subfolders <- list.files(file.path(path, paste0('Trial ', trial)))
#     for (folder in subfolders){
#         files <- list.files(file.path(path, paste0('Trial ', trial), folder))
#         for (file in files){
#             print(file)
#             pond <- substr(file, 18, 19)
#             load(file.path(path, paste0('Trial ', trial), folder, file))
#             ModDat$TimeNum <- NULL
#             # ModDat$Diel <- NA
#             # 
#             # # diel info; 1 is day, 0 is night
#             # SunRS <- sunrise.set(38.9122061924, -92.2795993947,
#             #                      paste0(year(min(ModDat$Time)), '/', month(min(ModDat$Time)), '/',
#             #                             day(min(ModDat$Time))-1), num.days = 8,
#             #                      timezone='America/North_Dakota/Center')
#             # SunRS[,1] <- as.POSIXct(SunRS[,1], origin="1970-01-01", tz = "America/Chicago")
#             # SunRS[,2] <- as.POSIXct(SunRS[,2], origin="1970-01-01", tz = "America/Chicago")
#             # SunRS$Date <- as.Date(SunRS[,1])
#             # 
#             # # Day is 0, night is 1
#             # for(j in 1:(nrow(SunRS)-1)){
#             #     ModDat$Diel[(ModDat$Time >= SunRS$sunrise[j]) & (ModDat$Time < SunRS$sunset[j])] <- 1
#             # }
#             # ModDat$Diel[is.na(ModDat$Diel)] <- 0
#         
# 
#             save(ModDat, file = file.path(path, paste0('Trial ', trial), folder, file))
#         }
#     }
# }
