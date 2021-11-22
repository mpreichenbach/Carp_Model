library(animation)
library(automap)
library(circular)
library(lubridate)
library(maps)
library(momentuHMM)
library(mvtnorm)
library(plyr)
library(raster)
library(rgdal)
library(sf)
library(sp)
library(tidyverse)
library(viridis)


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

plot.interp <- function(points, pond, sound, colours = c("blue", "red", "yellow")){
    plt <- ggplot(points, aes(x = x, y = y, fill = dB)) +
        geom_raster() +
        ggtitle(label = paste0(sound, ", Pond ", pond)) +
        xlab("Easting") +
        ylab("Northing") +
        scale_fill_gradientn(colours = colours) +
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


fit.crw <- function(sound_data, seconds_ba, timestep="6 sec", inits=c(2, 0.001),
                    ncores=detectCores(), retryParallel=TRUE, retryFits=100){
    # this function loads sound and processed telemetry data, and fits correlated random-walks to
    # the tracks. Instead of returning an object, this function saves .RDATA files to the out_path.
    
    for (trial in trials){
        # if no Trial X folder has been created in out_path, then create it
        if (!any(grepl(paste("Trial", trial), list.dirs(out_path, recursive=FALSE)))){
            dir.create(file.path(out_path, paste0("Trial ", trial)))
        }
        
        # get sound times (rounded to nearest 30 mins) for trial
        s.dat <- subset(sound_data, Trial==trial & Sound=="ON")
        s.DT <- unique(round_date(s.dat$locTimes, "30 mins"))
        for (i in 1:length(s.DT)){
            sound_time <- s.DT[i]
            if (hour(sound_time) == 0 & minute(sound_time) == 0){
                sound_time_str <-paste(as.character(sound_time), "00_00_00")
            }else{
                sound_time_str <- str_replace_all(as.character(sound_time), ":", "_")
            }
            
            # if sound time folder hasn't been created, then create it
            time_dirs <- list.dirs(file.path(out_path, paste("Trial", trial)), recursive=FALSE)
            if (!any(grepl(sound_time_str, time_dirs))){
                dir.create(file.path(out_path, paste("Trial", trial), sound_time_str))
            }
            for (pond in ponds){
                # store sound treatment for the particular trial and pond
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
                print(paste0("Data subsetting complete for Trial ", trial, ", Pond ", pond, 
                             ", sound-on time ", sound_time_str, "; fitting CRW."))
                
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
                ModDat <- prepData(data=tempDat0, covNames=c("Trial", "Pond", "Treatment", "Sound"),
                                   spatialCovs=raster_list)
                
                #output to directory
                if (out_path){
                    file_name <- paste0("CRW_Trial_", trial, "_Pond_", pond, "_SoundOn_", sound_time_str, ".RDATA")
                    save(ModDat, file=file.path(out_path, paste("Trial", trial), sound_time_str, file_name))
                    
                    print(paste0("Saved fitted random walk for Trial ", trial, 
                                 ", Pond ", pond, ", Sound-on time ", sound_time_str, ".")) 
                }
                
                return(ModDat)
            }
        }
    }
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
