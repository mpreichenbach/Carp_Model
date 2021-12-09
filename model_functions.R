library(animation)
library(automap)
library(circular)
library(ggplot2)
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


plot.tracks <- function(tel_data = NULL, crw_data = NULL, id = NULL, min_time, max_time, tel_col = "dodgerblue1", 
                        crw_col = "orange1", full_extent=FALSE, background = c("satellite", "sound")){
    # plots the telemetry tracks and the fitted CRW for the time period between min_time and
    # max_time. Returns a list of plots, one for each ID. Assumes that time values are POSIXct, and 
    # are in the columns tel_data$DT and crw_data$Time
    
    # time order check
    if (min_time > max_time){stop("Minimum time is greater than maximum time.")}
    
    # subset data
    if (!is.null(tel_data)){
        tel_data <- tel_data %>% 
            drop_na() %>%
            dplyr::select(DT, Easting, Northing, ID, Pond, Trial) %>%
            dplyr::filter((min_time < DT) & (DT < max_time))
        if (nrow(tel_data) == 0){stop("No telemetry data in specified time interval.")}
    }
    if (!is.null(crw_data)){
        crw_data <- crw_data %>%
            drop_na() %>%
            dplyr::select(ID, Time, x, y, Trial, Pond) %>%
            dplyr::filter((min_time < Time) & (Time < max_time))
        if (nrow(crw_data) == 0){stop("No CRW data in the specified time interval.")}
    }

    # make sure dataframes are sorted by times
    tel_data <- tel_data[order(tel_data$DT),]
    crw_data <- crw_data[order(crw_data$Time),]
    
    # ID checks
    tel_ids <- unique(tel_data$ID)
    crw_ids <- unique(crw_data$ID)
    if (is.null(id)){
        ids <- union(tel_ids, crw_ids)
        if (length(ids) != length(tel_ids)){stop("There are more CRW IDs than telemetry IDs.")}
        if (length(ids) != length(crw_ids)){stop("There are more telemetry IDs than CRW IDs.")}
        if (!(any(tel_ids %in% crw_ids))){stop("Telemetry IDs are not a subset of CRW IDs.")}
        if (!(any(crw_ids %in% tel_ids))){stop("CRW IDs are not a subset of telemetry IDs.")}
    }else{
        ids <- id
        if (!(any(ids %in% tel_ids))){stop("Some IDs supplied are not in the telemetry IDs.")}
        if (!(any(ids %in% crw_ids))){stop("Some IDs supplied are not in the CRW IDs.")}
    }
    
    # pond checks
    tel_pond <- unique(tel_data$Pond)
    crw_pond <- unique(crw_data$Pond)
    if (length(tel_pond) < 1){
        stop("No pond values in the telemetry data.")
    }else if (length(tel_pond) > 1){
        stop("More than one pond value in the telemetry data.")
    }
    
    if (length(crw_pond) < 1){
        stop("No pond values in the CRW data.")
    }else if (length(crw_pond) > 1){
        stop("More than one pond value in the CRW data.")
    }
    
    if (tel_pond != crw_pond){
        stop("Telemetry pond does not match the CRW pond.")
    }else{
        pond <- tel_pond
    }
    
    # trial checks
    tel_trial <- unique(tel_data$Trial)
    crw_trial <- unique(crw_data$Trial)
    if (length(tel_trial) < 1){
        stop("No trial values in the telemetry data.")
    }else if (length(tel_trial) > 1){
        stop("More than one trial value in the telemetry data.")
    }
    
    if (length(crw_trial) < 1){
        stop("No trial values in the CRW data.")
    }else if (length(crw_trial) > 1){
        stop("More than one trial value in the CRW data.")
    }
    
    if (tel_trial != crw_trial){
        stop("Telemetry pond does not match the CRW pond.")
    }else{
        trial <- tel_trial
    }
    
    # converts IDs to strings for later naming
    ids_str <- c()
    for (i in 1:length(ids)){
        ids_str <- append(ids_str, toString(ids[i]))
    }
    
    # get pond GPS coordinates
    bnd <- subset(pond.locations()$boundary, Pond == pond)
    min_x <- min(bnd$x)
    max_x <- max(bnd$x)
    min_y <- min(bnd$y)
    max_y <- max(bnd$y)
    
    if (background == satellite){
        sat <- read
    }
    
    
    # get time strings for plotting
    min_time_str <- time.to.str(min_time, sep=":")
    if (date(min_time) != date(max_time)){
        max_time_str <- time.to.str(max_time, sep=":")
    }else{
        max_time_str <- substr(time.to.str(max_time, sep=":"), 12, 19)
    }
    
    # create list of plots
    plot_list <- list()
    for (i in 1:length(ids)){
        tag <- ids[i]
        tel_sub <- subset(tel_data, ID == tag)
        crw_sub <- subset(crw_data, ID == tag)
        
        plt <- list(ggplot() + 
                        geom_point(data = tel_sub, mapping = aes(x = Easting, y = Northing),
                                   colour = tel_col, size=2) + 
                        geom_point(data = crw_sub, mapping = aes(x = x, y = y), colour = crw_col, 
                                   size = 2) +
                        geom_path(data = tel_sub, mapping = aes(x = Easting, y = Northing,
                                                                colour = "Telemetry"), size = 1) +
                        geom_path(data = crw_sub, mapping = aes(x = x, y = y, colour = "CRW"), 
                                  size = 1) +
                        scale_colour_manual("", breaks = c("Telemetry", "CRW"),
                                            values = c(tel_col, crw_col)) +
                        ggtitle(label = paste0("Trial ", trial, ", Pond ", pond, ", ID = ", tag, 
                                               ",\n", min_time_str, " to ", max_time_str)) +
                        {if (full_extent)xlim(min_x, max_x)} + 
                        {if (full_extent)ylim(min_y, max_y)} + 
                        xlab("Easting") +
                        ylab("Northing") +
                        theme_bw() + 
                        theme(
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            plot.title = element_text(size = 18, hjust = 0.5, 
                                                      margin=margin(t = 15)),
                            axis.title.x = element_text(size = 16, margin=margin(b = 0)),
                            axis.title.y = element_text(size = 16, margin=margin(l = 0)),
                            legend.position = c(0.9, 0.9),
                            legend.text = element_text(size=16)
                        )
        )
        
        plot_list <- append(plot_list, plt)
    }
    
    # make ID string the name of each respective plot
    names(plot_list) <- ids_str
    
    return(plot_list)
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
        dplyr::select(x, y, Pond)
    KET <- LocsDat %>%
        subset(., Type=="KET") %>%
        dplyr::select(x, y, Pond)
    HYD <- LocsDat %>%
        subset(., Type=="HYD") %>%
        dplyr::select(x, y, Pond)
    BND <- LocsDat %>%
        subset(., Type=="BND") %>%
        dplyr::select(x, y, Pond)

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
    
    if (nchar(hour(timeVal)) == 1){
        hr_str <- paste0(0, hour(timeVal))
    }else{
        hr_str <- as.character(hour(timeVal))
    }
    if (nchar(minute(timeVal)) == 1){
        min_str <- paste0(0, minute(timeVal))
    }else{
        min_str <- as.character(minute(timeVal))
    }
    if (nchar(second(timeVal)) == 1){
        sec_str <- paste0(0, second(timeVal))
    }else{
        sec_str <- as.character(minute(timeVal))
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

##### this code generates maps centered in the middle of the ponds
# for (pond in c(26, 27, 30, 31)){
#     bnd <- subset(pond_locations$boundary, Pond == pond)
#     center <- c((min(bnd$x) + max(bnd$x)) / 2, (min(bnd$y) + max(bnd$y)) / 2)
#     center_utm <- as.data.frame(t(center))
#     colnames(center_utm) <- c("x", "y")
#     coordinates(center_utm) <- ~x + y
#     proj4string(center_utm) <- CRS("+proj=utm +zone=15 +datum=WGS84 +units=m +ellps=WGS84")
#     center_ll <- spTransform(center_utm, CRS("+proj=longlat +datum=WGS84"))
# 
#     long <- center_ll$x
#     lat <- center_ll$y
# 
#     name <- paste("Pond_", pond)
#     plt <- get_googlemap(center = c(long, lat), zoom = 20, scale = 2, maptype = "satellite")
#     saveRDS(plt, file.path("~/Carp-Model/Supplementary Files/Pond Maps", paste0(name, ".RDS")))
# }

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
