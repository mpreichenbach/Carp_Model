library(automap)
library(lubridate)
library(sf)
library(sp)
library(tidyverse)


convert_coords <- function(df, 
                           input_names = c("Long", "Lat"),
                           input_crs = "+proj=longlat +datum=WGS84",
                           output_crs = "+proj=utm +zone=15 +datum=WGS84 +units=m +ellps=WGS84") {
    
    # takes a dataframe of x,y (Easting, Northing) points and converts them to output projection.
    
    # ensure columns are numeric
    for (name in input_names) {df[[name]] <- as.numeric(df[[name]])}
    
    # convert to new coordinate system
    df <- df[, input_names]
    colnames(df) <- c("x", "y")
    
    coordinates(df) <- ~x + y
    proj4string(df) <- CRS(input_crs)
    
    converted <- as.data.frame(spTransform(df, CRSobj=CRS(output_crs)))

    converted
}

fit_krig <- function(sound_data, pred_data, 
                     crs_string="+proj=utm +zone=15 +ellps=WGS84 +datum=WGS84 +units=m"){
    # this function performs an autoKriging on new_data. Assumes column labels are x, y, and dB.

    sf_sound <- st_as_sf(sound_data, coords = c("x", "y"), crs = CRS(crs_string))
    sp_pred_data <- SpatialPoints(as.data.frame(pred_data), proj4string = CRS(crs_string))
    
    krig <- automap::autoKrige(
        formula = dB ~ 1,
        input_data = as(sf_sound, "Spatial"),
        new_data = sp_pred_data
    ) %>%
        .$krige_output %>%
        as.data.frame() %>%
        dplyr::select(x, y, dB = var1.pred)
    
    krig
}


correct_tags <- function(trial, pond){
    # gives a numeric vector with the correct tag codes for each trial


    if (trial == 2 & pond == 26){
        return(c(2047.365, 2064.365, 2081.365, 2098.365, 2115.365, 2132.365, 2149.365, 2166.365,
                 2183.365, 2200.365))
    }else if (trial == 2 & pond == 27){
        return(c(2047.348, 2064.348, 2081.348, 2098.348, 2115.348, 2132.348, 2149.348, 2166.348,
               2183.348, 2200.348))
    }else if (trial == 2 & pond == 30){
        return(c(2047.330, 2064.330, 2081.330, 2098.330, 2115.330, 2132.330, 2149.330, 2166.330,
                 2183.330, 2200.330))
    }else if (trial == 2 & pond == 31){
        return(c(2047.311, 2064.311, 2081.311, 2098.311, 2115.311, 2132.311, 2149.311, 2166.311,
                 2183.311, 2200.311))
    }else if (trial == 3 & pond == 26){
        return(c(2047.423, 2064.423, 2081.423, 2098.423, 2115.423, 2132.423, 2149.423, 2166.423,
                 2183.423, 2200.423))
    }else if (trial == 3 & pond == 27){
        return(c(2047.41, 2064.41, 2081.41, 2098.41, 2115.41, 2132.41, 2149.41, 2166.41, 2183.41,
                 2200.41))
    }else if (trial == 3 & pond == 30){
        return(c(2047.396, 2064.396, 2081.396, 2098.396, 2115.396, 2132.396, 2149.396, 2166.396,
                 2183.396, 2200.396))
    }else if (trial == 3 & pond == 31){
        return(c(2047.381, 2064.381, 2081.381, 2098.381, 2115.381, 2132.381, 2149.381, 2166.381,
                 2183.381, 2200.381))
    }else if (trial == 4 & pond == 26){
        return(c(2047.512, 2064.512, 2081.512, 2098.512, 2115.512, 2132.512, 2149.512, 2166.512,
                 2183.512, 2200.512, 2217.512, 2234.512, 2251.512))
    }else if (trial == 4 & pond == 27){
        return(c(2047.456, 2064.456, 2081.456, 2098.456, 2115.456, 2132.456, 2149.456, 2166.456,
                 2183.456, 2200.456, 2217.456, 2234.456, 2251.456))
    }else if (trial == 4 & pond == 30){
        return(c(2047.446, 2064.446, 2081.446, 2098.446, 2115.446, 2132.446, 2149.446, 2166.446,
                 2183.446, 2200.446, 2217.446))
    }else if (trial == 4 & pond == 31){
        return(c(2047.435, 2064.435, 2081.435, 2098.435, 2115.435, 2132.435, 2149.435, 2166.435,
                 2183.435, 2200.435, 2217.435, 2234.435))
    }else if (trial == 5 & pond == 26){
        return(c(2047.562, 2064.562, 2081.562, 2098.562, 2115.562, 2132.562, 2149.562, 2166.562,
                 2183.562, 2200.562))
    }else if (trial == 5 & pond == 27){
        return(c(2047.548, 2064.548, 2081.548, 2098.548, 2115.548, 2132.548, 2149.548, 2166.548,
                 2183.548, 2200.548))
    }else if (trial == 5 & pond == 30){
        return(c(2047.535, 2064.535, 2081.535, 2098.535, 2115.535, 2132.535, 2149.535, 2166.535,
                 2183.535, 2200.535))
    }else if (trial == 5 & pond == 31){
        return(c(2047.523, 2064.523, 2081.523, 2098.523, 2115.523, 2132.523, 2149.523, 2166.523,
                 2183.523, 2200.523))
    }
}


get_formulas <- function(nCov, 
                         must_have_covs=NULL,
                         ignore_covs=NULL) {
    # outputs a list of formulas with the specified number of covariates; here we use the covariates
    # "Trial", "Temperature", "dB", "Treatment". Additionally, we allow the option of including
    # "Diel" (aka, day/night) as a covariate. The "must_have_covs" and "ignore_covs" arguments, if 
    # they're lists longer than 1, are assumed to mean formulas must have all entries.
    
    if (!(nCov %in% c(0, 1, 2, 3, 4, 5))) {
        stop("nCov must be 0, 1, 2, 3, 4, 5")
    }
    
    # yield a character vector of formulas with nCov covariates, including "Diel"
    if (nCov == 0) {
        frm_names <- c("~1")
    } else if (nCov == 1) {
        frm_names <- c("~Trial", "~Temperature", "~Diel", "~dB", "~Treatment")
    } else if (nCov == 2) {
        frm_names <- c("~Trial+Temperature", "~Trial+Diel", "~Trial+dB", "~Trial+Treatment",
                   "~Temperature+Diel", "~Temperature+dB", "~Temperature+Treatment", "~Diel+dB", 
                   "~Diel+Treatment", "~dB+Treatment")
    } else if (nCov == 3) {
        frm_names <- c("~Trial+Temperature+Diel", "~Trial+Temperature+dB", 
                   "~Trial+Temperature+Treatment", "~Trial+Diel+dB", "~Trial+Diel+Treatment",
                   "~Trial+dB+Treatment", "~Temperature+Diel+dB", "~Temperature+Diel+Treatment", 
                   "~Temperature+dB+Treatment", "~Diel+dB+Treatment")
    } else if (nCov == 4) {
        frm_names <- c("~Trial+Temperature+Diel+dB", "~Trial+Temperature+Diel+Treatment",
                   "~Trial+Temperature+dB+Treatment", "~Trial+Diel+dB+Treatment", 
                   "~Temperature+Diel+dB+Treatment")
    } else if (nCov == 5) {
        frm_names <- c("~Trial+Temperature+Diel+dB+Treatment")
    }
    
    # sequentially remove formulas which do not have the must_have_covs entries
    for (cvrt in must_have_covs) {
        frm_names <- frm_names[grep(cvrt, frm_names)]
    }
    
    # sequentially remove covariates which have the ignore_covs entries
    for (cvrt in ignore_covs) {
        frm_names <- frm_names[!grepl(cvrt, frm_names)]
    }
    
    # define a list to hold the formulas as formulas
    frm_list <- list()
    
    if (length(frm_names) == 0) {
        return(frm_list)
    } else {
        for (i in 1:length(frm_names)) {
            frm_list[[frm_names[i]]] <- formula(paste(frm_names[i], collapse=""))
        }
    }
    
    frm_list
}

pond_locations <- function(path=file.path(getwd(), "Supplementary Files"), 
                           bnd_corners_only=TRUE){
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
  
    loc_list = list("speakers"=SPK, "kettles"=KET, "hydrophones"=HYD, "boundary"=BND)
    loc_list
}


sound_data <- function(path=file.path(getwd(), "Supplementary Files"), 
                       round_str="30 min"){
    # loads the data frame with sound on/off times, pond, and sound type; also processes date/time
    # values to POSIXct format. Rounds values to nearest interval given in round_str; round_str=None
    # keeps unrounded values. Returns a dataframe ordered by local times.
    
    SoundDat <- read.csv(file.path(path, "Master_Sound_Tag_20200925.csv"), stringsAsFactors=FALSE)
    SoundDat[SoundDat$Sound=="ON ", c("Sound")] <- "ON"
    SoundDat <- SoundDat[SoundDat$Sound == "ON",]
    
    SoundDat$Time <- paste(SoundDat$DOY, SoundDat$LocalTime..CT.)
    SoundDat$Time <- as.POSIXct(SoundDat$Time, format="%m/%d/%Y %H:%M:%S", tz = "America/Chicago")
    if (typeof(round_str) == "character"){
        SoundDat$Time <- round_date(SoundDat$Time, round_str)
    }

    unique_times <- unique(SoundDat$Time)
    df <- data.frame(matrix(nrow=length(unique_times), ncol=2, 
                            dimnames=list(1:length(unique_times), c("Time", "Repetition"))))
    

    df$Time <- unique_times[order(unique_times)]
    df$Repetition <- (((as.numeric(row.names(df)) - 1) %% 24) + 1)

    df    
}



time_to_str <- function(timeVal, sep = "_", hasDate = TRUE, hasTime = TRUE){
    # turns a POSIXct value into a string
    
    if (!(hasDate | hasTime)) stop("At least one of hasDate or hasTime must be TRUE.")
    
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
    
    dt_str
}

treatment_key <- function(trial, pond){
    # this returns the sound treatment type ("Control", "ChirpSquare", "BoatMotor", "ChirpSaw"),
    # given the trial and pond numbers. Note that this function is only valid for the 2018 CERC
    # trials; it will need to change for any future experimental setup.
    
    if (trial == 1) {
        if (pond == 31) {
            treatment <- "Saw"
        } else if (pond == 30) {
            treatment <- '100Hp'
        } else if (pond == 27) {
            treatment <- 'Square'
        } else if (pond == 26) {
            treatment <- "Control"
        }
    } 
    if (trial == 2) {
        if (pond == 31) {
            treatment <- "Control"
        } else if (pond == 30) {
            treatment <- "Saw"
        } else if (pond == 27) {
            treatment <- "100Hp"
        } else if (pond == 26) {
            treatment <- 'Square'
        }
    }
    if (trial == 3) {
        if (pond == 31) {
            treatment <- 'Square'
        } else if (pond == 30) {
            treatment <- "Control"
        } else if (pond == 27) {
            treatment <- "Saw"
        } else if (pond == 26) {
            treatment <- "100Hp"
        }
    }
    if (trial == 4) {
        if (pond == 31) {
            treatment <- "100Hp"
        } else if (pond == 30) {
            treatment <- 'Square'
        } else if (pond == 27) {
            treatment <- 'Control'
        } else if (pond == 26) {
            treatment <- "Saw"
        }
    }
    if (trial == 5) {
        if (pond == 31) {
            treatment <- 'Square'
        } else if (pond == 30) {
            treatment <- "Saw"
        } else if (pond == 27) {
            treatment <- "100Hp"
        } else if (pond == 26) {
            treatment <- "Control"
        }
    }
    
    treatment
}