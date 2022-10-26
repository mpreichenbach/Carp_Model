library(hms)
library(momentuHMM)
library(parallel)
library(suncalc)
library(stats)
library(tidyverse)
source('model_functions.R')


predict_crw <- function(.data, 
                        tel_colnames=c("ID", "Time", "Easting", "Northing"),
                        crw_colnames=c("ID", "Time", "x", "y"),
                        telCovs=c("Trial", "Pond", "Repetition"),
                        timestep="6 sec",
                        id_batch_size=12,
                        inits=c(1, 0.01),  
                        retry_fits=100, 
                        attempts=100, 
                        doParallel = TRUE) {
    
    # this function loads sound and processed telemetry data, and fits correlated random-walks to
    # the tracks.
    
    # remove rows with missing locations, subset and rename columns
    tel <- .data %>%
        select(all_of(tel_colnames)) %>%
        drop_na()
        
    colnames(tel) <- crw_colnames
    
    # get unique IDs to define batches for crawlWrap
    ids <- unique(tel[,crw_colnames[1]])
    
    # crawlWrap does best with fewer tracks to fit; this generates smaller batches of tracks
    batch_counter <- 0
    predictions <- data.frame(matrix(nrow=0, ncol=ncol(tel) + 2))
    colnames(predictions) <- colnames(tel)
    
    # since all batches may not fit, this object will 
    error_list <- list()
    
    while (batch_counter < length(ids)){
        batch_ids <- ids[(1 + batch_counter):(id_batch_size + batch_counter)]
        batch_tel <- tel[tel[[crw_colnames[1]]] %in% batch_ids,]
        
        # with the try() function, this will output either a CRW or an error message
        crw_attempt <- try({
            batch_crw <- crawlWrap(obsData=batch_tel, 
                                   timeStep=timestep,
                                   theta=inits,
                                   fixPar=c(NA, NA),
                                   Time.name = crw_colnames[2],
                                   retryFits = retry_fits,
                                   attempts=attempts,
                                   doParallel = doParallel,
                                   ncores = min(id_batch_size, ceiling(0.75 * detectCores())))
            
        })
        
        # append CRWs to a list of batches, and an error to a list of errors
        if ("crwData" %in% class(crw_attempt)){
            batch_predictions <- prepData(crw_attempt)
        } else {
            error_list <- append(error_list, crw_attempt)
        }

        # extract the covariates from telemetry to merge with crw predictions
        covData <- .data[.data$ID %in% batch_ids, c(crw_colnames[1:2], telCovs)]

        # merge covariates with predictions
        batch_predictions <- merge(batch_predictions, covData, by=crw_colnames[1:2])
        
        predictions <- rbind(predictions, batch_predictions)
        batch_counter <- batch_counter + id_batch_size
    }
    
    out_list = list("pred"=predictions, "err"=error_list)
    
    return(out_list)    
}


add_diel <- function(.data, 
                     time_name="Time", 
                     day_night_values=c("Day", "Night"), 
                     latlong=c(38.9122061924, -92.2795993947), 
                     timezone="America/Chicago"){
    # adds a column of day/night values given in day_night_values vector
    
    .data$Diel <- day_night_values[2]
    all_dates <- unique(as.Date(.data[,time_name]))

    for (the_date in all_dates) {
        the_date <- as.Date(the_date, origin="1970-01-01")
        sun_times <- getSunlightTimes(date=the_date,
                                      lat=latlong[1],
                                      lon=latlong[2],
                                      keep=c("sunrise", "sunset"),
                                      tz=timezone)
        sunRise <- sun_times[, "sunrise"]
        sunSet <- sun_times[, "sunset"]
        
        .data[sunRise < .data[,time_name] & .data[, time_name] < sunSet, "Diel"] <- day_night_values[1]
    }
    
    return(.data)
}


add_intensity <- function(.data, 
                          sound_samples_path="C:/Users/RDGRLMPR/r_dev/Carp-Model/Supplementary Files/Sound Mapping/Data_Master_PosUp_Compiled.csv",
                          silent_treatments = c("Control", "Silence"),
                          input_names = c("Pond",
                                          "Long", 
                                          "Lat", 
                                          "Actual.Signal", 
                                          "RMS.SPL..dB.re.1uPa."),
                          convert_to_utm = TRUE,
                          crs_string = "+proj=utm +zone=15 +ellps=WGS84 +datum=WGS84 +units=m") {
    # This function performs autokriging on sound intensity data, predicts dB levels at the 
    # appropriate coordinates in .data, and outputs .data with a a column of those dB values.
    # ".data" should be the output (or subset thereof) of the fit.crw function, piped through 
    # add_sound and add_treatment. 
    
    # The argument sound_samples_path should point to a single CSV file with pond, location, 
    #treatment type, and dB levels. 
    
    # The argument input_names should be a 5-element vector with pond, x, y, treatment type, and dB
    # column names (in that exact order) from the dataframe at sound_samples_path.
    
    # create column and load the sound samples
    .data$dB <- 0 
    df_sound <- read.csv(sound_samples_path)
    df_sound[df_sound == "*"] <- NA
    df_sound <- drop_na(df_sound)
    
    ponds <- as.numeric(unique(.data$Pond))
    treatments <- setdiff(unique(.data$Treatment), silent_treatments)
    
    # convert coordinates if necessary
    if (convert_to_utm) {
        df_sound[, input_names[2:3]] <- convert_coords(df_sound[, input_names[2:3]],
                                                       output_crs=crs_string)
    }
    
    # rename various columns
    colnames(df_sound)[which(colnames(df_sound) == input_names[1])] <- "Pond"
    colnames(df_sound)[which(colnames(df_sound) == input_names[2])] <- "x"
    colnames(df_sound)[which(colnames(df_sound) == input_names[3])] <- "y"
    colnames(df_sound)[which(colnames(df_sound) == input_names[4])] <- "Treatment"
    colnames(df_sound)[which(colnames(df_sound) == input_names[5])] <- "dB"
    
    # loop through trials/ponds and update new column with kriging predictions
    for (pond in ponds) {
        
        # define the minimum intensity for "Silence" and "Control" treatments
        min_intensity <- min(df_sound[df_sound$Pond == pond, "dB"])
        .data[.data$Pond == pond & .data$Treatment %in% silent_treatments, "dB"] <- min_intensity
        
        for (tmnt in treatments) {
            # get only the relevant data for this iteration
            tel_tmnt_bool <- .data$Pond == pond & .data$Treatment == tmnt
            tel_sub <- .data[tel_tmnt_bool, ]
            sound_sub <- df_sound[df_sound$Pond == pond & df_sound$Treatment == tmnt, ]

            if (nrow(tel_sub) == 0){
                print(paste0("Treatment ", tmnt, ", Pond ", pond, " has no data."))
                next
            }
            
            # fit kriging model and update the sound-on dB values
            .data[tel_tmnt_bool, "dB"] <- fit_krig(sound_sub,
                                                   pred_data = tel_sub[, c("x", "y")],
                                                   crs_string = crs_string)$dB
        }
    }
    
    # no outputted dB values should be 0
    if (min(.data$dB == 0)) {
        stop("The minimum dB value of the dataset is 0, but this should not be the case.")
    }
    
    return(.data)
}


add_sound <- function(.data,
                      values = c("on", "off")) {
    
    # adds a column to indicate before/after the sound plays
    
    # initialize the column
    .data$Sound <- values[2]
    
    # get the repetition numbers
    rep_numbers <- unique(.data$Repetition)

    # load the sound sample data
    df_sound <- sound_data()
    
    # extract the relevant times
    on_times <- unique(as_hms(df_sound[df_sound$Repetition %in% rep_numbers, "Time"]))
    
    for (rep in rep_numbers) {
        on_time <- unique(as_hms(df_sound[df_sound$Repetition == rep, "Time"]))
        .data[(as_hms(.data$Time) >= on_time) & (.data$Repetition == rep), "Sound"] <- values[1]
    }

    return(.data)
}


add_temperature <- function(.data,
                            temperature_path="C:/Users/RDGRLMPR/r_dev/Carp-Model/Supplementary Files/2018CERC_WaterTemperature_AllTrialsPonds.csv",
                            input_colnames=c("DateTime", "Temp_C"),
                            timezone="America/Chicago",
                            trials=1:5,
                            ponds=c(26, 27, 30, 31)) {
    # adds a column with temperature values. Since the data I'm using has frequent measurements (15
    # minutes) and low variation, this performs linear interpolation. However, this may not be a 
    # good model for less frequent measurements (where a sinusoidal model may be better).
    
    temperature_data <- read.csv(temperature_path)
    temperature_data[[input_colnames[1]]] <- as.POSIXct(temperature_data[[input_colnames[1]]],
                                                        tz=timezone)
    temperature_data <- temperature_data[order(temperature_data[[input_colnames[1]]]),]
    
    
    # first check whether there is sufficient temperature data
    if (min(.data$Time) < min(temperature_data[[input_colnames[1]]]) || 
        max(.data$Time) > max(temperature_data[[input_colnames[1]]])){
        stop("Insufficient temperature data; check min/max times of the random walk.")
    }
    
    # create column
    .data$Temperature <- NA
    
    # predict temperature values
    for (trial in trials){
        for (pond in ponds){
            # subset the positions
            df_pos <- .data[.data$Trial == trial & .data$Pond == pond,]

            # subset the temperature data
            df_temp <- temperature_data[temperature_data$Trial == trial & 
                                            temperature_data$Pond == pond,]
            
            # fit a line to consecutive measurements, and predict temperatures
            for (i in 1:(nrow(df_temp) - 1)){
                y <- df_temp[i:(i + 1), input_colnames[[2]]]
                x <- df_temp[i:(i + 1), input_colnames[[1]]]
                
                # skip past times intervals which don't intersect with the position data
                if (x[2] < min(df_pos$Time) || max(df_pos$Time) < x[1]){next}
                
                fit.lm <- lm(y~x)
                
                # yields the rows to update
                subset_condition <- x[1] <= .data$Time & 
                    .data$Time < x[2] & 
                    .data$Trial == trial & 
                    .data$Pond == pond
                
                predicted <- predict(fit.lm, newdata=data.frame(x=.data[subset_condition, "Time"]))
                
                .data[subset_condition, "Temperature"] <- predicted
            }
        }
    }
    
    if (any(is.na(.data$Temperature))){stop("There are still NA's in the temperature column.")}
    
    return(.data)
}


add_treatment <- function(.data) {
    # this function adds a column for treatment type; requires add_sound to be run first
    
    .data$Treatment <- "Silence"
    
    trials <- unique(.data$Trial)
    ponds <- unique(.data$Pond)
    
    for (trial in trials){
        for (pond in ponds){
            tmnt <- treatment_key(trial, pond)
            .data[.data$Trial == trial
                  & .data$Pond == pond
                  & .data$Sound == "on", "Treatment"] <- tmnt
        }
    }
    
    .data$Treatment <- as.factor(.data$Treatment)
    
    return(.data)
}


concatenate_tracks <- function(df_before, df_after, xycols=c("x", "y")) {
    # this function concatenates the tracks in the before/after dataframes.
    
    x_name = xycols[1]
    y_name = xycols[2]
    
    df_before0 <- df_before
    common_ids <- intersect(drop_na(df_before)$ID, drop_na(df_after)$ID)
    
    for (id in common_ids) {
        before_sub <- df_before[df_before$ID == id, ]
        after_sub <- df_after[df_after$ID == id, ]
        
        # make sure everything is in chronological order
        before_sub <- before_sub[order(before_sub$Time), ]
        after_sub <- after_sub[order(after_sub$Time), ]
        
        # get the starting and ending positions
        before_sub_drop_na <- drop_na(before_sub)
        after_sub_drop_na <- drop_na(after_sub)
        
        # move on to the next ID if there are no records
        if (nrow(after_sub_drop_na) == 0) {
            next
        }
        
        x_ending <- before_sub_drop_na[nrow(before_sub_drop_na), x_name]
        y_ending <- before_sub_drop_na[nrow(before_sub_drop_na), y_name]
        
        x_starting <- after_sub_drop_na[nrow(after_sub_drop_na), x_name]
        y_starting <- after_sub_drop_na[nrow(after_sub_drop_na), y_name]
        
        # update positions
        after_sub[, x_name] <- after_sub[, x_name] - x_starting + x_ending
        after_sub[, y_name] <- after_sub[, y_name] - y_starting + y_ending
        
        # row-bind the repositioned track onto the before-data frame
        df_before0 <- rbind(df_before0, after_sub)
        }
    
    return(df_before0)
}
