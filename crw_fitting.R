library(hms)
library(momentuHMM)
library(parallel)
library(suncalc)
library(stats)
library(tidyverse)
source('model_functions.R')


predict_crw <- function(.data, 
                        tel_colnames=c('ID', 'Time', 'Easting', 'Northing'),
                        crw_colnames=c('ID', 'Time', 'x', 'y'),
                        telCovs=c('Trial', 'Pond'),
                        timestep="6 sec",
                        id_batch_size=10,
                        inits=c(2, 0.001),  
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
    
    while (batch_counter < length(ids)){
        batch_ids <- ids[(1 + batch_counter):(id_batch_size + batch_counter)]
        batch_tel <- tel[tel[[crw_colnames[1]]] %in% batch_ids,]
        
        batch_crw <- crawlWrap(obsData=batch_tel, 
                               timeStep=timestep,
                               theta=inits,
                               fixPar=c(NA, NA),
                               Time.name = crw_colnames[2],
                               retryFits = retry_fits,
                               attempts=attempts,
                               doParallel = doParallel,
                               ncores = min(id_batch_size, ceiling(0.75 * detectCores()))
                               )
        
        batch_predictions <- prepData(batch_crw)

        # extract the covariates from telemetry to merge with crw predictions
        covData <- .data[.data$ID %in% batch_ids, c(crw_colnames[1:2], telCovs)]

        # merge covariates with predictions
        batch_predictions <- merge(batch_predictions, covData, by=crw_colnames[1:2])
        
        predictions <- rbind(predictions, batch_predictions)
        batch_counter <- batch_counter + id_batch_size
    }
    
    return(predictions)
}


add_diel <- function(.data, 
                     colname="Diel",
                     time_name="Time", 
                     day_night_values=c("Day", "Night"), 
                     latlong=c(38.9122061924, -92.2795993947), 
                     timezone="America/Chicago"){
    # adds a column of day/night values given in dn_vals vector
    
    .data[[colname]] <- day_night_values[2]
    all_dates <- unique(as.Date(.data[,time_name]))
    
    for (the_date in all_dates){
        the_date <- as.Date(the_date)
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
                          new_colname = "dB",
                          sound_samples_path = "~/Carp-Model/Supplementary Files/Sound Mapping/Data_Master_PosUp_Compiled.csv",
                          treatments = list(c("ChirpSaw", "Saw"), 
                                          c("ChirpSquare", "Square"), 
                                          c("BoatMotor", "100Hp")),
                          input_names = c("Pond",
                                          "Long", 
                                          "Lat", 
                                          "Actual.Signal", 
                                          "RMS.SPL..dB.re.1uPa."),
                          convert_to_utm = TRUE,
                          crs_string = "+proj=utm +zone=15 +ellps=WGS84 +datum=WGS84 +units=m") {
    # this function performs autokriging on sound intensity data, predicts dB levels at the 
    # appropriate coordinates in .data, and outputs .data with a a column of those dB values.
    # ".data" should be the output (or subset thereof) of the fit.crw function, and the files in
    # db_data_path should be .cvs files with columns "x", "y", and "dB", with names like
    # PondXXTreatment, i.e., Pond27ChirpSquare.
    
    # create column and load the sound samples
    .data[[new_colname]] <- 0
    df_sound <- read.csv(sound_samples_path)
    df_sound[df_sound == "*"] <- NA
    df_sound <- drop_na(df_sound)
    
    # trials <- unique(.data$Trial)
    ponds <- unique(.data$Pond)
    
    # convert coordinates if necessary
    if (convert_to_utm) {
        df_sound[, input_names[2:3]] <- convert_coords(df_sound[, input_names[2:3]],
                                                       output_crs=crs_string)
    }
    
    # rename verious columns
    colnames(df_sound)[which(colnames(df_sound) == input_names[2])] <- "x"
    colnames(df_sound)[which(colnames(df_sound) == input_names[3])] <- "y"
    colnames(df_sound)[which(colnames(df_sound) == input_names[5])] <- "dB"
    
    # loop through trials/ponds and update new column with kriging predictions
    for (pond in ponds) {
        for (tmnt_pair in treatments) {
            tel_tmnt <- tmnt_pair[1]
            sound_tmnt <- tmnt_pair[2]
            
            # define telemetry sub-setting condition
            tel_subset_bool <- .data$Sound == "on" & 
                .data$Pond == pond & 
                .data$Treatment == tel_tmnt
            
            # define sound sub-setting condition
            sound_subset_bool <- df_sound[[input_names[1]]] == pond & 
                df_sound[[input_names[4]]] == sound_tmnt
            
            # get only the relevant data for this iteration
            tel_sub <- .data[tel_subset_bool, ]
            sound_sub <- df_sound[sound_subset_bool,]
            
            if (nrow(tel_sub) == 0){
                print(paste0("Treatment ", tel_tmnt, ", Pond ", pond, " has no data."))
                next
            }
            
            .data[tel_subset_bool, new_colname] <- fit_krig(sound_sub,
                                                            pred_data = tel_sub[, c("x", "y")],
                                                            crs_string = crs_string)$dB
            
        }
    }
    
    .data
    
    # for (trial in trials){
    #     for (pond in ponds){
    #         tmnt <- treatment.key(trial, pond)
    #         if (tmnt == "Control"){next}
    #         db_data <- read.csv(file.path(intensity_data_path, paste0("Pond", pond, tmnt, ".csv")))
    #         subset_condition <- .data$Sound == "on" & .data$Trial == trial & .data$Pond == pond
    #         sub_data <- .data[subset_condition, ]
    #         if (nrow(sub_data) == 0){
    #             print(paste0("Trial ", trial, ", Pond ", pond, " has no data."))
    #             next
    #         }
    #         .data[subset_condition, colname] <- fit.krig(db_data, sub_data[, c("x", "y")])$dB
    #     }
    # }
    # return(.data)
}


add_temperature <- function(.data,
                            temperature_data_path,
                            colname="Temperature",
                            input_colnames=c("DateTime", "Temp_C"),
                            timezone="America/Chicago",
                            trials=1:5,
                            ponds=c(26, 27, 30, 31)) {
    # adds a column with temperature values. Since the data I'm using has frequent measurements (15
    # minutes) and low variation, this performs linear interpolation. However, this may not be a 
    # good model for less frequent measurements (where a sinusoidal model may be better).
    
    temperature_data <- read.csv(temperature_data_path)
    temperature_data[[input_colnames[1]]] <- as.POSIXct(temperature_data[[input_colnames[1]]],
                                                        tz=timezone)
    temperature_data <- temperature_data[order(temperature_data[[input_colnames[1]]]),]
    
    
    # first check whether there is sufficient temperature data
    if (min(.data$Time) < min(temperature_data[[input_colnames[1]]]) || 
        max(.data$Time) > max(temperature_data[[input_colnames[1]]])){
        stop("Insufficient temperature data; check min/max times of the random walk.")
    }
    
    # create column
    .data[[colname]] <- NA
    
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
                
                .data[[colname]][subset_condition] <- predict(fit.lm, 
                                                              newdata=data.frame(x=.data[subset_condition, "Time"]))
            }
        }
    }
    
    if (any(is.na(.data[[colname]]))){stop("There are still NA's in the temperature column.")}
    
    return(.data)
}

add_sound <- function(.data,
                      sound_data,
                      rep_number,
                      colname = "Sound",
                      values = c("on", "off")) {
    
    # adds a column to indicate before/after the sound plays
    
    # extract the relevant time
    on_time <- unique(as_hms(sound_data[sound_data$Repetition == rep_number, "Time"]))
    
    # assign the on/off values
    .data[, colname] <- values[2]
    .data[as_hms(.data$Time) >= on_time, colname] <- values[1]
    
    return(.data)
}

add_treatment <- function(.data, 
                          colname="Treatment") {
    # this function adds a column for treatment type
    
    .data[[colname]] <- "placeholder"
    
    trials <- unique(.data$Trial)
    ponds <- unique(.data$Pond)
    
    for (trial in trials){
        for (pond in ponds){
            .data[.data$Trial == trial & .data$Pond == pond, colname] <- treatment.key(trial, pond)
        }
    }
    
    .data[[colname]] <- as.factor(.data[[colname]])
    
    return(.data)
}

