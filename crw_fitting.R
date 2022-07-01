library(momentuHMM)
library(parallel)
library(tidyverse)
source('model_function.R')


crw.prediction <- function(.data, 
                           tel_colnames=c('ID', 'Time', 'Easting', 'Northing'),
                           crw_colnames=c('ID', 'Time', 'x', 'y'),
                           telemetryCovs=c('Trial', 'Pond'),
                           timestep="6 sec",
                           id_batch_size=10,
                           inits=c(2, 0.001),  
                           retry_fits=100, 
                           attempts=100, 
                           doParallel = TRUE){
    # this function loads sound and processed telemetry data, and fits correlated random-walks to
    # the tracks.
    
    # remove rows with missing locations, subset and rename columns
    tel <- .data %>%
        select(tel_colnames) %>%
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
        covData <- telemetry[telemetry$ID %in% batch_ids, c(crw_colnames[1:2], telemetryCovs)]

        # merge covariates with predictions
        batch_predictions <- merge(batch_predictions, covData, by=crw_colnames[1:2])
        
        predictions <- rbind(predictions, batch_predictions)
        batch_counter <- batch_counter + id_batch_size
    }
    
    return(predictions)
}

add.treatment <- function(.data, 
                          colname="Treatment"){
    # this function adds a column for treatment type
    
    trials <- unique(.data$Trial)
    ponds <- unique(.data$Pond)
    .data[[colname]] <- "placeholder"
    
    for (trial in trials){
        for (pond in ponds){
            .data[.data$Trial == trial & .data$Pond == pond, colname] <- treatment.key(trial, pond)
        }
    }
    
    .data[[colname]] <- as.factor(.data[[colname]])
    
    return(.data)
}

add.db <- function(.data, 
                   db_data_path,
                   colname="dB",
                   crs_string="+proj=utm +zone=15 +ellps=WGS84 +datum=WGS84 +units=m"){
    # this function performs autokriging on sound intensity data, predicts dB levels at the 
    # appropriate coordinates in .data, and outputs .data with a a column of those dB values.
    # ".data" should be the output (or subset thereof) of the fit.crw function, and the files in
    # db_data_path should be .cvs files with columns "x", "y", and "dB", with names like
    # PondXXTreatment, i.e., Pond27ChirpSquare.
    
    # create column and determine trials/ponds in the provided data
    .data$dB <- 0
    trials <- unique(.data$Trial)
    ponds <- unique(.data$Pond)
    
    # loop through trials/ponds and update dB column with kriging predictions
    for (trial in trials){
        for (pond in ponds){
            tmnt <- treatment.key(trial, pond)
            if (tmnt == "Control"){next}
            db_data <- read.csv(file.path(db_data_path, paste0("Pond", pond, tmnt, ".csv")))
            sub_data <- .data[.data$Trial == trial & .data$Pond == pond,]
            if (nrow(sub_data) == 0){
                print(paste0("Trial ", trial, ", Pond ", pond, " has no data."))
                next
            }
            .data[.data$Trial == trial & .data$Pond == pond, "dB"] <- fit.krig(db_data, 
                                                                               sub_data[, c("x", "y")])$dB
        }
    }
    
    return(.data)
}


diel.column <- function(.data, time_name="Time", dn_vals=c(0, 1), 
                        latlong=c(38.9122061924, -92.2795993947), timezone="America/Chicago"){
    # this outputs a column of diel values (should be rewritten to output df with the Diel column).
    
    .data$Diel <- dn_vals[2]
    df_dates <- unique(as.Date(.data[,time_name]))
    
    for (the_date in df_dates){
        the_date <- as.Date(the_date)
        sunRS <- sunrise.set(latlong[1], latlong[2],
                             paste0(year(the_date), '/', month(the_date), '/', day(the_date)),
                             timezone=timezone)
        sunRise <- as.POSIXct(sunRS[,1], origin="1970-01-01", tz = timezone)
        sunSet <- as.POSIXct(sunRS[,2], origin="1970-01-01", tz = timezone)
        
        .data[sunRise < .data[,time_name] & .data[, time_name] < sunSet, "Diel"] <- dn_vals[1]
    }
    
    return(.data)
}