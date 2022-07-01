library(momentuHMM)
library(parallel)
library(tidyverse)


fit.crw <- function(telemetry, 
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
    tel <- telemetry %>%
        select(tel_colnames) %>%
        drop_na()
        
    colnames(tel) <- crw_colnames
    
    # get unique IDs to define batches for crawlWrap
    ids <- unique(tel[,crw_colnames[1]])
    
    # crawlWrap does best with fewer tracks to fit; this generates smaller batches of tracks
    batch_counter <- 0
    predictions <- data.frame(matrix(nrow=0, ncol=ncol(tel)))
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
        
        batch_predictions <- batch_crw$crwPredict
        
        # extract the covariates from telemetry to merge with crw predictions
        covData <- telemetry[, c(crw_colnames[1:2], telemetryCovs)]
        
        # merge covariates with predictions
        batch_predictions <- merge(batch_predictions, covData, by=crw_colnames[1:2])
        
        predictions <- rbind(predictions, batch_predictions)
        batch_counter <- batch_counter + id_batch_size
    }
    
    #prepare data for input into movement model and define covariates
    ModDat <- prepData(data=predictions, covNames=covariates)
    
    return(ModDat)
}

compile.crw <- function(on_time, trials = c(1, 2, 3, 4, 5), path="~/Carp-Model/Fitted CRWs"){
    # compiles various fitted random walks into a single dataset, to be used for HMM fitting.
    
    dt_str <- time.to.str(on_time)
    date_str <- as.character(date(on_time))
    
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
    rm(ModDat)
    for (i in 2:length(pond_files)){
        load(pond_files[i])
        df0 <- ModDat
        rm(ModDat)
        df <- rbind(df, df0)
    }
    
    return(df)
}

convert.coords <- function(df, 
                           input_crs = CRS("+proj=utm +zone=15 +datum=WGS84 +units=m +ellps=WGS84"),
                           output_crs = CRS("+proj=longlat +datum=WGS84")){
    # takes a dataframe of x,y (Easting, Northing) points and converts them to output projection.
    
    x <- colnames(df)[1]
    y <- colnames(df)[2]
    coordinates(df) <- ~x + y
    proj4string(df) <- input_crs
    
    df_out <- as.data.frame(spTransform(df, output_crs)@coords)
    
    return(df_out)
}

db.column <- function(.data, db_data_path, trials=1:5, ponds=c(26, 27, 30, 31),
                      crs_string="+proj=utm +zone=15 +ellps=WGS84 +datum=WGS84 +units=m"){
    # this function performs autokriging on sound intensity data, and predicts dB levels at the 
    # appropriate coordinates in crw.
    # "crw" should be the output (or subset thereof) of the fit.crw function, and the files in
    # db_data_path should be .cvs files with columns "x", "y", and "dB", with names like
    # PondXXTreatment, i.e., Pond27ChirpSquare.
    
    df_holder <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(df_holder) <- c("x", "y", "dB")
    
    .data$dB <- 0
    
    for (trial in trials){
        for (pond in ponds){
            tmnt <- treatment.key(trial, pond)
            if (tmnt == "Control"){next}
            db_data <- read.csv(file.path(db_data_path, paste0("Pond", pond, tmnt, ".csv")))
            sub_data <- .data[.data$Trial == trial & crw$Pond == pond,]
            if (nrow(sub_data) == 0){
                print(paste0("Trial ", trial, ", Pond ", pond, " have no data."))
                next
            }
            sub_data$dB <- fit.krig(db_data, sub_data[, c("x", "y")])$dB
            
            df_holder <- rbind(df_holder, sub_data)
        }
    }
    if (nrow(.data) != nrow(df_holder)){print("Number of rows of input/output data do not match, 
                                            when they should.")}
    
    return(df_holder)
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

fit.krig <- function(sound_data, pred_data, 
                     crs_string="+proj=utm +zone=15 +ellps=WGS84 +datum=WGS84 +units=m"){
    # this function performs an autoKriging on new_data, and extracts dataframe. Assumes labels of
    # x, y, and dB.
    
    sf_sound <- st_as_sf(sound_data, coords = c("x", "y"), crs = CRS(crs_string))
    sp_pred_data <- SpatialPoints(as.data.frame(pred_data), proj4string = CRS(crs_string))
    
    fit_KRIG <- automap::autoKrige(
        formula = dB ~ 1,
        input_data = as(sf_sound, "Spatial"),
        new_data = sp_pred_data
    ) %>%
        .$krige_output %>%
        as.data.frame() %>%
        dplyr::select(x, y, dB = var1.pred)
    
    return(fit_KRIG)
}