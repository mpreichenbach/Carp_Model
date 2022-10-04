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
    # this function performs an autoKriging on new_data, and extracts dataframe. Assumes labels of
    # x, y, and dB.
    
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
                         ignore_covs=c("Diel")) {
    # outputs a list of formulas with the specified number of covariates; here we use the covariates
    # "Trial", "Temperature", "dB", "Treatment". Additionally, we allow the option of includeing
    # "Diel" (aka, day/night) as a covariate.
    
    if (!(nCov %in% c(0, 1, 2, 3, 4, 5, 6))) {
        stop("nCov must be 0, 1, 2, 3, 4, 5, or 6.")
    }
    
    # yield a character vector of formulas with nCov covariates, including "Diel"
    if (nCov == 0) {
        frm_names <- c("~1")
    } else if (nCov == 1) {
        frm_names <- c("~Trial", "~Temperature", "~Diel", "~dB", "~Treatment")
    } else if (nCov == 2) {
        frm_names <- c("~Trial+Temperature", "~Trial+Diel", "~Trial+dB", "~Trial+Treatment",
                   "~Temperature+Diel", "~Temperature+Diel", "~Temperature+dB", 
                   "~Temperature+Treatment", "~Diel+dB", "~Diel+Treatment", "~dB+Treatment",
                   "~dB:Treatment")
    } else if (nCov == 3) {
        frm_names <- c("~Trial+Temperature+Diel", "~Trial+Temperature+dB", 
                   "~Trial+Temperature+Treatment", "~Trial+Diel+dB", "~Trial+Diel+Treatment",
                   "~Trial+dB+Treatment", "~Trial+dB:Treatment", "~Temperature+Diel+dB",
                   "~Temperature+Diel+Treatment", "~Temperature+dB+Treatment", 
                   "~Temperature+dB:Treatment", "~Diel+dB+Treatment", "~Diel+dB:Treatment",
                   "~dB*Treatment")
    } else if (nCov == 4) {
        frm_names <- c("~Trial+Temperature+Diel+dB", "~Trial+Temperature+Diel+Treatment",
                   "~Trial+Temperature+dB+Treatment", "~Trial+Temperature+dB:Treatment",
                   "~Trial+dB*Treatment", "~Temperature+Diel+dB+Treatment",
                   "~Temperature+Diel+dB:Treatment", "~Temperature+dB*Treatment",
                   "~Diel+dB*Treatment")
    } else if (nCov == 5) {
        frm_names <- c("~Trial+Temperature+Diel+dB+Treatment", "~Trial+Temperature+Diel+dB:Treatment",
                       "~Trial+Temperature+dB*Treatment", "~Temperature+Diel+dB*Treatment")
    } else if (nCov == 6) {
        frm_names <- c("~Trial+Temperature+Diel+dB*Treatment")
    }
    
    # sequentially remove covariates from the ignore_covs vector
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


plot_interp <- function(points, 
                        pond, 
                        sound, 
                        grad_min, 
                        grad_max, 
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
    
    plt    
}


plot_tracks <- function(tel_data = NULL, 
                        crw_data = NULL, 
                        id = NULL, 
                        min_time, 
                        max_time, 
                        tel_col = "dodgerblue1", 
                        crw_col = "orange1", full_extent=FALSE){
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

    plot_list    
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

##### interpolation script
# We will need some packages for (spatial) data processing
# library(tidyverse) # wrangling tabular data and plotting
# library(sf) # processing spatial vector data
# library(sp) # another vector data package necessary for continuity
# library(raster) # processing spatial raster data. !!!overwrites dplyr::select!!!
# 
# # And a lot of different packages to test their interpolation functions
# library(gstat)  # inverse distance weighted, Kriging
# library(fields) # Thin Plate Spline
# library(interp) # Triangulation
# library(mgcv)   # Spatial GAM
# library(automap)# Automatic approach to Kriging
# 
# # Finally, some packages to make pretty plots
# library(patchwork)
# library(viridis)
# 
# crs_string <- "+proj=utm +zone=15 +datum=WGS84 +units=m +ellps=WGS84"
# 
# pts_sound <- readr::read_csv(
#     "~/Carp-Model/Supplementary Files/Sound Mapping/UTM, Zone 15/Pond31BoatMotor.csv",
#     col_types = cols(dB = col_double(), 
#                      x = col_double(), y = col_double())
# ) %>% 
#     dplyr::select(x, y, dB)
# 
# bbox <- c(
#     "xmin" = min(pts_sound$x),
#     "ymin" = min(pts_sound$y),
#     "xmax" = max(pts_sound$x),
#     "ymax" = max(pts_sound$y)
# )
# 
# grd_template <- expand.grid(
#     x = seq(from = bbox["xmin"], to = bbox["xmax"], by = 0.05),
#     y = seq(from = bbox["ymin"], to = bbox["ymax"], by = 0.05) # 0.05 m resolution
# )
# 
# sf_sound <- st_as_sf(pts_sound, coords = c("x", "y"), crs = crs_string)
# 
# alt_grd_template_sf <- sf_sound %>% 
#     st_bbox() %>% 
#     st_as_sfc() %>% 
#     st_make_grid(
#         cellsize = c(0.05, 0.05),
#         what = "centers"
#     ) %>%
#     st_as_sf() %>%
#     cbind(., st_coordinates(.)) %>% 
#     st_drop_geometry() %>% 
#     mutate(dB = 0)
# 
# colnames(alt_grd_template_sf) <- c("x", "y",)
# 
# grd_template_raster <- grd_template %>% 
#     dplyr::mutate(dB = 0) %>% 
#     raster::rasterFromXYZ( 
#         crs = crs_string)
# 
# alt_grd_template_raster <- alt_grd_template_sf %>% 
#     raster::rasterFromXYZ(
#         crs = crs_string
#     )
# 
# fit_IDW <- gstat::gstat( # The setup here is quite similar to NN
#     formula = dB ~ 1,
#     data = as(sf_sound, "Spatial"),
#     nmax = 20, nmin = 3,
#     set = list(idp = 0.5) # inverse distance power
# )
# 
# fit_NN <- gstat::gstat( # using package {gstat} 
#     formula = dB ~ 1,    # The column `NH4` is what we are interested in
#     data = as(sf_sound, "Spatial"), # using {sf} and converting to {sp}, which is expected
#     nmax = 20, nmin = 3 # Number of neighboring observations used for the fit
# )
# 
# fit_TPS <- fields::Tps( # using {fields}
#     x = as.matrix(pts_sound[, c("x", "y")]), # accepts points but expects them as matrix
#     Y = pts_sound$dB,  # the dependent variable
#     miles = FALSE     # EPSG 25833 is based in meters
# )
# 
# fit_GAM <- mgcv::gam( # using {mgcv}
#     dB ~ s(x, y),      # here come our X/Y/Z data - straightforward enough
#     data = pts_sound      # specify in which object the data is stored
# )
# 
# fit_TIN <- interp::interp( # using {interp}
#     x = pts_sound$x,           # the function actually accepts coordinate vectors
#     y = pts_sound$y,
#     z = pts_sound$dB,
#     xo = grd_template$x,     # here we already define the target grid
#     yo = grd_template$y,
#     output = "points"
# ) %>% bind_cols()
# 
# sp_new_data <- SpatialPoints(alt_grd_template_sf, proj4string = CRS(crs_string))
# 
# fit_KRIG <- automap::autoKrige(
#     formula = dB ~ 1,
#     input_data = as(sf_sound, "Spatial"),
#     new_data = sp_new_data
# ) %>%
#     .$krige_output %>%
#     as.data.frame() %>%
#     dplyr::select(x, y, dB = var1.pred)
# 
# interp_TIN <- raster::rasterFromXYZ(fit_TIN, crs = crs_string)
# 
# interp_KRIG <- raster::rasterFromXYZ(fit_KRIG, crs = crs_string)
# 
# interp_NN <- interpolate(grd_template_raster, fit_NN)
# 
# interp_IDW <- interpolate(grd_template_raster, fit_IDW)
# 
# interp_TPS <- interpolate(grd_template_raster, fit_TPS)
# 
# interp_GAM <- grd_template %>% 
#     mutate(dB = predict(fit_GAM, .)) %>% 
#     rasterFromXYZ(crs = crs_string)
# 
# plot_my_rasters <- function(raster_object, raster_name){
#     
#     df <- rasterToPoints(raster_object) %>% as_tibble()
#     colnames(df) <- c("X", "Y", "Z")
#     
#     ggplot(df, aes(x = X, y = Y, fill = Z)) +
#         geom_raster() +
#         ggtitle(label = raster_name) +
#         scale_fill_viridis(option = "C") +
#         theme_bw() +
#         theme(
#             axis.text = element_blank(),
#             axis.title = element_blank(),
#             axis.ticks = element_blank()
#         )
# }
# 
# rasterlist <- list(
#     "Nearest Neighbor" = interp_NN, 
#     "Inverse Distance Weighted" = interp_IDW, 
#     "Kriging" = interp_KRIG, 
#     "Thin Plate Spline Regression" = interp_TPS,
#     "Triangular Irregular Surface" = interp_TIN, 
#     "Generalized Additive Model" = interp_GAM
# )
# 
# plotlist <- map2(
#     rasterlist,
#     names(rasterlist),
#     plot_my_rasters
# )
# 
# # Note that the next trick only works because of library(patchwork)
# (plotlist[[1]] + plotlist[[2]]) /
#     (plotlist[[3]] + plotlist[[4]]) /
#     (plotlist[[5]] + plotlist[[6]])
