library(momentuHMM)
library(readr)

##### These functions extract info from the fitted HMM models

compile.AIC <- function(path, verbose=FALSE){
    # extract all AIC scores from files in Repetition X folders
    tic <- Sys.time()
    
    files <- list.files(path)
    aic_holder <- data.frame(formula=character(0),
                                nCov=integer(0),
                                AIC=numeric(0))

    for (f in files){
        print(f)
        nCov <- parse_number(f)
        hmm_list <- readRDS(paste0(path, f))
        print(paste0(path, f))
        nHMM <- length(hmm_list)
        =
        # collect the formulas from the HMMs into a character vector
        formulas <- c()
        aic_values <- c()
        for (i in 1:nHMM){
            frm <- format(hmm_list[[i]]$conditions$formula)
            frm <- gsub(" ", "", frm, fixed=TRUE)
            
            formulas <- append(formulas, frm)
            aic_values <- append(aic_values, AIC(hmm_list[[i]]))
        }
        
        # dataframe of model formulas and AIC scores for a given number of covariates
        sub_aic_holder <- data.frame(formula=formulas,
                                        nCov=rep(nCov, nHMM),
                                        AIC=aic_values)
        
        aic_holder <- rbind(aic_holder, sub_aic_holder)
    }
    
    toc <- Sys.time()
    t_elapsed <- toc - tic
    
    if (verbose){
        print(paste0("Finished computing AIC values for HMMs; time elapsed: ", 
                     round(t_elapsed), " seconds."))
    }
    
    return(aic_holder)
}

mean.positions <- function(path, nCov=0:6){
    # takes the path of the folder containing all AIC Tables, and outputs a dataframe with mean
    # ranking over all reps for each formula (or those specified by nCov).
    
    files <- list.files(path)
    source_f <- read.csv(paste0(path, files[1]))
    source_f <- source_f[source_f$nCov %in% nCov,]
    
    formulas <- source_f$formula
        
    position_holder <- data.frame(matrix(0, nrow=length(formulas), ncol=length(files)),
                                  row.names=formulas)
    colnames(position_holder) <- files
    
    for (f in files){
        csv <- read.csv(paste0(path, f))
        csv <- csv[csv$nCov %in% nCov,]
        # check to make sure that the formulas are in the correct order
        if (!(any(csv$formula == formulas))){
            stop("Formulas have different orders in the AIC tables.")
        }
        position_holder[, f] <- as.numeric(rank(csv$AIC)) 
    }
    
    position_holder$Mean <- rowMeans(position_holder, na.rm=TRUE)
    
    return(position_holder)
}






























