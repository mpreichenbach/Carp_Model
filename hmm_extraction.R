library(momentuHMM)
libary(readr)

##### These functions extract info from the fitted HMM models

compile.AIC <- function(path){
    # extract all AIC scores from files in Repetition X folders
    
    files <- list.files(path)
    aic_holder <- as.data.frame(formula=character(0),
                                nCov=integer(0),
                                AIC=numeric(0))

    for (file in files){
        nCov <- parse_number(file)
        hmm_list <- readRDS(paste0(path, file))
        nHMM <- length(hmm_list)
        
        # collect the formulas from the HMMs into a character vector
        formulas <- c()
        aic_values <- c()
        for (i in 1:nHMM){
            frm <- format(hmm_list[[i]]$conditions$formula)
            formulas <- append(formulas, frm)
            aic_values <- append(aic_values, AIC(hmm_list[[i]]))
        }
        
        # dataframe of model formulas and AIC scores for a given number of covariates
        sub_aic_holder <- as.data.frame(formula=formulas,
                                        nCov=rep(nCov, nHMM),
                                        AIC=aic_values)
        
        aic_holder <- rbind(aic_holder, sub_aic_holder)
    }
    
    return(aic_holder)
}