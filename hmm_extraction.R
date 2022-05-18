library(momentuHMM)
library(readr)
library(dplyr)

##### These functions extract info from the fitted HMM models

compile.AIC <- function(path, verbose=FALSE){
    # extract all AIC scores from files in Repetition X folders
    tic <- Sys.time()
    
    files <- list.files(path)
    aic_holder <- data.frame(formula=character(0),
                                nCov=integer(0),
                                AIC=numeric(0))

    for (f in files){
        nCov <- parse_number(f)
        hmm_list <- readRDS(paste0(path, f))
        nHMM <- length(hmm_list)
        
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
    
    aic_holder$Rank <- rank(aic_holder$AIC)
    aic_holder <- aic_holder[order(aic_holder$Rank),]
    aic_holder$DeltaAIC <- 0
    aic_holder$DeltaAIC[2:nrow(aic_holder)] <- diff(aic_holder$AIC)
    aic_holder <- relocate(aic_holder, Rank)
    
    toc <- Sys.time()
    t_elapsed <- toc - tic
    
    if (verbose){
        print(paste0("Finished computing AIC values for HMMs; time elapsed: ", 
                     round(t_elapsed), " seconds."))
    }
    
    return(aic_holder)
}

top_models <- function(reps=1:24, aic_scores_path, fitted_hmm_path){
    ##### after running compile.AIC above, this function extracts the top model, and adds it to a
    ##### list indexed by the repetition number. The argument aic_scores_path should point to the
    ##### directory containing the output of compile.AIC, and fitted_hmm_path should point to the
    ##### directory which contains the "nCov covariates.RDS" files/
    
    for (rep in reps){
        aic_scores <- read.csv(paste0(aic_scores_path, rep, ".csv"))
        top_model_info <- aic_scores[aic_scores$Rank == 1,]
        nCov <- top_model_info$nCov
        top_formula <- top_model_info$formula
        fitted_models <- readRDS(paste0(fitted_hmm_path, rep, "/", nCov, " covariates.RDS"))
        n_models <- length(fitted_models)
        for (mod in 1:n_models){
            mod_formula <- format(fitted_models[[mod]]$conditions$formula)
            mod_formula <- gsub(" ", "", mod_formula, fixed=TRUE)
            if (mod_formula == top_formula){
                top_hmms[[rep]] <- fitted_models[[mod]]
            }
        }
        print(paste0("Repetition ", rep, " complete."))
    }
}
