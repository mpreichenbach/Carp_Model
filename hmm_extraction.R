library(momentuHMM)
library(readr)
library(dplyr)

##### These functions extract info from the fitted HMM models

compile.aic <- function(path, verbose=FALSE){
    # extract all AIC scores from files in Repetition X folders
    tic <- Sys.time()
    
    files <- list.files(path)
    aic_holder <- data.frame(formula=character(0),
                                nCov=integer(0),
                                AIC=numeric(0))

    for (f in files){
        nCov <- parse_number(f)
        hmm_list <- readRDS(file.path(path, f))
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


top.formulas <- function(model_lists, model_list_names=c("5min", "30min")){
    # takes the top formulas in different b/a time frames, and merges them by repetition. Allows for
    # comparison of which covariates are more effective in the long/short runs. For example, if
    # model_lists=c(top_hmms_5min, top_hmms_30), this will output a dataframe with a column for each
    # input list.
    
    # check whether the model_list entries have the name number of models
    n_model_lists <- length(model_lists)
    n_models <- c()
    for (n in 1:n_model_lists){
        n_models <- append(n_models, length(model_lists[[n]]))
    }
    if (length(unique(n_models)) != 1){stop("Each entry in model_lists must be the same length.")}
    
    n_reps <- length(model_lists[[1]])
    
    # create and fill comparison dataframe
    df <- data.frame(matrix(nrow=n_reps, ncol=length(model_lists) + 1))
    colnames(df) <- append("Repetition", model_list_names)
    df$Repetition <- 1:n_reps
    
    for (rep in 1:n_reps){
        for (i in 1:length(model_list_names)){
            name <- model_list_names[i]
            df[rep, name] <- gsub(" ", "", deparse1(model_lists[[i]][[rep]]$conditions$formula))
        }
    }
    return(df)
}


top.models <- function(aic_scores_path, fitted_hmm_path, reps=1:24, verbose=TRUE){
    ##### after running compile.AIC above, this function extracts the top model, and adds it to a
    ##### list indexed by the repetition number. The argument aic_scores_path should point to the
    ##### directory containing the output of compile.AIC, and fitted_hmm_path should point to the
    ##### directory which contains the "Repetition n" folders.
    
    top_hmms = list()
    for (rep in reps){
        aic_scores <- read.csv(paste0(aic_scores_path, "Repetition ", rep, ".csv"))
        top_model_info <- aic_scores[aic_scores$Rank == 1,]
        nCov <- top_model_info$nCov
        top_formula <- top_model_info$formula
        fitted_models <- readRDS(paste0(fitted_hmm_path, "Repetition ", rep, "/", nCov, 
                                        " covariates.RDS"))
        n_models <- length(fitted_models)
        for (mod in 1:n_models){
            mod_formula <- format(fitted_models[[mod]]$conditions$formula)
            mod_formula <- gsub(" ", "", mod_formula, fixed=TRUE)
            if (mod_formula == top_formula){
                top_hmms[[rep]] <- fitted_models[[mod]]
            }
        }
        if (verbose){print(paste0("Repetition ", rep, " complete."))}
    }
    return(top_hmms)
}
