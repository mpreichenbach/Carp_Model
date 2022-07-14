library(momentuHMM)
library(readr)
library(dplyr)

##### These functions extract info from the fitted HMM models

compile_aic_scores <- function(path, verbose=FALSE){
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


top_formulas <- function(model_lists, model_list_names=c("5min", "30min")){
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


top_models <- function(aic_scores_path, fitted_hmm_path, reps=1:24, verbose=TRUE){
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

get_mean_estimates <- function(hmm_list,
                               save_dir = NA,
                               parms = c("step", "angle"),
                               factor_covs = c("Trial", "Pond", "Treatment"),
                               state_names = c("exploratory", "encamped"),
                               state_colors = c("#E69F00", "#56B4E9"),
                               plotCI = TRUE,
                               verbose = TRUE){
    # plots the means of the fitted distributions for values in the parms arguments. It would be
    # nice to plot mean step lengths over all factor covariates; but this requires some fancy math
    # I don't yet know.
    
    n_models <- length(hmm_list)
    data <- hmm_list[[1]]$data
    
    # the entries of this list are vectors of the unique values for each factor covariate
    factor_values <- list()
    
    for (fac in factor_covs) {
        factor_values[[fac]] <- unique(data[[fac]])
    }
    
    
    # this dataframe holds every combination of parameter estimates, with full name
    df_colnames <- expand.grid(list(state = state_names, parm = parms, 
                                    val = c("lower", "est", "upper")))
    
    df_colnames$FullName <- do.call(paste, c(df_colnames[c("state", "parm", "val")]))
    
    # very slow loop which extracts parameter estimates for each model
    for (i in 1:n_models) {
        hmm <- hmm_list[[i]]
        
        # this dataframe holds every combination of the factor covariates
        df_factors <- expand.grid(factor_values)
        
        # initialize the parameter value columns
        for (fullname in df_colnames$FullName) {
            df_factors[[fullname]] <- 0.0
        }
        
        # extract the predicted values for each covariate combination (necessarily row-by-row; v. slow)
        for (j in 1:nrow(df_factors)) {
            # this step must proceed row-by-row
            estimates <- CIreal(hmm, covs = df_factors[j, factor_covs], parms = parms)
            for (fullname in df_colnames$FullName) {
                split_names <- df_colnames[df_colnames$FullName == fullname, ]
                state <- as.character(split_names[1, "state"])
                parm <- as.character(split_names[1, "parm"])
                val <- as.character(split_names[1, "val"])
                
                df_factors[j, fullname] <- as.data.frame(estimates[[parm]][[val]])[1, state]
            }
        }
        if (verbose) {print(paste0("Finished extracting parameters for model ", i, "/", n_models))}
    }
    
    df_factors
}