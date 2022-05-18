library(ggplot2)
library(dplyr)

aic.caterpillar <- function(aic_df, loc="mean", colors=c("blue", "red"), descending_rank=TRUE){
    # Generates a caterpillar plot of the AIC scores for the various models fitted to a given
    # dataset. The loc argument stands for "line of comparison", and can take values "mean",
    # "median", and "none". The colors argument should be a vector of strings.
    
    if (loc == "mean"){
        y_loc <- mean(aic_df$AIC)
    }else if (loc == "median"){
        y_loc <- median(aic_df$AIC)
    }
    
    if (loc == "none"){
        
    }
}

proportion.plots <- function(best_models, trials=1:5, on_times, plot=TRUE, save_path=NA, 
                             state_names=c("exploratory", "encamped"), include_sun_times=TRUE){
    # Generates plots for each trial with the relative proportions of behavioral states.
    
    # make sure that all trials have data
    n_models <- length(best_models)
    for (i in 1:length(models)){
        trial_comparison <- as.numeric(unique(best_models[[1]]$data$Trial))
        if (setequal(trials, trial_comparison)){
            next
        }else{
            stop(paste0("Repetition ", i, " has fewer trials than required."))
        }
    }
    
    # make sure that on_times matches n_models
    if (n_models != on_times){stop("Number of fitted models does not match number of on_times.")}
    
    # compile the relevant dataframe
    df_proportions <- data.frame(matrix(ncol=4, nrow=0))
    colnames(df_proportions) <- c("Rep", "Trial", "States", "Names")
    for (rep in 1:n_models){
        best_model <- best_models[[rep]]
        state_sequence <- viterbi(best_model)
        
        holder <- data.frame(matrix(ncol=4, nrow=length(state_sequence)))
        colnames(holder) <- c("Rep", "Trial", "States", "Names")
        data <- df_proportions[df_proportions$Trial == trial,]
        
        holder["Rep"] <- paste0("Rep. ", rep)
        holder["Trial"] <- best_model$data$Trial
        holder["States"] <- state_sequence
        for (i in 1:length(state_names)){
            holder[holder$States == i, "Names"] <- state_names[i]
        }

        df_proportions <- rbind(df_proportions, holder)
    }

    # make the plots
    for (trial in trials){
        data <- df_proportions[df_proportions$Trial == trial,]
        
        ggplot(data, aes(x=Rep, y=States, fill=Names)) +
            geom_col(position="fill") +
            # possibly backwards
            scale_fill_manual(values=c("#E69F00", "#56B4E9"))
    }
}
















