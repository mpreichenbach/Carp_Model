library(tidyverse)
library(momentuHMM)


proportion.plots <- function(models, trials=1:5, plot=TRUE, save_path=NA, 
                             state_names=c("exploratory", "encamped"), 
                             state_colors=c("#56B4E9", "#E69F00")){
    # Generates plots for each trial with the relative proportions of behavioral states.
    
    # make sure that all trials have data
    n_models <- length(models)
    for (i in 1:length(models)){
        trial_comparison <- as.numeric(unique(models[[1]]$data$Trial))
        if (setequal(trials, trial_comparison)){
            next
        }else{
            stop(paste0("Repetition ", i, " has fewer trials than required."))
        }
    }

    # compile the relevant dataframe
    df_proportions <- data.frame(matrix(ncol=6, nrow=0))
    colnames(df_proportions) <- c("Rep", "Trial", "States", "Names", "Diel", "RepDiel")
    for (rep in 1:n_models){
        model <- models[[rep]]
        state_sequence <- viterbi(model)
        
        holder <- data.frame(matrix(ncol=6, nrow=length(state_sequence)))
        colnames(holder) <- c("Rep", "Trial", "States", "Names", "Diel")
        
        holder$Rep <- rep
        holder$Trial <- model$data$Trial
        holder$States <- state_sequence
        holder["Diel"] <- as.numeric(model$data$Diel)
        
        # for reps which a mix of day/night data, assign to night (1) if more than 50% are night
        holder$RepDiel <- as.factor(ifelse(sum(holder$Diel) / nrow(holder) > 0.5, 1, 0))
        for (i in 1:length(state_names)){
            holder[holder$States == i, "Names"] <- state_names[i]
        }

        df_proportions <- rbind(df_proportions, holder)
    }

    # make the plots
    for (trial in trials){
        data <- df_proportions[df_proportions$Trial == trial,]
        # get sunrise/sunset times
        
        plt <- ggplot(data, aes(x=Rep, y=States, fill=Names, alpha=factor(RepDiel))) +
            geom_col(position="fill") +
            scale_fill_manual(values=state_colors) +
            coord_cartesian(xlim=c(1, 24), ylim=c(0, 1)) +
            scale_x_continuous(breaks=seq(from=1, to=24, by=2)) +
            scale_alpha_manual(values=c("0"=0.5, "1"=1.0), guide="none") +
            labs(title=paste0("Trial ", trial, " Behavioral States"), 
                 x="Repetition Number", 
                 y="Proportion", 
                 fill="States") +
            theme(plot.title=element_text(hjust = 0.5), 
                  legend.title=element_blank(),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  panel.background=element_blank(),
                  axis.line=element_line(colour="black"))
                  
        print(plt)
        
        if (!(is.na(save_path))){
            ggsave(paste0("D:/Carp-Model/Trial ", trial, " State Proportions (5 min).png"))
        }
    }
}

aic.plot <- function(df, rep, colors=c("blue", "red"), disp_cov=NA){
    # makes a line graph rescaled so that the best model is at 0, and subsequent points increase by
    # their DeltaAIC value. The argument disp_cov should be one of the strings"dB", "Temp", "Pond", 
    # "Trial", "Treatment". The colors vector specifies first the main plot color, and secondarily 
    # the covariate point colors.
    
    # creates values to plot, and a factor variable for which covariates to highlight.
    df$Values <- df$AIC - min(df$AIC)
    df$DispCov <- paste0("Without ", disp_cov)
    if (!is.na(disp_cov)){
        for (covariate in disp_cov){
            df[grepl(covariate, df$formula, fixed=TRUE), "DispCov"] <- paste0("With ", disp_cov)
        }
        df$DispCov <- factor(df$DispCov)
    }

    # create the plot
    plt <- ggplot(data=df, aes(x=reorder(formula, Rank), y=Values)) + 
        geom_line(aes(group=1), size=1, color="blue") + 
        geom_point(aes(colour=DispCov), size=3) +
        scale_color_manual(values=colors) +
        scale_x_discrete(labels=df$formula) +
        labs(title=paste0("Cumulative \u0394AIC Scores for Rep. ", rep), color="Formulas") +
        theme(axis.text.x = element_text(angle=90, hjust=0.95, vjust=0.2),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              plot.title=element_text(hjust=0.5),
              legend.title=element_text(hjust=0.5))
    
    print(plt)
}
