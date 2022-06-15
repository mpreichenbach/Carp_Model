library(tidyverse)
library(momentuHMM)


proportion.plots <- function(models, ba="5min", trials=1:5, show_plot=TRUE, save_path=NA, 
                             state_names=c("exploratory", "encamped"), 
                             state_colors=c("#E69F00", "#56B4E9")){
    # Generates plots for each trial with the relative proportions of behavioral states; the
    # argument "ba" is the time period considered before/after the sound.
    
    if (is.na(ba)){stop("Specify the before/after period in 'ba' argument.")}
    
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
    df_proportions <- data.frame(matrix(ncol=5, nrow=0))
    colnames(df_proportions) <- c("Rep", "Trial", "States", "Names", "Diel")
    
    for (rep in 1:n_models){
        model <- models[[rep]]
        holder <- data.frame(matrix(ncol=5, nrow=nrow(model$data)))
        colnames(holder) <- c("Rep", "Trial", "States", "Diel", "RepDiel")
        
        holder$Rep <- rep
        holder$Trial <- model$data$Trial
        holder$Diel <- as.numeric(as.character(model$data$Diel))
        holder$RepDiel <- as.factor(model$data$Diel)
        
        if ("States" %in% colnames(model$data)){
            holder$States <- model$data$States
        }else{
            holder$States <- viterbi(model)
            print(paste0("Finished decoding state sequence for Repetition ", rep, "."))
        }
        
        df_proportions <- rbind(df_proportions, holder)
    }

    # make the plots
    for (trial in trials){
        data <- df_proportions[df_proportions$Trial == trial,]
        
        # if there are both day/night 
        for (rep in 1:24){
            if (0 %in% unique(data[data$Rep == rep, "Diel"])){
                data[data$Rep == rep, "RepDiel"] <- 0
            }
        } 
        data$RepDiel <- as.factor(data$RepDiel)

        plt <- ggplot(data, aes(x=Rep, y=States, fill=factor(States))) + 
            geom_bar(aes(alpha=factor(RepDiel)), position="fill", stat="identity") +
            scale_fill_manual(values=state_colors, labels=state_names) +
            coord_cartesian(xlim=c(1, 24), ylim=c(0, 1)) +
            scale_x_continuous(breaks=seq(from=1, to=24, by=2)) +
            scale_alpha_manual(values=c("0"=0.5, "1"=1.0), guide="none") +
            labs(title=paste0("Trial ", trial, " Behavioral States (", ba, ")"), x="Repetition Number",
                 y="Proportion") +
            theme(plot.title=element_text(hjust = 0.5), 
              legend.title=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.background=element_blank(),
              axis.line=element_line(colour="black"))
        
        if (show_plot){
            print(plt)
        }     
        
        if (!is.na(save_path)){
            ggsave(paste0(save_path, "Trial ", trial, " State Proportions (", ba, ").png"))
        }
    }
}

aic.plot <- function(df, rep, include_ranks=1:5, colors=c("red", "black"), disp_cov=NA){
    # makes a line graph rescaled so that the best model is at 0, and subsequent points increase by
    # their DeltaAIC value. The argument disp_cov should be one of the strings"dB", "Temp", "Pond", 
    # "Trial", "Treatment". The colors vector specifies first the main plot color, and secondarily 
    # the covariate point colors.
    
    # creates values to plot, and a factor variable for which covariates to highlight.
    df <- df[df$Rank %in% include_ranks,]
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
        geom_line(aes(group=1), size=1) + 
        geom_point(aes(colour=DispCov), size=3, show.legend=!is.na(disp_cov)) +
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

db.histogram <- function(models, ba=NA, bins=NA, binwidth=1, include_means=TRUE, 
                         show_plot=TRUE, save_path=NA, state_names=c("exploratory", "encamped"), 
                         state_colors=c("#E69F00", "#56B4E9")){
    # plots overlapping histograms of the dB variable, colored by the different states
    
    if (is.na(ba)){stop("Specify the before/after period in 'ba' argument.")}
    
    if (sum(c(is.na(bins), is.na(binwidth))) != 1){
        stop("Only assign a value to one of bins/binwidth.")
    }
    
    n_reps <- length(models)
    min_db <- Inf
    max_db <- -Inf
    
    # find the min/max dB among all reps to use as consistent x-axis limits of the histograms
    for (rep in 1:n_reps){
        df <- models[[rep]]$data
        df <- df[df$dB > 0,]
        df_min_db <- min(df$dB, na.rm=TRUE)
        df_max_db <- max(df$dB, na.rm=TRUE)
        if (df_min_db < min_db){min_db <- df_min_db}
        if (df_max_db > max_db){max_db <- df_max_db}
    }
    
    xlims <- c(min_db, max_db)
    
    # plot the histograms
    for (rep in 1:n_reps){
        df <- models[[rep]]$data
        
        df <- df[df$dB > 0,]
        df$Names <- "exploratory"
        df[df$States == 2, "Names"] <- "encamped"
        df$Names <- factor(df$Names, levels=c("exploratory", "encamped"))
        
        plt_title <- paste0("Repetition ", rep, ": Sound Intensity Histogram at Fish Positions")
        
        if (include_means){
            mean_1 <- mean(df[df$States == 1, "dB"], na.rm=TRUE)
            mean_2 <- mean(df[df$States == 2, "dB"], na.rm=TRUE)
            plt_title <-paste0(plt_title, " (with mean lines)")
        }
        
        plt <- ggplot(data=df, aes(x=dB, fill=Names)) + 
            geom_histogram(aes(y=..density..), binwidth=binwidth, alpha=0.5, position="identity") + 
            scale_fill_manual(values=state_colors) + 
            {if(include_means)geom_vline(xintercept=mean_1, size=1.5, color=state_colors[1])} +
            {if(include_means)geom_vline(xintercept=mean_2, size=1.5, color=state_colors[2])} +
            xlim(xlims) + 
            labs(title=plt_title, x="Sound Intensity (dB)") + 
            theme(plot.title=element_text(hjust=0.5),
                  legend.title=element_blank(),
                  panel.background=element_blank())
        
        if (show_plot){
            print(plt)
        }     
        
        if (!is.na(save_path)){
            ggsave(paste0(save_path, "Repetition ", rep, " dB Histogram by State (", ba, ").png"))
        }
    }
}

db.means.plot <- function(models, ba=NA, state_colors=c("#E69F00", "#56B4E9"), show_plot=TRUE, 
                          save_path=NA){
    # plots the means of dB levels experienced in different states
    
    if (is.na(ba)){stop("Specify the before/after period in 'ba' argument.")}
    
    state_1_means <- c()
    state_2_means <- c()
    
    for (i in 1:length(models)){
        df <- models[[i]]$data
        state_1_means[i] <- mean(df[df$dB > 0 & df$States == 1, "dB"], na.rm=TRUE)
        state_2_means[i] <- mean(df[df$dB > 0 & df$States == 2, "dB"], na.rm=TRUE)
    }
    
    diff_vec <- state_1_means - state_2_means
    
    df_means <- data.frame(diff_vec, sign(diff_vec))
    colnames(df_means) <- c("Difference", "Sign")
    
    df_means$Names <- "encamped"
    df_means[df_means$Sign > 0, "Names"] <- "exploratory"
    df_means$Sign <- factor(df_means$Sign)
    df_means$States <- factor(df_means$Names, levels=c("exploratory", "encamped"))
    
    # makes the line plots
    
    plt <- ggplot(df_means, aes(x=1:length(models), y=Difference)) + 
        geom_col(aes(fill=States)) +
        scale_fill_manual(values=state_colors) +
        scale_x_continuous(breaks=seq(from=1, to=24, by=2)) +
        labs(title="Difference between Exploratory and Encamped Mean dB", x="Repetition") +
        theme(plot.title=element_text(hjust=0.5),
              axis.title.y=element_blank(),
              legend.title=element_blank(),
              panel.background=element_blank())
    
    if (show_plot){
        print(plt)
    }     
    
    if (!is.na(save_path)){
        ggsave(paste0(save_path, "Difference Between Mean dB by State.png"))
    }
}
    