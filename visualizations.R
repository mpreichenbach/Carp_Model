library(grid)
library(gtable)
library(tidyverse)
library(momentuHMM)


proportion.plots <- function(models, rep_times, ba=NA, treatment=NA, trials=1:5, show_plot=TRUE, 
                             save_path=NA, state_names=c("exploratory", "encamped"), 
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
    df_proportions <- data.frame(matrix(ncol=6, nrow=0))
    colnames(df_proportions) <- c("Time", "Rep", "Trial", "States", "Names", "Diel")
    
    for (rep in 1:n_models){
        model <- models[[rep]]
        df <- model$data
        
        # calculate behavioral states if not included already
        if (!("States" %in% colnames(df))){
            df$States <- viterbi(model)
            print(paste0("Finished decoding state sequence for Repetition ", rep, "."))
        }
        
        # subset to the specified treatment
        if (!is.na(treatment)){df <- df[df$Treatment == treatment,]}
        
        holder <- data.frame(matrix(ncol=6, nrow=nrow(df)))
        colnames(holder) <- c("Time", "Rep", "Trial", "States", "Diel", "RepDiel")
        
        holder$Time <- df$Time
        holder$Rep <- rep
        holder$Trial <- df$Trial
        holder$States <- df$States
        holder$Diel <- as.numeric(as.character(df$Diel))
        holder$RepDiel <- as.factor(df$Diel)
        
        df_proportions <- rbind(df_proportions, holder)
    }

    # make the plots
    for (trial in trials){
        data <- df_proportions[df_proportions$Trial == trial,]
        data$BA <- "before"
        
        for (rep in 1:24){
            data0 <- data[data$Rep == rep,]
            # if there are both day/night values, assign to night
            if (0 %in% unique(data0[data0$Rep == rep, "Diel"])){
                data0[data0$Rep == rep, "RepDiel"] <- 0
            }
            
            # assign before/after labels relative to on_time
            on_time <- unique(as_hms(rep_times[rep_times$Repetition == rep, "Time"]))
            
            # note that hms inequalities assume the minimum time is 00:00:00; this logic deals with that
            if (on_time == hms(0, 0, 0)){
                data0[on_time <= as_hms(data0$Time) &
                          as_hms(data0$Time) <= hms(0, ba, 0), "BA"] <- "after"
            }else{
                data0[on_time <= as_hms(data0$Time), "BA"] <- "after"
            }
            data[data$Rep == rep,] <- data0
        } 
        
        data$RepDiel <- as.factor(data$RepDiel)
        
        # subset into before/after on_time
        data_before <- data[data$BA == "before",]
        data_after <- data[data$BA == "after",]
        
        # first plot, before on_time
        if (is.na(treatment)){
            before_title <- paste0("Trial ", trial, " Behavioral States (", ba, 
                                  "min before)")
        }else{
            before_title <- paste0("Trial ", trial, " ", treatment, " Behavioral States (", ba, 
                                  "min before)")
        }
        
        plt_before <- ggplot(data_before, aes(x=Rep, y=States, fill=factor(States))) + 
            geom_bar(aes(alpha=factor(RepDiel)), position="fill", stat="identity") +
            scale_fill_manual(values=state_colors, labels=state_names) +
            coord_cartesian(xlim=c(1, 24), ylim=c(0, 1)) +
            scale_x_continuous(breaks=seq(from=1, to=24, by=2)) +
            scale_alpha_manual(values=c("0"=0.5, "1"=1.0), guide="none") +
            labs(title=before_title, x="Repetition Number", y="Proportion") +
            theme(plot.title=element_text(hjust = 0.5), 
              legend.title=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.background=element_blank(),
              axis.line=element_line(colour="black"))
        
        # second plot, after on_time
        if (is.na(treatment)){
            after_title <- paste0("Trial ", trial, " Behavioral States (", ba, 
                                  "min after)")
        }else{
            after_title <- paste0("Trial ", trial, " ", treatment, " Behavioral States (", ba, 
                                  "min after)")
        }
        
        plt_after <- ggplot(data_after, aes(x=Rep, y=States, fill=factor(States))) + 
            geom_bar(aes(alpha=factor(RepDiel)), position="fill", stat="identity") +
            scale_fill_manual(values=state_colors, labels=state_names) +
            coord_cartesian(xlim=c(1, 24), ylim=c(0, 1)) +
            scale_x_continuous(breaks=seq(from=1, to=24, by=2)) +
            scale_alpha_manual(values=c("0"=0.5, "1"=1.0), guide="none") +
            labs(title=after_title, x="Repetition Number", y="Proportion") +
            theme(plot.title=element_text(hjust = 0.5), 
                  legend.title=element_blank(),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  panel.background=element_blank(),
                  axis.line=element_line(colour="black"))
        
        # bring the plots together into one figure (or gtable)
        g_before <- ggplotGrob(plt_before)
        g_after <- ggplotGrob(plt_after)
        g <- rbind(g_before, g_after, size="first")
        g$widths <- unit.pmax(g_before$widths, g_after$widths)

        # possibly print the gtable
        if (show_plot){
            grid.newpage()
            grid.draw(g)
        }
        
        # possibly save the object g
        if (!is.na(save_path)){
            if (is.na(treatment)){
                ptitle <- paste0("Trial ", trial, " State Proportions (", ba, "min).png")
            }else{
                ptitle <- paste0(treatment, " ", "Trial ", trial, " State Proportions (",
                                 ba, "min).png")
            }
            ggsave(file.path(save_path, ptitle), plot=g)
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


db.histograms <- function(models, ba=NA, bins=NA, binwidth=1, include_means=TRUE, 
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
        
        plt_title <- paste0("Rep ", rep, ": dB Histogram at Fish Positions")
        
        if (include_means){
            mean_1 <- mean(df[df$States == 1, "dB"], na.rm=TRUE)
            mean_2 <- mean(df[df$States == 2, "dB"], na.rm=TRUE)
            plt_title <-paste0(plt_title, " (with mean lines)")
        }
        
        plt <- ggplot(data=df, aes(x=dB, fill=Names)) + 
            geom_histogram(aes(y=..density..), binwidth=binwidth, alpha=0.5, position="identity") + 
            scale_fill_manual(values=state_colors) + 
            {if(include_means)geom_vline(xintercept=mean_1, size=1, color=state_colors[1])} +
            {if(include_means)geom_vline(xintercept=mean_2, size=1, color=state_colors[2])} +
            xlim(xlims) + 
            labs(title=plt_title, x="Sound Intensity (dB)") + 
            theme(plot.title=element_text(hjust=0.5),
                  legend.title=element_blank(),
                  panel.background=element_blank())
        
        if (show_plot){
            print(plt)
        }     
        
        if (!is.na(save_path)){
            ggsave(paste0(save_path, "Rep ", rep, " dB Histogram by State (", ba, "min).png"))
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
        labs(title=paste0("Difference between Exploratory and Encamped Mean dB (", ba, "min)"),
             x="Repetition") +
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


plot_predicted_means <- function(df,
                                 x_col = "dB",
                                 y_cols = c("step_lower", "step_est", "step_upper"),
                                 y_label = "Mean Step Estimate (meters)",
                                 plot_title = "",
                                 state_names = c("exploratory", "encamped"),
                                 state_colors = c("#E69F00", "#56B4E9")) {
    # makes a plot of the state means with shaded confidence intervals
    
    # if the step means don't correspond to the exploratory/encamped intuition, then switch labels
    exp_mean_step <- mean(df[df$State == state_names[1], c("step_est")])
    enc_mean_step <- mean(df[df$State == state_names[2], c("step_est")])
    
    if (exp_mean_step < enc_mean_step) {
        df$State <- ifelse(df$State == state_names[1], state_names[2], state_names[1])
    }
    
    # ensure the state factors are in a consistent order
    df$State <- factor(df$State, levels = state_names)
    
    # generate the plot
    plt <- ggplot(df, aes_(x = as.name(x_col), y = as.name(y_cols[2]), col = as.name("State"))) + 
        geom_line() + 
        geom_ribbon(aes_(ymin = as.name(y_cols[1]), ymax = as.name(y_cols[3]), fill = as.name("State")),
                    alpha = 0.2,
                    linetype="dotted") +
        scale_fill_manual(values = state_colors) + 
        labs(x = x_col, y = y_label, title = plot_title) + 
        ylim(0, 1.2 * max(df[[y_cols[2]]])) + 
        theme(plot.title = element_text(hjust = 0.5))
    
    plt
}






























