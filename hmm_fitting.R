library(momentuHMM)
library(parallel)
source("model_functions.R")

fit_model <- function(.data, 
                      modelFormula,
                      factorCovs=c("Trial", "Pond", "Treatment", "Sound", "Diel"),
                      numericCovs=c("Temperature", "dB"),
                      stateNames=c("exploratory", "encamped"),
                      dist=list(step="gamma", angle="vm"),
                      initPar=list(step = c(2, 1, 1, 1), angle=c(0.004, 0.004, 0.002, 0.002))) {
    
    # the pipe functions in crw_fitting remove the momentuHMMData class; this puts it back in
    .data <- .data[, ! names(.data) %in% c("step", "angle")]
    .data <- prepData(.data)
    
    # this ensures that covariate columns have the correct numeric/factor types
    for (factor_name in factorCovs){
        .data[[factor_name]] <- as.factor(.data[[factor_name]])
    }
    
    for (numeric_name in numericCovs){
        .data[[numeric_name]] <- as.numeric(.data[[numeric_name]])
    }
    
    # ensure that modelFormula is actually of class "formula"
    modelFormula <- as.formula(modelFormula)

    # this logic adds zero-mass parameters, which are necessary when 0 is a step-length
    has_zero_step <- (0 %in% .data$step)
    if (has_zero_step){
        zero_proportion <- length(which(.data$step == 0)) / nrow(.data)
        for (state in 1:length(stateNames)){
            initPar$step <- append(initPar$step, zero_proportion)
        }
    }
    
    # sometimes the initPar values don't match the fitted parameters; the while loop checks for that
    incorrect_means <- TRUE
    
    while (incorrect_means){
        # this section fits an initial movement model to give better starting values when fitting the full HMM.
        model_1 <- fitHMM(data=.data,
                          nbStates=length(stateNames),
                          dist=dist,
                          Par0=initPar,
                          estAngleMean=list(angle=TRUE),
                          stateNames=stateNames)
        
        print("Finished fitting movement model (step 1/3).")
        
        # this fits a model to estimate good starting transition probabilities
        initPar1 <- getPar0(model=model_1, formula=modelFormula)
        
        model_2 <- fitHMM(data=model_1$data,
                          dist=dist,
                          nbStates=length(stateNames),
                          estAngleMean=list(angle = TRUE),
                          stateNames=stateNames,
                          Par0=initPar1$Par,
                          beta0=initPar1$beta,
                          formula=modelFormula)
        
        DM <- list(step=list(mean=modelFormula, sd=~1),
                   angle=list(mean=modelFormula, concentration=~1))
        
        if (has_zero_step){DM$step$zeromass <- ~1}
        
        print("Finished fitting model for estimating initial transition probabilities (step 2/3).")
        
        # this fits the full model
        
        initPar2 <- getPar0(model=model_2,
                            formula=modelFormula,
                            DM=DM)
        
        FullModel <- fitHMM(data=model_2$data,
                            nbStates=length(stateNames),
                            dist=dist,
                            Par0=initPar2$Par,
                            beta0=initPar2$beta,
                            DM=DM,
                            stateNames=stateNames,
                            estAngleMean=list(angle = TRUE),
                            formula=modelFormula)
        
        # Makes sure that the fitted models have mean steps which align with expectations
        state_1_mean <- FullModel$CIreal$step$est["mean", stateNames[1]]
        state_2_mean <- FullModel$CIreal$step$est["mean", stateNames[2]]
        
        incorrect_means <- state_1_mean < state_2_mean
    }
    
    print("Finished fitting full model (step 3/3).")
    
    return(FullModel)
}

fit_model_list <- function(list_element) {
    # this runs fit_model, but with a single element so that it can be entered as an argument in 
    # parallel::mclapply().
    
    hmm <- fit_model(list_element$data, 
                     modelFormula=list_element$formula,
                     stateNames=c("exploratory", "encamped"), 
                     dist=list(step="gamma", angle="vm"),
                     initPar=list(step=c(2, 1, 1, 1), angle=c(0.004, 0.004, 0.002, 0.002)))
    
    return(hmm)
}

hmm_parallel_fit <- function(data,
                          step_max = 30,
                          cluster_size = 10,
                          n_covs = c(0, 1, 2, 3, 4, 5, 6),
                          must_have_covs = NULL,
                          ignore_covs = NULL){
    # fits an HMM for many models at once, using the parallel::parLapplyLB() function
    
    # create a list to hold the lists of fitted model, indexed by n_covs
    model_holder <- list()
    
    # large step lengths messes up the fitting process; removes IDs which go above step_max
    bad_ids <- unique(data[data$step > step_max, "ID"])
    data <- data[! data$ID %in% bad_ids, ]
    
    for (n_cov in n_covs) {
        formulas <- get_formulas(n_cov, must_have_covs=must_have_covs, ignore_covs=ignore_covs)
        
        # return the models if there are no more formulas to fit
        if (length(formulas) == 0) {
            return(model_holder)
        }
        
        # generates a list of formula/data pairs to input to fit_model_list
        frm_list <- list()
        for (frm in formulas) {
            key <- deparse(frm)
            frm_list[[key]] <- list("formula" = frm, "data" = data)
        }
        
        # set up cluster for multiprocessing
        cl <- makeCluster(min(length(formulas), cluster_size))
        clusterExport(cl, c("fit_model_list", "fit_model"))
        clusterEvalQ(cl, library(momentuHMM))
        
        # fit all the HMMs (this can be a long process)
        tic = Sys.time()
        hmm <- parLapplyLB(cl, frm_list, fit_model_list)
        toc = Sys.time()
        
        print(paste0("Fitting models with ", n_cov, " covariates is complete;"))
        print(toc - tic)
        
        model_holder[[as.character(n_cov)]] <- hmm
        stopCluster(cl)
    }
    
    model_holder
}