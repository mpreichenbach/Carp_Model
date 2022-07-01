library(mometuHMM)


fit.model <- function(df, stateNames = c("exploratory", "encamped"), dist = list(step = "gamma", angle = "vm"),
                      initPar = list(step = c(2, 1, 1, 1), angle = c(0, 0, 0, 0)),
                      modelFormula = ~ Trial + Pond + Diel + Temp + dB * Treatment){
    
    # this ensures that covariate columns have the correct numeric/factor types; change the initPar
    # to better 
    
    df <- within(df, {
        Trial <- as.factor(Trial)
        Pond <- as.factor(Pond)
        Treatment <- as.factor(Treatment)
        Sound <- as.factor(Sound)
        Diel <- as.factor(Diel)
        Temp <- as.numeric(Temp)
        dB <- as.numeric(dB)
    })
    
    # this logic takes adds zeromass parameters, if 0 exists as a step-length in the data
    has_zero_step <- (0 %in% df$step)
    if (has_zero_step){
        zero_proportion <- length(which(df$step == 0)) / nrow(df)
        for (state in 1:length(stateNames)){
            initPar$step <- append(initPar$step, zero_proportion)
        }
    }
    
    # this section fits an initial movement model to give better starting values when fitting the full HMM.
    
    m1 <- fitHMM(data = df, 
                 nbStates = length(stateNames),
                 dist = dist,
                 Par0 = initPar,
                 estAngleMean = list(angle=TRUE),
                 stateNames = stateNames)
    
    print("Finished fitting movement model (step 1/3).")
    
    # this fits a model to estimate good starting transition probabilities
    initPar1 <- getPar0(model = m1, formula = modelFormula)
    
    m2 <- fitHMM(data = m1$data,
                 dist = dist,
                 nbStates = length(stateNames),
                 estAngleMean = list(angle = TRUE),
                 stateNames = stateNames,
                 Par0 = initPar1$Par,
                 beta0 = initPar1$beta,
                 formula = modelFormula)
    
    DM <- list(step = list(mean = modelFormula, sd = ~ 1),
               angle = list(mean = modelFormula, concentration = ~ 1))
    
    if (has_zero_step){DM$step$zeromass <- ~1}
    
    print("Finished fitting model for estimating initial transition probabilities (step 2/3).")
    
    # this fits the full model
    
    initPar2 <- getPar0(model = m2, formula = modelFormula, DM = DM)
    
    FullMod <- fitHMM(data = m2$data,
                      nbStates = length(stateNames),
                      dist = dist,
                      Par0 = initPar2$Par,
                      beta0 = initPar2$beta,
                      DM = DM,
                      stateNames = stateNames,
                      estAngleMean = list(angle = TRUE),
                      formula = modelFormula)
    
    print("Finished fitting full model (step 3/3).")
    
    return(FullMod)
}

fit.model.list <- function(list_element){
    # this runs fit.model, but with a single element so that it can be entered as an argument in 
    # parallel::mclapply().
    
    hmm <- fit.model(list_element$data, stateNames=c("exploratory", "encamped"), dist=list(step="gamma", angle="vm"),
                     initPar=list(step=c(2, 1, 2, 1), angle=c(0.004, 0.004, 0.002, 0.002)),
                     modelFormula=list_element$formula)
    
    return(hmm)
}