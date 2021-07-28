# run_models.r
# Created by John Trochta
# Date created:  11/06/2019
# Summary:
# This is a bulky function that runs BASA models in parallel.
# There are various options to select, from which covariates to include and how, 
# different model structures, and different algorithms

run_models <- function(hypothesis_names,
                       include.null=TRUE,
                       start.years,
                       H.select,
                       n_cores,
                       n_sims,
                       n_age0_covs,
                       n_mor_covs,
                       custom.settings,
                       mcmc.algorithm="NUTS",
                       nuts.inputs=data.frame(reps=3,iter=1000,adapt_delta=0.95,warm_up=100,nuts.mode=1),
                       run_year="2019",
                       run_month="11"){
  
  ### Setup directory for storing ADMB input & output files
  n_run_types <- custom.settings$n_run_types
  R.or.M <- custom.settings$R.or.M
  sim_cov <- custom.settings$sim_cov
  series_start_yr <- custom.settings$series_start_yr
  estimate.1993.mortality <- custom.settings$estimate.1993.mortality
  estimate.mean.mortality <- custom.settings$estimate.mean.mortality
  mean.mortality <- custom.settings$mean.mortality
  h.block <- custom.settings$h.block
  block.on <- custom.settings$block.on
  run_ID <- custom.settings$run_ID
  base.M <- custom.settings$base.M
  
  hypothesis.1 <- gsub(pattern=" |-|/",replacement="_",x=hypothesis_names)
  hypothesis.1 <- gsub(pattern="[()']",replacement="",x=hypothesis.1)
  hypothesis_names.1 <- hypothesis_names
  
  modelPath <- paste0(here::here("results"),"/",run_ID,"/")
  create_modelPath <- 1
  if(create_modelPath==1){
    dir.create(modelPath)
    if(mcmc.algorithm=="MH"){
      folder_1 <- paste0(modelPath,"MCMC_diagnostics/")
      folder_2 <- paste0(modelPath,"Posterior_distributions/")
      folder_3 <- paste0(modelPath,"Parameter_correlations/")
      
      dir.create(folder_1)
      dir.create(folder_2)
      dir.create(folder_3)
    }else if(mcmc.algorithm=="NUTS"){
      folder_1 <- paste0(modelPath,"NUTS_chain_diagnostics/")
      folder_2 <- paste0(modelPath,"NUTS_RDS_files/")
      
      dir.create(folder_1)
      dir.create(folder_2)
    }
  }
  
  filename <- paste0(here::here("bin"),"/",base.M,"/PWS_ASA(covariate).ctl")
  phase.file <- paste0(here::here("bin"),"/",base.M,"/PWS_ASA(phases).ctl")
  pin.file <- paste0(here::here("bin"),"/",base.M,"/PWS_ASA.PIN")
  
  ### Allows flexibility to select only fitting time series that start from specified year
  turn_on_hypotheses.1 <- which(start.years<=(series_start_yr+1)) # Plus one to account for series starting one year later & thus still valid
  if(!is.null(H.select)){
  turn_on_hypotheses <- turn_on_hypotheses.1[H.select$single]
  }else{
    turn_on_hypotheses <- turn_on_hypotheses.1
  }
  
  hypothesis <- hypothesis.1[turn_on_hypotheses]
  hypothesis_names <- hypothesis_names.1[turn_on_hypotheses]
  #h.block <- h.block[h.block%in%turn_on_hypotheses]
  
  if(include.null){
    hypothesis <- c(hypothesis,"Null_model")#c("Null_model",hypothesis[turn_on_hypotheses])
    hypothesis_names <- c(hypothesis_names,"Null model")
    #h.block <- h.block+1 # Offset by 1 to account for Null_model not existing in covariate matrix
  }
  
  if(length(H.select)>1){
    extra_hypothesis <- names(H.select)[names(H.select)!='single']
    extra_hypothesis_names <- extra_hypothesis
    
    extra_hypothesis <- gsub(pattern=" |-|/",replacement="_",x=extra_hypothesis)
    extra_hypothesis <- gsub(pattern="[()']",replacement="",x=extra_hypothesis)
    
    hypothesis <- c(hypothesis,extra_hypothesis)
    hypothesis_names <- c(hypothesis_names,extra_hypothesis_names)
  }
  
  ### Compile code in base folder before copying over so all I have to do is run the executable
  template_files <- paste0(here::here("bin"),"/",base.M,"/")  
  setwd(template_files)
  system('admb PWS_ASA')
  
  if(mcmc.algorithm=="NUTS"){
    summary_diagnost <- paste0(folder_1,'model_convergence_summary.csv')
    sum_dia <- data.frame('model',
                          'divergences.from.extract.function',
                          'divergences.from.fit.list',
                          'min.ESS',
                          'which.min.ESS',
                          'max.Rhat',
                          'which.max.Rha',
                          'time.elapsed',
                          'nonconvergence',
                          'iter.per.chain')
    write.table(sum_dia,file=summary_diagnost,sep=",",col.names=FALSE,row.names=FALSE)
  }
  
  #################################################################
  ### Run models in parallel
  library(doParallel)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  ptm <- proc.time()
  #BASA <- foreach(i=10:12) %dopar%{
  BASA <- foreach(i=1:length(hypothesis),.combine='rbind') %dopar%{
    
    if(R.or.M[i]=="M"){
      covar.matrix <- 10
      covar.matrix.2 <- 14 # For turning on certain years to estimate covariate - a time block
      covar.matrix.3 <- 12 # Actual matrix
    }else if(R.or.M[i]=="R"){
      covar.matrix <- 4
      covar.matrix.2 <- 6
      covar.matrix.3 <- 5 # Actual matrix
    }
    
    library(dplyr)
    library(tidyr)
    library(prodlim)
    library(R2admb)
    source(file=paste0(here::here("src"),"/data_reader.R"))
    source(file=paste0(here::here("src"),"/data_header_reader.R"))
    source(file=paste0(here::here("src"),"/color_ts.R"))
    source(file=paste0(here::here("src"),"/calc_beta.R"))
    source(file=paste0(here::here("src"),"/calc_beta_fast.R"))
    
    setwd(template_files)
    PWS_ASA.dat <- data_reader(filename="PWS_ASA.dat")
    nyr <- PWS_ASA.dat[[1]]
    nyr_tobefit <- PWS_ASA.dat[[2]]
    nyr_fit <- nyr_tobefit
    nage <- PWS_ASA.dat[[3]]
    
    covariates <- data_reader(filename="PWS_ASA(covariate).ctl")
    headers <- data_header_reader(filename="PWS_ASA(covariate).ctl")
    
    # Create a new directory for the hypothesis
    new.directory <- paste0(modelPath,hypothesis[i])
    
    dir.create(new.directory)
    # Change to that working directory
    setwd(new.directory)
    targetdir <- new.directory
    
    # Copy the files that will remain unchanged to this new directory
    filestocopy <- list.files(template_files)
    filestocopy <- paste0(rep(template_files,times=length(filestocopy)),filestocopy)
    file.copy(from=filestocopy, to=targetdir)
    
    # Turn on/off phases for two additional mortality parameters added to baseline mortality in 1993/1994 (parameterized after Melissa's model)
    phase.values <- data_reader(filename=phase.file)
    phase.headers <- data_header_reader(filename=phase.file)
    
    pin.values <- data_reader(filename=pin.file)
    pin.headers <- data_header_reader(filename=pin.file)
      
    if(estimate.1993.mortality[i]==1){
      phase.values[[21]] <- 3
      phase.values[[22]] <- 3
    }else if(estimate.1993.mortality[i]==0){
      phase.values[[21]] <- -3
      phase.values[[22]] <- -3
    }
    if(estimate.mean.mortality[i]==1){
      phase.values[[2]] <- 1
    }else if(estimate.mean.mortality[i]==0){
      phase.values[[2]] <- -1
    } 

    if(!is.null(block.on)){
      if(R.or.M[i]=="R"){
        phase.values[[26]] <- -1
        phase.values[[29]] <- 1
      }else if(R.or.M[i]=="M"){
        phase.values[[27]] <- -2
        phase.values[[30]] <- 2
      }
    }
    write.table(rbind(phase.headers[[1]][1],phase.headers[[1]][2],phase.headers[[1]][3],phase.values[[1]]," "),
                file = "PWS_ASA(phases).ctl", append = F, sep = " ",
                row.names=FALSE,col.names=FALSE,quote=F)
    for(j in 2:length(phase.values)) {write.table(phase.headers[[1]][j*2+1],
                                                  file = "PWS_ASA(phases).ctl", append = T, sep = " ",
                                                  row.names=FALSE,col.names=FALSE,quote=F)
      write.table(rbind(phase.values[[j]]," "),
                  file = "PWS_ASA(phases).ctl", append = T, sep = " ",
                  row.names=FALSE,col.names=FALSE,quote=F)
    }
      
    # Set baseline mortality and rewrite PWS_ASA.PIN file with it
    pin.values[[3]] <- mean.mortality[i]
    write.table(pin.values[[1]],
                file = "PWS_ASA.PIN", append = F, sep = " ",
                row.names=FALSE,col.names=FALSE,quote=F)
    for(j in 2:length(pin.values)) {
      write.table(rbind(pin.values[[j]]),
                  file = "PWS_ASA.PIN", append = T, sep = " ",
                  row.names=FALSE,col.names=FALSE,quote=F)
    }
    
    # Check to see which covariate(s) is/are being modeled and setup the proper indices within the PWS_ASA(covariate).ctl
    H.select.2 <- if(hypothesis[i]%in%names(H.select)) H.select[[hypothesis_names[i]]] else H.select$single[i]
    
    if(hypothesis[i]!="Null_model"){
      if(is.null(H.select)){
        hypothesis_indexing_foreach <- which(h.block==i)
      }else{
        hypothesis_indexing_foreach <- which(h.block%in%H.select.2)
      } 
      
      # Create .ctl file with the hypothesis to be tested
      covariates[[covar.matrix]][hypothesis_indexing_foreach] <- 1 # THIS ASSUMES DISEASE HYPOTHESIS IS IN FIRST 3 COLUMNS 
      
      if(!is.null(block.on[[i]])){
        covariates[[covar.matrix.2]][block.on-1979] <- 1
      }
      
      write.table(rbind(headers[[1]][1],headers[[1]][2],headers[[1]][3],covariates[[1]]," "),
                  file = "PWS_ASA(covariate).ctl", append = F, sep = " ",
                  row.names=FALSE,col.names=FALSE,quote=F)
      for(j in 2:length(covariates)) {write.table(headers[[1]][j*2+1],
                                                  file = "PWS_ASA(covariate).ctl", append = T, sep = " ",
                                                  row.names=FALSE,col.names=FALSE,quote=F)
        write.table(rbind(covariates[[j]]," "),
                    file = "PWS_ASA(covariate).ctl", append = T, sep = " ",
                    row.names=FALSE,col.names=FALSE,quote=F)
      }
    }
    # Modify the file with the covariates to turn on covariate to fit in model
    
    ######################################################
    # Calculate ESS for particular hypothesis
    
    # Read in measured age comps
    seine.ac <- PWS_ASA.dat[[20]][1:nyr_fit,]
    spawn.ac <- PWS_ASA.dat[[21]][1:nyr_fit,]
    
    PWS_ASA_ESS.ctl <- data_reader(filename="PWS_ASA(ESS).ctl")
    ESS_Se <- PWS_ASA_ESS.ctl[[1]][1:nyr_fit]
    ESS_Sp <- PWS_ASA_ESS.ctl[[2]][1:nyr_fit]
    
    # Read in the actual sample sizes
    seine.SS <- read.table("agecomp_samp_sizes.txt",header=FALSE,skip=4,nrows=nyr)[1:nyr_fit,1]
    spawn.SS <- read.table("agecomp_samp_sizes.txt",header=FALSE,skip=4+nyr+1,nrows=nyr)[1:nyr_fit,1]
    
    # Create empty matrices to fill estimated ESS and age comps
    #Seine
    Seine.ess.its <- matrix(0, nyr_fit, 1) # Matrix to hold all iterations of the routine
    Seine.ess.its <- seine.SS # fill in the first column with the recorded sample size 
    #Spawn
    Spawn.ess.its <- matrix(0, nyr_fit, 1)
    Spawn.ess.its <- spawn.SS
    
    # Change phases of the ESS in phases file to use the PWS_ASA(ESS_estimate)
    ph <- readLines("PWS_ASA(phases).ctl",-1)
    ph[4] <- 1
    writeLines(ph,"PWS_ASA(phases).ctl")
    
    # LOOP THROUGH AND ITERATIVELY CALCULATE ESS
    convergence <- 0
    
    seine.ESS <- seine.SS
    spawn.ESS <- spawn.SS
    # seine.ESS <- data_reader(filename="PWS_ASA(ESS_estimate).ctl")[[1]][1:nyr_fit]
    # spawn.ESS <- data_reader(filename="PWS_ASA(ESS_estimate).ctl")[[2]][1:nyr_fit]
    its <- 1
    
    for(j in 1:2){
      # Create "PWS_ASA(ESS_estimate).ctl" with sample sizes (the original sample sizes on the first iteration)
      write.table(rbind("# PWS age comp effective sample sizes","# Seine ESS", as.matrix(seine.ESS),
                        " ", "# Spawn ESS", as.matrix(spawn.ESS)),
                  file = "PWS_ASA(ESS_estimate).ctl", append = F, sep = " ",
                  row.names=FALSE,col.names=FALSE,quote=F)
      
      # Compile and Run PWS_ASA
      # shell("PWS_ASA -pinwrite") # For running on Windows
      system('./PWS_ASA -pinwrite') # For running on MAC
      
      clean_admb("PWS_ASA")
      
      # Read in the estimated seine and spawner age comps
      seine.ac.est <- read.table("SeAC_pd.rep", header = FALSE)[,-(1:3)]
      spawn.ac.est <- read.table("SpAC_pd.rep", header = FALSE)[,-(1:3)]
      
      # Calculate the ESS
      seine.ESS <- rowSums(seine.ac.est*(1-seine.ac.est))/rowSums((seine.ac[,-(1:3)]-seine.ac.est)^2)
      spawn.ESS <- rowSums(spawn.ac.est*(1-spawn.ac.est))/rowSums((spawn.ac[,-(1:3)]-spawn.ac.est)^2)
      
      # Remove the missing years of age comps
      seine.ESS.rem <- seine.ESS[!(seine.ac[,1]==-9)]
      spawn.ESS.rem <- spawn.ESS[!(spawn.ac[,1]==-9)]
      
      # Calculate the ratio of ESS to original sample sizes
      seine.ratio <- seine.ESS.rem/seine.SS[seine.ac[,1]!=-9]
      spawn.ratio <- spawn.ESS.rem/spawn.SS[spawn.ac[,1]!=-9]
      
      # Calculate the harmonic means
      seine.hm <- 1/mean(1/seine.ratio)
      spawn.hm <- 1/mean(1/spawn.ratio)
      
      # Compare this harmonic mean to the previous using a convergence criteria (WHAT AM I CONVERGING!!!!)
      if(its==1) {
        seine.hmS <- seine.hm
        spawn.hmS <- spawn.hm
      } else{
        seine.test <- abs(seine.hm - seine.hmS[its-1])/seine.hmS[its-1]*100
        spawn.test <- abs(spawn.hm - spawn.hmS[its-1])/spawn.hmS[its-1]*100
        seine.hmS <- rbind(seine.hmS,seine.hm)
        spawn.hmS <- rbind(spawn.hmS,spawn.hm)      
      }
      
      # Now multiply the harmonic mean by the sample size to get the new ESS 
      seine.ESS <- seine.hm*seine.SS
      spawn.ESS <- spawn.hm*spawn.SS
      
      seine.ESS[seine.ESS>0] <- round(seine.ESS[seine.ESS>0], digits=0)
      spawn.ESS[spawn.ESS>0] <- round(spawn.ESS[spawn.ESS>0], digits=0)
      
      # seine.ESS[seine.ESS>0] <- round(mean(seine.ESS[seine.ESS>0]), digits=0)
      # spawn.ESS[spawn.ESS>0] <- round(mean(spawn.ESS[spawn.ESS>0]), digits=0)
      # Denote the missing values
      seine.ESS[(seine.ac[,1]==-9)] <- -9
      spawn.ESS[(spawn.ac[,1]==-9)] <- -9
      # Fill in this iteration's ESS
      Seine.ess.its <- cbind(Seine.ess.its,round(seine.ESS,0))
      Spawn.ess.its <- cbind(Spawn.ess.its,round(spawn.ESS,0))
      
      its <- its+1
    }
    
    # Turn of the phases for the ESS calculation so it no longer recalculates ESS in future model runs
    ph <- readLines("PWS_ASA(phases).ctl",-1)
    ph[4] <- -1
    writeLines(ph,"PWS_ASA(phases).ctl")
    
    try.this <- data.frame(seine.ESS)
    
    # Now write the converged ESS to a ctl file to be used for model runs
    write.table(rbind("# PWS age comp effective sample sizes",paste0("# (",date(),")")," ",
                      "# Seine ESS", as.matrix(seine.ESS),
                      " ", "# Spawn ESS", as.matrix(spawn.ESS)),
                file = "PWS_ASA(ESS).ctl", append = F, sep = " ",
                row.names=FALSE,col.names=FALSE,quote=F)
    
    ######################################################
    ### Use Metropolis-Hastings or NUTS to sample posterior
    ### Jul 2021: MH was used much earlier in study, and has not been checked
    ### since
    if(mcmc.algorithm=="MH"){
      
      N_chain <- 7000000
      N_thin <- 3500
      N_burn <- 200 # The number of initial draws to remove from the MCMC chain
      # NEXT STEPS:  Change phases? (e.g. do not estimate natural mortality or average recruitment OR use mortality deviates)
      # Compile and Run PWS_ASA
      #shell('admb PWS_ASA')
      shell(paste0('PWS_ASA -mcmc ',N_chain,' -mcsave ',N_thin))
      shell('PWS_ASA -mceval')
      #system('export ADMB_HOME="/Users/johnt23/ADMBTerminal.app/admb"; export PATH="/Users/johnt23/ADMBTerminal.app/admb/bin":${PATH}; admb -s PWS_ASA; ./PWS_ASA') 
      
      # PLOTTING ===========================================================================
      variable_names <- c(  'VHSV_age3_4_mort_93',
                            'ICH_age5_8_mort_93',
                            'Z_0_8',   
                            'Z_9',
                            'Z_0_8offset',
                            'Z_9offset',
                            'matur_age3_per1',
                            'matur_age4_per1',
                            'matur_age3_per2',
                            'matur_age4_per2',
                            'alpha_maturity_prespawn',
                            'beta_maturity_prespawn',
                            'alpha_v', 
                            'beta_v',
                            'survey_vul_alpha',
                            'survey_vul_beta',
                            paste0('loginit_pop',1:5),
                            'egg_add',
                            'logmdm_c',
                            'm_add',
                            'hydADFG_q',
                            'hydADFG_add',
                            'hydPWSSC_q',
                            'hydPWSSC_add',
                            paste0('annual_age0devs',1:nyr_tobefit),
                            'log_MeanAge0',
                            'Mean_Age0offset',
                            'sigma_age0devs',
                            paste0('beta_age0',1:n_age0_covs),
                            paste0('beta_mortality',1:n_mor_covs),
                            paste0('annual_mortdevs',1:nyr_tobefit),
                            'sigma_mortdevs',
                            paste0('beta_age0_offset',1:n_age0_covs),
                            paste0('beta_mortality_offset',1:n_mor_covs),
                            'sigma_age0covar',
                            'sigma_morcovar')

      posterior.plots <- function(Outs){
        par.names <-colnames(Outs)
        Icnt <- length(Outs[,1]) # Rows are the draws-the length of the chain
        Ndim <- length(Outs[1,]) # Cols are the estimated parameter plus likelihood
        xx <- seq(1:Icnt)
        
        # Now plot parameters
        mat <- matrix(1:20,
                      nrow=5,ncol=4,byrow=T)
        par(oma=c(0,0,0,1),mar=c(2,1,2,1))
        layout(mat=mat)
        for (II in 1:Ndim) {
          yy <- Outs[,II][1:Icnt]
          lab1 <- par.names[II]
          plot(density(yy),main=lab1,yaxt="n",pch=16,cex=0.02)
          segments(x0=quantile(yy,probs=0.5), x1=quantile(yy,probs=0.5),
                   y0=0, y1=max(density(yy)$y), 
                   lty=1, lwd=3)
        }
      }
      library(coda)
      library(ggplot2)
      library(ggmcmc)
      
      #for(i in 1:length(hypothesis)){
      # setwd(paste0(modelPath,hypothesis[i]))
      targetdir <- targetdir
      Outs<-read.table(paste0(targetdir,"/iterations.csv"), header = FALSE, sep = ",", dec=".")
      names(Outs) <- c(variable_names,"objective_fn_val")
      Outs<-Outs[-seq(1,N_burn,by=1),] # Removing fixed parameters
      
      MCMC1 <- mcmc(Outs)
      ggmcmc(ggs(MCMC1),file=paste0(folder_1,hypothesis[i],"_",N_chain,"_diagnostics.pdf"),
             plot=c("ggs_traceplot()",
                    #"ggs_running()",
                    "ggs_compare_partial()",
                    "ggs_autocorrelation()",
                    "ggs_Rhat()"))
      
      pdf(file=paste0(folder_2,hypothesis[i],"_parameter_posterior_distributions.pdf"))
      posterior.plots(Outs)
      dev.off()
      
      library(corrplot)
      cor.struc <- cor(Outs,use="pairwise.complete.obs")
      trial <- cor.struc[apply(cor.struc,1,function(x) all(is.na(x)))==FALSE,
                         apply(cor.struc,2,function(x) all(is.na(x)))==FALSE]
      pdf(file=paste0(folder_3,hypothesis[i],"_parameter_correlations.pdf"))
      corrplot(trial,method="color",type="upper",
               tl.cex=0.5,tl.col="black",diag=FALSE)
      dev.off()
      #}
    }else if(mcmc.algorithm=="NUTS"){
      library(adnuts)
      library(snowfall)
      library(rstan)
      library(r4ss)
      library(shinystan)
      library(flock)

      # Tips
      # 1) If divergent transitions>0, increase target acceptance rate to reduce step size
      #    e.g. control = list(adapt_delta = 0.9)
      # 2) IF extreme global correlations, pass dense matrix estimated from previous run
      #    e.g. control-list(metric=M) where M is matrix in untransformed space
      #    for ADMB models, use MLE covairance with control=list(metric="mle")
      
      # With this adnuts, most important diagnostics are the:
      #   1) ESS (accounts for autocorrelation)-1000 ESS is sufficient for most quantities
      #   2) Potential Scale reduction (R hat)-R hat fails if >1.1
      #   3) No max tree depths exceeded (<12)
      #   4) <1% divergences (ideally 0)
      # reps <- parallel::detectCores()-1 # chains to run in parallel
      reps <- nuts.inputs$reps
      set.seed(63594)
      seeds <- sample(1:1e4, size=reps)
      
      # Jul 2021: Again, nuts.mode=2 used earlier as a check and may not work if 
      # set now
      if(nuts.inputs$nuts.mode==1){
        # shell('PWS_ASA -pinwrite')
        system('./PWS_ASA -pinwrite')
        
        pars <- read.admbFit('PWS_ASA')
        #pars <- read_admb('PWS_ASA')
        inits <- list()
        for(j in 1:reps){
          # This is a check on the CV of pars - basically high uncertainty
          # indicates parameter is fixed or uninformed attempts to avoid
          # boundary issues in further estimation
          if(pars$std[1]/pars$est[1]<100 & !is.infinite(pars$std[1]/pars$est[1])){
            inits[[j]] <- rnorm(1,pars$est[1],sd=pars$std[1])  
          }else{inits[[j]] <- pars$est[1]}
          
          for(k in 2:length(pars$est)){
            if(pars$std[k]/pars$est[k]<100 & !is.infinite(pars$std[k]/pars$est[k])){
              inits[[j]] <- c(inits[[j]],rnorm(1,pars$est[k],sd=pars$std[k]))  
            }else{inits[[j]] <- c(inits[[j]],pars$est[k])}
          }
          
          inits[[j]] <- inits[[j]][-length(inits[[j]])]
        }
        
        start_time <- proc.time()
        nuts.fit <- sample_admb(model='./PWS_ASA',path=targetdir, # change model='PWS_ASA' if using Windows
                                  iter=nuts.inputs$iter,
                                  #iter=1000,
                                  init=inits,
                                  algorithm='NUTS',
                                  seeds=seeds,
                                  chains=reps,
                                  parallel=TRUE,
                                  cores=reps,
                                  #duration=45,
                                  warmup=nuts.inputs$warm_up,
                                  #warmup=100,
                                  mceval=TRUE, 
                                  control=list(adapt_delta=nuts.inputs$adapt_delta))
        end_time <- proc.time()
        
        mon <- monitor(nuts.fit$samples, warmup=nuts.fit$warmup, print=FALSE)
        x <- extract_sampler_params(nuts.fit)
        
        if(sum(x$divergent__)/nrow(x)<=0.01 & max(mon[,'Rhat'])<=1.1){
          nonconvergence <- 0
        }else{nonconvergence <- 1}
        iter.multiplier <- 1 # Not used in this condition - needed here so code doesn't break
      }else if(nuts.inputs$nuts.mode==2){

        shell('PWS_ASA -pinwrite -hbf 1 -nox -iprint 200 -mcmc 15')

        pars <- read.admbFit('PWS_ASA')
        inits <- list()
        for(j in 1:reps){
          # This is a check on the CV of pars - basically high uncertainty indicates parameter is fixed or uninformed
          # attempts to avoid boundary issues in further estimation
          if(pars$std[1]/pars$est[1]<100 | !is.infinite(pars$std[1]/pars$est[1])){
            inits[[j]] <- rnorm(1,pars$est[1],sd=pars$std[1])  
          }else{inits[[j]] <- pars$est[1]}
          
          for(k in 2:length(pars$est)){
            if(pars$std[k]/pars$est[k]<100 | !is.infinite(pars$std[k]/pars$est[k])){
              inits[[j]] <- c(inits[[j]],rnorm(1,pars$est[k],sd=pars$std[k]))  
            }else{inits[[j]] <- c(inits[[j]],pars$est[k])}
          }
          
          inits[[j]] <- inits[[j]][-length(inits[[j]])]
        }
        # inits <- NULL
        start_time <- proc.time()
        nuts.trial <- sample_admb(model='PWS_ASA',path=targetdir,
                             iter=nuts.inputs$iter, 
                             init=inits,
                             algorithm='NUTS',
                             seeds=seeds,
                             chains=reps,
                             parallel=TRUE,
                             cores=reps,
                             #duration=45,
                             warmup=nuts.inputs$warm_up,
                             mceval=TRUE, 
                             control=list(adapt_delta=nuts.inputs$adapt_delta))
        end_time <- proc.time()
        # Launch Shiny App to check diagnostics online
        # launch_shinyadmb(nuts.trial)
        # Do my own monitoring
        mon <- monitor(nuts.trial$samples, warmup=nuts.trial$warmup, print=FALSE)
        x <- extract_sampler_params(nuts.trial)
      
        nonconvergence <- 1
        nuts.fit <- nuts.trial
        # Automatic check for converging chains - if so, use mass matrix from trial to run full chains
        # else, store warm-up results and come back to later
        if(sum(x$divergent__)/nrow(x)<=0.01 & max(mon[,'Rhat'])<=1.1){
          mass <- nuts.trial$covar.est
          inits_2 <- sample_inits(nuts.trial, reps)
          #start_time <- Sys.time()
          if(ceiling(nuts.inputs$iter/min(mon[,'n_eff']))<=1){
            iter.multiplier=2
          }else{
            #500 is the target ESS we want for the slowest mixing parameter
            iter.multiplier<-2.5#ceiling(500/min(mon[,'n_eff']))
          }
          
          nuts.fit <-
            sample_admb(model='PWS_ASA', path=targetdir, 
                        iter=nuts.inputs$iter*iter.multiplier, 
                        init=inits_2, 
                        algorithm='NUTS',
                        seeds=seeds,
                        parallel=TRUE, 
                        chains=reps, 
                        warmup=nuts.inputs$iter*iter.multiplier*0.2, # 20% warmup 
                        cores=reps,
                        mceval=TRUE, 
                        control=list(metric=mass, adapt_delta=nuts.inputs$adapt_delta))
          end_time <- proc.time() 
          #end_time_1 - start_time
          #launch_shinyadmb(nuts.updated)
          mon <- monitor(nuts.fit$samples, warmup=nuts.fit$warmup, print=FALSE)
          x <- extract_sampler_params(nuts.fit)
          nonconvergence <- 0
        }
    }
        
      #launch_shinyadmb(nuts.fit)
        write.csv(mon, file=paste0(folder_1,hypothesis[i],'_table_diagnostics.csv'))
        saveRDS(nuts.fit, file=paste0(folder_2,hypothesis[i],"_NUTS.RDS"))
        locked <- flock::lock(summary_diagnost)
        sum_dia <- data.frame(model=hypothesis[i],
                   divergences.from.extract.function=sum(x$divergent__)/nrow(x),
                   divergences.from.fit.list=sum(nuts.fit$sampler_params[[1]][,5])/length(nuts.fit$sampler_params[[1]][,5]),
                   min.ESS=min(mon[,'n_eff']),
                   which.min.ESS=names(which.min(mon[,'n_eff'])),
                   max.Rhat=max(mon[,'Rhat']),
                   which.max.Rhat=names(which.max(mon[,'Rhat'])),
                   time.elapsed=(end_time-start_time)[3]/60,
                   nonconvergence=nonconvergence,
                   iter.per.chain=ifelse(nonconvergence==1,nuts.inputs$iter,nuts.inputs$iter*iter.multiplier))
        write.table(sum_dia,file=summary_diagnost,append=TRUE,sep=",",col.names=FALSE,row.names=FALSE)
        flock::unlock(locked)
        clean_admb("PWS_ASA",which="sys")
    }
    
  }
  
  print(proc.time() - ptm)
  stopCluster(cl)
}

