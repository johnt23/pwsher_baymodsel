# calc_model_selection.R
# Formerly:  plot_hypothesis_test_results.r
# Created by John Trochta
# Date created:  11/06/2019
# Summary:
# This script calculates several model selection indices and summary outputs of 
# results to include in manuscript
library(ggplot2)
library(dplyr)
library(R2admb)
library(loo)

mod_sel <- function(what.to.calculate,run_ID,basePath,H_set,run_year,run_month){
  source(file=paste0(here::here("src"),"/data_reader.R"))
  source(file=paste0(here::here("src"),"/data_header_reader.R"))
  source(file=paste0(here::here("src"),"/function_effect_posteriors.R"))
  source(file=paste0(here::here("src"),"/function_DIC.R"))
  
  modelPath <-paste0(here::here("results"),"/",run_ID,"/") # paste0(here::here("output"),"/",run_year,"/",run_year,"_",run_month,"/",run_ID,"/")

  resultsPath <- paste0(here::here("results"),"/figs & tables/")
  
  N_burn <- 1 # The number of initial draws to remove from the MCMC chain
  
  
  # MAKE SURE TO CHANGE THESE IF CHANGED IN ADMB FILES
  PWS_ASA.dat <- data_reader(filename=paste0(basePath,"PWS_ASA.dat"))
  nyr <- PWS_ASA.dat[[1]]
  nyr_tobefit <- PWS_ASA.dat[[2]]
  nage <- PWS_ASA.dat[[3]]
  
  covariates <- data_reader(filename=paste0(basePath,"PWS_ASA(covariate).ctl"))
  n_age0_covs <- covariates[[2]]
  n_mor_covs <- covariates[[7]]
  
  variable.names <- c(  'VHSV_age3_4_mort_93',
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
                        'beta_age0',
                        'beta_mortality',
                        paste0('annual_mortdevs',1:nyr_tobefit),
                        'sigma_mortdevs',
                        'beta_age0_offset',
                        'beta_mortality_offset',
                        'sigma_age0covar',
                        'sigma_morcovar')
  
  #Below is a comprehensive list of all covariates/hypotheses I tested
  
  if(H_set=="M"){
    hypothesis_names <- c("I. hoferi pre-2007 (Ages 5+)",
                          "I. hoferi post-2007 (Ages 5+)",
                          "VHSV (Ages 3-4)",
                          "Summer upwelling (Ages 3+)",
                          "Summer PDO (Ages 3+)",
                          "Winter PDO (Ages 3+)",
                          'Summer NPGO (Ages 3+)',
                          'Winter NPGO (Ages 3+)',
                          'Total pink salmon returns (Ages 3+)',
                          'GOA walleye pollock SSB (Ages 3+)',
                          'GOA Pacific cod SB (Ages 3+)',
                          'GOA Arrowtooth flounder female SB (Ages 3+)',
                          'GOA Arrowtooth flounder total biomass (Ages 3+)',
                          "Summer humpback whale estimates (Age 3+)",
                          "Humpback whale counts (Age 3+)",
                          "Disease indices",
                          "Disease",
                          "Null",
                          "Full")
    
    hypothesis <- gsub(pattern=" |-|/",replacement="_",x=hypothesis_names)
    hypothesis <- gsub(pattern="[()']",replacement="",x=hypothesis)
    
    h.block <- c(1:14)
    
    H.within <- which(hypothesis %in% list.files(modelPath))
    h.block[!((h.block) %in% H.within)] <- NA
    hypothesis <- hypothesis[H.within]
    hypothesis_names <- hypothesis_names[H.within]
  }else if(H_set=="R"){
    hypothesis_names <- c("1989 regime shift",
                          "Hatchery-released juvenile pink salmon",
                          "GOA walleye pollock age 1",
                          'Summer NPGO',
                          'Summer PDO',
                          'Age 0 scale growth (lagged 1 year)',
                          'Freshwater discharge',
                          "PWS Avg. Zooplankton density '81-'99",
                          "Full")
    hypothesis <- gsub(pattern=" |-|/",replacement="_",x=hypothesis_names)
    hypothesis <- gsub(pattern="[()']",replacement="",x=hypothesis)
    h.block <- 1:8
    H.within <- 1:8
    
    H.within <- which(hypothesis %in% list.files(modelPath))
    h.block[!((h.block) %in% H.within)] <- NA
    hypothesis <- hypothesis[H.within]
    hypothesis_names <- hypothesis_names[H.within]
  }
  
  
  #####################################################################################
  # COMPARING POSTERIORS ON EFFECTS OF COVARIATES
  # Must loop through as I have to load in .csv files from MCMC results of ADMB model
  # Edit this loop because this is pretty hard coded
  
  if("Estimated Effects" %in% what.to.calculate){
    function_effect_posteriors(hypothesis,
                               modelPath,
                               variable.names,
                               N_burn,
                               resultsPath,
                               run_ID)
  }
  
  #####################################################################################
  # CALCULATE DIC
  # Edit because iterations file cannot be used to directly write PIN file
  if("DIC" %in% what.to.calculate){
    source(here::here("src/function_DIC.R"))
    targetdir <- paste0(resultsPath,"DIC_calculations/")
    if(!dir.exists(targetdir)){
      dir.create(paste0(resultsPath,"DIC_calculations/"))
    }
    
    Model.DIC<- vector(length=length(hypothesis_names))
    for(i in 1:length(hypothesis)){  
      # Change to that working directory
      setwd(paste0(modelPath,hypothesis[i]))
      
      # Copy the files that will remain unchanged to this new directory
      filestocopy <- list.files(paste0(modelPath,hypothesis[i]))
      file.copy(from=filestocopy, to=targetdir)
      Model.DIC[i] <- function_DIC(targetdir,N_burn=1)
      file.remove(list.files())
    }
    setwd(resultsPath)
    unlink(paste0(resultsPath,"DIC_calculations"),recursive = TRUE,force=TRUE)
  }
  
  #####################################################################################
  # CALCULATE WAIC
  if("WAIC" %in% what.to.calculate){
    WAIC <- function (log_lik,N_burn){
      lppd <- sum (log (colMeans(exp(log_lik[-c(1:N_burn),])))) # Checked on 08/28/2018 - looks right with the Hooten and Hobbs equation
      #p_waic_1 <- 2*sum (log(colMeans(exp(log_lik))) - colMeans(log_lik))
      p_waic_2 <- sum (colVars(log_lik[-c(1:N_burn),]))
      waic_2 <- -2*lppd + 2*p_waic_2
      #return (list (waic=waic_2, p_waic=p_waic_2, lppd=lppd, p_waic_1=p_waic_1))
      return (waic_2)
    }
    colVars <- function (a){
      diff <- a - matrix (colMeans(a), nrow(a), ncol(a), byrow=TRUE)
      vars <- colMeans (diff^2)*nrow(a)/(nrow(a)-1)
      return (vars)
    }
    Model.WAIC.1 <- vector(length=length(hypothesis))
    for(i in 1:length(hypothesis)){
      # Change to that working directory
      setwd(paste0(modelPath,hypothesis[i]))
      llk.data <- read.csv("llik_observations.csv",header=FALSE)
      log_lik <- -llk.data[,colSums(llk.data)!=0][,-(1:8)] # Don't want first column because it is priors & penllk and convert from NLL to LL
      Model.WAIC.1[i] <- WAIC(log_lik,N_burn)
    }
  }
  
  #####################################################################################
  # CALCULATE WAIC and POSTERIOR PREDICTIVE LOSS
  if("Predictive Loss" %in% what.to.calculate){
    
    # source(file=paste0(here::here("src"),"/function_posterior_predictive_loss.R"))
    ptm <- proc.time()
    Model.PPL <- matrix(nrow=length(hypothesis),ncol=54)
    Model.WAIC.2 <- vector(length=length(hypothesis))
    
    library(doParallel)
    n_cores <- 9
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    ptm <- proc.time()
    
    post.pred.output <- foreach(i=1:length(hypothesis)) %dopar%{
      source(file=paste0(here::here("src"),"/function_PPL_WAIC.R"))
      #for(i in 1:length(hypothesis)){  
      # Change to that working directory
      setwd(paste0(modelPath,hypothesis[i]))
      nburn=N_burn
      #function_posterior_predictive_loss(nburn=nburn)
      function_PPL_WAIC(nburn=nburn)
    }
    stopCluster(cl)
    proc.time() - ptm
    
    
    for(i in 1:length(hypothesis)){
      if(i==1){
        Model.PPL <- post.pred.output[[i]]$ppl
        Model.fit <- data.frame(model=hypothesis_names[i],post.pred.output[[i]]$model.fits)
        Model.WAIC.2 <- post.pred.output[[i]]$waic
      }else{
        Model.PPL <- rbind(Model.PPL,post.pred.output[[i]]$ppl)
        Model.fit <- rbind(Model.fit,data.frame(model=hypothesis_names[i],post.pred.output[[i]]$model.fits))
        Model.WAIC.2 <- rbind(Model.WAIC.2,post.pred.output[[i]]$waic)
      }
      
    }
    
    final.names <- hypothesis_names
    if(length(final.names)==1){
      Model.PPL = t(as.matrix(Model.PPL))
    }
    model.selection <- data.frame(final.names,
                                  Model.DIC[1:length(final.names)],
                                  delta.Model.DIC=Model.DIC[1:length(final.names)]-min(Model.DIC[1:length(final.names)]),
                                  Model.WAIC.1[1:length(final.names)],
                                  delta.Model.WAIC.1=Model.WAIC.1[1:length(final.names)]-min(Model.WAIC.1[1:length(final.names)]),
                                  Model.WAIC.2[1:length(final.names)],
                                  delta.Model.WAIC.2=Model.WAIC.2[1:length(final.names)]-min(Model.WAIC.2[1:length(final.names)]),
                                  Posterior.predict.loss=Model.PPL[1:length(final.names),1],
                                  delta.posterior.predict.loss=Model.PPL[1:length(final.names),1]-min(Model.PPL[1:length(final.names),1]),
                                  Error.loss.G=Model.PPL[1:length(final.names),2],
                                  Predictive.variance.P=Model.PPL[1:length(final.names),3])
    write.csv(model.selection,file=paste0(resultsPath,run_ID,"_model selection results.csv"))
    
    ppl.selection <- Model.PPL
    colnames(ppl.selection) <- c('Dinf.squ.err','G.squ.err','P.squ.err',
                                 'G.HYD_PWSSC.norm.squ.err','G.HYD_ADFG.norm.squ.err',
                                 'G.EGG.norm.squ.err','G.MDM.norm.squ.err',
                                 'G.SeAC.norm.squ.err','G.SpAC.norm.squ.err',
                                 'P.HYD_PWSSC.norm.squ.err','P.HYD_ADFG.norm.squ.err',
                                 'P.EGG.norm.squ.err','P.MDM.norm.squ.err',
                                 'P.SeAC.norm.squ.err','P.SpAC.norm.squ.err',
                                 'G.HYD_PWSSC.squ.err','G.HYD_ADFG.squ.err',
                                 'G.EGG.squ.err','G.MDM.squ.err','G.SeAC.squ.err',
                                 'G.SpAC.squ.err','P.HYD_PWSSC.squ.err',
                                 'P.HYD_ADFG.squ.err','P.EGG.squ.err',
                                 'P.MDM.squ.err','P.SeAC.squ.err','P.SpAC.squ.err')
    
    # colnames(ppl.selection) <- c('Dinf.squ.err','G.squ.err','P.squ.err',
    #                              'G.HYD_PWSSC.norm.squ.err','G.HYD_ADFG.norm.squ.err',
    #                              'G.EGG.norm.squ.err','G.MDM.norm.squ.err',
    #                              'G.SeAC.norm.squ.err','G.SpAC.norm.squ.err',
    #                              'P.HYD_PWSSC.norm.squ.err','P.HYD_ADFG.norm.squ.err',
    #                              'P.EGG.norm.squ.err','P.MDM.norm.squ.err',
    #                              'P.SeAC.norm.squ.err','P.SpAC.norm.squ.err',
    #                              'G.HYD_PWSSC.squ.err','G.HYD_ADFG.squ.err',
    #                              'G.EGG.squ.err','G.MDM.squ.err','G.SeAC.squ.err',
    #                              'G.SpAC.squ.err','P.HYD_PWSSC.squ.err',
    #                              'P.HYD_ADFG.squ.err','P.EGG.squ.err',
    #                              'P.MDM.squ.err','P.SeAC.squ.err','P.SpAC.squ.err',
    #                              'Dinf.abs.err','G.abs.err','P.abs.err',
    #                              'G.HYD_PWSSC.norm.abs.err','G.HYD_ADFG.norm.abs.err',
    #                              'G.EGG.norm.abs.err','G.MDM.norm.abs.err',
    #                              'G.SeAC.norm.abs.err','G.SpAC.norm.abs.err',
    #                              'P.HYD_PWSSC.norm.abs.err','P.HYD_ADFG.norm.abs.err',
    #                              'P.EGG.norm.abs.err','P.MDM.norm.abs.err','P.SeAC.norm.abs.err',
    #                              'P.SpAC.norm.abs.err','G.HYD_PWSSC.abs.err','G.HYD_ADFG.abs.err',
    #                              'G.EGG.abs.err','G.MDM.abs.err','G.SeAC.abs.err','G.SpAC.abs.err',
    #                              'P.HYD_PWSSC.abs.err','P.HYD_ADFG.abs.err','P.EGG.abs.err',
    #                              'P.MDM.abs.err','P.SeAC.abs.err','P.SpAC.abs.err')
    write.csv(ppl.selection,file=paste0(resultsPath,run_ID,"_posterior predictive loss.csv"))
    
    write.csv(Model.fit,file=paste0(resultsPath,run_ID,"_model fits to data.csv"))
  }
  
  #####################################################################################
  # CALCULATE PSIS LOO CV
  if("PSIS LOO" %in% what.to.calculate){
    ptm <- proc.time()
    library(doParallel)
    n_cores <- 9
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    ptm <- proc.time()
    
    loo.output <- foreach(i=1:length(hypothesis)) %dopar%{
      source(file=paste0(here::here("src"),"/function_PSIS_LOO.R"))
      # Change to that working directory
      setwd(paste0(modelPath,hypothesis[i]))
      function_PSIS_LOO(main_path=modelPath,mod_name=hypothesis[i])
    }
    stopCluster(cl)
    proc.time() - ptm
    
    
    Model.LOO <- matrix(nrow=length(hypothesis),ncol=5)
    Model.LOO.pareto.k <- vector(length=length(hypothesis))
    
    for(i in 1:length(hypothesis)){
      if(i==1){
        Model.LOO <- loo.output[[i]]$psis_loo
        Model.LOO.pareto.k <- data.frame(model=hypothesis_names[i],loo.output[[i]]$k_vals)
      }else{
        Model.LOO <- rbind(Model.LOO,loo.output[[i]]$psis_loo)
        
        if(nrow(loo.output[[i]]$k_vals)==0){
          Model.LOO.pareto.k <- rbind( Model.LOO.pareto.k,
                                       data.frame(model=hypothesis_names[i],bad_obs=NA,bad_k=NA))
        }else{
          Model.LOO.pareto.k <- rbind( Model.LOO.pareto.k,
                                       data.frame(model=hypothesis_names[i],loo.output[[i]]$k_vals))
        }
        
      }
      
    }
    Model.LOO <- data.frame(final.names=hypothesis_names,Model.LOO)
    write.csv(Model.LOO,file=paste0(resultsPath,run_ID,"_PSIS LOO estimates.csv"),row.names = FALSE)
    write.csv(Model.LOO.pareto.k,file=paste0(resultsPath,run_ID,"_PSIS LOO k diagnostics.csv"),row.names = FALSE)
  }
  
  }

OS <- "MAC"
run_year = '2021'
run_month = '07'

run_ID <- "2021_07_all fixed mortality effects"
what.to.calculate <- c("Estimated Effects","DIC","WAIC","Predictive Loss","PSIS LOO")
basePath <- paste0(here::here("bin"),"/PWS_ASA_all fixed effects/")
H_set <- "M"
ptm = proc.time()
mod_sel(what.to.calculate,run_ID,basePath,H_set,run_year,run_month)
proc.time() - ptm


run_ID <- "2021_07_all fixed mortality effects _ start from 1994"
basePath <- paste0(here::here("bin"),"/PWS_ASA_all fixed effects/")
H_set <- "M"
ptm = proc.time()
mod_sel(what.to.calculate,run_ID,basePath,H_set,run_year,run_month)
proc.time() - ptm


run_ID <- "2021_07_all fixed age0 effects"
#what.to.calculate <- c("Estimated Effects","DIC","WAIC","Predictive Loss")
basePath <- paste0(here::here("bin"),"/PWS_ASA_all fixed effects _ pen rec/")
H_set <- "R"
ptm = proc.time()
mod_sel(what.to.calculate,run_ID,basePath,H_set,run_year,run_month)
proc.time() - ptm


run_ID <- "2021_07_all fixed age0 effects _ start from 1994"
#what.to.calculate <- c("Estimated Effects","DIC","WAIC","Predictive Loss")
basePath <- paste0(here::here("bin"),"/PWS_ASA_all fixed effects _ pen rec/")
H_set <- "R"
ptm = proc.time()
mod_sel(what.to.calculate,run_ID,basePath,H_set,run_year,run_month)
proc.time() - ptm


run_ID <- "2021_07_latent mortality effects _ year specific SE"
what.to.calculate <- c("DIC","WAIC","Predictive Loss","PSIS LOO")
basePath <- paste0(here::here("bin"),"/PWS_ASA_latent vars for mortality _ year specific SE/")
H_set <- "M"
ptm = proc.time()
mod_sel(what.to.calculate,run_ID,basePath,H_set,run_year,run_month)
proc.time() - ptm


run_ID <- "2021_07_latent age0 effects _ year specific SE"
what.to.calculate <- c("DIC","WAIC","Predictive Loss","PSIS LOO")
basePath <- paste0(here::here("bin"),"/PWS_ASA_latent vars for age0 _ year specific SE/")
H_set <- "R"
ptm = proc.time()
mod_sel(what.to.calculate,run_ID,basePath,H_set,run_year,run_month)
proc.time() - ptm


run_ID <- "2021_07_all fixed age0 effects_low sigmaR"
run_year = '2021'
run_month = '07'
what.to.calculate <- c("Estimated Effects","DIC","WAIC","Predictive Loss","PSIS LOO")
basePath <- paste0(here::here("bin"),"/PWS_ASA_all fixed effects _ pen rec_low sigmaR/")
H_set <- "R"
ptm = proc.time()
mod_sel(what.to.calculate,run_ID,basePath,H_set,run_year,run_month)
proc.time() - ptm
