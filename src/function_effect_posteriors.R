# function_effect_posteriors.R
# Created by: John Trochta
# Date:01/14/2020
# Summary: Takes BASA runs, extracts estimated covariates, 
#          and plots/stores covariate posteriors (median & 95% credibility intervals)

function_effect_posteriors<- function(hypothesis,
         modelPath,
         variable.names,
         N_burn,
         resultsPath,
         run_ID){
  
  # Loop through hypothesis name
  
  start.index <- 1
  if("Null_model" %in% hypothesis){
    hypothesis <- hypothesis[1:(length(hypothesis)-1)]
    #hypothesis_names <- hypothesis_names[1:(length(hypothesis_names)-1)]
  }
  
  for(i in start.index:(length(hypothesis))){
    # Change to that working directory
    targetdir <- paste0(modelPath,hypothesis[i],"/")
    setwd(targetdir)
    
    covariates <- data_reader(filename=paste0(targetdir,"PWS_ASA(covariate).ctl"))
    age0_cov_index <- covariates[[4]]
    age0_offset_index <- covariates[[6]]
    mort_cov_index <- covariates[[10]]
    mort_offset_index <- covariates[[14]]
    
    phases <- data_reader(filename=paste0(targetdir,"PWS_ASA(phases).ctl"))
    set_beta <- c(phases[[26]],phases[[27]],phases[[29]],phases[[30]])
    set_beta <- which(set_beta>0)
    
    PWS_ASA.dat <- data_reader(filename=paste0(basePath,"PWS_ASA.dat"))
    nyr <- PWS_ASA.dat[[2]]
    
    # This next part modifies the variable_names vectors for the MCMC draws matrix for parameters
    # (ADMB code adaptively creates parameter vectors associated with covariates to dimension
    # that matches number of covariates included)
    beta = c("beta_age0","beta_mortality",
             "beta_age0_offset","beta_mortality_offset")
    beta_1 = beta[set_beta]
    beta_names = beta[set_beta]
    
    if(sum(age0_cov_index)>1){
      if(beta[1]%in%beta_1){
        beta_names = c(beta_names,paste0(beta[1],"_",2:sum(age0_cov_index)))
      }
      if(beta[3]%in%beta_1){
        beta_names = c(beta_names,paste0(beta[3],"_",2:sum(age0_cov_index)))
      }
      age0_pars = c("beta_age0","beta_age0_offset")
      for(j in 1:length(age0_pars)){
        variable.names = append(variable.names,
                                paste0(age0_pars[j],"_",2:sum(age0_cov_index)),
                                after = which(variable.names==age0_pars[j]))
      }
    }
    
    if(sum(mort_cov_index)>1){
      if(beta[2]%in%beta_1){
        beta_names = c(beta_names,paste0(beta[2],"_",2:sum(mort_cov_index)))
      }
      if(beta[4]%in%beta_1){
        beta_names = c(beta_names,paste0(beta[4],"_",2:sum(mort_cov_index)))
      }
      # If multiple covariates are included, need to also expand annual_mortdevs, beta_mort
      mort_pars = c("beta_mortality","beta_mortality_offset","sigma_morcovar")
      for(j in 1:length(mort_pars)){
        variable.names = append(variable.names,
                                paste0(mort_pars[j],"_",2:sum(mort_cov_index)),
                                after = which(variable.names==mort_pars[j]))
      }
      # START HERE - NEED TO INDEX UP FOR EACH SET OF MORTDEVS
      mortdev_names = paste0(rep(paste0('annual_mortdevs',1:nyr),sum(mort_cov_index)-1),"_",
                             rep(2:sum(mort_cov_index),each=nyr))
      variable.names = append(variable.names,
                              mortdev_names,
                              after = which(variable.names==paste0("annual_mortdevs",nyr)))
    }
  
    Outs<-read.table(paste0(targetdir,"iterations.csv"), header = FALSE, sep = ",", dec=".")
    names(Outs) <- c(variable.names,"objective_fn_val")
    if(NROW(Outs)!=0 ){
      Outs<-Outs[-seq(1,N_burn,by=1),] # Removing fixed parameters
      if(i==start.index){
        Outs <- data.frame(Outs[,names(Outs)%in%beta_names])
        names(Outs)[1] <- hypothesis[i]
        if(ncol(Outs)>1){
          names(Outs)[2:ncol(Outs)] <- paste0(hypothesis[i],'_',2:ncol(Outs))
        }
        hypothesis_posteriors <- data.frame(Outs)
      }else{
        Outs <- data.frame(Outs[,names(Outs)%in%beta_names])
        names(Outs)[1] <- hypothesis[i]
        if(ncol(Outs)>1){
          names(Outs)[2:ncol(Outs)] <- paste0(hypothesis[i],'_',2:ncol(Outs))
        }
        hypothesis_posteriors <- data.frame(hypothesis_posteriors,Outs)
      }
    }
  }
  
  library(reshape)
  hyp_post <- melt(hypothesis_posteriors)
  
  # Save all the posterior draws from each model
  write.csv(hyp_post,file=paste0(resultsPath,run_ID,"_effect posteriors.csv"),row.names=FALSE)
  
  hyp_post <- hyp_post %>% group_by(variable) %>% summarise(
    minim=min(value),
    lower.99.9th=quantile(value,probs=0.0005),
    lower.99th=quantile(value,probs=0.005),
    lower.95th=quantile(value,probs=0.025),
    median.50th=quantile(value,probs=0.5),
    upper.95th=quantile(value,probs=0.975),
    upper.99th=quantile(value,probs=0.995),
    upper.99.9th=quantile(value,probs=0.9995),
    maxim=max(value))

    hypothesis_names_2 <- hypothesis
  
  hyp_post <- data.frame(hyp_post,id=1:NROW(hyp_post),covariates=hypothesis_names_2[start.index:length(hypothesis_names_2)])
  
  covariate.effect.posteriors <- data.frame(hyp_post$covariates,
                                            min.posterior.effect=c(hyp_post$minim),
                                            lower.99.9th.credibility.effect=c(hyp_post$lower.99.9th),
                                            lower.99th.credibility.effect=c(hyp_post$lower.99th),
                                            lower.95th.credibility.effect=c(hyp_post$lower.95th),
                                            median.posterior.effect=c(hyp_post$median.50th),
                                            upper.95th.credibility.effect=c(hyp_post$upper.95th),
                                            upper.99th.credibility.effect=c(hyp_post$upper.99th),
                                            upper.99.9th.credibility.effect=c(hyp_post$upper.99.9th),
                                            max.posterior.effect=c(hyp_post$maxim))
  write.csv(covariate.effect.posteriors,file=paste0(resultsPath,run_ID,"_effect credibility intervals.csv"))
}
