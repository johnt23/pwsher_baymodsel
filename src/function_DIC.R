# function_DIC.R
# Created by: John Trochta
# Date:01/14/2020
# Summary: Takes BASA runs, extracts estimated covariates, 
#          and plots/stores covariate posteriors (median & 95% credibility intervals)

function_DIC <- function(targetdir,N_burn=1,if.null.model=FALSE,OS="MAC") {
  setwd(targetdir)
  # Model deviance
  llk <- read.table(file = paste0(targetdir,"/llikcomponents.csv"), header = FALSE, sep = ",", dec=".")  
  # Extract specific objective function components, namely the likelihoods & priors
  # modal <- rowSums(llk[,c(1:6,11:(NCOL(llk)-1))])
  # modal <- -llk[,1] -llk[,2] -rowSums(llk[,c(3:6,11:(NCOL(llk)-1))])
  modal <- llk[,1] +llk[,2] +rowSums(llk[,c(3:6)]) # Only take conditional likelihood of data
  # modal <- llk[,NCOL(llk)]
  D.bar <- mean(2*modal[-(1:N_burn)]) # llk is the NLL
  rm(llk)
  
  # Model complexity
  # params <- read.table(file = paste0(targetdir,"/iterations.csv"), header = TRUE, sep = ",", dec=".")
  # params <- read_psv('PWS_ASA') # Reads in MCMC draws for only active parameters
  params<-read.table(paste0(targetdir,"iterations.csv"), header = FALSE, sep = ",", dec=".")
  names(params) <- c(variable.names,"objective_fn_val")
  params <- params[-c(1:N_burn),]
  
  params <- colMeans(params,na.rm=TRUE)
  
  #new.pin <- append(params,rep(0,times=nyr-nyr_tobefit),after=which(names(params)==paste0('annual_age0devs',nyr_tobefit)))
  #new.pin <- append(new.pin,rep(0,times=nyr-nyr_tobefit),after=which(names(new.pin)==paste0('annual_mortdevs',nyr_tobefit)))
  new.pin <- params
  write("# Mean of par posteriors",file="NEW.PIN",append=FALSE)
  for(j in 1:(length(new.pin))){
    write(new.pin[j],file="NEW.PIN",append=TRUE)
  }
  
  # Create my new PIN file
  # file.remove("pws_asa.PIN")
  # file.rename(from = "NEW.PIN", to = "PWS_ASA.PIN") 
  if(OS=="MAC"){
    #system('admb -s PWS_ASA', ignore.stdout=TRUE)
    system('./PWS_ASA -ainp NEW.PIN -noest', ignore.stdout=TRUE)
  }else if(OS=="PC"){
    #shell('admb PWS_ASA')
    shell('PWS_ASA -ainp NEW.PIN -noest')
  }
  file.remove("NEW.PIN")
  
  # START HERE!!
  # f <- as.numeric(readLines("deterministic_run.rep")[6])
  f <- scan("deterministic_run.rep",skip=5,nlines=1)
  
  # We only want the data likelihood components, not the priors OR penalized likelihoods
  f <- sum(f[1:6])
  # f <- read_pars("PWS_ASA",warn_nonstd_rep=FALSE)$loglik
  # f <- read.csv("ParReport.csv")$f_llk
  D.theta.bar <- 2*f
  
  DIC <- D.bar-D.theta.bar+D.bar
  return(DIC)
}