function_PSIS_LOO <- function(main_path,mod_name){
  library(loo)
  read.ADMB.files <- function(){
    # Parameters that are fixed and thus included as data
    #Z_3_8 <- 0.25
    egg_add <- 0.4
    
    # Store the file names from which data is available
    filename <- vector(length=2)
    filename[1]="PWS_ASA.dat"
    filename[2]="PWS_ASA(ESS).ctl"
    
    source(file=paste0(here::here("src"),"/data_reader.R"))
    source(file=paste0(here::here("src"),"/data_header_reader.R"))
    
    PWS_ASA.dat <- data_reader(filename=filename[1]) # This is nyr - we want to start at nyr_tobefit
    PWS_ASA.dat <- PWS_ASA.dat[-1]
    nyr <- PWS_ASA.dat[[1]]
    nage <- PWS_ASA.dat[[2]]
    w_a_a <- PWS_ASA.dat[[3]] # From PWS_ASA.dat
    fecun <- PWS_ASA.dat[[4]]
    pc <- PWS_ASA.dat[[5]]
    pk <- PWS_ASA.dat[[6]]
    fbc <- PWS_ASA.dat[[7]]
    gc <- PWS_ASA.dat[[8]]
    sc <- PWS_ASA.dat[[9]]
    f_sp <- PWS_ASA.dat[[10]]
    
    mdm <- PWS_ASA.dat[[11]]
    egg <- PWS_ASA.dat[[12]]
    cv_egg <- PWS_ASA.dat[[13]]
    hydADFG <- PWS_ASA.dat[[15]]
    hydPWSSC <- PWS_ASA.dat[[17]]
    cv_hydPWSSC <- PWS_ASA.dat[[18]]
    seine <- PWS_ASA.dat[[19]]
    spac <- PWS_ASA.dat[[20]]
    
    PWS_ASA_ESS.ctl <- data_reader(filename=filename[2])
    ESS_Se <- PWS_ASA_ESS.ctl[[1]]
    ESS_Sp <- PWS_ASA_ESS.ctl[[2]]
  
    
    seine_indices <- which(rowSums(seine[1:nyr,])>0)
    spawnsurvey_indices <- which(rowSums(spac[1:nyr,])>0)
    egg_indices <- which(egg[1:nyr]>0)
    hydADFG_indices <- which(hydADFG[1:nyr]>0)
    hydPWSSC_indices <- which(hydPWSSC[1:nyr]>0)
    
    model.data <- list(nyr=nyr,
                       nage=nage,
                       w_a_a=w_a_a[1:nyr,],
                       fecun=fecun[1:nyr,],
                       pc=pc[1:nyr,],
                       pk=pk,
                       fbc=fbc[1:nyr,],
                       gc=gc[1:nyr,],
                       sc=sc[1:nyr],
                       f_sp=f_sp[1:nyr],
                       ESS_Se=ESS_Se,
                       ESS_Sp=ESS_Sp,
                       seine=round(seine[1:nyr,],3),
                       spac=round(spac[1:nyr,],3),
                       mdm=mdm[1:nyr],
                       egg=egg[1:nyr],
                       cv_egg=cv_egg[1:nyr],
                       hydADFG=hydADFG[1:nyr],
                       hydPWSSC=hydPWSSC[1:nyr],
                       cv_hydPWSSC=cv_hydPWSSC[1:nyr],
                       seine_indices=seine_indices,
                       spawnsurvey_indices=spawnsurvey_indices,
                       egg_indices=egg_indices,
                       hydADFG_indices=hydADFG_indices,
                       hydPWSSC_indices=hydPWSSC_indices,
                       #Z_3_8=Z_3_8,
                       egg_add=egg_add)
    return(model.data)
  }
  model.data <- read.ADMB.files()
  
  nYr <- model.data$nyr
  ncol <- model.data$nage
  nage <- ncol
  
  # The following reads in the model estimates (draws from MCMC chain of N length) for the different survey indices
  HYD_ADFG <-read.table("HYD_ADFG.csv", header = FALSE, sep = ",", dec=".")
  #HYD_ADFG <-HYD_ADFG[-c(1:nburn),] # Get rid of the burn-in draws
  HYD_PWSSC <-read.table("HYD_PWSSC.csv", header = FALSE, sep = ",", dec=".")
  #HYD_PWSSC <-HYD_PWSSC[-c(1:nburn),]
  MDM <-read.table("MDM.csv", header = FALSE, sep = ",", dec=".")
  #MDM <-MDM[-c(1:nburn),]
  EGG <-read.table("EGG.csv", header = FALSE, sep = ",", dec=".")
  #EGG <-EGG[-c(1:nburn),]
  EGG[EGG==0] <- NA # Replace years with no data (read in as zero to these output files) with NA so that they do not appear on the plot
  VARSReport<-read.table("VarsReport.csv", header = F, sep = ",", dec=".")
  #VARSReport<-VARSReport[-c(1:nburn),]
  
  # Order is m_add, egg_add, hydADFG_add, hydPWSSC_add.
  # m_add:  additional error on the mile-days of milt survey observations (estimated within model)
  # egg_add:  additional error on the egg deposition survey observations (fixed within model)
  # hydADFG_add:  additional error on the ADF&G hydroacoustic survey observations (estimated within model)
  # hydPWSSC_add:  additional error on the PWSSC hydroacoustic survey observations (estimated within model)
  Vars<-data.frame(VARSReport)
  
  # Calculate CV's for each survey type, using both estimated and provided variances (where appropriate)
  egg.se <- model.data$cv_egg
  PWSSC.se <- model.data$cv_hydPWSSC
  
  mediansAdd <- apply(Vars,2, median); #mediansAdd
  m.cv <- mediansAdd[1] #
  egg.cv <- as.vector(sqrt(egg.se^2 + mediansAdd[2]^2))
  ADFG.cv <- mediansAdd[3] # 
  PWSSC_cv <- sqrt((PWSSC.se^2) + (mediansAdd[4]^2))
    
  # Calculated 95% confidence intervals from survey CV (from Buckland 1992 in References)
  calc_buck_C <- function(cv=cv) {
    buck_c <- NULL
    cv <- cv
    buck_c <- exp(1.96*sqrt(log(1+cv*cv)))
    return(buck_c)
  }
  
  MDM.density <- matrix(NA,nrow=NROW(MDM),ncol=NCOL(MDM))
  colnames(MDM.density) <- paste0("MDM_",1:NCOL(MDM)+1979)
  data.indices <- which(model.data$mdm!=-9)
  for(i in 1:length(data.indices)){
    fits <- as.numeric(MDM[,data.indices[i]]) # Take the true value of the survey estimate
    # errors <- as.numeric(Vars[i,1]*fits) # Take the error based on the true value of the survey estimate
    # Need to reparameterize to calculate the log-normal mean and sd for the survey data
    # From: https://msalganik.wordpress.com/2017/01/21/making-sense-of-the-rlnorm-function-in-r/
    # E <- log(fits^2/sqrt(errors^2+fits^2))
    # SD <- sqrt(log(1+(errors^2/fits^2)))
    E <- log(fits)
    SD <- Vars[,1]
    MDM.density[,data.indices[i]] <- dnorm(x=log(model.data$mdm[data.indices[i]]),
                                         mean=E,
                                         sd=SD)
  }
  
  EGG.density <- matrix(NA,nrow=NROW(EGG),ncol=NCOL(EGG))
  colnames(EGG.density) <- paste0("EGG_",1:NCOL(EGG)+1979)
  data.indices <- which(model.data$egg!=-9)
  for(i in 1:length(data.indices)){
    fits <- as.numeric(EGG[,data.indices[i]])
    E <- log(fits)
    SD <- sqrt(egg.se[data.indices[i]]^2+Vars[,2]^2)
    EGG.density[,data.indices[i]] <- dnorm(x=log(model.data$egg[data.indices[i]]),mean=E, sd=SD)
  }
  
  HYD_ADFG.density <- matrix(NA,nrow=NROW(HYD_ADFG),ncol=NCOL(HYD_ADFG))
  colnames(HYD_ADFG.density) <- paste0("HYD_ADFG_",1:NCOL(HYD_ADFG)+1979)
  data.indices <- which(model.data$hydADFG!=-9)
  for(i in 1:length(data.indices)){
    fits <- as.numeric(HYD_ADFG[,data.indices[i]])
    E <- log(fits)
    SD <- Vars[,3]
    HYD_ADFG.density[,data.indices[i]] <- dnorm(x=log(model.data$hydADFG[data.indices[i]]),mean=E,sd=SD)
  }
  
  HYD_PWSSC.density <- matrix(NA,nrow=NROW(HYD_PWSSC),ncol=NCOL(HYD_PWSSC))
  colnames(HYD_PWSSC.density) <- paste0("HYD_PWSSC_",1:NCOL(HYD_PWSSC)+1979)
  data.indices <- which(model.data$hydPWSSC!=-9)
  for(i in 1:length(data.indices)){
    fits <- as.numeric(HYD_PWSSC[,data.indices[i]])
    E <- log(fits)
    SD <- sqrt(PWSSC.se[data.indices[i]]^2+Vars[,4]^2)
    HYD_PWSSC.density[,data.indices[i]] <- dnorm(x=log(model.data$hydPWSSC[data.indices[i]]),mean=E,sd=SD)
  }
  
  SeESS <- model.data$ESS_Se
  SpESS <- model.data$ESS_Sp
  
  SeESS <- SeESS[1:nYr]
  SpESS <- SpESS[1:nYr]
  
  SeESS[SeESS==-9] <- 0
  SpESS[SpESS==-9] <- 0
  
  seData <- model.data$seine
  spData <- model.data$spac 
  seData<-seData*SeESS
  spData<-spData*SpESS
  
  seData[10,]<-rep(0,ncol)
  seData[14:17,]<-rep(0,ncol)
  seData[20:nYr,]<-rep(0,ncol)
  
  DATA.SeAC <- as.vector(t(seData))
  DATA.SpAC <- as.vector(t(spData))
  
  # read in the model estimates for spawners
  ## SpAC is a # of saved draws X 231 matrix (33yrs * 7 age classes)
  SpAC<-read.csv("SpAC.csv", header = FALSE, dec=".") 
  #SpAC<-SpAC[-c(1:nburn),] # Just 'cause it's easier to work with this scale
  ###################
  ## read in the model estimates for seine catch
  SeAC<-read.csv("SeAC.csv", header = FALSE, dec=".") 
  #SeAC<-SeAC[-c(1:nburn),] # Just 'cause it's easier to work with this scale
  
  # Individual densities
  Density.SeAC<-matrix(NA, NROW(SeAC), nYr)
  Density.SpAC<-matrix(NA, NROW(SpAC), nYr)
  
  colnames(Density.SeAC) <- paste0("SeAC_",1:nYr+1979)
  colnames(Density.SpAC) <- paste0("SpAC_",1:nYr+1979)
  
  # sample using the MCMC draws 
  for(j in 1:length(SeESS)){ # Loop through each year
    Density.SeAC[,j] <- apply(SeAC[,(j*ncol-(ncol-1)):(j*ncol)],MARGIN=1,FUN=function(x) dmultinom(x=DATA.SeAC[(j*ncol-(ncol-1)):(j*ncol)],prob=as.numeric(x)))

    if(all(SpAC[, (j*ncol-(ncol-1)) : (j*ncol) ]==0)){
    }else{
      Density.SpAC[,j] <- apply(SpAC[,(j*ncol-(ncol-1)):(j*ncol)],MARGIN=1,FUN=function(x) dmultinom(x=DATA.SpAC[(j*ncol-(ncol-1)):(j*ncol)],prob=as.numeric(x)))
    }
  }
  
  rm(SpAC,SeAC)
  
  #################################### CALCULATE PSIS LOO ####################################
  
  Density.SeAC[Density.SeAC==1] = NA
  Density.SpAC[Density.SpAC==1] = NA
  
  total.density = cbind(MDM.density,EGG.density,HYD_ADFG.density,HYD_PWSSC.density,Density.SeAC,Density.SpAC)
  # rm(MDM.density,EGG.density,HYD_ADFG.density,HYD_PWSSC.density,Density.SeAC,Density.SpAC)
  
  # Remove missing years (columns)
  total.density = total.density[,apply(total.density,2,function(x) all(!is.na(x)))]
  
  #total.density = total.density[c(-2500,-7499),]
  
  fit <- readRDS(paste0(main_path,"NUTS_RDS_files/",mod_name,"_NUTS.RDS"))
  
  # Log-likelihoods
  log.density = log(total.density)
  
  # Compute relative efficiencies
  r_eff = relative_eff(total.density,chain_id = c(rep(1,times=2500),rep(2,times=2500),rep(3,times=2500)))
  x = loo(log.density,r_eff=r_eff)
  #print(x)
  #plot(x)
  
  # These represent influence of observations on model posterior distribution
  # pareto_k_influence_values(x)[pareto_k_influence_values(x)>0.7]
  
  # Observations that are above threshold 
  # colnames(total.density)[pareto_k_ids(x)]
   
  # Output:  elpd_loo, se_elpd_loo, p_loo, se_p_loo, total # parameters, # obs with k>0.7
  # 2nd output:  which obs with k>0.7, k values for those obs
  
  return(list(psis_loo=data.frame(elpd=x$estimates[1,1],se_elpd=x$estimates[1,2],
                                  eff_p=x$estimates[2,1],se_eff_p=x$estimates[2,2],
                                  tot_p=fit$mle$nopar, no_bad_k=sum(pareto_k_influence_values(x)>0.7)),
              k_vals=data.frame(bad_obs=names(pareto_k_influence_values(x)[pareto_k_influence_values(x)>0.7]),
                                bad_k=pareto_k_influence_values(x)[pareto_k_influence_values(x)>0.7])))
}




