function_WAIC <- function(nburn=500){
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
  
  # The following reads in the model estimates (draws from MCMC chain of N length) for the different survey indices
  HYD_ADFG <-read.table("HYD_ADFG.csv", header = FALSE, sep = ",", dec=".")
  HYD_ADFG <-HYD_ADFG[-c(1:nburn),] # Get rid of the burn-in draws
  HYD_PWSSC <-read.table("HYD_PWSSC.csv", header = FALSE, sep = ",", dec=".")
  HYD_PWSSC <-HYD_PWSSC[-c(1:nburn),]
  MDM <-read.table("MDM.csv", header = FALSE, sep = ",", dec=".")
  MDM <-MDM[-c(1:nburn),]
  EGG <-read.table("EGG.csv", header = FALSE, sep = ",", dec=".")
  EGG <-EGG[-c(1:nburn),]
  EGG[EGG==0] <- NA # Replace years with no data (read in as zero to these output files) with NA so that they do not appear on the plot
  VARSReport<-read.table("VarsReport.csv", header = F, sep = ",", dec=".")
  VARSReport<-VARSReport[-c(1:nburn),]
  
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
  data.indices <- model.data$mdm!=-9
  for(i in 1:NROW(MDM)){
    fits <- as.numeric(MDM[i,]) # Take the true value of the survey estimate
    # errors <- as.numeric(Vars[i,1]*fits) # Take the error based on the true value of the survey estimate
    # Need to reparameterize to calculate the log-normal mean and sd for the survey data
    # From: https://msalganik.wordpress.com/2017/01/21/making-sense-of-the-rlnorm-function-in-r/
    # E <- log(fits^2/sqrt(errors^2+fits^2))
    # SD <- sqrt(log(1+(errors^2/fits^2)))
    E <- log(fits)
    SD <- Vars[i,1]
    MDM.density[i,data.indices] <- dnorm(x=log(model.data$mdm[data.indices]),
                                         mean=E[data.indices],
                                         sd=SD)# Generate a random variable based on the expected value (predictive) using the assumed 
  }
  
  EGG.density <- matrix(NA,nrow=NROW(EGG),ncol=NCOL(EGG))
  data.indices <- model.data$egg!=-9
  for(i in 1:NROW(EGG)){
    fits <- as.numeric(EGG[i,])
    E <- log(fits)
    SD <- sqrt(egg.se[data.indices]^2+Vars[i,2]^2)
    EGG.density[i,data.indices] <- dnorm(x=log(model.data$egg[data.indices]),
                                         mean=E[data.indices],
                                         sd=SD)
  }
  
  HYD_ADFG.density <- matrix(NA,nrow=NROW(HYD_ADFG),ncol=NCOL(HYD_ADFG))
  data.indices <- model.data$hydADFG!=-9
  for(i in 1:NROW(HYD_ADFG)){
    fits <- as.numeric(HYD_ADFG[i,])
    E <- log(fits)
    SD <- Vars[i,3]
    HYD_ADFG.density[i,data.indices] <- dnorm(x=log(model.data$hydADFG[data.indices]),
                                  mean=E[data.indices],
                                  sd=SD)
  }
  
  HYD_PWSSC.density <- matrix(NA,nrow=NROW(HYD_PWSSC),ncol=NCOL(HYD_PWSSC))
  data.indices <- model.data$hydPWSSC!=-9
  for(i in 1:NROW(HYD_PWSSC)){
    fits <- as.numeric(HYD_PWSSC[i,])
    E <- log(fits)
    SD <- sqrt(PWSSC.se[data.indices]^2+Vars[i,4]^2)
    HYD_PWSSC.density[i,data.indices] <- dnorm(x=log(model.data$hydPWSSC[data.indices]),
                                    mean=E[data.indices],
                                    sd=SD)
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
  SpAC<-SpAC[-c(1:nburn),]*100 # Just 'cause it's easier to work with this scale
  ###################
  ## read in the model estimates for seine catch
  SeAC<-read.csv("SeAC.csv", header = FALSE, dec=".") 
  SeAC<-SeAC[-c(1:nburn),]*100 # Just 'cause it's easier to work with this scale
  
  ## SEINE
  Density.SeAC<-matrix(NA, NROW(SeAC), nYr) # CHANGE
  Density.SpAC<-matrix(NA, NROW(SpAC), nYr)
  
  # Calculate densities of each MCMC draw 
  for(i in 1:length(SeAC[,1])){ # Loop through the MCMC draws
    for(j in 1:length(SeESS)){ # Loop through each year
      Density.SeAC[i,j] <- dmultinom(x=DATA.SeAC[(j*ncol-(ncol-1)):(j*ncol)],
                                      prob=as.numeric(SeAC[i,(j*ncol-(ncol-1)):(j*ncol)]))
      if(all(SpAC[i,(j*ncol-(ncol-1)):(j*ncol)]==0)){}else{
        Density.SpAC[i,j] <- dmultinom(x=DATA.SpAC[(j*ncol-(ncol-1)):(j*ncol)],
                                       prob=as.numeric(SpAC[i,(j*ncol-(ncol-1)):(j*ncol)]))
      }
    }
  }
  rm(SpAC,SeAC)
  Density.SeAC[Density.SeAC==1] = NA
  Density.SpAC[Density.SpAC==1] = NA
  
  total.density = cbind(MDM.density,EGG.density,HYD_ADFG.density,HYD_PWSSC.density,Density.SeAC,Density.SpAC)
  # Remove missing years (columns)
  total.density = total.density[,apply(total.density,2,function(x) all(!is.na(x)))]
  
  lppd <- sum (log (colMeans(total.density))) # Checked on 08/28/2018 - looks right with the Hooten and Hobbs equation
  
  colVars <- function (a){
    diff <- a - matrix (colMeans(a), nrow(a), ncol(a), byrow=TRUE)
    vars <- colMeans (diff^2)*nrow(a)/(nrow(a)-1)
    return (vars)
  }
  p_waic_2 <- sum (colVars(log(total.density)))
  
  waic_2 <- -2*lppd + 2*p_waic_2
  
  return(waic_2)
}




