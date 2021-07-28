function_PPL_WAIC <- function(nburn=500){
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
  
  MDM.PP <- matrix(NA,nrow=NROW(MDM),ncol=NCOL(MDM))
  MDM.density <- matrix(NA,nrow=NROW(MDM),ncol=NCOL(MDM))
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
    # Generate a random variable based on the expected value (predictive) using the assumed 
    MDM.PP[,data.indices[i]] <- exp(rnorm(n=length(fits),mean=E,sd=SD))
  }
  
  EGG.PP <- matrix(NA,nrow=NROW(EGG),ncol=NCOL(EGG))
  EGG.density <- matrix(NA,nrow=NROW(EGG),ncol=NCOL(EGG))
  data.indices <- which(model.data$egg!=-9)
  for(i in 1:length(data.indices)){
    fits <- as.numeric(EGG[,data.indices[i]])
    E <- log(fits)
    SD <- sqrt(egg.se[data.indices[i]]^2+Vars[,2]^2)
    EGG.density[,data.indices[i]] <- dnorm(x=log(model.data$egg[data.indices[i]]),mean=E, sd=SD)
    # Generate a random variable based on the expected value (predictive) using the assumed 
    EGG.PP[,data.indices[i]] <- exp(rnorm(n=length(fits),mean=E,sd=SD))
  }
  
  HYD_ADFG.PP <- matrix(NA,nrow=NROW(HYD_ADFG),ncol=NCOL(HYD_ADFG))
  HYD_ADFG.density <- matrix(NA,nrow=NROW(HYD_ADFG),ncol=NCOL(HYD_ADFG))
  data.indices <- which(model.data$hydADFG!=-9)
  for(i in 1:length(data.indices)){
    fits <- as.numeric(HYD_ADFG[,data.indices[i]])
    E <- log(fits)
    SD <- Vars[,3]
    HYD_ADFG.density[,data.indices[i]] <- dnorm(x=log(model.data$hydADFG[data.indices[i]]),mean=E,sd=SD)
    # Generate a random variable based on the expected value (predictive) using the assumed 
    HYD_ADFG.PP[,data.indices[i]] <- exp(rnorm(n=length(fits),mean=E,sd=SD))
  }
  
  HYD_PWSSC.PP <- matrix(NA,nrow=NROW(HYD_PWSSC),ncol=NCOL(HYD_PWSSC))
  HYD_PWSSC.density <- matrix(NA,nrow=NROW(HYD_PWSSC),ncol=NCOL(HYD_PWSSC))
  data.indices <- which(model.data$hydPWSSC!=-9)
  for(i in 1:length(data.indices)){
    fits <- as.numeric(HYD_PWSSC[,data.indices[i]])
    E <- log(fits)
    SD <- sqrt(PWSSC.se[data.indices[i]]^2+Vars[,4]^2)
    HYD_PWSSC.density[,data.indices[i]] <- dnorm(x=log(model.data$hydPWSSC[data.indices[i]]),mean=E,sd=SD)
    # Generate a random variable based on the expected value (predictive) using the assumed 
    HYD_PWSSC.PP[,data.indices[i]] <- exp(rnorm(n=length(fits),mean=E,sd=SD))
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
  SpAC<-SpAC[-c(1:nburn),] # Just 'cause it's easier to work with this scale
  ###################
  ## read in the model estimates for seine catch
  SeAC<-read.csv("SeAC.csv", header = FALSE, dec=".") 
  SeAC<-SeAC[-c(1:nburn),] # Just 'cause it's easier to work with this scale
  
  # Posterior predictions
  PPSeAC<-matrix(0, NROW(SeAC), NCOL(SeAC)) 
  PPSpAC<-matrix(0, NROW(SpAC), NCOL(SpAC))
  
  # Individual densities
  Density.SeAC<-matrix(NA, NROW(SeAC), nYr)
  Density.SpAC<-matrix(NA, NROW(SpAC), nYr)
  
  # sample using the MCMC draws 
  for(j in 1:length(SeESS)){ # Loop through each year
    PPSeAC[,(j*ncol-(ncol-1)):(j*ncol)] <- t(apply(SeAC[,(j*ncol-(ncol-1)):(j*ncol)],MARGIN=1,FUN=function(x) rmultinom(1,size=as.integer(SeESS[j]),x)))
    Density.SeAC[,j] <- apply(SeAC[,(j*ncol-(ncol-1)):(j*ncol)],MARGIN=1,FUN=function(x) dmultinom(x=DATA.SeAC[(j*ncol-(ncol-1)):(j*ncol)],prob=as.numeric(x)))

    if(all(SpAC[, (j*ncol-(ncol-1)) : (j*ncol) ]==0)){
      PPSpAC[, (j*ncol-(ncol-1)) : (j*ncol) ] <- 0
    }else{
      PPSpAC[,(j*ncol-(ncol-1)):(j*ncol)] <- t(apply(SpAC[,(j*ncol-(ncol-1)):(j*ncol)],MARGIN=1,FUN=function(x) rmultinom(1,size=as.integer(SpESS[j]),x)))
      Density.SpAC[,j] <- apply(SpAC[,(j*ncol-(ncol-1)):(j*ncol)],MARGIN=1,FUN=function(x) dmultinom(x=DATA.SpAC[(j*ncol-(ncol-1)):(j*ncol)],prob=as.numeric(x)))
    }
  }
  
  # sample using the MCMC draws 
  # for(i in 1:length(SeAC[,1])){ # Loop through the MCMC draws
  #   for(j in 1:length(SeESS)){ # Loop through each year
  #     PPSeAC[i, (j*ncol-(ncol-1)):(j*ncol) ] <- t(rmultinom(1,size=as.integer(SeESS[j]),SeAC[i,(j*ncol-(ncol-1)):(j*ncol)]))
  #     
  #     Density.SeAC[i,j] <- dmultinom(x=DATA.SeAC[(j*ncol-(ncol-1)):(j*ncol)],
  #                                    prob=as.numeric(SeAC[i,(j*ncol-(ncol-1)):(j*ncol)]))
  #     if(all(SpAC[ i , (j*ncol-(ncol-1)) : (j*ncol) ]==0)){
  #       PPSpAC[i, (j*ncol-(ncol-1)) : (j*ncol) ] <- rep(0,times=ncol)
  #     }else{
  #       PPSpAC[i, (j*ncol-(ncol-1)) : (j*ncol) ] <- t(rmultinom(1,SpESS[j], SpAC[ i , (j*ncol-(ncol-1)) : (j*ncol) ] ) ) 
  #     
  #       Density.SpAC[i,j] <- dmultinom(x=DATA.SpAC[(j*ncol-(ncol-1)):(j*ncol)],
  #                                      prob=as.numeric(SpAC[i,(j*ncol-(ncol-1)):(j*ncol)]))
  #     }
  #   }
  # }
  rm(SpAC,SeAC)
  
  #################################### CALCULATE WAIC ####################################
  
  Density.SeAC[Density.SeAC==1] = NA
  Density.SpAC[Density.SpAC==1] = NA
  
  total.density = cbind(MDM.density,EGG.density,HYD_ADFG.density,HYD_PWSSC.density,Density.SeAC,Density.SpAC)
  rm(MDM.density,EGG.density,HYD_ADFG.density,HYD_PWSSC.density,Density.SeAC,Density.SpAC)
  
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
  
  rm(total.density)
   
  #################################### CALCULATE G ####################################
  PPSeAC <- signif(PPSeAC,3)
  PPSpAC <- signif(PPSpAC,3)
  
  # For normalization of age comp, can divide each age-specific sample size by 
  # mean of that age across all years such as below
  # seine_age_means = rep(colMeans(matrix(DATA.SeAC,ncol=10,byrow=TRUE)),times=nYr)
  # spawn_age_means = rep(colMeans(matrix(DATA.SpAC,ncol=10,byrow=TRUE)),times=nYr)
  
  # or use the mean of the sample size from each year AND each age
  # seine_age_means = mean(DATA.SeAC[DATA.SeAC>0])
  # spawn_age_means = mean(DATA.SpAC[DATA.SpAC>0])
  
  # OR instead I use the mean sample size across all YEARS (the sum of ages in each year)
  seine_age_means = mean(SeESS[SeESS>0])
  spawn_age_means = mean(SeESS[SpESS>0])
  
  #   HYD_PWSSC.G <- abs(apply(HYD_PWSSC.PP[,model.data$hydPWSSC!=-9],2,median)-model.data$hydPWSSC[model.data$hydPWSSC!=-9])
  #   HYD_ADFG.G <- abs(apply(HYD_ADFG.PP[,model.data$hydADFG!=-9],2,median)-model.data$hydADFG[model.data$hydADFG!=-9])
  #   EGG.G <- abs(apply(EGG.PP[,model.data$egg!=-9],2,median)-model.data$egg[model.data$egg!=-9])
  #   MDM.G <- abs(apply(MDM.PP[,model.data$mdm!=-9],2,median)-model.data$mdm[model.data$mdm!=-9])
  #   SeAC.G <- abs(apply(PPSeAC[,DATA.SeAC>0],2,median)-DATA.SeAC[DATA.SeAC>0])
  #   SpAC.G <- abs(apply(PPSpAC[,DATA.SpAC>0],2,median)-DATA.SpAC[DATA.SpAC>0])
  #   
  #   HYD_PWSSC.G.1 <- HYD_PWSSC.G/IQR(model.data$hydPWSSC[model.data$hydPWSSC!=-9])
  #   HYD_ADFG.G.1 <- HYD_ADFG.G/IQR(model.data$hydADFG[model.data$hydADFG!=-9])
  #   EGG.G.1 <- EGG.G/IQR(model.data$egg[model.data$egg!=-9])
  #   MDM.G.1 <- MDM.G/IQR(model.data$mdm[model.data$mdm!=-9])
  #   SeAC.G.1 <- SeAC.G/IQR(DATA.SeAC[DATA.SeAC>0])
  #   SpAC.G.1 <- SpAC.G/IQR(DATA.SpAC[DATA.SpAC>0])
  # 
  # G.data.AbsErr <- c(sum(HYD_PWSSC.G),sum(HYD_ADFG.G),sum(EGG.G),sum(MDM.G),sum(SeAC.G),sum(SpAC.G))
  # G.data.NormAbsErr <- c(sum(HYD_PWSSC.G.1),sum(HYD_ADFG.G.1),sum(EGG.G.1),sum(MDM.G.1),sum(SeAC.G.1),sum(SpAC.G.1))
  # G.AbsErr <- sum(HYD_PWSSC.G.1) + sum(HYD_ADFG.G.1) + sum(EGG.G.1) + sum(MDM.G.1) + sum(SeAC.G.1) + sum(SpAC.G.1) 
  # 
    HYD_PWSSC.G <- (apply(HYD_PWSSC.PP[,model.data$hydPWSSC!=-9],2,mean)-model.data$hydPWSSC[model.data$hydPWSSC!=-9])^2
    HYD_ADFG.G <- (apply(HYD_ADFG.PP[,model.data$hydADFG!=-9],2,mean)-model.data$hydADFG[model.data$hydADFG!=-9])^2
    EGG.G <- (apply(EGG.PP[,model.data$egg!=-9],2,mean)-model.data$egg[model.data$egg!=-9])^2
    MDM.G <- (apply(MDM.PP[,model.data$mdm!=-9],2,mean)-model.data$mdm[model.data$mdm!=-9])^2
    SeAC.G <- (apply(PPSeAC[,DATA.SeAC>0],2,mean)-DATA.SeAC[DATA.SeAC>0])^2
    SpAC.G <- (apply(PPSpAC[,DATA.SpAC>0],2,mean)-DATA.SpAC[DATA.SpAC>0])^2 
    
    HYD_PWSSC.G.1 <- HYD_PWSSC.G/mean(model.data$hydPWSSC[model.data$hydPWSSC!=-9])^2
    HYD_ADFG.G.1 <- HYD_ADFG.G/mean(model.data$hydADFG[model.data$hydADFG!=-9])^2
    EGG.G.1 <- EGG.G/mean(model.data$egg[model.data$egg!=-9])^2
    MDM.G.1 <- MDM.G/mean(model.data$mdm[model.data$mdm!=-9])^2
    SeAC.G.1 <- SeAC.G/seine_age_means^2
    SpAC.G.1 <- SpAC.G/spawn_age_means^2
    
  G.data.SquErr <- c(sum(HYD_PWSSC.G),sum(HYD_ADFG.G),sum(EGG.G),sum(MDM.G),sum(SeAC.G),sum(SpAC.G))
  G.data.NormSquErr <- c(sum(HYD_PWSSC.G.1),sum(HYD_ADFG.G.1),sum(EGG.G.1),sum(MDM.G.1),sum(SeAC.G.1),sum(SpAC.G.1))
  G.SquErr <- sum(HYD_PWSSC.G.1) + sum(HYD_ADFG.G.1) + sum(EGG.G.1) + sum(MDM.G.1) + sum(SeAC.G.1) + sum(SpAC.G.1) 
  
  #################################### CALCULATE P ####################################

    # HYD_PWSSC.P <- apply(HYD_PWSSC.PP[,model.data$hydPWSSC!=-9],2,FUN=function(X) median(abs(X-median(X))))
    # HYD_ADFG.P  <- apply(HYD_ADFG.PP[,model.data$hydADFG!=-9],2,FUN=function(X) median(abs(X-median(X))))
    # EGG.P       <- apply(EGG.PP[,model.data$egg!=-9],2,FUN=function(X) median(abs(X-median(X))))
    # MDM.P       <- apply(MDM.PP[,model.data$mdm!=-9],2,FUN=function(X) median(abs(X-median(X))))
    # SeAC.P      <- apply(PPSeAC[,DATA.SeAC>0],2,FUN=function(X) median(abs(X-median(X))))
    # SpAC.P      <- apply(PPSpAC[,DATA.SpAC>0],2,FUN=function(X) median(abs(X-median(X))))
    # 
    # HYD_PWSSC.P.1 <- HYD_PWSSC.P/IQR(model.data$hydPWSSC[model.data$hydPWSSC!=-9])
    # HYD_ADFG.P.1  <- HYD_ADFG.P/IQR(model.data$hydADFG[model.data$hydADFG!=-9])
    # EGG.P.1       <- EGG.P/IQR(model.data$egg[model.data$egg!=-9])
    # MDM.P.1       <- MDM.P/IQR(model.data$mdm[model.data$mdm!=-9])
    # SeAC.P.1      <- SeAC.P/IQR(DATA.SeAC[DATA.SeAC>0])
    # SpAC.P.1      <- SpAC.P/IQR(DATA.SpAC[DATA.SpAC>0])
  #}
    # P.data.AbsErr <- c(sum(HYD_PWSSC.P),sum(HYD_ADFG.P),sum(EGG.P),sum(MDM.P),sum(SeAC.P),sum(SpAC.P))
    # P.data.NormAbsErr <- c(sum(HYD_PWSSC.P.1),sum(HYD_ADFG.P.1),sum(EGG.P.1),sum(MDM.P.1),sum(SeAC.P.1),sum(SpAC.P.1))
    # P.AbsErr <- sum(HYD_PWSSC.P.1) + sum(HYD_ADFG.P.1) + sum(EGG.P.1) + sum(MDM.P.1) + sum(SeAC.P.1) + sum(SpAC.P.1)
    
  #if(loss.function=="squared.error"){
    HYD_PWSSC.P <- apply(HYD_PWSSC.PP[,model.data$hydPWSSC!=-9],2,FUN=function(X) sum((X-mean(X))^2)/length(X))
    HYD_ADFG.P  <- apply(HYD_ADFG.PP[,model.data$hydADFG!=-9],2,FUN=function(X) sum((X-mean(X))^2)/length(X))
    EGG.P       <- apply(EGG.PP[,model.data$egg!=-9],2,FUN=function(X) sum((X-mean(X))^2)/length(X))
    MDM.P       <- apply(MDM.PP[,model.data$mdm!=-9],2,FUN=function(X) sum((X-mean(X))^2)/length(X))
    SeAC.P      <- apply(PPSeAC[,DATA.SeAC>0],2,FUN=function(X) sum((X-mean(X))^2)/length(X))
    SpAC.P      <- apply(PPSpAC[,DATA.SpAC>0],2,FUN=function(X) sum((X-mean(X))^2)/length(X))
    
    HYD_PWSSC.P.1 <- HYD_PWSSC.P/mean(model.data$hydPWSSC[model.data$hydPWSSC!=-9])^2
    HYD_ADFG.P.1  <- HYD_ADFG.P/mean(model.data$hydADFG[model.data$hydADFG!=-9])^2
    EGG.P.1       <- EGG.P/mean(model.data$egg[model.data$egg!=-9])^2
    MDM.P.1       <- MDM.P/mean(model.data$mdm[model.data$mdm!=-9])^2
    SeAC.P.1      <- SeAC.P/seine_age_means^2
    SpAC.P.1      <- SpAC.P/spawn_age_means^2
    
  P.data.SquErr <- c(sum(HYD_PWSSC.P),sum(HYD_ADFG.P),sum(EGG.P),sum(MDM.P),sum(SeAC.P),sum(SpAC.P))
  P.data.NormSquErr <- c(sum(HYD_PWSSC.P.1),sum(HYD_ADFG.P.1),sum(EGG.P.1),sum(MDM.P.1),sum(SeAC.P.1),sum(SpAC.P.1))
  P.SquErr <- sum(HYD_PWSSC.P.1) + sum(HYD_ADFG.P.1) + sum(EGG.P.1) + sum(MDM.P.1) + sum(SeAC.P.1) + sum(SpAC.P.1)
  
  # Dinf.AbsErr <- P.AbsErr+G.AbsErr
  Dinf.SquErr <- P.SquErr+G.SquErr
  
  posterior.predictions <- function(obs,est,data.set.name,units,age.comp=FALSE){
    # Gonna have to create data frames
    # Columns:  Data set, Year, Obs, Lower 95th, Median, Upper 95th
    if(age.comp){
      out <- data.frame(data.type=data.set.name,
                        units=units,
                        year=rep(1980:(1980+nYr-1),each=nage),
                        age=rep(0:(nage-1),times=nYr),
                        obs=obs,
                        lower.95th=apply(est,2,quantile,probs=0.025,na.rm=TRUE)/100,
                        median = apply(est,2,median,na.rm=TRUE)/100,
                        upper.95th=apply(est,2,quantile,probs=0.975,na.rm=TRUE)/100) 
    }else{
      out <- data.frame(data.type=data.set.name,
                        units=units,
                        year=1980:(1980+length(apply(est,2,median))-1),
                        age=NA,
                        obs=obs,
                        lower.95th=apply(est,2,quantile,probs=0.025,na.rm=TRUE),
                        median = apply(est,2,median,na.rm=TRUE),
                        upper.95th=apply(est,2,quantile,probs=0.975,na.rm=TRUE)) 
    }
    
    return(out)
  }
  
  A1 <- posterior.predictions(obs=model.data$hydPWSSC,est=HYD_PWSSC.PP,data.set.name="PWSSC hydroacoustic",units="tons")
  B1 <- posterior.predictions(obs=model.data$hydADFG,est=HYD_ADFG.PP,data.set.name="ADFG hydroacoustic",units="tons")
  C1 <- posterior.predictions(obs=model.data$egg,est=EGG.PP,data.set.name="Egg deposition",units="trillion")
  D1 <- posterior.predictions(obs=model.data$mdm,est=MDM.PP,data.set.name="Milt",units="mile-days")
  E1 <- posterior.predictions(obs=as.vector(t(model.data$seine)),est=PPSeAC,data.set.name="Seine fishery",units="proportion",age.comp=TRUE)
  F1 <- posterior.predictions(obs=as.vector(t(model.data$spac)),est=PPSpAC,data.set.name="Spawner survey",units="proportion",age.comp=TRUE)
  
  # 05/23/2020:  3.8 minutes for a single model
  
  # return(list(waic=waic_2,ppl=c(Dinf.SquErr,G.SquErr,P.SquErr,G.data.NormSquErr,P.data.NormSquErr,G.data.SquErr,P.data.SquErr,
  #                               Dinf.AbsErr,G.AbsErr,P.AbsErr,G.data.NormAbsErr,P.data.NormAbsErr,G.data.AbsErr,P.data.AbsErr),
  #             model.fits=rbind(A1,B1,C1,D1,E1,F1)))
  return(list(waic=waic_2,
              ppl=c(Dinf.SquErr,G.SquErr,P.SquErr,G.data.NormSquErr,P.data.NormSquErr,G.data.SquErr,P.data.SquErr),
              model.fits=rbind(A1,B1,C1,D1,E1,F1)))
}




