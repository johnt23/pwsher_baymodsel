# model_run_mortality_effects.r
# Created by John Trochta
# Date started:  12/01/2016
# Summary:
# This script runs the Bayesian ASA model (BASA) for PWS using different
# covariates and extracts the output. Each section calls run_models to perform
# the following steps:
#     1) Compiles and extracts each version of BASA from bin directory.
#     2) Runs parallel processing where each covariate is individually turned on
#        within BASA. Iterative reweighting of eff. sample sizes for age comps
#        is done before running NUTS.
#     3) Runs 3 MCMC chains with NUTS for each BASA-one covariate model.
#        Importantly, parallel processsing is setup so that the max number of
#        processors running at any time is n_cores x nuts.inputs$reps.
#     4) Stores output files of each BASA model in their respective folders
#        (following naming of hypothesis_names) and summary files of NUTS runs
#        (i.e. parameter posterior quantiles & convergence diagnostics) within
#        NUTS_chain_diagnostics

library(here)
setwd(here::here("src"))
library(data.table)
source(file="run_models.R")

### Specifications for all model sets
run_year="2021"
run_month="07"
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
                      'GOA Arrowtooth flounder total biomass (Ages 3+)',
                      "Summer humpback whale estimates (Age 3+)",
                      "Humpback whale counts (Age 3+)")
h.block <- c(1:14) # Indices of above hypotheses in columns (matching length of this vector) of covariate matrix in PWS_ASA(covariate).ctl
start.years <- c(1994,2007,1994,rep(1980,9),1981,2006)
include.null=FALSE
n_run_types=length(hypothesis_names)
n_age0_covs=8
n_mor_covs=14

n_cores <- 3
mcmc.algorithm="NUTS"
nuts.inputs=data.frame(reps=3,iter=3000,adapt_delta=0.925,warm_up=500,nuts.mode=1)

### Specifications for each model set

##############################################
# 2021_07_all fixed mortality effects
# Run with
#  - all mortality covariates as fixed effects from 1980-2017
#  - rec devs are freely estimated (not centered or bias corrected)
#  - mean mortality is fixed at 0.25 (line 461 in PWS_ASA.TPL)

H.select <- list(single=4:13)
n_models <- 1+length(H.select$single)
custom.settings <- list(n_run_types=n_run_types,
                        R.or.M=rep('M',n_run_types),
                        series_start_yr=2007, # This is the most recent start year of any time series
                        estimate.1993.mortality=rep(1,n_models),
                        estimate.mean.mortality=rep(0,n_models),
                        mean.mortality=c(rep(0.25,length(hypothesis_names))),
                        h.block=h.block,
                        sim_cov=0,
                        block.on=NULL,
                        run_ID=paste0(run_year,'_',run_month,'_','all fixed mortality effects'),
                        base.M="PWS_ASA_all fixed effects")

run_models(hypothesis_names,
           include.null,
           start.years,
           H.select,
           n_cores,
           n_sims,
           n_age0_covs,
           n_mor_covs,
           custom.settings,
           mcmc.algorithm,
           nuts.inputs,
           run_year,
           run_month)
# 2020-04-23:  277.6 mins.


##############################################
# 2021_07_all fixed mortality effects & start from 1994
# Run with
#  - all covariates as fixed effects
#  - rec devs are freely estimated (not centered or bias corrected)
#  - mean mortality is fixed at 0.25 (line 461 in PWS_ASA.TPL)

H.select <- list(single=4:13,Disease=1:3)
n_models <- 1+length(H.select$single)
custom.settings <- list(n_run_types=n_run_types,
                        R.or.M=rep('M',n_run_types),
                        series_start_yr=2007, # This is the most recent start year of any time series
                        estimate.1993.mortality=rep(1,length(n_models)),
                        estimate.mean.mortality=rep(0,length(n_models)),
                        mean.mortality=c(rep(0.25,length(hypothesis_names))),
                        h.block=h.block,
                        sim_cov=0,
                        block.on=NULL,
                        run_ID=paste0(run_year,'_',run_month,'_','all fixed mortality effects & start from 1994'),
                        base.M="PWS_ASA_all fixed effects _ start from 1994")

run_models(hypothesis_names,
           include.null,
           start.years,
           H.select,
           n_cores,
           n_sims,
           n_age0_covs,
           n_mor_covs,
           custom.settings,
           mcmc.algorithm,
           nuts.inputs,
           run_year,
           run_month)
# 2020-04-24:  200 mins.


##############################################
# 2021_07_latent mortality effects & year specific SE
# Run with
#  - select covariates are latent variables
#  - rec devs are freely estimated (not centered or bias corrected)
#  - mean mortality is fixed at 0.25 (line 461 in PWS_ASA.TPL)
#  - Adapted from 2020_01_only latent mortality effects in years with obs, 
#  - CV of effects taken from assessment reports, not estimated

H.select <- list(single=4:14,Disease=1:3)
n_models <- 1+length(H.select$single)
# H.select <- list(single=4:7)
custom.settings <- list(n_run_types=n_run_types,
                        R.or.M=rep('M',n_run_types),
                        series_start_yr=2007, # This is the most recent start year of any time series
                        estimate.1993.mortality=c(rep(0,n_models-2),1,1),
                        estimate.mean.mortality=rep(0,n_models),
                        mean.mortality=rep(0.25,n_models),
                        h.block=h.block,
                        sim_cov=0,
                        block.on=NULL,
                        run_ID=paste0(run_year,'_',run_month,'_','latent mortality effects _ year specific SE'),
                        base.M="PWS_ASA_latent vars for mortality _ year specific SE")

run_models(hypothesis_names,
           include.null,
           start.years,
           H.select,
           n_cores,
           n_sims,
           n_age0_covs,
           n_mor_covs,
           custom.settings,
           mcmc.algorithm,
           nuts.inputs,
           run_year,
           run_month)
# 2021-07-06:  323 minutes
