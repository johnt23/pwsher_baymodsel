# model_run_recruit_effects.r
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
hypothesis_names <- c("1989 regime shift",
                      "Hatchery-released juvenile pink salmon",
                      "GOA walleye pollock age 1",
                      'Summer NPGO',
                      'Summer PDO',
                      'Age 0 scale growth (lagged 1 year)',
                      'Freshwater discharge',
                      "PWS Avg. Zooplankton density '81-'99")
h.block <- c(1:8) # Indices of above hypotheses in columns (matching length of this vector) of covariate matrix in PWS_ASA(covariate).ctl
start.years <- c(rep(1980,7),1981)
include.null=FALSE
n_run_types=length(hypothesis_names)
n_age0_covs=8
n_mor_covs=14

n_cores <- 3
mcmc.algorithm="NUTS"
nuts.inputs=data.frame(reps=3,iter=3000,adapt_delta=0.925,warm_up=500,nuts.mode=1)

#### 2020_04_all fixed age0 effects ####
# Run with
#  - all covariates as fixed effects from 1980-2017
#  - rec devs are centered around estimated mean & estimated sigma
#  - mean mortality is fixed at 0.25 (line 461 in PWS_ASA.TPL)

H.select <- list(single=1:length(hypothesis_names))
custom.settings <- list(n_run_types=n_run_types,
                        R.or.M=rep('R',n_run_types),
                        series_start_yr=2000,
                        estimate.1993.mortality=rep(1,length(hypothesis_names)),
                        estimate.mean.mortality=rep(0,length(hypothesis_names)),
                        mean.mortality=c(rep(0.3,length(hypothesis_names))),
                        h.block=h.block,
                        sim_cov=0,
                        block.on=NULL,
                        run_ID='2021_07_all fixed age0 effects',
                        base.M="PWS_ASA_all fixed effects _ pen rec")

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

# 2021-07-14: 183 mins.  


#### 2020_04_all fixed age0 effects & start from 1994 ####
# Run with
#  - all covariates as fixed effects
#  - only 1994-2017 of covariates are incorporated
#  - rec devs are centered around estimated mean & estimated sigma
#  - mean mortality is fixed at 0.25 (line 461 in PWS_ASA.TPL)

H.select <- list(single=1:7)
custom.settings <- list(n_run_types=n_run_types,
                        R.or.M=rep('R',n_run_types),
                        series_start_yr=1981,
                        estimate.1993.mortality=rep(1,length(hypothesis_names)),
                        estimate.mean.mortality=rep(0,length(hypothesis_names)),
                        mean.mortality=c(rep(0.25,length(hypothesis_names))),
                        h.block=h.block,
                        sim_cov=0,
                        block.on=NULL,
                        run_ID='2020_04_all fixed age0 effects _ start from 1994',
                        base.M="PWS_ASA_all fixed effects & pen rec _ start from 1994")

start_time <- proc.time()
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
print(proc.time() - start_time)
# 2020-01-25:  120 mins.


#### 2021_07_latent age0 effects & year specific SE ####
# Run with
#  - select covariates are latent variables
#  - rec devs are centered around estimated mean & estimated sigma
#  - mean mortality is fixed at 0.25 (line 461 in PWS_ASA.TPL)
#    This version ensures estimated recruit devs for those years without obs do not impact mortality

run_year="2021"
run_month="07"

H.select <- list(single=2:7)
custom.settings <- list(n_run_types=n_run_types,
                        R.or.M=rep('R',n_run_types),
                        series_start_yr=1980,
                        estimate.1993.mortality=rep(1,length(hypothesis_names)),
                        estimate.mean.mortality=rep(0,length(hypothesis_names)),
                        mean.mortality=c(rep(0.25,length(hypothesis_names))),
                        h.block=h.block,
                        sim_cov=0,
                        block.on=NULL,
                        run_ID=paste0(run_year,'_',run_month,'_','latent age0 effects _ year specific SE'),
                        base.M="PWS_ASA_latent vars for age0 _ year specific SE")

nuts.inputs=data.frame(reps=3,iter=3000,adapt_delta=0.95,warm_up=500,nuts.mode=1)

start_time <- proc.time()
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
print(proc.time() - start_time)

# 07/09/2021:  64.9 mins


#### 2021_07_all fixed age0 effects_low sigmaR ####
H.select <- list(single=1:length(hypothesis_names))
include.null=FALSE
n_run_types=length(hypothesis_names)
n_age0_covs=8
n_mor_covs=14
custom.settings <- list(n_run_types=n_run_types,
                        R.or.M=rep('R',n_run_types),
                        series_start_yr=2000,
                        estimate.1993.mortality=rep(1,length(hypothesis_names)),
                        estimate.mean.mortality=rep(0,length(hypothesis_names)),
                        mean.mortality=c(rep(0.3,length(hypothesis_names))),
                        h.block=h.block,
                        sim_cov=0,
                        block.on=NULL,
                        run_ID='2021_07_all fixed age0 effects_low sigmaR',
                        base.M="PWS_ASA_all fixed effects _ pen rec_low sigmaR")

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
# 07/14/2021 Run: ### minutes


