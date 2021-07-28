# PWS_ASA parameter phases

# Determines which ctl file for the age comps and ESS to use (1 uses ESS control to be iteratively estimated)
-1
 
# Ages 0-8 mortality phase (ph_Z_0_8)
-1
 
# Plus age group mortality phase (ph_Z_9)
2
 
# Ages 0-8 mortality offset phase (ph_Z_0_8offset)
-1
 
# Plus age group mortality offset phase (ph_Z_9offset)
-2
 
# Prop. of prop. age 4 mature that represents age 3 maturity - First period (1994 or 1997)
3
 
# Prop. age 4 mature - First period (1994 or 1997)
3
 
# Prop. of prop. age 4 mature that represents age 3 maturity - Second period (1994 or 1997)
-3
 
# Prop. age 4 mature - Second period (1994 or 1997)
-3
 
# Selectivity (alpha; ph_alpha_v)
4
 
# Selectivity (beta; ph_beta_v)
4
 
# Age-0 deviates/recruitment (ph_age0devs)
1
 
# Initial population (ph_init_pop) ages 1-5
1
 
# ph_eggAdd - This one should always be turned off
-2
 
# Mile-days coefficient (ph_mdm)
2
 
# ph_mdmAdd
5
 
# ADF&G Hydroacoustic coefficient (ph_hyd1)
2
 
# ADF&G Additional Hyd variance term (ph_hydAdd1) - I use this phase for all additional variance terms - they should be in last phase; are needed to get better est's of uncertainty -not needed for model fits.
5
 
# PWSSC Hydroacoustic coefficient (ph_hyd2)
2
 
# PWSSC Additional Hyd variance term (ph_hydAdd2)
5
 
# Additional mortality due to VHSV in 1992-1993 (Ages3-4)
2
 
# Additional mortality due to I. hoferi in 1992-1993 (Ages5-8)
2
 
# Mean age 0 over entire time series (ph_meanage0)
1
 
# Offset from Mean age 0 (ph_meanage0_offset)
-2
 
# Sigma for mean age 0 to allow for shrinkage (ph_sigmaage0) - set to -11 fixes sigma to 0 so rec_devs are freely estimated
-4
 
# Beta on covariate effects on age 0 values (ph_betaage0)
1
 
# Beta on covariate effects on mortality for ages 1-9+ (ph_betamortality)
-4
 
# Deviates on mortality covariate effects if included (ph_mortdevs)
-4
 
# Offset to the covariate effect on age 0 (ph_betaage0_offset)
-1
 
# Offset to the covariate effect on mortality (ph_betamortatliy_offset)
-1
 
# Survey Selectivity parameter 1 in logistic (alpha; ph_survey_vul_alpha)
-4
 
# Survey Selectivity parameter 2 in logistic (beta; ph_survey_vul_beta)
-4
 
# Sigma for penalty on annual mortality deviates (sigma_mortdevs)
-5
 
# Sigma for normal priors with recruitment covariates (sigma_age0covar)
-5
 
# Sigma for normal priors with mortality covariates (sigma_mortdevs)
-5
 
