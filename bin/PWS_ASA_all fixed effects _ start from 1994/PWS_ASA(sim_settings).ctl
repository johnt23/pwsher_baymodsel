# Specifies settings for running simulations with the PWS BASA model

# sim_catches - should catches be simulated by randomly resampling past years of estimated vulnerability (1) or not (0)
1

# nyr_resample_period
5

# resample_period - vector of indices specifying which time frame from which to resample for simulating catches
5
6
7
8
9

# exploitation_history
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1
0.1

# age0_dev_option - 0 if use devs specified in PIN file or 1 to randomly sample devs given mean and SD for recruitment
1

# data_avg_option - 0 if use historical values of all other data sets not fit (fixed in the model w/out error) from PWS_ASA.dat or 1 to calculate average of entire time frame to use for simulations
# This includes 3 data sets:  empirical weight-at-age, fecundity-at-age, % of female spawners
0
