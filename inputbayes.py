from inputsimobs import one_planet
init_ratio_gp = 0.1

one_planet_fit = True
if one_planet_fit:
    n_dim = 8
else:
    n_dim = 15
# 
test_ar_start = 150000 # donot check ar before first * iteration
n_iter_tunesig = 20 # not big than 30, otherwise will dead! 
target_ar = 0.24
ar_min_crit = 0.16
ar_max_crit = 0.33
tune_ar = 0.1 # reduce the frequency of tuning
# 
ar_diff_crit = 0.1
#
para_tune_min = -1.2 # 0.10 # in log 
#para_tune_min = -0.7 # 0.20 # in log 
#para_tune_min = -0.3 # 0.50 # in log 
para_tune_max = 1.2 # 5.01 times # in log 
#para_tune_max = 0.3 # 1.99 times # in log 
scale_min = 0.00005
scale_max = 0.5
# reach bound, scale to (bound to bound(*or/)scale_bound)
scale_bound = 300

# iterations of mcmc
n_iter = 1000000
n_iter_swap_av = 25
# n_iter_swap = n_iter_swap_av +/- dn_swap
dn_swap = 5
#
N_swap = 4 # test swap 4 adjencent chains

beta_PT = [0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0]

###############################
#Prior setting
ecc_max = 1 
ecc_min = 0
#
## in degree
#incl_max = 180 
#incl_min = 0
cos_incl_max = 1
cos_incl_min = -1
#
# in degree
an_Omega_max = 360
an_Omega_min = 0
#
# in degree
p_omega_max = 360
p_omega_min = 0
#
# in degree
M0_max = 360
M0_min = 0
#
# in days
period_max = 3650
period_min = 0.5
#
# in Earth mass
p_mass_max = 3000
p_mass_min = 0.1
#
# in mu_as
var_uk_err_max = 100
var_uk_err_min = 0
var_uk_err_a = 1



