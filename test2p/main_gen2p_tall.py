import numpy as np
import copy
import time
from scipy import constants
#import pandas as pd
#import pickle
#from matplotlib.colors import LogNorm

from constant import *
from input import *
from inputbayes import *
import class_mcmc
import func_orbit,func_ref,func_bayes 
import plot

debug = True
#ddebug = True
ddebug = False


print(time.ctime())

##########################
time_con = np.genfromtxt('t_all_2.dat')

################################
################################
##################################

time_dt = time_con[:-1]


##########################
# astrometry sequence
# as_cal_con shape: ra dec
as_cal_con = func_orbit.cal_t_radec(time_con,mp_Mearth,ms_Msun,a_AU,d_pc,e_orbit,periapsis_omega,i_orbit, ascend_node_Omega, M0)
# add 2nd planet
as_cal_con2 = func_orbit.cal_t_radec(time_con,mp_Mearth2,ms_Msun,a_AU2,d_pc,e_orbit2,periapsis_omega2,i_orbit2, ascend_node_Omega2, M02)
#as_mu = np.transpose(as_cal_con)
as_mu1 = np.transpose(as_cal_con)
as_mu2 = np.transpose(as_cal_con2)
as_mu = as_mu1 + as_mu2 

np.savetxt("as_mu_sim1.dat", np.transpose([as_mu1[:,0],as_mu1[:,1]]))
np.savetxt("as_mu_sim2.dat", np.transpose([as_mu2[:,0],as_mu2[:,1]]))
if debug:
    np.savetxt("as_mu_sim.dat", np.transpose([as_mu[:,0],as_mu[:,1]]))

delta_dx_dy = func_ref.create_delta_dxdy_pureAsMu(time_con, as_mu)
np.savetxt("delta_dx_dy_sim.dat", np.transpose([delta_dx_dy[:,0],delta_dx_dy[:,1]]))
quit()

##########################
##### simulation delta dx dy
delta_dx_dy_pm = func_ref.create_delta_dxdy(time_con, n_ref_stars, w_ra_dec, as_mu, noise_mean, noise_std, ref_stars_pm_ra, ref_stars_pm_dec)
if ddebug:
    np.savetxt("delta_dx_pm.dat", np.transpose([delta_dx_dy_pm[:,0,0],delta_dx_dy_pm[:,1,0]]))
    np.savetxt("delta_dy_pm.dat", np.transpose([delta_dx_dy_pm[:,0,1],delta_dx_dy_pm[:,1,1]]))

##### calc mean
delta_dx_dy_mean = np.zeros((int(N_time)-1, 2), dtype=np.float64)
delta_dx_dy_sig = np.zeros((int(N_time)-1, 2), dtype=np.float64)
for i in range(int(N_time)-1):
    delta_dx_dy_mean[i,0] = delta_dx_dy_pm[i,:,0].mean()
    delta_dx_dy_mean[i,1] = delta_dx_dy_pm[i,:,1].mean()
    delta_dx_dy_sig[i,0] = np.std(delta_dx_dy_pm[i,:,0])
    delta_dx_dy_sig[i,1] = np.std(delta_dx_dy_pm[i,:,1])
if debug:
    np.savetxt("delta_dx_dy_mean.dat", np.transpose([delta_dx_dy_mean[:,0], delta_dx_dy_mean[:,1]]))
    np.savetxt("delta_dx_dy_sig.dat", np.transpose([delta_dx_dy_sig[:,0], delta_dx_dy_sig[:,1]]))

#figname = 'dt_dx_dy.png'
#plot.plot_dx_dy_smltmean(time_dt, delta_dx_dy_mean, figname)


##########################
## parallel tempering
delta_dx_dy_onerefstar = delta_dx_dy_pm[:,0,:]
if debug:                                                                                      
    np.savetxt("delta_dx_dy_one.dat", np.transpose([delta_dx_dy_onerefstar[:,0], delta_dx_dy_onerefstar[:,1]]))
chains_list = func_bayes.mcmc_PT(time_con, delta_dx_dy_onerefstar, delta_dx_dy_sig)

n_PT = len(beta_PT)
for i in range(n_PT):
    fname = "%s%s" % ("chain", i)
    chains_list[i].dump(fname)

#plot.plot_chains_list(chains_list, n_PT, time_dt)

print(time.ctime())
print("end")




