import numpy as np
import copy
import time
from scipy import constants
import create_time_con
#import pandas as pd
#import pickle
#from matplotlib.colors import LogNorm

from constant import *
from inputsimobs import *
from inputbayes import *
import class_mcmc
import func_orbit,func_ref,func_bayes 
import plot

debug = True
#ddebug = True
ddebug = False


print(time.ctime())

# new time seq or not
new_time = False

# create new ref astrometry
create_sim_delta_dxdy = True #True

# ref star
i_tofit = 0

if new_time:
    ##########################
    time_con0 = np.linspace(t0,t1,int(N_time_init))
    np.savetxt("t_all0.dat", np.transpose([time_con0]))
    #
    retain_ind = create_time_con.create_time_seg()
    #
    time_con = time_con0[retain_ind]
    np.savetxt("t_all.dat", np.transpose([time_con]))
else:
    time_con = np.genfromtxt('t_all.dat')
    np.savetxt("t_all_clone.dat", np.transpose([time_con]))



time_dt = time_con[:-1]
if debug:
    np.savetxt('t_seq.out', np.transpose([time_dt]))


N_time = len(time_con)



#
if one_planet:
    as_cal_con = func_orbit.cal_t_radec(time_con,mp_Mearth,ms_Msun,a_AU,d_pc,e_orbit,periapsis_omega,cos_i_orbit, ascend_node_Omega, M0)
    as_mu = np.transpose(as_cal_con)
else:
    ####################################
    # two planet using dynamic-varied parameters generated from mercury6 package
    # add 2nd planet
    ##################
    # parameters
    # aei files from mercury6
    aei_file_name1 = "21.aei"
    aei_file_name2 = "22.aei"
    # M0 should be provided 
    M01 = 100
    M02 = 200
    ##################
    aei_file1 = np.genfromtxt(aei_file_name1)
    aei_file2 = np.genfromtxt(aei_file_name2)
    as_mu1 = func_orbit.cal_t_radec_fromAEI(M01, time_con,aei_file1,mp_Mearth,ms_Msun,d_pc)
    as_mu2 = func_orbit.cal_t_radec_fromAEI(M02, time_con,aei_file2,mp_Mearth,ms_Msun,d_pc)
    as_mu = as_mu1 + as_mu2 
    np.savetxt("as_mu.dat", np.transpose([as_mu[:,0],as_mu[:,1]]))
    np.savetxt("as_mu21.dat", np.transpose([as_mu1[:,0],as_mu1[:,1]]))
    np.savetxt("as_mu22.dat", np.transpose([as_mu2[:,0],as_mu2[:,1]]))

if debug:
    np.savetxt("as_mu.dat", np.transpose([as_mu[:,0],as_mu[:,1]]))

quit()

if create_sim_delta_dxdy:
    ##########################
    ##### simulation delta dx dy
    delta_dx_dy_pm = func_ref.create_delta_dxdy(time_con, n_ref_stars, w_ra_dec, as_mu, noise_mean, noise_std, ref_stars_pm_ra, ref_stars_pm_dec)
    if debug:
        for i in range(delta_dx_dy_pm.shape[1]):
            delta_dx_dy_name = "%s%s%s" % ("delta_dx_dy_", i, ".dat")
            np.savetxt(delta_dx_dy_name, np.transpose([delta_dx_dy_pm[:,i,0],delta_dx_dy_pm[:,i,1]]))
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
    

delta_dx_dy_sig = np.genfromtxt("delta_dx_dy_sig.dat")
np.savetxt("delta_dx_dy_sig.dat.clone", np.transpose([delta_dx_dy_sig[:,0], delta_dx_dy_sig[:,1]]))

##########################
## parallel tempering
i_dxdy_name =  "%s%s%s" % ("delta_dx_dy_", int(i_tofit), ".dat")
# genformtxt
delta_dx_dy_onerefstar = np.genfromtxt(i_dxdy_name)
if ddebug:
    np.savetxt("delta_dx_dy_one.dat", np.transpose([delta_dx_dy_onerefstar[:,0], delta_dx_dy_onerefstar[:,1]]))

chains_list = func_bayes.mcmc_PT(time_con, delta_dx_dy_onerefstar, delta_dx_dy_sig)

n_PT = len(beta_PT)
for i in range(n_PT):
    fname = "%s%s" % ("chain", i)
    chains_list[i].dump(fname)

#plot.plot_chains_list(chains_list, n_PT, time_dt)

print(time.ctime())
print("end")


