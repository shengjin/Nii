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

time_con = np.genfromtxt('t_all.dat')
np.savetxt("t_all_2.dat", np.transpose([time_con]))



delta_dx_dy_sig = np.genfromtxt("delta_dx_dy_sig.dat")

delta_dx_dy_one = np.genfromtxt("delta_dx_dy_one.dat")
delta_dx_dy_simone = np.genfromtxt("sim_fit_dx_dy.dat")

if debug:
    np.savetxt("delta_dx_dy_sig_2.dat", np.transpose([delta_dx_dy_sig[:,0], delta_dx_dy_sig[:,1]]))


##########################
## parallel tempering
delta_dx_dy_onerefstar = delta_dx_dy_one
delta_dx_dy_onerefstar[:,0] = delta_dx_dy_one[:,0] - delta_dx_dy_simone[:,0] 
delta_dx_dy_onerefstar[:,1] = delta_dx_dy_one[:,1] - delta_dx_dy_simone[:,1] 
if debug:
    np.savetxt("delta_dx_dy_left2.dat", np.transpose([delta_dx_dy_onerefstar[:,0], delta_dx_dy_onerefstar[:,1]]))


chains_list = func_bayes.mcmc_PT(time_con, delta_dx_dy_onerefstar, delta_dx_dy_sig)

n_PT = len(beta_PT)
for i in range(n_PT):
    fname = "%s%s" % ("chain", i)
    chains_list[i].dump(fname)

#plot.plot_chains_list(chains_list, n_PT, time_dt)

print(time.ctime())
print("end")




