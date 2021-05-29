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


time_dt = np.genfromtxt('t_all.dat')
print(len(time_dt))
#time_dt = np.genfromtxt('t_seq.out')
time_dt = time_dt[1:]
print(len(time_dt))
delta_dx_dy_mean = np.genfromtxt("delta_dx_dy_mean.dat")
sim_dx_dy = np.genfromtxt("sim_fit_dx_dy.dat")

figname = 'dt_dx_dy.png'
plot.plot_dx_dy_smltmean(time_dt, delta_dx_dy_mean, figname)
#plot.plot_dx_dy_smltmean_fit(time_dt, delta_dx_dy_mean, sim_dx_dy, figname)

n_PT = len(beta_PT)
chains_list = {x: class_mcmc.Chains(n_iter,n_dim) for x in range(n_PT)}

for i in range(n_PT):
    fname = "%s%s" % ("chain", i)
    fname_ch = "%s%s" % (fname, ".ch")
    fname_ar = "%s%s" % (fname, ".ar")
    fname_lp = "%s%s" % (fname, ".lp")
    chains_list[i].chain = np.genfromtxt(fname_ch)
    chains_list[i].lnprob = np.genfromtxt(fname_lp)
    chains_list[i].accept_rate = np.genfromtxt(fname_ar)

plot.plot_chains_list(chains_list, n_PT, time_dt)

