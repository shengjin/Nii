# coding: utf-8
import numpy as np
from input import ms_Msun,d_pc
import func_orbit  as fobt

au = fobt.period_2_au( 542.8360079168133 , ms_Msun)
print(au)
au = fobt.period_2_au( 882.9453243126169, ms_Msun)
print(au)
quit()
skipnum = 200000
skipnum_agl = 499990

chdat = np.genfromtxt("chain7.ch", skip_header=skipnum)
chdat.shape

x = np.zeros(8, dtype = float)
#print(x)
for i in range(chdat.shape[1]):
    print(i)
    x[i] = chdat[:,i].mean()
    print(i, chdat[:,i].mean())

chdat = np.genfromtxt("chain7.ch", skip_header=skipnum_agl)
x[2] = chdat[:,2].mean()
x[3] = chdat[:,3].mean()
x[4] = chdat[:,4].mean()
print(x[2])
print(x[3])
print(x[4])

time_con = np.genfromtxt('t_all.dat')
au = fobt.period_2_au(x[7], ms_Msun)
period = fobt.au_2_period(1.3, ms_Msun)
print(period)


sim_dx_dy = fobt.gen_dx_dy_mc(x, time_con)
np.savetxt("sim_fit_dx_dy.dat", np.transpose([sim_dx_dy[:,0], sim_dx_dy[:,1]]))


