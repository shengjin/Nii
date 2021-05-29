##################
# plot  delta_dx_dy of simulation signale (ref star mean) in a two-panel fig

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from cycler import cycler


font = {'family' : 'serif', #monospace
        'weight' : 'normal', #bold
        'size'   : 14,
       }

mpl.rc('font', **font)


#plt.style.use("bmh")
# Later in the code
plt.style.use('ggplot')


def plot_dx_dy_smltmean(dt_con, delta_dx_dy, figname):
    #import matplotlib
    import matplotlib.pyplot as plt
    from constant import year
    plt.clf()
    plt.subplot(2,1,1) 
    plt.errorbar(dt_con/year, delta_dx_dy[:,0], c='k', yerr=1, linestyle='', linewidth=0.6)
    plt.scatter(dt_con/year, delta_dx_dy[:,0], s=5, c='k', label = "simulation mean")
    plt.ylabel('d_ra (mu arcsec)')
    plt.subplot(2,1,2) 
    plt.errorbar(dt_con/year, delta_dx_dy[:,1], c='k', yerr=1, linestyle='', linewidth=0.6 )
    plt.scatter(dt_con/year, delta_dx_dy[:,1], s=5, c='k', label = "simulation mean")
    plt.ylabel('d_dec (mu arcsec)')
    plt.xlabel('dt [years]')
    plt.legend(loc='best')
    plt.savefig(figname, dpi=300)

# plot  delta_dx_dy of simulation signale (ref star mean) in a two-panel fig
def plot_dx_dy_smltmean_fit(dt_con, delta_dx_dy, sim_dx_dy, figname):
    #import matplotlib
    import matplotlib.pyplot as plt
    from constant import year
    plt.clf()
    plt.subplot(2,1,1) 
    plt.minorticks_on()
    plt.xlim(-0.0,5.01)
    plt.ylim(-4.0,8.0)
    plt.errorbar(dt_con/year, delta_dx_dy[:,0], c='k', yerr=1, linestyle='', linewidth=0.6)
    plt.plot(dt_con/year, sim_dx_dy[:,0], c='r', label = "Fitted Orbit Parameters")
    plt.scatter(dt_con/year, delta_dx_dy[:,0], s=5, c='k', label = "Simulated Observation")
    #plt.ylabel('d_ra (mu arcsec)')
    #plt.grid(False)
    plt.ylabel(r'$\Delta{\mathrm{ra}}$ $(\mu \mathrm{as})$') #, size=17)
    #plt.legend(loc=4, prop={'size': 7})
    plt.legend(loc=4, ncol=2, prop={'size': 8})
    plt.subplot(2,1,2) 
    plt.minorticks_on()
    plt.xlim(-0.0,5.01)
    plt.ylim(-8.0,4.0)
    plt.errorbar(dt_con/year, delta_dx_dy[:,1], c='k', yerr=1, linestyle='', linewidth=0.6 )
    plt.plot(dt_con/year, sim_dx_dy[:,1], c='r', label = "Fitted Orbit Parameters")
    plt.scatter(dt_con/year, delta_dx_dy[:,1], s=5, c='k', label = "Simulated Observation")
    #plt.ylabel('d_dec (mu arcsec)')
    plt.ylabel(r'$\Delta{\mathrm{dec}}$ $(\mu \mathrm{as})$') #, size=17)
    #plt.xlabel('dt [years]')
    #plt.xlabel('dt  [ years ]') #, size=17)
    plt.xlabel(r'$\mathrm{dt}$ $[ \mathrm{years} ]$') #, size=17)
    plt.legend(loc=4, ncol=2, prop={'size': 8})
    #plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    #plt.legend(loc=2)
    plt.savefig(figname, dpi=300)

def plot_dx_dy_mclast(dt_con, delta_dx_dy, figname, lab):
    #import matplotlib
    import matplotlib.pyplot as plt
    from constant import year
    plt.clf()
    plt.subplot(2,1,1) 
    plt.plot(dt_con/year, delta_dx_dy[:,0], label = lab)
    plt.ylabel('d_ra (mu arcsec)')
    plt.subplot(2,1,2) 
    plt.plot(dt_con/year, delta_dx_dy[:,1], label = lab)
    plt.ylabel('d_dec (mu arcsec)')
    plt.xlabel('dt [years]')
    plt.legend(loc='best')
    plt.savefig(figname, dpi=300)


# plot  chain
def plot_chains_list(chains_list, n_PT, time_con):
    #import matplotlib
    import matplotlib.pyplot as plt
    import func_orbit
    plt.clf()
    figname = "lnprob.png"
    for i in range(n_PT):
        #plt.subplot(n_PT+1,1,i+1) 
        plt.plot(chains_list[i].lnprob, label =i)
        plt.legend(loc='best')
    plt.ylabel('lnprob')
    plt.yscale('symlog')
    plt.xlabel('n')
    plt.savefig(figname, dpi=600)
    #
    plt.clf()
    figname = "ar.png"
    for i in range(n_PT):
        #plt.subplot(n_PT+1,1,i+1) 
        plt.plot(chains_list[i].accept_rate , label =i)
        plt.legend(loc='best')
    plt.ylabel('ar')
    plt.xlabel('n')
    plt.savefig(figname, dpi=600)
    #
    plt.clf()
    figname = "i.png"
    for i in range(n_PT):
        #plt.subplot(n_PT+1,1,i+1) 
        plt.plot(chains_list[i].chain[:,0] , label =i)
        plt.legend(loc='best')
    plt.ylabel('i')
    plt.xlabel('n')
    plt.savefig(figname, dpi=600)
    #
    plt.clf()
    figname = "e.png"
    for i in range(n_PT):
        #plt.subplot(n_PT+1,1,i+1) 
        plt.plot(chains_list[i].chain[:,1] , label =i)
        plt.legend(loc='best')
    plt.ylabel('e')
    plt.xlabel('n')
    plt.savefig(figname, dpi=600)
    #
    plt.clf()
    figname = "an0.png"
    for i in range(n_PT):
        #plt.subplot(n_PT+1,1,i+1) 
        plt.plot(chains_list[i].chain[:,2] , label =i)
        plt.legend(loc='best')
    plt.ylabel('an0')
    plt.xlabel('n')
    plt.savefig(figname, dpi=600)
    #
    plt.clf()
    figname = "po.png"
    for i in range(n_PT):
        #plt.subplot(n_PT+1,1,i+1) 
        plt.plot(chains_list[i].chain[:,3] , label =i)
        plt.legend(loc='best')
    plt.ylabel('po')
    plt.xlabel('n')
    plt.savefig(figname, dpi=600)
    #
    plt.clf()
    figname = "M0.png"
    for i in range(n_PT):
        #plt.subplot(n_PT+1,1,i+1) 
        plt.plot(chains_list[i].chain[:,4] , label =i)
        plt.legend(loc='best')
    plt.ylabel('M0')
    plt.xlabel('n')
    plt.savefig(figname, dpi=600)
    #
    plt.clf()
    figname = "mp.png"
    for i in range(n_PT):
        #plt.subplot(n_PT+1,1,i+1) 
        plt.plot(chains_list[i].chain[:,5] , label =i)
        plt.legend(loc='best')
    plt.ylabel('mp')
    plt.yscale('symlog')
    plt.xlabel('n')
    plt.savefig(figname, dpi=600)
    #
    plt.clf()
    figname = "var_uke.png"
    for i in range(n_PT):
        #plt.subplot(n_PT+1,1,i+1) 
        plt.plot(chains_list[i].chain[:,6] , label =i)
        plt.legend(loc='best')
    plt.ylabel('var_uke')
    plt.xlabel('n')
    plt.savefig(figname, dpi=600)
    #
    plt.clf()
    figname = "period.png"
    for i in range(n_PT):
        #plt.subplot(n_PT+1,1,i+1) 
        plt.plot(chains_list[i].chain[:,7] , label =i)
        plt.legend(loc='best')
    plt.ylabel('period')
    plt.yscale('symlog')
    plt.xlabel('n')
    plt.savefig(figname, dpi=600)
    #
    time_dt = time_con[:-1]
    for i in range(n_PT):
        ch_last = chains_list[i].chain[-1,:]
        dx_dy_mc = func_orbit.gen_dx_dy_mc(ch_last, time_con)
        figname = "%s%s" % ("delta_dx_dy_mc", i)
        plot_dx_dy_mclast(time_dt, dx_dy_mc, figname, i)
    
