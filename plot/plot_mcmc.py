import numpy as np
import os
import matplotlib.pyplot as plt
from input import one_planet

def plot_mcmc(fname, burn_skip, bname): 
    chain  = np.genfromtxt(fname, skip_header=burn_skip)
    plt.figure(figsize=(15, 8))
    plt.subplot(241)
    plt.hist(chain[:,0], 100);
    plt.xlabel(r'$i$', fontsize=15)
    #plt.ylabel(r'$x$', fontsize=15)
    #xx = np.linspace(-3.5, 3.5, 1000)
    #plt.plot(xx, scipy.stats.norm(loc=0, scale=1).pdf(xx), lw=2)
    plt.subplot(242)
    plt.hist(chain[:,1], 100);
    plt.xlabel(r'$e$', fontsize=15)
    plt.subplot(243)
    plt.hist(chain[:,2], 100);
    plt.xlabel(r'$Omega$', fontsize=15)
    plt.subplot(244)
    plt.hist(chain[:,3], 100);
    plt.xlabel(r'$omega$', fontsize=15)
    plt.subplot(245)
    plt.hist(chain[:,4], 100);
    plt.xlabel(r'$M0$', fontsize=15)
    plt.subplot(246)
    plt.hist(chain[:,5], 100);
    plt.xlabel(r'$mp$', fontsize=15)
    plt.subplot(247)
    plt.hist(chain[:,6], 100);
    plt.xlabel(r'$var-uke$', fontsize=15)
    plt.subplot(248)
    plt.hist(chain[:,7], 100);
    plt.xlabel(r'$period$', fontsize=15)
    plt.tight_layout()
    figname = "%s%s" % (fname, ".png")
    fname = "%s%s%s" % (i, " ", bname)
    plt.suptitle(fname)
    plt.savefig(figname)

def plot_mcmc_2p(fname, burn_skip, bname): 
    chain  = np.genfromtxt(fname, skip_header=burn_skip)
    plt.figure(figsize=(15, 8))
    plt.subplot(331)
    plt.hist(chain[:,0], 100);
    plt.xlabel(r'$i$', fontsize=15)
    #plt.ylabel(r'$x$', fontsize=15)
    #xx = np.linspace(-3.5, 3.5, 1000)
    #plt.plot(xx, scipy.stats.norm(loc=0, scale=1).pdf(xx), lw=2)
    plt.subplot(332)
    plt.hist(chain[:,1], 100);
    plt.xlabel(r'$e$', fontsize=15)
    plt.subplot(333)
    plt.hist(chain[:,5], 100);
    plt.xlabel(r'$mp$', fontsize=15)
    plt.subplot(334)
    plt.hist(chain[:,6], 100);
    plt.xlabel(r'$var-uke$', fontsize=15)
    plt.subplot(335)
    plt.hist(chain[:,7], 100);
    plt.xlabel(r'$period$', fontsize=15)
    plt.subplot(336)
    plt.hist(chain[:,8], 100);
    plt.xlabel(r'$i2$', fontsize=15)
    plt.subplot(337)
    plt.hist(chain[:,9], 100);
    plt.xlabel(r'$e2$', fontsize=15)
    plt.subplot(338)
    plt.hist(chain[:,13], 100);
    plt.xlabel(r'$mp2$', fontsize=15)
    plt.subplot(339)
    plt.hist(chain[:,14], 100);
    plt.xlabel(r'$period2$', fontsize=15)
    plt.tight_layout()
    figname = "%s%s" % (fname, ".png")
    fname = "%s%s%s" % (i, " ", bname)
    plt.suptitle(fname)
    plt.savefig(figname)

cwd = os.getcwd()
bname = os.path.basename(cwd)

burn_skip = 90000
for i in range(8):
    fname = "%s%s%s" % ("chain", i, ".ch")
    if one_planet:
        plot_mcmc(fname, burn_skip, bname)
    else:
        plot_mcmc_2p(fname, burn_skip, bname)


