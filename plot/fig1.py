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


plt.style.use("bmh")
# Later in the code
#plt.style.use('ggplot')
#plt.style.use("Solarize_Light2")

chfile1 = "chain7_1.ch"
chfile2 = "chain7.ch"
ini_burn = 200000
figname = "tmp.png"



dat_all1 = pd.read_table(chfile1, sep="\s+", skiprows=ini_burn)
dat_all2 = pd.read_table(chfile2, sep="\s+", skiprows=ini_burn)

fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(18, 12), dpi=160, sharex=False, sharey=False)
fig.subplots_adjust(hspace=0.30, wspace=0.05)

# i
dat_plt = dat_all1.iloc[:,0]
dat_plt.plot(ax=axes[0,0], kind = "hist", density = True, alpha=0.65, bins = 30) # change density to true, because KDE uses density
dat_plt.plot(ax=axes[0,0], kind = "kde") # Plot KDE
#axes[0,0].axvline(ii, c='w', alpha=0.65, linestyle = ":")
#axes[0,0].set_xlim(0, 20)
axes[0,0].set_xlim(101, 109)
axes[0,0].set_xlabel("$\mathrm{1st}$ $\mathrm{Planet:}$ i $(\mathrm{degree})$", style='italic') #, size=17)
axes[0,0].set_yticks([])
axes[0,0].set_ylabel("Marginal Posterior", labelpad=1.0)
axes[0,0].yaxis.set_label_coords(-0.06,0.5)
# bmh style add grid as default
axes[0,0].grid(False)
# Remove ticks and spines
axes[0,0].tick_params(left = False, bottom = False)
for axes[0,0], spine in axes[0,0].spines.items():
    spine.set_visible(False)

# e
dat_plt = dat_all1.iloc[:,1]
dat_plt.plot(ax=axes[0,1], kind = "hist", density = True, alpha=0.65, bins = 30) # change density to true, because KDE uses density
dat_plt.plot(ax=axes[0,1], kind = "kde") # Plot KDE
#axes[0,1].axvline(ee, c='w', alpha=0.65, linestyle = ":")
#axes[0,1].set_xlim(0, 1)
axes[0,1].set_xlim(0.06, 0.24)
axes[0,1].set_xlabel("$\mathrm{1st}$ $\mathrm{Planet:}$ e", style='italic') #, size=17)
        
axes[0,1].set_yticks([])
#axes[0,1].set_ylabel("Marginal Posterior")
axes[0,1].set_ylabel("")
# bmh style add grid as default
axes[0,1].grid(False)
# Remove ticks and spines
axes[0,1].tick_params(left = False, bottom = False)
for axes[0,1], spine in axes[0,1].spines.items():
    spine.set_visible(False)


# Omega
dat_plt = dat_all1.iloc[:,5]
dat_plt.plot(ax=axes[0,2], kind = "hist", density = True, alpha=0.65, bins = 30) # change density to true, because KDE uses density
dat_plt.plot(ax=axes[0,2], kind = "kde") # Plot KDE
#axes[0,2].axvline(omO, c='w', alpha=0.65, linestyle = ":")
#axes[0,1].set_xlim(0, 1)
axes[0,2].set_xlabel(r'$\mathrm{1st}$ $\mathrm{Planet:}$ $M_\mathrm{p}$ $(M_{\oplus})$') #, size=17)
axes[0,2].set_yticks([])
#axes[0,2].set_ylabel("Marginal Posterior")
axes[0,2].set_ylabel("")
#axes[0,2].set_xlim(0, 360)
# bmh style add grid as default
axes[0,2].grid(False)
# Remove ticks and spines
axes[0,2].tick_params(left = False, bottom = False)
for axes[0,2], spine in axes[0,2].spines.items():
    spine.set_visible(False)


# omega
dat_plt = dat_all1.iloc[:,7]
dat_plt.plot(ax=axes[1,0], kind = "hist", density = True, alpha=0.65, bins = 30) # change density to true, because KDE uses density
dat_plt.plot(ax=axes[1,0], kind = "kde") # Plot KDE
#axes[0,3].axvline(om, c='w', alpha=0.65, linestyle = ":")
#axes[0,1].set_xlim(0, 1)
#axes[1,0].set_xlim(0, 360)
axes[1,0].set_xlabel("$\mathrm{1st}$ $\mathrm{Planet:}$ P $(\mathrm{days})$", style='italic') #, size=17)
axes[1,0].set_yticks([])
axes[1,0].set_ylabel("Marginal Posterior", labelpad=1.0)
axes[1,0].yaxis.set_label_coords(-0.06,0.5)
#axes[0,3].set_ylabel("Marginal Posterior")
# bmh style add grid as default
axes[1,0].grid(False)
# Remove ticks and spines
axes[1,0].tick_params(left = False, bottom = False)
for axes[1,0], spine in axes[1,0].spines.items():
    spine.set_visible(False)

# M0
dat_plt = dat_all2.iloc[:,0]
dat_plt.plot(ax=axes[1,1], kind = "hist", density = True, alpha=0.65, bins = 30) # change density to true, because KDE uses density
dat_plt.plot(ax=axes[1,1], kind = "kde") # Plot KDE
#axes[1,0].axvline(m0, c='w', alpha=0.65, linestyle = ":")
#axes[0,1].set_xlim(0, 1)
axes[1,1].set_xlabel("$\mathrm{2nd}$ $\mathrm{Planet:}$ i $(\mathrm{degree})$", style='italic') #, size=17)
axes[1,1].set_yticks([])
axes[1,1].yaxis.set_label_coords(-0.06,0.5)
# bmh style add grid as default
axes[1,1].grid(False)
# Remove ticks and spines
axes[1,1].tick_params(left = False, bottom = False)
for axes[1,1], spine in axes[1,1].spines.items():
    spine.set_visible(False)

# var
dat_plt = dat_all2.iloc[:,1]
dat_plt.plot(ax=axes[1,2], kind = "hist", density = True, alpha=0.65, bins = 30) # change density to true, because KDE uses density
dat_plt.plot(ax=axes[1,2], kind = "kde") # Plot KDE
#axes[1,2].set_xlim(0, 2)
#axes[1,2].axvline(vke, c='w', alpha=0.65, linestyle = ":")
axes[1,2].set_xlabel("$\mathrm{2nd}$ $\mathrm{Planet:}$ e", style='italic') #, size=17)
axes[1,2].set_xlim(-0.05, 0.65)
#axes[1,2].set_xlabel(r'$\mathrm{\epsilon_{x}}$ $(\mathrm{\mu as})$') #, size=17)
#axes[0,3].set_xlabel(r'$\omega$') #, size=17)
axes[1,2].set_yticks([])
#axes[1,2].set_ylabel("Marginal Posterior")
axes[1,2].set_ylabel("")
# bmh style add grid as default
axes[1,2].grid(False)
# Remove ticks and spines
axes[1,2].tick_params(left = False, bottom = False)
for axes[1,2], spine in axes[1,2].spines.items():
    spine.set_visible(False)

# period
dat_plt = dat_all2.iloc[:,5]
dat_plt.plot(ax=axes[2,0], kind = "hist", density = True, alpha=0.65, bins = 30) # change density to true, because KDE uses density
dat_plt.plot(ax=axes[2,0], kind = "kde") # Plot KDE
#
#axes[1,3].axvline(ped, c='w', alpha=0.65, linestyle = ":")
#axes.axvline(539.97, alpha = i[1], ymax = i[2], linestyle = ":")
#
#axes[0,1].set_xlim(0, 1)
axes[2,0].set_xlabel(r'$\mathrm{2nd}$ $\mathrm{Planet:}$ $M_\mathrm{p}$ $(M_{\oplus})$') #, size=17)
#set_xlabel(r'$\epsilon_{x}$ $(\mu \mathrm{as})$') #, size=17)
axes[2,0].set_yticks([])
axes[2,0].set_ylabel("Marginal Posterior", labelpad=1.0)
axes[2,0].yaxis.set_label_coords(-0.06,0.5)
# bmh style add grid as default
axes[2,0].grid(False)
# Remove ticks and spines
axes[2,0].tick_params(left = False, bottom = False)
for axes[2,0], spine in axes[2,0].spines.items():
    spine.set_visible(False)

# Mp
dat_plt = dat_all2.iloc[:,7]
dat_plt.plot(ax=axes[2,1], kind = "hist", density = True, alpha=0.65, bins = 30) # change density to true, because KDE uses density
#axes[1,1].axvline(mp, c='w', alpha=0.65, linestyle = ":")
dat_plt.plot(ax=axes[2,1], kind = "kde") # Plot KDE
#axes[0,1].set_xlim(0, 1)
axes[2,1].set_xlabel("$\mathrm{2nd}$ $\mathrm{Planet:}$ P $(\mathrm{days})$", style='italic') #, size=17)
axes[2,1].set_yticks([])
#axes[1,1].set_ylabel("Marginal Posterior")
axes[2,1].set_ylabel("")
# bmh style add grid as default
axes[2,1].grid(False)
# Remove ticks and spines
axes[2,1].tick_params(left = False, bottom = False)
for axes[2,1], spine in axes[2,1].spines.items():
    spine.set_visible(False)

# period
dat_plt = dat_all2.iloc[:,6]
dat_plt.plot(ax=axes[2,2], kind = "hist", density = True, alpha=0.65, bins = 30) # change density to true, because KDE uses density
dat_plt.plot(ax=axes[2,2], kind = "kde") # Plot KDE
#
#axes[1,3].axvline(ped, c='w', alpha=0.65, linestyle = ":")
#axes.axvline(539.97, alpha = i[1], ymax = i[2], linestyle = ":")
#
#axes[0,1].set_xlim(0, 1)
axes[2,2].set_xlabel(r'$\epsilon_{x}$ $(\mu \mathrm{as})$') #, size=17)
#set_xlabel(r'$\epsilon_{x}$ $(\mu \mathrm{as})$') #, size=17)
axes[2,2].set_yticks([])
#axes[1,3].set_ylabel("Marginal Posterior")
axes[2,2].set_ylabel("")
# bmh style add grid as default
axes[2,2].grid(False)
# Remove ticks and spines
axes[2,2].tick_params(left = False, bottom = False)
for axes[2,2], spine in axes[2,2].spines.items():
    spine.set_visible(False)





#plt.show()
plt.savefig(figname)
