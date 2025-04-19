import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
for ii in range(1,13):
    n = int((ii-1)/2)

    dat = np.loadtxt('./jk_output_1/avg_dsigma.dat_%d'%ii)
    rdat = np.loadtxt('./jk_random_1/avg_dsigma.dat_%d'%ii)

    ax = plt.subplot(3,3,n+1)

    frac =  dat[:,-2]*(1.0/rdat[:,-2])
    fracerr = frac * np.sqrt((dat[:,-1]*1.0/dat[:,-2])**2 + (rdat[:,-1]*1.0/rdat[:,-2])**2)

    ax.errorbar(dat[:,0], frac, yerr=fracerr, fmt='.', label='bin=%d'%ii, capsize=3)
    if n>2:
        ax.set_xlabel(r'$R[{\rm h^{-1}Mpc}]$')
    
    ax.set_xscale('log')
    if n==1:
        ax.legend(markerscale=1.0,loc='upper right')
    
    if n==0 or n==3:
        ax.set_ylabel(r'$C(R)$')
                    
    if n<=2:
        ax.set_xticklabels([])

    ax.axhline(y=1.0,ls='--',color='grey')
    ax.legend()

plt.tight_layout()
plt.savefig('./jk_output_1/boost.png',dpi=300)


