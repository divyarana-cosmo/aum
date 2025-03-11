import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from glob import glob

lenspath = 'sn_preet_gama_output_cosmo_zu15_full_revised'
randpath = 'sn_preet_gama_rand_output_cosmo_zu15_full_revised'
nbins = 8

for ii,bb in enumerate(range(0,nbins)):
    n = int((bb)/2)
    
    filename = '%s/dsigma.dat_%d' % (lenspath,bb)
    dfdict = pd.read_csv(filename, delim_whitespace = True, usecols=([5,6,7,8,11,13,14]), header=None, names=(["dsigma","Sumwls_by_sumwl","r","dsigma_cross","sn_err","tot_pairs","r12_sel"]), comment='#')
    xjk = dfdict.r.values
    y1jk = dfdict.Sumwls_by_sumwl.values

    filename = '%s/avg_dsig_uni_%d.dat' % (randpath,bb)
    dfdict = np.loadtxt(filename)

    rrxjk = dfdict[:,0]
    rry1jk = dfdict[:,3]
    rry1errjk = dfdict[:,4]
    bst = y1jk*1.0/rry1jk  
    
    x  = xjk
    y1 = bst
    y1err = bst * (rry1errjk * 1.0/rry1jk)
    np.savetxt('./%s/boost.dat_uni_%d'%(lenspath,bb),np.transpose([x,y1,y1err]),header='rbins\tmean_boost\tstd_boost')   
    ax = plt.subplot(3,3,n+1)
    ax.errorbar(x[:10], y1[:10], yerr=y1err[:10], fmt='.', label='bin=%d'%(bb), capsize=3)
    #ax.errorbar(x, y1, yerr=y1err, fmt='.', label='bin=%d'%(bb), capsize=3)

    if n>2:
        ax.set_xlabel(r'$R[{\rm h^{-1}Mpc}]$')
    ax.set_xscale('log')
    if n==1:
        ax.legend(markerscale=1.0,loc='upper right')
    if n==0 or n==3:
        ax.set_ylabel(r'$C(R)$')
    if n<=2:
        ax.set_xticklabels([])
    #if n!=6:
    
    ax.set_ylim(0.96,1.02)
        
    ax.axhline(y=1.0,ls='--',color='grey')
    ax.legend()

plt.tight_layout()
plt.savefig('./%s/boost_uni.png'%lenspath,dpi=300)

