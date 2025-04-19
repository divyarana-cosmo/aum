import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from glob import glob

lenspath = sys.argv[1]
randpath = sys.argv[2]

flist = glob('%s/dsigma.dat_*'%lenspath)

nbins = len(flist)

for bb,ii in enumerate(range(1,nbins+1)):
    n = int((bb)/2)
    
    filename = '%s/dsigma.dat_%d' % (lenspath,ii)
    dfdict = pd.read_csv(filename, delim_whitespace = True, usecols=([5,6,7,8,11,13,14]), header=None, names=(["dsigma","Sumwls_by_sumwl","r","dsigma_cross","sn_err","tot_pairs","r12_sel"]), comment='#')
    xjk = dfdict.r.values
    y1jk = dfdict.Sumwls_by_sumwl.values

    filename = '%s/dsigma.dat_%d' % (randpath,ii)
    dfdict = pd.read_csv(filename, delim_whitespace = True, usecols=([5,6,7,8,11,13,14]), header=None, names=(["dsigma","Sumwls_by_sumwl","r","dsigma_cross","sn_err","tot_pairs","r12_sel"]), comment='#')
    rrxjk = dfdict.r.values
    rry1jk = dfdict.Sumwls_by_sumwl.values

    bst = y1jk*1.0/rry1jk  #jackknife regions vise boost param
<<<<<<< HEAD
    print(ii, np.shape(y1jk), np.shape(rry1jk))
=======
    print ii, np.shape(y1jk), np.shape(rry1jk)
>>>>>>> 2cb08947c4825a1e0e38de69838ff8bfaff35728
    x = np.unique(xjk)
    njacks = int(len(xjk)/len(x))
    y1 = 0.0*x
    y1err = 0.0*x
    
    for jj,rr in enumerate(x):
        idx = (xjk==rr)
        val = bst[idx]
        y1[jj] = np.mean(val)
        y1err[jj] = np.sqrt(njacks -1) * np.std(val)

    #np.savetxt('./jk_output_1/boost_1.dat_%d'%ii,np.transpose([x,y1,y1err]),header='rbins\tmean_boost\tstd_boost')   
    ax = plt.subplot(3,3,n+1)
    ax.errorbar(x, y1, yerr=y1err, fmt='.', label='bin=%d'%(bb+1), capsize=3)

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
    #    ax.set_ylim(0.98,1.08)
        
    ax.axhline(y=1.0,ls='--',color='grey')
    ax.legend()

plt.tight_layout()
plt.savefig('./%s/boost.png'%lenspath,dpi=300)

