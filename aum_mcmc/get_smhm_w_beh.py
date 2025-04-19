import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from behroz import Behroozi_2013a
bh = Behroozi_2013a()

import sys


chnpath = sys.argv[1]
predpath = sys.argv[2]

chn  = pd.read_csv(chnpath, delim_whitespace=True,header=None)   
pred = pd.read_csv(predpath, delim_whitespace=True,header=None)   

idx = (pred.values[:,-1]<200)
minchi = np.min(pred.values[:,-1])


chn = chn.values[idx,:]

mh = np.linspace(11,16,20)

mat = np.zeros((len(chn), len(mh)))
#mat = np.zeros((len(chn.values[:,0]), len(mh)))

#for i in range(len(chn.values[:,0])):
def ffunc(p,x):
    logeps, logM1, alpha, delta, gamma = p
    val = -np.log10(10**(alpha * x) + 1) + delta*(np.log10(1+np.exp(x))**gamma)/(1 + np.exp(10**(-x)))
    return val

def beh_2013(x,xm):
    logeps, logM1, alpha, delta, gamma = x
    xm = xm - np.log10(0.73)
    val = logeps + logM1 + ffunc(x,xm-logM1) - ffunc(x,0) + 2*np.log10(0.73)
    return val


for i in range(len(chn[:,0])):
    p = chn[i,:5]

    val = beh_2013(p,mh)
    mat[i,:] = val
    print(minchi,i)


ax =  plt.subplot(2,2,1)
zcl = 0.158
ax.plot(mh, bh.SHMRbeh(mh-np.log10(0.73), zcl)+2*np.log10(0.73), ls='--', color='black', label='behroozi-2013a')

ymin0 = np.percentile(mat,16,axis=0)
ymax0 = np.percentile(mat,84,axis=0)
 
ymin1 = np.percentile(mat,2.5,axis=0)
ymax1 = np.percentile(mat,97.5,axis=0)
 
ax.fill_between(mh, ymin1, ymax1, color='#52aae7', alpha = 0.45, zorder=1) # 95 percentile
ax.fill_between(mh, ymin0, ymax0, color='#1f77b4', alpha = 0.45, zorder=2) # 68 percentile

#nfw = pd.read_csv('chainfile.dat_nfw',delim_whitespace=True,header=None)
#logmh = nfw.values[:,:10]
#mstels = np.loadtxt('../../DataStore/preetish/mstel_bins.dat')
#mstels = mstels[2:-1,:]
#x = np.percentile(logmh,50, axis=0)
#xerr = [x - np.percentile(logmh,16,axis=0), np.percentile(logmh,84,axis=0)- x]
#
#ax.errorbar(x, mstels[:,3], xerr=xerr, fmt='.', color='C2', capsize=3, zorder=10, label='1h-nfw-fits')


yang = np.loadtxt('smhm_yang_2008_fig8.csv')
x = yang[:,0]
y = yang[:,2]
yerr = [yang[:,2]-yang[:,1], yang[:,3]-yang[:,2]]

ax.errorbar(x, y, yerr=yerr, fmt='.', color='C1', capsize=3, zorder=10, label='yang2008')

ax.legend(loc='lower right')
ax.set_xlabel(r'$\log[{\rm M_{\rm h} /(h^{-1} M_\odot)}]$')
ax.set_ylabel(r'$\log[{\rm M_{\rm *,c} /(h^{-2} M_\odot)}]$')

ax.set_ylim(6.5,13)
plt.savefig('%s_smhm.png'%chnpath ,dpi=250)    

