import numpy as np
import matplotlib.pyplot as plt

def ab2loglum(x):
    return -(x-4.67)/2.5

logM0  = 9.95
logM1  = 11.24
gamma1 = 3.18
gamma2 = 0.245
sig0   = 0.157
b0     = -1.17
b1     = 1.53
b2     = -0.217
alphas = -1.18


data = np.loadtxt('clf_esd_bin5_marc2013.csv')
rbin = 10**data[1::3,0]
esds = 10**data[1::3,1]
yerrl = esds - 10**data[0::3,1]
yerru = 10**data[2::3,1] - esds

plt.subplot(2,2,1)

m1 = -22.0
m2 = -21.5
rsft  = 0.17

plt.errorbar(rbin, esds, yerr=[yerrl,yerru], fmt='.',capsize=3, label='(%2.2f,%2.2f)'%(m1,m2))

logMb = ab2loglum(m1) 
logMa = ab2loglum(m2)

from sn_mcmc_csmf import *
print logMa,logMb
 #     logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alpsat, cfac, poff, roff     , ap = x

aum = init_aum()
x = [logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alphas, 1.0, 0.0,0.1,0.0]

plt.plot(rbin, esd(x, rbin, rsft, logMa, logMb, 10.0, aum=aum)[0], c='C1')

plt.yscale('log')
plt.xscale('log')

plt.ylabel(r'$\Delta \Sigma [{\rm h  M_\odot pc^{-2}}]$')
plt.xlabel(r'$R[{\rm h^{-1} Mpc}]$')
plt.legend()
plt.savefig('chk_marc.png',dpi=300)



