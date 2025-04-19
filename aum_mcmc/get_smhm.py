import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from behroz import Behroozi_2013a
bh = Behroozi_2013a()

import sys

chnpath = sys.argv[1]
predpath = sys.argv[2]


try :
    #mh, ymin0, ymax0, ymin1, ymax1 = np.loadtxt('%s_smhm.dat_fsat_leq_0.5'%chnpath, unpack=1)
    #mh, ymed, ymin0, ymax0, ymin1, ymax1 = np.loadtxt('%s_smhm.dat'%chnpath, unpack=1)
    mh, ymed, ymin0, ymax0, ymin1, ymax1 = np.loadtxt('/mnt/home/student/cdivya/github/weaklens_pipeline/preet/sn_preet_gama_output_cosmo_zu15_full_revised_fits/%s_smhm.dat'%chnpath.split('/')[-1], unpack=1)
except:
    pred  = pd.read_csv('%s'%predpath,delim_whitespace=True,header=None)
    chn  = pd.read_csv('%s'%chnpath,delim_whitespace=True,header=None)   
    idx  = (pred.values[:,-1]<300)
    chn = chn.values[idx,:]
    mh = np.linspace(11,16,20)
    mat = np.zeros((len(chn), len(mh)))
    
    for i in range(len(chn[:,0])):
        logM0, logM1, gamma1, gamma2 = chn[i,:4]
        val = logM0 + gamma1*(mh - logM1)-(gamma1 - gamma2)*np.log10(1+10**(mh-logM1))
        mat[i,:] = val
        print(i)
    
    ymed  = np.percentile(mat,50,axis=0)
    ymin0 = np.percentile(mat,16,axis=0)
    ymax0 = np.percentile(mat,84,axis=0)
    ymin1 = np.percentile(mat,2.5,axis=0)
    ymax1 = np.percentile(mat,97.5,axis=0)

    smhmdat = np.transpose([mh, ymed, ymin0, ymax0, ymin1, ymax1])
    
    np.savetxt('/mnt/home/student/cdivya/github/weaklens_pipeline/preet/sn_preet_gama_output_cosmo_zu15_full_revised_fits/%s_smhm.dat'%chnpath.split('/')[-1], smhmdat, header='logmh(h-1Msun)\t50logMcen_m(h-2Msun)\t68logMcen_m(h-2Msun)\t68logMcen_p(h-2Msun)\t95logMcen_m(h-2Msun)\t95logMcen_p(h-2Msun)')
    #np.savetxt('%s_smhm.dat'%chnpath, smhmdat, header='logmh(h-1Msun)\t68logMcen_m(h-2Msun)\t68logMcen_p(h-2Msun)\t95logMcen_m(h-2Msun)\t95logMcen_p(h-2Msun)')

ax =  plt.subplot(2,2,1)
zcl = 0.16
ax.plot(mh, bh.SHMRbeh(mh-np.log10(0.73), zcl)+2*np.log10(0.73), ls='--', color='black', label='behroozi-2013a', zorder=2)

ax.plot(mh, ymed, color='C0', ls='--') # 95 percentile
ax.fill_between(mh, ymin1, ymax1, color='#52aae7', alpha = 0.45, zorder=0) # 95 percentile
ax.fill_between(mh, ymin0, ymax0, color='#1f77b4', alpha = 0.45, zorder=1) # 68 percentile

#ax.fill_between(mh, ymin1, ymax1, color='#52aae7', alpha = 0.45, zorder=0) # 95 percentile
#ax.fill_between(mh, ymin0, ymax0, color='#1f77b4', alpha = 0.45, zorder=1) # 68 percentile

yang = np.loadtxt('smhm_yang_2008_fig8.csv')
x = yang[:,0]
y = yang[:,2]
yerr = [yang[:,2]-yang[:,1], yang[:,3]-yang[:,2]]
ax.errorbar(x, y, yerr=yerr, fmt='.', color='C1', capsize=3, zorder=1, label='yang2008')

zu = np.loadtxt('smhm_zu_2015_fig8.dat')
x = zu[:,0]
y = zu[:,1]
ax.plot(x[:13], y[:13], ':', c='C2', label='zu2015', zorder=3)
ax.plot(x[13:], y[13:], ':', c='C2', zorder=3)


zu = np.loadtxt('smhm_van_uitert_2016_fig3.dat')
x = zu[:,0]
y = zu[:,1]
ax.plot(x[:8], y[:8], '-.', ms=1.0, c='C3', label='van2016', zorder=4)
ax.plot(x[8:], y[8:], '-.', ms=1.0, c='C3', zorder=4)


zu = np.loadtxt('smhm_1d_dvornik_2020_fig8_.dat')
x = zu[:,0] + 2*np.log10(0.73)
y = zu[:,1] + np.log10(0.73)
ax.plot(y, x, '--', color='C4', label='dvornik2020-1d', zorder=5)

#zu = np.loadtxt('smhm_2d_dvornik_2020_fig8_.dat')
#x = zu[:,0] + 2*np.log10(0.73)
#y = zu[:,1] + np.log10(0.73)
#ax.plot(y, x, '--', color='C5', label='dvornik2020-2d', zorder=6)



ax.legend(loc='lower right')
ax.set_xlabel(r'$\log[{\rm M_{\rm h} /(h^{-1} M_\odot)}]$')
ax.set_ylabel(r'$\log[{\rm M_{\rm *,c} /(h^{-2} M_\odot)}]$')

ax.set_ylim(6.5,13)
plt.savefig('/mnt/home/student/cdivya/github/weaklens_pipeline/preet/sn_preet_gama_output_cosmo_zu15_full_revised_fits/%s_smhm.png'%chnpath.split('/')[-1], dpi=250)    
#plt.savefig('%s_smhm_fsat_leq_0.5.png'%chnpath,dpi=250)    

