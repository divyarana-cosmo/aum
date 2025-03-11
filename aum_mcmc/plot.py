import numpy as np
import matplotlib.pyplot as plt

mstel_bin = np.loadtxt('../DataStore/preetish/mstel_bins_gama.dat')

plt.figure(figsize=[8,6])
#for bb,ii in enumerate(range(3,13)):
for bb,ii in enumerate(range(1,12)):
    sel = (mstel_bin[:,0]==ii)
    logMa = float(mstel_bin[sel,1])
    logMb = float(mstel_bin[sel,2])

    dat0  = np.loadtxt('./jk_preet_gama_output_cosmo_y08/dsigma.dat_%d'%ii);
    rdat0 = np.loadtxt('./jk_preet_gama_rand_output_cosmo_y08/dsigma.dat_%d'%ii);

    ax = plt.subplot(3,4,bb+1)

    y0 = dat0[:,5]*(1.0/(1+dat0[:,14])) - rdat0[:,5]*(1.0/(1+rdat0[:,14]))
    uqr = np.unique(dat0[:,7])
    avgy0 = 0.0*uqr
    yerr0 = 0.0*uqr
 
    njacks = int(len(y0)/len(uqr))
    for jj,rr in enumerate(uqr):
        idx = (dat0[:,7]==rr)
        avgy0[jj] = np.mean(y0[idx])
        yerr0[jj] = np.sqrt(njacks - 1) * np.std(y0[idx])  
    
    cov0 = np.loadtxt('./jk_preet_gama_output_cosmo_y08/cov_w_rand_subs_no_boost.dat',dtype=float,unpack=True)
    x0 = 10*(bb)
    x1 = 10*(bb+1)
    cov0 = cov0[x0:x1,x0:x1]
    icov0 = np.linalg.inv(cov0)

    chisq0 = np.dot(avgy0, np.dot(icov0, avgy0))
 
    ax.errorbar(uqr, avgy0, yerr0, fmt='.', label=r'${\rm SNR} = %2.2f$'%(np.sqrt(chisq0)), capsize=3)
    
    ax.set_title(r'$\log M_* \in (%2.2f,%2.2f]$'%(logMa,logMb))
    ax.set_ylim(0.05,4500) 
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel(r'$\Delta\Sigma [{\rm h M_\odot pc^{-2}}]$')
    
    ax.legend()
   
     
    ax.set_xlabel(r'$R[{\rm h^{-1}Mpc}]$')

plt.tight_layout()
    
    
#plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.savefig('./jk_preet_gama_output_cosmo_y08/esds_no_boost.png',dpi=300)
#plt.savefig('./re_full_jk_output_1_cosmo_y08_deblend/esds_no_boost.png',dpi=300)
