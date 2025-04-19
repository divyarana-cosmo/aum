import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=[10,10])
#for bb,ii in enumerate(range(3,13)):
for ii,bb in enumerate(range(0,8)):

    dat0  = np.loadtxt('./sn_preet_gama_output_cosmo_zu15_full_revised/dsigma.dat_%d'%bb);
    rdat0 = np.loadtxt('./sn_preet_gama_rand_output_cosmo_zu15_full_revised/avg_dsig_uni_%d.dat'%bb);

    ax = plt.subplot(4,4,ii+1)

    y0 = dat0[:10,5]*(1.0/(1+dat0[:10,14])) - rdat0[:10,1]
    #y0 = dat0[:,5]*(1.0/(1+dat0[:,14])) - rdat0[:,1]
    uqr = np.unique(dat0[:10,7])
    avgy0 = y0
 
    
    cov0 = np.loadtxt('./sn_preet_gama_output_cosmo_zu15_full_revised/cov_no_boost.dat',dtype=float,unpack=True)
    x0 = len(uqr)*(ii)
    x1 = len(uqr)*(ii+1)

    cov0 = cov0[x0:x1,x0:x1]
    yerr0 = np.sqrt(np.diag(cov0))
    icov0 = np.linalg.inv(cov0)

    chisq0 = np.dot(avgy0, np.dot(icov0, avgy0))
 
    ax.errorbar(uqr, avgy0, yerr0, fmt='.', label=r'${\rm SNR} = %2.2f$'%(np.sqrt(chisq0)), capsize=3)
    
    ax.set_ylim(0.1,450) 
    ax.set_yscale('log')
    ax.set_xscale('log')
    if ii==0 or ii==4 :#or ii==6  :
       ax.set_ylabel(r'$\Delta\Sigma [{\rm h M_\odot pc^{-2}}]$')

    if not (ii==0 or ii==4): #or ii==6) :
        ax.set_yticklabels([])
    if not (ii==4 or ii==5 or ii==6 or ii==7):
        ax.set_xticklabels([])
    if ii==4 or ii==5 or ii==6 or ii==7:
        ax.set_xlabel(r'$R[{\rm h^{-1}Mpc}]$')

       

    #ax.set_ylabel(r'$\Delta\Sigma [{\rm h M_\odot pc^{-2}}]$')
    
    ax.legend()
   
     
    ax.set_xlabel(r'$R[{\rm h^{-1}Mpc}]$')

#plt.tight_layout()
plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.savefig('./sn_preet_gama_output_cosmo_zu15_full_revised/esds_no_boost_uni.png',dpi=300)
#plt.savefig('./sn_preet_gama_output_cosmo_y08_full_revised/esds_no_boost.png',dpi=300)
#plt.savefig('./re_full_jk_output_1_cosmo_y08_deblend/esds_no_boost.png',dpi=300)
