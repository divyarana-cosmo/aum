import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import healpy as hp

full_gal = fits.open('../DataStore/preetish/revised_gama_specz_sample_zleq_03.fits')[1].data
msk = hp.read_map("/mnt/home/student/cdivya/github/weaklens_pipeline/DataStore/data/S16A_mask_w_holes.fits")


plt.figure(figsize=[8,6])
#for bb,ii in enumerate(range(3,13)):
for bb,ii in enumerate(range(7)):
    ipix = hp.ang2pix(int(np.sqrt(msk.size/12)), full_gal['ra'], full_gal['dec'], lonlat=1)
    flg  = msk[ipix]
    idx  = (full_gal['bin']==bb) & (flg==1.0) & (full_gal['logmstar_h2'] < 12.0)

    logMa    = float(full_gal['logMstar_h2'][idx].min())
    logMb    = float(full_gal['logMstar_h2'][idx].max())
    logMstel = float(np.log10(np.mean(10**full_gal['logMstar_h2'][idx])))
    zred     = np.median(full_gal['z_tonry'][idx])
    
<<<<<<< HEAD
    print(('%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\n'%(logMa, logMb, logMstel, zred, sum(idx))))
=======
    print('%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\n'%(logMa, logMb, logMstel, zred, sum(idx)))
>>>>>>> 2cb08947c4825a1e0e38de69838ff8bfaff35728
    dat0  = np.loadtxt('./sn_preet_gama_output_cosmo_zu15_full_revised/dsigma.dat_%d'%ii);

    ax = plt.subplot(3,4,bb+1)
    y = dat0[:10,8]*(1.0/(1+dat0[:10,14]))
    uqr = np.unique(dat0[:10,7])
    yerr = 0.0*uqr
 
    cov = np.loadtxt('./sn_preet_gama_output_cosmo_zu15_full_revised/cov_no_boost_x.dat',dtype=float,unpack=True)
    x0 = 10*(bb)
    x1 = 10*(bb+1)

    yerr = np.sqrt(np.diag(cov))[x0:x1]
    icov = np.linalg.inv(cov[x0:x1,x0:x1])

    chisq0 = np.dot(y, np.dot(icov, y))
    from scipy.stats import chi2
    
    ax.errorbar(uqr, y*uqr, yerr*uqr, fmt='.', label=r"$[%2.2f,%2.2f)$"%(logMa,logMb), capsize=3)
    ax.text(0.1, 0.05, r'${\rm p value} = %2.2f$'%(chi2.sf(chisq0,df=len(uqr))), transform=ax.transAxes, fontsize=6, zorder=10)
    
    #ax.set_title(r'$bin=%d$'%(bb))
    #ax.set_ylim(0.05,4500) 
    #ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel(r'$R\times\Delta\Sigma_{\times}[{\rm 10^6 M_\odot pc^{-1}}]$')
    ax.legend(markerscale=1.0, ncol=1, fontsize='xx-small', loc='upper right', frameon=0)
    ax.set_ylim(-5,5) 
    ax.axhline(0.0, ls='--', color='grey') 
    ax.set_xlabel(r'$R[{\rm h^{-1}Mpc}]$')

plt.tight_layout()
    
    
#plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.savefig('./sn_preet_gama_output_cosmo_zu15_full_revised/esds_x.png',dpi=300)
#plt.savefig('./re_full_jk_output_1_cosmo_y08_deblend/esds_no_boost.png',dpi=300)
