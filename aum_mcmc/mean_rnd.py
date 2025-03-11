import matplotlib.pyplot as plt
import numpy as np
import glob
import sys
fpath = 'sn_preet_gama_rand_output_cosmo_zu15_full_revised' #sys.argv[1]
#fpath = 'sn_preet_gama_rand_output_cosmo_y08_full_revised' #sys.argv[1]
#nbin  = int(sys.argv[2])


from astropy.io import fits
import healpy as hp
full_gal = fits.open('../DataStore/preetish/revised_gama_specz_sample_zleq_03.fits')[1].data
msk = hp.read_map("/mnt/home/student/cdivya/github/weaklens_pipeline/DataStore/data/S16A_mask_w_holes.fits")

def avg_rand(nbin):
    njacks = 100


    ipix = hp.ang2pix(int(np.sqrt(msk.size/12)), full_gal['ra'], full_gal['dec'], lonlat=1)
    flg  = msk[ipix]
    idx  = (full_gal['bin']==nbin) & (flg==1.0) & (full_gal['logmstar_h2'] < 12.0)

    logMa    = float(full_gal['logMstar_h2'][idx].min())
    logMb    = float(full_gal['logMstar_h2'][idx].max())
    logMstel = float(np.log10(np.mean(10**full_gal['logMstar_h2'][idx])))
    zred     = np.median(full_gal['z_tonry'][idx])



    dat             = np.loadtxt('%s/dsigma_uni.dat%05d_%d'%(fpath,0,0))
    #dat             = np.loadtxt('%s/dsigma.dat%05d_%d'%(fpath,0,0))
    rbin            = 0.0*dat[:,7]
    dsig            = 0.0*dat[:,5]/(1.0 + dat[:,14])
    dsigerr         = 0.0*(dat[:,5]/(1.0 + dat[:,14]))**2
    avg_sumwls      = 0.0*dat[:,6]
    avg_sumwls_err  = 0.0*dat[:,6]**2
    
    
    plt.subplot(2,2,1)
    
    for cnt in range(njacks):
        if cnt==0:
            dat     = np.loadtxt('%s/dsigma_uni.dat%05d_%d'%(fpath,cnt,nbin))
            #dat     = np.loadtxt('%s/dsigma.dat%05d_%d'%(fpath,cnt,nbin))
            rbin    = dat[:,7]
            dsig    = dat[:,5]/(1.0 + dat[:,14])
            dsigerr = (dat[:,5]/(1.0 + dat[:,14]))**2
            avg_sumwls     = dat[:,6]
            avg_sumwls_err = dat[:,6]**2
        else:
            dat     = np.loadtxt('%s/dsigma_uni.dat%05d_%d'%(fpath,cnt,nbin))
            #dat     = np.loadtxt('%s/dsigma.dat%05d_%d'%(fpath,cnt,nbin))
            dsig    = dsig + dat[:,5]/(1.0 + dat[:,14])
            dsigerr = dsigerr + (dat[:,5]/(1.0 + dat[:,14]))**2
            avg_sumwls     = avg_sumwls + dat[:,6]
            avg_sumwls_err = avg_sumwls_err + dat[:,6]**2
    
    rbin = rbin
    dsig        = dsig*1.0/njacks
    dsigerr     = np.sqrt(dsigerr*1.0/njacks - (dsig)**2)
    avg_sumwls  = avg_sumwls*1.0/njacks
    avg_sumwls_err = np.sqrt(avg_sumwls_err*1.0/njacks - (avg_sumwls)**2)
    
    dat = np.transpose([rbin, dsig, dsigerr, avg_sumwls, avg_sumwls_err])
    np.savetxt('./%s/avg_dsig_uni_%d.dat'%(fpath, nbin), dat, header='0-rbin\t1-dsigma\t2-dsigmaerrt\t3-avgwls\t4-avgwlserr')
    #np.savetxt('./%s/avg_dsig_%d.dat'%(fpath, nbin), dat, header='0-rbin\t1-dsigma\t2-dsigmaerrt\t3-avgwls\t4-avgwlserr')
    
    plt.axhline(y=0.0, ls='--', color='grey', zorder=4)
    plt.errorbar(rbin[:10], dsig[:10]*rbin[:10], yerr=dsigerr[:10]*rbin[:10]/10.0, fmt='.', capsize=3, zorder=10, label=r"$[%2.2f,%2.2f)$"%(logMa,logMb))
    
    plt.ylim(-0.5,0.5)
    plt.xscale('log')
    plt.xlabel(r'$R [{\rm h^{-1}Mpc}]$')
    plt.ylabel(r'$R \times \Delta \Sigma [{\rm 10^6 M_\odot pc^{-1}}]$')
    plt.legend()
    plt.savefig('./%s/rand_uni_%d.pdf'%(fpath,nbin), dpi=300)
    plt.savefig('./%s/rand_uni_%d.png'%(fpath,nbin), dpi=300)
    plt.clf() 
    return 0



for nn in range(0,8):
    avg_rand(nn)
