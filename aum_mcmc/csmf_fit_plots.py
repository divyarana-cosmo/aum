import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import corner
from astropy.io import fits
import sys

plt.figure(figsize=[8.5,8.5])
#plt.figure(figsize=[12.5,10])

from math import log10
import numpy as np
from scipy.integrate import quad
#from sn_mcmc_csmf import *
from sn_mcmc_csmf_w_loglin_alpsat_fid import *

#chnname  = 'chainfile_bosch_2013_loglin_alpsat_uni_five_per_boost.dat'
#predname =  'predfile_bosch_2013_loglin_alpsat_uni_five_per_boost.dat'

#chnname  = 'chainfile_bosch_2013_loglin_alpsat_uni.dat'
#predname =  'predfile_bosch_2013_loglin_alpsat_uni.dat'

chnname  = 'chainfile_bosch_2013_loglin_alpsat_uni_const_alpha_no_off.dat'
predname = 'predfile_bosch_2013_loglin_alpsat_uni_const_alpha_no_off.dat'


chn  = pd.read_csv('/scratch/cdivya/%s'%chnname, delim_whitespace=True, header=None)
pred = pd.read_csv('/scratch/cdivya/%s'%predname, delim_whitespace=True, header=None) 

#chn  = pd.read_csv('./sn_preet_gama_output_cosmo_zu15_full_revised_fits/%s'%chnname, delim_whitespace=True, header=None)
#pred = pd.read_csv('./sn_preet_gama_output_cosmo_zu15_full_revised_fits/%s'%predname, delim_whitespace=True, header=None) 


idx = (pred.values[:,-1] < 300)
chn = chn.values[idx,:]
pred = pred.values[idx,:]

cov = np.loadtxt('./sn_preet_gama_output_cosmo_zu15_full_revised/cov_no_boost.dat',dtype=float,unpack=True)
#cov = cov *1.05**2 #boosted
from astropy.io import fits
import healpy as hp
full_gal = fits.open('../DataStore/preetish/revised_gama_specz_sample_zleq_03.fits')[1].data
msk = hp.read_map("/mnt/home/student/cdivya/github/weaklens_pipeline/DataStore/data/S16A_mask_w_holes.fits")


def flat_ele(x0, x1):
    ilim = np.transpose([x0,x1])
    cij = np.zeros(len(x0))
    cnt=0
    for n in range(len(x0)):
        i = ilim[n][0]
        j = ilim[n][1]
        a = ((j**3 - i**3)/3.0)/(j - i)
        b = ((j**2 - i**2)/2.0)/(j  - i)
        b = b**2
        cij[n] = a-b
    return cij


def n_eff(mat):
    #logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alpsat, eta, cfac, poff, roff , ap = x
    chnn       = 0.0*mat[:,:11]
    chnn[:,:9] = mat[:,:9]
    chnn[:,9]  = mat[:,10]
    chnn[:,10] = mat[:,13]

    #logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alpsat, cfac, ap = x
    p = len(chnn[0,:])
    cij = np.zeros(p)
    #flat prior limits
    x0 = [8.0, 9.5, 1.0, 0.01, 0.05, -2, -2, -2, -2, 0.5]
    x1 = [ 11.5, 15, 4.0, 2.0, 0.5, 2, 2, 2, -1, 5.0]

    ss = flat_ele(x0, x1)
    cij[:-2] = ss[:-1]
    cij[-2] = 0.2
    cij[-1] = ss[-1]

    ic_prior = np.diag(1.0/cij)

    cov = 0.0 * ic_prior
    size = len(cov[:,0])

    for i1 in range(size):
        for i2 in range(size):
            r0 = chnn[:,i1] - np.mean(chnn[:,i1])
            r1 = chnn[:,i2] - np.mean(chnn[:,i2])
            cov[i1][i2] = np.mean(r0*r1)

    mat = np.dot(ic_prior, cov)
    neff = p - np.trace(mat)
    print((p,neff))
    return neff




#logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alpsat, cfac, ap = x

lab = [r'${\rm \log M_{0}}$', r'${\rm \log M_1}$', r'${\rm \gamma_1}$', r'${\rm \gamma_2}$', r'${\rm \sigma_0}$', r'${\rm b_0$}', r'${\rm b_1}$', r'${\rm b_2}$',r'${\rm \alpha_s}$',r'${\rm c_{fac}}$',r'${\rm a_p}$']

def texit(chn=chn):
    dat       = 0.0*chn[:,:11]
    dat[:,:9] = chn[:,:9]
    dat[:,9]  = chn[:,10]
    dat[:,10] = chn[:,13]

    for ii in range(len(dat[0,:])):
        med = np.median(dat[:,ii])
        perr = np.percentile(dat[:,ii],84) - med
        merr = med - np.percentile(dat[:,ii],16)
        print(('\t %s & $%2.2f^{+%2.2f}_{-%2.2f} \\vspace{0.1cm} $\\\\'%(lab[ii], med, perr, merr)))
    return 0




def pltsigdata(ii,bb):
    ipix = hp.ang2pix(int(np.sqrt(msk.size/12)), full_gal['ra'], full_gal['dec'], lonlat=1)
    flg  = msk[ipix]
    idx  = (full_gal['bin']==bb) & (flg==1.0) & (full_gal['logmstar_h2'] < 12.0)

    logMa    = float(full_gal['logMstar_h2'][idx].min())
    logMb    = float(full_gal['logMstar_h2'][idx].max())
    logMstel = float(np.log10(np.mean(10**full_gal['logMstar_h2'][idx])))
    zred     = np.median(full_gal['z_tonry'][idx])

    dat  = np.loadtxt('./sn_preet_gama_output_cosmo_zu15_full_revised/dsigma.dat_%d'%bb);
    rdat = np.loadtxt('./sn_preet_gama_rand_output_cosmo_zu15_full_revised/avg_dsig_uni_%d.dat'%bb);

    y = dat[:,5]*(1.0/(1+dat[:,14])) - rdat[:,1]
    #y = dat[:,5]*(1.0/(1+dat[:,14]))*1.05 - rdat[:,1]
    uqr = np.unique(dat[:10,7])
    avgy = y[:10] 

    x0 = ii*len(uqr)
    x1 = (ii + 1)*len(uqr)
    yerr = np.sqrt(np.diag(cov)[x0:x1])
    icov = np.linalg.inv(cov[x0:x1,x0:x1])
    njacks = 500#int(len(y)/len(uqr)
    #print njacks
    #for jj,rr in enumerate(uqr):
    #    idx = (dat[:,7]==rr)
    #    avgy[jj] = np.mean(y[idx])
    #    yerr[jj] = np.sqrt(njacks - 1) * np.std(y[idx])  

    #x0 = len(uqr)*ii
    #x1 = len(uqr)*(ii+1)
    #if ii<5:
    #    ii=ii
    #else:
    #    ii=ii+5

    ax = plt.subplot(4,4,ii+1)
    mat = np.transpose([uqr, avgy, yerr])
    np.savetxt('/mnt/home/student/cdivya/github/size_mass_gama_hsc_s16a/dsigma_bin_%d.dat'%(ii+1), mat)
    np.savetxt('/mnt/home/student/cdivya/github/size_mass_gama_hsc_s16a/cov.dat', cov)

    #ax = plt.subplot(4,5,ii+1)
    ax.errorbar(uqr, avgy, yerr, fmt='.', label=r"$[%2.2f,%2.2f)$"%(logMa,logMb),capsize=3, zorder=10)
    #ax.errorbar(uqr, avgy, yerr, fmt='.', label='[%2.2f,%2.2f)'%(logMa,logMb), capsize=3, zorder=10)
    y    = pred[:,x0:x1]
    #icov = np.linalg.inv(cov[x0:x1,x0:x1])
    hartlap = (njacks-len(uqr)-2) /(njacks - 1.0)
    #sys.exit()
    icov = hartlap*icov
    chisq = pred[:,-1]
    chain = chn[:,:]
    ymin0 = np.percentile(y,16,axis=0)
    ymax0 = np.percentile(y,84,axis=0)
    
    ymin1 = np.percentile(y,2.5,axis=0)
    ymax1 = np.percentile(y,97.5,axis=0)
    
    print(len(uqr), len(ymin1), x0, x1, len(pred[0,:]))
    ax.fill_between(uqr, ymin1, ymax1, color='#52aae7', alpha = 0.45, zorder=1) # 95 percentile
    ax.fill_between(uqr, ymin0, ymax0, color='#1f77b4', alpha = 0.45, zorder=2) # 68 percentile
    
    #best fit model
    idx = np.argmin(chisq)
    del_sig = y[idx,:]
    #bestp = chain[idx,:14]
    #print bestp, len(bestp)
    bestp = chain[idx,:14]

    Delta = y[idx,:] - avgy
    res = Delta*1.0/yerr

    dof = 70 - n_eff(chain) # we have 70 datapoints

    #chisq = np.dot(Delta, np.dot(icov, Delta))
    snr   = (np.dot(avgy, np.dot(icov, avgy)))**0.5
    if ii==0:
        ax.plot(uqr, y[idx,:], '-', color='red', label=r'$\chi^2/{\rm dof} = %2.2f/%2.2f=%2.2f$'%(np.min(pred[:,-1]),dof,np.min(pred[:,-1])/dof))
        #ax.legend(markerscale=1.0, ncol=1, fontsize='xx-small', loc='upper right', frameon=0)
    else:    
        ax.plot(uqr, y[idx,:], '-', color='red')
    #ax.plot(uqr, y[idx,:], '-', color='red', label= r"$\chi^2 = %2.2f$"%(np.dot(Delta, np.dot(icov, Delta))))
    ax.legend(markerscale=1.0, ncol=1, fontsize='xx-small', loc='upper right', frameon=0)
    #ax.text(0.5, 0.05, r"$[%2.2f,%2.2f)$"%(logMa,logMb), transform=ax.transAxes, fontsize=6)
    ax.text(0.1, 0.05, r"${\rm SNR=%2.2f}$"%(snr), transform=ax.transAxes, fontsize=6, zorder=10)
    aum = init_aum()

    full = esd(bestp, uqr, zred, logMa=logMa, logMb=logMb, logMstel=logMstel, aum = aum, whichopt=0)
    ohc  = esd(bestp, uqr, zred, logMa=logMa, logMb=logMb, logMstel=logMstel, aum = aum, whichopt=1)
    ohs  = esd(bestp, uqr, zred, logMa=logMa, logMb=logMb, logMstel=logMstel, aum = aum, whichopt=2)
    twh  = full[0] - ohc - ohs - stellar(bestp[-1], logMstel, uqr)
    esdstel = stellar(bestp[-1], logMstel, uqr)
    ax.set_title(r'$\log\langle {\rm M_c} \rangle = %2.2f, f_{\rm sat} = %2.4f$'%(full[1],full[2]))
    ax.plot(uqr, ohc, ls='--', label = "1h-c")
    ax.plot(uqr, ohs, ls='-.', label = "1h-s")
    ax.plot(uqr, twh, ls=':', label = "2h")
    ax.plot(uqr, esdstel, ls=':', label = "bary")
    
    if ii<=3:
        ax.set_ylim(1,100)
    else:
        ax.set_ylim(4,400)

    ax.set_yscale('log')
    ax.set_xscale('log')
    #ax.legend(markerscale=1.0, ncol=1, fontsize='xx-small', loc='upper right')
    #ax.legend(markerscale=1.0, ncol=1, fontsize='xx-small', loc='upper right')

    #if ii<=4:
    #    ax.set_xticklabels([])
    #if ii==0 or ii==10 :
    if ii==0 or ii==4 :#or ii==6  :
        ax.set_ylabel(r'$\Delta\Sigma [{\rm h M_\odot pc^{-2}}]$')
 
    #if not (ii==0 or ii==10) :
    if not (ii==0 or ii==4): #or ii==6) :
        ax.set_yticklabels([])
    #if ii<10:
    #if not (ii==4 or ii==5 or ii==6 or ii==7):
    if ii<=2:
        ax.set_xticklabels([])


    #ax = plt.subplot(4,5,ii+6)
    #ax.plot(uqr,res)
    #ax.set_ylim(-5,5)
    #ax.axhline(0, ls='--', color='grey')
    #ax.set_xscale('log')
    #if ii>10:
    #if ii==4 or ii==5 or ii==6 or ii==7:
    if ii>=3:
        ax.set_xlabel(r'$R[{\rm h^{-1}Mpc}]$')
    #ax.set_xlabel(r'$R[{\rm h^{-1}Mpc}]$')
    #if ii==0 or ii==10 :
    #    ax.set_ylabel(r'$(\Delta \Sigma_{\rm best} - \Delta \Sigma_{\rm meas})/\sigma$')
    #if not (ii==0 or ii==10) :
    #    ax.set_yticklabels([])
    #if ii<10:
    #    ax.set_xticklabels([])
    print('logM0\tlogM1\tgamma1\tgamma2\tsig0\tb0\tb1\tb2\talpsat\tcfac\tap\tpoff\troff\tap\tchisq\n')    
    print(bestp, np.min(pred[:,-1]))

#def get_fsat(b_min,b_max):
#    #bins,fsat,fsaterr = np.loadtxt('/mnt/home/student/cdivya/github/weaklens_pipeline/DataStore/preetish/gama_equatorial_zleq_03_Mrpetro_cen_sat_iso_flags_gr_colour_reduced_sample.fits_fsat.dat', unpack=1)
#    #logMstelarr = -999*np.ones(len(np.unique(bins)))
#    #idx = (bins>=bmin) & (bins<bmax)
#    #fsat = fsat[idx]
#    #fsaterr = 2*fsaterr[idx]
#
#    #for bb in range(bmin, bmax):
#    #    idx = (full_gal['bin']==bb)
#    #    logMstelarr[bb] = float(np.log10(np.mean(10**full_gal['logmstar'][idx])) + 2*np.log10(0.7))
#    #idx =(logMstelarr !=-999) 
#    #logMstelarr = logMstelarr[idx]
#
#    #xx = logMstelarr
#    #yy = fsat
#    #yyerr = fsaterr
# 
#    fsatdat = pd.read_csv('./fsat_zu_2015_tab1.csv', delim_whitespace=1)
#    xx    = (fsatdat['bmin'].values + fsatdat['bmax'].values)/2.0
#    yy    = fsatdat['fsat'].values
#    yyerr = (fsatdat['fsaterr_p'].values + fsatdat['fsaterr_m'].values)/2.0
#
#
#    ii = 15
#    ax = plt.subplot(4,5,ii)
#    ax.errorbar(xx, yy ,yyerr, fmt='.', capsize=3) 
#    y    = pred[:,77:-1]
#
#    chisq = pred[:,-1]
#    ymin0 = np.percentile(y,16,axis=0)
#    ymax0 = np.percentile(y,84,axis=0)
#    ymin1 = np.percentile(y,2.5,axis=0)
#    ymax1 = np.percentile(y,97.5,axis=0)
#   
#    ax.fill_between(xx, ymin1, ymax1, color='#52aae7', alpha = 0.45, zorder=1) # 95 percentile
#    ax.fill_between(xx, ymin0, ymax0, color='#1f77b4', alpha = 0.45, zorder=2) # 68 percentile
#    
#    #best fit model
#    idx = np.argmin(chisq)
#    del_sig = y[idx,:]
#    Delta = y[idx,:] - yy
#    res  = Delta*1.0/yyerr
#    icov = np.diag(1/(yyerr)**2)
#
#    chisq = np.dot(Delta, np.dot(icov, Delta))
#    snr   = (np.dot(yy, np.dot(icov, yy)))**0.5
#    ax.plot(xx, y[idx,:], '-', color='red', label= r"$\chi^2 = %2.2f$"%(np.dot(Delta, np.dot(icov, Delta))))
#
#    ax.text(0.08, 0.05, r"${\rm SNR} = %2.2f$"%(snr), transform=ax.transAxes, fontsize=7)
#    ax.legend(markerscale=1.0, ncol=2, fontsize='xx-small', loc='upper right')
#    ax.set_ylabel(r'$f_{\rm sat}$')
#    ax.set_xticklabels([])
#
#    ax = plt.subplot(4,5,ii+5)
#    ax.plot(xx,res)
#    ax.set_ylim(-5,5)
#    ax.axhline(0, ls='--', color='grey')
#    ax.set_xlabel(r'$\log[M^*/h^{-2}M_\odot]$')
#    ax.set_ylabel(r'$(f^{\rm mod}_{\rm sat} - f^{\rm meas}_{\rm sat})/\sigma$')
#    ax.set_xticklabels([])
 
       
for ii,bb in enumerate(range(0,7)):
    pltsigdata(ii,bb)
    print(bb)

texit()
#get_fsat(0,7)
#plt.tight_layout()
#plt.savefig('./sn_preet_gama_output_cosmo_zu15_full_revised_fits/fits_%s.pdf'%chnname,dpi=300)
plt.savefig('./sn_preet_gama_output_cosmo_zu15_full_revised_fits/fits_full_%s.pdf'%chnname,dpi=300)
#plt.savefig('./sn_preet_gama_output_cosmo_zu15_full_revised_fits/fits_%s.png'%chnname,dpi=300)


