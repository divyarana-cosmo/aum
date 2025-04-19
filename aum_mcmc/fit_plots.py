import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import corner
from astropy.io import fits
import sys
#from mcmc import stellar

plt.figure(figsize=[10,7.5])
#plt.figure(figsize=[8.8,6.6])
sys.path.append("/mnt/home/student/cdivya/github/aum-master/mandel_2005_whichopt/lib/python2.7/site-packages/") #r200m in 10**9 - 10**16 halo mass

from math import log10
import cosmology as cc
import numpy as np
import hod as h
from scipy.integrate import quad

def getdblarr(r):
    temp=h.doubleArray(r.size)
    for i in range(r.size):
        temp[i]=r[i]
    return temp

def getnparr(r,n):
    temp=np.zeros(n)
    for i in range(n):
        temp[i]=r[i]
    return temp

p = h.cosmo()
p.Om0 = 0.31
p.w0 = -1
p.wa = 0
p.Omk = 0.0
p.hval = 0.67386
p.Omb = 0.04856
p.th = 2.726
p.s8 = 0.829
p.nspec = 0.9603
p.ximax = log10(8.0)
p.cfac = 1.0
q = h.hodpars()


def aum():
  q.Mmin = 13
  q.siglogM = 0.1
  q.alpsat = 3 # normalization factor for satellites
  ##satellite parameters
  q.Msat = 17.0
  q.Mcut = 16.5 #large mcut M > kappa_M_min zero contri from satellite
  q.csbycdm = 1.0 
  q.fac = 1.0 # what is this
  a = h.hod(p, q)
  return a

a = aum()

data = fits.open("../DataStore/preetish/data_table.fits")[1].data
full_gal = fits.open('../DataStore/preetish/pdr2_zleq_03_frm_16A_overlap_blend_mass_cut.fits')[1].data

def model(x, rbin, whichopt, binno, a=a):
    #logMmin, siglogM, alpha, poff, roff, cfac, logMstel = x
    logMmin, siglogM, alphasat, poff, roff, cfac = x
    #logMmin, siglogM, alphasat, poff, roff, cfac, ap = x

    sel = (data['bin']==binno) #& (z>0) & (dec<50) #also removing aegis
    idx = (full_gal['logM']> data['logM'][sel].min()) & (full_gal['logM']<= data['logM'][sel].max()) & (full_gal['dec']<50) & (full_gal['photoz_best']>0)
    
    rsft = np.mean(full_gal['photoz_best'][idx])


    q.Mmin = logMmin
    q.siglogM = siglogM
    #q.Mcut = alphacen
    q.alpsat = alphasat

    a.hod_renew(q)
    a.set_cfactor(cfac)
    a.set_cen_offset_params(poff,roff)

    nbin = int(np.size(rbin))
    dbin = nbin+4    
    a.set_cfactor(cfac)
    a.set_whichopt(whichopt)

    esd = getdblarr(np.zeros(nbin))
    a.ESD(rsft,nbin,getdblarr(rbin),esd,dbin)
    result = getnparr(esd,nbin)
    return result


def pltsigdata(ii):
    ax=plt.subplot(3,4,ii)

    dat= np.loadtxt('./full_jk_output_1/dsigma.dat_%d'%ii);
    rdat= np.loadtxt('./full_jk_random_1/dsigma.dat_%d'%ii);

    y = dat[:,5]*(1.0/(1+dat[:,14])) - rdat[:,5]*(1.0/(1+rdat[:,14]))
    #y = dat[:,5]*(1.0*dat[:,6]/(rdat[:,6]*(1+dat[:,14]))) - rdat[:,5]*(1.0/(1+rdat[:,14]))
    uqr = np.unique(dat[:,7])
    avgy = 0.0*uqr
    yerr = 0.0*uqr
    njacks = int(len(y)/len(uqr))
    for jj,rr in enumerate(uqr):
        idx = (dat[:,7]==rr)
        avgy[jj] = np.mean(y[idx])
        yerr[jj] = np.sqrt(njacks - 1) * np.std(y[idx])  
    print(ii,yerr)
    #if ii>=11:
    #    uqr = uqr[2:]
    #    avgy = avgy[2:]
    #    yerr = yerr[2:]
    #    ax.errorbar(uqr, avgy, yerr, fmt='.', label='bin=%d'%ii, capsize=3, zorder=10)
    #else:
    #    ax.errorbar(uqr, avgy, yerr, fmt='.', label='bin=%d'%ii, capsize=3, zorder=10)
    ax.errorbar(uqr, avgy, yerr, fmt='.', label='bin=%d'%ii, capsize=3, zorder=10)

    chn = pd.read_csv('./full_jk_fits_1/new_x9_16_off_chainfile_%s.dat'%ii,delim_whitespace=True,header=None)
    pred = pd.read_csv('./full_jk_fits_1/new_x9_16_off_predfile_%s.dat'%ii,delim_whitespace=True,header=None)
 
    if ii>=8:
        print("cutting chains")
        chisq_cut = np.percentile(pred.values[:,-1],70)
        idx  = (pred.values[:,-1] < chisq_cut)
        y = pred.values[idx,:len(uqr)]
        chisq = pred.values[idx,-1]
        chn = chn.values[idx,:]
    else:
        y = pred.values[:,:len(uqr)]
        chisq = pred.values[:,-1]
        chn = chn.values[:,:]

    dof = len(uqr) - 6
     
    ymin0 = np.percentile(y,16,axis=0)
    ymax0 = np.percentile(y,84,axis=0)
    
    ymin1 = np.percentile(y,2.5,axis=0)
    ymax1 = np.percentile(y,97.5,axis=0)
   
    ax.fill_between(uqr, ymin1, ymax1, color='#52aae7', alpha = 0.45, zorder=1) # 95 percentile
    ax.fill_between(uqr, ymin0, ymax0, color='#1f77b4', alpha = 0.45, zorder=2) # 68 percentile
    
    #best fit model
    idx = np.argmin(chisq)
    del_sig = y[idx,:]
    bestp = chn[idx,:6]

    #ax.plot(uqr, y[idx,:], '-', color='red', label= r"$\chi^2_{\rm min}/{\rm dof} = %2.2f/%d$"%(chisq[idx],dof))
    ax.plot(uqr, y[idx,:], '-', color='red', label= r"$\chi^2_{\rm red} = %2.2f/%d$"%(chisq[idx],dof))
    ax.plot(uqr, model(bestp, uqr, whichopt=1, binno=ii), ls='--', label= "1h-c")
    ax.plot(uqr, model(bestp, uqr, whichopt=2, binno=ii), ls='-.', label= "1h-s")
    ax.plot(uqr, model(bestp, uqr, whichopt=3, binno=ii), ls=':', label= "2h")

    from astropy.io import fits
    #samp = fits.open('../DataStore/preetish/data_table.fits')[1].data;
    #logMstel = np.mean(samp['logM'][(samp['bin']==ii)&(samp['photoz_best']>0)&(samp['logM']>0)])
    
    #data = fits.open("../DataStore/preetish/data_table.fits")[1].data
    #sel = (data['bin']==ii) #& (z>0) & (dec<50) #also removing aegis
                      
    #full_gal = fits.open('../DataStore/preetish/pdr2_zleq_03_frm_16A_overlap_blend_mass_cut.fits')[1].data
    #                           
    #fidx = (full_gal['logM']> data['logM'][sel].min()) & (full_gal['logM']<= data['logM'][sel].max()) & (full_gal['dec']<50) & (full_gal['photoz_best']>0)

    #rsft = np.mean(full_gal['photoz_best'][fidx])
    #logMstel = np.mean(full_gal['logM'][fidx])


    #logMstel = np.log10(10**logMstel)
    #ax.plot(uqr, stellar(1.0, logMstel, uqr), '--', color='black', label=r'$\log M_{*}$=(%2.2f)'%(logMstel))
    
    #if ii==2:
    #    ax.set_title(r'$\log M_{*}$ is in units of $h^{-1}M_{\odot}$')
    #if ii<=4:
    #    ax.set_ylim(0.05,100)
 
    #if 5<=ii<=8:
    #    ax.set_ylim(0.1,500)
 
    #if ii>8:
    #    ax.set_ylim(2.0,700)
    ax.set_ylim(0.05,4500)
    if ii>8:
        ax.set_xlabel(r'$R[{\rm h^{-1}Mpc}]$')
 
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.legend(markerscale=1.0, fontsize='x-small', loc='upper right',  ncol = 2 )

    #if not (ii==1 or ii==5 or ii==9):
    #    ax.set_yticklabels([])
    if ii==1 or ii==5 or ii==9:
        ax.set_ylabel(r'$\Delta\Sigma [{\rm h M_\odot pc^{-2}}]$')
 
for i in range(1,13):
    #for i in range(1,13):
    pltsigdata(i)
    print(i)

#plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.savefig('./full_jk_fits_1/new_fits_no_bary.png', dpi=600)


