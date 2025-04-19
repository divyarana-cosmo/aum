#!/usr/bin/env python3
import sys
sys.path.append("/mnt/csoft/tools/anaconda2/bin/python")
import numpy as np
import pandas as pd
import emcee
import time
from emcee.utils import MPIPool
import os
import numpy as np
from scipy import integrate
sys.path.append("/mnt/home/student/cdivya/github/aum-master/csmf_bosch_2013_loglin_alpsat_pvt_13_beh/lib/python2.7/site-packages/") #r200m in 10**9 - 10**16 halo mass
from math import log10
import cosmology as cc
import numpy as np
import hod as h
from scipy.integrate import quad
from astropy.io import fits

def gauss(x,mean,sig):
    val = np.exp(-(x-mean)**2/(2*sig**2))/np.sqrt(2*np.pi*sig**2)
    return val

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
p.Om0 = 0.238
p.w0 = -1
p.wa = 0
p.Omk = 0.0
p.hval = 0.73
p.Omb = 0.0418
p.th = 2.726
p.s8 = 0.75
p.nspec = 0.951
p.ximax = log10(8.0)
p.cfac = 1.0
q = h.hodpars()

def init_aum():
    q.Mmin = 13
    q.siglogM = 0.1
    ##other params not required but need to initialize
    q.Msat = 17.0
    q.Mcut = 16.5 #large mcut M > kappa_M_min zero contri from satellite
    q.csbycdm = 1.0
    q.fac = 1.0 # what is this
  
    #mod params for csmf
    #parametrization params behroz_2013a
    q.logM0 = -1.7538816143114042    # logeps
    q.logM1 = 11.535270078042663  # logM1
    q.alpha = -1.4224932155336807     # alpha
    q.beta = 3.469797805247476    # delta
    q.alpsat = 0.302032396110954    # gamma
    #central HOD
    q.sig0 = 0.17
    #satellite HOD
    q.alpha_15 = -1.3
    q.eta = 0.0
    q.b0 = 0.19
    q.b1 = 0.83
    q.b2 = -0.02
    #stellar bin limits
    q.logMa = 9
    q.logMb = 10
    a = h.hod(p, q)
    return a

def stellar(ap, logMstel, rbin): 
    ds = ap*(10**logMstel)/(np.pi*rbin**2)
    ds = ds/1e12
    return ds


def esd(x, rbin, zred, logMa, logMb, logMstel, aum, whichopt=0):
    logeps, logM1, alpha, delta, gamma, sig0, b0, b1, b2, alpha_15, eta, cfac, ap, poff, roff = x
    #beh_2013 SMHM params
    q.logM0  = logeps
    q.logM1  = logM1
    q.alpha  = alpha
    q.beta   = delta
    q.alpsat = gamma
    #scatter
    q.sig0 = sig0
    #satellite galaxies params
    q.b0 = b0
    q.b1 = b1
    q.b2 = b2
    q.alpha_15 = alpha_15
    q.eta = eta
    #bin limits
    q.logMa = logMa
    q.logMb = logMb

    aum.hod_renew(q)
    aum.set_cfactor(cfac)
    aum.set_cen_offset_params(poff, roff)
    nbin = int(np.size(rbin))
    dbin = nbin+4    
    esd = getdblarr(np.zeros(nbin))
    
    aum.set_whichopt(whichopt)
    aum.ESD(zred, nbin, getdblarr(rbin), esd, dbin);
    result = getnparr(esd,nbin) 
    halo_mass = aum.avmass_cen(zred)*1e12
    if whichopt==0:
        avgnsat = quad(lambda y: aum.nsat(y) * aum.nofm(10**y,zred) * 10**y, 9, 16)[0]
        avgncen = quad(lambda y: aum.ncen(y) * aum.nofm(10**y,zred) * 10**y, 9, 16)[0]
        try:
            fsat =  avgnsat*1.0/(avgnsat + avgncen)
        except ZeroDivisionError:
            result[:] = np.NaN
            fsat = -99
        return result + stellar(ap, logMstel, rbin), np.log10(halo_mass), fsat
   
    return result

def lnprior(x):
    logeps, logM1, alpha, delta, gamma, sig0, b0, b1, b2, alpha_15, eta, cfac, ap, poff, roff = x

    if -2<=logeps<=0.0 and 10.0<=logM1<=12 and -2<=alpha<=0 and 1<=delta<=4 and 0.0<=gamma<=1.0 and 0.01<=sig0<=1.5 and -2<=b0<=2 and -5<=b1<=5 and -2<=b2<=2 and -2<=alpha_15<=0 and  -1.0<=eta<=0.0 and 0.5<=ap<=5  and 0<=poff<=1.0 and 0.0<=roff<=0.2 and cfac>0:
        return 0.0 + np.log(gauss(cfac,1.0,0.2)) 
    return -np.inf

cnt=0
aumdict={}
def lnprob(x, data, icov, rbin, rsftarr, logMaarr, logMbarr, logMstelarr):
    global cnt, aumdict
    lp = lnprior(x)
    nbins = int(len(data)/len(rbin))
    if not np.isfinite(lp) or np.isnan(lp):
       dirt = 5*np.ones(len(data) + 2*nbins)
       return -np.inf,dirt
    mod   = np.zeros(nbins*len(rbin))
    logMh = np.zeros(nbins)
    fsat  = np.zeros(nbins)
    print("cnt  sdkjsahdkjhskas", cnt) 
    if cnt==0:
        for bb in range(len(rsftarr)):
            aumdict['%d'%bb] = init_aum()
        cnt = cnt + 1    
   
    for ii in range(len(rsftarr)):
        logMa    = logMaarr[ii]
        logMb    = logMbarr[ii]
        logMstel = logMstelarr[ii]
        zred     = rsftarr[ii]
        aum      = aumdict['%d'%ii]

        pred     = esd(x, rbin, zred, logMa, logMb, logMstel, aum=aum) 
        x0 = (ii)*len(rbin)
        x1 = (ii+1)*len(rbin)
        mod[x0:x1] = pred[0]          
        logMh[ii]  = pred[1]
        fsat[ii]   = pred[2]
        #print zred, logMa, logMb, logMstel, pred[0]
    Delta = mod - data
    chisq = np.dot(Delta, np.dot(icov, Delta))
    blob  = np.append(mod,np.append(logMh,fsat))
    blob  = np.append(blob,chisq)
  
    if np.isnan(chisq):
        res = -np.inf
        return res,blob
    else:
        res = lp-chisq*0.5

    if chisq<0:
        print('alert')
        sys.exit(0)
    # logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alpsat, cfac, poff, roff , ap = x
    print('logM0\tlogM1\tgamma1\tgamma2\tsig0\tb0\tb1\tb2\talpsat\teta\tcfac\tap\tpoff\troff\tap\tchisq')
    print((x,chisq))
    print('blob is Del_sig, logMh, fsat, chisq')
    print(( 'size of blob=%d+%d+%d+%d'%(len(mod),len(logMh),len(fsat),len([chisq]))))
     
    return res,blob

def runchain(Ntotal,sampler,chainf,blobf,pos):
    blnk=[""];
    fchain=open(chainf,"w");
    fblob=open(blobf,"w");
    iterno=1;
    # Store chainfile and prednfile in the same format as before
    for result in sampler.sample(pos, iterations=Ntotal, storechain=0):
        posn,probn,staten,blobsn = result;
        sys.stdout.flush()
        sys.stderr.flush()
 
        #posn,probn,staten = result;
        for i in range(nwalkers):
            np.savetxt(fchain,posn[i],newline=' ');
            np.savetxt(fchain,[sampler.acceptance_fraction[i],-2.*probn[i]],newline=' ');
            np.savetxt(fblob,blobsn[i],newline=' ');
            np.savetxt(fchain,blnk,fmt='%s');
            np.savetxt(fblob,blnk,fmt='%s');
        print(("Iteration number: %d of %d done"%(iterno,Ntotal)));
        iterno=iterno+1;
        posnew=result[0];

    fchain.close();
    fblob.close();
    return posnew;


if __name__ == "__main__":

    from astropy.io import fits
    full_gal = fits.open('../DataStore/preetish/gama_equatorial_zleq_03_Mrpetro_cen_sat_iso_flags_gr_colour_reduced_sample.fits')[1].data
    
    bins =  full_gal['bin']
    rsftarr     = -999*np.ones(len(np.unique(bins)))
    logMaarr    = -999*np.ones(len(np.unique(bins)))
    logMbarr    = -999*np.ones(len(np.unique(bins)))
    logMstelarr = -999*np.ones(len(np.unique(bins)))

    data = 0.0*rsftarr
    rbin = 0.0*rsftarr
    #joint fit first collecting the measurements
    for i,bb in enumerate(range(2,9)):
        idx = (bins==bb)
        rsftarr[i]     = np.median(full_gal['z_tonry'][idx])
        logMaarr[i]    = full_gal['logmstar'][idx].min() + 2*np.log10(0.7) 
        logMbarr[i]    = full_gal['logmstar'][idx].max() + 2*np.log10(0.7)
        logMstelarr[i] = np.log10(np.mean(10**full_gal['logmstar'][idx])) + 2*np.log10(0.7)

        dat  = np.loadtxt('./sn_preet_gama_output_cosmo_y08_full/dsigma.dat_%d'%bb);
        rdat = np.loadtxt('./sn_preet_gama_rand_output_cosmo_y08_full/avg_dsig_%d.dat'%bb);
        if i==0:
            rbin = np.unique(dat[:,7])
            yy = dat[:,5]*(1.0/(1+dat[:,14])) - rdat[:,1]
            data = yy    
        else:
            yy = dat[:,5]*(1.0/(1+dat[:,14])) - rdat[:,1]
            data = np.concatenate((data,yy))
    
    idx = (rsftarr>=0) 
    rsftarr     = rsftarr[idx]
    logMaarr    = logMaarr[idx]
    logMbarr    = logMbarr[idx]
    logMstelarr = logMstelarr[idx]

    #then get the full covariance
    cov = np.loadtxt('./sn_preet_gama_output_cosmo_y08_full/cov_no_boost.dat',dtype=float,unpack=True)
    #cov = cov[len(rbin):, len(rbin):] 


    icov = np.linalg.inv(cov)
    
    njacks = 500
    hartlap_factor = (njacks - len(data) - 2) * 1.0/(njacks - 1)

    icov = hartlap_factor*icov

    # Add in a MPI pool of workers
    pool = MPIPool()
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    ndim = 15
    nwalkers = 256
    np.random.seed(1996)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool, args=[data, icov, rbin, rsftarr, logMaarr, logMbarr, logMstelarr])


    #logeps, logM1, alpha, delta, gamma, sig0, b0, b1, b2, alpha_15, eta, cfac, ap, poff, roff = x

    p_logeps    = np.random.uniform(-2, 0.0, nwalkers)
    p_logM1     = np.random.uniform(10.0,12,nwalkers)
    p_alpha     = np.random.uniform(-2,0.0,nwalkers)
    p_delta     = np.random.uniform(1.0,4.0,nwalkers)
    p_gamma     = np.random.uniform(0.0,1.0,nwalkers)
    p_sig0      = np.random.uniform(0.05,1.5,nwalkers)
    p_b0        = np.random.uniform(-2,2,nwalkers)
    p_b1        = np.random.uniform(-5,5,nwalkers)
    p_b2        = np.random.uniform(-2,2,nwalkers)
    p_alpha_15  = np.random.uniform(-2,0.0,nwalkers)
    p_eta       = np.random.uniform(-1.0,0.0,nwalkers)
    p_cfac      = np.random.uniform(0.8,1.2,nwalkers)
    p_ap        = np.random.uniform(0.5,5.0, nwalkers)
    p_poff      = np.random.uniform(0.0,1.0,nwalkers)
    p_roff      = np.random.uniform(0.0,0.2,nwalkers)
 
    p_0 = list(zip(p_logeps, p_logM1, p_alpha, p_delta, p_gamma, p_sig0, p_b0, p_b1, p_b2, p_alpha_15, p_eta, p_cfac, p_ap, p_poff, p_roff))

    
    print(("Starting iteration, ndim is ",ndim));
    Nburn=5000
    pos = runchain(Nburn,sampler,"./sn_preet_gama_output_cosmo_y08_full_fits/burnfile.dat_bosch_2013_loglin_alpsat_pvt_13_beh","./sn_preet_gama_output_cosmo_y08_full_fits/burnpredfile.dat_bosch_2013_loglin_alpsat_pvt_13_beh",p_0);
    
    sampler.reset();
    
    Ntotal=10000
    posfinal = runchain(Ntotal,sampler,"./sn_preet_gama_output_cosmo_y08_full_fits/chainfile.dat_bosch_2013_loglin_alpsat_pvt_13_beh","./sn_preet_gama_output_cosmo_y08_full_fits/predfile.dat_bosch_2013_loglin_alpsat_pvt_13_beh", pos);
    pool.close()


