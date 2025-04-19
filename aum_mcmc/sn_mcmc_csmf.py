#!/usr/bin/env python3
import sys
sys.path.append("/mnt/csoft/tools/anaconda2/bin/python")
import numpy as np
import pandas as pd
import emcee
import time
from emcee.utils import MPIPool
import os
sys.path.append("/mnt/home/student/cdivya/github/aum-master/csmf_bosch_2013/lib/python2.7/site-packages")
#import cosmology as cc
import hod as h
from scipy.integrate import quad

from math import log10

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
#fixing cosmology
p = h.cosmo()
p.Om0 = 0.26
#p.Om0 = 0.238
p.w0 = -1
p.wa = 0
p.Omk = 0.0
p.hval = 0.72
p.Omb = 0.0418
p.th = 2.726
p.s8 = 0.77   #check zu and madelbaumm 2015
p.nspec = 0.951
p.ximax = log10(8.0)
p.cfac = 1.0
q = h.hodpars()


def init_aum():
    q.Mmin = 13
    q.siglogM = 0.1
    ##satellite parameters
    q.Msat = 17.0
    q.Mcut = 16.5 #large mcut M > kappa_M_min zero contri from satellite
    q.csbycdm = 1.0 
    q.fac = 1.0 # what is this
    #mod params for csmf
    q.logM0 = 12
    q.logM1 = 11.6
    q.alpha = 3  #gamma1
    q.beta = 0.3 #gamma2
    q.sig0 = 0.4
    q.b0 = -13
    q.b1 = 3
    q.b2 = 1.0
    q.alpsat = -1.5 # nor
    q.logMa = 9
    q.logMb = 10
    return h.hod(p, q)

def stellar(ap, logMstel, rbin): 
    logMstel = logMstel - np.log10(0.72) #stellar function needs inputs in h-1 Msun
    ds = ap*(10**logMstel)/(np.pi*rbin**2)
    ds = ds/1e12
    return ds

def esd(x, rbin, zred, logMa, logMb, logMstel, aum, whichopt=0):
    logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alpsat, cfac, poff, roff , ap = x
    q.logM0 = logM0
    q.logM1 = logM1
    q.alpha = gamma1
    q.beta  = gamma2
    q.sig0 = sig0
    q.b0 = b0
    q.b1 = b1
    q.b2 = b2
    q.alpsat = alpsat
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
    logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alpsat, cfac, poff, roff , ap = x
    if 8.0<=logM0<=11.5 and 9.5<=logM1<=15 and 2<=gamma1<=5 and 0.01<=gamma2<=2 and 0.05<=sig0<=0.5 and -2<=b0<=2 and -2<=b1<=2 and -2<=b2<=2 and -2<=alpsat<=-1  and 0.5<=ap<=5  and 0<=poff<=1.0 and 0.0<=roff<=0.2 and cfac>0:
        return 0.0 + np.log(gauss(cfac,1.0,0.2)) 
    return -np.inf

cnt=0
aumdict={}
def lnprob(x, data, icov, rbin, rsftarr, logMaarr, logMbarr, logMstelarr):
    #logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alpsat, cfac, poff, roff , ap = x
    #logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alpsat, cfac, poff, roff , ap = x
    global cnt, aumdict
    lp = lnprior(x)
    nbins = int(len(data)/len(rbin))
    if not np.isfinite(lp) or np.isnan(lp):
       dirt = 5*np.ones(len(data) + 2*nbins)
       return -np.inf,dirt
    mod   = np.zeros(len(data))
    logMh = np.zeros(nbins)
    fsat  = np.zeros(nbins)
    #print "cnt  sdkjsahdkjhskas", cnt 
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
<<<<<<< HEAD
        print('alert')
        sys.exit(0)
    print('logM0\tlogM1\tgamma1\tgamma2\tsig0\tb0\tb1\tb2\talpsat\tcfac\tap\tpoff\troff\tap\tchisq')
    print((x,chisq))
    print('blob is Del_sig, logMh, fsat, chisq')
    print(( 'size of blob=%d+%d+%d+%d'%(len(mod),len(logMh),len(fsat),len([chisq]))))
=======
        print 'alert'
        sys.exit(0)
    print('logM0\tlogM1\tgamma1\tgamma2\tsig0\tb0\tb1\tb2\talpsat\tcfac\tap\tpoff\troff\tap\tchisq')
    print(x,chisq)
    print('blob is Del_sig, logMh, fsat, chisq')
    print( 'size of blob=%d+%d+%d+%d'%(len(mod),len(logMh),len(fsat),len([chisq])))
>>>>>>> 2cb08947c4825a1e0e38de69838ff8bfaff35728
     
    return res,blob

def runchain(Ntotal,sampler,chainf,blobf,pos):
    blnk=[""];
    fchain=open(chainf,"w");
    fblob=open(blobf,"w");
    iterno=1;
    # Store chainfile and prednfile in the same format as before
    for result in sampler.sample(pos, iterations=Ntotal, storechain=0):
        posn,probn,staten,blobsn = result;
        #posn,probn,staten = result;
<<<<<<< HEAD
        for i in range(nwalkers):
=======
        for i in xrange(nwalkers):
>>>>>>> 2cb08947c4825a1e0e38de69838ff8bfaff35728
            np.savetxt(fchain,posn[i],newline=' ');
            np.savetxt(fchain,[sampler.acceptance_fraction[i],-2.*probn[i]],newline=' ');
            np.savetxt(fblob,blobsn[i],newline=' ');
            np.savetxt(fchain,blnk,fmt='%s');
            np.savetxt(fblob,blnk,fmt='%s');
<<<<<<< HEAD
        print(("Iteration number: %d of %d done"%(iterno,Ntotal)));
=======
        print ("Iteration number: %d of %d done"%(iterno,Ntotal));
>>>>>>> 2cb08947c4825a1e0e38de69838ff8bfaff35728
        iterno=iterno+1;
        posnew=result[0];

    fchain.close();
    fblob.close();
    return posnew;


if __name__ == "__main__":

    from astropy.io import fits
    import healpy as hp
    full_gal = fits.open('../DataStore/preetish/revised_gama_specz_sample_zleq_03.fits')[1].data
    msk = hp.read_map("/mnt/home/student/cdivya/github/weaklens_pipeline/DataStore/data/S16A_mask_w_holes.fits")
    ipix = hp.ang2pix(int(np.sqrt(msk.size/12)), full_gal['ra'], full_gal['dec'], lonlat=1)
    flg  = msk[ipix]
    sel  = (flg==1.0) & (full_gal['logmstar_h2'] < 12.0)
    sel  = sel & (full_gal['dec']>=-6)

    rsftarr     = np.zeros(8)
    logMaarr    = np.zeros(8)
    logMbarr    = np.zeros(8)
    logMstelarr = np.zeros(8)

    data = 0.0*rsftarr
    #joint fit first collecting the measurements
    for i,bb in enumerate(range(0,8)):
        idx            = sel & (full_gal['bin']==bb)
        rsftarr[i]     = np.median(full_gal['z_tonry'][idx])
        logMaarr[i]    = full_gal['logMstar_h2'][idx].min()  
        logMbarr[i]    = full_gal['logMstar_h2'][idx].max() 
        logMstelarr[i] = np.log10(np.mean(10**full_gal['logMstar_h2'][idx]))

        dat  = np.loadtxt('./sn_preet_gama_output_cosmo_zu15_full_revised/dsigma.dat_%d'%bb);
        rdat = np.loadtxt('./sn_preet_gama_rand_output_cosmo_zu15_full_revised/avg_dsig_%d.dat'%bb);
        yy   = dat[:,5]*(1.0/(1+dat[:,14])) - rdat[:,1]
        if i==0:
            data = yy[:10]    
        else:
            yy = yy[:10]            
            data = np.concatenate((data,yy))
    
    rdat = np.loadtxt('./sn_preet_gama_rand_output_cosmo_zu15_full_revised/avg_dsig_0.dat');
    rbin = rdat[:10,0]
    
    #then get the full covariance
    cov = np.loadtxt('./sn_preet_gama_output_cosmo_zu15_full_revised/cov_no_boost.dat',dtype=float,unpack=True)
    icov = np.linalg.inv(cov)
<<<<<<< HEAD
    print(rbin)
=======
    print rbin
>>>>>>> 2cb08947c4825a1e0e38de69838ff8bfaff35728
    njacks = 500
    hartlap_factor = (njacks - len(data) - 2) * 1.0/(njacks - 1)
    icov = hartlap_factor*icov

    # Add in a MPI pool of workers
    pool = MPIPool()
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    ndim = 13
    #nwalkers = 128
    nwalkers = 256
    np.random.seed(1991)
    #def lnprob(x, data, icov, rbin, rsftarr, logMaarr, logMbarr, logMstelarr, aumdict):
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool, args=[data, icov, rbin, rsftarr, logMaarr, logMbarr, logMstelarr])
 

    p_logM0    = np.random.uniform(8.0,11.5,nwalkers)
    p_logM1    = np.random.uniform(10.0,15,nwalkers)
    p_gamma1   = np.random.uniform(2,5,nwalkers)
    p_gamma2   = np.random.uniform(0.01,2.0,nwalkers)
    p_sig0     = np.random.uniform(0.05,0.5,nwalkers)
    p_b0       = np.random.uniform(-1,2,nwalkers)
    p_b1       = np.random.uniform(-2,2,nwalkers)
    p_b2       = np.random.uniform(-1,2,nwalkers)
    p_alpsat   = np.random.uniform(-2,-1,nwalkers)
    p_cfac     = np.random.uniform(0.8,1.2,nwalkers)
    p_ap       = np.random.uniform(0.5,5.0, nwalkers)
    p_poff     = np.random.uniform(0.0,1.0,nwalkers)
    p_roff     = np.random.uniform(0.0,0.2,nwalkers)
 
    #logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alpsat, cfac, poff, roff , ap = x
<<<<<<< HEAD
    p_0 = list(zip(p_logM0, p_logM1, p_gamma1, p_gamma2, p_sig0, p_b0, p_b1, p_b2, p_alpsat, p_cfac, p_poff, p_roff, p_ap))

    print(("Starting iteration, ndim is ",ndim));
=======
    p_0 = zip(p_logM0, p_logM1, p_gamma1, p_gamma2, p_sig0, p_b0, p_b1, p_b2, p_alpsat, p_cfac, p_poff, p_roff, p_ap)

    print ("Starting iteration, ndim is ",ndim);
>>>>>>> 2cb08947c4825a1e0e38de69838ff8bfaff35728
    Nburn=5000
    pos = runchain(Nburn,sampler,"./sn_preet_gama_output_cosmo_zu15_full_revised_fits/burnfile.dat","./sn_preet_gama_output_cosmo_zu15_full_revised_fits/burnpredfile.dat",p_0);
    
    sampler.reset();
    
    Ntotal=10000
    posfinal = runchain(Ntotal,sampler,"./sn_preet_gama_output_cosmo_zu15_full_revised_fits/chainfile.dat","./sn_preet_gama_output_cosmo_zu15_full_revised_fits/predfile.dat", pos);
    pool.close()


