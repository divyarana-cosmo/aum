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
#sys.path.append("/mnt/home/student/cdivya/github/aum-master/mandel_2005/lib/python2.7/site-packages/") #r200m in 10**9 - 10**16 halo mass
#sys.path.append("/mnt/home/student/cdivya/github/aum-master/mod_csmf_hod/lib/python2.7/site-packages/")
sys.path.append("/mnt/home/student/cdivya/github/aum-master/mod_csmf_hod_b0_pvt_13/lib/python2.7/site-packages/")
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

#p.Om0 = 0.31
#p.w0 = -1
#p.wa = 0
#p.Omk = 0.0
#p.hval = 0.67386
#p.Omb = 0.04856
#p.th = 2.726
#p.s8 = 0.829
#p.nspec = 0.9603
#p.ximax = log10(8.0)
#p.cfac = 1.0
q = h.hodpars()

def stellar(ap, logMstel, rbin): 
    ds = ap*(10**logMstel)/(np.pi*rbin**2)
    ds = ds/1e12
    return ds

def aum():
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
  q.alpha = 3
  q.beta = 0.3
  q.sig0 = 0.4
  q.b0 = -13
  q.b1 = 3
  q.b2 = 1.0
  q.alpsat = 3 # nor
  q.logMa = 9
  q.logMb = 10
  a = h.hod(p, q)
  return a

aux = aum()

def model(x, rbin, rsft, logMa, logMb, whichopt=0,  ina=aux):
    logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alphas, cfac, ap, poff, roff = x
    q.logM0 = logM0
    q.logM1 = logM1
    q.alpha = gamma1
    q.beta  = gamma2
    q.sig0 = sig0
    q.b0 = b0
    q.b1 = b1
    q.b2 = b2
    q.alpsat = alphas
    q.logMa = logMa
    q.logMb = logMb
    ina.hod_renew(q)
    ina.set_cfactor(cfac)
    ina.set_cen_offset_params(poff, roff)
    nbin = int(np.size(rbin))
    dbin = nbin+4    
    esd = getdblarr(np.zeros(nbin))
    
    ina.set_whichopt(whichopt)
    ina.ESD(rsft, nbin, getdblarr(rbin), esd, dbin);
    result = getnparr(esd,nbin) 
    halo_mass = ina.avmass_cen(rsft)*1e12
    if whichopt==0:
        avgnsat = quad(lambda y: ina.nsat(y) * ina.nofm(10**y,rsft) * 10**y, 9, 16)[0]
        avgncen = quad(lambda y: ina.ncen(y) * ina.nofm(10**y,rsft) * 10**y, 9, 16)[0]
        try:
            fsat =  avgnsat*1.0/(avgnsat + avgncen)
        except ZeroDivisionError:
            result[:] = np.NaN
            fsat = -99
        return result, np.log10(halo_mass), fsat
   
    return result

def lnprior(x):
    logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alphas, cfac, ap, poff, roff = x
    #logM0, logM1, alpha, beta, sig0, b0, b1, b2, alphas, cfac, poff, roff = x
    #if 9.0<=logM0<=11 and 11.0<=logM1<=13 and 1<=gamma1<=5 and 0<=gamma2<=1 and 0.05<=sig0<=1.5 and -5<=b0<=5 and -5<=b1<=5 and -1<=b2<=1 and -2<=alphas<=0:# and 0<=poff<=1.0 and 0.0<=roff<=0.5:
    if 9.0<=logM0<=12 and 11.0<=logM1<=13 and 1<=gamma1<=5 and 0<=gamma2<=1 and 0.05<=sig0<=1.5 and -5<=b0<=5 and -5<=b1<=5 and -1<=b2<=1 and -2<=alphas<=0 and 0.5<=ap<=5  and 0<=poff<=1.0 and 0.0<=roff<=0.5:
        return 0.0 + np.log(gauss(cfac,1.0,0.2))
    return -np.inf


def lnprob(x, data, icov, rbin, rsft):
    #logM0, logM1, alpha, beta, sig0, b0, b1, b2, alphas = x
    #logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alphas, cfac, ap = x
    lp = lnprior(x)
    nbins = int(len(data)/len(rbin))
    #print nbins
    if not np.isfinite(lp) or np.isnan(lp):
       dirt = 5*np.ones(len(data) + 2*nbins)
       return np.inf
       #return -np.inf,dirt
    mod   = np.zeros(nbins*len(rbin))
    logMh = np.zeros(nbins)
    fsat  = np.zeros(nbins)
    
    bins,mstel_min,mstel_max,mstel_avg,z_mean,z_med = np.loadtxt("../DataStore/preetish/mstel_bins.dat", unpack=True)
    #for ii,bb in enumerate(range(3,13)):
    for ii,bb in enumerate(range(6,13)):
        sel = (bins==bb) 
        logMa    = float(mstel_min[sel])
        logMb    = float(mstel_max[sel])
        logMstel = mstel_avg[sel] - np.log10(0.73)
        #logMstel = mstel_avg[sel] - np.log10(0.67386)
        pred = model(x, rbin, rsft, logMa, logMb) 
        x0 = (ii)*len(rbin)
        x1 = (ii+1)*len(rbin)
        mod[x0:x1] = pred[0] + stellar(x[10], logMstel, rbin)
        logMh[ii]  = pred[1]
        fsat[ii]   = pred[2]

    Delta = mod - data
    chisq = np.dot(Delta, np.dot(icov, Delta))
    #blob  = np.append(mod,np.append(logMh,fsat))
    #blob  = np.append(blob,chisq)
   
    if np.isnan(chisq):
        res = np.inf
        return res#,blob
    else:    
        res = chisq
        #res = lp-chisq*0.5
    if chisq<0:
        print 'alert'
        sys.exit(0)
    
    #print('logM0\tlogM1\talpha\tbeta\tsig0\tb0\tb1\tb2\talphas\tcfac\tpoff\troff\tchisq')
    #print('logM0\tlogM1\tgamma1\tgamma2\tsig0\tb0\tb1\tb2\talphas\tcfac\tap\tpoff\troff\tlogMh_cen\tfsat\tchisq')
    #print(x,logMh,fsat,chisq)
    #print('blob is Del_sig, logMh, fsat, chisq')
    #print( 'size of blob=%d+%d+%d+%d'%(len(mod),len(logMh),len(fsat),len([chisq])))
    print res, np.isscalar(res)
    return res#,blob

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
        for i in xrange(nwalkers):
            np.savetxt(fchain,posn[i],newline=' ');
            np.savetxt(fchain,[sampler.acceptance_fraction[i],-2.*probn[i]],newline=' ');
            np.savetxt(fblob,blobsn[i],newline=' ');
            np.savetxt(fchain,blnk,fmt='%s');
            np.savetxt(fblob,blnk,fmt='%s');
        print ("Iteration number: %d of %d done"%(iterno,Ntotal));
        iterno=iterno+1;
        posnew=result[0];

    fchain.close();
    fblob.close();
    return posnew;


if __name__ == "__main__":
    for bb,i in enumerate(range(6,13)):
        dat  = np.loadtxt('./re_full_jk_output_1_cosmo_y08/dsigma.dat_%d'%i)
        rdat = np.loadtxt('./re_full_jk_random_1_cosmo_y08/dsigma.dat_%d'%i)


        rbin = np.unique(dat[:,7])
        y = 0.0*rbin
        if bb==0:
            yy = dat[:,5]*(1.0/(1+dat[:,14])) - rdat[:,5]*(1.0/(1+rdat[:,14]))
            data = 0.0*rbin
            for jj,rr in enumerate(rbin):
                idx = (dat[:,7]==rr)
                y[jj] = np.mean(yy[idx])
            data = y    
        else:
            yy = dat[:,5]*(1.0/(1+dat[:,14])) - rdat[:,5]*(1.0/(1+rdat[:,14]))
            for jj,rr in enumerate(rbin):
                idx = (dat[:,7]==rr)
                y[jj] = np.mean(yy[idx])
            data = np.concatenate((data,y))
    
    rsft = 0.24 #fixing the redshift to median value for the sample
    rbin = np.unique(dat[:,7])

    cov = np.loadtxt('./re_full_jk_output_1_cosmo_y08/cov_w_rand_subs_no_boost.dat',dtype=float,unpack=True)
    cov = cov[30:, 30:] #just last 5 bins
    icov = np.linalg.inv(cov)
    

    def _lnprob(x):
        res = lnprob(x, data, icov, rbin, rsft)
        print res, np.isscalar(res)
        return res
        

    chn  = pd.read_csv('./re_full_jk_fits_csmf_1_cosmo_y08/chainfile.dat_lt_5_w_cfac_ap_poff_b0_pvt_13',delim_whitespace=True,header=None)
    pred = pd.read_csv('./re_full_jk_fits_csmf_1_cosmo_y08/predfile.dat_lt_5_w_cfac_ap_poff_b0_pvt_13',delim_whitespace=True,header=None)
    
    idx = np.argmin(pred.values[:,-1])
    x0 = [chn.values[idx,:13]]
    from scipy.optimize import minimize

    samp = minimize(_lnprob , x0=x0)

    print samp.x
 
   


