import sys
#sys.path.append("mandel_2005_whichopt/lib/python2.7/site-packages/")
sys.path.append("mandel_2005/lib/python2.7/site-packages/")
#sys.path.append("test_divya/lib/python2.7/site-packages/")
#sys.path.append("divya/lib/python2.7/site-packages/")
from math import log10
#import cosmology as cc
import numpy as np
import hod as h
import matplotlib.pyplot as plt

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
p.Om0 = 0.307115
p.w0 = -1
p.wa = 0
p.Omk = 0.0
p.hval = 0.6777
p.Omb = 0.048206
p.th = 2.726
p.s8 = 0.8228
p.nspec = 0.96
p.ximax = log10(8.0)
p.cfac = 1.0
q = h.hodpars()


#ni = cc.cosmology()

def aum():
#def aum(x,rp,nbin,M):
  #logMmin, siglogM, poff, roff_rs, cfac, alpsat, rsft=x

  q.Mmin = 13
  q.siglogM = 0.1
  
  q.alpsat = 3 # normalization factor for satellites

  ##satellite parameters
  q.Msat = 17.0
  #q.alpsat = 0.0
  q.Mcut = 16.5 #large mcut M > kappa_M_min zero contri from satellite
  q.csbycdm = 1.0 
  q.fac = 1.0 # what is this
  a = h.hod(p, q)
  print 'here'  
  return a
nbin=12
rp = getdblarr(np.linspace(0.2,2,nbin))

m = np.linspace(13,16,10)

x = [12, 0.81030222, 0.5, 0.5, 3.75151809, 1.5, 0.197079]
a = aum()
for cnt,p in enumerate([1,2]):
    #x = [i, 0.81030222, 0.5, 0.5, 0.5, 0.5, 1.5, 0.197079]
    x = [11.73240861,  1.39641369,  0.93582792,  1.0,  0.31538348,
                    0.36457677,  p]

    logMmin, siglogM, alphacen, alphasat, poff, roff_rs, cfac=x
    rsft= 0.1
    q.Mmin = logMmin
    q.siglogM = siglogM
    
    #satellite parameters
    #q.Msat = 14.0
    q.alpsat = alphasat
    q.Mcut = alphacen #large mcut M > kappa_M_min zero contri from satellite
    #q.csbycdm = 1.0 
    #q.fac = 1.0 # what is this
    a.hod_renew(q)
    a.set_cfactor(cfac)
    a.set_cen_offset_params(poff,roff_rs)
    #a.set_whichopt(3)
    esd = getdblarr(np.zeros(nbin))
    a.ESD(rsft,12,rp,esd,16);

    esd = getnparr(esd,nbin)
    mh = np.linspace(9,16,100)
    n_cen = 0.0*mh
    n_sat = 0.0*mh
    avg_n_cen = 0.0
    avg_n_sat = 0.0
 
    for i in range(100):
        n_cen[i] = a.ncen(mh[i])#*a.nofm(10**mh[i],rsft) * 10**mh[i]
        n_sat[i] = a.nsat(mh[i])#*a.nofm(10**mh[i],rsft) * 10**mh[i]
        avg_n_cen += a.ncen(mh[i])*a.nofm(10**mh[i],rsft) * 10**mh[i]
        avg_n_sat += a.nsat(mh[i])*a.nofm(10**mh[i],rsft) * 10**mh[i]
    
    #print avg_n_sat/(avg_n_sat+avg_n_cen) 
    plt.subplot(2,2,1)
    plt.plot(10**mh,n_cen,label='cfac=%s'%p)
    plt.xscale('log')
    plt.subplot(2,2,2)
    plt.plot(10**mh,n_sat,label='cfac=%s'%p)
    plt.xscale('log')
    plt.yscale('log')
    plt.subplot(2,2,3)
    plt.plot(10**mh,n_sat+n_cen,label='cfac%s'%p)
    #plt.plot(10**mh,n_sat,ls='--',label='cfac=%s'%i)
    plt.xscale('log')
    plt.yscale('log')
    plt.subplot(2,2,4)
    plt.plot(getnparr(rp,12),esd,label='cfac=%s'%p)
    plt.xscale('log')
    plt.yscale('log')

plt.legend()
plt.savefig('spooky.png', dpi=300)
plt.show()
