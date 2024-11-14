import sys
sys.path.append("mandle_2005_whichopt/lib/python2.7/site-packages/")
#sys.path.append("test_divya/lib/python2.7/site-packages/")
#sys.path.append("divya/lib/python2.7/site-packages/")
from math import log10
import cosmology as cc
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

ni = cc.cosmology()
a =  ni.nofm(10**(14.1953),0.3329355)
b =  ni.nofm(10**(14.53),0.35)
print a, b, a/b

#sys.exit()

def aum(x,rp,nbin,M):
  logMmin,siglogM,poff,roff_rs,cfac,rsft=x

  q.Mmin = logMmin
  q.siglogM = siglogM
  
  #satellite parameters
  q.Msat = 14.0
  q.alpsat = 0.0
  q.Mcut = 16.5 #large mcut M > kappa_M_min zero contri from satellite
  q.csbycdm = 1.0 
  q.fac = 1.0 # what is this
  a = h.hod(p, q)
  print 'here'  
  #poff and roff/rs initialization
  a.set_cen_offset_params(poff,roff_rs)
  #multiplicative constant to multiply all dark matter concentration
  a.set_cfactor(cfac)
  
  #set incompleteness params
  # inc_alp: Slope for the incompleteness
  # inc_xM: Logarithm of mass above which sample is complete, below a log-linear form with slope inc_alp
  inc_alp,inc_xM = [0.0,12.0]
  a.set_inc_params(inc_alp,inc_xM)
  esd = getdblarr(np.zeros(nbin))
  a.ESD(rsft,12,rp,esd,16);
  
  result = getnparr(esd,nbin)
  value = a.ncen(M)
  print result
  return result

nbin=12
rp = getdblarr(np.linspace(0.2,2,nbin))

m = np.linspace(13,16,10)
#m = [12.34456488]#,12.0,13.0]
x = [14,0.81030222,0.04510043,0.33656621,3.75151809,0.27079]
ncen=0.0*m

#for i in range(len(m)):
#    print 'blah'
#    ncen[i] = aum(x,rp,nbin,m[i])	
#
#print m
#print ncen
#plt.plot(m, ncen)
off = [0.0,0.5,1.0]
for i in off:
    x = [12,0.81030222,i,0.5,3.75151809,0.197079]
    print i
    plt.plot(np.logspace(0.2,2,12), aum(x,rp,nbin,12))
plt.xscale('log')
plt.savefig('spooky.png')
#plt.show()
