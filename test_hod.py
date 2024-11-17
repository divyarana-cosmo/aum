import sys
import numpy as np
sys.path.append('/Users/divyarana/github/aum/install/lib/python3.13/site-packages/')
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
q = h.hodpars()
p.Om0 = 0.307115
p.w0 = -1
p.wa = 0
p.Omk = 0.0
p.hval = 0.6777
p.Omb = 0.048206
p.th = 2.726
p.s8 = 0.8228
p.nspec = 0.96
p.ximax = np.log10(8.0)
p.cfac = 1.0
q.Mmin = 13.0
q.siglogM = 0.5
q.Msat = 14.0
q.alpsat = 1.0
q.Mcut = 13.5
q.csbycdm = 1.0
q.fac = 1.0
a = h.hod(p, q)


zred=0.1
nbin=10
esd = getdblarr(np.zeros(nbin))
rbin = np.logspace(-1,1,nbin)
a.ESD(zred, nbin, getdblarr(rbin), esd, nbin+4)
result = getnparr(esd,nbin)

plt.subplot(2,2,1)
plt.plot(rbin,result)
plt.xscale('log')
plt.yscale('log')
plt.savefig('test.png', dpi=300)
#help(a)
