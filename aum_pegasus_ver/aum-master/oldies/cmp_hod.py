import sys
import matplotlib.pyplot as plt
import numpy as np

def run_hod(seln, dirpath):
    if seln==0:
        sys.path.append(dirpath)
        from math import log10
        import hod as h
        
    else:
        print 'newone'
        del sys.modules["hod"]
        sys.path.append(dirpath)
        from math import log10
        import hod as h
        

    print h.__file__ 
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
    
    if seln==0:
        x = [12, 0.101682050, 1e2, 0.613934100, 0.0738630841, 1.61797865]
        logMmin, siglogM, normnsat, poff, roff_rs, cfac=x
        rsft=0.1
        q.Mmin = logMmin
        q.siglogM = siglogM
        q.alpsat = normnsat

        q.Msat = 14.0
        q.Msat = 17.0
        q.Mcut = 16.5 #large mcut M > kappa_M_min zero contri from satellite
        q.csbycdm = 1.0 
        q.fac = 1.0 # what is this

        a = h.hod(p,q)
        a.set_cfactor(cfac)
        a.set_cen_offset_params(poff,roff_rs)
        label = 'mandel2005'
    else:
        #logMmin, siglogM, alphasat, kappa, msat = x
        x = [10.49, 0.1017, 4.0, 1.0, 14]
        logMmin, siglogM, alphasat, kappa, msat = x
        rsft= 0.1
        q.Mmin = logMmin
        q.siglogM = siglogM
        q.Msat = msat
        q.alpsat = alphasat
        q.Mcut = np.log10(kappa) + q.Mmin 

        q.csbycdm = 1.0 
        q.fac = 1.0 # what is this
        a = h.hod(p,q)
 
        #a.hod_renew(q)
        a.set_cfactor(1.0)
        a.set_cen_offset_params(0.0,0.1)
        label = 'modified std HOD'

    mh = np.linspace(9,16,100)
    n_cen = 0.0*mh
    n_sat = 0.0*mh
    
    for i in range(100):
        n_cen[i] = a.ncen(mh[i])*a.nofm(10**mh[i],0.1)
        n_sat[i] = a.nsat(mh[i])*a.nofm(10**mh[i],0.1)

    plt.subplot(2,2,1)
    plt.plot(10**mh, n_cen, label=label)
    plt.xlabel(r'$M$')
    plt.ylabel(r'$N_{\rm cen} * n(M)$')
    plt.xscale('log')
    plt.legend()

    plt.subplot(2,2,2)
    plt.plot(10**mh, n_sat)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$M$')
    plt.ylabel(r'$N_{\rm sat} * n(M)$')
    
    plt.subplot(2,2,3)
    plt.plot(10**mh, n_sat+n_cen)
    plt.xlabel(r'$M$')
    plt.ylabel(r'$N_{\rm tot} * n(M)$')
    plt.xscale('log')
    plt.yscale('log')
    plt.subplot(2,2,4)
    plt.plot(10**mh, n_sat/(n_sat + n_cen))
    plt.ylabel(r'$f_{\rm off}$')
    plt.xlabel(r'$M$')
    plt.xscale('log')
    #plt.yscale('log')
    sys.path.remove(dirpath)




run_hod(0,"mandel_2005/lib/python2.7/site-packages/")
run_hod(1,"mod_std_hod/lib/python2.7/site-packages/")

plt.tight_layout()
plt.savefig('spooky.png', dpi=300)
plt.show()
