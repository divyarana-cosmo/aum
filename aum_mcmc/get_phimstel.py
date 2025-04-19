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
from scipy.integrate import quad,simps,fixed_quad

from math import log10
import matplotlib.pyplot as plt

hval = 0.73
#fixing cosmology
p = h.cosmo()
p.Om0 = 0.238
p.w0 = -1
p.wa = 0
p.Omk = 0.0
p.hval = hval
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

def smhm(xh, logM0, logM1, gamma1, gamma2):
    ans = logM0 + gamma1*(xh - logM1) - (gamma1 - gamma2)*np.log10(1 + 10**(xh - logM1))
    return ans

def phi_cen(x, logmstel, xh):
    #logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alpsat, cfac, poff, roff , ap = x
    logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alpsat, eta, cfac, poff, roff, ap = x
    logMc = smhm(xh, logM0, logM1, gamma1, gamma2) 
    #ans = 1/(np.sqrt(2*np.pi)*sig0) * np.exp(-(logmstel - logMc)**2/(2*sig0**2))#*1/10**logmstel
    ans = np.log10(np.e)/(np.sqrt(2*np.pi)*sig0) * np.exp(-(logmstel - logMc)**2/(2*sig0**2))*1/10**logmstel
    return ans

def phi_sat(x, logmstel, xh):
    #logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alpsat, cfac, poff, roff , ap = x
    logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alpsat, eta, cfac, poff, roff, ap = x
    logMs = smhm(xh, logM0, logM1, gamma1, gamma2) + np.log10(0.562)
    phis  = 10**(b0 + b1 * (xh-12) + b2 * (xh-12)**2)
    ans   = phis*10**((logmstel - logMs)*(alpsat+1)) * np.exp(-10**((logmstel - logMs)*2)) * 1/10**logmstel
    return ans


aum = init_aum()

def phi_mstel(x, logMstelarr, zred):
    logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alpsat, cfac, poff, roff , ap = x
    anscen = 0.0*logMstelarr
    anssat = 0.0*logMstelarr
    for ii in range(len(logMstelarr)):
        logmstel = logMstelarr[ii]

        anscen[ii] = np.log(10)**2* 10**logmstel * quad(lambda y: phi_cen(x, logmstel, y) * aum.nofm(10**y,zred) * 10**y, 9, 16)[0]
        anssat[ii] = np.log(10)**2* 10**logmstel * quad(lambda y: phi_sat(x, logmstel, y) * aum.nofm(10**y,zred) * 10**y, 9, 16)[0]
 
    return anscen , anssat 



if __name__ == "__main__":
    chnpath = '/scratch/cdivya/chainfile_bosch_2013_loglin_alpsat_uni_const_alpha_no_off.dat'  #chainfile.dat'
    predpath = '/scratch/cdivya/predfile_bosch_2013_loglin_alpsat_uni_const_alpha_no_off.dat' #predfile.dat'

    pred    = pd.read_csv('%s'%predpath,  delim_whitespace=True)
    chn     = pd.read_csv('%s'%chnpath,    delim_whitespace=True)
    idx     = (pred.values[:,-1]<300)
    chn     = chn.values[idx,:]
    mh      = np.linspace(11,16,20)
    mat     = np.zeros((len(chn), len(mh)))
    
    for i in range(len(chn[:,0])):
        #logM0, logM1, gamma1, gamma2 = chn[i,:4]
        logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alpsat, cfac, poff, roff , ap = chn[i,:13]
        val  = b0 + b1 * (mh-12) + b2 * (mh-12)**2
        mat[i,:] = val
        print(i)
 
    ymed  = np.percentile(mat,50,axis=0)
    ymin0 = np.percentile(mat,16,axis=0)
    ymax0 = np.percentile(mat,84,axis=0)
    ymin1 = np.percentile(mat,2.5,axis=0)
    ymax1 = np.percentile(mat,97.5,axis=0)

    ax =  plt.subplot(2,2,1)

    xx,yy = np.loadtxt('phis_yang_2009_fig6.csv',unpack=1)
    ax.plot(xx,yy,'.',color='C1',label='yang2009', zorder=10)
    ax.plot(mh, ymed, color='C0', ls='--', label='This work') # 95 percentile
    ax.fill_between(mh, ymin1, ymax1, color='#52aae7', alpha = 0.45, zorder=0) #     95 percentile
    ax.fill_between(mh, ymin0, ymax0, color='#1f77b4', alpha = 0.45, zorder=1) #     68 percentile
    ax.plot(mh, 0.18 +0.83*(mh-13), color='C3', ls='--', label='Uitert 2016') # 95 percentile
    ax.legend()
    ax.set_xlabel(r'$\log M_{200m}$')
    ax.set_ylabel(r'$\log \phi_{\rm s}^*$')
    plt.savefig('/mnt/home/student/cdivya/github/weaklens_pipeline/preet/sn_preet_gama_output_cosmo_zu15_full_revised_fits/%s_phis_csmf_sat.png'%chnpath.split('/')[-1], dpi=250)  
      

    #from astropy.io import fits
    #full_gal = fits.open('../DataStore/preetish/revised_gama_specz_sample_zleq_03.fits')[1].data
    #bins =  full_gal['bin']
    #rsftarr     = np.zeros(7)
    #logMaarr    = np.zeros(7)
    #logMbarr    = np.zeros(7)
    #data = 0.0*rsftarr
    ##joint fit first collecting the measurements
    #for i,bb in enumerate(range(0,7)):
    #    idx = (bins==bb)
    #    rsftarr[i]     = np.median(full_gal['z_tonry'][idx])
    #    logMaarr[i]    = full_gal['logMstar_h2'][idx].min()
    #    logMbarr[i]    = full_gal['logMstar_h2'][idx].max()
    #print rsftarr
    #
    #nofm = {}
    #xh = 0

    #from scipy.interpolate import interp1d
    #xx, yy = np.loadtxt('gama_phistel.dat', unpack=1)
    #gama = interp1d(xx,yy, kind='cubic')
    #
    #try :
    #    for kk in range(len(rsftarr)):
    #        plt.subplot(3,3,kk+1)

    #        logMstelarr, cymin0, cymax0, cymin1, cymax1 = np.loadtxt('%s_zeq_%2.2f_phimstel_cen.dat'%(chnpath,rsftarr[kk]), unpack=1)
    #        logMstelarr, symin0, symax0, symin1, symax1 = np.loadtxt('%s_zeq_%2.2f_phimstel_sat.dat'%(chnpath,rsftarr[kk]), unpack=1)
    #        logMstelarr, tymin0, tymax0, tymin1, tymax1 = np.loadtxt('%s_zeq_%2.2f_phimstel_tot.dat'%(chnpath,rsftarr[kk]), unpack=1)

    #       #p1 = plt.fill_between(logMstelarr, cymin0, cymax0, color='C0', alpha = 0.45, zorder=1) 
    #       #p2 = plt.fill_between(logMstelarr, cymin1, cymax1, color='C0', alpha = 0.25, zorder=2) 
    #       #p3 = plt.fill_between(logMstelarr, symin0, symax0, color='C1', alpha = 0.45, zorder=3) 
    #       #p4 = plt.fill_between(logMstelarr, symin1, symax1, color='C1', alpha = 0.25, zorder=4) 
    #        p5 = plt.fill_between(logMstelarr, tymin0, tymax0, color='C2', alpha = 0.45, zorder=5) 
    #        p6 = plt.fill_between(logMstelarr, tymin1, tymax1, color='C2', alpha = 0.25, zorder=6) 


    #        plt.legend([(p5,p6)], ["tot"], loc='lower left')
    #        #plt.legend([(p1, p2), (p3, p4), (p5,p6)], ["cen", "sat", "tot"], loc='lower left')
    #        plt.title('$z=%2.2f$'%rsftarr[kk])
    #        idx = ((logMstelarr - 2*np.log10(hval))> xx.min()) & ((logMstelarr - 2*np.log10(hval))< xx.max())
    #        plt.plot(logMstelarr[idx], 10**gama(logMstelarr[idx] - 2*np.log10(hval)), '-', c='black')

    #        plt.yscale('log')
    #        plt.ylim(10**-4,0.2)
    #        plt.xlabel(r'$\log[M_{*}/(h^{-2}M_\odot)]$')
    #        plt.ylabel(r'$\phi(\log M_{*})$')

    #except:
    #    chn   = pd.read_csv('%s'%chnpath,  delim_whitespace=True, header= None)
    #    pred  = pd.read_csv('%s'%predpath,  delim_whitespace=True, header= None)
    #    np.random.seed(1123)
    #    idx = np.random.choice(np.arange(len(chn.values[:,0])), size=100, replace=0)
    #    #idx = (np.random.uniform(size=len(chn.values[:,0]))<0.01) 
    #    chn = chn.values[idx,:]
    #    #pred = pd.read_csv('%s'%predpath, delim_whitespace=True, header= None)
    #    logMstelarr = np.linspace(8,12,20) 
    #    cmat = np.zeros((len(chn[:,0]), len(logMstelarr), len(rsftarr)))
    #    smat = np.zeros((len(chn[:,0]), len(logMstelarr), len(rsftarr)))
    #    tmat = np.zeros((len(chn[:,0]), len(logMstelarr), len(rsftarr)))
    #    
    #    for ii in range(len(chn[:,0])):
    #        pp = chn[ii,:13]
    #        for jj in range(len(rsftarr)):
    #            yy    = phi_mstel(pp, logMstelarr, zred=rsftarr[jj])
    #            cmat[ii, :, jj] = yy[0] 
    #            smat[ii, :, jj] = yy[1] 
    #            tmat[ii, :, jj] = yy[0] + yy[1]  
    #        print ii

    #    for kk in range(len(rsftarr)):
    #        plt.subplot(3,3,kk+1)
    #        cymin0 = np.percentile(cmat[:,:,kk], 16,   axis=0)
    #        cymax0 = np.percentile(cmat[:,:,kk], 84,   axis=0)
    #        cymin1 = np.percentile(cmat[:,:,kk], 2.5,  axis=0)
    #        cymax1 = np.percentile(cmat[:,:,kk], 97.5, axis=0)

    #        symin0 = np.percentile(smat[:,:,kk], 16,   axis=0)
    #        symax0 = np.percentile(smat[:,:,kk], 84,   axis=0)
    #        symin1 = np.percentile(smat[:,:,kk], 2.5,  axis=0)
    #        symax1 = np.percentile(smat[:,:,kk], 97.5, axis=0)

    #        tymin0 = np.percentile(tmat[:,:,kk], 16,   axis=0)
    #        tymax0 = np.percentile(tmat[:,:,kk], 84,   axis=0)
    #        tymin1 = np.percentile(tmat[:,:,kk], 2.5,  axis=0)
    #        tymax1 = np.percentile(tmat[:,:,kk], 97.5, axis=0)


    #        dat = np.transpose([logMstelarr, cymin0, cymax0, cymin1, cymax1])
    #        np.savetxt('%s_zeq_%2.2f_phimstel_cen.dat'%(chnpath,rsftarr[kk]), dat, header='logmh(h-1Msun)\t68logMcen_m(h-2Msun)\t68logMcen_p(h-2Msun)\t95logMcen_m(h-2Msun)\t95logMcen_p(h-2Msun)')

    #        dat = np.transpose([logMstelarr, symin0, symax0, symin1, symax1])
    #        np.savetxt('%s_zeq_%2.2f_phimstel_sat.dat'%(chnpath,rsftarr[kk]), dat, header='logmh(h-1Msun)\t68logMcen_m(h-2Msun)\t68logMcen_p(h-2Msun)\t95logMcen_m(h-2Msun)\t95logMcen_p(h-2Msun)')

    #        dat = np.transpose([logMstelarr, tymin0, tymax0, tymin1, tymax1])
    #        np.savetxt('%s_zeq_%2.2f_phimstel_tot.dat'%(chnpath,rsftarr[kk]), dat, header='logmh(h-1Msun)\t68logMcen_m(h-2Msun)\t68logMcen_p(h-2Msun)\t95logMcen_m(h-2Msun)\t95logMcen_p(h-2Msun)')


    #        p1 = plt.fill_between(logMstelarr, cymin0, cymax0, color='C0', alpha = 0.45, zorder=1) 
    #        p2 = plt.fill_between(logMstelarr, cymin1, cymax1, color='C0', alpha = 0.25, zorder=2) 
    #        p3 = plt.fill_between(logMstelarr, symin0, symax0, color='C1', alpha = 0.45, zorder=3) 
    #        p4 = plt.fill_between(logMstelarr, symin1, symax1, color='C1', alpha = 0.25, zorder=4) 
    #        p5 = plt.fill_between(logMstelarr, tymin0, tymax0, color='C2', alpha = 0.45, zorder=5) 
    #        p6 = plt.fill_between(logMstelarr, tymin1, tymax1, color='C2', alpha = 0.25, zorder=6) 


    #        plt.legend([(p1, p2), (p3, p4), (p5,p6)], ["cen", "sat", "tot"], loc='lower right')
    #        plt.yscale('log')
    #        plt.ylim(10**-4,0.2)
    #        plt.xlabel(r'$\log[M_{*}]$')
    #        plt.ylabel(r'$\phi(M_{*})$')


    #plt.tight_layout()        
    #plt.savefig('%s_phistar.png'%chnpath, dpi=300)

    #plt.title('%2.2f,%2.2f'%(logMaarr[zz],logMbarr[zz]))
    #plt.plot(logMstelarr, yy[0], label='cen')
    #plt.plot(logMstelarr, yy[1], label='sat')

    #plt.subplot(2,2,1)
    #for zz in range(len(rsftarr)):
    #    plt.subplot(3,3,zz+1)
    #    zred = rsftarr[zz]
    #    yy =  phi_mstel(bestp, logMstelarr, rsftarr[zz])
    #    plt.title('%2.2f,%2.2f'%(logMaarr[zz],logMbarr[zz]))
    #    plt.plot(logMstelarr, yy[0], label='cen')
    #    plt.plot(logMstelarr, yy[1], label='sat')
    #    plt.yscale('log')
    #
    #    plt.xlabel(r'$\log[M_{*}/(h^{-2}M_\odot)]$')
    #    plt.ylabel(r'$\phi(\log M_{*})$')
    #    plt.legend()
    #plt.tight_layout()    



    ##logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alpsat, cfac, poff, roff , ap = x
    #x = [9.95, 11.24, 3.18, 0.245, 0.157, -1.17, 1.53, -0.217, -1.18, 1.0, 1.0, 0.1, 1.0]

    #logMstelarr = np.linspace(8,12,20)
    #yy    = phi_mstel(x, logMstelarr, zred=0.1)

    #yylogc = np.log(10) * 10**logMstelarr * yy[0]
    #yylogs = np.log(10) * 10**logMstelarr * yy[1]

    #plt.subplot(2,2,1)
    ##plt.plot(-2.5*logMstelarr+4.67, np.log10(yy[0]))
    ##plt.plot(-2.5*logMstelarr+4.67, np.log10(yy[1]))
    #plt.plot(logMstelarr , np.log10(yy[0] + yy[1]), label='total')
    #plt.plot(logMstelarr , np.log10(yy[0]), ls ='--', label='cen')
    #plt.plot(logMstelarr , np.log10(yy[1]), ls='--', label='sat')

    ##xx,yy = np.loadtxt('clf_phimstel_mrac2013.dat', unpack=1)
    #plt.xlabel(r'$\log[M^*]$')
    #plt.ylabel(r'$\log[\phi(M^*)]$')
    ##plt.plot(xx,yy, '.')
    ##plt.xlim(-17,-24)
    ##plt.ylim(-20, -10)
    #plt.legend()
    ##plt.yscale('log')
    #plt.savefig('test.png', dpi=300)
    #sys.exit()
 


