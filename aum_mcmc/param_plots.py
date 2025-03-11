import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import corner
import sys
from astropy.io import fits
plt.figure(figsize=[10,8])

chnfil = sys.argv[1]
predfil = sys.argv[2]
def get_scat():
    darkdata = pd.read_csv('%s'%chnfil, delim_whitespace=True, header=None)
    darkdata1 = pd.read_csv('%s'%predfil, delim_whitespace=True, header=None)
    print('data reading completed')
    darkdat = 0.0 * darkdata1.values[:,:20]
    darkdat[:,:13] = darkdata.values[:,:13]
    darkdat[:,13:20] = darkdata1.values[:,77:84]
    
    labels = [r"$\log M_0$", r"$\log M_1$", r"$\gamma_1$", r"$\gamma_2$", r"$\sigma_0$", r"$b_0$", r"$b_1$", r"$b_2$", r"$\alpha_{\rm sat}$", r"$c_{\rm fac}$", r"$p_{\rm off}$", r"$r_{\rm off}$", r"$a_{\rm p}$", r"$f^1_{sat}$", r"$f^2_{sat}$", r"$f^3_{sat}$", r"$f^4_{sat}$", r"$f^5_{sat}$", r"$f^6_{sat}$", r"$f^7_{sat}$"]    

    for i in range(20):
        ax = plt.subplot(4,5,i+1)
        ax.plot(darkdat[:,i] , darkdata1.values[:,-1],'.',ms=0.1,lw=0.0, label = labels[i])
        ax.legend()
        if i%5==0:
            ax.set_ylabel(r'$\chi^{2}$')
        if i%5!=0:
            ax.set_yticklabels([])
        ax.set_ylim(30,60)        
    plt.tight_layout()            
    plt.savefig('%s_param_scat.png'%chnfil,dpi=300)


get_scat()

#data = fits.open("../DataStore/preetish/data_table.fits")[1].data
#full_gal = fits.open('../DataStore/preetish/pdr2_zleq_03_frm_16A_overlap_blend_mass_cut.fits')[1].data
#
#
#def plotc(n):
#    #plt.rc('font', size=18)
#    darkdata = pd.read_csv('./full_jk_fits_1/new_x9_16_off_chainfile_%s.dat'%n,delim_whitespace=True,header=None)
#    darkdata1 = pd.read_csv('./full_jk_fits_1/new_x9_16_off_predfile_%s.dat'%n,delim_whitespace=True,header=None)
#    mcen = np.loadtxt('./full_jk_fits_1/new_avg_mhcen_%s.dat'%n) 
#    print('data reading completed')
#
#    chisq_cut = np.percentile(darkdata1.values[:,-1],70)
#    
#    idx  = (darkdata1.values[:,-1] < chisq_cut)
#    darkdat = 0.0 * darkdata1.values[:,:9]
#
#    darkdat[:,:6] = darkdata.values[:,:6]
#    darkdat[:,6] = mcen[:]
#    darkdat[:,7] = darkdata1.values[:,-3]
#    darkdat[:,8] = darkdata1.values[:,-2]
#    chisq =  darkdata1.values[:,-1]
#
#
#    if n>=8:
#        print 'cutting chains'
#        darkdat = darkdat[idx,:]
#        chisq = chisq[idx]
#    #logMmin, siglogM, alphasat, poff, roff, cfac, ap
#
#    #dnames = [r"$\log M_{\rm min}$", r"$a_{\rm p}$", r"$\log[\langle M \rangle]$", r'$f_{\rm sat}$']
#    dnames = [r"$\log M_{min}$", r"$\sigma_{\log M}$", r'$\alpha_{\rm sat}$', r"$f_{\rm off}$", r"$r_{\rm off}$", r"$c_{\rm fac}$", r"$\log[\langle M_{\rm cen} \rangle] $", r"$\log[\langle M \rangle]$", r'$f_{\rm sat}$']
#    #dnames = [r"$\log M_{min}$", r"$\sigma_{\log M}$", r'$\alpha_{\rm sat}$', r"$f_{\rm off}$", r"$r_{\rm off}$", r"$c_{\rm fac}$", r"$a_{\rm p}$", r"$\log[\langle M \rangle]$", r'$f_{\rm sat}$']
#    
#    #mat = 0.0*darkdat[:,:4]
#    #mat[:,0] = darkdat[:,0]
#    #mat[:,1:] = darkdat[:,-3:]
#    #corner.corner(mat, labels=dnames, show_titles=True)
#    
#    bstpar = np.argmin(chisq)
#    
#    corner_kwargs = {}
#    corner_kwargs['plot_datapoints'] = False
#    corner_kwargs['plot_density'] = False
#    corner_kwargs['smooth'] = True
#    corner_kwargs['fill_contours'] = True
#    corner_kwargs['levels'] = [0.68,0.95]
#
#    
#    corner.corner(darkdat, labels=dnames, color='#1f77b4',  show_titles=True, truths=darkdat[bstpar,:], truth_color='red', **corner_kwargs)
#
#    sel = (data['bin']==n) #& (z>0) & (dec<50) #also removing aegis
#   
#    idx = (full_gal['logM']> data['logM'][sel].min()) & (full_gal['logM']<= data['logM'][sel].max()) & (full_gal['dec']<50) & (full_gal['photoz_best']>0)
#
#    rsft = np.mean(full_gal['photoz_best'][idx])
#    logMstel = np.mean(full_gal['logM'][idx])
#
#
#    plt.suptitle(r'Binno=%d, $\langle\log[M_*/(M_{\odot})]\rangle = %2.2f$, $z_{\rm bin} = %2.2f$, Halo-Masses in units of $h^{-1} M_{\odot}$'%(n, logMstel,rsft ),ha='left')
#    plt.savefig('./full_jk_fits_1/new_param_cor_%s.png'%n,dpi=400)
#    plt.clf()
#for i in range(1,13):
#    plt.rcParams.update({'font.size': 20})
#    plt.rcParams.update({'xtick.labelsize': 16})
#    plt.rcParams.update({'ytick.labelsize': 16})
#    plotc(i)
#    #get_scat(i)
#    print i













 
