impo
import numpy as np
import matplotlib.pyplot as plt
import corner
import sys
from astropy.io import fits
samp = fits.open('../DataStore/preetish/data_table.fits')[1].data;

from behroz import Behroozi_2013a

bh = Behroozi_2013a()

def get_yvals(x):
    y = np.median(x)
    ymin = y - np.percentile(x,16),
    ymax = np.percentile(x,84) - y
    return y, ymin, ymax


xbh = np.zeros(12)
ybh = 0.0*xbh
def plotc(n):
    plt.rc('font', size=14)
    darkdata = pd.read_csv('./full_jk_fits_1/new_x9_16_off_chainfile_%s.dat'%n,delim_whitespace=True,header=None)
    darkdata1 = pd.read_csv('./full_jk_fits_1/new_x9_16_off_predfile_%s.dat'%n,delim_whitespace=True,header=None)
    
       
    print('data reading completed')

    darkdat = 0.0 * darkdata1.values[:,:9]

    darkdat[:,:7] = darkdata.values[:,:7]
    darkdat[:,7] = darkdata1.values[:,-3]
    darkdat[:,8] = darkdata1.values[:,-2]

      
    #logMmin, siglogM, alphasat, poff, roff, cfac, ap

    #dnames = [r"$\log M_{\rm min}$", r"$a_{\rm p}$", r"$\log[\langle M \rangle]$", r'$f_{\rm sat}$']
    
    mat = 0.0*darkdat[:,:4]
    mat[:,0] = darkdat[:,0]
    mat[:,1:] = darkdat[:,-3:]

    mcen = np.loadtxt('./full_jk_fits_1/new_avg_mhcen_%d.dat'%n)
 
    if n>=8:
<<<<<<< HEAD
       print("cutting chains")
=======
       print "cutting chains"
>>>>>>> 2cb08947c4825a1e0e38de69838ff8bfaff35728
       chisq_cut = np.percentile(darkdata1.values[:,-1],70)
       chiidx  = (darkdata1.values[:,-1] < chisq_cut)
       mat = mat[chiidx,:]
       mcen = mcen[chiidx]

    mcen = mcen - np.log10(0.67386) #removing h-1

    data = fits.open("../DataStore/preetish/data_table.fits")[1].data
    sel = (data['bin']==n) #& (z>0) & (dec<50) #also removing aegis

    full_gal = fits.open('../DataStore/preetish/pdr2_zleq_03_frm_16A_overlap_blend_mass_cut.fits')[1].data

    fidx = (full_gal['logM']> data['logM'][sel].min()) & (full_gal['logM']<=data['logM'][sel].max()) & (full_gal['dec']<50) & (full_gal['photoz_best']>0)

    rsft = np.mean(full_gal['photoz_best'][fidx])
    logMstel = full_gal['logM'][fidx]
    #logMstel = np.log10(10**logMstel)

    x = [np.median(logMstel)]
    xerr = [np.median(logMstel)-np.percentile(logMstel,16) , np.percentile(logMstel,84) - np.median(logMstel)]

    ax = plt.subplot(2,2,1)
    y = [np.median(mat[:,0])]
    yerr = [np.median(mat[:,0])-np.percentile(mat[:,0],16),np.percentile(mat[:,0],84) - np.median(mat[:,0])]
    ax.errorbar(np.array(x), np.array(y), xerr=np.array(xerr)[:,None], yerr=np.array(yerr)[:,None], color='blue', capsize=3)
    if n==1:
        ax.set_xlabel(r'$\log[M_*/ ( M_{\odot})]$')
        ax.set_ylabel(r'$\log[M_{\rm min} /(h^{-1} M_{\odot})]$')
 
    ax = plt.subplot(2,2,2)
    y = np.median(mcen)

    #ybh[n-1] = bh.SHMRbeh(np.array(y),rsft)
    #xbh[n-1] =  np.array(y)       
    #if n==12:
    #    ax.legend(fontsize='small')

    yerr = [y-np.percentile(mcen,16),np.percentile(mcen,84) - y]
 
    ax.errorbar(np.array(y), np.array(x), yerr=np.array(xerr)[:,None], xerr=np.array(yerr)[:,None], color='blue', capsize=3)
    if n==1:
        ax.set_ylabel(r'$\log[M_*/(M_{\odot})]$')
        ax.set_xlabel(r"$\log[\langle M_{\rm cen} \rangle / (M_{\odot})]$")
        #ax.set_xlabel(r"$\log[\langle M_{\rm cen} \rangle / (M_{\odot})]$")


    ax = plt.subplot(2,2,3)
    y = np.median(mat[:,-1])
    yerr = [y-np.percentile(mat[:,-1],16),np.percentile(mat[:,-1],84) - y]
    ax.errorbar(np.array(x), np.array(y), xerr=np.array(xerr)[:,None], yerr=np.array(yerr)[:,None], color='blue', capsize=3)
    if n==1:
        ax.set_xlabel(r'$\log[M_*/ (M_{\odot})]$')
        ax.set_ylabel(r"$f_{\rm sat}$")
 
        
    
for i in range(1,13):
    plotc(i)
<<<<<<< HEAD
    print(i)
=======
    print i
>>>>>>> 2cb08947c4825a1e0e38de69838ff8bfaff35728

plt.tight_layout()
plt.savefig('full_jk_fits_1/new_mstel_rels.png', dpi=300)












 
