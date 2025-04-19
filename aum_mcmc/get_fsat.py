import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

#predpath = 'predfile_bosch_2013_loglin_alpsat_uni_const_alpha_no_off.dat'
predpath = 'predfile_bosch_2013_loglin_alpsat_uni_const_alpha_no_off.dat_wide_gamma1'
pred = pd.read_csv('./sn_preet_gama_output_cosmo_zu15_full_revised_fits/%s'%predpath, delim_whitespace=1,header=None)
#pred = pd.read_csv('/scratch/cdivya/%s'%predpath, delim_whitespace=1,header=None)
idx = (pred.values[:,-1]<300)
#idx = (pred.values[:,77]<0.5)
fsat = pred.values[:, 77:-1]
#fsat = pred.values[idx, 77:-1]

#bins, mstel_min, mstel_max, mstel_avg, z_mean, z_med = np.loadtxt("../DataStore/preetish/mstel_bins_gama.dat", unpack=True)
#
xx    = np.zeros(7)
yy    = np.zeros(7)
yerr0 = 0.0 *xx
yerr1 = 0.0 *xx

from astropy.io import fits
full_gal = fits.open('../DataStore/preetish/revised_gama_specz_sample_zleq_03.fits')[1].data
#full_gal = fits.open('../DataStore/preetish/gama_equatorial_zleq_03_Mrpetro_cen_sat_iso_flags_gr_colour_reduced_sample.fits')[1].data
bins =  full_gal['bin']

for ii,bb in enumerate(range(0,7)):
    sel    = (bins==bb)
    xx[ii] = float(np.log10(np.mean(10**full_gal['logMstar_h2'][sel])) )
    yy[ii]    = np.percentile(fsat[:,ii],50)
    yerr0[ii] = yy[ii] - np.percentile(fsat[:,ii],16)
    yerr1[ii] = np.percentile(fsat[:,ii],84)- yy[ii]

plt.subplot(2,1,1)

plt.errorbar(xx, yy, yerr=[yerr0,yerr1], fmt='.', capsize=3, label='This work')

yangfsat = pd.read_csv('./fsat_yang_2008_tab3.csv', delim_whitespace=1)
plt.errorbar(yangfsat['bin'], yangfsat['all'], yerr=yangfsat['allerr'] , fmt='.', label='yang2008', capsize=3)   

zufsat = pd.read_csv('./fsat_zu_2015_tab1.csv', delim_whitespace=1)
plt.errorbar(np.log10((10**zufsat['bmin'] + 10**zufsat['bmax'])/2.0), zufsat['fsat'], yerr=[zufsat['fsaterr_m'],zufsat['fsaterr_p']] , fmt='.', label='zu2015', capsize=3)   


plt.xlabel(r'$\log[{\rm M_{\rm *} /(h^{-2} M_\odot)}]$')
plt.ylabel(r'$f_{\rm sat}$')
plt.ylim(0.0,1.0)
plt.legend()

plt.savefig('/mnt/home/student/cdivya/github/weaklens_pipeline/preet/sn_preet_gama_output_cosmo_zu15_full_revised_fits/%s_fsat.png'%predpath, dpi=200)
#plt.savefig('%s_fsat_leq_0.5.png'%predpath, dpi=200)
