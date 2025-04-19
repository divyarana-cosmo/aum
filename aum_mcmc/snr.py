#!/usr/bin/env python3
import sys
sys.path.append("/mnt/csoft/tools/anaconda2/bin/python")
import numpy as np
import pandas as pd
import time
import os
#import cosmology as cc
from scipy.integrate import quad

from math import log10


bins, mstel_min, mstel_max, mstel_avg, z_mean, z_med = np.loadtxt("../DataStore/preetish/mstel_bins_gama.dat",     unpack=True)


#joint fit first collecting the measurements
for i,bb in enumerate(range(5,12)):
    idx = (bins==bb)
    dat  = np.loadtxt('./jk_preet_gama_output_cosmo_y08/dsigma.dat_%d'%bb)
    rdat = np.loadtxt('./jk_preet_gama_rand_output_cosmo_y08/dsigma.dat_%d'%bb)

    rbin = np.unique(dat[:,7])
    y = 0.0*rbin
    if i==0:
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

#then get the full covariance
cov = np.loadtxt('./jk_preet_gama_output_cosmo_y08/cov_w_rand_subs_no_boost.dat',dtype=float,unpack=True)
cov = cov[40:, 40:] #removing first 4 bins
icov = np.linalg.inv(cov)

njacks = int(len(dat[:,7])/len(rbin))
hartlap_factor = 1.0#(njacks - len(data) - 2) * 1.0/(njacks - 1)

icov = hartlap_factor*icov


Delta = data
chisq = np.dot(Delta, np.dot(icov, Delta))
<<<<<<< HEAD
print(np.sqrt(chisq))
=======
print np.sqrt(chisq)
>>>>>>> 2cb08947c4825a1e0e38de69838ff8bfaff35728

