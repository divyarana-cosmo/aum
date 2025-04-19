import numpy as np
import pandas as pd
import csv
import sys
import matplotlib.pyplot as plt
from glob import glob

#ONLY TAKING FIRST 10 RADIAL BINS

lenspath = './sn_preet_gama_output_cosmo_zu15_full_revised'
#lenspath = './sn_preet_gama_output_cosmo_y08_full_revised'
njacks = 500 
nbins = 7
#input file reading
df_ds = {}
for ii in range(njacks):
    for jj,bb in enumerate(range(0,nbins)):
        df  = np.loadtxt('%s/dsigma.dat%05d_%d'%(lenspath,ii,bb))
        #if bb<=3:
        #    ds  = df[:7,5]
        #    r12 = df[:7,14]
        #else:
        #    ds  = df[:10,5]
        #    r12 = df[:10,14]
        #ds  = df[:10,5] # tangential
        ds  = df[:10,8]  # cross
        r12 = df[:10,14]


        ds = ds/(1.0+r12)

        if jj==0:
            df_ds['%d'%ii] = ds
        else:
            df_ds['%d'%ii] = np.concatenate((df_ds['%d'%ii],ds))
        

print('reading done')
size = len(df_ds['2'])

mean_m = np.zeros(size)
#segment for computing mean
for i1 in range(size):
    sm = 0.0
    for i2 in range(njacks):
        sm += df_ds['%d'%i2][i1]
    mean_m[i1] =  sm/(1.0*njacks)
    
print('mean computed')
#segment for computing convariance matrix

cov = np.zeros((size,size))
for i1 in range(size):
    for i2 in range(size):
        s1=0.0
        for i3 in range(njacks):
            r0 =  df_ds['%d'%i3][i1] - mean_m[i1]
            r1 =  df_ds['%d'%i3][i2] - mean_m[i2] 
            s1 += r0 * r1

        cov[i1][i2] = s1/(1.0*njacks) #covariance for jackknifes
        
print('writing to the file')
cross_cof_matrix = np.matrix(cov)#cov

fname ='%s/cov_no_boost_x.dat'%lenspath
#fname ='%s/cov_no_boost_7_rbins_for_first_4_bins.dat'%lenspath
np.savetxt(fname,cross_cof_matrix)

