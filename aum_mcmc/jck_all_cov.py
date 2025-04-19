import numpy as np
import pandas as pd
import csv
import sys
import matplotlib.pyplot as plt
from glob import glob

lenspath = sys.argv[1]
randpath = sys.argv[2]
flist = glob('%s/dsigma.dat_*'%lenspath)

nbins = len(flist)-1
print(nbins)
df={}
rdf={}

for i in range(nbins):
    df['%d'%i]  = np.loadtxt('%s/dsigma.dat_%d'%(lenspath,i))
    rdf['%d'%i] = np.loadtxt('%s/dsigma.dat_%d'%(randpath,i))


njacks = len(np.unique(df['0'][:,-1]))

print(njacks)
#input file reading
df_ds = {};df_r12 = {}
rdf_ds = {};rdf_r12 = {}

for ii in range(njacks):
    for jj in range(nbins):
        idx = (df['%d'%jj][:,-1] == ii)   
        ds = df['%d'%jj][idx,5]
        r12 = df['%d'%jj][idx,14]
        ds = ds/(1.0+r12)

        ridx = (rdf['%d'%jj][:,-1] == ii) 
        rds = rdf['%d'%jj][ridx,5]
        rr12 = rdf['%d'%jj][ridx,14]
        rds = rds/(1.0+rr12)
        
        print(ii, jj, np.shape(ds))
        ds = ds - rds #no boosted signal
        ds = ds[:10] # taking only first 10 bins
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

        cov[i1][i2] = (njacks - 1)*s1/(1.0*njacks) #covariance for jackknifes
        
print('writing to the file')
cross_cof_matrix = np.matrix(cov)#cov

fname ='%s/cov_w_rand_subs_no_boost_removed_last_bin.dat'%lenspath
np.savetxt(fname,cross_cof_matrix)

