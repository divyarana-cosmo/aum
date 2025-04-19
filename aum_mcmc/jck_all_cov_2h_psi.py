#removing radial bins from 0.2-0.8 h-1 Mpc
import numpy as np
import pandas as pd
import csv
import sys
import matplotlib.pyplot as plt

df={}
rdf={}
for i in range(3,13):
    df['%d'%i]  = np.loadtxt('./re_full_jk_output_1_cosmo_y08/dsigma.dat_%d'%i)
    rdf['%d'%i] = np.loadtxt('./re_full_jk_random_1_cosmo_y08/dsigma.dat_%d'%i)
    #df['%d'%i]  = np.loadtxt('./re_full_jk_output_1/dsigma.dat_%d'%i)
    #rdf['%d'%i] = np.loadtxt('./re_full_jk_random_1/dsigma.dat_%d'%i)


njacks = len(np.unique(df['5'][:,-1]))
print(njacks)
#input file readingi
df_ds = {};df_r12 = {}
rdf_ds = {};rdf_r12 = {}

for ii in range(njacks):
    for jj in range(3,13):
        idx = (df['%d'%jj][:,-1] == ii) & ~((df['%d'%jj][:,7]>0.2) & (df['%d'%jj][:,7]<0.8))
        
        #idx = (df['%d'%jj][:,-1] == ii)
        ds = df['%d'%jj][idx,5]
        r12 = df['%d'%jj][idx,14]
        avg_wls = df['%d'%jj][idx,6]
        ds = ds/(1.0+r12)

        ridx = (rdf['%d'%jj][:,-1] == ii) & ~((rdf['%d'%jj][:,7]>0.2) & (rdf['%d'%jj][:,7]<0.8))
        #ridx = (rdf['%d'%jj][:,-1] == ii)
        rds = rdf['%d'%jj][ridx,5]
        rr12 = rdf['%d'%jj][ridx,14]
        ravg_wls = rdf['%d'%jj][ridx,6]
        rds = rds/(1.0+rr12)

        ds = ds - rds #no boosted signal
        #ds = ds*1.0*avg_wls/ravg_wls - rds #boosted signal
        if jj==3:
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

#fname ='./re_full_jk_output_1/cov_w_rand_subs_no_boost.dat'
fname ='./re_full_jk_output_1_cosmo_y08/cov_w_rand_subs_no_boost.dat_2h_psi'
np.savetxt(fname,cross_cof_matrix)

#with open(fname,'wb') as f:
#    for line in cross_cof_matrix:
#        np.savetxt(f, line, fmt='%2.5f')
