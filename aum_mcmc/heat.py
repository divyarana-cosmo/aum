import numpy as np
import matplotlib.pyplot as plt
import sys

def tlab(i,j):
    a = r'$\Delta\Sigma_{%d,%d}$'%(i,j)
    return a

#print ticklab(0,1)
plt.subplot(2,2,1)
#plt.figure(figsize=[6.0,6.0])
fname = sys.argv[1] 

cov  = np.loadtxt('%s'%fname)
#cov  = np.loadtxt('./gama_v10_lumb/%s.dat'%fname)

size = np.size(cov[:,0])

cross_coeff = 0.0*cov
for i1 in range(size):
    for i2 in range(size):
        cross_coeff[i1][i2] = cov[i1][i2]/np.sqrt(cov[i1][i1]*cov[i2][i2])

#print cross_coeff
#plt.imshow(cov,cmap='PuOr_r',origin='lower')
plt.imshow(cross_coeff,cmap='PuOr_r',vmin=-1,vmax=1,origin='lower',aspect='equal')
#plt.imshow(cov,cmap='PuOr_r',origin='lower')
#rbin=np.logspace(np.log10(0.02),np.log10(2),10)
#plt.pcolormesh(rbin,rbin,cross_coeff,cmap='PuOr_r',vmin=-1,vmax=1)
#plt.xlabel(r'$\Delta\Sigma_{i,j}$')
#plt.ylabel(r'$\Delta\Sigma_{i,j}$')

#plt.xscale('log')
#plt.yscale('log')
#for i in range(1,13):
for i in range(1,8):
    if i==1:
        lticks =[tlab(i,1),tlab(i,6)]
    else:
        lticks = np.concatenate((lticks,[tlab(i,1),tlab(i,6)]))

cb = plt.colorbar(fraction=0.046,pad=0.04)
#cb = plt.colorbar(label=r"$r_{\rm ij}$",fraction=0.046,pad=0.04)
cb.set_label(label=r"$r_{\rm ij}$",size=15)
#cb.set_label(label=r"$r^{\rm JK}_{\rm ij}$",size=15)
cb.ax.tick_params(direction='in',length=2.5,labelsize=10)
plt.xticks(np.arange(0,70,5),lticks,fontsize=8,rotation=90)
plt.yticks(np.arange(0,70,5),lticks,fontsize=8)


#plt.xticks(np.arange(0,120,5),lticks,fontsize=8,rotation=90)
#plt.yticks(np.arange(0,120,5),lticks,fontsize=8)
plt.tick_params(axis='both',length=2.5,direction='in')

plt.savefig('%s.pdf'%fname, dpi=300)
#plt.savefig('%s.png'%fname, dpi=300)
#plt.savefig('./re_full_jk_output_1_cosmo_y08/cov_w_rand_subs_no_boost_nfw.png', dpi=300)
#plt.savefig('./re_full_jk_output_1_cosmo_y08/cov_w_rand_subs_no_boost.png', dpi=300)
#plt.savefig('./re_full_jk_output_1/cov_w_rand_subs_no_boost.png', dpi=300)
#plt.savefig('./jk_output_1/cov_w_rand_subs.png', dpi=300)

