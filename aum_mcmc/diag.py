import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

def pltshpnos(n):
    m = n+1
    ax = plt.subplot(3,3,m)

    for i in range(320):
        filename = '~/github/weaklens_pipeline/gama/gama_v10_lumb/output/%s_%s/dsigma.dat%05d' % (lum[n],lum[m],i)
        dfdict = pd.read_csv(filename, delim_whitespace = True, usecols=([5,6,7,8,11,13,14]), header=None, names=(["dsigma","Sumwls_by_sumwl","r","dsigma_cross","sn_err","tot_pairs","r12_sel"]), comment='#')

        ax.plot(dfdict.r.values,dfdict.dsigma.values,'.',lw=0.0)

    ax.text(0.1, 0.1, r"$[%s - %s]$"%(lum[n],lum[m]), transform=ax.transAxes)
    ax.axhline(y=0.0,color='grey',alpha=0.5)
    #if n<=2:
    #    ax.set_ylim(0,1.5)
    #else:
    #    ax.set_ylim(0,1.5)
    if n>2:
        ax.set_xlabel(r'$R[{\rm h^{-1}Mpc}]$')
    
    ax.set_xscale('log')
    #if not (n==0 or n==3):
    #    ax.set_yticklabels([])
    if n==0 or n==3:
        ax.set_ylabel(r'$\Delta\Sigma$')
                             
    if n<=2:
        ax.set_xticklabels([])

def pltrandnos(n):
    for i in range(32):
        filename = './random_1/dsigma.dat%05d_%d' % (i,n)
        dfdict = pd.read_csv(filename, delim_whitespace = True, usecols=([5,6,7,8,11,13,14]), header=None, names=(["dsigma","Sumwls_by_sumwl","r","dsigma_cross","sn_err","tot_pairs","r12_sel"]), comment='#')

        plt.plot(dfdict.r.values,dfdict.dsigma.values*(1.0/(1+dfdict.r12_sel.values)),'.',lw=0.0)

    plt.title( "binno = %d"%(n))
    plt.axhline(y=0.0,color='grey',alpha=0.5)
    #if n<=2:
    #    ax.set_ylim(0,1.5)
    #else:
    #    ax.set_ylim(0,1.5)
    plt.xlabel(r'$R$')
    
    plt.xscale('log')
    plt.ylabel(r'$\Delta\Sigma$')
                             
    plt.show()
    plt.clf()

def pltjk(n):
    ii = int((n-1)/2) + 1
    ax=plt.subplot(3,3, ii)
    #plt.subplot(2,2,)
    #filename = './full_jk_random_1/dsigma.dat_%d' % (n)
    filename = './re_full_jk_random_1_cosmo_y08/dsigma.dat_%d' % (n)
    dfdict = pd.read_csv(filename, delim_whitespace = True, usecols=([5,6,7,8,11,13,14]), header=None, names=(["dsigma","Sumwls_by_sumwl","r","dsigma_cross","sn_err","tot_pairs","r12_sel"]), comment='#')

    xjk = dfdict.r.values
    yjk = dfdict.dsigma.values/(1.0 + dfdict.r12_sel.values)
    y1jk = dfdict.Sumwls_by_sumwl.values

    x = np.unique(xjk)
    njacks = int(len(xjk)/len(x))
    #print njacks
    y = 0.0*x
    yerr = 0.0*y
    y1 = 0.0*x
    y1err = 0.0*y

    for ii1,rr in enumerate(x):
        idx = (xjk==rr)
        y[ii1] = np.mean(yjk[idx])
        yerr[ii1] = np.sqrt(njacks -1) * np.std(yjk[idx])

        y1[ii1] = np.mean(y1jk[idx])
        y1err[ii1] = np.sqrt(njacks -1) * np.std(y1jk[idx])

    ax.errorbar(x, y, yerr, fmt = '.', capsize=3, label= "binno = %d"%(n))

    ax.axhline(y=0.0,color='grey', ls='--', alpha=0.5)
    
    if ii<=3:
        ax.set_ylim(-2,2)
    else:
        ax.set_ylim(-4,4)

    if ii>3:
        ax.set_xlabel(r'$R$')
    ax.set_xscale('log')
    if ii==1 or ii==4:
        ax.set_ylabel(r'$\Delta\Sigma_{\rm rand}$')
    ax.legend()


#plt.figure(figsize=[8.0,8.0])

#for i in range(1,13):
for i in range(3,13):
    pltjk(i)
    #pltrandnos(i)
    #print i
plt.tight_layout()
plt.savefig('./re_full_jk_random_1_cosmo_y08/rand.png', dpi=300)
#plt.savefig('./full_jk_random_1/rand.png', dpi=300)
