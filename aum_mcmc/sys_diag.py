import matplotlib as mpl
mpl.use('Qt4Agg')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import seaborn as sns

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
    plt.subplot(1,1,1)
    filename = './jk_output_1/dsigma.dat_%d' % (n)
    dfdict = pd.read_csv(filename, delim_whitespace = True, usecols=([5,6,7,8,11,13,14]), header=None, names=(["dsigma","Sumwls_by_sumwl","r","dsigma_cross","sn_err","tot_pairs","r12_sel"]), comment='#')

    xjk = dfdict.r.values
    yjk = dfdict.dsigma.values/(1.0 + dfdict.r12_sel.values)
    y1jk = dfdict.Sumwls_by_sumwl.values


    #randrel = np.loadtxt('./random_1/dsigma.dat_%d' % (n))
    x = np.unique(xjk)
    njacks = int(len(xjk)/len(x))
    print(njacks)
    y = 0.0*x
    yerr = 0.0*y
    y1 = 0.0*x
    y1err = 0.0*y


    for ii,rr in enumerate(x):
        idx = (xjk==rr)
        y = yjk[idx]
        
        sns.violinplot(x=rr, y=y)
    #plt.plot(dfdict.r.values,dfdict.dsigma.values*(1.0/(1+dfdict.r12_sel.values)),'.',lw=0.0)
    
    #plt.plot(randrel[:,0], randrel[:,1],'.', lw=0.0)

    plt.title( "binno = %d"%(n))
    #plt.axhline(y=0.0,color='grey',alpha=0.5)
    #if n<=2:
    #    ax.set_ylim(0,1.5)
    #else:
    #    ax.set_ylim(0,1.5)
    plt.xlabel(r'$R$')
    
    #plt.ylim(0.05,)
    plt.xscale('log')
    #plt.yscale('log')
    plt.ylabel(r'$\Delta\Sigma$')
                            
    plt.show()
    plt.clf()



for i in range(1,13):
    pltjk(i)
    #pltrandnos(i)
    #print i

