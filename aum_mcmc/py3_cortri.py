import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import corner

import sys
import pygtc

chnfil  = sys.argv[1]
predfil = sys.argv[2]
chn  = pd.read_csv('%s'%chnfil, delim_whitespace=True,header=None)
nrows = len(chn.values[:,0])
pred = pd.read_csv('%s'%predfil,delim_whitespace=True,header=None)
#pred = pd.read_csv('%s'%predfil, usecols=[84], delim_whitespace=True,header=None)
idx = np.argmin(pred.values[:,-1])
bestp = chn.values[idx,:14]

idx = (pred.values[:,-1] < 300)
chn = chn.values[idx,:]

dfname = '%s_corner.png'%chnfil

#logM0, logM1, gamma1, gamma2, sig0, b0, b1, b2, alpsat, cfac, poff, roff , ap = x
def plotc():

    dnames = [r"$\log M_0$", r"$\log M_1$", r"$\gamma_1$", r"$\gamma_2$", r"$\sigma_0$", r"$b_0$", r"$b_1$", r"$b_2$", r"$\alpha_{\rm sat}$", r"$\eta$", r"$c_{\rm fac}$", r"$p_{\rm off}$", r"$r_{\rm off}$", r"$a_{\rm p}$"]
    #dnames = [r"$\log M_0$", r"$\log M_1$", r"$\gamma_1$", r"$\gamma_2$", r"$\sigma_0$", r"$b_0$", r"$b_1$", r"$b_2$", r"$\alpha_{\rm sat}$", r"$c_{\rm fac}$", r"$p_{\rm off}$", r"$r_{\rm off}$", r"$a_{\rm p}$"]
    GTC = pygtc.plotGTC(chains=chn[:,:14], paramNames=dnames, truths=bestp, figureSize='MNRAS_page', panelSpacing = 'loose', customTickFont={'family':'Arial', 'size':5},customLabelFont={'family':'Arial', 'size':6} )
    plt.savefig(dfname, dpi=600)

plotc()




 
