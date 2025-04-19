from mcmc import model
import matplotlib.pyplot as plt
import numpy as np
#print model([12.94356702,  0.30769993,  0.76663147,  0.38569895,  0.07206964,0.98933412],np.linspace(0.01,1.0,10),0.25,11.3)
plt.subplot(2,2,1)
#logMmin, siglogM, alphasat, poff, roff, cfac = x
rbin = np.linspace(0.01,1.0,10)
for cc in np.linspace(1,10,5):
    x = [12.94356702,  0.30769993,  0.76663147,  0.38569895,  0.07206964, cc]

    plt.plot(rbin,model(x,rbin,0.1,9)[:len(rbin)], label = r'$c_{\rm fac} = %2.2f$'%cc)

plt.title(r'$\log_{\rm Mmin}=%2.2f, \sigma_{\log M}=%2.2f, \alpha_{\rm sat}=%2.2f, {\rm p_{off}}=%2.2f, {\rm r_{off}}=%2.2f, z_{\rm red}=%s, \log M^*=%s $'%(x[0], x[1], x[2], x[3], x[4], 0.1, 9))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$R [{\rm h^{-1} Mpc}] $')
plt.ylabel(r'$\Delta\Sigma [{\rm h M_\odot pc^{-2}}]$')
plt.legend()
plt.savefig('test.png',dpi=300)

