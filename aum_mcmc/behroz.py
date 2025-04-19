import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import pandas as pd
class Behroozi_2013a:
    def __init__(self):
        self.eps0=-1.777;
        self.epsa=-0.006;
        self.epsz=-0.000;
        self.epsa2=-0.119;

        self.xM10=11.514;
        self.xM1a=-1.793;
        self.xM1z=-0.251;

        self.alp0=1.412;
        self.alpa=-0.731;

        self.delta0=3.508;
        self.deltaa=2.608;
        self.deltaz=-0.043;

        self.gamma0=0.316;
        self.gammaa=1.319;
        self.gammaz=0.279;


    def evolved_factors(self,zz):
        a=1./(1.+zz);
        nu=np.exp(-4*a*a);
        alp=self.alp0+self.alpa*(a-1.0)*nu;
        delta=self.delta0+(self.deltaa*(a-1.)+self.deltaz*zz)*nu;
        gamma=self.gamma0+(self.gammaa*(a-1.)+self.gammaz*zz)*nu;
        eps=self.eps0+(self.epsa*(a-1.0)+self.epsz*zz)*nu+self.epsa2*(a-1.)
        xM1=self.xM10+(self.xM1a*(a-1.0)+self.xM1z*zz)*nu;
        return nu, alp, delta, gamma, eps, xM1

    def fbeh12(self, x, nu, alp, delta, gamma):
        return -np.log10(10.**(-alp*x)+1)+delta*(np.log10(1.+np.exp(x)))**gamma/( 1. + np.exp(10.**-x) )

    def SHMRbeh(self, xmh,zz):
        nu, alp, delta, gamma, eps, xM1 = self.evolved_factors(zz)
<<<<<<< HEAD
        print(self.evolved_factors(zz))
=======
        print self.evolved_factors(zz)
>>>>>>> 2cb08947c4825a1e0e38de69838ff8bfaff35728

        xmstel=eps+xM1+self.fbeh12(xmh-xM1, nu, alp, delta, gamma)-self.fbeh12(0.0, nu, alp, delta, gamma);
        return xmstel;

#def sfe(logM,zref):
#    # (1+z), Log10(Halo Mass), Log10(SFR)
#    data = pd.read_csv('release-sfh_z0_z8_052913/sfe/sfe.dat',delimiter=' ')
#    z = data.values[:,0]-1
#    ids = (z<=5)
#    z = z[ids]
#    log_M = data.values[ids,1]
#    sfe = data.values[ids,2]
#
#    srt = np.argmin((z - zref)**2)
#    suck = (z==z[srt])
#
#    srt_logM,srt_sfe = log_M[suck],sfe[suck]
#
#    fun = interpolate.interp1d(srt_logM,srt_sfe)
#    result = 10**(fun(logM))
#
#    if np.isnan(result):
#        result=0.0
#
#    return result
#

if __name__ == "__main__":
    bh = Behroozi_2013a()
    mh = np.linspace(9,16)
    plt.plot(mh, bh.SHMRbeh(mh-np.log10(0.73),0.24)+2*np.log10(0.73))

    plt.savefig('test.png', dpi=300)

