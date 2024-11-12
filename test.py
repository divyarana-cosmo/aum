import sys
sys.path.append('/Users/divyarana/github/aum/install/lib/python3.13/site-packages/')
import cosmology as cc

# This is the default constructor with some basic cosmological parameters
a=cc.cosmology()
# Prints out the comoving distance in the fiducial cosmology
print(a.Dcofz(2.0))

# Prints the abundance of 1e9 Msun halos at z=0.0
print(a.nofm(1e9,0.0))

# Print all functions
#help(a)
