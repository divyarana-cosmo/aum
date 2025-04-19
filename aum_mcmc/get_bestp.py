import numpy as np
import pandas as pd
import sys

chnfil = sys.argv[1]
predfil = sys.argv[2]

pred = pd.read_csv(predfil, delim_whitespace=1, header=None)
idx = np.argmin(pred.values[:,-1])
<<<<<<< HEAD
print(np.min(pred.values[:,-1]))
chn = pd.read_csv(chnfil, delim_whitespace=1, header=None)

print(chn.values[idx,:-2])
=======
print np.min(pred.values[:,-1])
chn = pd.read_csv(chnfil, delim_whitespace=1, header=None)

print chn.values[idx,:-2]
>>>>>>> 2cb08947c4825a1e0e38de69838ff8bfaff35728
