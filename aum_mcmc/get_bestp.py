import numpy as np
import pandas as pd
import sys

chnfil = sys.argv[1]
predfil = sys.argv[2]

pred = pd.read_csv(predfil, delim_whitespace=1, header=None)
idx = np.argmin(pred.values[:,-1])
print(np.min(pred.values[:,-1]))
chn = pd.read_csv(chnfil, delim_whitespace=1, header=None)

print(chn.values[idx,:-2])
