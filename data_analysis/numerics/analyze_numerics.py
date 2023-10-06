"""
Analyze the output data of the BlumeCappelKagome code to provide data
for the "numerics" part of Fig. 3a.

Author: Alexis Poncet <alexis.poncet@ens-lyon.fr>
"""

import numpy as np
import pandas as pd
from scipy.stats import binned_statistic

dir_numerics = "data_numerics/"
dir_analyze = "data_numerics/"
N = 30
K = 2
skip = 10000
Js = [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]
Jcs = [-4*J/K for J in Js]

data_len = []
data_Rg = []
data_nloops = []

for J in Js:
    fname = dir_numerics + f"BCK_N{N}_K{K}_J{J}_skip{skip}_len.dat.gz"
    data_len.append(pd.read_csv(fname, dtype=int, skiprows=2).values[:, 0])
    print(fname, "loaded")
    fname = dir_numerics + f"BCK_N{N}_K{K}_J{J}_skip{skip}_Rg.dat.gz"
    data_Rg.append(pd.read_csv(fname, skiprows=2).values[:, 0])
    print(fname, "loaded")
    fname = dir_numerics + f"BCK_N{N}_K{K}_J{J}_skip{skip}_nloops.dat.gz"
    print("Loading", fname)
    data_nloops.append(pd.read_csv(fname, skiprows=2, delimiter=" ").values)

meanPara = [np.mean(d[:, 3]) for d in data_nloops]
meanHeight = [np.mean(d[:, 4]) for d in data_nloops]
bins_l = np.logspace(0.75, 3.25, 50)
bin_centers_l = np.sqrt(bins_l[1:] * bins_l[:-1])
bin_centers_log_num = []
means_logRg_num = []
std_logRg_num = []

for i in range(len(Js)):
    print(i, end='\r')
    cRg, _, _ = binned_statistic(data_len[i], np.log(data_Rg[i]),
                                 statistic='count', bins=bins_l)
    m = (cRg > 1)
    bin_centers_log_num.append(np.log(bin_centers_l[m]))
    mRg, _, _ = binned_statistic(data_len[i], np.log(data_Rg[i]),
                                 statistic='mean', bins=bins_l)
    means_logRg_num.append(mRg[m])
    sRg, _, _ = binned_statistic(data_len[i], np.log(data_Rg[i]),
                                 statistic='std', bins=bins_l)
    std_logRg_num.append(sRg[m] / np.sqrt(cRg)[m])

expRg2 = []
prefRg2 = []
err2 = []
minBins = [15, 15, 15, 15, 20, 20, 20, 20, 20]
maxBins = [2e2, 2e2, 2e2, 2e2, 4e2, 5e2, 5e2, 5e2, 5e2]

for i, mRg in enumerate(means_logRg_num):
    print(i)
    m = (
            (np.exp(bin_centers_log_num[i]) > minBins[i])
            & (np.exp(bin_centers_log_num[i]) < maxBins[i])
        )
    (A, B) = np.polyfit(bin_centers_log_num[i][m], mRg[m], 1)
    expRg2.append(A)
    prefRg2.append(np.exp(B))


skip = 1000
data_correl = []
for J in Js:
    fname = dir_numerics + f"correl/BCK_N{N}_K{K}_J{J}_skip{skip}_correl.dat.gz"
    print("Loading", fname)
    data_correl.append(pd.read_csv(fname, skiprows=2, delimiter=" ").values)

# Num Rg
for i in reversed(range(len(Js))):
    dd = np.column_stack((np.exp(bin_centers_log_num[i]),
                          np.exp(means_logRg_num[i])))
    np.savetxt(dir_analyze + f"gyration_radius_Jc{Jcs[i]}.txt", dd,
               header="L Rg")

# Num nu
dd = np.column_stack((Jcs, expRg2))
np.savetxt(dir_analyze + "exponent_nu.txt", dd, header="Jc nu")

# Num height
dd = np.column_stack((Jcs, meanHeight))
np.savetxt(dir_analyze + "nesting_level.txt", dd, header="Jc h")

# Num correl
for i in range(len(Js)):
    q = len(Js)-i-1
    m = data_correl[q][:, 1] > 0
    dd = np.column_stack((data_correl[q][m, 0],
                          data_correl[q][m, 1]/data_correl[q][m, 2]/(1000*50)))
    np.savetxt(dir_analyze + f"pair_correlations_Jc{Jcs[i]}.txt", dd,
               header="r C")
